#!/Users/privateprivate/envs/bin/python

import networkx as nx
from geopy.distance import lonlat, distance
from sklearn.metrics.pairwise import euclidean_distances
from osgeo import gdal, ogr, osr,gdalconst,gdal_array
import os
import glob
import matplotlib.pyplot as plt
import time
import pathlib 
from pathlib import Path
import numpy as np
import sys
import geopandas as gpd
import pandas as pd
from shapely.wkt import loads
from shapely.geometry import shape, mapping,LineString, Point,MultiLineString,MultiPoint
from shapely.ops import nearest_points
import fiona
import math
import csv
import shapelytools as shptools
import progress as progress

gdal.UseExceptions()
gdal.PushErrorHandler('CPLQuietErrorHandler')
gdal.SetConfigOption("GDAL_DATA", "gdal-data/")	

'''
Saeon Bioenergy Atlas model 
@ Main authors: Marc Pienaar (marc@saeon.ac.za, marc.pienaar@gmail.com) and Hayden Wilson (hayden@saeon.ac.za)
'''

#Main entry function for performing the clustering and routing algorithms
def run(inshape,attrs1,attrs1_ids,feedstock,technology,facility_start_point,capacity1,outputfile,roads,inpoly,buffer):
	'''
	@ Author Marc Pienaar 
	1) Iteratively clusters an input centroid layer using the max volume and it's closest neighbours using the generate facilities function and a facility threshold
	2) Determines if the new Facility location and its feedstock locations (clusstered centroids) need to snap to the closest road or closet virtual road, in which case it maps the distance from a virtual road path to the closest road.
	3) Clips various inputs files to the extent of the current facility catchment (these include roads, irtual road start and end points, and proximity rasters for the Virtual roads and main roads 
	4) add the new points that are generated into the road clip to perform routing
	5) fetches the polygoin areas associated with each facility feedstock location
	6) genereates geopackage outputs of the facility, its feeedstock locations and the catchement area
	7) Performs a one to many routing using NetworkX Dijkstra algorithm
	8) writes out the paths layer and creates a technology csv file with all the paramters (distances and areas etc)
	9) Processes the technology csv into 'real' costs. 

	Parameters
	----------
	inshape : String
		Filename of the input centroid layer to read.
		
	attrs1 : Dictionary of attribute names from the inshape file			
		
	attrs1_ids : Array[]
		A list of indices to use listed the attrs1 dictionary [id, area, volume in total tonnes]
	
	feedstock :  String
			The name of the feedstock e.g. Alien invasives

	technology : String
		The fname of the current technology e.g. CBP-Chipping
	
	facility_start_point : int
		The current facility number to start with (1 if starting a new set of scenarios)
	
	capacity1 : int 
		The max capacity the facility to process feedstock volumes (in tonnes)
		
	outputfile : String
		Filename of the csv file to update.
			
	roads : String
		Filename of the vector road layer too use for routing. 
	
	inpoly : String
		Filename of the polygon layer with the same id values as the inshape centroid layer in order to caluclate and generate a catcment for each facility. 

	buffer : real / double
		a buffer amount in degrees (WGS84) to add to the clipping area

	Returns
	-------
	
	1) Geopackage outputs of the facilities, their feeedstock locations catchement area, and paths if routing occured
	2) A csv output oof all the paramters (faciity id, feeedstocks, percent filled and distances need to travel according to road type)
	3) the last facility generated as a return statement.
	'''
	start_time = time.time()
	#Create a road features array to mmap distances travelled to
	features=['Arterial Road (m)','Interchange (m)', 'Main Road (m)', 'National Freeway (m)','National Road (m)','None (m)','On/OffRamp (m)','Other Road (m)','Secondary Road (m)','Street (m)','Virtual road (m)']
	feature_values = [0] * len(features)#zero initialise the values 
	#Declare some output files
	VR_vertices_in_wgs84_temp=os.path.sep.join(['outputs','temp', 'start_vertices.gpkg'])
	VR_vertices_out_wgs84_temp=os.path.sep.join(['outputs','temp', 'end_vertices.gpkg'])
	roads_temp=os.path.sep.join(['outputs', 'temp','temp_roads.shp'])
	roads_temp2=os.path.sep.join(['outputs', 'temp','temp_roads2.shp'])
	outShapefileextent=os.path.sep.join(['outputs', 'temp','extent.gpkg'])
	#Declare some input data 	
	VR_vertices_in_wgs84=os.path.sep.join(['data', 'start_vertices.shp'])
	VR_vertices_end_wgs84=os.path.sep.join(['data', 'end_vertices.shp'])
	VR_proximity=os.path.sep.join(['data', 'VR_roads_100m_proximity.tif'])
	ngi_proximity=os.path.sep.join(['data', 'NGI_roads_100m_proximity.tif'])
	#get all the start points, ids and vlaues as a dictionary object
	start_points=get_start_point_Dict(attrs1,attrs1_ids,inshape,'GPKG')
	start_points_new2=start_points.copy()#the object to itterate over until the volume has been depleted
	#this is how many runs are required
	valll=np.sum(start_points_new2['values'])/capacity1	
	#create a directory to store outputs from the current run
	outputdir=os.path.sep.join(['outputs', technology+"_"+str(capacity1)])
	mkdir_recursive(outputdir)		
	count=-100
	endpointnumber=facility_start_point
	fsp_temp=1
	while len(start_points_new2['id']) >0:
		endpoint_id_temp=[]
#	for kk in range(0,265): #in case we want to test a certain facility
		#perform a sort [max to min], then Euclidean distance from max to other points and sort again		
		start_points_new2= sorting(start_points_new2)
		fix_line_ends=False
		#generate facilities
		facilities=generate_facilities(capacity1,start_points_new2,facility_start_point)
		#Adjust the start point arrays	
		start_points_new2['id']=facilities['id']
		start_points_new2['points']=facilities['points']
		start_points_new2['values']=facilities['values']
		#initiate some variables
		points_temp=[];values_temp=[];id_temp=[];iterr=0;sp_id=[]
		sp_id.append(facilities['start_point_id'][0])
		tons=facilities['tonnes contributed'][0]
		#increment the facility_start_point by 1	
		count=max(len(start_points_new2['id']),count)
		term=[];term.append('Run ');term.append(str(fsp_temp) + " of ");term.append(str(math.ceil(valll)) + ': ')
		percent= round((fsp_temp/valll)*100,3)
		if percent >100:
			percent=100			
		term.append("%0.2f" % percent);	term.append('% complete: ')
		term.append(str(len(facilities['id'])));term.append(' points remaining: ')
		term.append('processing volume from ');term.append(str(len(sp_id[0])));term.append(' point(s): ')
		#increment the facility start point number
		facility_start_point=facility_start_point+len(facilities['End_point_id'])
		fsp_temp=fsp_temp+len(facilities['End_point_id'])
		#adjust the time printout
		runtime=(time.time() - start_time);runtime2=(time.time() - start_time)
		if runtime < 60:
			term.append("%0.2f" % runtime);term.append(' seconds runtime')
		elif runtime >= 60 and runtime < 3600 : 
			runtime2=runtime2/60
			term.append("%0.2f" % runtime2);term.append(' minutes runtime')
		elif runtime >= 3600 : 
			runtime2=runtime2/3600
			term.append("%0.2f" % runtime2);term.append(' hours runtime')
		print(''.join(term))
		#create some shapefiles to write values to. i.e start points, end points, paths, aqnd catchement 
		outpaths,outpaths_sp_original,outpaths_sp,outpaths_ep,outpaths_area=\
				create_output_shapefiles(outputdir,facilities,capacity1)
		#initialise the Dict object
		Dict={}		
		#check for VR or national road proximty and return the results to the proximity variable			
		proximity=Check_VR_or_NGI_roads_proximity(VR_proximity,ngi_proximity,facilities['start_points'][0],facilities['start_point_id'][0],False)
		#combine all the points
		combined_points=get_combined_points(facilities['start_points'][0],facilities['start_point_id'][0],start_points_new2['id'],facilities['end_points'][0])
		#Creat a temporary extent to use for cropping various layers
		create_extent("GPKG",combined_points['combined_points'],inshape,outShapefileextent,buffer)#
		#####crop the roads layer, and VR roads centtroid layers (faster than doing a Vr mapping on the bigger files
		clip_vector_by_shape('ESRI Shapefile','GPKG','ESRI Shapefile',roads,outShapefileextent,roads_temp,"line")		
		clip_vector_by_shape('ESRI Shapefile','GPKG','GPKG',VR_vertices_in_wgs84,outShapefileextent,VR_vertices_in_wgs84_temp,"point")
		clip_vector_by_shape('ESRI Shapefile','GPKG','GPKG',VR_vertices_end_wgs84,outShapefileextent,VR_vertices_out_wgs84_temp,"point")#		
		#adjust and remap any virtual roads points, a longish function as it has to go through the whole of SA
		start_points_temp,ids_temp,vr_distance= adjust_points('GPKG',proximity,VR_vertices_in_wgs84_temp,VR_vertices_out_wgs84_temp,False)
		#re-combine all the points using the adjused values
		combined_points=get_combined_points(start_points_temp,ids_temp,start_points_new2['id'],facilities['end_points'][0])
		#make sure the new values are on the road, if not snap thm to the road		
		start,end,multi_point= get_road_aligned_points(roads_temp,combined_points)
		#create some shapefiles
		create_new_point_layer2('GPKG', inshape, outpaths_sp,start,'point_ids',sp_id[0],technology,capacity1)
		create_new_point_layer2('GPKG', inshape, outpaths_sp_original,facilities['start_points'][0],"point_ids",sp_id[0],technology,capacity1)		
		end2=[];end2.append(facilities['End_point_id'][iterr])
		create_new_point_layer2('GPKG', inshape, outpaths_ep,end,'point_ids',end2,technology,capacity1)
		area=create_field_scenario('ESRI Shapefile','GPKG',facilities['start_point_id'][0],attrs1,attrs1_ids,inpoly,outpaths_area,attrs1['attributes'][attrs1_ids[1]][1],technology,capacity1)		
		if len(facilities['start_point_id'][0])<=1:
			Dict = create_output_dictionary_single(Dict,feature_values,features,vr_distance,facilities,capacity1,area)
		else:
			#update start points
			updatenode(roads_temp,outpaths_sp,roads_temp)
			#update end point
			updatenode(roads_temp,outpaths_ep,roads_temp)
			exlusions=['bad road','Dirt Road']
			#do the routing
			a=routing(roads_temp,start,end,outpaths,False,'w2',False,'FEAT_TYPE',\
			exlusions,1,True,sp_id[0],facilities['End_point_id'][iterr],tons,area,technology,capacity1,vr_distance)
			Dict=create_output_dictionary_paths(Dict,facilities,capacity1,area,a)
		#Update the csv file, then we are done!
		update_csv(outputfile,Dict,feedstock,technology)
		endpointnumber=Dict['End_point_id']
	return(endpointnumber)



def print_runtime(start_time):
	'''
	@ Author Marc Pienaar
	
	A small function to print the runtime in seconds, minutes or hours

	Parameters
	----------
	start_time : time
		The start time of a proccess

	Returns
	-------
	a String representation of the runtime taken 
	
	'''
	runtime=(time.time() - start_time);runtime2=(time.time() - start_time)
	term=[]
	if runtime < 60:
		term.append("%0.2f" % runtime);term.append(' seconds total runtime')
	elif runtime >= 60 and runtime < 3600 : 
		runtime2=runtime2/60
		term.append("%0.2f" % runtime2);term.append(' minutes total runtime')
	elif runtime >= 3600 : 
		runtime2=runtime2/3600
		term.append("%0.2f" % runtime2);term.append(' hours total runtime')
	return(''.join(term))
	
def rmdir(directory):
	'''
	A small pathlib function recursively remove a directory

	Parameters
	----------
	directory : String
		The directory path to remove
	'''
	
	directory = Path(directory)
	for item in directory.iterdir():
		if item.is_dir():
			rmdir(item)
		else:
			item.unlink()
	directory.rmdir()
	


def create_outputs(technology,feedstock,start_new,create_file):
	'''
	@ Author Marc Pienaar
	A small utility function to create the ourput directory, and csv output file 

	Parameters
	----------
	technology : String
		The fname of the current technology e.g. CBP-Chipping
	
	feedstock :  String
			The name of the feedstock e.g. Alien invasives	
	
	start_new: bool
		if True: first tries to remove the exisiting directories
		if False: does nothing

	create_file: bool
		if True: Creates a new csv file
		if False: does nothing
	
	Returns
	-------
	the name of the output csv file 
	'''
	
	#Make an outputs folder for temporary files in the working directory
	if start_new:
		try:
			rmdir(Path(os.path.sep.join(['outputs'])))
			rmdir(Path(os.path.sep.join(['Processed'])))
		except:
			pass
	mkdir_recursive(os.path.sep.join(['outputs', 'temp']))
	mkdir_recursive(os.path.sep.join(['Processed']))
	outputfile=os.path.sep.join(['outputs', technology+'_'+feedstock+'.csv'])
	if create_file:
		with open(outputfile, 'w', newline='') as file:
			writer = csv.writer(file)
			writer.writerow(['hits', 'feedstock','technology','end_point_id', 'X_coords','Y_coords','capacity','percent_full','original_feedstock_tonnes','tonnes_contributed','tonnes_still_needed','remaining_feedstock_tonnes','start_point_id','service_area_(ha)','Arterial_Road_(m)', 'Interchange_(m)', 'Main_Road_(m)', 'National_Freeway_(m)', 'National_Road_(m)', 'None_(m)', 'On/OffRamp_(m)', 'Other_Road_(m)', 'Secondary_Road_(m)', 'Street_(m)', 'Virtual_road_(m)'])
	
	return outputfile


def sorting(start_points_new2):
	'''
	@ Author Marc Pienaar 
	1) Sorts the start_points_new2 values [ids,vlaues,coordinates] from largest to smalles (according to value), the reorders according to the nearest neighbours (coordinates) of the max value

	Parameters
	----------
	start_points_new2 : Dictionary object with 3 values [ids, values, coords]		

	Returns
	-------
	
	1) A sorted Dictionary object from max value to min, then ordered according to the euclidean distance from the max to its closests neighbours .
	'''
	spv=np.array(start_points_new2['values'])
	ids = np.flip(np.argsort(spv))
	spid2=[];spv2=[];spp2=[]
	for i in ids:
		spid2.append(start_points_new2['id'][i])
		spv2.append(start_points_new2['values'][i])
		spp2.append(start_points_new2['points'][i])
	start_points_new2['id']=spid2
	start_points_new2['values']=spv2
	start_points_new2['points']=spp2
	#now perform an Euclidean distance
	sps1=start_points_new2['points'];sps=sps1.copy()
	X=np.asarray(sps) 
	xx=np.asarray(sps[0]).reshape(1, -1)
	dist=euclidean_distances(xx,X)
	dist2 = np.argsort(dist.flatten())
	spid2=[];spv2=[];spp2=[]
	for i in dist2:
		spid2.append(start_points_new2['id'][i])
		spv2.append(start_points_new2['values'][i])
		spp2.append(start_points_new2['points'][i])
	start_points_new2['id']=spid2
	start_points_new2['values']=spv2
	start_points_new2['points']=spp2
	return start_points_new2

def fix_lines(shpin,tempshp,drivername):
	'''
	@ Author Marc Pienaar 
	A function to attempt to fix line endings that don't join in a 5km buffer in the road layer (uses the library from https://github.com/ojdo/python-tools/blob/master/shapelytools.py to do the line snapping)
	This fucntion is currently not implemented in this test example

	Parameters
	----------
	shpin : String
		The name of the road layer to try and fix

	tempshp : String
		The name of a temporary shapefile layer to uses to copy values to and from 

	drivername:  String
		Ogr driver name to handle the file e.g. 'ESRI_Shapefile' or 'GPKG'

	Returns
	-------
	An updated shpin file with fixed line endings if there where any broken within a 5km buffer
	'''	
	a,b,c=shptools.read_shp(shpin)
	# use a 5km buffer e.i 0.00005
	d=shptools.snappy_endings(a,0.00005)
	shptools.write_shp(tempshp, d, records=b, fields=c)
	spatialRef = osr.SpatialReference()
	spatialRef.ImportFromEPSG(4326)#wgs84
	spatialRef.MorphToESRI()
	file = open(os.path.splitext(tempshp)[0]+'.prj', 'w')
	file.write(spatialRef.ExportToWkt())
	file.close()
	#now rewrite to shpin
	driver = ogr.GetDriverByName(drivername)
	inDataSet = driver.Open(tempshp)
	inLayer = inDataSet.GetLayer()
	srs = inLayer.GetSpatialRef()
	outDataSet = driver.CreateDataSource(shpin)
	outLayer = outDataSet.CreateLayer('roads',srs,geom_type=ogr.wkbMultiLineString)
	featureCount = inLayer.GetFeatureCount()
	inLayerDefn = inLayer.GetLayerDefn()
	for i in range(0, inLayerDefn.GetFieldCount()):
		fieldDefn = inLayerDefn.GetFieldDefn(i)
		outLayer.CreateField(fieldDefn)
	# get the output layer's feature definition
	outLayerDefn = outLayer.GetLayerDefn()
	for feature in inLayer:
		geom = feature.GetGeometryRef()
		outFeature = ogr.Feature(outLayerDefn)
		outFeature.SetGeometry(geom)
		for i in range(0, outLayerDefn.GetFieldCount()):
			outFeature.SetField(inLayerDefn.GetFieldDefn(i).GetNameRef(), feature.GetField(i))
		outLayer.CreateFeature(outFeature)
		# dereference the features and get the next input feature
		outFeature = None
	inDataSet = None
	outDataSet=None
	

def write_geometry(geometry, out_path, srs_id=4326):
	'''
	Saves the geometry in an ogr.Geometry object to a shapefile.

	Parameters
	----------
	geometry
		An ogr.Geometry object
	out_path
		The location to save the output shapefile
	srs_id
		The projection of the output shapefile. Can be an EPSG number or a WKT string.

	Notes
	-----
	The shapefile consists of one layer named 'geometry'.

	'''
	# TODO: Fix this needing an extra filepath on the end
	driver = ogr.GetDriverByName("ESRI Shapefile")
	data_source = driver.CreateDataSource(out_path)
	srs = osr.SpatialReference()
	if type(srs_id) is int:
		srs.ImportFromEPSG(srs_id)
	if type(srs_id) is str:
		srs.ImportFromWkt(srs_id)
	layer = data_source.CreateLayer(
		"geometry",
		srs,
		geom_type=geometry.GetGeometryType())
	feature_def = layer.GetLayerDefn()
	feature = ogr.Feature(feature_def)
	feature.SetGeometry(geometry)
	layer.CreateFeature(feature)
	data_source.FlushCache()
	data_source = None 

def update_csv(outputfile,Dict,feedstock,technology):
	'''
	@ Author Marc Pienaar 
	A function to update the csv file for each technology
	
	Parameters
	----------
	outputfile : String
		The name of the rcsv file to update

	Dict : Dictionary object
		The dictionary object from the main run function 

	feedstock:  String
		The name of the current feedstock

	technology:  String
		The name of the current technology

	Returns
	-------
	An updated csv file (i.e outputfile)
	'''
	with open(outputfile, 'a', newline='') as file:
		writer = csv.writer(file)
		for i in range(0,Dict['hits']):
			outout=[]	
			for j in range(0,len(Dict['start_point_id'])):
				outout.append(Dict['hits'])
				outout.append(feedstock)
				outout.append(technology)
				outout.append(Dict['End_point_id'])
				for k in Dict['end_point']:
					outout.append(k[0])
					outout.append(k[1])
				outout.append(Dict['Capacity'])
				for k in Dict['percent_full']:
					outout.append(k[i])
				for k in Dict['original_feedstock_tonnes']:
					outout.append(k[i])
				for k in Dict['tonnes_contributed']:
					outout.append(k[i])
				for k in Dict['tonnes_still_needed']:
					outout.append(k[i])
				for k in Dict['remaining_feedstock_tonnes']:
					outout.append(k[i])
				for k in Dict['start_point_id']:
					outout.append(k[i])
				for k in Dict['service_area_(ha)']:
					outout.append(k[i])
				for k in Dict['feature_type_distance'][i]:
					outout.append(k)
				writer.writerow(outout)


def create_output_dictionary_paths(Dict,facilities,capacity1,area,a):
	'''
	@ Author Marc Pienaar 
	A utility function to update the Dictionary from the main run function
	
	Parameters
	----------
	Dict : Dictionary object
		The dictionary object from the main run function 

	facilities:  Dictionary
		The factilities dictionary

	capacity1:  int
		the capacity of the current technology
	
	area:  Real or Double 
		the feedstock catchement area of the current facility

	a : Dictionary object
		The Dictionary output from the road routing with distancees per road feature type

	Returns
	-------
	Dict - the updated dictionary object
	'''
	Dict['hits']=a['hits']
	Dict['End_point_id']=facilities['End_point_id'][0]
	temp=[];temp.append(facilities['end_points'][0])
	Dict['end_point']=temp
	Dict['Capacity']=capacity1
	temp=[];temp.append(facilities['percent facility full'][0])
	Dict['percent_full']=temp
	temp=[];temp.append(facilities['original values'][0])
	Dict['original_feedstock_tonnes']=temp
	temp=[];temp.append(facilities['tonnes contributed'][0])
	Dict['tonnes_contributed']=temp
	temp=[];temp.append(facilities['tonnes still needed'][0])
	Dict['tonnes_still_needed']=temp
	temp=[];temp.append(facilities['excess'][0])
	Dict['remaining_feedstock_tonnes']=temp
	temp=[];temp.append(facilities['start_point_id'][0])
	Dict['start_point_id']=temp
	temp=[];temp.append(area)
	Dict['service_area_(ha)']=temp
	Dict['feature_type']=a['feature_type']
	Dict['feature_type_distance']=a['feature_type_distance']
	return Dict

def create_output_dictionary_single(Dict,feature_values,features,vr_distance,facilities,capacity1,area):
	'''
	@ Author Marc Pienaar 
	A utility function to update the Dictionary from the main run function without routing information
	
	Parameters
	----------
	Dict : Dictionary object
		The dictionary object from the main run function 

	feature_values :  String Array[]
		An array of feature names in the road

	features :  Array[]
		An array of distance values (zero, unless there was a virtual road travelled
	
	vr_distance:  Real or Double 
		The distance travlled along a virtual road (could be zero if no VR was travelled)

	facilities:  Dictionary
		The factilities dictionary

	capacity1:  int
		the capacity of the current technology
	
	area:  Real or Double 
		the feedstock catchement area of the current facility

	Returns
	-------
	Dict - the updated dictionary object
	'''
	feature_values[len(feature_values)-1]=vr_distance[0]
	Dict['hits']=1
	Dict['End_point_id']=facilities['End_point_id'][0]
	temp=[];temp.append(facilities['end_points'][0])
	Dict['end_point']=temp
	Dict['Capacity']=capacity1
	temp=[];temp.append(facilities['percent facility full'][0])
	Dict['percent_full']=temp
	temp=[];temp.append(facilities['original values'][0])
	Dict['original_feedstock_tonnes']=temp
	temp=[];temp.append(facilities['tonnes contributed'][0])
	Dict['tonnes_contributed']=temp
	temp=[];temp.append(facilities['tonnes still needed'][0])
	Dict['tonnes_still_needed']=temp
	temp=[];temp.append(facilities['excess'][0])
	Dict['remaining_feedstock_tonnes']=temp
	temp=[];temp.append(facilities['start_point_id'][0])
	Dict['start_point_id']=temp
	temp=[];temp.append(area)
	Dict['service_area_(ha)']=temp
	temp=[];temp.append(features)
	Dict['feature_type']=temp
	temp=[];temp.append(feature_values)
	Dict['feature_type_distance']=temp
	return Dict

def get_stats_of_attr(drivername,shp_in,attr):
	'''
	@ Author Marc Pienaar
	Get basic stats of an attribute from a shapefile
	Parameters
	----------
	drivername:  String
		Ogr driver name to handle the file e.g. ESRI_Shapefile or GPKG
		
	shp_in : File or string
		File, directory, or filename to read.

	attr : String
		name of a numeric-based attribute
	
	Returns
	-------
	Dictionary object:  'min','1st quantile','median','3rd quantile','max','mean','sum'and 'stdev'		
	'''
	start_time = time.time()
	shpDriver = ogr.GetDriverByName(drivername)
	source_ds = ogr.Open(shp_in)  # open the original shp
	layer = source_ds.GetLayer()
	featureCount = layer.GetFeatureCount()
	Dict={}
	temp=[]
	for feature in layer:
		temp.append(feature.GetField(attr))
	source_ds=None
	try:
		values=np.asarray(temp)
		Dict['min']=min(values)
		Dict['1st quantile']=np.quantile(values, .25)
		Dict['median']=np.median(values)
		Dict['3rd quantile']=np.quantile(values, .75)
		Dict['max']=max(values)
		Dict['mean']=np.mean(values)
		Dict['sum']=np.sum(values)
		Dict['stdev']=np.std(values)
	except:
		Dict=[]
	#print("Attribute statistics took --- %s seconds ---" % (time.time() - start_time))
	return(Dict)

def get_attributes(drivername, shp_in):
	'''
	@ Author Marc Pienaar
	Get the field names and feature count from a shapefile
	Parameters
	----------
	drivername:  String
		Ogr driver name to handle the file e.g. ESRI_Shapefile or GPKG
		
	shp_in : file or string
		File, directory, or filename to read.
	
	Returns
	-------
	Dictionary object:  'attributes' and 'number of features'		
	'''
	start_time = time.time()
	attrs=[]
	Dict={}
	shpDriver = ogr.GetDriverByName(drivername)
	source_ds = ogr.Open(shp_in)  # open the resource
	layer = source_ds.GetLayer()
	featureCount = layer.GetFeatureCount()
	layerDefinition = layer.GetLayerDefn()
	accum=0	
	for i in range(layerDefinition.GetFieldCount()):
		temp=[]
		temp.append(accum)
		temp.append(layerDefinition.GetFieldDefn(i).GetName())
		accum+=1
		attrs.append(temp)
	Dict['attributes']=attrs
	Dict['number of features']=featureCount
	source_ds=None
	#print("Get attributes from file took --- %s seconds ---" % (time.time() - start_time))
	return(Dict)

def mkdir_recursive(path):
	'''
	create a nested directory
	Parameters
	----------
	path:  String
		directory, and/or sub directories to create
	Returns
	-------
	Null, Creates a nested directory if it doesn't exists		
	'''
	if not os.path.exists(path):
		pathlib.Path(path).mkdir(parents=True, exist_ok=True)
		
def get_centroids(inshp,outshp,drivername,drivername2,overwrite):
	
	'''
	@ Author Marc Pienaar
	create a centroid shapefile from a polygon shapefile 
	Parameters
	----------
	inshp:  String
		filename for the input shape

	outshp:  String
		filename for the output shape

	drivername:  String
		Driver name to write e.g. 'ESRI Shapefile', or 'GPKG'

	drivername2:  String
		Driver name to write e.g. 'ESRI Shapefile', or 'GPKG'

	overwrite: bool
		if True: overwrites an existing file
		if False: does nothing if file exists

	Returns
	-------
	a new output vector (point) file		
	'''
	
	start_time = time.time()
	def get_centroids2(inshp,outshp,drivername,drivername2):
		driver=ogr.GetDriverByName(drivername)
		driver2=ogr.GetDriverByName(drivername2)
		inDataSource = driver.Open(inshp, 0)	
		inLayer = inDataSource.GetLayer()
		srs = inLayer.GetSpatialRef()
		outDataSet = driver2.CreateDataSource(outshp)
		outLayer = outDataSet.CreateLayer('centroids',srs,geom_type=ogr.wkbPoint)
		featureCount = inLayer.GetFeatureCount()
		inLayerDefn = inLayer.GetLayerDefn()
		for i in range(0, inLayerDefn.GetFieldCount()):
			fieldDefn = inLayerDefn.GetFieldDefn(i)
			outLayer.CreateField(fieldDefn)
		#get the output layer's feature definition
		outLayerDefn = outLayer.GetLayerDefn()
		for feature in inLayer:
			geom = feature.GetGeometryRef()
			geom2=geom.Centroid()
			outFeature = ogr.Feature(outLayerDefn)
			for j in range(0, inLayerDefn.GetFieldCount()):
				outFeature.SetField(inLayerDefn.GetFieldDefn(j).GetNameRef(), feature.GetField(j))
			outFeature.SetGeometry(geom2)
			outLayer.CreateFeature(outFeature)
			# dereference the features and get the next input feature
			outFeature = None
		inDataSet = None
		outDataSet=None
		
	if os.path.exists(outshp):
		if overwrite:
			get_centroids2(inshp,outshp,drivername,drivername2)
		else:
			print(outshp, "already exists, moving on")
	else:
		get_centroids2(inshp,outshp,drivername,drivername2)
	print("creation of centroid shp file --- %s seconds ---" % (time.time() - start_time))
	
def get_points_tuple(drivername,shp_in,attrs,N):
	'''
	@ Author Marc Pienaar
	Get points and other attributes
	Parameters
	----------
	drivername:  String
		Ogr driver name to handle the file e.g. ESRI_Shapefile or GPKG
		
	shp_in : File or string
		File, directory, or filename to read.

	attrs : String list[]
		names of field attributes

	N : Numeric
		number of attributes in each field
	
	Returns
	-------
	Dictionary object:  'coordinates as piont tupple','attribute values array'		
	'''
	start_time = time.time()
	points2=[]
	shpDriver = ogr.GetDriverByName(drivername)
	source_ds = ogr.Open(shp_in)  # open the original shp
	layer = source_ds.GetLayer()
	featureCount = layer.GetFeatureCount()
	Dict={};ids=[];	vals=[];arr = []
	rows, cols = (N,len(attrs))	
	accum=0
	for feature in layer:
		points=[]
		point = feature.GetGeometryRef()
		points.append(point.GetX())
		points.append(point.GetY())
		aa=tuple(points)
		points2.append(aa)
		temp=[]
		for i in range(0,cols):
			temp.append(feature.GetField(attrs[i]))
		arr.append(temp)
	source_ds=None
	for j in range(0,cols):
		temp=[]
		for i in range(0,rows):
			temp.append(arr[i][j])
		Dict[attrs[j]]=temp	
	Dict['points']=points2
	#print("point creation took --- %s seconds ---" % (time.time() - start_time))
	return(Dict)

def get_start_point_Dict(attrs1,attrs1_ids,inshape,drivername):
	'''
	@ Author Marc Pienaar
	get a dictionary of ids,values,and coords with all zero values removed
	Parameters
	----------
	attrs1:  Dictionary
		dictionary of attribute names
		
	attrs1_ids : [] array
		the index of the the attributes to extract

	inshape : String
		name of the shapefile to extract values from

	drivername:  String
		Ogr driver name to handle the file e.g. ESRI_Shapefile or GPKG
	
	Returns
	-------
	Dictionary object:  'ids', 'values', and 'points'		
	'''
	attributes=[]
	for i in attrs1_ids:
		attributes.append(attrs1['attributes'][i][1])#	
	####### Get the start points, IDS, and some attribute fields for weighting ###
	start_points2=get_points_tuple(drivername,inshape,attributes,attrs1['number of features'])
	start_points={}
	start_points['id']=start_points2[attrs1['attributes'][attrs1_ids[0]][1]]
	start_points['values']=start_points2[attrs1['attributes'][attrs1_ids[len(attrs1_ids)-1]][1]]
	start_points['points']=start_points2['points']
	#do a filter to get rid of all values that are zero or less
	unwanted=[]
	for i in range(0,len(start_points['values'])):
		if start_points['values'][i]<=0:			
			unwanted.append(i)
	for ele in sorted(unwanted, reverse = True):  
		del start_points['values'][ele]
		del start_points['id'][ele]
		del start_points['points'][ele]
		
	spv=np.array(start_points['values'])
	ids = np.flip(np.argsort(spv))
	spid2=[]
	spv2=[]
	spp2=[]
	
	for i in ids:
		spid2.append(start_points['id'][i])
		spv2.append(start_points['values'][i])
		spp2.append(start_points['points'][i])

	start_points2={}
	start_points2['id']=spid2
	start_points2['values']=spv2
	start_points2['points']=spp2
	return start_points


def generate_facilities(threhold,start_points,facility_startnumer):
	'''
	@ Author Marc Pienaar
	Generates a clustring of propossed facilitis realtive to volumes in a collection of input points
	Parameters
	----------
	threhold:  number
		the maximum threhold for a processing facility
		
	start_points : dictionary object {}
		IDs, values, coordinates

	facility_startnumer : int
		the facility number
	
	Returns
	-------
	Dictionary object:  values, percentages and assocaited input ids for each facility		
	'''
	sps1=start_points['points'];sps=sps1.copy()
	sps_ids1=start_points['id'];sps_ids=sps_ids1.copy()
	sps_vals1=start_points['values'];sps_vals=sps_vals1.copy()
	tempperecent=[];t_contributed=0
	threshold=threhold
	
	temp_tons_contributed=[];temp_tons_still_needed=[];tons_still_needed=[];tons_contributed=[]
	tempexcess=[];excess=[];excess=[]
	temppercentfull=[];perecent=[];perecent2=[];percentfull=[];tempperecent2=[]
	tempexcess2=[];excess2=[]
	
	temp_startpointsnew=[];startpointsnew=[]
	ep=[];ep_ids=[];sp_ids=[];sps_ids_etmp=[];points1=[];ids1=[]
	accumid=facility_startnumer	
	weightss=[]
	Dict={}	
	
	start_p = sps[0]
	X=np.asarray(sps) 
	xx=np.asarray(start_p).reshape(1, -1)
	dist=euclidean_distances(xx,X)
	dist2 = np.argsort(dist.flatten())#a sorted array of euclidean distancees from the start point coord with the max value
	dist2=range(0,len(sps))
	new_threhold=threshold
	summm=0
	for i in range(0,len(dist2)):
		#the values associated with the closet points
		t2=sps_vals[dist2[i]]
		a=(new_threhold/t2)*100	
		#tempperecent is the % of t2 that was contributed or needs to be contributed to match the threshold
		#t_contributed is how many units from t2, were contributed to the threshold
		if a>= 100:
			tempperecent.append((new_threhold/t2)*100)
			t_contributed=round(1*t2,3)
		else:
			tempperecent.append((new_threhold/t2)*100)
			t_contributed=round((new_threhold/t2)*t2,3)
			
		summm=summm+t_contributed
		temp_tons_contributed.append(t_contributed)# an array to use outsid the loop
		tempexcess.append(t2-new_threhold)#the excess from t2, if positive there is still something to use
		ids1.append(sps_ids[dist2[i]])
		temppercentfull.append(round((summm/threshold)*100.,3))#an accumulation valriable for each id
		temp_tons_still_needed.append(round(new_threhold-((t_contributed/new_threhold)*new_threhold),3))
		tempexcess2.append(t2)#original value
		weightss.append(sps_vals[dist2[i]])#to use for the location of the end point
		points1.append(sps[dist2[i]])
		if t2>=new_threhold:
			break
		new_threhold=new_threhold-t_contributed
		
	startpointsnew.append(points1)	
	percentfull.append(temppercentfull)
	tons_still_needed.append(temp_tons_still_needed)
	tons_contributed.append(temp_tons_contributed)	
	excess.append(tempexcess)
	excess2.append(tempexcess2)
	perecent.append(tempperecent)
	aa=tuple(np.average(points1, axis=0,weights=weightss))
	ep_ids.append(accumid)
	accumid=accumid+1
	sps_ids_etmp.append(ids1)
	ep.append(aa)
	sp_ids2=sps_ids_etmp.copy()
	excess4=excess.copy()
	for i in range(0,len(sp_ids2)):
		for j in range(0,len(sp_ids2[i])):
			index=sps_ids.index(sp_ids2[i][j])
			if excess[i][j] <= 0:
				excess4[i][j]=0
				sps_vals.remove(sps_vals[index])
				sps.remove(sps[index])
				sps_ids.remove(sps_ids[index])
			else:
				sps_vals[index]=excess[i][j]

	Dict['End_point_id']=ep_ids
	Dict['end_points']=ep	
	Dict['start_point_id']=sps_ids_etmp
	Dict['start_points']=startpointsnew
	Dict['original values']=excess2
	Dict['excess']=excess4
	Dict['tonnes contributed']=tons_contributed
	Dict['tonnes still needed']=tons_still_needed
	Dict['percent facility full']=percentfull
	Dict['percent feedstock needed from current start_point']=perecent
	Dict['id']=sps_ids
	Dict['values']=sps_vals
	Dict['points']=sps
	return(Dict)

def Check_VR_or_NGI_roads_proximity(raster,raster2,points,ids,showtime):
	'''
	@ Author Marc Pienaar
	Checks if the points are closer to an exisitng road or virtual road using pre computed proximity rasters (@ 100m)
	Parameters
	----------
	raster:  String
		Filename of the VR proximity raster to read.

	raster2:  String
		Filename of the NGI or exisiting roads proximity raster to read.
		
	points : array[]
		list of coordintates to chck for proximtiy

	ids : array[]
		the point array ids

	showtime :  bool
		If True, prints the time this function took to run.
		If False, prints nothing.
	
	Returns
	-------
	Dictionary object:  values, percentages and assocaited with either VR or NGI road proximity		
	'''
	if showtime:
		start_time = time.time()
	ds = gdal.Open(raster)
	ds2 = gdal.Open(raster2)
	target = osr.SpatialReference(wkt=ds.GetProjection())
	source = osr.SpatialReference()
	source.ImportFromEPSG(4326)#wgs84
	transformproj = osr.CoordinateTransformation(source,target)
	
	transform = ds.GetGeoTransform()
	xOrigin = transform[0] 
	yOrigin = transform[3]
	pixelWidth = transform[1] 
	pixelHeight = transform[5] 
	band = ds.GetRasterBand(1) 
	data = band.ReadAsArray()
	#get the bands from each raster	
	band2 = ds2.GetRasterBand(1) # 1-based index
	data2 = band2.ReadAsArray()	
	# loop through the coordinates
	values=[]
	values2=[]
	Dict = {}
	Dict['Points']=points
	vr_points=[]
	NGI_points=[]
	vr_ids=[]
	NGI_ids=[]
	i=0
	for point in points:
		newpoints=[]
		x2 = point[1]
		y2 = point[0]
		newpoints.append(y2)
		newpoints.append(x2)
		aa=tuple(newpoints)
		point = ogr.Geometry(ogr.wkbPoint)
		point.AddPoint(x2,y2)
		point.Transform(transformproj)
		x=point.GetX()
		y=point.GetY()
		xOffset = int((x - xOrigin) / pixelWidth)
		yOffset = int((y - yOrigin) / pixelHeight)		
		# get individual pixel values from roads and VR proximity layers
		value = data2[yOffset][xOffset]#roads
		value2 = data[yOffset][xOffset]#VR roads
		if value <= value2:#do normal roads point
			NGI_points.append(aa)
			NGI_ids.append(ids[i])
		else:
			vr_points.append(aa)
			vr_ids.append(ids[i])
		i=i+1
	Dict['VR points']=vr_points
	Dict['VR ids']=vr_ids
	Dict['NGI points']=NGI_points
	Dict['NGI ids']=NGI_ids
	#close resources
	ds=None
	ds2=None
	if showtime:
		print("Proximity checking took --- %s seconds ---" % (time.time() - start_time))
		print()#print a new line
	return(Dict)

def get_vr_road(drivername,centoids_in,centoids_in2,points,point_ids,showtime):
	'''
	@ Author Marc Pienaar
	Maps the VR start point to a VR end point (i.e along a road), and get the distance travelled
	----------
	drivername:  String
		Ogr driver name to handle the file e.g. ESRI_Shapefile or GPKG

	centoids_in:  String
		Filename of the VR start point shapefile.

	centoids_in2:  String
		Filename of the VR end point shapefile.
		
	points : array[]
		list of coordintates to check for VR mapping

	point_ids : array[]
		the point array ids

	showtime :  bool
		If True, prints the time this function took to run.
		If False, prints nothing.
	
	Returns
	-------
	Dictionary object:  values, percentages and assocaited with either VR or NGI road proximity		
	'''
	if showtime:
		start_time = time.time()
	shpDriver = ogr.GetDriverByName(drivername)
	source_ds = ogr.Open(centoids_in)  # open the original shp
	layer = source_ds.GetLayer()
	featureCount = layer.GetFeatureCount()
	ldefn = layer.GetLayerDefn()
	#get source points
	x1=[];y1=[];ids=[];ids_new=[];lengths=[]
	#get the points and ids of the start
	
	for feature in layer:#this takes long
		ids.append(feature.GetField("fid2"))
		lengths.append(feature.GetField("length"))
		point = feature.GetGeometryRef()
		x1.append(point.GetX())
		y1.append(point.GetY())
	source_ds=None
	X=np.column_stack((x1,y1)) 
	xx=np.asarray(points)
	dist=euclidean_distances(xx,X)
	fid=[]
	s_points=[]
	Dict={}
	for i in range(len(dist)):
		start_points=[]
		dist2 = np.argsort(dist[i].flatten())
		check=[]
		for j in range(4):
			check.append(round(dist[i][dist2[j]],6))
		idx = np.where(np.asarray(check)==np.asarray(check).min())
		if len(idx[0]) == 1:
			fid.append(ids[dist2[idx[0][0]]])
			start_points.append(x1[dist2[idx[0][0]]])
			start_points.append(y1[dist2[idx[0][0]]])
		else:
			check2=[]
			for j in idx[0]:
				check2.append(lengths[dist2[idx[0][j]]])
			idx2 = np.where(np.asarray(check2)==np.asarray(check2).min())
			fid.append(ids[dist2[idx[0][idx2[0][0]]]])
			start_points.append(x1[dist2[idx[0][idx2[0][0]]]])
			start_points.append(y1[dist2[idx[0][idx2[0][0]]]])
		aa=tuple(start_points)
		s_points.append(aa)
	Dict['start points']=s_points
	Dict['ids']=point_ids
	shpDriver = ogr.GetDriverByName(drivername)
	source_ds = ogr.Open(centoids_in2)  # open the original shp
	layer = source_ds.GetLayer()
	#featureCount = layer.GetFeatureCount()
	accum=len(fid)
	e_points=[]
	e_lengths=[]
	for i in fid:
		end_points=[]
		for feature in layer:
			if feature.GetField("fid2")==i:
				accum=accum-1
				point = feature.GetGeometryRef()
				end_points.append(point.GetX())
				end_points.append(point.GetY())
				aa=tuple(end_points)
				e_points.append(aa)
				e_lengths.append(feature.GetField("length"))
				break# add a break function to speed up
		if accum==0:
			break	
	source_ds=None
	Dict['end points']=e_points
	Dict['length (m)']=e_lengths
	if showtime:
		print("VR mapping took --- %s seconds ---" % (time.time() - start_time))
		print()#print a new line
	return Dict

def adjust_points(drivername, proximity,VR_vertices_in_wgs84,VR_vertices_end_wgs84,showtime):
	'''
	@ Author Marc Pienaar
	Maps any VR start points to a VR end point (i.e along a road), and get the distance travelled and the new poitn coordinates
	----------
	drivername:  String
		Ogr driver name to handle the file e.g. ESRI_Shapefile or GPKG

	proximity:  Dictionary object of VR values from Check_VR_or_NGI_roads_proximity function		

	VR_vertices_in_wgs84:  String
		Filename of the VR start point shapefile.

	VR_vertices_end_wgs84:  String
		Filename of the VR end point shapefile.

	showtime :  bool
		If True, prints the time this function took to run.
		If False, prints nothing.
	
	Returns
	-------
	Three arrays: new start points, new start point ids, and any Virtual road distances travelled
	'''
	if showtime:
		start_time = time.time()
	start_points_temp=[]
	ids_temp=[]
	vr_distance=[]
	#first do a proximity check
	if len(proximity['NGI points'])>0:
		for i in proximity['NGI points']:
			start_points_temp.append(i)
			vr_distance.append(0)
		for i in proximity['NGI ids']:
			ids_temp.append(i)
	#no do a VR mapping
	if len(proximity['VR points'])>0:
		vr_out=get_vr_road(drivername,VR_vertices_in_wgs84,VR_vertices_end_wgs84,proximity['VR points'],proximity['VR ids'],False)
		for i in vr_out['end points']:
			start_points_temp.append(i)
		for i in proximity['VR ids']:
			ids_temp.append(i)
		for i in vr_out['length (m)']:
			vr_distance.append(i)
	if showtime:
		print("Adjusting VR points took --- %s seconds ---" % (time.time() - start_time))
		print()#print a new line			
	return start_points_temp,ids_temp,vr_distance

def get_combined_points(start_points,new_start_ids,all_start_ids,end_points):	
	'''
	@ Author Marc Pienaar
	combine points
	Parameters
	----------
	start_points:  point tuple array
	start_ids:   ids of start point tuple array to extract
	all_start_ids: 
	end_points:  point tuple array
	Returns
	-------
	a Dictionary of combined points		
	'''
	start_time = time.time()
	sp_temp=start_points
	ep_temp=[]
	ep_temp.append(end_points)	
	points2=[]
	for i in sp_temp:
		points2.append(i)
	for i in ep_temp:
		points2.append(i)
	Dict={}	
	Dict['start_points']=sp_temp
	Dict['end_points']=ep_temp
	Dict['combined_points']=points2
	return(Dict)

def create_extent(drivername,points,shapein,outShapefile,offset):
	'''
	@ Author Marc Pienaar
	Creates a shapefile reprsenting an extent from a bunch of points
	----------
	drivername:  String
		Ogr driver name to handle the file e.g. ESRI_Shapefile or GPKG

	points : array[]
		list of coordintates 		

	shapein:  String
		Filename of the the input shapefile.

	outShapefile:  String
		Filename of the temp polygon extent (a rectangle polygon) to write to

	offset :  Numeric Real
		a value in degrees to buffer the extent by 0.2 in WGS84 is roughly 20km
	
	Returns
	-------
	A temporay polygion shapefile representing and extent of the input points + the offset buffer
	'''
	start_time = time.time()
	Dict={}
	NEW_Points=[]
	minx=sys.float_info.max
	maxx=sys.float_info.min
	miny=sys.float_info.max
	maxy=-1000
	for i in points:
		minx=min(minx,i[0])
		maxx=max(maxx,i[0])
		miny=min(miny,i[1])
		maxy=max(maxy,i[1])
	
	inDriver = ogr.GetDriverByName(drivername)
	inDataSource = inDriver.Open(shapein, 0)
	inLayer = inDataSource.GetLayer()
	srs = inLayer.GetSpatialRef()
	ring = ogr.Geometry(ogr.wkbLinearRing)
	ring.AddPoint(minx-offset,miny-offset)
	ring.AddPoint(maxx+offset,miny-offset)
	ring.AddPoint(maxx+offset, maxy+offset)
	ring.AddPoint(minx-offset, maxy+offset)
	ring.AddPoint(minx-offset,miny-offset)
	poly = ogr.Geometry(ogr.wkbPolygon)
	poly.AddGeometry(ring)
	# Save extent to a new Shapefile
	outDriver = ogr.GetDriverByName(drivername)
	# Create the output shapefile
	outDataSource = outDriver.CreateDataSource(outShapefile)
	outLayer = outDataSource.CreateLayer("extent", srs,geom_type=ogr.wkbPolygon)
	# Add an ID field
	idField = ogr.FieldDefn("id", ogr.OFTInteger)
	outLayer.CreateField(idField)
	# Create the feature and set values
	featureDefn = outLayer.GetLayerDefn()
	feature = ogr.Feature(featureDefn)
	feature.SetGeometry(poly)
	feature.SetField("id", 1)
	outLayer.CreateFeature(feature)
	# Close DataSource
	inDataSource.Destroy()
	outDataSource.Destroy()
	
def clip_vector_by_shape(driverName,driverName2,driverName3,shapein,clipsrc,outshape,typeout):
	'''
	@ Author Marc Pienaar
	Clips a vector by another vector
	----------
	drivername:  String
		Ogr driver name to handle the input file e.g. ESRI_Shapefile or GPKG

	Drivername2:  String
			Ogr driver name to handle the clipsrc file e.g. ESRI_Shapefile or GPKG

	Drivername3:  String
			Ogr driver name to handle the output file e.g. ESRI_Shapefile or GPKG
	
	shapein:  String
		Filename of the the input shapefile.

	clipsrc:  String
		Filename of the the temporat extent of clipping boundary shapefile.

	outshape:  String
		Filename of the temp clipped shapefile to write to

	typeout :  String
		'line' or 'point' to deternine whether to clip a input linestring or point shapefile
	
	Returns
	-------
	A clipped version of the input file
	'''
	
	driver = ogr.GetDriverByName(driverName)
	driver2=ogr.GetDriverByName(driverName2)
	driver3=ogr.GetDriverByName(driverName3)
	inDataSource = driver.Open(shapein, 0)
	inLayer = inDataSource.GetLayer()
	srs = inLayer.GetSpatialRef()
	## Clip
	inClipSource = driver2.Open(clipsrc, 0)
	inClipLayer = inClipSource.GetLayer()
	outDataSource = driver3.CreateDataSource(outshape)
	if typeout=="line":
		outLayer = outDataSource.CreateLayer('clip', srs,geom_type=ogr.wkbMultiLineString)
	if typeout=="point":
		outLayer = outDataSource.CreateLayer('clip', srs,geom_type=ogr.wkbPoint)
		
	ogr.Layer.Intersection(inLayer, inClipLayer, outLayer,["SKIP_FAILURES=YES","PROMOTE_TO_MULTI=YES"])
	#write and release resources
	inDataSource.Destroy()
	inClipSource.Destroy()
	outDataSource.Destroy()
	
def create_new_point_layer2(drivername, shapein, shpout,points,attrname,attrval,attr2,attr3):
	'''
	@Author Marc Pienaar
	Create a new point shapefile
	Parameters
	----------
	drivername:  String
		Ogr driver name to handle the file e.g. ESRI_Shapefile or GPKG

	shapein : File or string
		File, directory, or filename to read.

	shpout : File or string
		File, directory, or filename to write to.
		
	points : list of points (tuple)
		x,y, coordinates

	attrname : String
		Filed name to create

	attrval : Numeric 
		the values to etie to the attribute
	
	Returns
	-------
	a shapefile 		
	'''
	start_time = time.time()
	# set up the shapefile driver
	driver = ogr.GetDriverByName(drivername)
	inDataSource = driver.Open(shapein, 0)
	inLayer = inDataSource.GetLayer()
	srs = inLayer.GetSpatialRef()
	inDataSource.Destroy()
	driver = ogr.GetDriverByName(drivername)
	## create the data source
	data_source = driver.CreateDataSource(shpout)
	## create the layer
	layer = data_source.CreateLayer("points", srs, ogr.wkbPoint)
	layer.CreateField(ogr.FieldDefn(attrname, ogr.OFTInteger64))
	layer.CreateField(ogr.FieldDefn("technology", ogr.OFTString))
	layer.CreateField(ogr.FieldDefn("capacity", ogr.OFTInteger64))
#	layer.CreateField(ogr.FieldDefn("length", ogr.OFTReal))
	for i in range(0,len(points)):
		feature = ogr.Feature(layer.GetLayerDefn())
		feature.SetField(attrname, attrval[i])
		feature.SetField("technology", attr2)
		feature.SetField("capacity", attr3)
		point1 = Point(points[i][0],points[i][1])
		# Create the point from the Well Known Txt
		point = ogr.CreateGeometryFromWkt(point1.wkt)
		# Set the feature geometry using the point
		feature.SetGeometry(point)
		# Create the feature in the layer (shapefile)
		layer.CreateFeature(feature)
		# Dereference the feature
		feature = None
	# Save and close the data source
	data_source = None
	
def get_road_aligned_points(roads_temp,combined_points):
	'''
	@Author Marc Pienaar
	snaps point to an exiting linestring using shapely.ops
	Parameters
	----------
	roads_temp:  String
		he name of the LineString shapefile

	combined_points : array []
		Array of start and end points (the last coordinate is the end point)
	
	Returns
	-------
	Adjusted start, end and combined point arrays snapped to the input linestring file		
	'''
	gdf_segments = gpd.read_file(roads_temp)
	shply_line = gdf_segments.geometry.unary_union
	#declare some empty variables
	start_points_adjusted=[];end_point_adjusted=[];sp=[];ep=[];end_point_adjusted=[];points=[]
	accum=0
	length=len(combined_points['combined_points'])
	for i in combined_points['combined_points']:
		point = Point(i)
		new_point = nearest_points(shply_line, point)[0]
		if accum<length-1:
			start_points_adjusted.append(new_point.coords[0])
			sp.append(i)
		if accum==length-1:
			end_point_adjusted.append(new_point.coords[0])
			ep.append(i)
		accum=accum+1
		points.append(new_point)
	multi_point = MultiPoint(points)
	start=start_points_adjusted
	end=end_point_adjusted
	return start,end,multi_point

def create_output_shapefiles(outputdir,facilities,capacity1):
	'''
	@Author Marc Pienaar
	creates shapefile names
	Parameters
	----------
	outputdir:  String
		The name of the current output directory

	facilities : String
		Name of the current technology

	capacity1 : Numeric
		value of the current capacity of the facility
	
	Returns
	-------
	shapefile paths to create		
	'''
	outpaths = os.path.sep.join([outputdir, "facility" + str(facilities['End_point_id'][0])+"_paths_"+str(capacity1)+".gpkg"])	
	outpaths_sp_original = os.path.sep.join([outputdir, "facility" + str(facilities['End_point_id'][0])+"_start_points_original_capacity_"+str(capacity1)+".gpkg"])
	
	outpaths_sp = os.path.sep.join([outputdir, "facility" + str(facilities['End_point_id'][0])+"_start_points_capacity_"+str(capacity1)+".gpkg"])	
	
	outpaths_ep = os.path.sep.join([outputdir, "facility" + str(facilities['End_point_id'][0])+"_end_point_capacity_"+str(capacity1)+".gpkg"])
	
	outpaths_area =  os.path.sep.join([outputdir,"facility" + str(facilities['End_point_id'][0])+"_catchment_capacity_"+str(capacity1)+".gpkg"])	
	
	return outpaths,outpaths_sp_original,outpaths_sp,outpaths_ep,outpaths_area

def geom_type(shapein):
	'''
	Gete the geometry type of a shapefile
	Parameters
	----------
	shapein : File or string
		File, directory, or filename to read.
	
	Returns
	-------
	int as a gdal geometry type		
	'''
	
	shapefile = ogr.Open(shapein)
	layer = shapefile.GetLayer()
	layer_defn = layer.GetLayerDefn()
	out=layer_defn.GetGeomType()
	shapefile=None
	return(out)


def create_field_scenario(drivername,drivername2,facilities,attrs1,attrs1_ids,field_yield,outhsape,attr,attr2,attr3):
	'''
	@Author Marc Pienaar
	Gets the catchment area associated with a facility
	Parameters
	----------
	drivername:  String
		Ogr driver name to handle the file e.g. ESRI_Shapefile or GPKG

	drivername2:  String
		Ogr driver name to handle the file e.g. ESRI_Shapefile or GPKG

	facilities : Array []
		Facility IDs

	attrs1 : Dictionary
		Attribute names of the input shapefile
		
	attrs1_ids : Array []
		indecices of the Attributes to read

	field_yield : String
		File name of the polygon shapefil to read
	
	outhsape : String
		File name of the shapefile to write

	attr : String 
		The name of the attribute from the field_yield to read

	attr2 : String 
		The name of Facility technology

	attr3 : String 
		The capacity of the Facility technology
	
	Returns
	-------
	Writes a shapefile out and returns the area value associated withe the qttibue field		
	'''
	
	start_time = time.time()
	select=[]
	geomtype=geom_type(field_yield)	
	shpDriver = ogr.GetDriverByName(drivername)
	shpDriver2 = ogr.GetDriverByName(drivername2)
	source_ds = ogr.Open(field_yield)  # open the resource
	inLayer = source_ds.GetLayer()
	srs = inLayer.GetSpatialRef()
	outDataSet = shpDriver2.CreateDataSource(outhsape)
	outLayer = outDataSet.CreateLayer("catchment",srs,geom_type=geomtype)
	featureCount = inLayer.GetFeatureCount()
	inLayerDefn = inLayer.GetLayerDefn()
	outLayer.CreateField(ogr.FieldDefn("technology", ogr.OFTString))
	outLayer.CreateField(ogr.FieldDefn("capacity", ogr.OFTInteger64))
	for i in range(0, inLayerDefn.GetFieldCount()):
		fieldDefn = inLayerDefn.GetFieldDefn(i)
		outLayer.CreateField(fieldDefn)
	# get the output layer's feature definition
	outLayerDefn = outLayer.GetLayerDefn()
	temp=[]
	for ii in facilities:
		temp2=[]
		for feature in inLayer:
			if feature.GetField(attrs1['attributes'][attrs1_ids[0]][1])==ii:
				temp2.append(feature.GetField(attr))
				geom = feature.GetGeometryRef()
				outFeature = ogr.Feature(outLayerDefn)
				outFeature.SetGeometry(geom)
				outFeature.SetField("technology", attr2)
				outFeature.SetField("capacity", attr3)
				for i in range(0, inLayerDefn.GetFieldCount()):
					outFeature.SetField(inLayerDefn.GetFieldDefn(i).GetNameRef(), feature.GetField(i))
				outLayer.CreateFeature(outFeature)
				# dereference the features and get the next input feature
				outFeature = None
		values=np.asarray(temp2)
		temp.append(np.sum(values))
	inDataSet = None
	outDataSet=None
	return(temp)#return the area value


def add_node2(lines, point, tolerance=1e-7):
	'''
	library with functions to edit and perform spatial analysis on vector data.
	build around shapely and fiona libraries
	
	Created on 2016-July-16
	@author: Dirk Eilander (dirk.eilander@deltares.nl)
	
	# build on library from https://github.com/ojdo/python-tools/blob/master/shapelytools.py

	Split line at a given point

	Args:
		lines: a list of shapely LineStrings or a MultiLineString
		point: shapely Point
		tolerance: required to check if segment intersects with line

	Returns:
		a list of LineStrings
	'''
#		# find line which intersects with point, but not with start of end point
#		idx = [i for i, line in enumerate(lines) if
#				(line.distance(point) <= tolerance) & (line.boundary.distance(point) > tolerance)]
#		
#		if len(idx) == 0:
#			raise Warning('line does not intersect with point with given tolerance.')
#		elif len(idx) > 1:
#			raise Warning('more than one line within given tolerance')
#		else:
#			idx = idx[0]
	
	# for intersecting line, find intersecting segment and split line
	coords = list(lines.coords)
	segments = [LineString(s) for s in pairs(coords)]
	n = len(segments)
	val=-1
	
	for i, segment in enumerate(segments):
		# find intersecting segment
		if segment.distance(point) <= tolerance:
			
			# break if at point at existing node (do nothing)
			if Point(coords[i]).distance(point) <= tolerance:
				break
			# otherwise et node
			else:
				val=i+1
				lines = LineString(coords[:i+1] + [(point.x, point.y)] + coords[i+1:])
				break
	return lines

def updatenode(linesin,pointsin,outshp):
	#from: https://medium.com/@brendan_ward/how-to-leverage-geopandas-for-faster-snapping-of-points-to-lines-6113c94e59aa
	lines = gpd.read_file(linesin)
	points=gpd.read_file(pointsin)
	lines.to_crs(epsg=32662,inplace=True)
	points.to_crs(epsg=32662,inplace=True)
	offset = 1#0.000005
	tolerance=offset
	bbox = points.bounds + [-offset, -offset, offset, offset]
	hits = bbox.apply(lambda row: list(lines.sindex.intersection(row)), axis=1)
	tmp = pd.DataFrame({
		# index of points table
		"pt_idx": np.repeat(hits.index, hits.apply(len)),
		# ordinal position of line - access via iloc later
		"line_i": np.concatenate(hits.values)
	})
	# Join back to the lines on line_i; we use reset_index() to 
	# give us the ordinal position of each line
	tmp = tmp.join(lines.reset_index(drop=True), on="line_i")
	# Join back to the original points to get their geometry
	# rename the point geometry as "point"
	tmp = tmp.join(points.geometry.rename("point"), on="pt_idx")
	# Convert back to a GeoDataFrame, so we can do spatial ops
	tmp = gpd.GeoDataFrame(tmp, geometry="geometry", crs=points.crs)
	tmp["snap_dist"] = tmp.geometry.distance(gpd.GeoSeries(tmp.point))
	#print(tmp)
	## Discard any lines that are greater than tolerance from points
	tmp = tmp.loc[tmp.snap_dist <= tolerance]
	### Sort on ascending snap distance, so that closest goes to top
	tmp = tmp.sort_values(by=["snap_dist"])
	### group by the index of the points and take the first, which is the
	### closest line 
	closest = tmp.groupby("pt_idx").first()
	### construct a GeoDataFrame of the closest lines
	closest = gpd.GeoDataFrame(closest, geometry="geometry")
	
	
#	start_time = time.time()
	#filter
	while len(closest) >0:
		index_1=[]
		for idx,row in closest.iterrows():
			index_1.append(idx)
		temp2=closest[closest['OBJECTID']==closest.loc[index_1[0]]['OBJECTID']]
		index=[]
		for idx,row in temp2.iterrows():
			index.append(idx)
		for i in index:
			linest=add_node2(closest['geometry'].geometry[index[0]],temp2['point'].geometry[i])
			closest.loc[index[0], 'geometry']=linest
		#now update the dataframe and remoce the indices from closest and snapped
		temp_0=lines[lines['OBJECTID']==closest.loc[index_1[0]]['OBJECTID']]
		index_0=[]
		for idx,row in temp_0.iterrows():
			index_0.append(idx)
		lines.loc[index_0[0], 'geometry']=closest.loc[index[0], 'geometry']
		closest=closest[closest['OBJECTID']!=closest.loc[index_1[0]]['OBJECTID']]
#		print(len(closest))
		
	lines.to_crs(epsg=4326,inplace=True)
	lines.to_file(driver = 'ESRI Shapefile', filename= outshp)
	

	
def remove_redundant_nodes(lines, tolerance=1e-7):
	'''remove vertices with length smaller than tolerance'''
	lines_out = []
	for line in lines:
		coords = line.coords
		l_segments = np.array([Point(s[0]).distance(Point(s[1])) for s in pairs(coords)])
		idx = np.where(l_segments < tolerance)[0]
		lines_out.append(LineString([c for i, c in enumerate(coords) if i not in idx]))
	lines2 = MultiLineString(lines_out)
	return lines2

def remove_redundant_nodes2(lines, tolerance=1e-7):
	'''remove vertices with length smaller than tolerance'''
	lines_out = []
#	for line in lines:
	coords = lines.coords
	l_segments = np.array([Point(s[0]).distance(Point(s[1])) for s in pairs(coords)])
	idx = np.where(l_segments < tolerance)[0]
	return(LineString([c for i, c in enumerate(coords) if i not in idx]))
#	lines_out

def add_node(lines, point, tolerance=1e-7):
	'''
	library with functions to edit and perform spatial analysis on vector data.
	build around shapely and fiona libraries
	
	Created on 2016-July-16
	@author: Dirk Eilander (dirk.eilander@deltares.nl)
	
	# build on library from https://github.com/ojdo/python-tools/blob/master/shapelytools.py

	Split line at a given point

	Args:
		lines: a list of shapely LineStrings or a MultiLineString
		point: shapely Point
		tolerance: required to check if segment intersects with line

	Returns:
		a list of LineStrings
	'''
#		# find line which intersects with point, but not with start of end point
#		idx = [i for i, line in enumerate(lines) if
#				(line.distance(point) <= tolerance) & (line.boundary.distance(point) > tolerance)]
#		
#		if len(idx) == 0:
#			raise Warning('line does not intersect with point with given tolerance.')
#		elif len(idx) > 1:
#			raise Warning('more than one line within given tolerance')
#		else:
#			idx = idx[0]
	
	# for intersecting line, find intersecting segment and split line
	coords = list(lines.coords)
	segments = [LineString(s) for s in pairs(coords)]
	n = len(segments)
	val=-1
	for i, segment in enumerate(segments):
		# find intersecting segment
		if segment.distance(point) <= tolerance:
			
			# break if at point at existing node (do nothing)
			if Point(coords[i]).distance(point) <= tolerance:
				break
			# otherwise et node
			else:
				val=i+1
				lines = LineString(coords[:i+1] + [(point.x, point.y)] + coords[i+1:])
				break
	return val
	
def pairs(lst):
	'''
	library with functions to edit and perform spatial analysis on vector data.
	build around shapely and fiona libraries
	
	Created on 2016-July-16
	@author: Dirk Eilander (dirk.eilander@deltares.nl)
	
	# build on library from https://github.com/ojdo/python-tools/blob/master/shapelytools.py
	
	Iterate over a list in overlapping pairs.

	Args:
		lst: an iterable/list

	Returns:
		Yields a pair of consecutive elements (lst[k], lst[k+1]) of lst. Last
		call yields (lst[-2], lst[-1]).

	Example:
		lst = [4, 7, 11, 2]
		pairs(lst) yields (4, 7), (7, 11), (11, 2)

	Source:
		http://stackoverflow.com/questions/1257413/1257446#1257446
	'''
	i = iter(lst)
	prev = next(i)
	for item in i:
		yield prev, item
		prev = item

	
def routing_check(shapefile,start,end):
	'''
	@ Author Marc Pienaar 
	Generates a networkx.DiGraph to ceheck if the road network in navigable 

	Parameters
	----------
	shapefile : file or string
		File, directory, or filename to read.
		
	start and end : tuple
		e.g. [(25.5789695,-32.7209891)], or  [(25.578967,-32.721001),(25.578658,-32.718445)]
		Can't have multiple start and end arrays	
		
	shpOut : file or string
		File, directory, or filename to write to.

	usefieldweight :  bool
			If True, uses a fieldname from the input shapefile for additional weighting, if that field doesn't exist it defaults to 1.
			If False, no additional weighting.
	
	fieldname : String
		The fieldname in a shapefile with numeric values to use for additional edge weighting e.g. "weight"
	
	useexclusions : bool
		If True, uses a fieldname from the input shapefile for exluding edges / nodes
		If False, exclusions.
	
	exclusion_field : String 
		The fieldname in a shapefile with values to use for exclusions
		
	exclusion_vals : [String]
		attribute values in the field name e.g ["bad road", "dirt road"]
			
	ave_travel_time : float
		The average travel time is in Km/h e.g. 100.0.
		this value is used to calcualte average travel time in the function as distance / travel time. 
		if you assume an average travel time for a good road, you can use the fieldname variable to weight "worse"
		roads. i.e if your average travel time is 120 km/h assumed for a tar road and you weight all segments of 
		sand raods as 2 (the others 1), then it will double the avarege travel for all segaments that are weighted 2 
	
	returnshape : bool
			If True, writes out a shapefile specified by shpOut as a linestring, with length in meters as an attribute as well as the distandce in meters as a output
			If False, Does not write out a shapefile, but does still return the distance in meters of the path 	

	s_point_ids : Array []
		Start point IDs

	e_point_id : Array []
		End point IDs

	tons : Numeric
		The number of tons biomass being transported

	service_area : Numeric
		The area being serviced 

	attr2 : String
		he name of the facility technology

	attr3 : String
		the capacity of the facility technology


	Returns
	-------
	
	A bool of true or false 
	'''	
	returnval=False
	start_time = time.time()
	#create the graph from a shapefile
	G =nx.read_shp(shapefile,simplify=False)
	G=G.to_undirected()
	start_time = time.time()
	pos = {k: v for k,v in enumerate(G.nodes())}
	edg=[]
	for (u, v,c) in G.edges.data():
		a=[]
		a.append(u)
		a.append(v)
		aa=tuple(a)
		edg.append(aa)
	#print("--- %s total seconds greating new edge array---" % (time.time() - start_time))	
	start_time = time.time()
	temp=[]#
	for i in pos.keys():
		temp.append(pos[i])
	start2=[]
	end2=[]
	startcoords=[]
	endcoords=[]
	for i in start:
		dist=euclidean_distances(np.asarray(i).reshape(1,-1),np.asarray(temp))
		dist2 = np.argsort(dist.flatten())
		startcoords.append(temp[dist2[0]])
		start2.append(temp[dist2[0]])
	for i in end:
		dist=euclidean_distances(np.asarray(i).reshape(1,-1),np.asarray(temp))
		dist2 = np.argsort(dist.flatten())
		endcoords.append(temp[dist2[0]])
		end2.append(temp[dist2[0]])	
	start_time = time.time()
	dist=0
	length, path = nx.single_source_dijkstra(G, end2[0])
	try:
		for i in range(0,len(start2)):
			pathGraph = nx.path_graph(path[start2[i]]) 
	except:
		returnval =True
	return returnval
		
def routing(shapefile,start,end,shpOut,usefieldweight,fieldname,useexclusions,exclusion_field,exclusion_vals,ave_travel_time,returnshape,s_point_ids,e_point_id,tons,service_area,attr2,attr3,vr_distance):
	'''
	@ Author Marc Pienaar 
	Generates a networkx.DiGraph from shapefiles and finds the shortest path between two points using Dijkstra's algorithm. Point geometries are  translated into nodes, lines into edges. Coordinate tuples are used as keys. Attributes are preserved. The deafult is projection is wgs84. The distance in meters between nodes is automatically calcualted in this function and used as the default edge weighting, additional weightings can be applied through a field name in the shapefile with weighting values and / or a constant weight, incase of batch processing. All additional weightings are multiplied by the default weighting in the graph"

	Parameters
	----------
	shapefile : file or string
		File, directory, or filename to read.
		
	start and end : tuple
		e.g. [(25.5789695,-32.7209891)], or  [(25.578967,-32.721001),(25.578658,-32.718445)]
		Can't have multiple start and end arrays	
		
	shpOut : file or string
		File, directory, or filename to write to.

	usefieldweight :  bool
			If True, uses a fieldname from the input shapefile for additional weighting, if that field doesn't exist it defaults to 1.
			If False, no additional weighting.
	
	fieldname : String
		The fieldname in a shapefile with numeric values to use for additional edge weighting e.g. "weight"
	
	useexclusions : bool
		If True, uses a fieldname from the input shapefile for exluding edges / nodes
		If False, exclusions.
	
	exclusion_field : String 
		The fieldname in a shapefile with values to use for exclusions
		
	exclusion_vals : [String]
		attribute values in the field name e.g ["bad road", "dirt road"]
			
	ave_travel_time : float
		The average travel time is in Km/h e.g. 100.0.
		this value is used to calcualte average travel time in the function as distance / travel time. 
		if you assume an average travel time for a good road, you can use the fieldname variable to weight "worse"
		roads. i.e if your average travel time is 120 km/h assumed for a tar road and you weight all segments of 
		sand raods as 2 (the others 1), then it will double the avarege travel for all segaments that are weighted 2 
	
	returnshape : bool
			If True, writes out a shapefile specified by shpOut as a linestring, with length in meters as an attribute as well as the distandce in meters as a output
			If False, Does not write out a shapefile, but does still return the distance in meters of the path 	

	s_point_ids : Array []
		Start point IDs

	e_point_id : Array []
		End point IDs

	tons : Numeric
		The number of tons biomass being transported

	service_area : Numeric
		The area being serviced 

	attr2 : String
		he name of the facility technology

	attr3 : String
		the capacity of the facility technology

	vr_distance: Real or Double
		The distance travelled on a virtual road


	Returns
	-------
	
	A dictionary object of with various information 
	If return shape is True - It will also save a linestring shapefile with length of the path(s) in meters and 
	travel time in hours as attribute values
	'''	
	start_time = time.time()
	#create the graph from a shapefile
	G =nx.read_shp(shapefile,simplify=False)
	G=G.to_undirected()
	start_time = time.time()
	pos = {k: v for k,v in enumerate(G.nodes())}
	edg=[]
	for (u, v,c) in G.edges.data():
		a=[]
		a.append(u)
		a.append(v)
		aa=tuple(a)
		edg.append(aa)
	#print("--- %s total seconds greating new edge array---" % (time.time() - start_time))	
	start_time = time.time()
	temp=[]#
	for i in pos.keys():
		temp.append(pos[i])
	start2=[]
	end2=[]
	startcoords=[]
	endcoords=[]
	for i in start:
		dist=euclidean_distances(np.asarray(i).reshape(1,-1),np.asarray(temp))
		dist2 = np.argsort(dist.flatten())
		startcoords.append(temp[dist2[0]])
		start2.append(temp[dist2[0]])
	for i in end:
		dist=euclidean_distances(np.asarray(i).reshape(1,-1),np.asarray(temp))
		dist2 = np.argsort(dist.flatten())
		endcoords.append(temp[dist2[0]])
		end2.append(temp[dist2[0]])	
		
	start_time = time.time()
	dist=0
	distance_weights=[]	
	#if we use weights
	if usefieldweight:
		for (u, v, c) in G.edges.data():
			try:
				dist=distance(lonlat(*u), lonlat(*(v))).meters
				G[u][v]['weight']=dist*float(c[fieldname])
			except:
				G[u][v]['weight']=dist
	else:
		for (u, v, c) in G.edges.data():
			dist=distance(lonlat(*u), lonlat(*(v))).meters
			G[u][v]['weight']=dist
	start_time = time.time()
	#if we use exclusions
	if useexclusions:
		exclusions=[]
		for (u, v, c) in G.edges.data():
			temp=0
			for i in exclusion_vals:
				if c[exclusion_field]==i:
					temp=temp+1
				else:
					pass
			if temp>0:
				exclusions.append(1)
			else:
				exclusions.append(0)
		ebunch=[]
		for j in range(0,len(exclusions)):
			if exclusions[j]==1:
				ebunch.append(edg[j])
		G.remove_edges_from(ebunch)

	Dict = {}
	distout=[]
	pathouot=[]
	timeout=[]
	if len(end2)>1 and len(start2)>1:
		print("can't have multiple source and destination points")
		return
	if len(end2)==1 and len(start2)==1:#one start and end point
		#calculate the optimal path base on the weights in the graph
		start_time = time.time()
		length, path = nx.single_source_dijkstra(G, start2[0],end2[0])
		#print("--- %s total seconds to calcualte the path---" % (time.time() - start_time))	
		temp=[]
		for k in path:
			temp.append(k)
		for j in range(1,len(temp)):
			dist+=distance(lonlat(*temp[j-1]), lonlat(*(temp[j]))).meters
		distout.append(dist)
		tout=length/(ave_travel_time*1000)
		timeout.append(tout)
		Dict['hits'] = len(end2)
		Dict['start'] = startcoords
		Dict['end'] = endcoords
		Dict['time'] = timeout
		Dict['distance'] = distout			
		Dict['paths']=temp
		if returnshape:
			srs = osr.SpatialReference()       ###
			srs.SetFromUserInput("EPSG:4326")  ###
			wgs84 = srs.ExportToProj4() 
			schema = { 'geometry': 'LineString','properties': {'dist': 'float','time': 'float'}}
			## Create shp file
			with fiona.open(shpOut, 'w', crs=wgs84, driver='ESRI Shapefile',schema=schema) as output:
				line = LineString(temp)	
				output.write({'geometry': mapping(line),'properties': { 'dist': dist, 'time': tout }
				})
			return Dict
		else:
			return Dict
		
	if len(end2)>1:
		length, path = nx.single_source_dijkstra(G, start2[0])
		if returnshape:#writes a shapefile
			srs = osr.SpatialReference()       ###
			srs.SetFromUserInput("EPSG:4326")  ###
			wgs84 = srs.ExportToProj4() 
			schema = { 'geometry': 'LineString','properties': {'dist': 'float','time': 'float'}}
			## Create shp file
			with fiona.open(shpOut, 'w', crs=wgs84, driver='GPKG',
					schema=schema) as output:
						for i in range(0,len(end2)):
							temp=[]	
							for k in path[end2[i]]:
								temp.append(k)
							pathouot.append(temp)
							dist=0
							for j in range(1,len(temp)):
								dist+=distance(lonlat(*temp[j-1]), lonlat(*(temp[j]))).meters
							distout.append(dist)
							tout=length[end2[i]]/(ave_travel_time*1000)
							timeout.append(tout)
							line = LineString(temp)	
							output.write({'geometry': mapping(line),'properties': { 'dist': dist, 'time': tout }
					})
			Dict['hits'] = len(end2)
			Dict['start'] = startcoords
			Dict['end'] = endcoords
			Dict['time'] = timeout
			Dict['distance'] = distout			
			Dict['paths']=pathouot
			return Dict
		else:
			for i in range(0,len(end2)):
				temp=[]
				for k in path[end2[i]]:
					temp.append(k)
				pathouot.append(temp)
				dist=0
				for j in range(1,len(temp)):
					dist+=distance(lonlat(*temp[j-1]), lonlat(*(temp[j]))).meters
				distout.append(dist)
				tout=length[end2[i]]/(ave_travel_time*1000)
				timeout.append(tout)
			Dict['hits'] = len(end2)
			Dict['start'] = startcoords
			Dict['end'] = endcoords
			Dict['time'] = timeout
			Dict['distance'] = distout			
			Dict['paths']=pathouot
			return Dict
	if len(start2)>1:#this is the function we will be using in the BEA project
		length, path = nx.single_source_dijkstra(G, end2[0])
#		print(G)
		#print(length,path)
		if returnshape:
			srs = osr.SpatialReference()       ###
			srs.SetFromUserInput("EPSG:4326")  ###
			wgs84 = srs.ExportToProj4() 
			schema = { 'geometry': 'LineString','properties': {'start_point_id': 'int','end_point_id': 'int','technology':'str', 'capactiy':'float','tonnes': 'float','service area (ha)':'float','Total dist (m)': 'float',
				'Arterial Road': 'float','Interchange': 'float', 'Main Road': 'float', 'National Freeway': 'float','National Road': 'float','None': 'float','On/OffRamp': 'float','Other Road': 'float','Secondary Road': 'float','Street': 'float','Virtual road': 'float','time': 'float'}}
			feat_type1=[]
			feat_type_distance1=[]
			features=['Arterial Road','Interchange', 'Main Road', 'National Freeway','National Road','None','On/OffRamp','Other Road','Secondary Road','Street','Virtual road']
			#create the path shapefile
			
			with fiona.open(shpOut, 'w', crs=wgs84, driver="GPKG",
					schema=schema) as output:
						for i in range(0,len(start2)):
							pathGraph = nx.path_graph(path[start2[i]])  # does not pass edges attributes
							path3=[]
							summ=0
							feat_type=[]
							feat_type_distance=[]
							w=[]
							for ea in pathGraph.edges():
								try:
									feat_type.append(G.edges[ea[0], ea[1]]['FEAT_TYPE'])
									w.append(G.edges[ea[0], ea[1]]['weight'])
								except:
									feat_type.append("Other")
									w.append(G.edges[ea[0], ea[1]]['weight'])
							feature_values=[0,0,0,0,0,0,0,0,0,0,vr_distance[i]]
							for ii in set(feat_type):
								summ=0
								for j in range(0,len(feat_type)):
									if feat_type[j]==ii:
										summ=summ+w[j]
								feat_type_distance.append(summ)
							accum=0
							for j in set(feat_type):
								for k in range(0,len(features)):
									if j==features[k]:
										feature_values[k]=feat_type_distance[accum]
										break
								accum=accum+1
							feat_type1.append(features)
							feat_type_distance1.append(feature_values)
							temp=[]	
							for k in path[start2[i]]:
								temp.append(k)
							temp = temp[::-1]
							if len(temp) <2:
								temp2=temp[0]
								temp.append(temp2)
							pathouot.append(temp)
							dist=0
							for j in range(1,len(temp)):
								dist+=distance(lonlat(*temp[j-1]), lonlat(*(temp[j]))).meters
							dist=dist+vr_distance[i]
							distout.append(dist)
							tout=length[start2[i]]/(ave_travel_time*1000)
							timeout.append(tout)
							line = LineString(temp)								
							mills2=s_point_ids[i]
							output.write({'geometry': mapping(line),'properties': { 'start_point_id':mills2, 'end_point_id': e_point_id,'technology':attr2, 'capactiy':attr3,'tonnes':tons[i]  ,'service area (ha)':service_area[i], 'Total dist (m)':dist, 'Arterial Road': feature_values[0],'Interchange': feature_values[1], 'Main Road': feature_values[2], 'National Freeway': feature_values[3],'National Road': feature_values[4],'None':feature_values[5],'On/OffRamp': feature_values[6],'Other Road': feature_values[7],'Secondary Road': feature_values[8],'Street': feature_values[9],'Virtual road':feature_values[10], 'time': tout}})	
			
			for i in range(0,len(start2)):
				pathGraph = nx.path_graph(path[start2[i]])  # does not pass edges attributes
				# Read attributes from each edge
				path3=[]
				summ=0
				feat_type=[]
				feat_type_distance=[]
				w=[]
				for ea in pathGraph.edges():
					try:
						feat_type.append(G.edges[ea[0], ea[1]]['FEAT_TYPE'])
						w.append(G.edges[ea[0], ea[1]]['weight'])
					except:
						feat_type.append("Other")
						w.append(G.edges[ea[0], ea[1]]['weight'])
						
				feature_values=[0,0,0,0,0,0,0,0,0,0,vr_distance[i]]
				for ii in set(feat_type):
					summ=0
					for j in range(0,len(feat_type)):
						if feat_type[j]==ii:
							summ=summ+w[j]
					feat_type_distance.append(summ)
				accum=0
				for j in set(feat_type):
					for k in range(0,len(features)):
						if j==features[k]:
							feature_values[k]=feat_type_distance[accum]
							break
					accum=accum+1
				feat_type1.append(features)
				feat_type_distance1.append(feature_values)
				
			Dict['hits'] = len(start2)
			Dict['feature_type']=feat_type1
			Dict['feature_type_distance']=feat_type_distance1
			Dict['start'] = startcoords
			Dict['end'] = endcoords
			Dict['time'] = timeout
			Dict['distance'] = distout			
			Dict['paths']=pathouot			
			return Dict
		else:
			for i in range(0,len(start2)):
				pathGraph = nx.path_graph(path[start2[i]])  # does not pass edges attributes
				# Read attributes from each edge
				path3=[]
				summ=0
				feat_type=[]
				feat_type_distance=[]
				w=[]
				for ea in pathGraph.edges():
					try:
						feat_type.append(G.edges[ea[0], ea[1]]['FEAT_TYPE'])
						w.append(G.edges[ea[0], ea[1]]['weight'])
					except:
						feat_type.append("Other")
						w.append(G.edges[ea[0], ea[1]]['weight'])
						
				for ii in set(feat_type):
					summ=0
					for j in range(0,len(feat_type)):
						if feat_type[j]==ii:
							summ=summ+w[j]
					feat_type_distance.append(summ)					
				temp=[]
				for k in path[start2[i]]:
					temp.append(k)
				temp = temp[::-1]
				pathouot.append(temp)
				dist=0
				for j in range(1,len(temp)):
					dist+=distance(lonlat(*temp[j-1]), lonlat(*(temp[j]))).meters
				dist=dist+vr_distance[i]
				distout.append(dist)
				tout=length[start2[i]]/(ave_travel_time*1000)
				timeout.append(tout)
				
			Dict['feature_type']=set(feat_type)
			Dict['feature_type_distance']=set(feat_type_distance)
			Dict['hits'] = len(start2)
			Dict['start'] = startcoords
			Dict['end'] = endcoords
			Dict['time'] = timeout
			Dict['distance'] = distout			
			Dict['paths']=pathouot
			return Dict
	
