# BioEnergy Atlas Modelling Script test example
 
@ Main authors: Marc Pienaar (marc@saeon.ac.za, marc.pienaar@gmail.com) and Hayden Wilson (hayden@saeon.ac.za)
Run the BEA_model_test.py file as it to see what the model does.
Outputs are generated in an ouput folder and a processed folder 

Quick description
-----------------

The BEA_model_test.py demonstrates a test case example for the SAEON Bioenergy Atlas feasibility model. It uses theoretical biomass values for cleared invasive areas for a small subset in KwaZulu Natal, South Africa as the example.

The inputs (in the input directory) are: 

	1) The road shapefile (BEAroads1.shp) subset from South Africa's National Geo-spatial Information (NGI), a component of Department of Rural Development and Land Reform (DRDLR) (http://www.ngi.gov.za)
	
	2) A polygon subset of alien invasives (KZN_resampled.shp) - an input centroid layer is generated in the script from this

In addition to the input files are some standard data files (i the data directory) used by the model (also subsets in this example), these are:

	1) Start and end vertices representing virtual roads at 1km intervals, and proximity rasters@100m for the road layer and the virtual road start points

Parameter files used to process the outputs of the feasibility model are included in the parameters directory, these are:

	1) BioEnergy_impacts.csv
	
	2) Biofuel energy densities.csv
	
	3) Feedstock_energy_densities.csv 
	
	4) Technology cost curves.csv

The model performs the following operations:

	1) Using the max volume from a volume field (in total tonnes) from the centroid input layer, it begins filling a facility capacity; if the volume from the max volume centroid is used up before the capacity of a facility is met, it will fetch more volume from the next closest neighbour (in coordinate space) until the capacity of a facility is met. Using all the input points (whose volume was used to fill a facility), the facility location is determined as a weighted average in coordinate space). This process is repeated until all volumes are depleted from the input centroid layer. This is done using the generate facilities function
	
	2) Determines if the new Facility location and its feedstock locations (the clustered centroids) need to snap to the closest road or closet virtual road start point, in which case it maps the distance from a virtual road to the closest road.
	
	3) Clips various inputs files to the extent of the points created in step 1 and 2 above (these include roads, virtual road start and end points, and proximity rasters for the Virtual roads and main roads)
	
	4) Adds these points as new vertices into the clipped road layer in order to perform routing later on  
	
	5) Fetches the polygon areas associated with each facility feedstock location using an unique ID attribute in each layer 
	
	6) Generates geopackage outputs of the facility (x,y point), its feedstock locations (x,y coordinates) and the catchment area (polygon layer)
	
	7) Performs a one to many routing using NetworkX Dijkstra algorithm from all the start locations to the facility (or end point)
	
	8) writes out a paths geopackage layer and creates and an associated csv file with all the parameters from the model (distances, volumes, areas etc.)
	
	9) Finally, once a full set of scenarios with different technology capacities has been run, it processes the technology csv from step 8 above into a processed csv file. 
