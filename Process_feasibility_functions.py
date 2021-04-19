#!/Users/privateprivate/envs/bin/python
import os
import pandas as pd
import numpy as np
import math
import csv

'''
Saeon Bioenergy Atlas model cost utility script
@ Main authors: Marc Pienaar (marc@saeon.ac.za, marc.pienaar@gmail.com) and Hayden Wilson (hayden@saeon.ac.za)
'''

def munging_csv(wd):
	print(wd)
	#Combine .csv files into a single CSV with a unique ID
	li = []
	for subdir, dirs, files in os.walk(wd):
		for filename in files:
			
			filepath = subdir + os.sep + filename
			if filepath.endswith(".csv"):
				df = pd.read_csv(filepath, index_col=None, header=0)
				li.append(df)
				frame = pd.concat(li, axis=0, ignore_index=True)
				frame['unique_id'] = frame['feedstock'] + "_" + frame['technology'] + "_" +  frame['capacity'].astype(str) + '_'+ frame['end_point_id'].astype(str)
				frame.to_csv(wd + os.sep + 'combined.csv',index=False)
				
def do_aggreagation_csv(outputfile,csv2,cost_curves,product_energy,load_cost,unload_cost,road_v_capacity,offroad_v_capacity,time_horizon,inflation,road_speed,offroad_speed,truck_CPK,tractor_CPK,km_road_factor,km_offroad_factor,truck_working_days):
	#Read in the cost curve data - if not using the info provided in the params at the start
	cost_df = pd.read_csv(cost_curves, skiprows = 0)
	productE_df=pd.read_csv(product_energy, skiprows = 0)	
	df = pd.read_csv(csv2)
	#get unique_ids
	unique_id = df.unique_id.unique()
	#get colum headers
	Dict=aggregate_facility(cost_df,productE_df,df,unique_id[0],
		load_cost,unload_cost,road_v_capacity,offroad_v_capacity,time_horizon,
		inflation,road_speed,offroad_speed,truck_CPK,tractor_CPK,km_road_factor,km_offroad_factor,truck_working_days)
		#write headers to a csv
	with open(outputfile, 'w', newline='') as file:
		writer = csv.writer(file)
		writer.writerow(list(Dict.keys()))
		#now loop through and append the values to the csv
	for i in range(0,len(unique_id)):
		print(i+1,'of', len(unique_id),'done')
		Dict=aggregate_facility(cost_df,productE_df,df,
			unique_id[i],load_cost,unload_cost,road_v_capacity,offroad_v_capacity,time_horizon,
			inflation,road_speed,offroad_speed,truck_CPK,tractor_CPK,km_road_factor,km_offroad_factor,truck_working_days)
		update_csv(outputfile,Dict)
		
def aggregate_facility(cost_df,productE_df,df,id_no,load_cost,unload_cost,road_v_capacity,offroad_v_capacity,time_horizon,inflation,road_speed,offroad_speed,truck_CPK,tractor_CPK,km_road_factor,km_offroad_factor,truck_working_days):
	'''
	@ Author Marc Pienaar and Hayden Wilson
	A small utility function to create a aggreagate production cost output csv. it parses numerous user defined paramters
		
	Returns
	-------
	an aggregated production cost csv file 
	'''
	
	temp=pd.DataFrame()#create an empty dataframe to fill
	#get unique id rows
	unique=df[df.unique_id==id_no]
	#get paramters
	curve = cost_df[cost_df.ID.str.match(unique.technology.iloc[0],case=False)]
#	print(curve , unique.technology.iloc[0])
	time_horizon=curve.lifespan_years.iloc[0]	
	output_product=curve.Prod_id.iloc[0]
	Product_Energy_info = productE_df[productE_df.Prod_id.str.match(curve['Prod_id'].iloc[0],case=False)]
	conversion_efficiency = (curve['Efficiency'].iloc[0])/100
	Energy_Density = Product_Energy_info['Energy_(MJ/kg)'].iloc[0]
	#get capex and opex
	check=[curve.capacity1.iloc[0],curve.capacity2.iloc[0],curve.capacity3.iloc[0]]
	capex_array=[curve.capex1.iloc[0],curve.capex1.iloc[0],curve.capex1.iloc[0]]
	opex_array=[curve.opex1.iloc[0],curve.opex1.iloc[0],curve.opex1.iloc[0]]
	index_min = np.argmin(min(check/unique.capacity.iloc[0]))
	capex=capex_array[index_min]
	opex=opex_array[index_min]
	#get some values into a temp dataframe
	temp['tonnes_contributed'] = unique['tonnes_contributed'] 
	temp['capacity'] = unique['capacity'] 
	temp['load_cost'] = unique['tonnes_contributed'] * load_cost
	temp['unload_cost'] = unique['tonnes_contributed'] * unload_cost
	temp['road_trips'] = (unique['tonnes_contributed']/road_v_capacity)*2#we need a round trip 
	temp['VR_road_trips'] = (unique['tonnes_contributed']/offroad_v_capacity)*2#we need a round trip 
	#distances travelled
	temp['road_trip_dist'] = (temp['road_trips'] * (unique['Arterial_Road_(m)'] + unique['Interchange_(m)'] + unique['Main_Road_(m)'] + unique['National_Freeway_(m)']+unique['National_Road_(m)']+unique['None_(m)']+unique['On/OffRamp_(m)']+unique['Other_Road_(m)']+unique['Secondary_Road_(m)']+unique['Street_(m)']))/1000#convert to km
	#individual distances
	temp['Arterial_Road_(m)'] = (temp['road_trips'] * (unique['Arterial_Road_(m)']))/1000
	temp['Interchange_(m)'] = (temp['road_trips'] * (unique['Interchange_(m)']))/1000
	temp['Main_Road_(m)'] = (temp['road_trips'] * (unique['Main_Road_(m)']))/1000
	temp['National_Freeway_(m)'] = (temp['road_trips'] * (unique['National_Freeway_(m)']))/1000
	temp['National_Road_(m)'] = (temp['road_trips'] * (unique['National_Road_(m)']))/1000
	temp['None_(m)'] = (temp['road_trips'] * (unique['None_(m)']))/1000
	temp['On/OffRamp_(m)'] = (temp['road_trips'] * (unique['On/OffRamp_(m)']))/1000
	temp['Other_Road_(m)'] = (temp['road_trips'] * (unique['Other_Road_(m)']))/1000
	temp['Secondary_Road_(m)'] = (temp['road_trips'] * (unique['Secondary_Road_(m)']))/1000
	temp['Street_(m)'] = (temp['road_trips'] * (unique['Street_(m)']))/1000	
	temp['virtual_trip_dist'] = (temp['VR_road_trips']* unique['Virtual_road_(m)'])/1000
	#sum values
	a_costs=temp.sum(axis = 0, skipna = True)
	a_costs2=temp.mean(axis = 0, skipna = True)	
	#create dictionary object to so some calculations on
	my_dict=a_costs.to_dict() #get a dict
	my_dict_temp=a_costs2.to_dict() #get a dict
	my_dict2={}
	#fill with paramter values
	my_dict2['unique_id']=unique.unique_id.iloc[0]	
	my_dict2['end_point_id']=unique.end_point_id.iloc[0]
	my_dict2['feedstock']=unique.feedstock.iloc[0]
	my_dict2['technology']=unique.technology.iloc[0]
	my_dict2['capacity']=my_dict_temp['capacity']	
	my_dict2['output_product']=curve.Prod_id.iloc[0]
	my_dict2['X_coords']=unique.X_coords.iloc[0]
	my_dict2['Y_coords']=unique.Y_coords.iloc[0]
	my_dict2['time_horizon']=time_horizon
	my_dict2['inflation']=inflation	
	my_dict2['road_speed']=road_speed
	my_dict2['offroad_speed']=offroad_speed
	my_dict2['truck_CPK']=truck_CPK
	my_dict2['tractor_CPK']=tractor_CPK
	my_dict2['number_of_road_trucks_anum']=math.ceil(my_dict['road_trip_dist']/km_road_factor)
	my_dict2['number_of_vr_trucks_anum']=math.ceil(my_dict['virtual_trip_dist']/km_offroad_factor)	
	my_dict2['number_of_potential_road_trips_anum']=my_dict['road_trips']
	my_dict2['number_of_potential_vr_trips_anum']=my_dict['VR_road_trips']
	if my_dict['road_trip_dist']==0:
		my_dict2['number_of_actual_road_trips_anum']=0
		my_dict2['number_of_actual_road_trips_wd']=0
		my_dict2['km_wd_truck_road']=0
	else:
		my_dict2['number_of_actual_road_trips_anum']=my_dict['road_trips']
		my_dict2['number_of_actual_road_trips_wd']=my_dict['road_trips']/truck_working_days
		my_dict2['km_wd_truck_road']=(my_dict['road_trip_dist']/my_dict2['number_of_road_trucks_anum'])/truck_working_days
		
	if my_dict['virtual_trip_dist']==0:
		my_dict2['number_of_actual_vr_trips_anum']=0
		my_dict2['number_of_actual_vr_trips_wd']=0
		my_dict2['km_wd_truck_offroad']=0
	else:
		my_dict2['number_of_actual_vr_trips_anum']=my_dict['VR_road_trips']
		my_dict2['number_of_actual_vr_trips_wd']=my_dict['VR_road_trips']/truck_working_days
		my_dict2['km_wd_truck_offroad']=(my_dict['virtual_trip_dist']/my_dict2['number_of_vr_trucks_anum'])/truck_working_days
	#########################update with average lifetime costs
	my_dict2['Load_cost']=AVE_LTC(my_dict['load_cost'],inflation,time_horizon)
	my_dict2['Unload_cost']=AVE_LTC(my_dict['unload_cost'],inflation,time_horizon)	
	my_dict2['Load_cost_anum']=my_dict2['Load_cost']/time_horizon
	my_dict2['Unload_cost_anum']=my_dict2['Unload_cost']/time_horizon	
	my_dict2['load_cost_anum_pv']=pv(my_dict2['Load_cost_anum'],inflation,time_horizon)
	my_dict2['unload_cost_anum_pv']=pv(my_dict2['Unload_cost_anum'],inflation,time_horizon)	
	my_dict2['load_cost_tonne_anum']=AVE_LTC(load_cost,inflation,time_horizon)/time_horizon
	my_dict2['unload_cost_tonne_anum']=AVE_LTC(unload_cost,inflation,time_horizon)/time_horizon	
	my_dict2['load_cost_tonne_anum_pv']=pv(my_dict2['load_cost_tonne_anum'],inflation,time_horizon)
	my_dict2['unload_cost_tonne_anum_pv']=pv(my_dict2['unload_cost_tonne_anum'],inflation,time_horizon)
	#do some calculations
	my_dict2['tonnes_contributed']=my_dict['tonnes_contributed']*time_horizon
	my_dict2['tonnes_contributed_anum']=my_dict['tonnes_contributed']
	my_dict2['product_volume'] = (my_dict2['tonnes_contributed']*conversion_efficiency) #Volume OF bioEnergy product per lifetime 
	my_dict2['product_volume_anum'] = my_dict2['product_volume']/time_horizon #Volume OF bioEnergy product per anum 
	my_dict2['Capex']=my_dict['tonnes_contributed']*capex
	my_dict2['Opex']=(my_dict['tonnes_contributed']*opex)*time_horizon#here we assume opex is already a mean valu
	my_dict2['Opex_anum']=(my_dict['tonnes_contributed']*opex)
	my_dict2['Opex_anum_pv']=pv(my_dict2['Opex_anum'],inflation,time_horizon)
	###cost of transport
	###distances
	my_dict2['Arterial_Road_(km_anum)'] = my_dict['Arterial_Road_(m)']
	my_dict2['Interchange_(km_anum)'] = my_dict['Interchange_(m)']
	my_dict2['Main_Road_(km_anum)'] = my_dict['Main_Road_(m)']
	my_dict2['National_Freeway_(km_anum)'] = my_dict['National_Freeway_(m)']
	my_dict2['National_Road_(km_anum)'] = my_dict['National_Road_(m)']
	my_dict2['Unspecified_(km_anum)'] =my_dict['None_(m)']
	my_dict2['On/OffRamp_(km_anum)'] = my_dict['On/OffRamp_(m)']
	my_dict2['Other_Road_(km_anum)'] = my_dict['Other_Road_(m)']
	my_dict2['Secondary_Road_(km_anum)'] = my_dict['Secondary_Road_(m)']
	my_dict2['Street_(km_anum)'] = my_dict['Street_(m)']
	my_dict2['Road_dist_(km_anum)']=my_dict['road_trip_dist']
	my_dict2['VR_dist_(km_anum)']=my_dict['virtual_trip_dist']
	my_dict2['Total_dist_(km_anum)']=my_dict['virtual_trip_dist']+my_dict['road_trip_dist']
	##
	my_dict2['road_costs']=AVE_LTC((my_dict2['number_of_road_trucks_anum']*(truck_CPK*my_dict['road_trip_dist'])),inflation,time_horizon)
	my_dict2['road_costs_anum']=my_dict2['road_costs']/time_horizon
	my_dict2['road_costs_anum_pv']=pv(my_dict2['road_costs_anum'],inflation,time_horizon)	
	my_dict2['VR_costs']=AVE_LTC((my_dict2['number_of_vr_trucks_anum']*(tractor_CPK*my_dict['virtual_trip_dist'])),inflation,time_horizon)
	my_dict2['VR_costs_anum']=my_dict2['VR_costs']/time_horizon
	my_dict2['VR_costs_anum_pv']=pv(my_dict2['VR_costs_anum'],inflation,time_horizon)
	#total costs etc
	my_dict2['Total_costs_raw']=my_dict2['Capex']+my_dict2['Opex']+my_dict2['Load_cost']+my_dict2['Unload_cost']+my_dict2['road_costs']+my_dict2['VR_costs']	
	my_dict2['Total_costs_raw_anum']=my_dict2['Capex']+(my_dict2['Opex']+my_dict2['Load_cost_anum']+my_dict2['Unload_cost_anum']+my_dict2['road_costs_anum']+my_dict2['VR_costs_anum'])	
	my_dict2['Total_costs_raw_anum_pv']=my_dict2['Capex']+(my_dict2['Opex_anum_pv']+my_dict2['load_cost_anum_pv']+my_dict2['unload_cost_anum_pv']+my_dict2['road_costs_anum_pv']+my_dict2['VR_costs_anum_pv'])
	my_dict2['Total_costs_per_raw_ton']=my_dict2['Total_costs_raw_anum']/(my_dict['tonnes_contributed'])
	my_dict2['Total_costs_per_raw_ton_pv']=my_dict2['Total_costs_raw_anum_pv']/(my_dict['tonnes_contributed'])	
	my_dict2['energy_MJ_kg'] = ((my_dict2['product_volume']*1000)*Energy_Density) #Energy Equivalent 
	my_dict2['energy_mj_kg_anum'] = ((my_dict2['product_volume_anum']*1000)*Energy_Density) #Energy Equivalent 	
	my_dict2['energy_per_raw_tonne']=my_dict2['energy_mj_kg_anum']/my_dict2['tonnes_contributed_anum']
	my_dict2['energy_generation_Rand_kwh']=my_dict2['Total_costs_per_raw_ton']/my_dict2['energy_per_raw_tonne']
	return my_dict2

def AVE_LTC(pv,r,n):
	#pv = present value
	#r = annual interest rate
	#n= number of terms
	fv=pv*(1+r)**n
	total_payment=((fv+pv)/2)*n
	return total_payment

def pv(fv,r,n):
	return (2*fv)/(1+(1+r)**n)

def update_csv(outputfile,Dict):
	with open(outputfile, 'a', newline='') as file:
		writer = csv.writer(file)
		writer.writerow(list(Dict.values()))