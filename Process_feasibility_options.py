#!/Users/privateprivate/envs/bin/python
import os
import Process_feasibility_functions as func

time_horizon=20 #years of operation - updated from csv paramter files
inflation=0.06 #percent
#truck assumptions from here https://fleetwatch.co.za/wp-content/uploads/2019/09/pg-61-65.pdf
truck_working_days=260 #number of days a truck works in a year	
road_speed=50
offroad_speed=8
hours_per_day = 8
annual_km = 100000 #number of km the truck travells in a year
road_v_capacity = 26 #tons
truck_capital_cost=2000000	#R this value is not actually needed as it's included in the CPK
truck_standing_cost = 923352/ annual_km #R/km from https://fleetwatch.co.za/wp-content/uploads/2019/09/pg-61-65.pdf
truck_variable_cost = 1093536 /annual_km #R/km from https://fleetwatch.co.za/wp-content/uploads/2019/09/pg-61-65.pdf
truck_CPK = truck_standing_cost+truck_variable_cost #R/km or cost per kilometer
#need proper tracktor / offroad CPK values at some point
offroad_v_capacity = 4	#tons  https://sasri.org.za/storage/2019/07/Poster-truck-legal-payload.pdf
offroad_haulage=453.31 #rand per hour https://sasri.org.za/storage/2019/07/Mech-Report-No-1.pdf example 4 assumption
tractor_CPK=offroad_haulage/offroad_speed
#need to rename these to proper names at some point
km_road_factor = (road_speed*hours_per_day*truck_working_days)		#km/year@50km/hour # should be 100000 for 1 truck
km_offroad_factor = (offroad_speed*hours_per_day*truck_working_days)	#km/year@8km/hour
#cost parameters
capex = 0                   #user defined Parameter in R/ton of dry biomass - calculated from cost curve data
opex = 0                    #user defined parameter in R/ton of dry biomass - calculated from cost curve data
load_cost = 14.4			#ZAR/ton
unload_cost = 12.32			#ZAR/ton
harvest_cost = 0			#ZAR/ha 

#change accordingly 
cost_curves = '/Users/privateprivate/SAEON/Python_code/model_outputs/parameters/Technology_cost_curves.csv'
product_energy = '/Users/privateprivate/SAEON/Python_code/model_outputs/parameters/Biofuel energy densities.csv'
wd2 = '/Users/privateprivate/SAEON/Python_code/BEA_Feasibility_model_test/outputs/'
os.chdir(wd2)

feedstock='Aliens'
technology=['CBP-Pelleting']
capacity1=[[24000,48000,84000]]


for i in range(0,len(technology)):
	capa = capacity1[i]
	csv2 = wd2 + 'combined.csv'
	outputfile = wd2 + feedstock + '_' + str(technology[i]) + '_combined_aggregation.csv'
	#remove the files first in case
	if os.path.exists(outputfile):
		os.remove(outputfile)
	else:
		print("The file does not exist")
	if os.path.exists(csv2):
		os.remove(csv2)
	else:
		print("The file does not exist")
	#Combine .csv files into a single CSV with a unique ID
	func.munging_csv(wd2)
	#Combine .csv aggregations files into a single CSV with a unique ID
	func.do_aggreagation_csv(outputfile,csv2,cost_curves,product_energy,		load_cost,unload_cost,road_v_capacity,offroad_v_capacity,time_horizon,
		inflation,road_speed,offroad_speed,truck_CPK,tractor_CPK,km_road_factor,km_offroad_factor,truck_working_days)