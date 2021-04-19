# BioEnergy Atlas Modelling Script test example
 
@ Main authors: Marc Pienaar (marc@saeon.ac.za, marc.pienaar@gmail.com) and Hayden Wilson (hayden@saeon.ac.za)
Run the BEA_model_test.py file to see what the model does.
Outputs are generated in an output folder

Quick description
-----------------

The BEA_model_test.py demonstrates a test case example for the SAEON Bioenergy Atlas generalisable feasibility model. It uses theoretical biomass values for cleared invasive areas for a small subset in KwaZulu Natal, South Africa as an example.

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

![Overview_thread_marc_pienaar](https://user-images.githubusercontent.com/50328370/115220452-d9c04980-a108-11eb-9031-653f38ad51e5.png)
**Figure 1**. Overview of the BioEnergy Atlas generalisable feasibility model. The model is designed so that each model component or ‘thread’ can produce a set of analytical outputs that can be used as standalone products or can be connected to each other to produce further output products. The advantage of this approach is that not all model components need to be re-run in order to produce different model scenarios. The flowchart on the left visually summarises the flow of information through a series of connected model threads collectively called the BioEnergy Atlas generalisable feasibility model. The schematic on the right provides a conceptual overview of how the model works. In this schematic example, various inputs (scenario parameters, such as feedstock data and conversion technology data) along with transportation data (road networks and other transportation  parameters) are parsed to a facility allocation location component (the biomass conversion facility allocation thread) to determine the optimal location of a biomass process facilities. This is done by first finding the feedstock location with the most biomass volume, then depleting its volume followed by its nearest neighbours volume until a facility capacity has been met. Once a facility has enough feedstock, its location is calculated using a weighted geographic coordinate where the weights are feedstock volume. The outputs from this thread provide the start and end points that can be used for each biomass processing facility to calculate the transportation costs (the biomass acquisition thread figure) using a custom coded version of the python open source Network package. The outputs from this thread are saved as both shapefiles and tabular csv files, which are used by the biomass conversion cost thread and the comparisons thread to provide additional outputs, comparisons and analysis. 

![generate_facilities_thread_marc_pienaar](https://user-images.githubusercontent.com/50328370/115220932-6f5bd900-a109-11eb-80ec-d747dd33e8ea.png)

**Figure 2**. The biomass conversion facility location allocation thread. This algorithm is used to cluster feedstock locations according to the location and volume to generate an optimal location of a biomass processing technology facility. Using the max volume from a volume field (in total tonnes) from a centroid input layer representing feedstock locations, it begins filling a facility capacity; if the volume from the max volume centroid location is used up before the capacity of a facility is met, it will fetch more volume from the next closest neighbour (in coordinate space) until the capacity of a facility is met. Using all the input points (whose volume was used to fill a facility), the facility location is determined as a weighted average in coordinate space). This process is repeated until all volumes are depleted from the input centroid layer. Once all the volume in the feedstock input layer has been depleted the function ends. During the loop operation it parses various inputs to the biomass acquisition thread (see next diagram) to perform the transportation routing 

![Biomass_aquisition_thread_marc_piennaar](https://user-images.githubusercontent.com/50328370/115221401-eabd8a80-a109-11eb-9ace-4275af93857c.png)

**Figure 3**. The biomass acquisition thread. This thread forms part of a subprocess of the biomass conversion facility location allocation thread. It determines if the new Facility location and its feedstock locations need to snap to the closest road or closet virtual road start point, in which case it maps the distance from a virtual road to the closest road. During the operation it clips various inputs files to the extent of the points created in the biomass conversion facility location allocation thread (these include roads, virtual road start and end points, and proximity raster’s for the Virtual roads and main roads). It then Adds the temporary start and end points from the the biomass conversion facility location allocation thread as new vertices in the clipped road layer in order to create a graph and perform the transport routing. Transport modelling is done using a one to many routing algorithm (Dijkstra’s algorithm) using a custom coded version of the python NetworkX package. During the run of a full set of scenarios with different technology capacities from the biomass conversion facility location allocation thread, it appends information a technology csv and writes out various shapefiles. These outputs are used in the SAEON BEA but not described in further detail here.
