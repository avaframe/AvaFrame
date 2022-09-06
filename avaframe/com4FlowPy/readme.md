
# Flow-Py

You need to install GDAL and rasterio manually, it is not installed via the general AvaFrame requirements. 

If you have trouble installing GDAL or rasterio on Windows use these links to 
get the required version directly from their website, first install *GDAL* and then *rasterio*.

*GDAL*: https://www.lfd.uci.edu/~gohlke/pythonlibs/#gdal

*rasterio*: https://www.lfd.uci.edu/~gohlke/pythonlibs/#rasterio

#### MacOS 

com4FlowPy has not been tested it on MacOS. If you are able to run it there, please give us feedback.

<!-- ## Back-tracking extension -->

<!-- The back-tracking extension is an example of a custom built model extension used to identify the release areas, paths and deposition areas of GMF directly endangering infrastructure. -->

<!-- To activate the back-tracking extension an additional raster layer describing the spatial extent of the infrastructure must be called by the model. The GUI version of the model has a field where an infrastructure layer can be inserted. For the terminal version the “infra= path_to_infrastructure_raster” must be included as an additional argument (see command below). -->

<!-- ```markup -->
<!-- python3 main.py alpha_angle exponent working_directory path_to_dem path_to_release infra=path_to_infrastructure(Optional) flux_threshold=positiv_number(Optional) max_z_delta=positiv_number(Optional) -->
<!-- ``` -->

<!-- #### Example: -->

<!-- ```markup -->
<!-- python3 main.py 25 8 ./examples/dam/ ./examples/dam/dam_010m_standard_cr100_sw250_f2500.20.6_n0.asc ./examples/dam/release_dam.tif infra=./examples/dam/infra.tif flux=0.003 max_z=270 -->
<!-- ``` -->

<!-- The infrastructure layer must be in the same extent and resolution as the other input layers (DEM and release raster). Raster cells that contain infrastructure must have values > zero, raster cells with values = 0 represent locations without infrastructure (see infrastructure.tif in example folder). Different integers can be used to differentiate types of infrastructure, where higher valued infrastructure should correspond to higher integers . When a raster cell is associated with endangering infrastructure the integer associated with the infrastructure type is saved in the back-tracking output layer. When a raster cell is associated with endangering more than one type of infrastructure the larger integer is saved to the back tracking output.  -->


<!-- ### Back-tracking output: -->

<!-- - z_delta: the maximum z_delta of all paths for every raster cell (geometric measure of process magnitude, can be associated to kinetic energy/velocity) -->
<!-- - Flux: The maximum routing flux of all paths for every raster cell -->
<!-- - Flow Path Travel Angle, FP_TA: the gamma angle along the flow path -->
<!-- - Straight Line Travel Angle, SL_TA: Saves the travel angle along a straight line, i.e. distances are calculated via a direct line from the release cell to the current cell -->
<!-- - Back-tracking: Areas identified as endangering infrastructure.  -->


