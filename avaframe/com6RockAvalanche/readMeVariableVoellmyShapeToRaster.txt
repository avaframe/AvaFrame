The VariableVoellmyShapeToRaster.py script allows the user to define spatially different values for the voellmy Parameters mu and xsi, with the use of polygon shapefiles. For the extent of a DEM raster, all the areas that are not covered by a polygon get assigned a default mu or xsi value. The script then converts this Information into a raster mu and a raster xsi file, which can then be used in Avaframe Simulation runs, using the "spatialVoellmy" friction model. 

First, set up the Config File and provide inputs:
oInputs:
	dem: Path to the DEM Raster that is later on used for the avaframe simulation. This is needed, because the mu and xsi output rasters need to be the 			exact same size. Usually this file lies in \AvaFrame\avaframe\data\ *yourAvalancheDir*\Inputs
	mu_shapefile: Path to the mu shapefile, that is then converted to a raster file. Be aware, that the attribute has to be named “xsi”. 
	xsi_shapefile: Path to the xsi shapefile, that is then converted to a raster file. Be aware, that the attribute has to be named “mu”. 
oDefaults:
	default_mu: this is the default mu value, that gets assigned to all areas in the raster, that are not covered by shapefile-polygons
	default_xsi: this is the default xsi value, that gets assigned to all areas in the raster, that are not covered by shapefile-polygons
	
RunScript:
	oOnce everything is set up, run the script “runVariableVoellmyShapeToRaster.py”
	oIf libraries are missing use: pip install *name of missing library
