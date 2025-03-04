Prepare files and config file (avaframeCfg.ini):
- Set path to avalancheDir in avaframeCfg.ini and create it
- Inputs:
	o All input files are automatically read from the set avalancheDir. No file paths need to be specified
	o elevation: DEM (ASCII), which serves as the basis for calculating the scarps. Must be in avalancheDir/Inputs.
	o perimeter: Perimeter file (ASCII). This specifies the boundary area within which scarps are calculated. Must be located in avalancheDir/Inputs/RASTERS and file name must end with “_per”.
	o shapefile (optional): A shapefile containing point geometries can be specified here. These points represent the centers of the ellipsoids or planes. If a shapefile is specified, the coordinates (x,y) of these points are used. Must be located in avalancheDir/Inputs/POINTS and file name must end with “_feat”.
	o perimeter_shapefile (optional): ): A shapefile can be specified here, which specifies a boundary area instead of the perimeter ASCII file. Must be located in avalancheDir/Inputs/POLYGONS and file name must end with “_per”.

-Output:
	o elevscarp: Output DGM (ASCII or GeoTIFF), which maps the input DGM minus the calculated scarp. Is saved under “raster_scarp.asc” in avalancheDir/Inputs/RASTERS.
	o hrelease: File path to the output DGM (ASCII or GeoTIFF), which represents the calculated scarp volumes.Is saved under “raster_rel.asc” in avalancheDir/Inputs/REL.

Prepare the config file (scarpCfg.ini):
- Input:
	o features: A list of parameters that define the planes and ellipsoids is specified here.The list is only used if no shapefile is specified.Depending on whether the planes or ellipsoid method is used, the list must have a different form:
		 Planes: features = x1,y1,z1,dip1,slope1,x2,y2,z2,dip2,...
		 Ellipsoid: features = x1,y1,depth1,a1,b1,x2,y2,...
- Settings:
	o method: Here you specify whether the plane or the ellipsoid method should be used
	o If a shapefile is specified, the following parameters must also be defined:
		 For the plane method:
			- Planes_zseed: z1,z2,...
			- Planes_dip: dip1, dip2,...
			- Planes_slope: slope1, slope2...
		 For the ellipsoid method: 
			- ellipsoids_max_depth: max_depth1, max_depth2,...
			- ellipsoids_semi_major: semi_major1, semi_major,...
			- ellipsoids_semi_minor: semi_minor1, semi_minor2,...
		 It is important to note that the number of parameters specified must correspond to the number of point features specified in the shapefile. The order must also be observed. For example, Planes_zseed z1 is assigned to the point with ID = 1. 

If all the data is provided successfully, start the script by running runCom6Scarp.py

