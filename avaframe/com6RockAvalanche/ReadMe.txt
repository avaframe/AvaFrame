How to use com6scarp:

In the main configuration file avaframeCfg.ini: 
-set path to the scarpCfg configuration file in the MAIN section
-e.g: configFile = com6RockAvalanche/scarpCfg.ini

Prepare scarp configuration file:

- Input:
o elevation: File path to the DEM (ASCII or GeoTIFF), which serves as the basis for calculating the scarp. 
o perimeter: File path to the perimeter file (ASCII or GeoTIFF). This specifies the boundary area within which scarps are calculated. 
o shapefile (optional): The path to a shapefile containing point geometries can be specified here. These points represent the centers of the ellipsoids or 			planes. If a shapefile is specified, the coordinates (x,y) of these points are used. 
o perimeter_shapefile (optional): ): The path to a shapefile can be specified here, which specifies a boundary area instead of the perimeter ASCII file.
o features: A list of parameters that define the planes and ellipsoids is specified here. The list is only used if no shapefile is specified. Depending on 	    whether the planes or ellipsoid method is used, the list must have a different form:
	 Planes: features = x1,y1,z1,dip1,slope1,x2,y2,z2,dip2,...
	 Ellipsoid: features = x1,y1,depth1,a1,b1,x2,y2,...
- Output:
o elevscarp: File path to the output DEM (ASCII or GeoTIFF), which maps the input DGM minus the calculated scarps
o hrelease: File path to the output DEM (ASCII or GeoTIFF), which represents the calculated scarp volumes 

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
	 It is important to note that the number of parameters specified must correspond to the number of point features specified in the shapefile. The 	  order must also be observed. For example, Planes_zseed z1 is assigned to the point with ID = 1. 

If all the data is provided successfully, start the script by running runCom6Scarp.py