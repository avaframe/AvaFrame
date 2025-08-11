Prepare files and config file (avaframeCfg.ini):
- Set path to avalancheDir in avaframeCfg.ini and create it
- Inputs:
	o All input files are automatically read from the set avalancheDir. No file paths need to be specified
	o elevation: DEM (ASCII), which serves as the basis for calculating the scarps. Must be in avalancheDir/Inputs.
	o shapefile: A shapefile containing point geometries. These points represent the centers of the ellipsoids or planes. The coordinates (x,y) of these points are used. If the plane method is used, the shape file must contain the Attributes "zseed", "dip" and "slopeangle" as float values. If the ellipsoid method is used, the shape file must contain the Attributes "maxdepth", "semimajor", "semiminor", "tilt", "direc", "dip", "offset". The file must be located in avalancheDir/Inputs/POINTS and file name must end with “_coordinates”.
	o perimeter_shapefile: A shapefile that specifies a boundary area. Must be located in avalancheDir/Inputs/POLYGONS and file name must end with “_perimeter”.

-Output:
	o elevscarp: Output DGM (ASCII or GeoTIFF), which maps the input DGM minus the calculated scarp. Is saved under “raster_scarp.asc” in avalancheDir/Inputs/RASTERS.
	o hrelease: File path to the output DGM (ASCII or GeoTIFF), which represents the calculated scarp volumes.Is saved under “raster_rel.asc” in avalancheDir/Inputs/REL.

Prepare the config file (scarpCfg.ini):
- Input:
	o set useShapefiles = True 
- Settings:
	o method: Here you specify whether the plane or the ellipsoid method should be used
	
If all the data is provided successfully, start the script by running runCom6Scarp.py

Attribute meanings:
zseed: defines z coordinate of plane Center (m)
dip: direction in which the plane/slope is facing (degree)
slopeangle: steepness/angle of the slope (degree)

maxdepth: maximum depth of the ellipsoid (m)
semimajor: length of the major axis (m)
semiminor: length of the minor axis (m)
tilt: steepness/angle of the slope (degree)
direc: direction in which the slope is facing (degree)
dip: direction in which the ellipsoid is facing (degree)
Offset: offset, normal to the DEM slope (m)