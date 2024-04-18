
# Flow-Py

You need to install GDAL and rasterio manually, it is not installed via the general AvaFrame requirements. 

If you have trouble installing GDAL or rasterio on Windows use these links to 
get the required version directly from their website, first install *GDAL* and then *rasterio*.

*GDAL*: https://www.lfd.uci.edu/~gohlke/pythonlibs/#gdal

*rasterio*: https://www.lfd.uci.edu/~gohlke/pythonlibs/#rasterio

#### MacOS 

com4FlowPy has not been tested it on MacOS. If you are able to run it there, please give us feedback.



#### Branch: forest interaction

- Source: more/less actual master branch
- Modified ini - File to switch on (forestInteraction = 1) or off (forestInteraction = 0) the calculation
and additional output of a forest interaction layer
- **If forestInteraction is switched on:**
	- Additional input: forest layer (binary: 0: no forest, 1: forest) is required in same extent and resolution as DEM
(further extensions for forest density necessary???)
	- in *com4FlowPy.py* and *flowCore.py*: Read in and hand over forest & get forest_interaction variable from 
	*flowClass* and hand out minimum from all paths 
	- in *flowClass.py*: calculate forest interaction for each cell (in each path)
- The output *Forest_interaction.tif* shows how many forested cells a path passed previously. In the case of 
overlying paths, the minimum of the interaction value is taken

- Test cases:
	- 'avaInclinedPlane' (generic)
	- Davos_Jakobshorn_einHnag_test (real)
