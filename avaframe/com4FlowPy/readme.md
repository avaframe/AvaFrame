
# Flow-Py

You need to install GDAL and rasterio manually, it is not installed via the general AvaFrame requirements. 

If you have trouble installing GDAL or rasterio on Windows use these links to 
get the required version directly from their website, first install *GDAL* and then *rasterio*.

*GDAL*: https://www.lfd.uci.edu/~gohlke/pythonlibs/#gdal

*rasterio*: https://www.lfd.uci.edu/~gohlke/pythonlibs/#rasterio

#### MacOS 

com4FlowPy has not been tested it on MacOS. If you are able to run it there, please give us feedback.

## Branch: delineation of paths and thalwegs

- Source: master of AvaFrame version 1.8
- Modified ini - File to Switch on if Pathanalyses (PathAnalysis = 1) and plots of thalwegs are saved (PathPlots = 1)
- order, when cells are calculated (*flowCore.py*) differs: remember generation of every parentcell and calculate first all parentcells of one generation

- **If PathAnalysis is switched on:**
    - new File: *flowPath.py*
    - calculate for every path a thalweg (center of Energy and center of Flux in every generation). 
    - save several data of every thalweg (for every startcell/path) and every path in a csv File (e.g. runout length of every thalweg)
- **If PathPlot is switched on:**
    - save a plot for every path/thalweg containing information about location, elevation and zdelta off thalweg,..