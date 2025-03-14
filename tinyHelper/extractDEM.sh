#!/bin/bash

# Steps to take to generate a new topo manually (instead of script below):
# 1. go to https://www.data.gv.at/katalog/dataset/land-tirol_tirolgelnde#resources and get DGM 5m Tirol
# 2. Load zip into QGis, make sure the correct projection is used
# 3. Use Raster - Warp to reproject to epsg31287 (use resampling method: cubic) -> save as tiff
# 4. Draw the extend in a shapefile with epsg31287
# 5. Use Raster - Clip raster by mask layer to cut out the new DGM -> set .asc as output right away


# location of the folder containing the contours (AlrExtend.shp file for each avalanche in a seperate folder naved after the avalanche: Alr for exmple)
folder="tmp/topo/"
# location of the initial GEOtiff data
DEM="${folder}DGM_Tirol_5m_epsg31287.tif"

# move files to correct location
dest="tmp/topo/"

# Note that both .shp and geotiff should be in the same coord sys: epsg31287

# loop and extract tif from initial GEOtiff data and contours
for d in ${folder}**/*Extend.shp;do
	echo "extracting DEM from file: $d"
	filename=${d##*/}
	avaname=${filename%Extend*}
	# echo "$avaname"
	savename=${d%$filename}ava${avaname}.tif
        echo "to : $savename"
        # if output already exists remove it
        rm ${savename}
        gdalwarp -of GTiff -cutline "$d" -cl "${avaname}Extend" -crop_to_cutline "$DEM" "$savename"
        
    done
    
    
# loop & convert tif to asc recursively,
for j in ${folder}**/*.tif;do
	echo "converting file: $j"
        echo "to : ${j%.tif}.asc"
        # if output already exists remove it
        rm ${j%.tif}.asc
        # gdal_translate -of AAIGrid "$j" "${j%.tif}.asc"
        gdal_translate -co force_cellsize=true -of AAIGrid "$j" "${j%.tif}.asc"
        # remove the spurious file produced
        rm "${j%.tif}.asc.aux.xml"
        # rm "${j%.tif}.tif.aux.xml"
    done

# copy files to ava directory,
for j in ${folder}**/*.asc;do
	echo "moving file: $j"
	filename=${j##*/}
	avaname=${filename%.asc}
	destination=${dest}${avaname}/Inputs/${filename}
        echo "to : ${destination}"
        # if output already exists remove it
        rm ${destination}
        cp ${j} ${destination}
    done 



