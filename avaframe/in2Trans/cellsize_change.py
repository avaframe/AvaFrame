from osgeo import gdal
from osgeo import ogr
import os
import pathlib

#Change cell size 

def Cellsize_change(avaDir, cellsize):
    
    #1. Adapt cellsize for release area  
    input_file_path = os.path.join(avaDir, 'Avaflow_Input', 'rel_rst.tif')
    output_file_path = os.path.join(avaDir, 'Avaflow_Input', 'rel_rst.tif')
    gdal.Warp(output_file_path, input_file_path, xRes=cellsize, yRes=cellsize)
    
     ### .tif to ascii
    input_file_path = os.path.join(avaDir, 'Avaflow_Input', 'rel_rst.tif')
    output_file_path = os.path.join(avaDir, 'Avaflow_Input', 'rel_rst.asc')
    
    # Generate string of process.
    gdal_calc_str = 'gdal_translate -of AAIGrid {0} {1}'
    gdal_calc_process = gdal_calc_str.format(input_file_path, output_file_path)
    # Call process.
    os.system(gdal_calc_process)
    
    
    #2. Adapt cellsize for DEM
    input_file_path = os.path.join(avaDir, 'Avaflow_Input', 'DEM.tif')
    output_file_path = os.path.join(avaDir, 'Avaflow_Input', 'DEM.tif')
    gdal.Warp(output_file_path, input_file_path, xRes=cellsize, yRes=cellsize)
    
     ### .tif to ascii
    input_file_path = os.path.join(avaDir, 'Avaflow_Input', 'DEM.tif')
    output_file_path = os.path.join(avaDir, 'Avaflow_Input', 'DEM.asc')
    
    # Generate string of process.
    gdal_calc_str = 'gdal_translate -of AAIGrid {0} {1}'
    gdal_calc_process = gdal_calc_str.format(input_file_path, output_file_path)
    # Call process.
    os.system(gdal_calc_process)
    
    
    #3. Adapt cellsize for secondary release area 
    inputDir_rel = pathlib.Path(avaDir, 'Inputs/SECREL')
    shpFile = list(inputDir_rel.glob('*.shp'))

    if len(shpFile) == 1:
        input_file_path = os.path.join(avaDir, 'Avaflow_Input', 'rel2_rst.tif')
        output_file_path = os.path.join(avaDir, 'Avaflow_Input', 'rel2_rst.tif')
        gdal.Warp(output_file_path, input_file_path, xRes=cellsize, yRes=cellsize)
        
         ### .tif to ascii
        input_file_path = os.path.join(avaDir, 'Avaflow_Input', 'rel2_rst.tif')
        output_file_path = os.path.join(avaDir, 'Avaflow_Input', 'rel2_rst.asc')
        
        # Generate string of process.
        gdal_calc_str = 'gdal_translate -of AAIGrid {0} {1}'
        gdal_calc_process = gdal_calc_str.format(input_file_path, output_file_path)
        # Call process.
        os.system(gdal_calc_process)
    
    #4. Adapt cellsize for entrainment area 
    inputDir_rel = pathlib.Path(avaDir, 'Inputs/ENT')
    shpFile = list(inputDir_rel.glob('*.shp'))

    if len(shpFile) == 1:
        input_file_path = os.path.join(avaDir, 'Avaflow_Input', 'ent_rst.tif')
        output_file_path = os.path.join(avaDir, 'Avaflow_Input', 'ent_rst.tif')
        gdal.Warp(output_file_path, input_file_path, xRes=cellsize, yRes=cellsize)
        
         ### .tif to ascii
        input_file_path = os.path.join(avaDir, 'Avaflow_Input', 'ent_rst.tif')
        output_file_path = os.path.join(avaDir, 'Avaflow_Input', 'ent_rst.asc')
        
        # Generate string of process.
        gdal_calc_str = 'gdal_translate -of AAIGrid {0} {1}'
        gdal_calc_process = gdal_calc_str.format(input_file_path, output_file_path)
        # Call process.
        os.system(gdal_calc_process)
