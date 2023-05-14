from osgeo import gdal
import os
import pathlib

#Change cell size 

def Cellsize_change(avaDir, cellsize):
    """
    Changes the cellsize which can be defined in the r.avaflow config files,
    for all the files relevant for a r.avaflow simulations.

    Parameters
    ----------
    avaDir : str or pathlib Path
        path to avalanche data
    cellsize : int
        cellsize derived through the config file for r.avaflow

    """
    
    #1. Adapt cellsize for release area      
    input_file_path_tif = os.path.join(avaDir, 'Ravaflow_Input', 'rel_rst.tif')

    if os.path.exists(input_file_path_tif):
        # Adapt cellsize for release area
        output_file_path_tif = os.path.join(avaDir, 'Ravaflow_Input', 'rel_rst.tif')
        gdal.Warp(output_file_path_tif, input_file_path_tif, xRes=cellsize, yRes=cellsize)
    
        # .tif to ascii
        input_file_path_asc = os.path.join(avaDir, 'Ravaflow_Input', 'rel_rst.tif')
        output_file_path_asc = os.path.join(avaDir, 'Ravaflow_Input', 'rel_rst.asc')
    
        # Generate string of process.
        gdal_translate_str = 'gdal_translate -of AAIGrid {0} {1}'
        gdal_translate_process = gdal_translate_str.format(input_file_path_asc, output_file_path_asc)
    
        # Call process.
        os.system(gdal_translate_process)
    
    
    #2. Adapt cellsize for DEM
    input_file_path = os.path.join(avaDir, 'Ravaflow_Input', 'DEM.tif')
    
    if os.path.exists(input_file_path):
    
        output_file_path = os.path.join(avaDir, 'Ravaflow_Input', 'DEM.tif')
        gdal.Warp(output_file_path, input_file_path, xRes=cellsize, yRes=cellsize)
        
         ### .tif to ascii
        input_file_path = os.path.join(avaDir, 'Ravaflow_Input', 'DEM.tif')
        output_file_path = os.path.join(avaDir, 'Ravaflow_Input', 'DEM.asc')
        
        # Generate string of process.
        gdal_calc_str = 'gdal_translate -of AAIGrid {0} {1}'
        gdal_calc_process = gdal_calc_str.format(input_file_path, output_file_path)
        # Call process.
        os.system(gdal_calc_process)
        
    
    #3. Adapt cellsize for secondary release area 
    inputDir_rel = pathlib.Path(avaDir, 'Inputs/SECREL')
    shpFile = list(inputDir_rel.glob('*.shp'))

    if len(shpFile) == 1:
        input_file_path = os.path.join(avaDir, 'Ravaflow_Input', 'rel2_rst.tif')
        output_file_path = os.path.join(avaDir, 'Ravaflow_Input', 'rel2_rst.tif')
        gdal.Warp(output_file_path, input_file_path, xRes=cellsize, yRes=cellsize)
        
         ### .tif to ascii
        input_file_path = os.path.join(avaDir, 'Ravaflow_Input', 'rel2_rst.tif')
        output_file_path = os.path.join(avaDir, 'Ravaflow_Input', 'rel2_rst.asc')
        
        # Generate string of process.
        gdal_calc_str = 'gdal_translate -of AAIGrid {0} {1}'
        gdal_calc_process = gdal_calc_str.format(input_file_path, output_file_path)
        # Call process.
        os.system(gdal_calc_process)
    
    #4. Adapt cellsize for entrainment area 
    inputDir_rel = pathlib.Path(avaDir, 'Inputs/ENT')
    shpFile = list(inputDir_rel.glob('*.shp'))

    if len(shpFile) == 1:
        input_file_path = os.path.join(avaDir, 'Ravaflow_Input', 'ent_rst.tif')
        output_file_path = os.path.join(avaDir, 'Ravaflow_Input', 'ent_rst.tif')
        gdal.Warp(output_file_path, input_file_path, xRes=cellsize, yRes=cellsize)
        
         ### .tif to ascii
        input_file_path = os.path.join(avaDir, 'Ravaflow_Input', 'ent_rst.tif')
        output_file_path = os.path.join(avaDir, 'Ravaflow_Input', 'ent_rst.asc')
        
        # Generate string of process.
        gdal_calc_str = 'gdal_translate -of AAIGrid {0} {1}'
        gdal_calc_process = gdal_calc_str.format(input_file_path, output_file_path)
        # Call process.
        os.system(gdal_calc_process)
