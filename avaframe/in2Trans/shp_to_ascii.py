#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from osgeo import gdal
from osgeo import ogr
import os
import pathlib

###---------------------
### Shapefile of the release area needs the same size as the DEM
### 1. Make .asc DEM to perfect rectangel -> output is .tif
### 2. Write .tif DEM to ascii
### 3. Rasterize release area shapefile -> output is .tif
### 4. Write .tif release to ascii

def shptoasc(avaDir):
    ## Make .asc DEM to perfect rectangel
    gdal_path = 'gdal_translate'
    
    search_path = os.path.join(avaDir, 'Avaflow_Input')
    for root, dir, files in os.walk(search_path):
        for file in files:
            if file.endswith(".asc"):
                ascFile = os.path.join(root, file)
    input_file_path = ascFile
    output_file_path = os.path.join(avaDir, 'Avaflow_Input', 'DEM.tif')
    # Generate string of process.
    gdal_str = '{0} -of GTIFF -projwin ulx uly lrx lry {1} {2}'
    gdal_calc_process = gdal_str.format(gdal_path, input_file_path, 
        output_file_path)
    # Call process.
    os.system(gdal_calc_process)
    
    ##tiff to ascii
    input_file_path = os.path.join(avaDir, 'Avaflow_Input', 'DEM.tif')
    output_file_path = os.path.join(avaDir, 'Avaflow_Input', 'DEM.asc')
    # Generate string of process.
    gdal_calc_str = 'gdal_translate -of AAIGrid {0} {1}'
    gdal_calc_process = gdal_calc_str.format(input_file_path, output_file_path)
    # Call process.
    os.system(gdal_calc_process)
    
    
    ## Rasterize release area shapefile
    #fn_ras = ascFile 
    
    search_path = os.path.join(avaDir, 'Avaflow_Input')
    for root, dir, files in os.walk(search_path):
        for file in files:
            if file.endswith(".shp"):
                shpFile = os.path.join(root, file)
                
    fn_ras = os.path.join(avaDir, 'Avaflow_Input', 'DEM.asc')            
             
    fn_vec = shpFile 
    ras_ds = gdal.Open(fn_ras) 
    vec_ds = ogr.Open(fn_vec) 
    lyr = vec_ds.GetLayer() 
    geot = ras_ds.GetGeoTransform()
    
    ##new raster
    new_out = os.path.join(avaDir, 'Avaflow_Input', 'rel_rst.tif')
    drv_tiff = gdal.GetDriverByName("GTiff") 
    chn_ras_ds = drv_tiff.Create(new_out, ras_ds.RasterXSize, ras_ds.RasterYSize, 1, gdal.GDT_Float32)
    chn_ras_ds.SetGeoTransform(geot)
    init_value = 0
    chn_ras_ds.GetRasterBand(1).Fill(init_value)
    
    
    gdal.RasterizeLayer(chn_ras_ds, [1], lyr) 
    chn_ras_ds.GetRasterBand(1).SetNoDataValue(-99) 
    chn_ras_ds = None
    
    source = ogr.Open(fn_vec)
    layer = source.GetLayer()
    schema = []
    ldefn = layer.GetLayerDefn()
    for n in range(ldefn.GetFieldCount()):
        fdefn = ldefn.GetFieldDefn(n)
        schema.append(fdefn.name)
    print(schema)
    
    if "thickness" in schema:
        #####Get thickness
        shapefile = shpFile
        driver = ogr.GetDriverByName("ESRI Shapefile")
        ds = driver.Open(shapefile, 0)
        layer = ds.GetLayer()    
        thickness_value = layer[0].GetField("thickness")
        print("This is thickness release", thickness_value)
        # Thickness value needs to transfered to the thickness value in r.avaflow
        
        ## Raster calculation, set value 255 to the thickness value
        gdal_path = '/bin'
        gdal_calc_path = os.path.join(gdal_path, 'gdal_calc.py')
        
        # Arguements
        input_file_path = os.path.join(avaDir, 'Avaflow_Input', 'rel_rst.tif')
        output_file_path = os.path.join(avaDir, 'Avaflow_Input', 'rel_rst.tif')
        calc_expr = '"{0}*(A>=255)"'
        calc_expr_final = calc_expr.format(thickness_value)
        print(calc_expr_final)
        nodata = '0'
        typeof = '"Float32"'
        
        # Generate string of process.
        gdal_calc_str = 'python {0} -A {1} --outfile={2} --calc={3} --NoDataValue={4} --type={5} --overwrite '
        gdal_calc_process = gdal_calc_str.format(gdal_calc_path, input_file_path, 
            output_file_path, calc_expr_final, nodata, typeof)
        
        # Call process.
        os.system(gdal_calc_process)
        
    else:
        ## Raster calculation, set value 255 to the thickness value
        gdal_path = '/bin'
        gdal_calc_path = os.path.join(gdal_path, 'gdal_calc.py')
        
        # Arguements
        input_file_path = os.path.join(avaDir, 'Avaflow_Input', 'rel_rst.tif')
        output_file_path = os.path.join(avaDir, 'Avaflow_Input', 'rel_rst.tif')
        calc_expr = '"1*(A>=255)"'
        nodata = '0'
        typeof = '"Float32"'
        
        # Generate string of process.
        gdal_calc_str = 'python {0} -A {1} --outfile={2} --calc={3} --NoDataValue={4} --type={5} --overwrite '
        gdal_calc_process = gdal_calc_str.format(gdal_calc_path, input_file_path, 
            output_file_path, calc_expr, nodata, typeof)
        
        # Call process.
        os.system(gdal_calc_process)
        
        
    ### .tif to ascii
    input_file_path = os.path.join(avaDir, 'Avaflow_Input', 'rel_rst.tif')
    output_file_path = os.path.join(avaDir, 'Avaflow_Input', 'rel_rst.asc')
    
    # Generate string of process.
    gdal_calc_str = 'gdal_translate -of AAIGrid {0} {1}'
    gdal_calc_process = gdal_calc_str.format(input_file_path, output_file_path)
    # Call process.
    os.system(gdal_calc_process)
    
    ###-------------------------------------------------------------------------------
    ##Same process for secondary release area
    inputDir_rel = pathlib.Path(avaDir, 'Inputs/SECREL')
    shpFile = list(inputDir_rel.glob('*.shp'))
    
    if len(shpFile) == 1:
        ## Rasterize shapefile 
        
        search_path = os.path.join(avaDir, 'Avaflow_Input')
        for root, dir, files in os.walk(search_path):
            for file in files:
                if file.startswith("secondary") and file.endswith(".shp"):
                    shpFile = os.path.join(root, file)
        
                    
        fn_ras = os.path.join(avaDir, 'Avaflow_Input', 'DEM.asc')            
                 
        fn_vec = shpFile 
        ras_ds = gdal.Open(fn_ras) 
        vec_ds = ogr.Open(fn_vec) 
        lyr = vec_ds.GetLayer() 
        geot = ras_ds.GetGeoTransform()
        
        ##new raster
        new_out = os.path.join(avaDir, 'Avaflow_Input', 'rel2_rst.tif')
        drv_tiff = gdal.GetDriverByName("GTiff") 
        chn_ras_ds = drv_tiff.Create(new_out, ras_ds.RasterXSize, ras_ds.RasterYSize, 1, gdal.GDT_Float32)
        chn_ras_ds.SetGeoTransform(geot)
        init_value = 0
        chn_ras_ds.GetRasterBand(1).Fill(init_value)
        
        
        gdal.RasterizeLayer(chn_ras_ds, [1], lyr) 
        chn_ras_ds.GetRasterBand(1).SetNoDataValue(-99) 
        chn_ras_ds = None
        
        source = ogr.Open(fn_vec)
        layer = source.GetLayer()
        schema = []
        ldefn = layer.GetLayerDefn()
        for n in range(ldefn.GetFieldCount()):
            fdefn = ldefn.GetFieldDefn(n)
            schema.append(fdefn.name)
        print(schema)
        
        if "thickness" in schema:
            #####Get thickness
            shapefile = shpFile
            driver = ogr.GetDriverByName("ESRI Shapefile")
            ds = driver.Open(shapefile, 0)
            layer = ds.GetLayer()    
            thickness_value = layer[0].GetField("thickness")
            print("This is thickness release 2", thickness_value)
            # Thickness value needs to transfered to the thickness value in r.avaflow
            
            ## Raster calculation, set value 255 to the thickness value
            gdal_path = '/bin'
            gdal_calc_path = os.path.join(gdal_path, 'gdal_calc.py')
            
            # Arguements
            input_file_path = os.path.join(avaDir, 'Avaflow_Input', 'rel2_rst.tif')
            output_file_path = os.path.join(avaDir, 'Avaflow_Input', 'rel2_rst.tif')
            calc_expr = '"{0}*(A>=255)"'
            calc_expr_final = calc_expr.format(thickness_value)
            print(calc_expr_final)
            nodata = '0'
            typeof = '"Float32"'
            
            # Generate string of process.
            gdal_calc_str = 'python {0} -A {1} --outfile={2} --calc={3} --NoDataValue={4} --type={5} --overwrite '
            gdal_calc_process = gdal_calc_str.format(gdal_calc_path, input_file_path, 
                output_file_path, calc_expr_final, nodata, typeof)
            
            # Call process.
            os.system(gdal_calc_process)
            
        else:
            ## Raster calculation, set value 255 to the thickness value
            gdal_path = '/bin'
            gdal_calc_path = os.path.join(gdal_path, 'gdal_calc.py')
            
            # Arguements
            input_file_path = os.path.join(avaDir, 'Avaflow_Input', 'rel2_rst.tif')
            output_file_path = os.path.join(avaDir, 'Avaflow_Input', 'rel2_rst.tif')
            calc_expr = '"1*(A>=255)"'
            nodata = '0'
            typeof = '"Float32"'
            
            # Generate string of process.
            gdal_calc_str = 'python {0} -A {1} --outfile={2} --calc={3} --NoDataValue={4} --type={5} --overwrite '
            gdal_calc_process = gdal_calc_str.format(gdal_calc_path, input_file_path, 
                output_file_path, calc_expr, nodata, typeof)
            
            # Call process.
            os.system(gdal_calc_process)
            
            
        ### .tif to ascii
        input_file_path = os.path.join(avaDir, 'Avaflow_Input', 'rel2_rst.tif')
        output_file_path = os.path.join(avaDir, 'Avaflow_Input', 'rel2_rst.asc')
        
        # Generate string of process.
        gdal_calc_str = 'gdal_translate -of AAIGrid {0} {1}'
        gdal_calc_process = gdal_calc_str.format(input_file_path, output_file_path)
        # Call process.
        os.system(gdal_calc_process)
    
    ###-------------------------------------------------------------------------------
    ##Same process for entrainment area
    
    ## Rasterize release area shapefile 
    inputDir_rel = pathlib.Path(avaDir, 'Inputs/ENT')
    shpFile = list(inputDir_rel.glob('*.shp'))
    
    if len(shpFile) == 1:
        search_path = os.path.join(avaDir, 'Avaflow_Input')
        for root, dir, files in os.walk(search_path):
            for file in files:
                if file.startswith("ent") and file.endswith(".shp"):
                    shpFile = os.path.join(search_path, file)
    
                    
        fn_ras = os.path.join(avaDir, 'Avaflow_Input', 'DEM.asc')            
                 
        fn_vec = shpFile 
        ras_ds = gdal.Open(fn_ras) 
        vec_ds = ogr.Open(fn_vec) 
        lyr = vec_ds.GetLayer() 
        geot = ras_ds.GetGeoTransform()
        
        ##new raster
        new_out = os.path.join(avaDir, 'Avaflow_Input', 'ent_rst.tif')
        drv_tiff = gdal.GetDriverByName("GTiff") 
        chn_ras_ds = drv_tiff.Create(new_out, ras_ds.RasterXSize, ras_ds.RasterYSize, 1, gdal.GDT_Float32)
        chn_ras_ds.SetGeoTransform(geot)
        init_value = 0
        chn_ras_ds.GetRasterBand(1).Fill(init_value)
        
        
        gdal.RasterizeLayer(chn_ras_ds, [1], lyr) 
        chn_ras_ds.GetRasterBand(1).SetNoDataValue(-99) 
        chn_ras_ds = None
        
        source = ogr.Open(fn_vec)
        layer = source.GetLayer()
        schema = []
        ldefn = layer.GetLayerDefn()
        for n in range(ldefn.GetFieldCount()):
            fdefn = ldefn.GetFieldDefn(n)
            schema.append(fdefn.name)
        print(schema)
        
        if "thickness" in schema:
            #####Get thickness
            shapefile = shpFile
            driver = ogr.GetDriverByName("ESRI Shapefile")
            ds = driver.Open(shapefile, 0)
            layer = ds.GetLayer()    
            thickness_value = layer[0].GetField("thickness")
            print("This is thickness entrainment", thickness_value)
            # Thickness value needs to transfered to the thickness value in r.avaflow
            
            ## Raster calculation, set value 255 to the thickness value
            gdal_path = '/bin'
            gdal_calc_path = os.path.join(gdal_path, 'gdal_calc.py')
            
            # Arguements
            input_file_path = os.path.join(avaDir, 'Avaflow_Input', 'ent_rst.tif')
            output_file_path = os.path.join(avaDir, 'Avaflow_Input', 'ent_rst.tif')
            calc_expr = '"{0}*(A>=255)"'
            calc_expr_final = calc_expr.format(thickness_value)
            print(calc_expr_final)
            nodata = '0'
            typeof = '"Float32"'
            
            # Generate string of process.
            gdal_calc_str = 'python {0} -A {1} --outfile={2} --calc={3} --NoDataValue={4} --type={5} --overwrite '
            gdal_calc_process = gdal_calc_str.format(gdal_calc_path, input_file_path, 
                output_file_path, calc_expr_final, nodata, typeof)
            
            # Call process.
            os.system(gdal_calc_process)
            
        else:
            ## Raster calculation, set value 255 to the thickness value
            gdal_path = '/bin'
            gdal_calc_path = os.path.join(gdal_path, 'gdal_calc.py')
            
            # Arguements
            input_file_path = os.path.join(avaDir, 'Avaflow_Input', 'ent_rst.tif')
            output_file_path = os.path.join(avaDir, 'Avaflow_Input', 'ent_rst.tif')
            calc_expr = '"1*(A>=255)"'
            nodata = '0'
            typeof = '"Float32"'
            
            # Generate string of process.
            gdal_calc_str = 'python {0} -A {1} --outfile={2} --calc={3} --NoDataValue={4} --type={5} --overwrite '
            gdal_calc_process = gdal_calc_str.format(gdal_calc_path, input_file_path, 
                output_file_path, calc_expr, nodata, typeof)
            
            # Call process.
            os.system(gdal_calc_process)
            
            
        ### .tif to ascii
        input_file_path = os.path.join(avaDir, 'Avaflow_Input', 'ent_rst.tif')
        output_file_path = os.path.join(avaDir, 'Avaflow_Input', 'ent_rst.asc')
        
        # Generate string of process.
        gdal_calc_str = 'gdal_translate -of AAIGrid {0} {1}'
        gdal_calc_process = gdal_calc_str.format(input_file_path, output_file_path)
        # Call process.
        os.system(gdal_calc_process)
    
def cleandirec(avaDir):
    
    # Delete all tif files
    inputDir = pathlib.Path(avaDir, 'Avaflow_Input')        
    for f in os.listdir(inputDir):
        if not f.endswith(".tif"):
            continue
        os.remove(os.path.join(inputDir, f))