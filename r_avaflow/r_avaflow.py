"""
    Main functions to execute the r.avaflow simulation
"""

import logging
import subprocess
import os
import shutil
import glob

# Local imports
from avaframe.in3Utils import cfgUtils
import avaframe.in3Utils.initialiseDirs as inDirs
from avaframe.r_avaflow import generate_paramfile
from avaframe.in1Data.getInput import getInputDataravaflow as getData
import avaframe.in2Trans.shp_to_ascii as shp_to_ascii
from avaframe.out1Peak import outPlotAllPeak as oP
import avaframe.in2Trans.ascUtils as IOf

# create local logger
log = logging.getLogger(__name__)
cfgAVA = cfgUtils.getGeneralConfig()
debugPlot = cfgAVA['FLAGS'].getboolean('debugPlot')


def r_avaflow_Main(avalancheDir, cfgFile=''):
    """
    Executes the main file of r.avaflow (r.avaflow.main.c) through a linux executable.
    Input data for the executable is provided through a link leading to the folder with the asc.-raster files
    Afterwards small preparations are made for a possible AIMEC comparision run 
    
    Parameters
    ----------
    avalancheDir : str or pathlib Path
        path to avalanche data
    cfgFile : str or pathlib path
        path to cfgFile to read overall configuration - optional if not provided the local or default config is used

    """
        
    modName = 'r_avaflow'
    modNameCom1DFA = 'com1DFA'
    

    cfgAvaflow = cfgUtils.getModuleConfig_ravaflow(modName, fileOverride=cfgFile, toPrint=False)
    cfgCom1DFA = cfgUtils.getModuleConfig_ravaflow2(modNameCom1DFA, fileOverride=cfgFile, toPrint=False)
    ## Create output and work directories
    workDir, outDir = inDirs.initialiseRunDirs(avalancheDir, modName, cleanDEMremeshed=False)
    
    input_file_path_asc = os.path.join(avalancheDir, 'Ravaflow_Input', 'DEM.asc')

    if os.path.exists(input_file_path_asc):
        dem = input_file_path_asc
    else: 
        demFile = getData(avalancheDir, cfgFile)
    
        dem = IOf.readRaster(demFile[0])
         
    
    #Delete unnecessary input data
    shp_to_ascii.cleandirec(avalancheDir)
    
    ##start r.avaflow
    ##generate the parameter file for r.avaflow out of the input data given
    
    #Check if there is already a parameter file
    param_path = os.path.join(avalancheDir, "Ravaflow_Input", "param1.txt")
    if os.path.exists(param_path) == True:
        print("Old parameter file is going to be used")
        
        ##starting r.avaflow.main.c through linux executable + folder asc.-rasters as input for executlabe 
        subprocess.call(['./r_avaflow/ravaflow_main_160323', f'./{avalancheDir}/Ravaflow_Input'])
        
    else:
        generate_paramfile.generate_param(cfgAvaflow,cfgCom1DFA, avalancheDir)
    
        ##starting r.avaflow.main.c through linux executable + folder asc.-rasters as input for executlabe 
        subprocess.call(['./r_avaflow/ravaflow_main_160323', f'./{avalancheDir}/Ravaflow_Input'])
    
    ##Transfering local_r_avaflow.ini to the results  
    workdir_orig = os.getcwd()
    ravaflowIni_path = os.path.join(workdir_orig ,'r_avaflow','local_r_avaflowCfg.ini')
    dest_dir = os.path.join(workdir_orig ,avalancheDir, 'Outputs', 'r_avaflow')
    shutil.copy(ravaflowIni_path, dest_dir)
    
    #Prepare results for Ouptut folder and AIMEC (only if one phase is calculated)
    cfgGen = cfgAvaflow["GENERAL"]
    phases = cfgGen["phases"]
    phases = list(phases.split(","))
    process_forAIMEC(avalancheDir,cfgAvaflow)
            
    #Plot peak files to reports-directory
    #Check if com1DFA results are available as they are needed as reference
    workdir_orig = os.getcwd()
    folder_path = os.path.join(workdir_orig,"Outputs","com1DFA")
    if os.path.exists(folder_path):   
        oP.plotAllPeakFields_ravaflow(workdir_orig, cfgAVA['FLAGS'], modName, demData=dem)

    #Delete Ravaflow_Input-folder
    folder_path = os.path.join(workdir_orig,"Ravaflow_Input")
    # Iterate over all files in the folder and delete them
    for file_name in os.listdir(folder_path):
        file_path = os.path.join(folder_path, file_name)
        os.remove(file_path)
        

def process_forAIMEC(avalancheDir,cfgAvaflow):
    """
    Result values of r.avaflow are being renamed and stored seperately to then be used 
    for the AIMEC comparision

    Parameters
    ----------
    avalancheDir : str or pathlib Path
        path to avalanche data
    cfgAvaflow : configparser object
        full configuration settings for r.avaflow

    """
    cfgGen = cfgAvaflow["GENERAL"]
    prefix = cfgGen["prefix"]
    phases = cfgGen["phases"]
    cellsize = cfgGen["cellsize"]
    outname_basis = "ravaflow_" + prefix + phases + "_" + "cell" + cellsize
    
    #Make new folder for r.avaflow comparision
    workdir_orig = os.getcwd()
    os.chdir(avalancheDir)
    new_folder_path = os.path.join(workdir_orig ,avalancheDir, 'Outputs', "r_avaflow_compare")
    if not os.path.exists(new_folder_path):
        os.makedirs(new_folder_path)
        
    peakFiles_path = os.path.join(workdir_orig ,avalancheDir, 'Outputs', "r_avaflow_compare", "peakFiles")
    if not os.path.exists(peakFiles_path):
        os.makedirs(peakFiles_path)
    

    #Move max_files from r.avaflow results and rename to avaframe peak files
    phases = list(map(str, phases.split(",")))
    if len(phases) == 3 or phases == list("m"):
        
        #Transfer vflow_max and combine outputfiles for each phase
        results_glob1 = glob.glob("Ravaflow_Input/*results/*ascii/*vflow1_max.asc")
        src_dir1 = os.path.join(workdir_orig, avalancheDir ,results_glob1[0])
        results_glob2 = glob.glob("Ravaflow_Input/*results/*ascii/*vflow2_max.asc")
        src_dir2 = os.path.join(workdir_orig, avalancheDir ,results_glob2[0])
        results_glob3 = glob.glob("Ravaflow_Input/*results/*ascii/*vflow3_max.asc")
        src_dir3 = os.path.join(workdir_orig, avalancheDir ,results_glob3[0])
        
        gdal_path = '/bin'
        
        # ascii to tiff
        input_file_path = src_dir1
        output_file_path1 = os.path.join(workdir_orig, avalancheDir, "Outputs","r_avaflow_compare", "peakFiles", "vflow1_max.tif")
        # Generate string of process.
        gdal_calc_str = 'gdal_translate -of GTIFF {0} {1}'
        gdal_calc_process = gdal_calc_str.format(input_file_path, output_file_path1)
        # Call process.
        os.system(gdal_calc_process)
        
        # ascii to tiff
        input_file_path = src_dir2
        output_file_path2 = os.path.join(workdir_orig, avalancheDir, "Outputs","r_avaflow_compare", "peakFiles", "vflow2_max.tif")
        # Generate string of process.
        gdal_calc_str = 'gdal_translate -of GTIFF {0} {1}'
        gdal_calc_process = gdal_calc_str.format(input_file_path, output_file_path2)
        # Call process.
        os.system(gdal_calc_process)
        
        # ascii to tiff
        input_file_path = src_dir3
        output_file_path3 = os.path.join(workdir_orig, avalancheDir, "Outputs","r_avaflow_compare", "peakFiles", "vflow3_max.tif")
        # Generate string of process.
        gdal_calc_str = 'gdal_translate -of GTIFF {0} {1}'
        gdal_calc_process = gdal_calc_str.format(input_file_path, output_file_path3)
        # Call process.
        os.system(gdal_calc_process)
        
        # Raster calculator -> if overlap -> choose the raster with the highest value
        gdal_calc_path = os.path.join(gdal_path, 'gdal_calc.py')

        # Arguements
        output_file_path = os.path.join(workdir_orig ,avalancheDir, 'Outputs', "r_avaflow_compare", "peakFiles", "vflow_max.tif")
        calc_expr = '"maximum(A,B,C)"'
        nodata = '0'
        typeof = '"Float32"'

        # Generate string of process.
        gdal_calc_str = 'python {0} -A {1} -B {2} -C {3} --outfile={4} --calc={5} --NoDataValue={6} --type={7} --overwrite '
        gdal_calc_process = gdal_calc_str.format(gdal_calc_path, output_file_path1, output_file_path2, output_file_path3,
            output_file_path, calc_expr, nodata, typeof)

        # Call process.
        os.system(gdal_calc_process)
                
        ##tif to ascii
        input_file_path = os.path.join(workdir_orig ,avalancheDir, 'Outputs', "r_avaflow_compare", "peakFiles", "vflow_max.tif")
        output_file_path = os.path.join(workdir_orig, avalancheDir, "Outputs","r_avaflow_compare", "peakFiles", "vflow_max.asc")
        # Generate string of process.
        gdal_calc_str = 'gdal_translate -of AAIGrid {0} {1}'
        gdal_calc_process = gdal_calc_str.format(input_file_path, output_file_path)
        # Call process.
        os.system(gdal_calc_process)
        
        # # #Renaming
        dest_dir = os.path.join(workdir_orig ,avalancheDir, 'Outputs', "r_avaflow_compare", "peakFiles")
        old_name = glob.glob(dest_dir + "/vflow_max.asc" )
        old_name = str(old_name[0])
        outname = outname_basis + "_pfv.asc"
        new_name = os.path.join(dest_dir, outname)
        os.renames(old_name, new_name)
       
    else:  # one-phase
        results_glob = glob.glob("Ravaflow_Input/*results/*ascii/*vflow_max.asc")
        src_dir = os.path.join(workdir_orig, avalancheDir ,results_glob[0])
        dest_dir = os.path.join(workdir_orig ,avalancheDir, 'Outputs', "r_avaflow_compare", "peakFiles")
        shutil.copy(src_dir, dest_dir)
        # print("This is dest_dir:", dest_dir)
        
        #Renaming
        old_name = glob.glob(dest_dir + "/*vflow_max.asc" )
        old_name = str(old_name[0])
        outname = outname_basis + "_pfv.asc"
        new_name = os.path.join(dest_dir, outname)
        os.renames(old_name, new_name)
    

    # #Transfer pflow_max
    results_glob = glob.glob("Ravaflow_Input/*results/*ascii/*pflow_max.asc")
    src_dir = os.path.join(workdir_orig, avalancheDir ,results_glob[0])
    dest_dir = os.path.join(workdir_orig ,avalancheDir, 'Outputs', "r_avaflow_compare", "peakFiles")
    shutil.copy(src_dir, dest_dir)
    #Renaming
    old_name = glob.glob(dest_dir + "/*pflow_max.asc" )
    old_name = str(old_name[0])
    outname = outname_basis + "_ppr.asc"
    new_name = os.path.join(dest_dir, outname)
    os.renames(old_name, new_name)
    
    #Transfer tflow_max
    results_glob = glob.glob("Ravaflow_Input/*results/*ascii/*tflow_max.asc")
    src_dir = os.path.join(workdir_orig, avalancheDir ,results_glob[0])
    dest_dir = os.path.join(workdir_orig ,avalancheDir, 'Outputs', "r_avaflow_compare", "peakFiles")
    shutil.copy(src_dir, dest_dir)
    #Renaming
    old_name = glob.glob(dest_dir + "/*tflow_max.asc" )
    old_name = str(old_name[0])
    outname = outname_basis + "_pke.asc"
    new_name = os.path.join(dest_dir, outname)
    os.renames(old_name, new_name)
    
    #Transfer hflow_max
    results_glob = glob.glob("Ravaflow_Input/*results/*ascii/*hflow_max.asc")
    src_dir = os.path.join(workdir_orig, avalancheDir ,results_glob[0])
    dest_dir = os.path.join(workdir_orig ,avalancheDir, 'Outputs', "r_avaflow_compare", "peakFiles")
    shutil.copy(src_dir, dest_dir)
    #Renaming
    old_name = glob.glob(dest_dir + "/*hflow_max.asc" )
    old_name = str(old_name[0])
    outname = outname_basis + "_pft.asc"
    new_name = os.path.join(dest_dir, outname)
    os.renames(old_name, new_name)
    
    #Move param1.txt to output folder
    results_glob = glob.glob("Ravaflow_Input/param1.txt")
    src_dir = os.path.join(workdir_orig, avalancheDir ,results_glob[0])
    dest_dir = os.path.join(workdir_orig ,avalancheDir, 'Outputs', "r_avaflow")
    shutil.move(src_dir, dest_dir)
    
    #Transfering the output of r.avaflow to the main designated output folder
    results_glob = glob.glob("Ravaflow_Input/*results")
    src_dir = os.path.join(workdir_orig, avalancheDir ,results_glob[0])
    dest_dir = os.path.join(workdir_orig ,avalancheDir, 'Outputs', 'r_avaflow')
    shutil.move(src_dir, dest_dir)
    
    #---------------------------
    #Raster calculation - Correct unit from J to kJ from ravaflow_pke   
    ##ascii to tiff
    dest_dir = os.path.join(workdir_orig ,avalancheDir, 'Outputs', "r_avaflow_compare", "peakFiles")
    outname = outname_basis + "_pke.asc"
    new_name = os.path.join(dest_dir, outname)
    input_file_path = new_name
    output_file_path = os.path.join(dest_dir, "J_pke.tif")
    # Generate string of process.
    gdal_calc_str = 'gdal_translate -of GTIFF {0} {1}'
    gdal_calc_process = gdal_calc_str.format(input_file_path, output_file_path)
    # Call process.
    os.system(gdal_calc_process)
        
    ## Raster calc
    gdal_path = '/bin'
    gdal_calc_path = os.path.join(gdal_path, 'gdal_calc.py')
    
    # Arguements
    input_file_path = os.path.join(dest_dir, "J_pke.tif")
    output_file_path = os.path.join(dest_dir, 'kJ_pke.tif')
    calc_expr = '"1*(A*0.001)"'
    nodata = '0'
    typeof = '"Float32"'
    
    # Generate string of process.
    gdal_calc_str = 'python {0} -A {1} --outfile={2} --calc={3} --NoDataValue={4} --type={5} --overwrite '
    gdal_calc_process = gdal_calc_str.format(gdal_calc_path, input_file_path, 
        output_file_path, calc_expr, nodata, typeof)
    
    # Call process.
    os.system(gdal_calc_process)
    
    ##tiff to ascii
    input_file_path = os.path.join(dest_dir, "kJ_pke.tif")
    output_file_path = os.path.join(dest_dir, "kJ_pke.asc")
    # Generate string of process.
    gdal_calc_str = 'gdal_translate -of AAIGrid {0} {1}'
    gdal_calc_process = gdal_calc_str.format(input_file_path, output_file_path)
    # Call process.
    os.system(gdal_calc_process)
    
    # Delete all tif files
    inputDir = dest_dir        
    for f in os.listdir(inputDir):
        if not f.endswith(".tif") or f.endswith("_pke.asc") or f.endswith("xml") or f.endswith("prj"):
            continue
        os.remove(os.path.join(inputDir, f))
      
    #Renaming
    old_name = os.path.join(dest_dir, "kJ_pke.asc")
    # old_name = str(old_name[0])
    outname = outname_basis + "_pke.asc"
    new_name = os.path.join(dest_dir, outname)
    os.renames(old_name, new_name)
    
    ## Repeat
    #Raster calculation - Correct unit from Pa to kPa from ravaflow_ppr   
    ##ascii to tiff
    dest_dir = os.path.join(workdir_orig ,avalancheDir, 'Outputs', "r_avaflow_compare", "peakFiles")
    outname = outname_basis + "_ppr.asc"
    new_name = os.path.join(dest_dir, outname)
    input_file_path = new_name
    output_file_path = os.path.join(dest_dir, "Pa_ppr.tif")
    # Generate string of process.
    gdal_calc_str = 'gdal_translate -of GTIFF {0} {1}'
    gdal_calc_process = gdal_calc_str.format(input_file_path, output_file_path)
    # Call process.

    os.system(gdal_calc_process)
        
    ## Raster calc
    gdal_path = '/bin'
    gdal_calc_path = os.path.join(gdal_path, 'gdal_calc.py')
    
    # Arguements
    input_file_path = os.path.join(dest_dir, "Pa_ppr.tif")
    output_file_path = os.path.join(dest_dir, "kPa_ppr.tif")
    calc_expr = '"1*(A*0.001)"'
    nodata = '0'
    typeof = '"Float32"'
    
    # Generate string of process.
    gdal_calc_str = 'python {0} -A {1} --outfile={2} --calc={3} --NoDataValue={4} --type={5} --overwrite '
    gdal_calc_process = gdal_calc_str.format(gdal_calc_path, input_file_path, 
        output_file_path, calc_expr, nodata, typeof)
    
    # Call process.
    os.system(gdal_calc_process)
    
    ##tiff to ascii
    input_file_path = os.path.join(dest_dir, "kPa_ppr.tif")
    output_file_path = os.path.join(dest_dir, "kPa_ppr.asc")
    # Generate string of process.
    gdal_calc_str = 'gdal_translate -of AAIGrid {0} {1}'
    gdal_calc_process = gdal_calc_str.format(input_file_path, output_file_path)
    # Call process.
    os.system(gdal_calc_process)
    
    # Delete all tif files
    inputDir = dest_dir        
    for f in os.listdir(inputDir):
        if not f.endswith(".tif") or f.endswith("_ppr.asc") or f.endswith("xml") or f.endswith("prj"):
            continue
        os.remove(os.path.join(inputDir, f))
      
    #Renaming
    old_name = os.path.join(dest_dir, "kPa_ppr.asc")
    # old_name = str(old_name[0])
    outname = outname_basis + "_ppr.asc"
    new_name = os.path.join(dest_dir, outname)
    os.renames(old_name, new_name) 