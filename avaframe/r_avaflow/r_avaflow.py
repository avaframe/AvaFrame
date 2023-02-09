"""
    Main functions to execute the r.avaflow simulation
"""

import logging
import subprocess
import os
import shutil
import glob
import pathlib

# Local imports
from avaframe.in3Utils import cfgUtils
import avaframe.in3Utils.initialiseDirs as inDirs
from avaframe.r_avaflow import generate_paramfile
from avaframe.in1Data.getInput import getInputDataravaflow as getData
from avaframe.in1Data.getInput import getDEMPath
from avaframe.in1Data.getInput import readDEM
import avaframe.in2Trans.shp_to_ascii as shp_to_ascii
from avaframe.out1Peak import outPlotAllPeak as oP
from avaframe.com1DFA.com1DFA import com1DFAMain as CM
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
    

    cfgAvaflow = cfgUtils.getModuleConfig_avaflow(modName, fileOverride=cfgFile, toPrint=False)
    cfgCom1DFA = cfgUtils.getModuleConfig_avaflow2(modNameCom1DFA, fileOverride=cfgFile, toPrint=False)
    ## Create output and work directories
    workDir, outDir = inDirs.initialiseRunDirs(avalancheDir, modName, cleanDEMremeshed=False)
    
    demFile = getData(avalancheDir, cfgFile)

    dem = IOf.readRaster(demFile[0])
    
    
    #Delete unnecessary input data
    shp_to_ascii.cleandirec(avalancheDir)
    
    ##start r.avaflow
    
    ##generate the parameter file for r.avaflow out of the input data given
    
    #Check if there is already a parameter file
    param_path = os.path.join(avalancheDir, "Avaflow_Input", "param1.txt")
    if os.path.exists(param_path) == True:
        print("Old parameter file is going to be used")
        
        ##starting r.avaflow.main.c through linux executable + folder asc.-rasters as input for executlabe 
        #subprocess.call(['./r_avaflow/r_avaflow_mainc', f'./{avalancheDir}/Avaflow_Input'])
        subprocess.call(['./r_avaflow/r.avaflow.main.exe', f'./{avalancheDir}/Avaflow_Input'])
        
    else:
        generate_paramfile.generate_param(cfgAvaflow,cfgCom1DFA, avalancheDir)
    
        ##starting r.avaflow.main.c through linux executable + folder asc.-rasters as input for executlabe 
        #subprocess.call(['./r_avaflow/r_avaflow_mainc', f'./{avalancheDir}/Avaflow_Input'])
        subprocess.call(['./r_avaflow/r.avaflow.main.exe', f'./{avalancheDir}/Avaflow_Input'])
    
    #Prepare results for Ouptut folder and AIMEC (only if one phase is calculated)
    cfgGen = cfgAvaflow["GENERAL"]
    phases = cfgGen["phases"]
    phases = list(phases.split(","))
    #if len(phases) == 1 and phases[0] != "m" and phases[0] != "x":
    if len(phases) == 1 and phases[0] != "m":
        process_forAIMEC(avalancheDir,cfgAvaflow)
        
    else:
        
        workdir_orig = os.getcwd()
        os.chdir(avalancheDir)
        #Move param1.txt to output folder
        results_glob = glob.glob("Avaflow_Input/param1.txt")
        src_dir = os.path.join(workdir_orig, avalancheDir ,results_glob[0])
        dest_dir = os.path.join(workdir_orig ,avalancheDir, 'Outputs', "r_avaflow")
        shutil.move(src_dir, dest_dir)
        
        #Transfering the output of r.avaflow to the main designated output folder
        results_glob = glob.glob("Avaflow_Input/*results")
        src_dir = os.path.join(workdir_orig, avalancheDir ,results_glob[0])
        dest_dir = os.path.join(workdir_orig ,avalancheDir, 'Outputs', 'r_avaflow')
        shutil.move(src_dir, dest_dir)
    
    #Plot peak files to reports-directory
    #Check if com1DFA results are available as they are needed as reference
    folder_path = os.path.join(avalancheDir,"Outputs","com1DFA")
    if not os.path.exists(folder_path) == True:
        workdir_orig = os.getcwd()
        oP.plotAllPeakFields_ravaflow(workdir_orig, cfgAVA['FLAGS'], modName, demData=dem)

    #Prepare mass file
    #prepareMassfile(avalancheDir,cfgAvaflow)
    
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
    #Transfer vflow_max 
    results_glob = glob.glob("Avaflow_Input/*results/*ascii/*vflow_max.asc")
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
        
    #Transfer pflow_max
    results_glob = glob.glob("Avaflow_Input/*results/*ascii/*pflow_max.asc")
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
    results_glob = glob.glob("Avaflow_Input/*results/*ascii/*tflow_max.asc")
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
    results_glob = glob.glob("Avaflow_Input/*results/*ascii/*hflow_max.asc")
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
    results_glob = glob.glob("Avaflow_Input/param1.txt")
    src_dir = os.path.join(workdir_orig, avalancheDir ,results_glob[0])
    dest_dir = os.path.join(workdir_orig ,avalancheDir, 'Outputs', "r_avaflow")
    shutil.move(src_dir, dest_dir)
    
    #Transfering the output of r.avaflow to the main designated output folder
    results_glob = glob.glob("Avaflow_Input/*results")
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
    
#def writeReport(avalancheDir,cfgAvaflow,modName,dem):
    #reportDir = pathlib.Path(avalancheDir, 'Outputs', 'r_avaflow', 'reports')
    # Generate plots for all peakFiles
    #oP.plotAllPeakFields(avalancheDir, cfgAVA['FLAGS'], modName, demData=dem)
    

# def prepareMassfile(avalancheDir,cfgAvaflow):
    
#     cfgGen = cfgAvaflow["GENERAL"]
    
#     workdir_orig = os.getcwd()
#     workdir2 = os.path.join(workdir_orig, 'Outputs', "r_avaflow_compare/")    
#     os.chdir(workdir2)
#     file = glob.glob(workdir2 + "*.txt")
    
    #fp = os.path.join(workdir2, file[0])
    
    # with open(file[0]) as f:
    #     lines = f.readlines()
    
    # os.chdir(avalancheDir)
