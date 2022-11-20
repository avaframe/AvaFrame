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
import avaframe.in2Trans.shp_to_ascii as shp_to_ascii

# create local logger
log = logging.getLogger(__name__)
cfgAVA = cfgUtils.getGeneralConfig()
debugPlot = cfgAVA['FLAGS'].getboolean('debugPlot')


def r_avaflow_Main(avalancheDir, cfgFile=''):
    """
    Parameters
    ----------
    avalancheDir : str or pathlib Path
        path to avalanche data
    cfgFile : tr or pathlib path
        path to cfgFile to read overall configuration - optional if not provided the local or default config is used

    Executes the main file of r.avaflow (r.avaflow.main.c) through a linux executable.
    Input data for the executable in provided through a link leading to the folder with the asc.-raster files

    """
        
    modName = 'r_avaflow'
    modNameCom1DFA = 'com1DFA'
    

    cfgAvaflow = cfgUtils.getModuleConfig_avaflow(modName, fileOverride=cfgFile, toPrint=False)
    cfgCom1DFA = cfgUtils.getModuleConfig_avaflow2(modNameCom1DFA, fileOverride=cfgFile, toPrint=False)
    ## Create output and work directories
    workDir, outDir = inDirs.initialiseRunDirs(avalancheDir, modName,   cleanDEMremeshed=False)
    
    getData(avalancheDir, cfgFile)
    
    #Delete unnecessary input data
    shp_to_ascii.cleandirec(avalancheDir)
    
    ##start r.avaflow
    
    ##generate the parameter file for r.avaflow out of the input data given
    
    #Check if there is already a parameter file
    param_path = os.path.join(avalancheDir, "Avaflow_Input", "param1.txt")
    if os.path.exists(param_path) == True:
        print("Old parameter file is going to be used")
        
        ##starting r.avaflow.main.c through linux executable + folder asc.-rasters as input for executlabe 
        subprocess.call(['./r_avaflow/r_avaflow_mainc', f'./{avalancheDir}/Avaflow_Input'])
        
    else:
        generate_paramfile.generate_param(cfgAvaflow,cfgCom1DFA, avalancheDir)
    
        ##starting r.avaflow.main.c through linux executable + folder asc.-rasters as input for executlabe 
        subprocess.call(['./r_avaflow/r_avaflow_mainc', f'./{avalancheDir}/Avaflow_Input'])
    
    # ---- This works but not sure if really necessary
    #Transfering the output of r.avaflow to the main designated output folder
    workdir_orig = os.getcwd()
    os.chdir(avalancheDir)
    results_glob = glob.glob("Avaflow_Input/*results")
    src_dir = os.path.join(workdir_orig, avalancheDir ,results_glob[0])
    dest_dir = os.path.join(workdir_orig ,avalancheDir, 'Outputs', 'r_avaflow')
    shutil.move(src_dir, dest_dir)
    
