# -*- coding: utf-8 -*-
"""
Created on Mon Wed 12 2022

@author: dicko

Tools to extract information on the avalanche simulations run in the Output files  

"""

# Local imports
from avaframe.com1DFA import com1DFA
from avaframe.in3Utils import cfgUtils


# Function to extract the different avalanche simulations of the same output file 
def postProcess(avalancheDir, cfgMain, simDF):
    """ loop on all DFA simulations
    (avalancheDir, cfgMain, cfgSimi, simDF, solSimi, outDirTest)

    Parameters
    -----------
    avalancheDir: str or pathlib path
        avalanche directory
    cfgMain: confiparser
        avaframeCfg configuration
    simDF: pandas dataFrame
        configuration DF

    Returns
    --------
    simDF: pandas dataFrame
        configuration DF appended with the analysis results
    """
    # loop on all the simulations 
    for simHash, simDFrow in simDF.iterrows():
        simName = simDFrow['simName']
        # fetch the simulation results
        fieldsList, fieldHeader, timeList = com1DFA.readFields(avalancheDir, ['FT', 'FV', 'pft', 'pfv', 'ppr'],
                                                               simName=simName, flagAvaDir=True, comModule='com1DFA', atol=1.e-6)
    return simDF

    

# Function to extract the different avalanche simulations using postProcess 
def FindAvaSimu(): 
    # Load avalanche directory from general configuration file
    cfgMain = cfgUtils.getGeneralConfig()
    avalancheDir = cfgMain['MAIN']['avalancheDir']
    
    # Load configuration info of all com1DFA simulations
    simDF, _ = cfgUtils.readAllConfigurationInfo(avalancheDir) 
    
    # Extracting the different avalanche simulations of the output file 
    F = postProcess(avalancheDir, cfgMain, simDF)
    
    return F 