import numpy as np
import logging
import pathlib
from scipy.stats import qmc
from SALib.sample import morris
from SALib.analyze import morris as morris_analyze
import pickle

import pandas as pd
import configparser
import os
import matplotlib.pyplot as plt
import seaborn as sns
import rasterio
from rasterio.mask import mask
import geopandas as gpd
import json
import re

import avaframe.out3Plot.plotUtils as pU
from avaframe.ana6Optimisation.helper import optimiseNonSeqV1, EINextPoint
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils
from avaframe.in3Utils import cfgHandling
from avaframe.in3Utils import fileHandlerUtils as fU
import avaframe.in2Trans.rasterUtils as IOf
import avaframe.in1Data.computeFromDistribution as cP
import avaframe.com1DFA.deriveParameterSet as dP
from avaframe.in3Utils import geoTrans as gT
from avaframe.out3Plot import statsPlots as sP
from avaframe.in1Data import getInput as gI

import avaframe.out1Peak.outPlotAllPeakDiffs as oPD
import avaframe.in2Trans.shpConversion as shpConv
import avaframe.com1DFA.DFAtools as DFAtls


from avaframe.ana4Stats import probAna
from avaframe.com1DFA import com1DFA
from avaframe.com8MoTPSA import com8MoTPSA
from avaframe.ana4Stats import probAna

import helper
import plotSAResults
from avaframe.runScripts.runPlotAreaRefDiffs import runPlotAreaRefDiffs







# Determine if config should be read from current ini files or previous ini files (stored in cfgData folder)
# -----------------------------------------------------------------------------------------

# Ini file from current model run
#cfgProb = cfgUtils.getModuleConfig(probAna)
#cfgStart = cfgUtils.getModuleConfig(com8MoTPSA)

# -----------------------------------------------------------------------------------------
# Stored ini files, from,first morris run

# Stored ini files, from,first morris run
cfgProb = configparser.ConfigParser()
cfgStart = configparser.ConfigParser()

cfgProb.read('/home/lawine_naturgefahren/Fischbacher_Roland/AvaFrame/avaframe/data/cfgData/MorrisAllParamsCom8MoTPSA/morris_probAnaCfg.ini')
cfgStart.read('/home/lawine_naturgefahren/Fischbacher_Roland/AvaFrame/avaframe/data/cfgData/MorrisAllParamsCom8MoTPSA/morris_com8MoTPSACfg.ini')
varParList = cfgProb['PROBRUN']['varParList'].split('|')

# Load avalanche directory from general configuration file
cfgMain = cfgUtils.getGeneralConfig()
avalancheDir = cfgMain['MAIN']['avalancheDir']
avaName = pathlib.Path(avalancheDir).name


paramSelected = ['Dry-friction coefficient (-)', 'Density (kg/m^3)', 'Deposition rate 21  (m/s)', 'Basal drag coeff. 1-2 (-)',
              'Top drag coeff. (-)', 'Avalanche shear strength (Pa)', 'Turbulent drag coefficient (-)']
paramBounds = {
    'Dry-friction coefficient (-)': (0.15, 0.4),
    'Density (kg/m^3)': (50, 300),
    'Deposition rate 21  (m/s)': (0.15, 0.5),
    'Basal drag coeff. 1-2 (-)': (0.02, 0.05),
    'Top drag coeff. (-)': (1e-5, 1e-3),
    'Avalanche shear strength (Pa)': (0.0, 5.0),
    'Turbulent drag coefficient (-)': (0.0015, 0.005),
}

# plot results
outDir = pathlib.Path(avalancheDir, "Outputs", "ana6Morris")
# create folder if not already existing
fU.makeADir(outDir)



# calculate areal indicators and save pickle
resType = "ppr"
thresholdValueSimulation = 1
modName = 'com8MoTPSA'
runPlotAreaRefDiffs(resType, thresholdValueSimulation, modName)

# Read and merge results from parametersets for simulations, areal indicators
finalDF = helper.buildFinalDF(avalancheDir, cfgProb, paramSelected)

optimisationType = 'nonSeq'

if optimisationType == 'nonSeq':

    # save sim with currently best y
    helper.saveBestRow(finalDF, 'optimisationVariable',
                       csv_path=(avalancheDir + "/Outputs" + "/ana6Morris" + "/bestNonSeq.csv"))

    # ToDo Chose either all simulations or only LH simulations, currently --> all
    # create df with most important parameters and the loss function
    emulatorDF, emulatorScaledDF = helper.createDFParameterLoss(finalDF, varParList, paramSelected)
    # remove all entries where run out is not reached bzw where loss function is 1
    emulatorDF = emulatorDF[emulatorDF['Loss'] != 1]

    # train surrogate
    X, y, gp_pipe, etr_pipe = helper.fitSurrogate(emulatorDF)

    # K fold cross validation
    # helper.KFoldCV(X, y, gp_pipe, "Gaussian Process Matern 5/2 Kernel")
    # helper.KFoldCV(X,y, etr_pipe, "Extra Trees Surrogate")

    # fit final pipline
    gp_pipe.fit(X, y)

    # ----------------------------------------------------------------------------------------------------------------------
    # optimize
    # v1, non sequential (only use pipe once to find best param)
    topNStat = helper.optimiseNonSeqV1(gp_pipe, paramBounds)

    # run com8 with best x from N best surrogate mean
    simName = helper.runCom8MoTPSA(avalancheDir, topNStat['TopNBest']['mean_params'], cfgMain,
                                   optimisationType='nonSeq')

    # das gleiche für Best Surrogate




    # calculate areal indicators and save pickle
    runPlotAreaRefDiffs(resType, thresholdValueSimulation, modName)
    # Read and merge results from parametersets for simulations, areal indicators
    finalDF = helper.buildFinalDF(avalancheDir, cfgProb, paramSelected)

    # save best surrogate input params
    helper.saveTopCandidates(topNStat, finalDF, paramSelected,
                             out_path=avalancheDir + "/Outputs" + "/ana6Morris" + f"/bestNonSeq{avaName}.png",
                             title=f"{avaName} NonSeq-Analysis: Best Surrogate vs Best Real Model Run",
                             simName=simName)

    # save latest real sim
    helper.saveBestRow(df=finalDF, y='optimisationVariable', simName=simName,
                       csv_path=(avalancheDir + "/Outputs" + "/ana6Morris" + "/bestNonSeq.csv"))



else:

    # save sim with currently best y
    helper.saveBestRow(finalDF, 'optimisationVariable',
                       csv_path=(avalancheDir + "/Outputs" + "/ana6Morris" + "/bestBORows.csv"))

    for i in range(40):
        # ToDo Chose either all simulations or only LH simulations, currently --> all
        # create df with most important parameters and the loss function
        emulatorDF, emulatorScaledDF = helper.createDFParameterLoss(finalDF, varParList, paramSelected)
        # remove all entries where run out is not reached bzw where loss function is 1
        emulatorDF = emulatorDF[emulatorDF['Loss'] != 1]

        # train surrogate
        X, y, gp_pipe, etr_pipe = helper.fitSurrogate(emulatorDF)
        # K fold cross validation
        # helper.KFoldCV(X, y, gp_pipe, "Gaussian Process Matern 5/2 Kernel")
        # helper.KFoldCV(X,y, etr_pipe, "Extra Trees Surrogate")

        # fit final pipline
        gp_pipe.fit(X, y)

        # ---------------------------------------------------------------------------------------------------------------------
        # v2, sequential
        # get next input parameters with EI
        xBest, xBestDict, ei, lcb = helper.EINextPoint(gp_pipe, y, paramBounds)

        # run com8 with best x
        simName = helper.runCom8MoTPSA(avalancheDir, xBestDict, cfgMain, i)

        # calculate areal indicators and save pickle
        runPlotAreaRefDiffs(resType, thresholdValueSimulation, modName)

        # Read and merge results from parametersets for simulations, areal indicators
        finalDF = helper.buildFinalDF(avalancheDir, cfgProb, paramSelected)

        # save latest sim
        helper.saveBestRow(finalDF, 'optimisationVariable', ei, lcb, simName,
                           csv_path=(avalancheDir + "/Outputs" + "/ana6Morris" + "/bestBORows.csv"))

        # if ei < 0.002:
        #    break

    #
    #
    #
    #
    #
    #
    # ToDo Define how to optimize hyperparameters









print(gp_pipe)
print('sss')


















