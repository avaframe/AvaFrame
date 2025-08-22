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
from datetime import datetime

import avaframe.out3Plot.plotUtils as pU
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


from avaframe.com1DFA import com1DFA
from avaframe.com8MoTPSA import com8MoTPSA
from avaframe.ana4Stats import probAna
import helper
import plotSAResults










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


avaName = 'avaFleisskar'
avaDir = 'data/' + avaName


# Read results from AIMEC, parametersets for simulations, areal indicators and merge
finalDF = helper.buildFinalDF(avaName, cfgProb)


# use only morris samples
morrisDF = finalDF[finalDF['sampleMethod'] == 'morris'].copy(deep=True)
# set order as index
morrisDF.set_index('order', inplace=True)
# order df based on order (which is the index)
morrisDF.sort_index(inplace=True)


# Load variation data
valVariationValue = cfgProb['PROBRUN']['variationValue'].split('|')
varType = cfgProb['PROBRUN']['variationType'].split('|')
varParList = cfgProb['PROBRUN']['varParList'].split('|')

# define input for SA
paramValuesD = probAna.createSampleWithVariationStandardParameters(cfgProb, cfgStart, varParList, valVariationValue, varType)
problem = {
    'num_vars': len(varParList),
    'names': varParList,
    'bounds': paramValuesD['bounds']
}
samples = np.vstack(morrisDF['parameterSet'].values).astype(float)
Y = morrisDF['optimisationVariable'].values

Si = morris_analyze.analyze(
    problem,
    samples,
    Y,
    conf_level=0.95,
    print_to_console=True,
    num_levels=20
)

# Rank Parameters
SiData = {
    "Parameter": Si['names'],
    "mu_star": Si['mu_star'],
    "sigma": Si['sigma'],
    "mu_star_conf": Si['mu_star_conf']}
# In DataFrame umwandeln
SiDF = pd.DataFrame(SiData)

# plot results
outDir = pathlib.Path('../' + avaDir +  "/Outputs" + "/ana6Morris")
# create folder if not already existing
fU.makeADir(outDir)
plotSAResults.barplot(SiDF, avaName, outDir)
# ToDo make scatter plot better, visibility better, and maybe merge scatterplot functions
plotSAResults.scatterplot(SiDF, avaName, outDir)
plotSAResults.scatterplotUncertainty(SiDF, avaName, outDir)


# Get mean values of best fit parameter set of morris samples and create new dataframe for LHS
# Chosen input parameters are
paramSelected = ['Dry-friction coefficient (-)', 'Density (kg/m^3)', 'Deposition rate 21  (m/s)', 'Basal drag coeff. 1-2 (-)',
              'Top drag coeff. (-)', 'Avalanche shear strength (Pa)', 'Turbulent drag coefficient (-)']


# create df with most important parameters and the loss function
paramLossDF, paramLossScaledDF = helper.createDFParameterLoss(morrisDF, varParList, paramSelected)

N = 50 # for mean, number of morris simulation to take into account for LHS start
paramLossScaledSubsetDF = paramLossScaledDF.sort_values(by='Loss', ascending=True)[:N]
paramLossSubsetDF = paramLossDF.sort_values(by='Loss', ascending=True)[:N]

#save overview as csv, need in next step (LH sampling)
# include date, format: YYYYMMDD
date = datetime.now().strftime("%Y%m%d")
csvPath = f"{outDir}/{avaName}_LHSStartValues_{date}.csv"
paramLossSubsetDF.describe().to_csv(csvPath)







































