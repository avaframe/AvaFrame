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
from avaframe.ana6Optimisation.helper import optimiseNonSeqV1
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

avaName = 'avaFleisskar'
avaDir = 'data/' + avaName


# Read and merge results from AIMEC, parametersets for simulations, areal indicators
finalDF = helper.buildFinalDF(avaName, cfgProb)


paramSelected = ['Dry-friction coefficient (-)', 'Density (kg/m^3)', 'Deposition rate 21  (m/s)', 'Basal drag coeff. 1-2 (-)',
              'Top drag coeff. (-)', 'Avalanche shear strength (Pa)', 'Turbulent drag coefficient (-)']

# ToDo Chose either all simulations or only LH simulations, currently --> all
# create df with most important parameters and the loss function
emulatorDF, emulatorScaledDF = helper.createDFParameterLoss(finalDF, varParList, paramSelected)
# remove all entries where run out is not reached bzw where loss function is 1
emulatorDF = emulatorDF[emulatorDF['Loss'] != 1]




# now bayesian optimisation
# train surrogate
X, y, gp_pipe, etr_pipe = helper.fitSurrogate(emulatorDF)

# K fold cross validation
# helper.KFoldCV(X,y, gp_pipe, "Gaussian Process Matern 5/2 Kernel")
# helper.KFoldCV(X,y, etr_pipe, "Extra Trees Surrogate")

# fit final pipline ToDo: check if fit pipe can be done before CV
gp_pipe.fit(X, y)


# optimize
# v1, non sequential (only use pipe once to find best param) ToDo create LH samples with avaframe functions
optimiseNonSeqV1(gp_pipe, paramSelected)







print(gp_pipe)
print('sss')


















