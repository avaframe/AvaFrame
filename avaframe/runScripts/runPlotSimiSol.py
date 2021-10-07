"""
    Run com1DFA kernel and compare it tothe similarity solution
    This script computes the similarity solution for a gliding avalanche on
    a inclined plane according to similarity solution from :
    Hutter, K., Siegel, M., Savage, S.B. et al.
    Two-dimensional spreading of a granular avalanche down an inclined plane
    Part I. theory. Acta Mechanica 100, 37–68 (1993).
    https://doi.org/10.1007/BF01176861
    and compares it to the DFA kernel com1DFA
"""

import pathlib
import pickle

# Local imports
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils
import avaframe.out3Plot.outAna1Plots as outAna1Plots


# +++++++++SETUP CONFIGURATION++++++++++++++++++++++++
# log file name; leave empty to use default runLog.log
logName = 'runSimilarityTest'

# Load general configuration
cfgMain = cfgUtils.getGeneralConfig()
avalancheDir = 'data/avaSimilaritySol'

# Start logging
log = logUtils.initiateLogger(avalancheDir, logName)

# create output directory for test result plots
outDirTest = pathlib.Path(avalancheDir, 'Outputs', 'ana1Tests')
with open(outDirTest / 'results10.p', 'rb') as file:
    simDF10 = pickle.load(file)
with open(outDirTest / 'results5.p', 'rb') as file:
    simDF5 = pickle.load(file)
with open(outDirTest / 'results2.p', 'rb') as file:
    simDF2 = pickle.load(file)
simDF10 = simDF10.append(simDF5)
simDF10 = simDF10.append(simDF2)
print(simDF10[['Npart', 'hErrorL2', 'dt', 'sphKernelRadius']])
outAna1Plots.plotErrorLog(simDF10, outDirTest)
