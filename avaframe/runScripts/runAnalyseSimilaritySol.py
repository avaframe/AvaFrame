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
import pandas as pd

# Local imports
import avaframe.in3Utils.initializeProject as initProj
import avaframe.in3Utils.fileHandlerUtils as fU
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils
import avaframe.com1DFA.com1DFA as com1DFA
from avaframe.in1Data import getInput as gI
import avaframe.ana1Tests.simiSolTest as simiSolTest
import avaframe.out3Plot.outAna1Plots as outAna1Plots


# +++++++++SETUP CONFIGURATION++++++++++++++++++++++++
# log file name; leave empty to use default runLog.log
logName = 'runSimilarityTest'

# Load general configuration
cfgMain = cfgUtils.getGeneralConfig()
avalancheDir = 'data/avaSimilaritySol'

# Clean input directory(ies) of old work and output files
initProj.cleanSingleAvaDir(avalancheDir, keep=logName, deleteOutput=False)

# Start logging
log = logUtils.initiateLogger(avalancheDir, logName)
log.info('MAIN SCRIPT')
log.info('Current avalanche: %s', avalancheDir)

# Load configuration for similarity solution test
simiSolCfg = pathlib.Path(avalancheDir, 'Inputs', 'simiSol_com1DFACfg.ini')

# create output directory for test result plots
outDirTest = pathlib.Path(avalancheDir, 'Outputs', 'ana1Tests')
fU.makeADir(outDirTest)

cfg = cfgUtils.getModuleConfig(com1DFA, simiSolCfg)

# Define release thickness distribution
demFile = gI.getDEMPath(avalancheDir)
simDF, _ = cfgUtils.readAllConfigurationInfo(avalancheDir)

# compute the similartiy solution (this corresponds to our reference)
log.info('Computing similarity solution')
solSimi = simiSolTest.mainSimilaritySol(simiSolCfg)

# if the analysis already exists and you only want to replot uncoment this
# pathToResults = pathlib.Path(avalancheDir, 'Outputs', 'ana1Tests', 'results10.p')
# if pathToResults.is_file():
#     simDF = pd.read_pickle(pathToResults)
simDF = simiSolTest.postProcessSimiSol(avalancheDir, cfgMain, cfg['SIMISOL'], simDF, solSimi, outDirTest)

# select the simulations you want to plot
simDF = simDF[simDF['subgridMixingFactor'].isin([10])]
simDF = simDF[simDF['cMax'].isin([0.01])]
# simDF = simDF[simDF['aPPK']==-2]
simDF = simDF[simDF['sphKernelRadius']<=10]
# simDF = simDF[simDF['nPPK']>10]#, 0.1, 0.05
# simDF['refN'] = round(simDF['cPPK'] * 5 **(simDF['aPPK']))
# simDF = simDF[simDF['nPPK0']==15]
# now do some plotting
# compare the simulations to the reference
# outAna1Plots.plotErrorRef(simDF, outDirTest, cfg['SIMISOL'], 'subgridMixingFactor', ['hErrorL2', 'vhErrorL2'],
#                           'deltaTh', 'dt', logScale=False)

# make convergence plot
fig1, ax1, ax2, slopeU, slopeH = outAna1Plots.plotErrorConvergence(simDF, outDirTest, cfg['SIMISOL'], 'nPart', ['hErrorL2', 'vhErrorL2'],
                          'aPPK', 'nPPK0', logScale=True, fit=True)

outAna1Plots.plotTimeCPULog(simDF, outDirTest, cfg['SIMISOL'], 'nPart', 'aPPK', 'nPPK0')
