"""
    Run script for running Standalone DFA
    This file is part of Avaframe.
"""

# Load modules
import os

# Local imports
from avaframe.com1DFA import com1DFA
from avaframe.log2Report import generateReport as gR
from avaframe.in3Utils import fileHandlerUtils as fU
from avaframe.out3SimpPlot import outReportPlot as oP
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils
import time

# log file name; leave empty to use default runLog.log
logName = 'runCom1DFA'

# Load avalanche directory from general configuration file
cfgMain = cfgUtils.getGeneralConfig()
avalancheDir = cfgMain['MAIN']['avalancheDir']

# Start logging
log = logUtils.initiateLogger(avalancheDir, logName)
log.info('MAIN SCRIPT')
log.info('Current avalanche: %s', avalancheDir)

# Load input parameters from configuration file
cfg = cfgUtils.getModuleConfig(com1DFA)

startTime = time.time()
# Run Standalone DFA
reportDictList = com1DFA.runSamos(cfg, avalancheDir)

# Print time needed
endTime = time.time()
log.info(('Took %s seconds to calculate.' % (endTime - startTime)))

# Generate markdown reports for each simulation
for simDict in reportDictList:

    # add parameters collected from logFile
    simDict['simParameters'].update(fU.extractParameterInfo(avalancheDir, simDict['simName']))

    # add genereal simulation parameters from .ini file
    simDict['releaseArea'].update([('release densitiy [kgm-3]', float(cfg['REP']['rhoRelease']))])
    simDict['entrainmentArea'].update([('entrainment density [kgm-3]', float(cfg['REP']['rhoEntrainment'])),
                               ('entrainment thickness [m]', float(cfg['REP']['entH']))])

    # set which parameters shall be plotted
    resultParams = cfg['REP']['plotFields']
    resPar = resultParams.split('_')
    simDict['images'] = {}
    for var in resPar:
        fileName = oP.plotPeakField(avalancheDir, simDict['simName'], var)
        # example to include figures
        imagePath = os.path.join(os.getcwd(), fileName)
        simDict['images'].update({'peak field: %s' % var : imagePath})

    # example to add addtional information
    # simDict['text'] = {'Images' : 'the images have been produced with', 'Simulations' : 'the simulations have been produced with'}

    # write report
    gR.writeReport(avalancheDir, simDict)
