"""
    Generate a markdown report for data provided in dictionary with format to compare two simulations

    This file is part of Avaframe.
"""

# Load modules
import os
import glob
import logging
import numpy as np
import shutil
from avaframe.in3Utils import fileHandlerUtils as fU

# create local logger
# change log level in calling module to DEBUG to see log messages
log = logging.getLogger(__name__)


def copyPlots(avaName, outDir, plotListRep):
    """ copy the quick plots to report directory """

    plotPaths = {}
    for key in plotListRep:
        shutil.copy2(plotListRep[key], os.path.join(outDir, '%s_%s.png' % (avaName, key)))
        log.info('Copied: %s to %s' % (plotListRep[key], os.path.join(outDir, '%s_%s.png' % (avaName, key))))
        plotPaths[key] = os.path.join(outDir, '%s_%s.png' % (avaName, key))

    return plotPaths


def makeLists(simDict, benchDict):
    """ Create a list for table creation """

    parameterList = []
    valuesSim = []
    valuesBench = []
    for key in simDict:
        if key != 'type':
            if simDict[key] == '' and benchDict.get(key) is None:
                log.info('this parameter is not used: %s' % key)
            else:
                parameterList.append(key)
                valuesSim.append(simDict[key])
                if key in benchDict:
                    valuesBench.append(benchDict[key])
                else:
                    valuesBench.append('non existent')

    return parameterList, valuesSim, valuesBench


def writeCompareReport(reportFile, reportD, benchD, avaName, cfgRep):
    """ Write a report for simulation """

    # Set limit to produce warning for differences in peak fields
    diffLim = float(cfgRep['GENERAL']['diffLim'])

    # Start writing markdown style report
    with open(reportFile, 'a') as pfile:

        # HEADER BLOCK
        # Simulation name
        pfile.write('## *%s* \n' % avaName)
        simName = reportD['simName']['name']
        pfile.write('### Simulation name: *%s* \n' % reportD['simName']['name'])
        if benchD['simName'] != simName:
            textString = '<span style="color:red"> Reference simulation name is different: \
                         %s  </span>' % benchD['simName']
            pfile.write('#### %s \n' % textString)
        pfile.write(' \n')

        # TEST Info
        pfile.write('#### Test Info \n')
        pfile.write(' \n')
        for value in benchD['Test Info']:
            if value != 'type':
                # pfile.write('##### Topic:  %s \n' % value)
                pfile.write('%s \n' % (benchD['Test Info'][value]))

        # Create lists to write tables
        parameterList, valuesSim, valuesBench = makeLists(reportD['Simulation Parameters'], benchD['Simulation Parameters'])

        # PARAMETER BLOCK
        # Table listing all the simulation parameters
        pfile.write('#### Simulation Parameters \n')
        pfile.write('| Parameter | Reference | Simulation | \n')
        pfile.write('| --------- | --------- | ---------- | \n')
        countValue = 0
        for value in parameterList:
            # if reference and simulation parameter value do not match
            # mark red
            if valuesBench[countValue] != valuesSim[countValue]:
                valueRed = '<span style="color:red"> %s </span>' % valuesSim[countValue]
                pfile.write('| %s | %s | %s | \n' % (value, valuesBench[countValue], valueRed))
                countValue = countValue + 1
            else:
                pfile.write('| %s | %s | %s | \n' % (value, valuesBench[countValue], valuesSim[countValue]))
                countValue = countValue + 1

        # Add time needed for simulation to table
        pfile.write('| Run time [s] |  | %.2f | \n' % (reportD['runTime']))
        pfile.write(' \n')
        pfile.write(' \n')

        # IMAGE BLOCK
        pfile.write('#### Comparison Plots \n')
        for value in reportD['Simulation Results']:
            if value != 'type':
                pfile.write('##### Figure:   %s \n' % value)
                minDiff = float(reportD['Simulation Difference'][value][2])
                maxVal = float(reportD['Simulation Stats'][value][0])
                if minDiff < (-1.0*diffLim*maxVal) or minDiff > (diffLim*maxVal):
                    textString = '<span style="color:red"> Warning difference in %s is Max: %0.2f, \
                                 Mean %.02f and Min %.02f </span>' % (value, reportD['Simulation Difference'][value][0],
                                 reportD['Simulation Difference'][value][1], reportD['Simulation Difference'][value][2])
                    pfile.write(' %s \n' % textString)
                pfile.write('![%s](%s) \n' % (value, reportD['Simulation Results'][value]))

        pfile.write(' \n')
        pfile.write('------------------------------------------------------------------------- \n')
        pfile.write(' \n')
        pfile.write(' \n')

    log.info('standard tests performed - Report saved to: %s ' % reportFile)
