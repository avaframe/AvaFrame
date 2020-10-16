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


def fetchBenchParameters(avaDir):
    """ Collect simulation parameter info from all standard tests """

    # get name of avalanche
    avaName = os.path.basename(avaDir)
    avaDictName = avaName + 'Dict'

    # set desired benchmark simulation info dictionary
    if avaDictName == 'avaBowlDict':
        avaDictName = {'simName': 'release1BL_null_dfa_0.155',
                    'Simulation Parameters': {
                        'type': 'list',
                        'Release Area': 'release1BL',
                        'Entrainment Area': '',
                        'Resistance Area': '',
                        'Mu': '0.155',
                        'Release thickness [m]': '1.000',
                        'Release Mass [kg]': '20196.2',
                        'Final Mass [kg]': '20196.2'},
                    'Release area': {'type': 'columns', 'Release area scenario': 'release1BL'},
                    'Test Info': {'type': 'text', 'Test Info': 'This test uses a bowl-shaped geometry.'}}

    elif avaDictName == 'avaFlatPlaneDict':
        avaDictName = {'simName': 'release1FP_null_dfa_0.155',
                    'Simulation Parameters': {
                        'type': 'list',
                        'Release Area': 'release1FP',
                        'Entrainment Area': '',
                        'Resistance Area': '',
                        'Mu': '0.155',
                        'Release thickness [m]': '1.000',
                        'Release Mass [kg]': '20000.',
                        'Final Mass [kg]': '20000.'},
                    'Release area': {'type': 'columns', 'Release area scenario': 'release1FP'},
                    'Test Info': {'type': 'text', 'Test Info': 'This test runs on a flat plane geometry.'}}

    elif avaDictName == 'avaHelixDict':
        avaDictName = {'simName': 'release1HX_null_dfa_0.155',
                    'Simulation Parameters': {
                        'type': 'list',
                        'Release Area': 'release1HX',
                        'Entrainment Area': '',
                        'Resistance Area': '',
                        'Mu': '0.155',
                        'Release thickness [m]': '1.000',
                        'Release Mass [kg]': '20522.4',
                        'Final Mass [kg]': '20522.4'},
                    'Release area': {'type': 'columns', 'Release area scenario': 'release1HX'},
                    'Test Info': {'type': 'text', 'Test Info': 'This test uses a helix-shaped geometry.'}}

    elif avaDictName == 'avaHelixChannelDict':

        avaDictName = {'simName': 'release1HX_entres_dfa_0.155',
                    'Simulation Parameters': {
                        'type': 'list',
                        'Release Area': 'release1HX',
                        'Entrainment Area': 'entrainment1HX',
                        'Resistance Area': '',
                        'Mu': '0.155',
                        'Release thickness [m]': '1.000',
                        'Release Mass [kg]': '22398.8',
                        'Final Mass [kg]': '23117.6'},
                    'Release area': {'type': 'columns', 'Release area scenario': 'release1HX'},
                    'Entrainment area': {'type': 'columns', 'Entrainment area scenario': 'entrainment1HX'},
                    'Resistance area': {'type': 'columns', 'Resistance area scenario': ''},
                    'Test Info': {'type': 'text',
                    'Test Info': 'This test uses a helix-shaped geometry with a channel with a channel.'}}

    elif avaDictName == 'avaHockeyDict':

        avaDictName = {'simName': 'release1HS_null_dfa_0.155',
                    'Simulation Parameters': {
                        'type': 'list',
                        'Release Area': 'release1HS',
                        'Mu': '0.155',
                        'Release thickness [m]': '1.000',
                        'Release Mass [kg]': '20657.1',
                        'Final Mass [kg]': '20657.1'},
                    'Release area': {'type': 'columns', 'Release area scenario': 'release1HS'},
                    'Test Info': {'type': 'text',
                    'Test Info': 'This test runs on a parabolically sloping surface with a flat foreland.'}}

    elif avaDictName == 'avaHockeySmoothChannelDict':

        avaDictName = {'simName': 'release1HS2_entres_dfa_0.155',
                    'Simulation Parameters': {
                        'type': 'list',
                        'Release Area': 'release1HS2',
                        'Entrainment Area': 'entrainment1HS2',
                        'Resistance Area': '',
                        'Mu': '0.155',
                        'Release thickness [m]': '1.000',
                        'Release Mass [kg]': '20967.3',
                        'Final Mass [kg]': '21306.'},
                    'Release area': {'type': 'columns', 'Release area scenario': 'release1HS2'},
                    'Entrainment area': {'type': 'columns', 'Entrainment area scenario': 'entrainment1HS2'},
                    'Resistance area': {'type': 'columns', 'Resistance area scenario': ''},
                    'Test Info': {'type': 'text', 'Test Info': 'This test uses a hockey stick-shaped geometry, \
                     where a linearly sloping surface transitions smoothly into a flat foreland.'}}

    elif avaDictName == 'avaHockeySmoothSmallDict':
        avaDictName = {'simName': 'release1HS2_null_dfa_0.155',
                    'Simulation Parameters': {
                        'type': 'list',
                        'Release Area': 'release1HS2',
                        'Entrainment Area': '',
                        'Resistance Area': '',
                        'Mu': '0.155',
                        'Release thickness [m]': '1.000',
                        'Release Mass [kg]': '10000.',
                        'Final Mass [kg]': '10000.'},
                    'Release area': {'type': 'columns', 'Release area scenario': 'release1HS2'},
                    'Test Info': {'type': 'text', 'Test Info': 'This test uses a hockey stick-shaped geometry, \
                     where a linearly sloping surface transitions smoothly into a flat foreland. \
                     This geometry also includes a channel.'}}

    elif avaDictName == 'avaInclinedPlaneDict':
        avaDictName = {'simName': 'release1IP_entres_dfa_0.155',
                    'Simulation Parameters': {
                        'type': 'list',
                        'Release Area': 'release1IP',
                        'Entrainment Area': 'entrainment1IP',
                        'Resistance Area': '',
                        'Mu': '0.155',
                        'Release thickness [m]': '1.000',
                        'Release Mass [kg]': '20000.',
                        'Final Mass [kg]': '21735.1'},
                    'Release area': {'type': 'columns', 'Release area scenario': 'release1IP'},
                    'Entrainment area': {'type': 'columns', 'Entrainment area scenario': 'entrainment1IP'},
                    'Resistance area': {'type': 'columns', 'Resistance area scenario': ''},
                    'Test Info': {'type': 'text', 'Test Info': 'This test runs on a linearly sloping surface.'}}

    return avaDictName


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


def writeCompareReport(reportFile, reportD, benchD, avaName):
    """ Write a report for simulation """

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

        # IMAGE BLOCK
        pfile.write('#### Comparison Plots \n')
        for value in reportD['Simulation Results']:
            if value != 'type':
                pfile.write('##### Figure:   %s \n' % value)
                pfile.write('![%s](%s) \n' % (value, reportD['Simulation Results'][value]))

        pfile.write(' \n')
        pfile.write('------------------------------------------------------------------------- \n')
        pfile.write(' \n')
        pfile.write(' \n')

    log.info('standard tests performed - Report saved to: %s ' % reportFile)
