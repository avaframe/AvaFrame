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
                        'Release Area': 'releasehh1BL',
                        'Entrainment Area': '',
                        'Resistance Area': '',
                        'Mu': '0.155',
                        'Release thickness [m]': '1.000',
                        'Release Mass [kg]': '20196.2',
                        'Final Mass [kg]': '20196.2'},
                    'Release area': {'type': 'columns', 'Release area scenario': 'release1BL'}}

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
                    'Release area': {'type': 'columns', 'Release area scenario': 'release1FP'}}

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
                    'Release area': {'type': 'columns', 'Release area scenario': 'release1HX'}}

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
                    'Resistance area': {'type': 'columns', 'Resistance area scenario': ''}}

    elif avaDictName == 'avaHockeyDict':

        avaDictName = {'simName': 'release1HS_null_dfa_0.155',
                    'Simulation Parameters': {
                        'type': 'list',
                        'Release Area': 'release1HS',
                        'Mu': '0.155',
                        'Release thickness [m]': '1.000',
                        'Release Mass [kg]': '20657.1',
                        'Final Mass [kg]': '20657.1'},
                    'Release area': {'type': 'columns', 'Release area scenario': 'release1HS'}}

    elif avaDictName == 'avaHockeySmoothChannelDict':

        avaDictName = {'simName': 'release1HS2_entres_dfa_0.155',
                    'Simulation Parameters': {
                        'type': 'list',
                        'Release Area': 'release1HX',
                        'Entrainment Area': 'entrainment1HS2',
                        'Resistance Area': '',
                        'Mu': '0.155',
                        'Release thickness [m]': '1.000',
                        'Release Mass [kg]': '20967.3',
                        'Final Mass [kg]': '21306.'},
                    'Release area': {'type': 'columns', 'Release area scenario': 'release1HS2'},
                    'Entrainment area': {'type': 'columns', 'Entrainment area scenario': 'entrainment1HS2'},
                    'Resistance area': {'type': 'columns', 'Resistance area scenario': ''}}

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
                    'Release area': {'type': 'columns', 'Release area scenario': 'release1HS2'}}

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
                    'Resistance area': {'type': 'columns', 'Resistance area scenario': ''}}

    return avaDictName


def makeLists(simDict, benchDict):
    """ Create a list for table creation """

    parameterList = []
    valuesSim = []
    valuesBench = []
    for key in simDict:
        if key != 'type':
            if simDict[key] == '' and benchDict.get(key) is None:
                log.info('this parameter is not used: %s' % simDict[key])
            else:
                parameterList.append(key)
                valuesSim.append(simDict[key])
                if key in benchDict:
                    valuesBench.append(benchDict[key])
                else:
                    valuesBench.append('non existent')

    return parameterList, valuesSim, valuesBench


def writeCompareReport(outDir, reportD, benchD):
    """ Write a report for simulation """

    # Start writing markdown style report
    with open(os.path.join(outDir, 'fullSimulationReport.md'), 'a') as pfile:

        # HEADER BLOCK
        # Simulation name
        simName = reportD['simName']['name']
        pfile.write('### Simulation name: *%s* \n' % reportD['simName']['name'])
        if benchD['simName'] != simName:
            textString = '<span style="color:red"> Reference simulation name is different: %s  </span>' % benchD['simName']
            pfile.write('#### %s \n' % textString)

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
