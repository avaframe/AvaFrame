"""
    dictionary with simulation info for benchmarks

    This file is part of Avaframe.
"""

# Load modules
import os
from avaframe.in3Utils import fileHandlerUtils as fU


def fetchBenchParameters(avaDir):
    """ Collect simulation parameter info from standard tests """

    # get name of avalanche
    avaName = os.path.basename(avaDir)
    avaDictName = avaName + 'Dict'
    avaDictList = []

    # set desired benchmark simulation info dictionary
    if avaDictName == 'avaFlatPlaneDict':
        avaDictName = {'simName': {'type': 'simName', 'name': 'release1FP_null_dfa_2.00000'},
        		'testName': 'avaFlatPlaneVarParTest',
                    'Simulation Parameters': {
                        'type': 'list',
                        'Release Area Scenario': 'release1FP',
                        'Release Area': ['Rel_Example'],
                        'Entrainment Area': '',
                        'Resistance Area': '',
                        'Parameter variation on': 'RelTh',
                        'Parameter value': '2.00000',
                        'Mu': '0.15500',
                        'Release thickness [m]': '2.00000',
                        'Release Mass [kg]': '40000.',
                        'Final Mass [kg]': '40000.',
                        'Entrainment thickness [m]': 0.3},
                    'Release area': {'type': 'columns', 'Release area scenario': 'release1FP'},
                    'Test Info': {'type': 'text', 'Test Info': 'This test runs on a flat plane geometry.'}}
        avaDictList.append(avaDictName)

    elif avaDictName == 'avaHelixChannelDict':
        avaDictName = {'simName': {'type': 'simName', 'name': 'release1HX_entres_dfa_0.05500'},
        		'testName': 'avaHelixChannelVarParTest',
                    'Simulation Parameters': {
                        'type': 'list',
                        'Release Area Scenario': 'release1HX',
                        'Release Area': ['Rel_Example', 'Rel_Two'],
                        'Entrainment Area': 'entrainment1HX',
                        'Resistance Area': '',
                        'Parameter variation on': 'Mu',
                        'Parameter value': '0.05500',
                        'Mu': '0.05500',
                        'Release thickness [m]': ['1', '1'],
                        'Entrainment thickness [m]': 0.3},
                    'Release area': {'type': 'columns', 'Release area scenario': 'release1HX'},
                    'Entrainment area': {'type': 'columns', 'Entrainment area scenario': 'entrainment1HX'},
                    'Resistance area': {'type': 'columns', 'Resistance area scenario': ''},
                    'Test Info': {'type': 'text',
                    'Test Info': 'This test uses a helix-shaped geometry with a channel with a channel.'}}
        avaDictList.append(avaDictName)


        avaDictName = {'simName': {'type': 'simName', 'name': 'release1HX_entres_dfa_0.50000'},
        		'testName': 'avaHelixChannelEnt1mVarParTest',
                    'Simulation Parameters': {
                        'type': 'list',
                        'Release Area Scenario': 'release1HX',
                        'Release Area': ['Rel_Example', 'Rel_Two'],
                        'Entrainment Area': 'entrainment1HX',
                        'Resistance Area': '',
                        'Parameter variation on': 'RelTh',
                        'Parameter value': '0.50000',
                        'Mu': '0.15500',
                        'Release thickness [m]': ['0.5', '0.5'],
                        'Entrainment thickness [m]': 1.0},
                    'Release area': {'type': 'columns', 'Release area scenario': 'release1HX'},
                    'Entrainment area': {'type': 'columns', 'Entrainment area scenario': 'entrainment1HX'},
                    'Resistance area': {'type': 'columns', 'Resistance area scenario': ''},
                    'Test Info': {'type': 'text',
                    'Test Info': 'This test uses a helix-shaped geometry with a channel with a channel.'}}

        avaDictList.append(avaDictName)

    elif avaDictName == 'avaParabolaDict':
        avaDictName = {'simName': {'type': 'simName', 'name': 'release1PF_entres_dfa_0.50000'},
        		'testName': 'avaParabolaVarParTest',
                    'Simulation Parameters': {
                        'type': 'list',
                        'Release Area Scenario': 'release1PF',
                        'Release Area': ['Rel_Example'],
                        'Entrainment Area': '',
                        'Resistance Area': 'resistance1PF',
                        'Parameter variation on': 'RelTh',
                        'Parameter value': '0.50000',
                        'Mu': '0.15500',
                        'Release thickness [m]': '0.50000',
                        'Entrainment thickness [m]': 0.3},
                    'Release area': {'type': 'columns', 'Release area scenario': 'release1PF'},
                    'Test Info': {'type': 'text',
                    'Test Info': 'This test runs on a parabolically sloping surface with a flat foreland.'}}

        avaDictList.append(avaDictName)

    return avaDictList
