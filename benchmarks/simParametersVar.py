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

    # set desired benchmark simulation info dictionary
    if avaDictName == 'avaFlatPlaneDict':
        avaDictName = {'simName': 'release1FP_null_dfa_2.00000',
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
                        'Final Mass [kg]': '40000.'},
                    'Release area': {'type': 'columns', 'Release area scenario': 'release1FP'},
                    'Test Info': {'type': 'text', 'Test Info': 'This test runs on a flat plane geometry.'}}

    elif avaDictName == 'avaHelixChannelDict':
        avaDictName = {'simName': 'release1HX_entres_dfa_0.05500',
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
                        'Release Mass [kg]': '36865.6',
                        'Final Mass [kg]': '37584.5'},
                    'Release area': {'type': 'columns', 'Release area scenario': 'release1HX'},
                    'Entrainment area': {'type': 'columns', 'Entrainment area scenario': 'entrainment1HX'},
                    'Resistance area': {'type': 'columns', 'Resistance area scenario': ''},
                    'Test Info': {'type': 'text',
                    'Test Info': 'This test uses a helix-shaped geometry with a channel with a channel.'}}

    elif avaDictName == 'avaHockeyDict':
        avaDictName = {'simName': 'release1HS_entres_dfa_0.50000',
                    'Simulation Parameters': {
                        'type': 'list',
                        'Release Area Scenario': 'release1HS',
                        'Release Area': ['Rel_Example'],
                        'Entrainment Area': '',
                        'Resistance Area': 'resistance1HS',
                        'Parameter variation on': 'RelTh',
                        'Parameter value': '0.50000',
                        'Mu': '0.15500',
                        'Release thickness [m]': '0.50000',
                        'Release Mass [kg]': '10328.6',
                        'Final Mass [kg]': '10328.6'},
                    'Release area': {'type': 'columns', 'Release area scenario': 'release1HS'},
                    'Test Info': {'type': 'text',
                    'Test Info': 'This test runs on a parabolically sloping surface with a flat foreland.'}}


    return avaDictName
