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
        avaDictName = {'simName': 'release1FP_null_dfa_2.000',
                    'Simulation Parameters': {
                        'type': 'list',
                        'Release Area': 'release1FP',
                        'Entrainment Area': '',
                        'Resistance Area': '',
                        'Parameter variation on': 'RelTh',
                        'Parameter value': '2.000',
                        'Mu': '0.155',
                        'Release thickness [m]': '2.000',
                        'Release Mass [kg]': '40000.',
                        'Final Mass [kg]': '40000.'},
                    'Release area': {'type': 'columns', 'Release area scenario': 'release1FP'},
                    'Test Info': {'type': 'text', 'Test Info': 'This test runs on a flat plane geometry.'}}

    elif avaDictName == 'avaHelixChannelDict':
        avaDictName = {'simName': 'release1HX_entres_dfa_0.055',
                    'Simulation Parameters': {
                        'type': 'list',
                        'Release Area': 'release1HX',
                        'Entrainment Area': 'entrainment1HX',
                        'Resistance Area': '',
                        'Parameter variation on': 'Mu',
                        'Parameter value': '0.055',
                        'Mu': '0.055',
                        'Release thickness [m]': '1.000',
                        'Release Mass [kg]': '22398.8',
                        'Final Mass [kg]': '23117.6'},
                    'Release area': {'type': 'columns', 'Release area scenario': 'release1HX'},
                    'Entrainment area': {'type': 'columns', 'Entrainment area scenario': 'entrainment1HX'},
                    'Resistance area': {'type': 'columns', 'Resistance area scenario': ''},
                    'Test Info': {'type': 'text',
                    'Test Info': 'This test uses a helix-shaped geometry with a channel with a channel.'}}

    elif avaDictName == 'avaHockeyDict':
        avaDictName = {'simName': 'release1HS_null_dfa_0.500',
                    'Simulation Parameters': {
                        'type': 'list',
                        'Release Area': 'release1HS',
                        'Entrainment Area': '',
                        'Resistance Area': '',
                        'Parameter variation on': 'RelTh',
                        'Parameter value': '0.500',
                        'Mu': '0.155',
                        'Release thickness [m]': '0.500',
                        'Release Mass [kg]': '10328.6',
                        'Final Mass [kg]': '10328.6'},
                    'Release area': {'type': 'columns', 'Release area scenario': 'release1HS'},
                    'Test Info': {'type': 'text',
                    'Test Info': 'This test runs on a parabolically sloping surface with a flat foreland.'}}


    return avaDictName
