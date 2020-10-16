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
