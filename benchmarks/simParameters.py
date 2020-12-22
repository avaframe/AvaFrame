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
    if avaDictName == 'avaBowlDict':
        avaDictName = {'simName': 'release1BL_null_dfa_0.15500',
                    'Simulation Parameters': {
                        'type': 'list',
                        'Release Area': 'release1BL',
                        'Entrainment Area': '',
                        'Resistance Area': '',
                        'Mu': '0.15500',
                        'Release thickness [m]': '1.00000',
                        'Release Mass [kg]': '20196.2',
                        'Final Mass [kg]': '20196.2'},
                    'Release area': {'type': 'columns', 'Release area scenario': 'release1BL'},
                    'Test Info': {'type': 'text', 'Test Info': 'This test uses a bowl-shaped geometry.'}}

        avaDictList.append(avaDictName)

    elif avaDictName == 'avaFlatPlaneDict':
        avaDictName = {'simName': 'release1FP_null_dfa_0.15500',
                    'Simulation Parameters': {
                        'type': 'list',
                        'Release Area': 'release1FP',
                        'Entrainment Area': '',
                        'Resistance Area': '',
                        'Mu': '0.15500',
                        'Release thickness [m]': '1.00000',
                        'Release Mass [kg]': '20000.',
                        'Final Mass [kg]': '20000.'},
                    'Release area': {'type': 'columns', 'Release area scenario': 'release1FP'},
                    'Test Info': {'type': 'text', 'Test Info': 'This test runs on a flat plane geometry.'}}

        avaDictList.append(avaDictName)

    elif avaDictName == 'avaHelixDict':
        avaDictName = {'simName': 'release1HX_null_dfa_0.15500',
                    'Simulation Parameters': {
                        'type': 'list',
                        'Release Area': 'release1HX',
                        'Entrainment Area': '',
                        'Resistance Area': '',
                        'Mu': '0.15500',
                        'Release thickness [m]': '1.00000',
                        'Release Mass [kg]': '20522.4',
                        'Final Mass [kg]': '20522.4'},
                    'Release area': {'type': 'columns', 'Release area scenario': 'release1HX'},
                    'Test Info': {'type': 'text', 'Test Info': 'This test uses a helix-shaped geometry.'}}

        avaDictList.append(avaDictName)

    elif avaDictName == 'avaHelixChannelDict':

        avaDictName = {'simName': 'release1HX_entres_dfa_0.15500',
                    'Simulation Parameters': {
                        'type': 'list',
                        'Release Area': 'release1HX',
                        'Entrainment Area': 'entrainment1HX',
                        'Resistance Area': '',
                        'Mu': '0.15500',
                        'Release thickness [m]': '1.00000',
                        'Release Mass [kg]': '22398.8',
                        'Final Mass [kg]': '23117.6'},
                    'Release area': {'type': 'columns', 'Release area scenario': 'release1HX'},
                    'Entrainment area': {'type': 'columns', 'Entrainment area scenario': 'entrainment1HX'},
                    'Resistance area': {'type': 'columns', 'Resistance area scenario': ''},
                    'Test Info': {'type': 'text',
                    'Test Info': 'This test uses a helix-shaped geometry with a channel with a channel.'}}

        avaDictList.append(avaDictName)

    elif avaDictName == 'avaHockeyDict':

        avaDictName = {'simName': 'release1HS_entres_dfa_0.15500',
                    'Simulation Parameters': {
                        'type': 'list',
                        'Release Area': 'release1HS',
                        'Entrainment Area': '',
                        'Resistance Area': 'resistance1HS',
                        'Mu': '0.15500',
                        'Release thickness [m]': '1.00000',
                        'Release Mass [kg]': '20657.1',
                        'Final Mass [kg]': '20657.1'},
                    'Release area': {'type': 'columns', 'Release area scenario': 'release1HS'},
                    'Test Info': {'type': 'text',
                    'Test Info': 'This test runs on a parabolically sloping surface with a flat foreland.'}}

        avaDictList.append(avaDictName)


    elif avaDictName == 'avaHockeySmoothChannelDict':

        avaDictName = {'simName': 'release1HS2_entres_dfa_0.15500',
                    'Simulation Parameters': {
                        'type': 'list',
                        'Release Area': 'release1HS2',
                        'Entrainment Area': 'entrainment1HS2',
                        'Resistance Area': '',
                        'Mu': '0.15500',
                        'Release thickness [m]': '1.00000',
                        'Release Mass [kg]': '20967.3',
                        'Final Mass [kg]': '21306.'},
                    'Release area': {'type': 'columns', 'Release area scenario': 'release1HS2'},
                    'Entrainment area': {'type': 'columns', 'Entrainment area scenario': 'entrainment1HS2'},
                    'Resistance area': {'type': 'columns', 'Resistance area scenario': ''},
                    'Test Info': {'type': 'text', 'Test Info': 'This test uses a hockey stick-shaped geometry, \
                     where a linearly sloping surface transitions smoothly into a flat foreland.'}}

        avaDictList.append(avaDictName)

        avaDictName = {'simName': 'release2HS2_entres_dfa_0.15500',
                    'Simulation Parameters': {
                        'type': 'list',
                        'Release Area': 'release2HS2',
                        'Entrainment Area': 'entrainment1HS2',
                        'Resistance Area': '',
                        'Mu': '0.15500',
                        'Release thickness [m]': '1.00000',
                        'Release Mass [kg]': '26627.4',
                        'Final Mass [kg]': '26627.4'},
                    'Release area': {'type': 'columns', 'Release area scenario': 'release1HS2'},
                    'Entrainment area': {'type': 'columns', 'Entrainment area scenario': 'entrainment1HS2'},
                    'Resistance area': {'type': 'columns', 'Resistance area scenario': ''},
                    'Test Info': {'type': 'text', 'Test Info': 'This test uses a hockey stick-shaped geometry, \
                     where a linearly sloping surface transitions smoothly into a flat foreland.'}}

        avaDictList.append(avaDictName)

    elif avaDictName == 'avaHockeySmoothSmallDict':
        avaDictName = {'simName': 'release1HS2_null_dfa_0.15500',
                    'Simulation Parameters': {
                        'type': 'list',
                        'Release Area': 'release1HS2',
                        'Entrainment Area': '',
                        'Resistance Area': '',
                        'Mu': '0.15500',
                        'Release thickness [m]': '1.00000',
                        'Release Mass [kg]': '10000.',
                        'Final Mass [kg]': '10000.'},
                    'Release area': {'type': 'columns', 'Release area scenario': 'release1HS2'},
                    'Test Info': {'type': 'text', 'Test Info': 'This test uses a hockey stick-shaped geometry, \
                     where a linearly sloping surface transitions smoothly into a flat foreland. \
                     This geometry also includes a channel.'}}

        avaDictList.append(avaDictName)

    elif avaDictName == 'avaInclinedPlaneDict':
        avaDictName = {'simName': 'release1IP_entres_dfa_0.15500',
                    'Simulation Parameters': {
                        'type': 'list',
                        'Release Area': 'release1IP',
                        'Entrainment Area': 'entrainment1IP',
                        'Resistance Area': 'resistance1IP',
                        'Mu': '0.15500',
                        'Release thickness [m]': '1.00000',
                        'Release Mass [kg]': '20000.',
                        'Final Mass [kg]': '21735.1'},
                    'Release area': {'type': 'columns', 'Release area scenario': 'release1IP'},
                    'Entrainment area': {'type': 'columns', 'Entrainment area scenario': 'entrainment1IP'},
                    'Resistance area': {'type': 'columns', 'Resistance area scenario': ''},
                    'Test Info': {'type': 'text', 'Test Info': 'This test runs on a linearly sloping surface.'}}

        avaDictList.append(avaDictName)

    elif avaDictName == 'avaAlrDict':
        avaDictName = {'simName': 'relAlr_null_dfa_0.15500',
                    'Simulation Parameters': {
                        'type': 'list',
                        'Release Area': 'relAlr',
                        'Entrainment Area': '',
                        'Resistance Area': '',
                        'Mu': '0.15500',
                        'Release thickness [m]': '1.00000',
                        'Release Mass [kg]': '	10586.3',
                        'Final Mass [kg]': '10586.3'},
                    'Release area': {'type': 'columns', 'Release area scenario': 'relAlr'},
                    'Entrainment area': {'type': 'columns', 'Entrainment area scenario': ''},
                    'Resistance area': {'type': 'columns', 'Resistance area scenario': ''},
                    'Test Info': {'type': 'text', 'Test Info': 'This test runs on a Alr DEM.'}}

        avaDictList.append(avaDictName)

    elif avaDictName == 'avaKotDict':
        avaDictName = {'simName': 'relKot_null_dfa_0.15500',
                    'Simulation Parameters': {
                        'type': 'list',
                        'Release Area': 'relKot',
                        'Entrainment Area': '',
                        'Resistance Area': '',
                        'Mu': '0.15500',
                        'Release thickness [m]': '1.00000'},
                    'Release area': {'type': 'columns', 'Release area scenario': 'relKot'},
                    'Test Info': {'type': 'text', 'Test Info': 'This test runs the Kot test.'}}

        avaDictList.append(avaDictName)

    elif avaDictName == 'avaMalDict':
        avaDictName = {'simName': 'relMal_null_dfa_0.15500',
                    'Simulation Parameters': {
                        'type': 'list',
                        'Release Area': 'relMal',
                        'Entrainment Area': '',
                        'Resistance Area': '',
                        'Mu': '0.15500',
                        'Release thickness [m]': '1.00000'},
                    'Release area': {'type': 'columns', 'Release area scenario': 'relMal'},
                    'Test Info': {'type': 'text', 'Test Info': 'This test runs the Mal test.'}}

        avaDictList.append(avaDictName)

    elif avaDictName == 'avaWogDict':
        avaDictName = {'simName': 'relWog_null_dfa_0.15500',
                    'Simulation Parameters': {
                        'type': 'list',
                        'Release Area': 'relWog',
                        'Entrainment Area': '',
                        'Resistance Area': '',
                        'Mu': '0.15500',
                        'Release thickness [m]': '1.00000'},
                    'Release area': {'type': 'columns', 'Release area scenario': 'relWog'},
                    'Test Info': {'type': 'text', 'Test Info': 'This test runs the Wog test.'}}

        avaDictList.append(avaDictName)


    return avaDictList
