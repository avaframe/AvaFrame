"""
    dictionary with simulation info for benchmarks
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
                        'Release Area Scenario': 'release1BL',
                        'Release Area': ['Rel_Example'],
                        'Entrainment Area': '',
                        'Resistance Area': '',
                        'Mu': '0.15500',
                        'Release thickness [m]': ['1'],
                        'Release Mass [kg]': '20196.2',
                        'Final Mass [kg]': '20196.2'},
                    'Release area': {'type': 'columns', 'Release area scenario': 'release1BL'},
                    'Test Info': {'type': 'text', 'Test Info': 'This test uses a bowl-shaped geometry.'}}

        avaDictList.append(avaDictName)

    elif avaDictName == 'avaFlatPlaneDict':
        avaDictName = {'simName': 'release1FP_null_dfa_0.15500',
                    'Simulation Parameters': {
                        'type': 'list',
                        'Release Area Scenario': 'release1FP',
                        'Release Area': ['Rel_Example'],
                        'Entrainment Area': '',
                        'Resistance Area': '',
                        'Mu': '0.15500',
                        'Release thickness [m]': ['1.0'],
                        'Release Mass [kg]': '20000.',
                        'Final Mass [kg]': '20000.'},
                    'Release area': {'type': 'columns', 'Release area scenario': 'release1FP'},
                    'Test Info': {'type': 'text', 'Test Info': 'This test runs on a flat plane geometry.'}}

        avaDictList.append(avaDictName)

    elif avaDictName == 'avaHelixDict':
        avaDictName = {'simName': 'release1HX_null_dfa_0.15500',
                    'Simulation Parameters': {
                        'type': 'list',
                        'Release Area Scenario': 'release1HX',
                        'Release Area': ['Rel_Example'],
                        'Entrainment Area': '',
                        'Resistance Area': '',
                        'Mu': '0.15500',
                        'Release thickness [m]': ['1'],
                        'Release Mass [kg]': '20522.4',
                        'Final Mass [kg]': '20522.4'},
                    'Release area': {'type': 'columns', 'Release area scenario': 'release1HX'},
                    'Test Info': {'type': 'text', 'Test Info': 'This test uses a helix-shaped geometry.'}}

        avaDictList.append(avaDictName)

    elif avaDictName == 'avaHelixChannelDict':

        avaDictName = {'simName': 'release1HX_entres_dfa_0.15500',
                    'Simulation Parameters': {
                        'type': 'list',
                        'Release Area Scenario': 'release1HX',
                        'Release Area': ['Rel_Example', 'Rel_Two'],
                        'Entrainment Area': 'entrainment1HX',
                        'Resistance Area': '',
                        'Mu': '0.15500',
                        'Release thickness [m]': ['1', '1'],
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
                        'Release Area Scenario': 'release1HS',
                        'Release Area': ['Rel_Example'],
                        'Entrainment Area': '',
                        'Resistance Area': 'resistance1HS',
                        'Mu': '0.15500',
                        'Release thickness [m]': ['1.0'],
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
                        'Release Area Scenario': 'release1HS2',
                        'Release Area': ['Rel_Example'],
                        'Entrainment Area': 'entrainment1HS2',
                        'Resistance Area': '',
                        'Mu': '0.15500',
                        'Release thickness [m]': ['1.0'],
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
                        'Release Area Scenario': 'release2HS2',
                        'Release Area': ['Rel_Example'],
                        'Entrainment Area': 'entrainment1HS2',
                        'Resistance Area': '',
                        'Mu': '0.15500',
                        'Release thickness [m]': ['1.0'],
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
                        'Release Area Scenario': 'release1HS2',
                        'Release Area': ['Rel_Example'],
                        'Entrainment Area': '',
                        'Resistance Area': '',
                        'Mu': '0.15500',
                        'Release thickness [m]': ['1.0'],
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
                        'Release Area Scenario': 'release1IP',
                        'Release Area': ['Rel_Example'],
                        'Entrainment Area': 'entrainment1IP',
                        'Resistance Area': 'resistance1IP',
                        'Mu': '0.15500',
                        'Release thickness [m]': ['1.0'],
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
                        'Release Area Scenario': 'relAlr',
                        'Release Area': ['AlR'],
                        'Entrainment Area': '',
                        'Resistance Area': '',
                        'Mu': '0.15500',
                        'Release thickness [m]': ['1.0'],
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
                        'Release Area Scenario': 'relKot',
                        'Release Area': ['KoT'],
                        'Entrainment Area': '',
                        'Resistance Area': '',
                        'Mu': '0.15500',
                        'Release thickness [m]': ['1.0']},
                    'Release area': {'type': 'columns', 'Release area scenario': 'relKot'},
                    'Test Info': {'type': 'text', 'Test Info': 'This test runs the Kot test.'}}

        avaDictList.append(avaDictName)

    elif avaDictName == 'avaMalDict':
        avaDictName = {'simName': 'relMal1to3_null_dfa_0.15500',
                    'Simulation Parameters': {
                        'type': 'list',
                        'Release Area Scenario': 'relMal1to3',
                        'Release Area': ['MaL2', 'MaL3','MaL1'],
                        'Entrainment Area': '',
                        'Resistance Area': '',
                        'Mu': '0.15500',
                        'Release thickness [m]': ['1.0', '1.0', '1.0']},
                    'Release area': {'type': 'columns', 'Release area scenario': 'relMal1to3'},
                    'Test Info': {'type': 'text', 'Test Info': 'This test runs the Mal test.'}}

        avaDictList.append(avaDictName)

    elif avaDictName == 'avaWogDict':
        avaDictName = {'simName': 'relWog_null_dfa_0.15500',
                    'Simulation Parameters': {
                        'type': 'list',
                        'Release Area Scenario': 'relWog',
                        'Release Area': ['WoG'],
                        'Entrainment Area': '',
                        'Resistance Area': '',
                        'Mu': '0.15500',
                        'Release thickness [m]': ['1.0']},
                    'Release area': {'type': 'columns', 'Release area scenario': 'relWog'},
                    'Test Info': {'type': 'text', 'Test Info': 'This test runs the Wog test.'}}

        avaDictList.append(avaDictName)

    elif avaDictName == 'avaGarDict':
        avaDictName = {'simName': 'relGar_null_dfa_0.15500',
                    'Simulation Parameters': {
                        'type': 'list',
                        'Release Area Scenario': 'relGar',
                        'Release Area': ['GaR1'],
                        'Entrainment Area': '',
                        'Resistance Area': '',
                        'Mu': '0.15500',
                        'Release thickness [m]': ['1.2']},
                    'Release area': {'type': 'columns', 'Release area scenario': 'relGar'},
                    'Test Info': {'type': 'text', 'Test Info': 'This test runs the Gar test.'}}

        avaDictList.append(avaDictName)

        avaDictName = {'simName': 'relGar23_null_dfa_0.15500',
                    'Simulation Parameters': {
                        'type': 'list',
                        'Release Area Scenario': 'relGar23',
                        'Release Area': ['GaR2', 'GaR3'],
                        'Entrainment Area': '',
                        'Resistance Area': '',
                        'Mu': '0.15500',
                        'Release thickness [m]': ['1.5','1.35']},
                    'Release area': {'type': 'columns', 'Release area scenario': 'relGar23'},
                    'Test Info': {'type': 'text', 'Test Info': 'This test runs the Gar test.'}}

        avaDictList.append(avaDictName)

        avaDictName = {'simName': 'relGar6_null_dfa_0.15500',
                    'Simulation Parameters': {
                        'type': 'list',
                        'Release Area Scenario': 'relGar6',
                        'Release Area': ['GaR6'],
                        'Entrainment Area': '',
                        'Resistance Area': '',
                        'Mu': '0.15500',
                        'Release thickness [m]': ['1.0']},
                    'Release area': {'type': 'columns', 'Release area scenario': 'relGar6'},
                    'Test Info': {'type': 'text', 'Test Info': 'This test runs the Gar test.'}}

        avaDictList.append(avaDictName)


    elif avaDictName == 'avaHitDict':
        avaDictName = {'simName': 'relHit_null_dfa_0.15500',
                    'Simulation Parameters': {
                        'type': 'list',
                        'Release Area Scenario': 'relHit',
                        'Release Area': ['HiT'],
                        'Entrainment Area': '',
                        'Resistance Area': '',
                        'Mu': '0.15500',
                        'Release thickness [m]': ['1.0']},
                    'Release area': {'type': 'columns', 'Release area scenario': 'relHit'},
                    'Test Info': {'type': 'text', 'Test Info': 'This test runs the Hit test.'}}

        avaDictList.append(avaDictName)

    return avaDictList
