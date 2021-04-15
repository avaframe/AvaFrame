"""
    dictionary with simulation info for benchmarks
"""

# Load modules
import os
from avaframe.in3Utils import fileHandlerUtils as fU


def fetchBenchParameters(avaName):
    """ Collect simulation parameter info from standard tests """

    # get name of avalanche
    avaDictName = avaName + 'Dict'
    avaDictList = []

    # set desired benchmark simulation info dictionary
    if avaDictName == 'avaBowlDict':
        avaDictName = {'simName': {'type': 'simName', 'name': 'release1BL_null_dfa_0.15500'},
        		'testName': 'avaBowlNullTest',
                    'Simulation Parameters': {
                        'type': 'list',
                        'Release Area Scenario': 'release1BL',
                        'Release Area': ['Rel_Example'],
                        'Entrainment Area': '',
                        'Resistance Area': '',
                        'Mu': '0.15500',
                        'Release thickness [m]': [1],
                        'Release Mass [kg]': '20196.2',
                        'Final Mass [kg]': '20196.2',
                        'Entrainment thickness [m]': 0.3,
                        'run time [s]': ''},
                    'Release area': {'type': 'columns', 'Release area scenario': 'release1BL'},
                    'Test Info': {'type': 'text', 'Test Info': 'This test uses a bowl-shaped geometry.'}}

        avaDictList.append(avaDictName)

    elif avaDictName == 'avaFlatPlaneDict':
        avaDictName = {'simName': {'type': 'simName', 'name': 'release1FP_null_dfa_0.15500'},
        		'testName': 'avaFlatPlaneNullTest',
                    'Simulation Parameters': {
                        'type': 'list',
                        'Release Area Scenario': 'release1FP',
                        'Release Area': ['Rel_Example'],
                        'Entrainment Area': '',
                        'Resistance Area': '',
                        'Mu': '0.15500',
                        'Release thickness [m]': [1.0],
                        'Release Mass [kg]': '20000.',
                        'Final Mass [kg]': '20000.',
                        'Entrainment thickness [m]': 0.3,
                        'run time [s]': ''},
                    'Release area': {'type': 'columns', 'Release area scenario': 'release1FP'},
                    'Test Info': {'type': 'text', 'Test Info': 'This test runs on a flat plane geometry.'}}

        avaDictList.append(avaDictName)

    elif avaDictName == 'avaHelixDict':
        avaDictName = {'simName': {'type': 'simName', 'name': 'release1HX_null_dfa_0.15500'},
        		'testName': 'avaHelixNullTest',
                    'Simulation Parameters': {
                        'type': 'list',
                        'Release Area Scenario': 'release1HX',
                        'Release Area': ['Rel_Example'],
                        'Entrainment Area': '',
                        'Resistance Area': '',
                        'Mu': '0.15500',
                        'Release thickness [m]': [1],
                        'Release Mass [kg]': '20522.4',
                        'Final Mass [kg]': '20522.4',
                        'Entrainment thickness [m]': 0.3,
                        'run time [s]': ''},
                    'Release area': {'type': 'columns', 'Release area scenario': 'release1HX'},
                    'Test Info': {'type': 'text', 'Test Info': 'This test uses a helix-shaped geometry.'}}

        avaDictList.append(avaDictName)

    elif avaDictName == 'avaHelixChannelDict':

        avaDictName = {'simName': {'type': 'simName', 'name': 'release1HX_entres_dfa_0.15500'},
        		'testName': 'avaHelixChannelEntTest',
                    'Simulation Parameters': {
                        'type': 'list',
                        'Release Area Scenario': 'release1HX',
                        'Release Area': ['Rel_Example', 'Rel_Two'],
                        'Entrainment Area': 'entrainment1HX',
                        'Resistance Area': '',
                        'Mu': '0.15500',
                        'Release thickness [m]': [1, 1],
                        'Entrainment thickness [m]': 0.3,
                        'run time [s]': ''},
                    'Release area': {'type': 'columns', 'Release area scenario': 'release1HX'},
                    'Entrainment area': {'type': 'columns', 'Entrainment area scenario': 'entrainment1HX'},
                    'Resistance area': {'type': 'columns', 'Resistance area scenario': ''},
                    'Test Info': {'type': 'text',
                    'Test Info': 'This test uses a helix-shaped geometry with a channel with a channel.'}}

        avaDictList.append(avaDictName)

        avaDictName = {'simName': {'type': 'simName', 'name': 'release1HX_entres_dfa_0.15500'},
        		'testName': 'avaHelixChannelEnt1mTest',
                    'Simulation Parameters': {
                        'type': 'list',
                        'Release Area Scenario': 'release1HX',
                        'Release Area': ['Rel_Example', 'Rel_Two'],
                        'Entrainment Area': 'entrainment1HX',
                        'Resistance Area': '',
                        'Mu': '0.15500',
                        'Release thickness [m]': [1, 1],
                        'Entrainment thickness [m]': 1.0,
                        'run time [s]': ''},
                    'Release area': {'type': 'columns', 'Release area scenario': 'release1HX'},
                    'Entrainment area': {'type': 'columns', 'Entrainment area scenario': 'entrainment1HX', 'Entrainment Thickness': '1'},
                    'Resistance area': {'type': 'columns', 'Resistance area scenario': ''},
                    'Test Info': {'type': 'text',
                    'Test Info': 'This test uses a helix-shaped geometry with a channel with a channel.'}}

        avaDictList.append(avaDictName)

    elif avaDictName == 'avaParabolaDict':

        avaDictName = {'simName': {'type': 'simName', 'name': 'release1PF_entres_dfa_0.15500'},
        		'testName': 'avaParabolaResTest',
                    'Simulation Parameters': {
                        'type': 'list',
                        'Release Area Scenario': 'release1PF',
                        'Release Area': ['Rel_Example'],
                        'Entrainment Area': '',
                        'Resistance Area': 'resistance1PF',
                        'Mu': '0.15500',
                        'Release thickness [m]': [1.0],
                        'Release Mass [kg]': '20657.1',
                        'Final Mass [kg]': '20657.1',
                        'Entrainment thickness [m]': 0.3,
                        'run time [s]': ''},
                    'Release area': {'type': 'columns', 'Release area scenario': 'release1PF'},
                    'Test Info': {'type': 'text',
                    'Test Info': 'This test runs on a parabolically sloping surface with a flat foreland.'}}

        avaDictList.append(avaDictName)


    elif avaDictName == 'avaHockeyChannelDict':

        avaDictName = {'simName': {'type': 'simName', 'name': 'release1HS_entres_dfa_0.15500'},
        		'testName': 'avaHockeyChannelEntTest',
                    'Simulation Parameters': {
                        'type': 'list',
                        'Release Area Scenario': 'release1HS',
                        'Release Area': ['Rel_Example'],
                        'Entrainment Area': 'entrainment1HS',
                        'Resistance Area': '',
                        'Mu': '0.15500',
                        'Release thickness [m]': [1.0],
                        'Entrainment thickness [m]': 0.3,
                        'run time [s]': ''},
                    'Release area': {'type': 'columns', 'Release area scenario': 'release1HS'},
                    'Entrainment area': {'type': 'columns', 'Entrainment area scenario': 'entrainment1HS'},
                    'Resistance area': {'type': 'columns', 'Resistance area scenario': ''},
                    'Test Info': {'type': 'text', 'Test Info': 'This test uses a hockey stick-shaped geometry, \
                     where a linearly sloping surface transitions smoothly into a flat foreland.'}}

        avaDictList.append(avaDictName)

        avaDictName = {'simName': {'type': 'simName', 'name': 'release2HS_entres_dfa_0.15500'},
        		'testName': 'avaHockeyChannelEntTest',
                    'Simulation Parameters': {
                        'type': 'list',
                        'Release Area Scenario': 'release2HS',
                        'Release Area': ['Rel_one', 'Rel_two'],
                        'Entrainment Area': 'entrainment1HS',
                        'Resistance Area': '',
                        'Mu': '0.15500',
                        'Release thickness [m]': [1.0, 1.0],
                        'Entrainment thickness [m]': 0.3,
                        'run time [s]': ''},
                    'Release area': {'type': 'columns', 'Release area scenario': 'release1HS'},
                    'Entrainment area': {'type': 'columns', 'Entrainment area scenario': 'entrainment1HS'},
                    'Resistance area': {'type': 'columns', 'Resistance area scenario': ''},
                    'Test Info': {'type': 'text', 'Test Info': 'This test uses a hockey stick-shaped geometry, \
                     where a linearly sloping surface transitions smoothly into a flat foreland.'}}

        avaDictList.append(avaDictName)

    elif avaDictName == 'avaHockeySmallDict':
        avaDictName = {'simName': {'type': 'simName', 'name': 'release1HS_null_dfa_0.15500'},
       			'testName': 'avaHockeySmallNullTest',
                    'Simulation Parameters': {
                        'type': 'list',
                        'Release Area Scenario': 'release1HS',
                        'Release Area': ['Rel_Example'],
                        'Entrainment Area': '',
                        'Resistance Area': '',
                        'Mu': '0.15500',
                        'Release thickness [m]': [1.0],
                        'Final Mass [kg]': '10000.',
                        'Entrainment thickness [m]': 0.3,
                        'run time [s]': ''},
                    'Release area': {'type': 'columns', 'Release area scenario': 'release1HS'},
                    'Test Info': {'type': 'text', 'Test Info': 'This test uses a hockey stick-shaped geometry, \
                     where a linearly sloping surface transitions smoothly into a flat foreland. \
                     This geometry also includes a channel.'}}

        avaDictList.append(avaDictName)

    elif avaDictName == 'avaInclinedPlaneDict':
        avaDictName = {'simName': {'type': 'simName', 'name': 'release1IP_entres_dfa_0.15500'},
        		'testName': 'avaInclinedPlaneEntresTest',
                    'Simulation Parameters': {
                        'type': 'list',
                        'Release Area Scenario': 'release1IP',
                        'Release Area': ['Rel_Example'],
                        'Entrainment Area': 'entrainment1IP',
                        'Resistance Area': 'resistance1IP',
                        'Mu': '0.15500',
                        'Release thickness [m]': [1.0],
                        'Release Mass [kg]': '20000.',
                        'Final Mass [kg]': '21735.1',
                        'Entrainment thickness [m]': 0.3,
                        'run time [s]': ''},
                    'Release area': {'type': 'columns', 'Release area scenario': 'release1IP'},
                    'Entrainment area': {'type': 'columns', 'Entrainment area scenario': 'entrainment1IP'},
                    'Resistance area': {'type': 'columns', 'Resistance area scenario': ''},
                    'Test Info': {'type': 'text', 'Test Info': 'This test runs on a linearly sloping surface.'}}

        avaDictList.append(avaDictName)

    elif avaDictName == 'avaAlrDict':
        avaDictName = {'simName': {'type': 'simName', 'name': 'relAlr_null_dfa_0.15500'},
        		'testName': 'avaAlrNullTest',
                    'Simulation Parameters': {
                        'type': 'list',
                        'Release Area Scenario': 'relAlr',
                        'Release Area': ['AlR'],
                        'Entrainment Area': '',
                        'Resistance Area': '',
                        'Mu': '0.15500',
                        'Release thickness [m]': [1.0],
                        'Release Mass [kg]': '	10586.3',
                        'Final Mass [kg]': '10586.3',
                        'Entrainment thickness [m]': 0.3,
                        'run time [s]': ''},
                    'Release area': {'type': 'columns', 'Release area scenario': 'relAlr'},
                    'Entrainment area': {'type': 'columns', 'Entrainment area scenario': ''},
                    'Resistance area': {'type': 'columns', 'Resistance area scenario': ''},
                    'Test Info': {'type': 'text', 'Test Info': 'This test runs on a Alr DEM.'}}

        avaDictList.append(avaDictName)

    elif avaDictName == 'avaKotDict':
        avaDictName = {'simName': {'type': 'simName', 'name': 'relKot_null_dfa_0.15500'},
        		'testName': 'avaKotNullTest',
                    'Simulation Parameters': {
                        'type': 'list',
                        'Release Area Scenario': 'relKot',
                        'Release Area': ['KoT'],
                        'Entrainment Area': '',
                        'Resistance Area': '',
                        'Mu': '0.15500',
                        'Release thickness [m]': [1.0],
                        'Entrainment thickness [m]': 0.3,
                        'run time [s]': ''},
                    'Release area': {'type': 'columns', 'Release area scenario': 'relKot'},
                    'Test Info': {'type': 'text', 'Test Info': 'This test runs the Kot test.'}}

        avaDictList.append(avaDictName)

    elif avaDictName == 'avaMalDict':
        avaDictName = {'simName': {'type': 'simName', 'name': 'relMal1to3_null_dfa_0.15500'},
        		'testName': 'avaMalNullTest',
                    'Simulation Parameters': {
                        'type': 'list',
                        'Release Area Scenario': 'relMal1to3',
                        'Release Area': ['MaL2', 'MaL3','MaL1'],
                        'Entrainment Area': '',
                        'Resistance Area': '',
                        'Mu': '0.15500',
                        'Release thickness [m]': [1.0, 1.0, 1.0],
                        'Entrainment thickness [m]': 0.3,
                        'run time [s]': ''},
                    'Release area': {'type': 'columns', 'Release area scenario': 'relMal1to3'},
                    'Test Info': {'type': 'text', 'Test Info': 'This test runs the Mal test.'}}

        avaDictList.append(avaDictName)

    elif avaDictName == 'avaWogDict':
        avaDictName = {'simName': {'type': 'simName', 'name': 'relWog_null_dfa_0.15500'},
        		'testName': 'avaWogNullTest',
                    'Simulation Parameters': {
                        'type': 'list',
                        'Release Area Scenario': 'relWog',
                        'Release Area': ['WoG'],
                        'Entrainment Area': '',
                        'Resistance Area': '',
                        'Mu': '0.15500',
                        'Release thickness [m]': [1.0],
                        'Entrainment thickness [m]': 0.3,
                        'run time [s]': ''},
                    'Release area': {'type': 'columns', 'Release area scenario': 'relWog'},
                    'Test Info': {'type': 'text', 'Test Info': 'This test runs the Wog test.'}}

        avaDictList.append(avaDictName)

    elif avaDictName == 'avaGarDict':
        avaDictName = {'simName': {'type': 'simName', 'name': 'relGar_null_dfa_0.15500'},
        		'testName': 'avaGarNullTest',
                    'Simulation Parameters': {
                        'type': 'list',
                        'Release Area Scenario': 'relGar',
                        'Release Area': ['GaR1'],
                        'Entrainment Area': '',
                        'Resistance Area': '',
                        'Mu': '0.15500',
                        'Release thickness [m]': [1.2],
                        'Entrainment thickness [m]': 0.3,
                        'run time [s]': ''},
                    'Release area': {'type': 'columns', 'Release area scenario': 'relGar'},
                    'Test Info': {'type': 'text', 'Test Info': 'This test runs the Gar test.'}}

        avaDictList.append(avaDictName)

        avaDictName = {'simName': {'type': 'simName', 'name': 'relGar2_null_dfa_0.15500'},
        		'testName': 'avaGarNullTest',
                    'Simulation Parameters': {
                        'type': 'list',
                        'Release Area Scenario': 'relGar2',
                        'Release Area': ['GaR2'],
                        'Entrainment Area': '',
                        'Resistance Area': '',
                        'Mu': '0.15500',
                        'Release thickness [m]': [1.5],
                        'Entrainment thickness [m]': 0.3,
                        'run time [s]': ''},
                    'Release area': {'type': 'columns', 'Release area scenario': 'relGar2'},
                    'Test Info': {'type': 'text', 'Test Info': 'This test runs the Gar test.'}}

        avaDictList.append(avaDictName)

        avaDictName = {'simName': {'type': 'simName', 'name': 'relGar6_null_dfa_0.15500'},
        		'testName': 'avaGarNullTest',
                    'Simulation Parameters': {
                        'type': 'list',
                        'Release Area Scenario': 'relGar6',
                        'Release Area': ['GaR6'],
                        'Entrainment Area': '',
                        'Resistance Area': '',
                        'Mu': '0.15500',
                        'Release thickness [m]': [1.0],
                        'Entrainment thickness [m]': 0.3,
                        'run time [s]': ''},
                    'Release area': {'type': 'columns', 'Release area scenario': 'relGar6'},
                    'Test Info': {'type': 'text', 'Test Info': 'This test runs the Gar test.'}}

        avaDictList.append(avaDictName)


    elif avaDictName == 'avaHitDict':
        avaDictName = {'simName': {'type': 'simName', 'name': 'relHit_null_dfa_0.15500'},
        		'testName': 'avaHitNullTest',
                    'Simulation Parameters': {
                        'type': 'list',
                        'Release Area Scenario': 'relHit',
                        'Release Area': ['HiT'],
                        'Entrainment Area': '',
                        'Resistance Area': '',
                        'Mu': '0.15500',
                        'Release thickness [m]': [1.0],
                        'Entrainment thickness [m]': 0.3,
                        'run time [s]': ''},
                    'Release area': {'type': 'columns', 'Release area scenario': 'relHit'},
                    'Test Info': {'type': 'text', 'Test Info': 'This test runs the Hit test.'}}

        avaDictList.append(avaDictName)


    elif avaDictName == 'avaPyramidDict':
        avaDictName = {'simName': {'type': 'simName', 'name': 'release1PY_null_dfa_0.15500'},
        		'testName': 'avaPyramidNullTest',
                    'Simulation Parameters': {
                        'type': 'list',
                        'Release Area Scenario': 'release1PY',
                        'Release Area': ['rel1'],
                        'Entrainment Area': '',
                        'Resistance Area': '',
                        'Mu': '0.15500',
                        'Release thickness [m]': [1.0],
                        'Entrainment thickness [m]': 0.3,
                        'run time [s]': ''},
                    'Release area': {'type': 'columns', 'Release area scenario': 'release1PY'},
                    'Test Info': {'type': 'text', 'Test Info': 'This test runs the pyramid test.'}}

        avaDictList.append(avaDictName)


    elif avaDictName == 'avaPyramid45Dict':
        avaDictName = {'simName': {'type': 'simName', 'name': 'release1PY_null_dfa_0.15500'},
        		'testName': 'avaPyramid45NullTest',
                    'Simulation Parameters': {
                        'type': 'list',
                        'Release Area Scenario': 'release1PY',
                        'Release Area': ['rel'],
                        'Entrainment Area': '',
                        'Resistance Area': '',
                        'Mu': '0.15500',
                        'Release thickness [m]': [1.0],
                        'Entrainment thickness [m]': 0.3,
                        'run time [s]': ''},
                    'Release area': {'type': 'columns', 'Release area scenario': 'release1PY'},
                    'Test Info': {'type': 'text', 'Test Info': 'This test runs the rotation pyramid test.'}}

        avaDictList.append(avaDictName)

    return avaDictList
