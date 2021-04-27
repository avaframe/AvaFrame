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
        		'simType': 'null',
                    'Simulation Parameters': {
                        'type': 'list',
                        'Program version': 'macro',
                        'Friction model': 'samosAT',
                        'Release Area Scenario': 'release1BL',
                        'Release Area': ['Rel_Example'],
                        'Entrainment': 'No',
                        'Resistance': 'No',
                        'Mu': '0.15500',
                        'Density [kgm-3]': '200',
                        'Release thickness [m]': [1],
                        'Initial mass [kg]': '20196151.98',
                        'Final mass [kg]': '20196151.98',
                        'Entrained mass [kg]': '0.00',
                        'Entrainment thickness [m]': 0.3,
                        'Stop criterion': '< 1.00 percent of PKE',
                        'Avalanche run time [s]': '190.1'},
                    'Release area': {'type': 'columns', 'Release area scenario': 'release1BL'},
                    'Test Info': {'type': 'text', 'Test Info': 'This test uses a bowl-shaped geometry.'}}

        avaDictList.append(avaDictName)

    elif avaDictName == 'avaFlatPlaneDict':
        avaDictName = {'simName': {'type': 'simName', 'name': 'release1FP_null_dfa_0.15500'},
        		'testName': 'avaFlatPlaneNullTest',
        		'simType': 'null',
                    'Simulation Parameters': {
                        'type': 'list',
                        'Program version': 'macro',
                        'Friction model': 'samosAT',
                        'Release Area Scenario': 'release1FP',
                        'Release Area': ['Rel_Example'],
                        'Entrainment': 'No',
                        'Resistance': 'No',
                        'Mu': '0.15500',
                        'Density [kgm-3]': '200',
                        'Release thickness [m]': [1.0],
                        'Initial mass [kg]': '20000000.',
                        'Final mass [kg]': '20000000.',
                        'Entrained mass [kg]': '0.00',
                        'Entrainment thickness [m]': 0.3,
                        'Stop criterion': '< 1.00 percent of PKE',
                        'Avalanche run time [s]': '18.5'},
                    'Release area': {'type': 'columns', 'Release area scenario': 'release1FP'},
                    'Test Info': {'type': 'text', 'Test Info': 'This test runs on a flat plane geometry.'}}

        avaDictList.append(avaDictName)

    elif avaDictName == 'avaHelixDict':
        avaDictName = {'simName': {'type': 'simName', 'name': 'release1HX_null_dfa_0.15500'},
        		'testName': 'avaHelixNullTest',
        		'simType': 'null',
                    'Simulation Parameters': {
                        'type': 'list',
                        'Program version': 'macro',
                        'Friction model': 'samosAT',
                        'Release Area Scenario': 'release1HX',
                        'Release Area': ['Rel_Example'],
                        'Entrainment': 'No',
                        'Resistance': 'No',
                        'Mu': '0.15500',
                        'Release thickness [m]': [1],
                        'Density [kgm-3]': '200',
                        'Initial mass [kg]': '20522418.00',
                        'Final mass [kg]': '20522418.00',
                        'Entrained mass [kg]': '0.00',
                        'Entrainment thickness [m]': 0.3,
                        'Stop criterion': '< 1.00 percent of PKE',
                        'Avalanche run time [s]': '144.5'},
                    'Release area': {'type': 'columns', 'Release area scenario': 'release1HX'},
                    'Test Info': {'type': 'text', 'Test Info': 'This test uses a helix-shaped geometry.'}}

        avaDictList.append(avaDictName)

    elif avaDictName == 'avaHelixChannelDict':

        avaDictName = {'simName': {'type': 'simName', 'name': 'release1HX_ent_dfa_0.15500'},
        		'testName': 'avaHelixChannelEntTest',
        		'simType': 'ent',
                    'Simulation Parameters': {
                        'type': 'list',
                        'Program version': 'macro',
                        'Friction model': 'samosAT',
                        'Release Area Scenario': 'release1HX',
                        'Release Area': ['Rel_Example', 'Rel_Two'],
                        'Entrainment': 'Yes',
                        'Resistance': 'No',
                        'Mu': '0.15500',
                        'Density [kgm-3]': '200',
                        'Release thickness [m]': [1, 1],
                        'Entrainment thickness [m]': 0.3,
                        'Initial mass [kg]': '34354423.99',
                        'Final mass [kg]': '35018418.98',
                        'Entrained mass [kg]': '-na-',
                        'Stop criterion': '< 1.00 percent of PKE',
                        'Avalanche run time [s]': '141.5',
                        'run time [s]': ''},
                    'Release area': {'type': 'columns', 'Release area scenario': 'release1HX'},
                    'Entrainment area': {'type': 'columns', 'Entrainment area scenario': 'entrainment1HX'},
                    'Resistance area': {'type': 'columns', 'Resistance area scenario': ''},
                    'Test Info': {'type': 'text',
                    'Test Info': 'This test uses a helix-shaped geometry with a channel with a channel.'}}

        avaDictList.append(avaDictName)

        avaDictName = {'simName': {'type': 'simName', 'name': 'release1HX_ent_dfa_0.15500'},
        		'testName': 'avaHelixChannelEnt1mTest',
        		'simType': 'ent',
                    'Simulation Parameters': {
                        'type': 'list',
                        'Program version': 'macro',
                        'Friction model': 'samosAT',
                        'Release Area Scenario': 'release1HX',
                        'Release Area': ['Rel_Example', 'Rel_Two'],
                        'Entrainment': 'Yes',
                        'Resistance': 'No',
                        'Mu': '0.15500',
                        'Density [kgm-3]': '200',
                        'Release thickness [m]': [1, 1],
                        'Entrainment thickness [m]': 1.0,
                        'Initial mass [kg]': '34354423.99',
                        'Final mass [kg]': '36353187.28',
                        'Entrained mass [kg]': '-na-',
                        'Stop criterion': '< 1.00 percent of PKE',
                        'Avalanche run time [s]': '139.1',
                        'run time [s]': ''},
                    'Release area': {'type': 'columns', 'Release area scenario': 'release1HX'},
                    'Entrainment area': {'type': 'columns', 'Entrainment area scenario': 'entrainment1HX', 'Entrainment Thickness': '1'},
                    'Resistance area': {'type': 'columns', 'Resistance area scenario': ''},
                    'Test Info': {'type': 'text',
                    'Test Info': 'This test uses a helix-shaped geometry with a channel with a channel.'}}

        avaDictList.append(avaDictName)

    elif avaDictName == 'avaParabolaDict':

        avaDictName = {'simName': {'type': 'simName', 'name': 'release1PF_res_dfa_0.15500'},
        		'testName': 'avaParabolaResTest',
        		'simType': 'res',
                    'Simulation Parameters': {
                        'type': 'list',
                        'Program version': 'macro',
                        'Friction model': 'samosAT',
                        'Release Area Scenario': 'release1PF',
                        'Release Area': ['Rel_Example'],
                       'Entrainment': 'No',
                        'Resistance': 'Yes',
                        'Mu': '0.15500',
                        'Release thickness [m]': [1.0],
                        'Density [kgm-3]': '200',
                        'Initial mass [kg]': '20657132.08',
                        'Final mass [kg]': '20657132.08',
                        'Entrainment thickness [m]': 0.3,
                        'Entrained mass [kg]': '0.00',
                        'Stop criterion': '< 1.00 percent of PKE',
                        'Avalanche run time [s]': '159.1',
                        'run time [s]': ''},
                    'Release area': {'type': 'columns', 'Release area scenario': 'release1PF'},
                    'Test Info': {'type': 'text',
                    'Test Info': 'This test runs on a parabolically sloping surface with a flat foreland.'}}

        avaDictList.append(avaDictName)


    elif avaDictName == 'avaHockeyChannelDict':

        avaDictName = {'simName': {'type': 'simName', 'name': 'release1HS_ent_dfa_0.15500'},
        		'testName': 'avaHockeyChannelEntTest',
        		'simType': 'ent',
                    'Simulation Parameters': {
                        'type': 'list',
                        'Program version': 'macro',
                        'Friction model': 'samosAT',
                        'Release Area Scenario': 'release1HS',
                        'Release Area': ['Rel_Example'],
                        'Entrainment': 'Yes',
                        'Resistance': 'No',
                        'Mu': '0.15500',
                        'Density [kgm-3]': '200',
                        'Release thickness [m]': [1.0],
                        'Entrained mass [kg]': '-na-',
                        'Entrainment thickness [m]': 0.3,
                        'Initial mass [kg]': ' 19939294.00',
                        'Final mass [kg]': '20287572.96',
                        'Stop criterion': '< 1.00 percent of PKE',
                        'Avalanche run time [s]': '301.0',
                        'run time [s]': ''},
                    'Release area': {'type': 'columns', 'Release area scenario': 'release1HS'},
                    'Entrainment area': {'type': 'columns', 'Entrainment area scenario': 'entrainment1HS'},
                    'Resistance area': {'type': 'columns', 'Resistance area scenario': ''},
                    'Test Info': {'type': 'text', 'Test Info': 'This test uses a hockey stick-shaped geometry, \
                     where a linearly sloping surface transitions smoothly into a flat foreland.'}}

        avaDictList.append(avaDictName)

        avaDictName = {'simName': {'type': 'simName', 'name': 'release2HS_ent_dfa_0.15500'},
        		'testName': 'avaHockeyChannelEntTest',
        		'simType': 'ent',
                    'Simulation Parameters': {
                        'type': 'list',
                        'Program version': 'macro',
                        'Friction model': 'samosAT',
                        'Release Area Scenario': 'release2HS',
                        'Release Area': ['Rel_one', 'Rel_two'],
                        'Entrainment': 'Yes',
                        'Resistance': 'No',
                        'Mu': '0.15500',
                        'Density [kgm-3]': '200',
                        'Release thickness [m]': [1.0, 1.0],
                        'Entrainment thickness [m]': 0.3,
                        'Entrained mass [kg]': '-na-',
                        'Initial mass [kg]': '  24569688.015',
                        'Final mass [kg]': '24569688.01',
                        'Stop criterion': 'end Time reached: 400.00',
                        'Avalanche run time [s]': '400.00',
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
       			'simType': 'null',
                    'Simulation Parameters': {
                        'type': 'list',
                        'Program version': 'macro',
                        'Friction model': 'samosAT',
                        'Release Area Scenario': 'release1HS',
                        'Release Area': ['Rel_Example'],
                        'Entrainment': 'No',
                        'Resistance': 'No',
                        'Mu': '0.15500',
                        'Release thickness [m]': [1.0],
                        'Density [kgm-3]': '200',
                        'Initial mass [kg]': ' 10000000.04',
                        'Final mass [kg]': ' 10000000.04',
                        'Entrainment thickness [m]': 0.3,
                        'Entrained mass [kg]': '0.0',
                        'Stop criterion': 'end Time reached: 400.00',
                        'Avalanche run time [s]': '400.00',
                        'run time [s]': ''},
                    'Release area': {'type': 'columns', 'Release area scenario': 'release1HS'},
                    'Test Info': {'type': 'text', 'Test Info': 'This test uses a hockey stick-shaped geometry, \
                     where a linearly sloping surface transitions smoothly into a flat foreland. \
                     This geometry also includes a channel.'}}

        avaDictList.append(avaDictName)

    elif avaDictName == 'avaInclinedPlaneDict':
        avaDictName = {'simName': {'type': 'simName', 'name': 'release1IP_entres_dfa_0.15500'},
        		'testName': 'avaInclinedPlaneEntresTest',
        		'simType': 'entres',
                    'Simulation Parameters': {
                        'type': 'list',
                        'Program version': 'macro',
                        'Friction model': 'samosAT',
                        'Release Area Scenario': 'release1IP',
                        'Release Area': ['Rel_Example'],
                        'Entrainment': 'Yes',
                        'Resistance': 'Yes',
                        'Mu': '0.15500',
                        'Release thickness [m]': [1.0],
                        'Density [kgm-3]': '200',
                        'Initial mass [kg]': '20000000.09',
                        'Final mass [kg]': '21735111.13',
                        'Entrained mass [kg]': '-na-',
                        'Entrainment thickness [m]': 0.3,
                        'Stop criterion': 'end Time reached: 400.00',
                        'Avalanche run time [s]': '400.00',
                        'run time [s]': ''},
                    'Release area': {'type': 'columns', 'Release area scenario': 'release1IP'},
                    'Entrainment area': {'type': 'columns', 'Entrainment area scenario': 'entrainment1IP'},
                    'Resistance area': {'type': 'columns', 'Resistance area scenario': ''},
                    'Test Info': {'type': 'text', 'Test Info': 'This test runs on a linearly sloping surface.'}}

        avaDictList.append(avaDictName)

    elif avaDictName == 'avaAlrDict':
        avaDictName = {'simName': {'type': 'simName', 'name': 'relAlr_null_dfa_0.15500'},
        		'testName': 'avaAlrNullTest',
        		'simType': 'null',
                    'Simulation Parameters': {
                        'type': 'list',
                        'Program version': 'macro',
                        'Friction model': 'samosAT',
                        'Release Area Scenario': 'relAlr',
                        'Release Area': ['AlR'],
                        'Entrainment': 'No',
                        'Resistance': 'No',
                        'Mu': '0.15500',
                        'Release thickness [m]': [1.0],
                        'Density [kgm-3]': '200',
                        'Initial mass [kg]': '10634661.97',
                        'Final mass [kg]': '10634661.97',
                        'Entrained mass [kg]': '0.00',
                        'Entrainment thickness [m]': 0.3,
                        'Stop criterion': '< 1.00 percent of PKE',
                        'Avalanche run time [s]': '189.2',
                        'run time [s]': ''},
                    'Release area': {'type': 'columns', 'Release area scenario': 'relAlr'},
                    'Entrainment area': {'type': 'columns', 'Entrainment area scenario': ''},
                    'Resistance area': {'type': 'columns', 'Resistance area scenario': ''},
                    'Test Info': {'type': 'text', 'Test Info': 'This test runs on a Alr DEM.'}}

        avaDictList.append(avaDictName)

    elif avaDictName == 'avaKotDict':
        avaDictName = {'simName': {'type': 'simName', 'name': 'relKot_null_dfa_0.15500'},
        		'testName': 'avaKotNullTest',
        		'simType': 'null',
                    'Simulation Parameters': {
                        'type': 'list',
                        'Program version': 'macro',
                        'Friction model': 'samosAT',
                        'Release Area Scenario': 'relKot',
                        'Release Area': ['KoT'],
                        'Entrainment': 'No',
                        'Resistance': 'No',
                        'Mu': '0.15500',
                        'Density [kgm-3]': '200',
                        'Initial mass [kg]': '4325818.99',
                        'Final mass [kg]': '4325818.99',
                        'Entrained mass [kg]': '0.00',
                        'Release thickness [m]': [1.0],
                        'Entrainment thickness [m]': 0.3,
                        'Stop criterion': '< 1.00 percent of PKE',
                        'Avalanche run time [s]': '112.2',
                        'run time [s]': ''},
                    'Release area': {'type': 'columns', 'Release area scenario': 'relKot'},
                    'Test Info': {'type': 'text', 'Test Info': 'This test runs the Kot test.'}}

        avaDictList.append(avaDictName)

    elif avaDictName == 'avaMalDict':
        avaDictName = {'simName': {'type': 'simName', 'name': 'relMal1to3_null_dfa_0.15500'},
        		'testName': 'avaMalNullTest',
        		'simType': 'null',
                    'Simulation Parameters': {
                        'type': 'list',
                        'Program version': 'macro',
                        'Friction model': 'samosAT',
                        'Release Area Scenario': 'relMal1to3',
                        'Release Area': ['MaL2', 'MaL3','MaL1'],
                       'Entrainment': 'No',
                        'Resistance': 'No',
                        'Mu': '0.15500',
                        'Density [kgm-3]': '200',
                        'Initial mass [kg]': '87904128.03',
                        'Final mass [kg]': '87904128.03',
                        'Entrained mass [kg]': '0.00',
                        'Release thickness [m]': [1.0, 1.0, 1.0],
                        'Entrainment thickness [m]': 0.3,
                        'Stop criterion': '< 1.00 percent of PKE',
                        'Avalanche run time [s]': '176.1',
                        'run time [s]': ''},
                    'Release area': {'type': 'columns', 'Release area scenario': 'relMal1to3'},
                    'Test Info': {'type': 'text', 'Test Info': 'This test runs the Mal test.'}}

        avaDictList.append(avaDictName)

    elif avaDictName == 'avaWogDict':
        avaDictName = {'simName': {'type': 'simName', 'name': 'relWog_null_dfa_0.15500'},
        		'testName': 'avaWogNullTest',
        		'simType': 'null',
                    'Simulation Parameters': {
                        'type': 'list',
                        'Program version': 'macro',
                        'Friction model': 'samosAT',
                        'Release Area Scenario': 'relWog',
                        'Release Area': ['WoG'],
                        'Entrainment': 'No',
                        'Resistance': 'No',
                        'Mu': '0.15500',
                        'Density [kgm-3]': '200',
                        'Initial mass [kg]': '33090744.01',
                        'Final mass [kg]': '33090744.01',
                        'Entrained mass [kg]': '0.00',
                        'Release thickness [m]': [1.0],
                        'Entrainment thickness [m]': 0.3,
                        'Stop criterion': '< 1.00 percent of PKE',
                        'Avalanche run time [s]': '165.5',
                        'run time [s]': ''},
                    'Release area': {'type': 'columns', 'Release area scenario': 'relWog'},
                    'Test Info': {'type': 'text', 'Test Info': 'This test runs the Wog test.'}}

        avaDictList.append(avaDictName)

    elif avaDictName == 'avaGarDict':
        avaDictName = {'simName': {'type': 'simName', 'name': 'relGar_null_dfa_0.15500'},
        		'testName': 'avaGarNullTest',
        		'simType': 'null',
                    'Simulation Parameters': {
                        'type': 'list',
                        'Program version': 'macro',
                        'Friction model': 'samosAT',
                        'Release Area Scenario': 'relGar',
                        'Release Area': ['GaR1'],
                        'Entrainment': 'No',
                        'Resistance': 'No',
                        'Mu': '0.15500',
                        'Density [kgm-3]': '200',
                        'Initial mass [kg]': '12355490.01',
                        'Final mass [kg]': '12355490.01',
                        'Entrained mass [kg]': '0.00',
                        'Release thickness [m]': [1.2],
                        'Entrainment thickness [m]': 0.3,
                        'Stop criterion':'< 1.00 percent of PKE',
                        'Avalanche run time [s]': '123.9',
                        'run time [s]': ''},
                    'Release area': {'type': 'columns', 'Release area scenario': 'relGar'},
                    'Test Info': {'type': 'text', 'Test Info': 'This test runs the Gar test.'}}

        avaDictList.append(avaDictName)

        avaDictName = {'simName': {'type': 'simName', 'name': 'relGar2_null_dfa_0.15500'},
        		'testName': 'avaGarNullTest',
        		'simType': 'null',
                    'Simulation Parameters': {
                        'type': 'list',
                        'Program version': 'macro',
                        'Friction model': 'samosAT',
                        'Release Area Scenario': 'relGar2',
                        'Release Area': ['GaR2'],
                        'Entrainment': 'No',
                        'Resistance': 'No',
                        'Mu': '0.15500',
                        'Density [kgm-3]': '200',
                        'Initial mass [kg]': '67116520.00',
                        'Final mass [kg]': '67116520.00',
                        'Entrained mass [kg]': '0.00',
                        'Release thickness [m]': [1.5],
                        'Entrainment thickness [m]': 0.3,
                        'Stop criterion': '< 1.00 percent of PKE',
                        'Avalanche run time [s]': '104.9',
                        'run time [s]': ''},
                    'Release area': {'type': 'columns', 'Release area scenario': 'relGar2'},
                    'Test Info': {'type': 'text', 'Test Info': 'This test runs the Gar test.'}}

        avaDictList.append(avaDictName)

        avaDictName = {'simName': {'type': 'simName', 'name': 'relGar6_null_dfa_0.15500'},
        		'testName': 'avaGarNullTest',
        		'simType': 'null',
                    'Simulation Parameters': {
                        'type': 'list',
                        'Program version': 'macro',
                        'Friction model': 'samosAT',
                        'Release Area Scenario': 'relGar6',
                        'Release Area': ['GaR6'],
                        'Entrainment': 'No',
                        'Resistance': 'No',
                        'Mu': '0.15500',
                        'Density [kgm-3]': '200',
                        'Initial mass [kg]': '22167257.98',
                        'Final mass [kg]': '22167257.98',
                        'Entrained mass [kg]': '0.00',
                        'Release thickness [m]': [1.0],
                        'Entrainment thickness [m]': 0.3,
                        'Stop criterion': '< 1.00 percent of PKE',
                        'Avalanche run time [s]': '96.7',
                        'run time [s]': ''},
                    'Release area': {'type': 'columns', 'Release area scenario': 'relGar6'},
                    'Test Info': {'type': 'text', 'Test Info': 'This test runs the Gar test.'}}

        avaDictList.append(avaDictName)


    elif avaDictName == 'avaHitDict':
        avaDictName = {'simName': {'type': 'simName', 'name': 'relHit_null_dfa_0.15500'},
        		'testName': 'avaHitNullTest',
        		'simType': 'null',
                    'Simulation Parameters': {
                        'type': 'list',
                        'Program version': 'macro',
                        'Friction model': 'samosAT',
                        'Release Area Scenario': 'relHit',
                        'Release Area': ['HiT'],
                       'Entrainment': 'No',
                        'Resistance': 'No',
                        'Mu': '0.15500',
                        'Density [kgm-3]': '200',
                        'Initial mass [kg]': '53584431.99',
                        'Final mass [kg]': '53584431.99',
                        'Entrained mass [kg]': '0.00',
                        'Release thickness [m]': [1.0],
                        'Entrainment thickness [m]': 0.3,
                        'Stop criterion': '< 1.00 percent of PKE',
                        'Avalanche run time [s]': '158.0',
                        'run time [s]': ''},
                    'Release area': {'type': 'columns', 'Release area scenario': 'relHit'},
                    'Test Info': {'type': 'text', 'Test Info': 'This test runs the Hit test.'}}

        avaDictList.append(avaDictName)


    elif avaDictName == 'avaPyramidDict':
        avaDictName = {'simName': {'type': 'simName', 'name': 'release1PY_null_dfa_0.15500'},
        		'testName': 'avaPyramidNullTest',
        		'simType': 'null',
                    'Simulation Parameters': {
                        'type': 'list',
                        'Program version': 'macro',
                        'Friction model': 'samosAT',
                        'Release Area Scenario': 'release1PY',
                        'Release Area': ['rel1'],
                        'Entrainment': 'No',
                        'Resistance': 'No',
                        'Mu': '0.15500',
                        'Density [kgm-3]': '200',
                        'Release thickness [m]': [1.0],
                        'Entrainment thickness [m]': 0.3,
                        'run time [s]': ''},
                    'Release area': {'type': 'columns', 'Release area scenario': 'release1PY'},
                    'Test Info': {'type': 'text', 'Test Info': 'This test runs the pyramid test.'}}

        avaDictList.append(avaDictName)


    elif avaDictName == 'avaPyramid45Dict':
        avaDictName = {'simName': {'type': 'simName', 'name': 'release1PY_null_dfa_0.15500'},
        		'testName': 'avaPyramid45NullTest',
        		'simType': 'null',
                    'Simulation Parameters': {
                        'type': 'list',
                        'Program version': 'macro',
                        'Friction model': 'samosAT',
                        'Release Area Scenario': 'release1PY',
                        'Release Area': ['rel'],
                       'Entrainment': 'No',
                        'Resistance': 'No',
                        'Mu': '0.15500',
                        'Density [kgm-3]': '200',
                        'Release thickness [m]': [1.0],
                        'Entrainment thickness [m]': 0.3,
                        'run time [s]': ''},
                    'Release area': {'type': 'columns', 'Release area scenario': 'release1PY'},
                    'Test Info': {'type': 'text', 'Test Info': 'This test runs the rotation pyramid test.'}}

        avaDictList.append(avaDictName)

    return avaDictList
