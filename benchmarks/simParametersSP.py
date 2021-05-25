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
    if avaDictName == 'avaFlatPlaneDict':
        avaDictName = {'simName': {'type': 'simName', 'name': 'release1FP_null_dfa_f8fca36398'},
        	       'testName': 'avaFlatPlaneNullTestPy',
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
                        'Initial mass [kg]': '19500000.00',
                        'Final mass [kg]': '19500000.00',
                        'Entrained mass [kg]': '0.00',
                        'Entrainment thickness [m]': 0.3,
                        'Stop criterion': '< 1.00 percent of PKE',
                        'Avalanche run time [s]': '17.80'},
                    'Release area': {'type': 'columns', 'Release area scenario': 'release1FP'},
                    'Test Info': {'type': 'text', 'Test Info': 'This test runs on a flat plane geometry.'}}

        avaDictList.append(avaDictName)

    elif avaDictName == 'avaHelixChannelDict':

        avaDictName = {'simName': {'type': 'simName', 'name': 'release1HX_ent_dfa_f005e4e9c5'},
        	       'testName': 'avaHelixChannelEntTestPy',
        	       'simType': 'ent',
                    'Simulation Parameters': {
                        'type': 'list',
                        'Program version': 'com1DFAPy',
                        'Friction model': 'samosAT',
                        'Release Area Scenario': 'release1HX',
                        'Release Area': ['Rel_Example', 'Rel_Two'],
                        'Entrainment': 'Yes',
                        'Resistance': 'No',
                        'Mu': '0.15500',
                        'Density [kgm-3]': '200',
                        'Release thickness [m]': [1, 1],
                        'Entrainment thickness [m]': 0.3,
                        'Initial mass [kg]': '34264950.92',
                        'Final mass [kg]': '34921704.69',
                        'Entrained mass [kg]': '656753.77',
                        'Stop criterion': '< 1.00 percent of PKE',
                        'Avalanche run time [s]': '140.80',
                        'run time [s]': ''},
                    'Release area': {'type': 'columns', 'Release area scenario': 'release1HX'},
                    'Entrainment area': {'type': 'columns', 'Entrainment area scenario': 'entrainment1HX'},
                    'Resistance area': {'type': 'columns', 'Resistance area scenario': ''},
                    'Test Info': {'type': 'text',
                    'Test Info': 'This test uses a helix-shaped geometry with a channel with a channel.'}}

        avaDictList.append(avaDictName)

   


    elif avaDictName == 'avaHockeyChannelDict':

        avaDictName = {'simName': {'type': 'simName', 'name': 'release1HS_ent_dfa_5afbc252c7'},
        	       'testName': 'avaHockeyChannelEntTestPy',
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
                        'Entrained mass [kg]': '366651.77',
                        'Entrainment thickness [m]': 0.3,
                        'Initial mass [kg]': '19512767.70',
                        'Final mass [kg]': '19879419.47',
                        'Stop criterion': '< 1.00 percent of PKE',
                        'Avalanche run time [s]': '309.40',
                        'run time [s]': ''},
                    'Release area': {'type': 'columns', 'Release area scenario': 'release1HS'},
                    'Entrainment area': {'type': 'columns', 'Entrainment area scenario': 'entrainment1HS'},
                    'Resistance area': {'type': 'columns', 'Resistance area scenario': ''},
                    'Test Info': {'type': 'text', 'Test Info': 'This test uses a hockey stick-shaped geometry, \
                     where a linearly sloping surface transitions smoothly into a flat foreland.'}}

        avaDictList.append(avaDictName)
        avaDictName = {'simName': {'type': 'simName', 'name': 'release2HS_ent_dfa_4adbe5324d'},
        		'testName': 'avaHockeyChannelEntTestPy',
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
                        'Entrained mass [kg]': '0.0',
                        'Initial mass [kg]': '24580273.69',
                        'Final mass [kg]': '24580273.69',
                        'Stop criterion': 'end Time reached: 400.00',
                        'Avalanche run time [s]': '400.00',
                        'run time [s]': ''},
                    'Release area': {'type': 'columns', 'Release area scenario': 'release1HS'},
                    'Entrainment area': {'type': 'columns', 'Entrainment area scenario': 'entrainment1HS'},
                    'Resistance area': {'type': 'columns', 'Resistance area scenario': ''},
                    'Test Info': {'type': 'text', 'Test Info': 'This test uses a hockey stick-shaped geometry, \
                     where a linearly sloping surface transitions smoothly into a flat foreland.'}}

        avaDictList.append(avaDictName)
        
    elif avaDictName == 'avaInclinedPlaneDict':
        
        avaDictName = {'simName': {'type': 'simName', 'name': 'release1IP_entres_dfa_eeaf249f5c'},
        	       'testName': 'avaInclinedPlaneEntresTestPy',
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
                        'Initial mass [kg]': '19570886.22',
                        'Final mass [kg]': '21313398.02',
                        'Entrained mass [kg]': '1742511.80',
                        'Entrainment thickness [m]': 0.3,
                        'Stop criterion': 'end Time reached: 400.00',
                        'Avalanche run time [s]': '400.00',
                        'run time [s]': ''},
                    'Release area': {'type': 'columns', 'Release area scenario': 'release1IP'},
                    'Entrainment area': {'type': 'columns', 'Entrainment area scenario': 'entrainment1IP'},
                    'Resistance area': {'type': 'columns', 'Resistance area scenario': ''},
                    'Test Info': {'type': 'text', 'Test Info': 'This test runs on a linearly sloping surface.'}}

        avaDictList.append(avaDictName) 

    return avaDictList
