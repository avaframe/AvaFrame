import sys
import os
import logging
import glob

def readAIMECinputs(avalancheDir):
    pathPressure = avalancheDir + '/Outputs/dfa_pressure'
    pathFlowHeight = avalancheDir + '/Outputs/dfa_depth'

    ProfileLayer = glob.glob(avalancheDir + '/Inputs/LINES/*aimec*.shp')
    cfgPath['ProfileLayer'] = = ''.join(ProfileLayer)

    demSource = glob.glob(avalancheDir + '/Inputs/*.asc')
    try:
        assert len(demSource) == 1, 'There should be only one and only one DEM .asc file in ' + \
            avalancheDir + '/Inputs/'
    except AssertionError as e:
        raise
    cfgPath['demSource'] = ''.join(demSource)
    pressurefileList = [pathPressure +
                        '/' +
                        str(name) for name in
                        sorted(os.listdir(pathPressure)) if os.path.isfile(os.path.join(pathPressure, name))]
    cfgPath['pressurefileList'] = pressurefileList

    depthfileList = [str(pathFlowHeight) +
                     '/' +
                     str(name) for name in
                     sorted(os.listdir(pathFlowHeight)) if os.path.isfile(os.path.join(pathFlowHeight, name))]
    cfgPath['depthfileList'] = depthfileList

    pathResult = avalancheDir + '/Outputs/AimecResults'
    cfgPath['pathResult'] = pathResult

    defaultName = str(avalancheDir).split('/')[-1]
    cfgPath['defaultName'] = defaultName

    set_name = pressurefileList[0].split('/')[-3]
    cfgPath['set_name'] = set_name
    project_name = str(ProfileLayer).split('/')[-4]
    cfgPath['project_name'] = project_name
    path_name = str(ProfileLayer).split('/')[-1]
    cfgPath['path_name'] = path_name

    print(set_name, project_name, path_name)



    return cfgPath
