import aimecdata as aimecdata
from aimecrunner import *
import aimectaskcoordtrans
from coordtrans_process import processData
from coordtrans_process import analyze
from coordtrans_process import aimecOut

import string
import sys
import os
from matplotlib import pyplot as plt
import matplotlib
import numpy as np


def Run_aimec(fname):
    f = open(fname, 'r')
    aD = aimecdata.fromString(f.read())
    runner = AimecRunner()
    tasklist = [runner.tasks[0]]
    runner.runall(aD, tasklist)


def Run_part_aimec():
    pathPressure = '/home/matthiastonnel/Documents/github/AvaFrame/avaframe/aimec_2020/NameOfAvalanche/Outputs/dfa_pressure/'
    pathAvalanchePath = '/home/matthiastonnel/Documents/github/AvaFrame/avaframe/aimec_2020/NameOfAvalanche/Inputs/avalanche_path.xyz'
    pathResult = '/home/matthiastonnel/Documents/github/AvaFrame/avaframe/aimec_2020/NameOfAvalanche/Outputs/AimecResults/coordtrans/'
    depthfile = '/home/matthiastonnel/Documents/github/AvaFrame/avaframe/aimec_2020/NameOfAvalanche/Outputs/dfa_depth/'
    pathDHM = '/home/matthiastonnel/Documents/github/AvaFrame/avaframe/aimec_2020/NameOfAvalanche/Inputs/GBG_DGM.asc'
    calcPressureLimit = 1
    domainWidth = 600

    pressurefileList = [str(pathPressure) +
                        '/' +
                        str(name) for name in
                        sorted(os.listdir(pathPressure)) if os.path.isfile(os.path.join(pathPressure, name))]
    set_name = pressurefileList[0].split('/')[-3]
    project_name = pathAvalanchePath.split('/')[-4]
    path_name = pathAvalanchePath.split('/')[-1]

    # -----------------------------------------------------------
    # create new raster + preparing new raster assignment function
    # -----------------------------------------------------------
    deskewedRasterInd = processData.processDataInd(pressurefileList,
                                                   pathAvalanchePath,
                                                   domainWidth,
                                                   False,
                                                   '',
                                                   False,
                                                   '',
                                                   True,  # str(data.pathMapResult),
                                                   pathResult,
                                                   True,
                                                   calcPressureLimit)

    # -----------------------------------------------------------
    # transform pressure_data and depth_data in new raster
    # -----------------------------------------------------------

    deskewedRasterPressure = processData.assignData(pressurefileList, deskewedRasterInd)

    depthfileList = [str(depthfile) +
                     '/' +
                     str(name) for name in
                     sorted(os.listdir(depthfile)) if os.path.isfile(os.path.join(depthfile, name))]

    deskewedRasterDepth = processData.assignData(depthfileList, deskewedRasterInd)

    deskewedRasterDHM = processData.assignData([pathDHM], deskewedRasterInd)

    deskewedRasterDoku = None
    # -----------------------------------------------------------
    # analyze
    # -----------------------------------------------------------
    # analyze doku
    print('[MAIN] analyzing dokumentation data')
    doku, runout_doku, delta_h, elevRel = analyze.analyzeDocu(calcPressureLimit,
                                                              pressurefileList,
                                                              deskewedRasterInd,
                                                              deskewedRasterPressure,
                                                              deskewedRasterDepth,
                                                              deskewedRasterDoku,
                                                              deskewedRasterDHM,
                                                              True,
                                                              True,
                                                              False)

    [runout,
     runout_mean,
     AMPP,
     MMPP,
     AMD,
     MMD] = analyze.analyzeDataWithDepth(deskewedRasterInd,
                                         calcPressureLimit,
                                         pressurefileList, deskewedRasterPressure,
                                         deskewedRasterDepth,
                                         False,  # data.with__(data.pathMapResult),
                                         pathResult,
                                         True)

    # -----------------------------------------------------------
    # result visualisation + report
    # -----------------------------------------------------------
    relMass = 0.
    entMass = np.zeros((len(pressurefileList)))
    gr_index = np.zeros((len(pressurefileList)))
    gr_grad = np.zeros((len(pressurefileList)))
    aimecOut.result_visu(pressurefileList, pathAvalanchePath,
                         runout, AMPP, doku, gr_index, calcPressureLimit,
                         pathResult)

    # fig0, ax0 = plt.subplots()
    # ax0.imshow(deskewedRasterPressure[0], cmap=plt.cm.jet, interpolation='nearest')
    #
    # fig1, ax1 = plt.subplots()
    # ax1.imshow(deskewedRasterDepth[0], cmap=plt.cm.jet, interpolation='nearest')
    #
    plt.show()


if __name__ == '__main__':
    dirname = os.path.dirname(__file__)
    fname = os.path.join(dirname, 'NameOfAvalanche/Outputs/aimecPathFile.txt')

    # Run_aimec(fname)
    Run_part_aimec()
