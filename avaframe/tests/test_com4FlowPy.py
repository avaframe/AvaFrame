"""
    Pytest for module com4FlowPy
"""

#  Load modules
import numpy as np
import pytest

from avaframe.com4FlowPy import flowClass
import avaframe.com4FlowPy.flowCore as flowCore

def test_add_os():
    cell = flowClass.Cell(1,1,
                          np.array([[10,10,10], [10,10,10], [10,10,10]]), 10,
                          1,0, None,
                          20, 8, 3e-4, 270,
                          startcell=True)
    cell.add_os(0.2)
    assert cell.flux == 1.2

def test_reverseTopology():
    '''
    testing flowCore.reverseTopology() for different
    examples of dir graphs
    '''

    testGraph = {0: [1,2,3],
                      1: [2,4],
                      2: [4,5],
                      3: [2],
                      4: [],
                      5: [6],
                      6: []}
    
    testGraphReverse = {0: [],
                        1: [0],
                        2: [0,1,3],
                        3: [0],
                        4: [2,1],
                        5: [2],
                        6: [5]}
    
    reverseGraphCalc = flowCore.reverseTopology(testGraph)

    for key, item in reverseGraphCalc.items():
        assert key in testGraphReverse.keys()
        setTestChildren = set(item)
        setCalcChildren = set(testGraphReverse[key])
        assert setTestChildren == setCalcChildren
    
    testGraph = {0:[1,2,3],
                 1:[4],
                 2:[5,6],
                 3:[7,8],
                 4:[9],
                 5:[9,10],
                 6:[10],
                 7:[11],
                 8:[],
                 9:[],
                 10:[12],
                 11:[],
                 12:[]}
    
    testGraphReverse = {0:[],
                        1:[0],
                        2:[0],
                        3:[0],
                        4:[1],
                        5:[2],
                        6:[2],
                        7:[3],
                        8:[3],
                        9:[4,5],
                        10:[5,6],
                        11:[7],
                        12:[10]}
    
    reverseGraphCalc = flowCore.reverseTopology(testGraph)

    for key, item in reverseGraphCalc.items():
        assert key in testGraphReverse.keys()
        setTestChildren = set(item)
        setCalcChildren = set(testGraphReverse[key])
        assert setTestChildren == setCalcChildren
    

def test_backTracking():
    '''
    testing flowCore.backTracking() for different
    examples of dir graphs - basic graphs are the same as
    in test_reverseTopology() with added 'infra' values
    as valueDicts
    '''
    
    testGraph = {0: [1,2,3],
                 1: [2,4],
                 2: [4,5],
                 3: [2],
                 4: [],
                 5: [6],
                 6: []}
    
    testValsIn = {0: 0,
                  1: 0,
                  2: 0,
                  3: 0,
                  4: 3,
                  5: 0,
                  6: 2}
    
    testValsBT = {0: 3,
                  1: 3,
                  2: 3,
                  3: 3,
                  4: 3,
                  5: 2,
                  6: 2}
    
    calcValsBT = flowCore.backTracking(testGraph, testValsIn)

    for key,item in calcValsBT.items():
        assert calcValsBT[key] == testValsBT[key]
    
    testGraph = {0:[1,2,3],
                 1:[4],
                 2:[5,6],
                 3:[7,8],
                 4:[9],
                 5:[9,10],
                 6:[10],
                 7:[11],
                 8:[],
                 9:[],
                 10:[12],
                 11:[],
                 12:[]}
    
    testValsIn = {0:0,
                  1:0,
                  2:0,
                  3:0,
                  4:0,
                  5:0,
                  6:0,
                  7:0,
                  8:0,
                  9:1,
                  10:0,
                  11:3,
                  12:2}
    
    testValsBT = {0:3,
                  1:1,
                  2:2,
                  3:3,
                  4:1,
                  5:2,
                  6:2,
                  7:3,
                  8:0,
                  9:1,
                  10:2,
                  11:3,
                  12:2}
    
    calcValsBT = flowCore.backTracking(testGraph, testValsIn)

    for key,item in calcValsBT.items():
        assert calcValsBT[key] == testValsBT[key]


def test_calculation():
    dem = np.array([[20, 20, 20], [10, 10, 10], [0, 0, 0]])
    infra = None
    pra = np.array([[0, 1, 0], [0, 0, 0], [0, 0, 0]])
    alpha = 10
    exp = 99
    fluxTh = 1 * e - 3
    zDeltaMax = 270
    nodata = -9999
    cellsize = 10
    infraBool = False
    forestBool = False
    variableParameters = {
        "varUmaxBool": False,
        "varUmaxArray": None,
        "varAlphaBool": False,
        "varAlphaArray": None,
        "varExponentBool": False,
        "varExponentArray": None,
    }
    fluxDistOldVersionBool = False
    previewMode = False
    outputs = ['travelLengthMin']
    forestArray = None
    forestParams = None

    args[0] = dem
    args[1] = infra
    args[2] = pra
    args[3] = alpha
    args[4] = exp
    args[5] = fluxTh
    args[6] = zDeltaMax
    args[7] = nodata
    args[8] = cellsize
    args[9] = infraBool
    args[10] = forestBool
    args[11] = variableParameters
    ars[12] = fluxDistOldVersionBool
    args[13] = previewMode
    args[14] = forestArray
    args[15] = forestParams
    args[16] = outputs

    flux = np.array([-9999, 1, -9999], [-9999, 1, -9999], [-9999, 1, -9999])
    routFluxSum = flux.copy()
    depFluxSum = np.array([-9999, 0, -9999], [-9999, 0, -9999], [-9999, 0, -9999])
    travelLengthMin = np.array([-9999, 0, -9999], [-9999, np.sqrt(200), -9999], [-9999, np.sqrt(200) * 2, -9999])
    results = flowCore.calculation(args)

    assert len(results) == 12
    assert results[1] == flux
    assert results[10] == routFluxSum
    assert results[11] == depFluxSum
    assert results[9] == travelLengthMin
    assert results[8] == np.ones_like(flux) * -9999


if __name__=='__main__':
    test_add_os()
    test_reverseTopology()
    test_backTracking()