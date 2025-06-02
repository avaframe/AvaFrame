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


if __name__=='__main__':
    test_add_os()
    test_reverseTopology()
    test_backTracking()