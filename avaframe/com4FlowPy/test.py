from avaframe.com4FlowPy.flowClass import Cell
import avaframe.com4FlowPy.flowCore as fc

import numpy as np
#import pytest


# Test flowClass.py

def test_add_parent():
    '''
    Test method of class Cell in flowClass.py:
    Test the computation of forest interaction with different parent cells.
    '''

    # define Flow-Py input parameters (as same attributes for all cells)
    cellsize = 10
    alpha = 25
    exp = 8
    flux_threshold = 3e-4
    max_z_delta = 9000

    # define an example startcell
    row_idx = 0
    col_idx = 5
    dem_ng = np.array([[1100,1100,1100],[1000, 1000, 1000], [990, 990, 990]])
    forest_local = 1
    
    test_startcell = Cell(row_idx, col_idx, dem_ng, cellsize, 1,0, None,
                            alpha, exp, flux_threshold, max_z_delta,forest_local, startcell=True)


    # define an example child cell
    row_idx = 1
    col_idx = 5
    dem_ng = np.array([[1000, 1000, 1000], [990, 990, 990], [980,980,980]])
    forest_local = 1

    # test if hitted forest from parent cell (1) is added to test child cell (1)
    test_cell = Cell(row_idx, col_idx, dem_ng, cellsize, 0.33,2, test_startcell,
                            alpha, exp, flux_threshold, max_z_delta,forest_local, startcell=False)
    assert test_cell.hitted_forest == 2

    # test if the hitted_forest value is lowered when adding a parent with a lower 
    # hitted forest value
    parent_cell2 = Cell(row_idx, col_idx, dem_ng, cellsize, 0.33,2, None,
                            alpha, exp, flux_threshold, max_z_delta, 0, startcell=True)
    test_cell.add_parent(parent_cell2)
    assert len(test_cell.parent) == 2
    assert test_cell.hitted_forest == 1

    # test if value stays at 1 when adding a parent with a higher hitted_forest value
    parent_cell3 = Cell(row_idx, col_idx, dem_ng, cellsize, 0.33,2, test_startcell,
                            alpha, exp, flux_threshold, max_z_delta, 1, startcell=False)
    test_cell.add_parent(parent_cell2)
    assert test_cell.hitted_forest == 1
    assert len(test_cell.parent) == 3


def test_calc_distribution():
    '''
    Test (main) method of class Cell in flowClass.py:
    Test the distribution function for a simple example at an inclined plane
    without lateral spreading.
    The total flux should be distributed in 1 cell downslope of the parent cell.
    Additionally, the calculation of z_delta is tested for this example 
    (since the output includes calculated z_delta).
    '''
    # define Flow-Py input parameters (as same attributes for all cells)
    cellsize = 10
    alpha = 25
    exp = 99
    flux_threshold = 3e-4
    max_z_delta = 9000

    # define an example startcell
    row_idx = 0
    col_idx = 5
    dem_ng = np.array([[1100,1100,1100],[1000, 1000, 1000], [990, 990, 990]])
    forest_local = 1
    
    test_startcell = Cell(row_idx, col_idx, dem_ng, cellsize, 1,0, None,
                            alpha, exp, flux_threshold, max_z_delta,forest_local, startcell=True)


    # define an example child cell
    row_idx = 1
    col_idx = 5
    dem_ng = np.array([[1000, 1000, 1000], [990, 990, 990], [980,980,980]])
    forest_local = 1

    
    test_cell = Cell(row_idx, col_idx, dem_ng, cellsize, 1. ,2, test_startcell,
                            alpha, exp, flux_threshold, max_z_delta,forest_local, startcell=False)

    
    # altitude difference = 10
    # 10 = z_alpha + z_delta
    # tan(alpha) = z_alpha / cellsize

    test_cell.z_delta = 10 - cellsize * np.tan(np.deg2rad(alpha))
    result_distribution = test_startcell.calc_distribution()
    assert np.isclose(result_distribution[0], np.array(test_cell.rowindex))
    assert np.isclose(result_distribution[1], np.array(test_cell.colindex))
    assert np.isclose(result_distribution[2], test_cell.flux)
    assert np.isclose(result_distribution[3], test_cell.z_delta)


# Test flowCore.py

def test_calculation():

    #define arguments
    dem = np.array([[1000, 1000, 1000], [990, 990, 990], [980,980,980], [970,970,970], [960,960,960]])
    infra = np.zeros_like(dem)
    release = np.array([[0,0,0], [0,1,0], [0,0,0], [0,0,0], [0,0,0]])
    alpha = 5
    exp = 99
    flux_threshold = 3e-4
    max_z_delta=9000
    nodata = -9999
    cellsize = 10
    infraBool = False
    forest = None

    args = [dem,infra, release, alpha, exp,flux_threshold, max_z_delta, nodata, cellsize,
    infraBool, forest]
    
    

     # altitude difference = 10
    # 10 = z_alpha + z_delta
    # tan(alpha) = z_alpha / cellsize

    z_delta1 = 10 - cellsize * np.tan(np.deg2rad(alpha))
    z_delta2 = 20 - 2*cellsize * np.tan(np.deg2rad(alpha))
    gamma = np.rad2deg(np.arctan(10 / cellsize))
    
    z_delta_array = np.array([[0, 0, 0], [0, 0 , 0], [0, z_delta1 , 0], [0, z_delta2 , 0], [0, 0 , 0]])
    flux_array = np.array([[0, 0, 0], [0, 1, 0], [0, 1, 0], [0, 1, 0], [0, 0 , 0]])
    count_array = np.array([[0, 0, 0], [0, 1, 0], [0, 1, 0], [0, 1, 0], [0, 0 , 0]])
    z_delta_sum = z_delta_array
    backcalc = np.zeros_like(dem, dtype=np.int32)
    fp_travelangle_array = np.array([[0, 0, 0], [0, 0 , 0], [0, gamma , 0], [0, gamma , 0], [0, 0 , 0]])
    sl_travelangle_array = fp_travelangle_array
    forest_interaction_array = None
   
    assert np.allclose(fc.calculation(args)[0], z_delta_array) 
    assert np.allclose(fc.calculation(args)[1], flux_array)
    assert np.allclose(fc.calculation(args)[2],  count_array)
    assert np.allclose(fc.calculation(args)[3], z_delta_sum)
    assert np.allclose(fc.calculation(args)[4], backcalc)
    assert np.allclose(fc.calculation(args)[5], fp_travelangle_array)
    assert np.allclose(fc.calculation(args)[6], sl_travelangle_array)
    assert fc.calculation(args)[7] is None


test_add_parent()
test_calc_distribution()
test_calculation()