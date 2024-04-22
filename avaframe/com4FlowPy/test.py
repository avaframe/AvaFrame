from avaframe.com4FlowPy.flowClass import Cell
import numpy as np
#import pytest


def test_add_parent():

    # define an example startcell
    row_idx = 0
    col_idx = 5
    dem_ng = np.array([[1100,1100,1100],[1000, 1000, 1000], [990, 990, 990]])
    cellsize = 10
    alpha = 25
    exp = 8
    flux_threshold = 3e-4
    max_z_delta = 9000
    forest_local = 1
    
    test_startcell = Cell(row_idx, col_idx, dem_ng, cellsize, 1,0, None,
                            alpha, exp, flux_threshold, max_z_delta,forest_local, startcell=True)



    # define an example child cell
    row_idx = 1
    col_idx = 5
    dem_ng = np.array([[1000, 1000, 1000], [990, 990, 990], [980,980,980]])
    cellsize = 10
    alpha = 25
    exp = 8
    flux_threshold = 3e-4
    max_z_delta = 9000
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

test_add_parent()