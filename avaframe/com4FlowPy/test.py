from avaframe.com4FlowPy.flowClass import Cell
import numpy as np

row_idx = 1
col_idx = 5
dem_ng = np.array([[1000, 1000, 1000], [990, 990, 990], [980,980,980]])
cellsize = 10
alpha = 25
exp = 8
flux_threshold = 3e-4
max_z_delta = 9000
forest_local = [[0,0,0], [0,1,0], [1,0,0]]
# startcell:
test_cell = Cell(row_idx, col_idx, dem_ng, cellsize, 1,0, None,
                         alpha, exp, flux_threshold, max_z_delta,forest_local, startcell=False)

def test_forest_int():
    #assert test_cell.hitted_forest == 1

    #advanced
    test_cell.parent[0].hitted_forest = 2
    assert test_cell.hitted_forest == 3


test_forest_int()