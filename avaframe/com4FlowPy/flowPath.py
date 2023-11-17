import numpy as np
import matplotlib.pyplot as plt
import math

class Path:
    '''Class which contains a path, containing one startcell  and corresponding child cells'''
    
    def __init__(self, dem, start_row, start_col, gen_list):
        self.dem = dem
        self.gen_list = gen_list
        self.start_row = start_row
        self.start_col = start_col

        self.z_delta_array = np.zeros_like(dem, dtype=np.float32)
        self.flux_array = np.zeros_like(dem, dtype=np.float32)
        self.travel_length_array = np.zeros_like(dem, dtype=np.float32)
        self.flow_energy_array = np.zeros_like(dem, dtype=np.float32)

    def get_path_arrays(self):
        for gen, cell_list in enumerate(self.gen_list):
            for cell in cell_list:
                self.z_delta_array[cell.rowindex, cell.colindex] = max(self.z_delta_array[cell.rowindex, cell.colindex], cell.z_delta)
                self.flux_array[cell.rowindex, cell.colindex] = max(self.flux_array[cell.rowindex, cell.colindex], cell.flux)
                self.travel_length_array[cell.rowindex, cell.colindex] = max(self.travel_length_array[cell.rowindex, cell.colindex], cell.min_distance)
                self.flow_energy_array[cell.rowindex, cell.colindex] = max(self.flow_energy_array[cell.rowindex, cell.colindex], cell.flow_energy)

    def plot_test(self):
        self.get_path_arrays()
        fig, axs = plt.subplots()
        axs.imshow(self.dem, cmap ='Greys', alpha=0.8)
        axs.contour(self.dem, levels = 10, colors ='k',linewidths=0.5)
        axs.imshow(self.z_delta_array, cmap = 'Blues', alpha = 0.6)
        fig.savefig(f'/home/paula/data/Flowpy_test/plane/output_1cell_PRA/plots/path_analysis_TEST.png')
        plt.close(fig)


