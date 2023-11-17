import numpy as np
import matplotlib.pyplot as plt
import math

class Path:
    '''Class contains a path, containing one startcell  and corresponding child cells'''
    
    def __init__(self, dem, start_row, start_col, gen_list):
        self.dem = dem
        self.gen_list = gen_list
        self.start_row = start_row
        self.start_col = start_col

        self.z_delta_array = np.zeros_like(self.dem, dtype=np.float32)
        self.flux_array = np.zeros_like(self.dem, dtype=np.float32)
        self.travel_length_array = np.zeros_like(self.dem, dtype=np.float32)
        self.flow_energy_array = np.zeros_like(self.dem, dtype=np.float32)
        self.generation_array = np.full_like(self.dem, np.nan, dtype=np.float32)

        self.z_delta_generation = []
        self.flux_generation = []
        self.travel_length_generation = []
        self.flow_energy_generation = []

    def get_path_arrays(self):
        for gen, cell_list in enumerate(self.gen_list):
            for cell in cell_list:
                self.z_delta_array[cell.rowindex, cell.colindex] = max(self.z_delta_array[cell.rowindex, cell.colindex], cell.z_delta)
                self.flux_array[cell.rowindex, cell.colindex] = max(self.flux_array[cell.rowindex, cell.colindex], cell.flux)
                self.travel_length_array[cell.rowindex, cell.colindex] = max(self.travel_length_array[cell.rowindex, cell.colindex], cell.min_distance)
                self.flow_energy_array[cell.rowindex, cell.colindex] = max(self.flow_energy_array[cell.rowindex, cell.colindex], cell.flow_energy)
                self.generation_array[cell.rowindex, cell.colindex] = gen

    def get_variables_generation(self):
        for gen, cell_list in enumerate(self.gen_list):
            cell_list_z_delta = []
            cell_list_flux = []
            cell_list_min_distance = []
            cell_list_flow_energy = []

            for cell in cell_list:
                cell_list_z_delta.append(cell.z_delta)
                cell_list_flux.append(cell.flux)
                cell_list_min_distance.append(cell.min_distance)
                cell_list_flow_energy.append(cell.flow_energy)

            self.z_delta_generation.append(cell_list_z_delta)
            self.flux_generation.append(cell_list_flux)
            self.travel_length_generation.append(cell_list_min_distance)
            self.flow_energy_generation.append(cell_list_flow_energy)

    def get_coords_in_genlist_format(self):
        '''
        get coords in the format as the gen_list
        '''
        row = []
        col = []
        for gen in range(len(self.gen_list)):
            row_gen,col_gen = np.where(self.generation_array == gen)
            row.append(list(row_gen))
            col.append(list(col_gen))
        return row,col


    def calc_thalweg_centerof(self, variable, variable_co):
        '''
        variable: variable, which is centered (in format gen_list)
        variable_co: center of variable_so is calculated (variable is weighted) (in format gen_list)
        '''
        co_var = np.zeros(len(self.gen_list))
        for gen in range(1,len(self.gen_list)):
            try:
                var = np.array(variable[gen])
                co = np.array(variable_co[gen])
                variable_co_sum = np.sum(co)
                co_var[gen] = 1 / variable_co_sum * np.sum(var * co)
            except:
                continue

        return co_var

    def plot_test(self):
        self.get_path_arrays()
        self.get_variables_generation()

        row_generation, col_generation = self.get_coords_in_genlist_format()
        energy_coE = self.calc_thalweg_centerof(self.flow_energy_generation, self.flow_energy_generation)
        row_coE = self.calc_thalweg_centerof(row_generation, self.flow_energy_generation)
        col_coE = self.calc_thalweg_centerof(col_generation, self.flow_energy_generation)

        fig, axs = plt.subplots()
        axs.imshow(self.dem, cmap ='Greys', alpha=0.8)
        #axs.contour(self.dem, levels = 10, colors ='k',linewidths=0.5)
        #f = axs.imshow(self.generation_array)
        #fig.colorbar(f, ax = axs, label = 'flow energy')
        axs.scatter(col_coE[::1], row_coE[::1], c = 'r', s = 0.4, label = 'center of energy')
        axs.imshow(self.z_delta_array, cmap = 'Blues', alpha = 0.6)
        fig.savefig(f'/home/paula/data/Flowpy_test/plane/output_1cell_PRA/plots/path_analysis_TEST.png')
        plt.close(fig)

    

    

        