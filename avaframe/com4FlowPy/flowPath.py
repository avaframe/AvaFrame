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
        self.row_generation = []
        self.col_generation = []
        self.altitude_generation = []

    def get_path_arrays(self):
        ''' 
            write ARRAYS with size of dem containing the variable values of every path
            value 0 means, the path does not hit the cell
        '''
        for gen, cell_list in enumerate(self.gen_list):
            for cell in cell_list:
                self.z_delta_array[cell.rowindex, cell.colindex] = max(self.z_delta_array[cell.rowindex, cell.colindex], cell.z_delta)
                self.flux_array[cell.rowindex, cell.colindex] = max(self.flux_array[cell.rowindex, cell.colindex], cell.flux)
                self.travel_length_array[cell.rowindex, cell.colindex] = max(self.travel_length_array[cell.rowindex, cell.colindex], cell.min_distance)
                self.flow_energy_array[cell.rowindex, cell.colindex] = max(self.flow_energy_array[cell.rowindex, cell.colindex], cell.flow_energy)
                self.generation_array[cell.rowindex, cell.colindex] = gen

    def get_variables_generation(self):
        '''
            write LISTS with size and format of gen_list 
            (the main list contains lists for every generation)
        '''
        for gen, cell_list in enumerate(self.gen_list):
            cell_list_z_delta = []
            cell_list_flux = []
            cell_list_min_distance = []
            cell_list_flow_energy = []
            cell_list_row = []
            cell_list_col = []
            cell_list_alt = []

            for cell in cell_list:
                cell_list_z_delta.append(cell.z_delta)
                cell_list_flux.append(cell.flux)
                cell_list_min_distance.append(cell.min_distance)
                cell_list_flow_energy.append(cell.flow_energy)
                cell_list_row.append(cell.rowindex)
                cell_list_col.append(cell.colindex)
                cell_list_alt.append(cell.altitude)

            self.z_delta_generation.append(cell_list_z_delta)
            self.flux_generation.append(cell_list_flux)
            self.travel_length_generation.append(cell_list_min_distance)
            self.flow_energy_generation.append(cell_list_flow_energy)
            self.row_generation.append(cell_list_row)
            self.col_generation.append(cell_list_col)
            self.altitude_generation.append(cell_list_alt)

    def get_coords_in_genlist_format(self):
        '''
        get coords in the format as the gen_list
        UNNÖTIG????
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
        variable_sum = np.zeros(len(self.gen_list))
        for gen in range(1,len(self.gen_list)):
            var = np.array(variable[gen])
            co = np.array(variable_co[gen])
            variable_sum[gen] = np.sum(var)
            variable_co_sum = np.sum(co)
            co_var[gen] = 1 / variable_co_sum * np.sum(var * co)
        return variable_sum, co_var

    def get_centerofs(self):
        '''
        calculate sum of variable for every iteration step/ generation and
        center of energy & flux for the following variables:
        '''

        self.get_path_arrays()
        self.get_variables_generation()

        variables = {'col': self.col_generation, 'row': self.row_generation, 'flux':self.flux_generation, 
        'energy':self.flow_energy_generation, 'altitude': self.altitude_generation, 
        's' : self.travel_length_generation, 'z_delta': self.z_delta_generation}
        for var_name, var in variables.items():
            sumF, coF = self.calc_thalweg_centerof(var, self.flux_generation) # center of flux of every variable
            sumE, coE = self.calc_thalweg_centerof(var, self.flow_energy_generation) # center of energy of every variable

            setattr(self, f'{var_name}_sumF', sumF)
            setattr(self, f'{var_name}_coF', coF)
            setattr(self, f'{var_name}_sumE', sumE)
            setattr(self, f'{var_name}_coE', coE)


    def plot(self):
        fig, axs = plt.subplots()
        axs.imshow(self.dem, cmap ='Greys', alpha=0.8)
        axs.contour(self.dem, levels = 10, colors ='k',linewidths=0.5)
        #fig.colorbar(f, ax = axs, label = 'flow energy')
        axs.scatter(self.col_coE, self.row_coE, c = 'r', s = 0.4, label = 'center of energy')
        fig.savefig(f'/home/paula/data/Flowpy_test/plane/output_1cell_PRA/plots/path_analysis_TEST.png')
        plt.close(fig)

    

    

        