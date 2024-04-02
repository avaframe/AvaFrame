#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import math

class Cell:
    """ This is the com4FlowPy 'Cell ' class
        This class handles the calculation at the 'cell level' (vgl. D'Amboise et al., 2022)
    """
    def __init__(self, rowindex, colindex, dem_ng, cellsize, flux, 
                 z_delta, parent, alpha, exp, flux_threshold, 
                 max_z_delta, startcell):
        """ constructor for the Cell class
            the constructor function is called every time a new instance of type 'Cell' is
            initialized.
        """
        self.rowindex = rowindex #index of the Cell in row-direction (i.e. local y-index in the calculation domain)
        self.colindex = colindex #index of the Cell in column-direction (i.e. local x-index in the calculation domain)
        self.dem_ng = dem_ng #elevation values in the 3x3 neigbourhood around the Cell
        self.altitude = dem_ng[1, 1] #elevation value of the cell (central cell of 3x3 neighbourhood)
        self.cellsize = cellsize #cellsize in meters
        
        self.tan_beta = np.zeros_like(self.dem_ng)
        self.dist = np.zeros_like(self.dem_ng)
        self.persistence = np.zeros_like(self.dem_ng)
        self.r_t = np.zeros_like(self.dem_ng)
        self.no_flow = np.ones_like(self.dem_ng)

        self.flux = flux
        self.z_delta = z_delta

        self.alpha = float(alpha)
        self.exp = int(exp)
        self.max_z_delta = float(max_z_delta)
        self.flux_threshold = float(flux_threshold)
        
        self.min_distance = 0 #minimal distance to start-cell (i.e. along shortest path) min_distance >= 
        self.max_distance = 0 #NOTE: self.max_distance is never used - maybe remove!?
        self.min_gamma = 0  #NOTE: self.min_gamma (assumingly minimal travel angle to cell) is never used - maybe remove!?
        self.max_gamma = 0
        self.sl_gamma = 0             

        if type(startcell) == bool:     # if a boolean variable (i.e.'True') is passed to the constructor
            self.is_start = True        # set is_start to True
        else:            
            self.startcell = startcell  # give startcell to cell
            self.is_start = False       # set is_start to False

        self.parent = []                #NOTE: maybe renameto 'parents' or 'lOfParents' for clarity ...
        if type(parent) == Cell:
            self.parent.append(parent)

    def add_os(self, flux):
        self.flux += flux

    def add_parent(self, parent):
        self.parent.append(parent)

    def calc_fp_travelangle(self):
        """ function calculates the travel-angle along the shortest flow-path from the start-cell to the current cell
            the trave-angle along the shortest flow-path is equivalent to the maximum travel angle along all paths from
            the startcell to this cell.
        """
        dist_min = [] #
        dh = self.startcell.altitude - self.altitude #elevation difference from cell to start-cell
        for parent in self.parent:
            dx = abs(parent.colindex - self.colindex)
            dy = abs(parent.rowindex - self.rowindex)
            dist_min.append(math.sqrt(dx ** 2 + dy ** 2) * self.cellsize + parent.min_distance)
        self.min_distance = np.amin(dist_min)
        self.max_gamma = np.rad2deg(np.arctan(dh / self.min_distance))

    def calc_sl_travelangle(self):
        dx = abs(self.startcell.colindex - self.colindex)
        dy = abs(self.startcell.rowindex - self.rowindex)
        dh = self.startcell.altitude - self.altitude

        ds = math.sqrt(dx ** 2 + dy ** 2) * self.cellsize
        self.sl_gamma = np.rad2deg(np.arctan(dh / ds))

    def calc_z_delta(self):
        self.z_delta_neighbour = np.zeros((3, 3))
        self.z_gamma = self.altitude - self.dem_ng
        ds = np.array([[np.sqrt(2), 1, np.sqrt(2)], [1, 0, 1], [np.sqrt(2), 1, np.sqrt(2)]])
        tan_alpha = np.tan(np.deg2rad(self.alpha))
        self.z_alpha = ds * self.cellsize * tan_alpha
        self.z_delta_neighbour = self.z_delta + self.z_gamma - self.z_alpha
        self.z_delta_neighbour[self.z_delta_neighbour < 0] = 0
        self.z_delta_neighbour[self.z_delta_neighbour > self.max_z_delta] = self.max_z_delta
           
    def calc_tanbeta(self):
        ds = np.array([[np.sqrt(2), 1, np.sqrt(2)], [1, 1, 1], [np.sqrt(2), 1, np.sqrt(2)]])
        distance = ds * self.cellsize
        
        beta = np.arctan((self.altitude - self.dem_ng) / distance) + np.deg2rad(90)
        self.tan_beta = np.tan(beta/2)

        self.tan_beta[self.z_delta_neighbour <= 0] = 0
        self.tan_beta[self.persistence <= 0] = 0
        self.tan_beta[1, 1] = 0
        if abs(np.sum(self.tan_beta)) > 0:
            self.r_t = self.tan_beta ** self.exp / np.sum(self.tan_beta ** self.exp)

    def calc_persistence(self):
        self.persistence = np.zeros_like(self.dem_ng)
        if self.is_start:
            self.persistence += 1
        elif self.parent[0].is_start:
            self.persistence += 1
        else:
            for parent in self.parent:
                dx = (parent.colindex - self.colindex) 
                dy = (parent.rowindex - self.rowindex)

                self.no_flow[dy + 1,dx + 1] = 0  # 3x3 Matrix of ones, every parent gets a 0, so no flow to a parent field.
                
                maxweight = parent.z_delta
                # Old Calculation
                if dx == -1:
                    if dy == -1:
                        self.persistence[2, 2] += maxweight
                        self.persistence[2, 1] += 0.707 * maxweight
                        self.persistence[1, 2] += 0.707 * maxweight
                    if dy == 0:
                        self.persistence[1, 2] += maxweight
                        self.persistence[2, 2] += 0.707 * maxweight
                        self.persistence[0, 2] += 0.707 * maxweight
                    if dy == 1:
                        self.persistence[0, 2] += maxweight
                        self.persistence[0, 1] += 0.707 * maxweight
                        self.persistence[1, 2] += 0.707 * maxweight

                if dx == 0:
                    if dy == -1:
                        self.persistence[2, 1] += maxweight
                        self.persistence[2, 0] += 0.707 * maxweight
                        self.persistence[2, 2] += 0.707 * maxweight
                    if dy == 1:
                        self.persistence[0, 1] += maxweight
                        self.persistence[0, 0] += 0.707 * maxweight
                        self.persistence[0, 2] += 0.707 * maxweight

                if dx == 1:
                    if dy == -1:
                        self.persistence[2, 0] += maxweight
                        self.persistence[1, 0] += 0.707 * maxweight
                        self.persistence[2, 1] += 0.707 * maxweight
                    if dy == 0:
                        self.persistence[1, 0] += maxweight
                        self.persistence[0, 0] += 0.707 * maxweight
                        self.persistence[2, 0] += 0.707 * maxweight
                    if dy == 1:
                        self.persistence[0, 0] += maxweight
                        self.persistence[0, 1] += 0.707 * maxweight
                        self.persistence[1, 0] += 0.707 * maxweight
                        

                    
    def calc_distribution(self):

        self.calc_z_delta()
        self.calc_persistence()
        self.persistence *= self.no_flow
        self.calc_tanbeta()
        #print(self.persistence)

        if not self.is_start:
            self.calc_fp_travelangle()
            self.calc_sl_travelangle()

        threshold = self.flux_threshold
        if np.sum(self.r_t) > 0:
            self.dist = (self.persistence * self.r_t) / np.sum(self.persistence * self.r_t) * self.flux
        # This lines handle if a distribution to a neighbour cell is lower then the threshold, so we don´t lose
        # flux.
        # The flux of this cells will then spread equally to all neighbour cells
        count = ((0 < self.dist) & (self.dist < threshold)).sum()
        mass_to_distribute = np.sum(self.dist[self.dist < threshold])
        '''Checking if flux is distributed to a field that isn't taking in account, when then distribute it equally to
         the other fields'''
        if mass_to_distribute > 0 and count > 0:
            self.dist[self.dist > threshold] += mass_to_distribute / count
            self.dist[self.dist < threshold] = 0
        if np.sum(self.dist) < self.flux and count > 0:
            self.dist[self.dist > threshold] += (self.flux - np.sum(self.dist))/count

        row_local, col_local = np.where(self.dist > threshold)

        return self.rowindex - 1 + row_local, self.colindex - 1 + col_local, self.dist[row_local, col_local], self.z_delta_neighbour[row_local, col_local]
