#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 16 15:14:39 2019

@author: Michael Neuhauser
This is the flow class
"""

import numpy as np
import math


class Cell:
    
    def __init__(self, process,  rowindex, colindex, dem_ng, cellsize, susceptibility, z_delta, parent, alpha, exp, startcell):
        '''This class handles the spreading over the DEM!
        Depending on the process different alpha angles are used for energy dissipation.'''
        self.rowindex = rowindex
        self.colindex = colindex
        self.altitude = dem_ng[1, 1]
        self.dem_ng = dem_ng
        self.cellsize = cellsize
        self.tan_beta = np.zeros_like(self.dem_ng)
        self.dist = np.zeros_like(self.dem_ng)
        self.persistence = np.zeros_like(self.dem_ng)
        self.p_fd = np.zeros_like(self.dem_ng)
        self.susceptibility = susceptibility
        self.z_delta = z_delta
        self.alpha = float(alpha)
        self.exp = int(exp)
        self.min_distance = 0
        self.max_distance = 0
        self.min_gamma = 0
        self.max_gamma = 0
        self.sl_gamma = 0

        if process == 'Avalanche':
            self.p_threshold = 3 * 10 ** -4
            self.max_z_delta = 270  # maximum velocity this process can reach
        if process == 'Rockfall':
            self.p_threshold = 3 * 10 ** -4
            self.max_z_delta = 50  # maximum velocity this process can reach
        if process == 'Soil Slides':
            self.p_threshold = 3 * 10 ** -4
            self.max_z_delta = 12  # maximum velocity this process can reach

        if type(startcell) == bool:  # check, if start cell exist (start cell is release point)
            self.is_start = True  # set is_start to True
        else:            
            self.startcell = startcell  # give startcell to cell
            self.is_start = False  # set is_start to False

        self.parent = []
        if type(parent) == Cell:
            self.parent.append(parent)

    def add_os(self, susceptibility):
        self.susceptibility += susceptibility

    def add_parent(self, parent):
        self.parent.append(parent)

    def calc_fp_travelangle(self):
        dist_min = []
        dh = self.startcell.altitude - self.altitude
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
            self.p_fd = self.tan_beta ** self.exp / np.sum(self.tan_beta ** self.exp)

    def calc_persistence(self):
        self.persistence = np.zeros_like(self.dem_ng)
        if self.is_start:
            self.persistence += 1
        elif self.parent[0].is_start:
            self.persistence += 1
        else:
            for parent in self.parent:
                dx = (self.colindex - parent.colindex) + 1 # plus 1 to bring it from range [-1,0,1] to [0,1,2] = index of neighbour array
                dy = (self.rowindex - parent.rowindex) + 1
                maxweight = parent.z_delta
                if dx == 0:
                    if dy == 0:
                        self.persistence[0, 0] += maxweight
                        self.persistence[1, 0] += 0.707 * maxweight
                        self.persistence[0, 1] += 0.707 * maxweight
                    if dy == 1:
                        self.persistence[1, 0] += maxweight
                        self.persistence[2, 0] += 0.707 * maxweight
                        self.persistence[0, 0] += 0.707 * maxweight
                    if dy == 2:
                        self.persistence[2, 0] += maxweight
                        self.persistence[1, 0] += 0.707 * maxweight
                        self.persistence[2, 1] += 0.707 * maxweight

                if dx == 1:
                    if dy == 0:
                        self.persistence[0, 1] += maxweight
                        self.persistence[0, 0] += 0.707 * maxweight
                        self.persistence[0, 2] += 0.707 * maxweight
                    if dy == 2:
                        self.persistence[2, 1] += maxweight
                        self.persistence[2, 0] += 0.707 * maxweight
                        self.persistence[2, 2] += 0.707 * maxweight

                if dx == 2:
                    if dy == 0:
                        self.persistence[0, 2] += maxweight
                        self.persistence[0, 1] += 0.707 * maxweight
                        self.persistence[1, 2] += 0.707 * maxweight
                    if dy == 1:
                        self.persistence[1, 2] += maxweight
                        self.persistence[0, 2] += 0.707 * maxweight
                        self.persistence[2, 2] += 0.707 * maxweight
                    if dy == 2:
                        self.persistence[2, 2] += maxweight
                        self.persistence[2, 1] += 0.707 * maxweight
                        self.persistence[1, 2] += 0.707 * maxweight
                    
    def calc_distribution(self):

        self.calc_z_delta()
        self.calc_persistence()
        self.calc_tanbeta()

        if not self.is_start:
            self.calc_fp_travelangle()
            self.calc_sl_travelangle()

        threshold = self.p_threshold
        if np.sum(self.p_fd) > 0:
            self.dist = (self.persistence * self.p_fd) / np.sum(self.persistence * self.p_fd) * self.susceptibility
        # This lines handle if a distribution to a neighbour cell is lower then the threshold, so we don´t lose
        # susceptibility.
        # The susceptibility of this cells will then spread equally to all neighbour cells
        count = ((0 < self.dist) & (self.dist < threshold)).sum()
        mass_to_distribute = np.sum(self.dist[self.dist < threshold])
        '''Checking if susceptibility is distributed to a field that isn't taking in account, when then distribute it to
         the other fields'''
        if mass_to_distribute > 0 and count > 0:
            self.dist[self.dist > threshold] += mass_to_distribute / count
            self.dist[self.dist < threshold] = 0
        if np.sum(self.dist) < self.susceptibility and count > 0:
            self.dist[self.dist > threshold] += (self.susceptibility - np.sum(self.dist))/count

        row_local, col_local = np.where(self.dist > threshold)

        return self.rowindex - 1 + row_local, self.colindex - 1 + col_local, self.dist[row_local, col_local], self.z_delta_neighbour[row_local, col_local]
