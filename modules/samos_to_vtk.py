# -*- coding: utf-8 -*-
"""
Created on Wed Apr 28 14:31:36 2021

@author: neuhauser
"""

from pyevtk.hl import pointsToVTK
import numpy as np
npoints = 100

for i in range(100):
    x = np.random.rand(npoints)
    y = np.random.rand(npoints)
    z = np.random.rand(npoints)
    pressure = np.random.rand(npoints)
    temp = np.random.rand(npoints)
    pointsToVTK("./points_{}".format(i), x, y, z, data = {"temp" : temp, "pressure" : pressure})