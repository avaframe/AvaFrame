# -*- coding: utf-8 -*-
"""
Created on Wed Apr 28 11:23:58 2021

@author: neuhauser
"""



import h5py
import os
import numpy as np





X = [0, 1, 2, 3]
Y = [0, 1, 2, 3]
Z = [0, 1, 2, 3]

fnout = 'particles.h5'

if os.path.isfile(fnout):
    os.remove(fnout)
with h5py.File(fnout, 'a') as h5f:
    #mesh = h5f.create_group("mesh")
    for i in range(len(X)):
        time = h5f.create_group("{}".format(i))
        time.create_dataset("X", data=X)
        time.create_dataset("Y", data=Y)
        time.create_dataset("Z", data=Z)
    # save arrays in different datasets
    #x = h5f["mesh"]["x"][X]
    #y = h5f.create_dataset("Y", data=Y)
    #z = h5f.create_dataset("Z", data=Z)
h5f.close()