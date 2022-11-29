import glob
import os
from paraview.simple import *
import shutil

hcont = []
hmin = 1
hmax = 100
hint = 2
for i in range(hmin, hmax, hint): hcont.append(i)

elevcont = []
elevmin = -11000
elevmax = 9000
elevint = 100
for i in range(elevmin, elevmax, elevint): elevcont.append(i)

if os.path.exists('surface'): shutil.rmtree('surface')
if os.path.exists('contoursh'): shutil.rmtree('contoursh')
if os.path.exists('contoursz'): shutil.rmtree('contoursz')

intab = paraview.simple.CSVReader(FileName=glob.glob('data/pv*.csv'))
points = paraview.simple.TableToPoints(intab, XColumn='x', YColumn='y', ZColumn='z', KeepAllDataArrays=True)
surface = paraview.simple.Delaunay2D(points)
surfcalc = paraview.simple.Calculator(Input=surface, ResultArrayName='rgb', Function='r * iHat + g * jHat + b * kHat')
contoursh = paraview.simple.Contour(Input=surface, ContourBy='h', Isosurfaces=hcont)
contoursz = paraview.simple.Contour(Input=surface, ContourBy='z', Isosurfaces=elevcont)

print('Writing surfaces ...')
paraview.simple.SaveData('surface.pvd', surfcalc, WriteTimeSteps=True)

print('Writing flow contours ...')
paraview.simple.SaveData('contoursh.pvd', contoursh, WriteTimeSteps=True)

print('Writing surface contours ...')
paraview.simple.SaveData('contoursz.pvd', contoursz, WriteTimeSteps=True)

print ('completed.')
