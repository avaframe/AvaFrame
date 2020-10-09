"""
    trackedPath
    get information from Outputs/com1DFAPy/particles/particlesxxxx.xxxx.p
    create shape file as Inputs/LINES/pathAB for com2AB
"""

import shapefile
import numpy as np
import matplotlib.pyplot as plt
from pyevtk.hl import pointsToVTK

# Local imports
import avaframe.com1DFAPy.com1DFA as com1DFA
import avaframe.in3Utils.geoTrans as geoTrans
import avaframe.in2Trans.shpConversion as shpConv
import avaframe.in2Trans.ascUtils as IOf
import avaframe.com1DFAPy.DFAtools as DFAtls
from avaframe.out3Plot.plotUtils import *
from avaframe.out3Plot.makePalette import *


def export_to_paraview():

    # inDir = '/home/marie/ava0/AvaFrame/avaframe/data/avaParabola/Outputs/com1DFAPy/particles'
    # inDir = '/home/marie/ava0/AvaFrame/avaframe/data/avaHelixChannel/Outputs/com1DFAPy/particles'
    # inDir = '/home/marie/ava0/AvaFrame/avaframe/data/avaAlr/Outputs/com1DFAPy/particles'
    inDir = "C:/Users/neuhauser/Documents/ParaView/particles_avatest/"
    inDir = "C:/Users/neuhauser/OneDrive - Bundesforschungszentrum fuer Wald/20210225_Literature_AvaFrameMasters/ParaView/avatest/Outputs/com1DFAPy/particles/"
    # inDEM = '/home/marie/ava0/AvaFrame/avaframe/data/avaParabola/Inputs/DEM_PF_Topo.asc'
    # inDEM = '/home/marie/ava0/AvaFrame/avaframe/data/avaHelixChannel/Inputs/DEM_HX_Topo.asc'
    # inDEM = '/home/marie/ava0/AvaFrame/avaframe/data/avaAlr/Inputs/avaAlr.asc'
    inDEM = "C:/Users/neuhauser/Documents/ParaView/Nordkette.asc"

    #pathAB_p = 'C:/Users/neuhauser/Documents/ParaView/Output/pathAB_p'
    #pathAB_m = 'C:/Users/neuhauser/Documents/ParaView/Output/pathAB_m'
    #pathAB_kE = 'C:/Users/neuhauser/Documents/ParaView/Output/pathAB_kE'
    # pathAB = '/home/marie/ava0/AvaFrame/avaframe/data/avaAlr/Inputs/LINES/pathAB'
    # splitPoint = '/home/marie/ava0/AvaFrame/avaframe/data/avaHelixChannel/Inputs/POINTS/splitPoint'
    header = IOf.readASCheader(inDEM)
    dem = IOf.readRaster(inDEM)
    ncols = header.ncols
    nrows = header.nrows
    xllc = header.xllcenter
    yllc = header.yllcenter
    csz = header.cellsize
    xgrid = np.linspace(xllc, xllc+(ncols-1)*csz, ncols)
    ygrid = np.linspace(yllc, yllc+(nrows-1)*csz, nrows)
    PointsX, PointsY = np.meshgrid(xgrid, ygrid)
    XX = PointsX[0, :]
    YY = PointsY[:, 0]
    ZZ = dem['rasterData']

    Particles, TimeStepInfo = com1DFA.readPartFromPickle(inDir, flagAvaDir=False)

    count = 0

    for t in TimeStepInfo:
        particles = Particles[count]
        m = particles['m']
        X = particles['x'] + xllc
        Y = particles['y'] + yllc
        Z = particles['z']
        ux = particles['ux']
        uy = particles['uy']
        uz = particles['uz']
        u = DFAtls.norm(ux, uy, uz)
        U2 = u*u
        Npart = particles['Npart']
        S = particles['s']
        # L = particles['l']

        pointsToVTK("./Output_vtk_Seilbahn_t05/points_{}".format(count), X, Y, Z, data = {"u" : u})

        count = count + 1


if __name__ == "__main__":
    export_to_paraview()

    
