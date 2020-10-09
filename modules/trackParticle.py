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


def writeLine2SHPfile(part, lineName, fileName):
    """copied from
    https://pypi.org/project/pyshp/#writing-shapefiles
    section: Adding a LineString shape
    """

    w = shapefile.Writer(fileName)
    w.field('name', 'C')
    w.line([part])
    w.record(lineName)
    w.close()


def trackParticle():

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
    
    x = []
    y = []
    z = []
    vel = []
    
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

        x.append(X[0])
        y.append(Y[0])
        z.append(Z[0])
        vel.append(u[0])
    
    x = np.array(x-x[0])
    y = np.array(y-y[0])
    h = np.sqrt(x**2 + y**2)    
    #plt.plot(h, z)
    fig, ax = plt.subplots()
    ax.plot(h, z, label='trajectory', color='blue')
    ax.set_xlabel('Horizontal Distance [m]')
    ax.set_ylabel('Altitude [m]')
    #ax.legend()
    ax2 = ax.twinx()
    ax2.scatter(h, vel, label='Velocity', color='red')
    ax2.set_ylabel('Velocity [m/s]')
    #ax2.legend()
    ax.set_title('trajectory and velocity of one particle')
    # ask matplotlib for the plotted objects and their labels
    lines, labels = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax2.legend(lines + lines2, labels + labels2, loc=0)
    
        
    
    
def write_hdf5(fnout, X, Y, Z):
    # save raw data to hdf5
    if os.path.isfile(fnout):
        os.remove(fnout)
    with h5py.File(fnout, 'a') as h5f:
        # save arrays in different datasets
        x = h5f.create_dataset("X", data=X)
        y = h5f.create_dataset("Y", data=Y)
        z = h5f.create_dataset("Z", data=Z)
        
def write_to_vtk(fnout):
    from pyevtk.hl import pointsToVTK
    import numpy as np
    npoints = 100
    x = np.random.rand(npoints)
    y = np.random.rand(npoints)
    z = np.random.rand(npoints)
    pressure = np.random.rand(npoints)
    temp = np.random.rand(npoints)
    pointsToVTK("./points", x, y, z, data = {"temp" : temp, "pressure" : pressure})
        


if __name__ == "__main__":
    trackParticle()
    #write_hdf5('particles.h5', x, y, z)
    
