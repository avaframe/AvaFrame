"""
    Get release area corner coordinates for a rectangle of an area of approx. 100 000 m2

    This file is part of Avaframe.
"""

# Load modules
import numpy as np
import matplotlib.pyplot as plt
from avaframe.in3Utils import generateTopo as gT
import os
import logging
import shutil

# create local logger
# change log level in calling module to DEBUG to see log messages
log = logging.getLogger(__name__)

def makexyPoints(x1, x2, y1, cfgGen):
    """ Make xy Points of release area from given extent and start and end points """

    lenp = int(cfgGen['lenP'])

    # define corner coordinates
    xyPoints = np.zeros((lenp, 2))
    xyPoints[0, 0] = x1 + x2          # xlr
    xyPoints[0, 1] = -0.5 * y1        # ylr
    xyPoints[1, 0] = x1 + x2          # xur
    xyPoints[1, 1] = 0.5 * y1         # yur
    xyPoints[2, 0] = x1               # xul
    xyPoints[2, 1] = 0.5 * y1         # yul
    xyPoints[3, 0] = x1               # xll
    xyPoints[3, 1] = -0.5 * y1        # yll

    # list all points
    for m in range(lenp):
        log.info('Point %d: x %f y %f' % (m, xyPoints[m, 0], xyPoints[m, 1]))

    return xyPoints


def getCorners_FP(xv, yv, dx, cfgFP, cfgGen):

    # Terrain parameters
    hr = float(cfgGen['hr'])
    vol = float(cfgGen['vol'])
    dh = float(cfgGen['dh'])
    xStart = float(cfgGen['xStart'])

    # Compute release area ---------------------------------------------------------------------------------------
    xr = float(cfgFP['xExtent'])
    yr = vol / (xr * dh)

    log.info('volume is: %f' % (yr*xr*dh))

    # make corner coordinates of release area
    xyPoints = makexyPoints(xStart, xr, yr, cfgGen)

    return xv, yv, xyPoints


def getCorners_IP(xv, yv, dx, mean_alpha, cfgGen):

    # Terrain parameters
    hr = float(cfgGen['hr'])
    vol = float(cfgGen['vol'])
    dh = float(cfgGen['dh'])
    xStart = float(cfgGen['xStart'])

    # Compute release area ---------------------------------------------------------------------------------------
    xr = hr / np.sin(np.radians(mean_alpha))        # along valley extent of release area
    xrp = hr / np.tan(np.radians(mean_alpha))       # projected to the flat plane
    yr = vol / (xr * dh)

    log.info('volume is: %f, xr: %f, xrp: %f' % (yr*xr*dh, xr, xrp))

    # make corner coordinates of release area
    xyPoints = makexyPoints(xStart, xrp, yr, cfgGen)

    return xv, yv, xyPoints


def getCorners_HS(xv, yv, dx, mean_alpha, cfgHS, cfgGen, cfgTopo):

    # Terrain parameters
    hr = float(cfgGen['hr'])
    vol = float(cfgGen['vol'])
    dh = float(cfgGen['dh'])
    lenp = int(cfgGen['lenP'])
    alpha_stop = float(cfgHS['alpha_stop'])

    # get A, B from HS
    [A, B, f_len] = gT.getParams(cfgTopo)

    # Compute release area ---------------------------------------------------------------------------------------
    # along valley margin of release area at alpha_stop° point
    xStop = (np.tan(np.radians(-alpha_stop)) - B) / (2. * A)
    xr = hr / np.sin(np.radians(alpha_stop))
    xrp = hr / np.tan(np.radians(alpha_stop))
    yr = vol / (xr * dh)

    xStart = xStop - xrp
    # make corner coordinates of release area
    xyPoints = makexyPoints(xStart, xrp, yr, cfgGen)

    return xv, yv, xyPoints


def correctOrigin(xv, yv, xyPoints, cfgT, y_end):
    """ Move the points so that the correct origin is set """

    xl = float(cfgT['DEMDATA']['xl'])
    yl = float(cfgT['DEMDATA']['yl'])

    xv = xv + xl
    yv = yv + yl + 0.5 * y_end

    xyPoints[:,0] = xyPoints[:,0] + xl
    xyPoints[:,1] = xyPoints[:,1] + yl + 0.5 * y_end

    return xv, yv, xyPoints


def writeNXYZ(xyPoints, DEM_type, cfgFile, cfgGen, outDir):
    """ Write topography information to file """

    lenp = len(xyPoints)
    rel_no = int(cfgFile['rel_no'])
    rel_h = float(cfgGen['dh'])
    rel_name = cfgFile['rel_name']

    # Add vertical coordinate
    z = np.zeros((lenp, 1))
    p_mat = np.matrix(np.append(xyPoints, z, axis=1))

    # Save elevation data to .asc file and add header lines
    with open(os.path.join(outDir, 'release_%d%s.nxyz' % (rel_no, DEM_type)), 'w') as f:
        f.write('name=%s\n' % (rel_name))
        f.write('d0=%.2f\n' % (rel_h))
        f.write('rho=None\n')
        f.write('%d\n' % (lenp))
        for line in p_mat:
            np.savetxt(f, line, fmt='%f')

    # Log info here
    log.info('Release Area written to: %s/release_%d%s.nxyz' % (outDir, rel_no, DEM_type))
    if cfgGen.getboolean('outputtxt'):
        shutil.copyfile(os.path.join(outDir, 'release_%d%s.nxyz' % (rel_no, DEM_type)),
        os.path.join(outDir, 'release_%d%s.txt' % (rel_no, DEM_type)))


def getReleaseArea(cfgT, cfgR, avalancheDir):
    """ Main function to compute release areas """

    cfgTopo = cfgT['TOPO']
    cfgChannel = cfgT['CHANNELS']
    cfgDEM = cfgT['DEMDATA']
    cfgGen = cfgR['GENERAL']
    cfgFP = cfgR['FP']
    cfgHS = cfgR['HS']
    cfgFile = cfgR['FILE']

    # Which DEM type
    DEM_type = cfgTopo['DEM_type']

    log.info('DEM type is set to: %s' % DEM_type)

    # Set Output directory
    outDir = os.path.join(avalancheDir, 'Inputs', 'REL')
    if os.path.isdir(outDir):
        log.info('The new release area is saved to %s' % (outDir))
    else:
        log.error('Required folder structure: NameOfAvalanche/Inputs missing! \
                    Run runInitializeProject first!')

    # determine number of rows and columns to define domain
    dx = float(cfgTopo['dx'])
    x_end = float(cfgTopo['x_end']) + dx
    y_end = float(cfgTopo['y_end']) + dx
    xv = np.arange(0, x_end, dx)
    yv = np.arange(-0.5 * y_end, 0.5 * y_end, dx)
    nrows = int(y_end / dx)                    # number of rows
    ncols = int(x_end / dx)                    # number of columns

    flagCont = False
    # Get release area
    if DEM_type == 'FP':
        [xv, yv, xyPoints] = getCorners_FP(xv, yv, dx, cfgFP, cfgGen)
        flagCont = True

    elif DEM_type == 'IP':
        [xv, yv, xyPoints] = getCorners_IP(
            xv, yv, dx, float(cfgTopo['mean_alpha']), cfgGen)
        flagCont = True

    elif DEM_type == 'HS':
        [xv, yv, xyPoints] = getCorners_HS(
            xv, yv, dx, float(cfgTopo['mean_alpha']), cfgHS, cfgGen, cfgTopo)
        flagCont = True

    elif DEM_type == 'HS2':
        [xv, yv, xyPoints] = getCorners_IP(
            xv, yv, dx, float(cfgTopo['mean_alpha']), cfgGen)
        flagCont = True

    elif DEM_type == 'HX':
        log.info('no release area available for DEM_type: %s' % (DEM_type))

    elif DEM_type == 'BL':
        log.warning('no release area available for DEM_type: %s' % (DEM_type))

    if flagCont:

        # Move to correct correctOrigin
        [xv, yv, xyPoints] = correctOrigin(xv, yv, xyPoints, cfgT, y_end)

        # Write release area
        writeNXYZ(xyPoints, DEM_type, cfgFile, cfgGen, outDir)

    return xv, yv, xyPoints
