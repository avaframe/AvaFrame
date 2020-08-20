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

def makexyPoints(x1, x2, y1, cfgR):
    """ Make xy Points of release area from given extent and start and end points """

    lenp = int(cfgR['GENERAL']['lenP'])

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


def getCorners_FP(cfgR, cfgT):

    # Terrain parameters
    hr = float(cfgR['GENERAL']['hr'])
    vol = float(cfgR['GENERAL']['vol'])
    dh = float(cfgR['GENERAL']['dh'])
    xStart = float(cfgR['GENERAL']['xStart'])

    # Compute release area ---------------------------------------------------------------------------------------
    xr = float(cfgR['FP']['xExtent'])
    yr = vol / (xr * dh)

    log.info('volume is: %f' % (yr*xr*dh))

    # make corner coordinates of release area
    xyPoints = makexyPoints(xStart, xr, yr, cfgR)

    return xyPoints


def getCorners_IP(cfgR, cfgT):

    # Terrain parameters
    hr = float(cfgR['GENERAL']['hr'])
    vol = float(cfgR['GENERAL']['vol'])
    dh = float(cfgR['GENERAL']['dh'])
    xStart = float(cfgR['GENERAL']['xStart'])
    meanAlpha = float(cfgT['TOPO']['meanAlpha'])

    # Compute release area ---------------------------------------------------------------------------------------
    xr = hr / np.sin(np.radians(meanAlpha))        # along valley extent of release area
    xrp = hr / np.tan(np.radians(meanAlpha))       # projected to the flat plane
    yr = vol / (xr * dh)

    log.info('volume is: %f, xr: %f, xrp: %f' % (yr*xr*dh, xr, xrp))

    # make corner coordinates of release area
    xyPoints = makexyPoints(xStart, xrp, yr, cfgR)

    return xyPoints


def getCorners_HS(cfgR, cfgT):

    # Terrain parameters
    hr = float(cfgR['GENERAL']['hr'])
    vol = float(cfgR['GENERAL']['vol'])
    dh = float(cfgR['GENERAL']['dh'])
    lenp = int(cfgR['GENERAL']['lenP'])
    alpha_stop = float(cfgR['HS']['alpha_stop'])

    # get A, B from HS
    [A, B, fLen] = gT.getParabolaParams(cfgT)

    # Compute release area ---------------------------------------------------------------------------------------
    # along valley margin of release area at alpha_stop° point
    xStop = (np.tan(np.radians(-alpha_stop)) - B) / (2. * A)
    xr = hr / np.sin(np.radians(alpha_stop))
    xrp = hr / np.tan(np.radians(alpha_stop))
    yr = vol / (xr * dh)

    xStart = xStop - xrp
    # make corner coordinates of release area
    xyPoints = makexyPoints(xStart, xrp, yr, cfgR)

    return xyPoints


def correctOrigin(xyPoints, cfgT):
    """ Move the points so that the correct origin is set """

    # determine number of rows and columns to define domain
    dx = float(cfgT['TOPO']['dx'])
    xEnd = float(cfgT['TOPO']['xEnd']) + dx
    yEnd = float(cfgT['TOPO']['yEnd']) + dx

    # Compute coordinate grid
    xv = np.arange(0, xEnd, dx)
    yv = np.arange(-0.5 * yEnd, 0.5 * yEnd, dx)

    xl = float(cfgT['DEMDATA']['xl'])
    yl = float(cfgT['DEMDATA']['yl'])

    xv = xv + xl
    yv = yv + yl + 0.5 * yEnd

    xyPoints[:,0] = xyPoints[:,0] + xl
    xyPoints[:,1] = xyPoints[:,1] + yl + 0.5 * yEnd

    return xv, yv, xyPoints


def writeNXYZ(xyPoints, DEM_type, cfgR, outDir):
    """ Write topography information to file """

    lenp = len(xyPoints)
    relNo = int(cfgR['FILE']['relNo'])
    relH = float(cfgR['GENERAL']['dh'])
    relName = cfgR['FILE']['relName']

    # Add vertical coordinate
    z = np.zeros((lenp, 1))
    p_mat = np.matrix(np.append(xyPoints, z, axis=1))

    # Save elevation data to .asc file and add header lines
    with open(os.path.join(outDir, 'release_%d%s.nxyz' % (relNo, DEM_type)), 'w') as f:
        f.write('name=%s\n' % (relName))
        f.write('d0=%.2f\n' % (relH))
        f.write('rho=None\n')
        f.write('%d\n' % (lenp))
        for line in p_mat:
            np.savetxt(f, line, fmt='%f')

    # Log info here
    log.info('Release Area written to: %s/release_%d%s.nxyz' % (outDir, relNo, DEM_type))
    if cfgR.getboolean('GENERAL','outputtxt'):
        shutil.copyfile(os.path.join(outDir, 'release_%d%s.nxyz' % (relNo, DEM_type)),
        os.path.join(outDir, 'release_%d%s.txt' % (relNo, DEM_type)))


def getReleaseArea(cfgT, cfgR, avalancheDir):
    """ Main function to compute release areas """

    # Which DEM type
    DEM_type = cfgT['TOPO']['DEM_type']

    log.info('DEM type is set to: %s' % DEM_type)

    # Set Output directory
    outDir = os.path.join(avalancheDir, 'Inputs', 'REL')
    if os.path.isdir(outDir):
        log.info('The new release area is saved to %s' % (outDir))
    else:
        log.error('Required folder structure: NameOfAvalanche/Inputs missing! \
                    Run runInitializeProject first!')


    flagCont = False
    # Get release area
    if DEM_type == 'FP':
        xyPoints = getCorners_FP(cfgR, cfgT)
        flagCont = True

    elif DEM_type == 'IP':
        xyPoints = getCorners_IP(cfgR, cfgT)
        flagCont = True

    elif DEM_type == 'HS':
        xyPoints = getCorners_HS(cfgR, cfgT)
        flagCont = True

    elif DEM_type == 'HS2':
        xyPoints = getCorners_IP(cfgR, cfgT)
        flagCont = True

    elif DEM_type == 'HX':
        log.info('no release area available for DEM_type: %s' % (DEM_type))

    elif DEM_type == 'BL':
        log.warning('no release area available for DEM_type: %s' % (DEM_type))

    if flagCont:

        # Move to correct correctOrigin
        [xv, yv, xyPoints] = correctOrigin(xyPoints, cfgT)

        # Write release area
        writeNXYZ(xyPoints, DEM_type, cfgR, outDir)

    return xv, yv, xyPoints
