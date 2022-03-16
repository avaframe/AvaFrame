"""

Simple python script to reproduce analytic solution for a Riemann problem,
following the derivations in Faccanoni and Mangeney (2012), Test 2, Case 1.2.
but scaled up in size.

Here the instantanous release of fluid from rest is described using incompressible,
depth-avaeraged mass and momentum conservation equations and a Coulomb-tpye friction law.


"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import os
import logging

# local imports
from avaframe.in3Utils import fileHandlerUtils as fU
import avaframe.out3Plot.plotUtils as pU
import avaframe.in2Trans.ascUtils as IOf
from avaframe.com1DFA import DFAtools

# create local logger
# change log level in calling module to DEBUG to see log messages
log = logging.getLogger(__name__)


def _plotVariable(var, x, dtInd, dtStep, label):

    fig = plt.figure(figsize=(pU.figW, pU.figH))
    plt.title('Dry-Bed')
    plt.plot(x, var[:,0], 'k--', label='t_init')
    plt.plot(x, var[:,dtInd], label='t = %.1fs' % dtStep)
    plt.xlabel('x-coordinate [m]')
    plt.ylabel(label)
    plt.legend()

    return fig


def plotResults(x, h, u, dtStep, cfg):
    """ Create plots of the analytical solution for the given settings,
        including an animation

    """

    # index of time steps
    dtInd = int(100. * dtStep)
    # output directory
    avaDir = cfg['MAIN']['avalancheDir']
    outDir = os.path.join(avaDir, 'Outputs', 'ana1Tests')
    fU.makeADir(outDir)

    fig = _plotVariable(h, x, dtInd, dtStep, 'Flow depth [m]')
    fig.savefig(os.path.join(outDir, 'damBreakFlowDepth.%s' % (pU.outputFormat)))

    if cfg['FLAGS'].getboolean('showPlot'):
        plt.show()
    else:
        plt.close(fig)

    fig = _plotVariable(u, x, dtInd, dtStep, 'Flow velocity [ms-1]')
    fig.savefig(os.path.join(outDir, 'damBreakFlowVelocity.%s' % (pU.outputFormat)))
    if cfg['FLAGS'].getboolean('showPlot'):
        plt.show()
    else:
        plt.close(fig)

    if cfg['FLAGS'].getboolean('showPlot'):
        plt.show()

        # now start visualizing results using animation
        fig, ax = plt.subplots()

        def make_step(step):
            ax.clear()
            ax.plot(x, h[:,0], 'k--', label='t_init')
            ax.plot(x, h[:,step], label='t=%d' % (step))
            plt.title('Animation of dry-bed test')
            ax.set_xlabel('x-coordinate [m]')
            ax.set_ylabel('Flow Depth [m]')
            plt.legend()

        anim = animation.FuncAnimation(fig, make_step, interval=0.1, frames=1000)
        plt.show()


def _plotMultVariables(x, y, nx_loc, dtAnalysis, data1, data2, xR, dataR, tR, label, unit):
    """ generate plots """

    fig, ax = plt.subplots(nrows=1, sharex=True)
    ax.plot(x, y, 'grey', linestyle='--')
    ax.plot(x, data1[nx_loc, :], 'k--', label='init')
    ax.plot(x, data2[nx_loc, :], 'b', label='com1DFAPy')
    ax.plot(xR, dataR[:,tR], 'r-', label='analyt')
    ax.set_xlabel('Along track [ncols]')
    ax.set_ylabel('%s [%s]' % (label, unit))
    plt.legend()
    ax.set_title('%s at time step %.02f s' % (label, dtAnalysis))

    return fig


def plotComparison(cfg, fields0, fieldsT, header, hL, xR, hR, uR, dtAnalysis, cfgMain):
    """ Generate plots that compare the simulation results to the analytical solution

        Parameters
        -----------

        dataComSol: dataFrame
            dataframe of simulation results (including name, file path, result type, etc.)
        hL: float
            initial release thickness
        xR: numpy array
            x extent of domain
        hR: numpy array
            flow depth of analytical solution
        uR: numpy array
            flow velocity of analytical solution
        dtAnalysis: float
            time step of analaysis
        cfgMain: dict
            main configuration for AvaFrame

    """
    phi = cfg['DAMBREAK'].getfloat('phi')
    # Load data
    dataIniFD = fields0['FD']
    dataAnaFD = fieldsT['FD']
    dataIniVx = fields0['Vx']
    dataIniVy = fields0['Vy']
    dataIniVz = fields0['Vz']
    dataAnaVx = fieldsT['Vx']
    dataAnaVy = fieldsT['Vy']
    dataAnaVz = fieldsT['Vz']
    # project velocity on inclined plane
    dataIniV = DFAtools.scalProd(dataIniVx, dataIniVy, dataIniVz, np.cos(np.radians(phi)), 0, -np.sin(np.radians(phi)))
    dataAnaV = DFAtools.scalProd(dataAnaVx, dataAnaVy, dataAnaVz, np.cos(np.radians(phi)), 0, -np.sin(np.radians(phi)))

    # Location of Profiles
    cellSize = header['cellsize']
    ny = dataAnaFD.shape[0]
    nx = dataAnaFD.shape[1]
    xllc = header['xllcenter']
    nx_loc = int(ny *0.5)

    # set x Vector
    x = np.arange(xllc, xllc + nx*cellSize, cellSize)
    y = np.zeros(len(x))
    y[x<0] = hL
    y[x>=0] = 0.0
    y[x<-120] = 0.0

    # setup index for time of analyitcal solution
    tR = int(dtAnalysis * 100.0)

    # setup output directory
    outDir = os.path.join(cfgMain['MAIN']['avalancheDir'], 'Outputs', 'ana1Tests')
    fU.makeADir(outDir)

    fig = _plotMultVariables(x, y, nx_loc, dtAnalysis, dataIniFD, dataAnaFD, xR, hR, tR, 'Flow depth', 'm')
    fig.savefig(os.path.join(outDir, 'CompareDamBreakH.%s' % (pU.outputFormat)))

    y = np.zeros(len(x))
    fig = _plotMultVariables(x, y, nx_loc, dtAnalysis, dataIniV, dataAnaV, xR, uR, tR,
                             'Flow velocity in slope directrion', 'm/s')
    fig.savefig(os.path.join(outDir, 'CompareDamBreakVel.%s' % (pU.outputFormat)))
    plt.show()
    if cfgMain['FLAGS'].getboolean('showPlot'):
        plt.show()
    else:
        plt.close(fig)


def damBreakSol(avaDir, cfg, cfgC):
    """ Compute flow depth for dam break for granular flow over a dry rough sloping bed with the Savage Hutter model

        Parameters
        -----------
        avaDir: str
            path to avalanche directory
        cfg: configParser object
            main avaframe configuration - here showPlot flag used
        cfgC: configParser object
            configuration setting for avalanche simulation including DAMBREAK section

        Returns
        -------
        hL: float
            release thickness
        h: numpy array
            flow depth
        u: numpy array
            flow velocity
        phi: float
            slope angle
        x: numpy array
            extent of domain in the local plane coordinate system
    """

    # Set Parameters
    # Coordinate system chosen in the direction of the inclined plane
    g = cfgC['GENERAL'].getfloat('gravAcc')       # acceleration due to gravity [ms-2]
    phiDeg = cfgC['DAMBREAK'].getfloat('phi')       # acceleration due to gravity [ms-2]
    deltaDeg = cfgC['DAMBREAK'].getfloat('delta')       # acceleration due to gravity [ms-2]
    phi = np.radians(phiDeg)                          # slope angle [°]
    delta = np.radians(deltaDeg)                        # bed friction angle [°]
    gz = g * np.cos(phi)                          # projection of g perpendicular to the inclined plane
    m0 = gz * (np.tan(phi) - np.tan(delta))       # constant x-acceleration resulting from gravity and friction force
    hL = cfgC['GENERAL'].getfloat('relTh')        # initial height [m] in Riemann problem in state 1 (x<0), hR (x>0)=0
    cL = np.sqrt(gz * hL)
    hR = cfgC['GENERAL'].getfloat('relTh')        # initial height [m] in Riemann problem in state 1 (x<0), hR (x>0)=0
    cR = np.sqrt(gz * hR)
    x0R = -120 / np.cos(phi)                           # wave celeritiy
    dtStep = cfgC['DAMBREAK'].getfloat('dtStep')
    # Define time [0-1] seconds and space [-2,2] meters domains multiplied times 100
    t = np.linspace(0, 20, 2000)
    x = np.linspace(-200, 200, 1000)
    nt = len(t)
    nx = len(x)
    # Initialise flow depth solution and velocity
    h = np.zeros((nx, nt))
    u = np.zeros((nx, nt))

    # Compute exact solution for case: 'dry bed' - including three different states
    for m in range(nt):
        for k in range(nx):
            cond1R = x0R + ((m0*t[m]) / 2.0 - 2*cR) * t[m]
            cond2R = x0R + ((m0*t[m]) / 2.0 + cR) * t[m]
            cond1 = ((m0*t[m]) / 2.0 - cL) * t[m]
            cond2 = (2.0 *cL + ((m0*t[m]) / 2.0)) * t[m]
            if x[k] <= cond1R:
                h[k,m] = 0
            elif cond1R < x[k] <= cond2R:
                h[k,m] = ((2.* cR + ((x[k]-x0R) / t[m]) - ((m0 * t[m]) / 2.))**2) / (9. * gz)
            elif cond2R < x[k] <= cond1:
                h[k,m] = hL
            elif cond1 < x[k] <= cond2:
                h[k,m] = ((2.* cL - (x[k] / t[m]) + ((m0 * t[m]) / 2.))**2) / (9. * gz)
            elif x[k] > cond2:
                h[k,m] = 0.0

    # Compute exact solution for velocity
    for m in range(nt):
        for k in range(nx):
            cond1R = x0R + ((m0*t[m]) / 2.0 - 2*cR) * t[m]
            cond2R = x0R + ((m0*t[m]) / 2.0 + cR) * t[m]
            cond1 = ((m0*t[m]) / 2.0 - cL) * t[m]
            cond2 = (2.0 *cL + ((m0*t[m]) / 2.0)) * t[m]
            if x[k] <= cond1R:
                u[k,m] = 0
            elif cond1R < x[k] <= cond2R:
                u[k,m] = (2./3.) * (-cR + ((x[k]-x0R) / t[m]) + m0 * t[m])
            elif cond2R < x[k] <= cond1:
                u[k,m] = m0 * t[m]
            elif cond1 < x[k] <= cond2:
                u[k,m] = (2./3.) * (cL + (x[k] / t[m]) + m0 * t[m])

    #-----------------------------Plot results --------------
    # Reproduce figure 6, case 1.2.1 - Test 2
    plotResults(x, h, u, dtStep, cfg)
    return hL, h, u, phi, x
