"""

Simple python script to reproduce analytic solution for a Riemann problem,
following the derivations in Faccanoni and Mangeney (2012), Test 2, Case 1.2.
but scaled up in size.

Here the instantanous release of fluid from rest is described using incompressible,
depth-avaeraged mass and momentum conservation equations and a Coulomb-tpye friction law.


"""

import numpy as np
import matplotlib.pyplot as plt
import seaborn
import matplotlib.animation as animation
import os

# local imports
from avaframe.in3Utils import fileHandlerUtils as fU
import avaframe.out3Plot.plotUtils as pU


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

    fig = plt.figure(figsize=(pU.figW, pU.figH))
    plt.title('Dry-Bed, $\delta=21, \phi$=22')
    plt.plot(x, h[:,0], 'k--', label='t_init')
    plt.plot(x, h[:,dtInd], label='t = %.1fs' % dtStep)
    plt.xlabel('x-coordinate [m]')
    plt.ylabel('Flow depth [m]')
    plt.legend()

    fig.savefig(os.path.join(outDir, 'damBreakFlowDepth.%s' % (pU.outputFormat)))

    if cfg['FLAGS'].getboolean('showPlot'):
        plt.show()

    fig = plt.figure(figsize=(pU.figW, pU.figH))
    plt.title('Dry-Bed, $\delta=21, \phi$=22')
    plt.plot(x, u[:,0], 'k--', label='t_init')
    plt.plot(x, u[:,dtInd], label='t = %.1fs' % dtStep)
    plt.xlabel('x-coordinate [m]')
    plt.ylabel('velocity [ms-1]')
    plt.legend()

    fig.savefig(os.path.join(outDir, 'damBreakVelocity.%s' % (pU.outputFormat)))


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


def damBreakSol(avaDir, cfg, cfgC):
    """ Compute flow depth for dam break for granular flow over a dry rough sloping bed with the Savage Hutter model """

    # Set Parameters
    # Coordinate system chosen in the direction of the inclined plane
    g = cfgC['GENERAL'].getfloat('gravAcc')       # acceleration due to gravity [ms-2]
    phi = np.radians(22)                          # slope angle [°]
    delta = np.radians(21)                        # bed friction angle [°]
    gz = g * np.cos(phi)                          # projection of g perpendicular to the inclined plane
    m0 = gz * (np.tan(phi) - np.tan(delta))       # constant x-acceleration resulting from gravity and friction force
    F = - gz * np.tan(delta)                      # Friction force, Coulomb friction
    hL = cfgC['GENERAL'].getfloat('relTh')        # initial height [m] in Riemann problem in state 1 (x<0), hR (x>0)=0
    cL = np.sqrt(gz * hL)                         # wave celeritiy
    dtStep = cfgC['DAMBREAK'].getfloat('dtStep')

    # Define time [0-1] seconds and space [-2,2] meters domains multiplied times 100
    t = np.linspace(0, 20, 2000)
    x = np.linspace(-200, 200, 1000)
    y = np.linspace(0, 1, 1000)
    nt = len(t)
    nx = len(x)
    # Initialise flow depth solution and velocity
    h = np.zeros((nx, nt))
    u = np.zeros((nx, nt))

    # Compute exact solution for case: 'dry bed' - including three different states
    for m in range(nt):
        for k in range(nx):
            cond1 = ((m0*t[m]) / 2.0 - cL) * t[m]
            cond2 = (2.0 *cL + ((m0*t[m]) / 2.0)) * t[m]
            if x[k] <= cond1:
                h[k,m] = hL
            elif cond1 < x[k] <= cond2:
                h[k,m] = ((2.* cL - (x[k] / t[m]) + ((m0 * t[m]) / 2.))**2) / (9. * gz)
            elif x[k] > cond2:
                h[k,m] = 0.0

    # Compute exact solution for velocity
    for m in range(nt):
        for k in range(nx):
            cond1 = ((m0*t[m]) / 2.0 - cL) * t[m]
            cond2 = (2.0 *cL + ((m0*t[m]) / 2.0)) * t[m]
            if x[k] <= cond1:
                u[k,m] = m0 * t[m]
            elif cond1 < x[k] <= cond2:
                u[k,m] = (2./3.) * (cL + (x[k] / t[m]) + m0 * t[m])

    #-----------------------------Plot results --------------
    # Reproduce figure 6, case 1.2.1 - Test 2
    plotResults(x, h, u, dtStep, cfg)

    return hL, h, u, phi
