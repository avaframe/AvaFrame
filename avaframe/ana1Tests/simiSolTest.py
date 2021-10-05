"""
Similarity solution module

This module contains functions that compute the similarity solution
for a gliding avalanche on a inclined plane according to similarity solution from :
Hutter, K., Siegel, M., Savage, S.B. et al.
Two-dimensional spreading of a granular avalanche down an inclined plane
Part I. theory. Acta Mechanica 100, 37–68 (1993).
https://doi.org/10.1007/BF01176861
"""

# imports
import numpy as np
from scipy.integrate import ode
import math
import pathlib
import logging

# local imports
from avaframe.in3Utils import cfgUtils
from avaframe.in1Data import getInput as gI
import avaframe.in3Utils.geoTrans as geoTrans
import avaframe.com1DFA.com1DFA as com1DFA
import avaframe.com1DFA.DFAtools as DFAtls
import avaframe.in2Trans.ascUtils as IOf
import avaframe.ana1Tests.analysisTools as anaTools
import avaframe.out3Plot.outAna1Plots as outAna1Plots


# create local logger
# change log level in calling module to DEBUG to see log messages
log = logging.getLogger(__name__)


def mainCompareSimSolCom1DFA(avalancheDir, cfgMain, simiSolCfg, outDirTest):
    """ Compute com1DFA sol, similarity solution and compare
    Main script for the comparison of DFA simulations to the analytic similarity solution
    Parameters
    -----------
    avalancheDir: str or pathlib path
        avalanche directory
    cfgMain: confiparser
        avaframeCfg configuration
    simiSolCfg: pathlib path
        path to simiSol configuration file
    outDirTest: pathlib path
        path to output directory

    """

    cfg = cfgUtils.getModuleConfig(com1DFA, simiSolCfg)
    # Define release thickness distribution
    demFile = gI.getDEMPath(avalancheDir)
    relDict = getReleaseThickness(avalancheDir, cfg, demFile)
    relTh = relDict['relTh']

    # call com1DFA to perform simulations - provide configuration file and release thickness function
    # (may be multiple sims)
    _, _, Tsave, dem, _, _ = com1DFA.com1DFAMain(avalancheDir, cfgMain, cfgFile=simiSolCfg, relThField=relTh)
    relDict['dem'] = dem

    # compute the similartiy solution (this corresponds to our reference)
    log.info('Computing similarity solution')
    solSimi = mainSimilaritySol()

    # now compare the simulations to the reference
    # first fetch info about all the simulations performed (and maybe order them)
    varParList = cfg['ANALYSIS']['varParList'].split('|')
    ascendingOrder = cfg['ANALYSIS']['ascendingOrder']
    # load info for all configurations and order them
    simDF = cfgUtils.orderSimFiles(avalancheDir, '', varParList, ascendingOrder)
    # loop on all the simulations and make the comparison to reference
    for index, row in simDF.iterrows():
        simName = row['simName']
        # fetch the simulation results
        particlesList, _ = com1DFA.readPartFromPickle(avalancheDir, simName=simName, flagAvaDir=True, comModule='com1DFA')
        fieldsList, fieldHeader = com1DFA.readFields(avalancheDir, ['FD', 'FV', 'Vx', 'Vy', 'Vz'], simName=simName, flagAvaDir=True, comModule='com1DFA')
        simDF.loc[index, 'Npart'] = particlesList[0]['Npart']
        # analyze and compare results
        hEL2Array, hELMaxArray, vEL2Array, vELMaxArray = analyzeResults(particlesList, fieldsList, solSimi, relDict,
                                                                        cfg, outDirTest, index)
        # add result of error analysis
        # save results in the simDF
        simDF.loc[index, 'hErrorL2'] = hEL2Array[-1]
        simDF.loc[index, 'vErrorL2'] = vEL2Array[-1]
        simDF.loc[index, 'hErrorLMax'] = hELMaxArray[-1]
        simDF.loc[index, 'vErrorLMax'] = vELMaxArray[-1]
        # +++++++++POSTPROCESS++++++++++++++++++++++++
        # -------------------------------
        if cfgMain['FLAGS'].getboolean('showPlot'):
            outAna1Plots.plotContoursSimiSol(particlesList, fieldsList, solSimi, relDict, cfg, outDirTest)

        outAna1Plots.showSaveTimeSteps(cfgMain, cfg, particlesList, fieldsList, solSimi, Tsave, relDict, outDirTest, index)

    simDF.to_pickle(outDirTest / 'results.p')
    outAna1Plots.plotError(simDF, outDirTest)


def mainSimilaritySol():
    """ Compute similarity solution"""

    # Load configuration
    simiSolCfg = pathlib.Path('data/avaSimilaritySol', 'Inputs', 'simiSol_com1DFACfg.ini')
    cfg = cfgUtils.getModuleConfig(com1DFA, simiSolCfg)
    cfgGen = cfg['GENERAL']
    cfgSimi = cfg['SIMISOL']
    bedFrictionAngleDeg = cfgSimi.getfloat('bedFrictionAngle')
    planeinclinationAngleDeg = cfgSimi.getfloat('planeinclinationAngle')
    internalFrictionAngleDeg = cfgSimi.getfloat('internalFrictionAngle')
    # Dimensioning parameters L
    L_x = cfgSimi.getfloat('L_x')
    L_y = cfgSimi.getfloat('L_y')
    H = cfgGen.getfloat('relTh')

    # Set parameters
    Pi = math.pi
    gravAcc = cfgGen.getfloat('gravAcc')
    zeta = planeinclinationAngleDeg * Pi /180       # plane inclination
    delta = bedFrictionAngleDeg * Pi /180           # basal angle of friction
    phi = internalFrictionAngleDeg * Pi /180        # internal angle of friction phi>delta

    # Dimensioning parameters
    U = np.sqrt(gravAcc*L_x)
    V = np.sqrt(gravAcc*L_y)
    T = np.sqrt(L_x/gravAcc)

    # calculate aspect ratios
    eps_x = H/L_x
    eps_y = H/L_y
    eps_xy = L_y/L_x

    # Full scale end time
    T_end = cfgGen.getfloat('tEnd') + cfgGen.getfloat('maxdT')

    # Non dimensional time for similarity sim calculation
    t_1 = 0.1/T         # start time for ode solvers
    t_end = T_end/T     # end time
    dt_early = 0.01/T   # time step for early sol
    dt = 0.01/T         # time step for early sol

    # Initial conditions [g0 g_p0 f0 f_p0]
    x_0 = [1.0, 0.0, 1.0, 0.0]  # here a circle as start point

    # compute earth pressure coefficients
    flagEarth = cfgSimi.getboolean('flagEarth')
    if flagEarth:
        earthPressureCoefficients = defineEarthPressCoeff(phi, delta)
    else:
        earthPressureCoefficients = np.ones((6, 1))

    # Early time solution
    t_early = np.arange(0, t_1, dt_early)
    solSimi = calcEarlySol(t_early, earthPressureCoefficients, x_0, zeta, delta, eps_x, eps_xy, eps_y)

    # Runge-Kutta integration away from the singularity
    # initial conditions
    t_start = t_1 - dt_early
    x_1 = np.empty((4, 1))
    x_1[0] = solSimi['g_sol'][-1]
    x_1[1] = solSimi['g_p_sol'][-1]
    x_1[2] = solSimi['f_sol'][-1]
    x_1[3] = solSimi['f_p_sol'][-1]

    # Create an `ode` instance to solve the system of differential
    # equations defined by `fun`, and set the solver method to'dopri5' 'dop853'.
    solver = ode(Ffunction)
    solver.set_integrator('dopri5')
    solver.set_f_params(earthPressureCoefficients, zeta, delta, eps_x, eps_xy, eps_y)
    solver.set_initial_value(x_1, t_start)
    solSimi = odeSolver(solver, dt, t_end, solSimi)

    Time = solSimi['time']*T
    solSimi['Time'] = Time

    return solSimi

#####################################
# Compute similarity solution
#####################################


def defineEarthPressCoeff(phi, delta):
    """ Define earth pressure coefficients

        Parameters
        -----------
        phi: float
            internal friction angle
        delta: float
            bed friction angle

        Returns
        --------
        earthPressureCoefficients: numpy array
            [Kxact Kxpass Kyact(Kxact) Kypass(Kxact) Kyact(Kxpass) Kypass(Kxpass)]
    """

    cos2phi = np.cos(phi)**2
    cos2delta = np.cos(delta)**2
    tan2delta = np.tan(delta)**2
    root1 = np.sqrt(1.0 - cos2phi / cos2delta)
    earthPressureCoefficients = np.zeros((6, 1))
    # kx active / passive
    earthPressureCoefficients[0] = 2 / cos2phi * (1.0 - root1) - 1.0
    earthPressureCoefficients[1] = 2 / cos2phi * (1.0 + root1) - 1.0
    # ky active / passive for kx active
    kx = earthPressureCoefficients[0]
    root2 = np.sqrt((1.0 - kx) * (1.0 - kx) + 4.0 * tan2delta)
    earthPressureCoefficients[2] = 0.5 * (kx + 1.0 - root2)
    earthPressureCoefficients[3] = 0.5 * (kx + 1.0 + root2)
    # ky active / passive for kx passive
    kx = earthPressureCoefficients[1]
    root2 = np.sqrt((1.0 - kx) * (1.0 - kx) + 4.0 * tan2delta)
    earthPressureCoefficients[4] = 0.5 * (kx + 1.0 - root2)
    earthPressureCoefficients[5] = 0.5 * (kx + 1.0 + root2)

    return earthPressureCoefficients


def computeEarthPressCoeff(x, earthPressureCoefficients):
    """ Computes the earth pressure coefficients function of sng of f and g
        i.e depending on if we are in the active or passive case

        Parameters
        -----------
        x: float
            internal friction angle
        earthPressureCoefficients: numpy array
            [Kxact Kxpass Kyact(Kxact) Kypass(Kxact) Kyact(Kxpass) Kypass(Kxpass)]

        Returns
        --------
        K_x: float
            earth pressure coefficient in x direction
        K_y: float
            earth pressure coefficient in y direction
    """

    g_p = x[1]
    f_p = x[3]
    if g_p >= 0:
        K_x = earthPressureCoefficients[0]
        if f_p >= 0:
            K_y = earthPressureCoefficients[2]
        else:
            log.info('ky passive')
            K_y = earthPressureCoefficients[3]
    else:
        log.info('kx passive')
        K_x = earthPressureCoefficients[1]
        if f_p >= 0:
            K_y = earthPressureCoefficients[4]
        else:
            log.info('ky passive')
            K_y = earthPressureCoefficients[5]

    return K_x, K_y


def computeFCoeff(x, K_x, K_y, zeta, delta, eps_x, eps_xy, eps_y):
    """ Compute coefficients eq 3.2 for the function F
    """

    A = np.sin(zeta)
    B = eps_x * np.cos(zeta) * K_x
    C = np.cos(zeta) * np.tan(delta)
    D = eps_y / eps_xy * np.cos(zeta) * K_y
    if A == 0:
        E = 1
        C = 0
    else:
        E = (A-C)/A
        C = np.cos(zeta) * np.tan(delta)

    return A, B, C, D, E


def calcEarlySol(t, earthPressureCoefficients, x_0, zeta, delta, eps_x, eps_xy, eps_y):
    """ Compute the early solution for 0<t<t_1 to avoid singularity in the
        Runge-Kutta integration process
    """

    # early solution exists only if first derivative of f at t=0 is zero
    assert x_0[3] == 0, "f'(t=0)=f_p0 must be equal to 0"
    K_x, K_y = computeEarthPressCoeff(x_0, earthPressureCoefficients)
    A, B, C, D, E = computeFCoeff(x_0, K_x, K_y, zeta, delta, eps_x, eps_xy, eps_y)
    g0 = x_0[0]
    g_p0 = x_0[1]
    f0 = x_0[2]
    f_p0 = x_0[3]

    g = g0 + g_p0*t + B/(f0*g0**2)*t**2
    g_p = g_p0 + 2*B/(f0*g0**2)*t
    f = f0 + f_p0*t + D*E/(g0*f0**2)*t**2
    f_p = f_p0 + 2*D*E/(g0*f0**2)*t

    solSimi = {}
    solSimi['time'] = t
    solSimi['g_sol'] = g
    solSimi['g_p_sol'] = g_p
    solSimi['f_sol'] = f
    solSimi['f_p_sol'] = f_p

    return solSimi


def Ffunction(t, x, earthPressureCoefficients, zeta, delta, eps_x, eps_xy, eps_y):
    """ Calculate right hand side of the differential equation :
        dx/dt = F(x,t) F is discribed in Hutter 1993.

        Parameters
        -----------
        t: float
            curent time
        x: numpy array
            column vector of size 4

        Returns:
        F: numpy array
            column vector of size 4
    """

    global A, C
    K_x, K_y = computeEarthPressCoeff(x, earthPressureCoefficients)
    A, B, C, D, E = computeFCoeff(x, K_x, K_y, zeta, delta, eps_x, eps_xy, eps_y)
    u_c = (A - C)*t
    g = x[0]
    g_p = x[1]
    f = x[2]
    f_p = x[3]

    dx0 = g_p
    dx1 = 2*B/(g**2*f)
    dx2 = f_p
    if C == 0:
        dx3 = 2*D/(g*f**2)
    else:
        dx3 = 2*D/(g*f**2)-C*f_p/u_c
    F = [dx0, dx1, dx2, dx3]

    return F


def odeSolver(solver, dt, t_end, solSimi):
    """ Solve the ODE using a Runge-Kutta method
    """
    time = solSimi['time']
    g_sol = solSimi['g_sol']
    g_p_sol = solSimi['g_p_sol']
    f_sol = solSimi['f_sol']
    f_p_sol = solSimi['f_p_sol']

    while solver.successful() and solver.t < t_end:
        solver.integrate(solver.t+dt, step=True)
        x_sol = solver.y
        t_sol = solver.t
        time = np.append(time, t_sol)
        g_sol = np.append(g_sol, x_sol[0])
        g_p_sol = np.append(g_p_sol, x_sol[1])
        f_sol = np.append(f_sol, x_sol[2])
        f_p_sol = np.append(f_p_sol, x_sol[3])

    solSimi['time'] = time
    solSimi['g_sol'] = g_sol
    solSimi['g_p_sol'] = g_p_sol
    solSimi['f_sol'] = f_sol
    solSimi['f_p_sol'] = f_p_sol

    return solSimi


def h(solSimi, x1, y1, i, L_y, L_x, H):
    """ get flow depth from f and g solutions
    """
    time = solSimi['time']
    g_sol = solSimi['g_sol']
    f_sol = solSimi['f_sol']
    y1 = -(y1/L_y)**2/(f_sol[i])**2
    x1 = -(x1/L_x-(A-C)/2*(time[i])**2)**2/(g_sol[i])**2
    z = H*(1+x1+y1)/(f_sol[i]*g_sol[i])

    return z


def u(solSimi, x1, y1, i, L_x, U):
    """ get flow velocity in x direction from f and g solutions
    """
    time = solSimi['time']
    g_sol = solSimi['g_sol']
    g_p_sol = solSimi['g_p_sol']
    z = U*((A-C)*time[i]+(x1/L_x-(A-C)/2*(time[i])**2)*g_p_sol[i]/g_sol[i])

    return z


def v(solSimi, x1, y1, i, L_y, V):
    """ get flow velocity in y direction from f and g solutions
    """
    f_sol = solSimi['f_sol']
    f_p_sol = solSimi['f_p_sol']
    z = V*y1/L_y*f_p_sol[i]/f_sol[i]

    return z


def xc(solSimi, x1, y1, i, L_x):
    """ get center of mass location
    """
    time = solSimi['time']
    z = L_x*(A-C)/2*(time[i])**2

    return z


def getSimiSolParameters(solSimi, relDict, ind_time, cfg):
    """ get flow depth, flow velocity and center location of flow mass of similarity solution
        for required time step

        Parameters
        -----------
        solSimi: dict
            similarity solution
        relDict: dict
            dictionary with info of release
        ind_time: int
            index for required time step in similarity solution
        cfg: dict
            configuration

        Returns
        --------
        simiDict: dict
            dictionary of similiarty solution with flow depth, flow velocity,
            and center location in x for required time step
        """

    cfgSimi = cfg['SIMISOL']
    L_x = cfgSimi.getfloat('L_x')
    L_y = cfgSimi.getfloat('L_y')
    Hini = cfg['GENERAL'].getfloat('relTh')
    gravAcc = cfg['GENERAL'].getfloat('gravAcc')
    cos = relDict['cos']
    sin = relDict['sin']
    X1 = relDict['X1']
    Y1 = relDict['Y1']

    # Dimensioning parameters
    U = np.sqrt(gravAcc*L_x)
    V = np.sqrt(gravAcc*L_y)

    # get simi sol
    hSimi = h(solSimi, X1, Y1, ind_time, L_y, L_x, Hini)
    hSimi = np.where(hSimi <= 0, 0, hSimi)
    uxSimi = u(solSimi, X1, Y1, ind_time, L_x, U)
    uxSimi = np.where(hSimi <= 0, 0, uxSimi)
    uySimi = v(solSimi, X1, Y1, ind_time, L_y, V)
    uySimi = np.where(hSimi <= 0, 0, uySimi)
    vSimi = DFAtls.norm(uxSimi, uySimi, 0*uySimi)
    xCenter = xc(solSimi, X1, Y1, ind_time, L_x)*cos

    simiDict = {'hSimi': hSimi, 'vSimi': vSimi, 'vxSimi': uxSimi*cos, 'vySimi': uySimi, 'vzSimi': -uxSimi*sin,
                'xCenter': xCenter}

    return simiDict

##########################
# Analyze and compare analytic to numerical solution
#########################


def analyzeResults(particlesList, fieldsList, solSimi, relDict, cfg, outDirTest, index):
    """Compare analytical and com1DFA results"""

    hErrorL2Array = np.zeros((len(particlesList)))
    vErrorL2Array = np.zeros((len(particlesList)))
    hErrorLMaxArray = np.zeros((len(particlesList)))
    vErrorLMaxArray = np.zeros((len(particlesList)))
    count = 0
    time = []
    # run the comparison routine for each saved time step
    for particles, field in zip(particlesList, fieldsList):
        t = particles['t']
        time.append(t)
        # get similartiy solution h, u at required time step
        ind_time = np.searchsorted(solSimi['Time'], t)
        simiDict = getSimiSolParameters(solSimi, relDict, ind_time, cfg)
        cellSize = relDict['dem']['header']['cellsize']
        cosAngle = relDict['cos']
        hSimi = simiDict['hSimi']
        hNumerical = field['FD']
        vSimi = {'fx': simiDict['vxSimi'], 'fy': simiDict['vySimi'], 'fz': simiDict['vzSimi']}
        vNumerical = {'fx': field['Vx'], 'fy': field['Vy'], 'fz': field['Vz']}
        hErrorL2, hErrorLmax = anaTools.normL2Scal(hSimi, hNumerical, cellSize, cosAngle)
        hErrorL2Array[count] = hErrorL2
        hErrorLMaxArray[count] = hErrorLmax
        log.debug("L2 error on the Flow Depth at t=%.2f s is : %.4f" % (t, hErrorL2))
        vErrorL2, vErrorLmax = anaTools.normL2Vect(vSimi, vNumerical, cellSize, cosAngle)
        vErrorL2Array[count] = vErrorL2
        vErrorLMaxArray[count] = vErrorLmax
        log.debug("L2 error on the Flow velocity at t=%.2f s is : %.4f" % (t, vErrorL2))
        count = count + 1
    outAna1Plots.plotErrorTime(time, hErrorL2Array, hErrorLMaxArray, vErrorL2Array, vErrorLMaxArray, outDirTest, index)

    return hErrorL2Array, hErrorLMaxArray, vErrorL2Array, vErrorLMaxArray


def getReleaseThickness(avaDir, cfg, demFile):
    """ Define release thickness for the similarity solution test

    Release area is defined as an elipse or main radius Lx and Ly.
    Release thickness has a parabolic shape from relTh in the
    center to 0 on the edges

    Parameters
    -----------
    avaDir: str
        path to avalanche directory
    cfg: dict
        confguration settings
    demFile: str
        path to DEM file

    Returns
    --------
    relDict: dict
        dictionary with info on release thickness distribution

    """

    # Read dem
    demOri = IOf.readRaster(demFile)
    cszNew = cfg.getfloat('GENERAL', 'meshCellSize')
    x, y, xNew, yNew, diffExtentX, diffExtentY = geoTrans.getMeshXY(demOri, cellSizeNew=cszNew)

    nrows = len(yNew)
    ncols = len(xNew)
    xllc = demOri['header']['xllcenter']
    yllc = demOri['header']['yllcenter']
    csz = cszNew

    # define release thickness distribution
    cfgSimi = cfg['SIMISOL']
    L_x = cfgSimi.getfloat('L_x')
    L_y = cfgSimi.getfloat('L_y')
    Hini = cfg['GENERAL'].getfloat('relTh')
    planeinclinationAngleDeg = cfgSimi.getfloat('planeinclinationAngle')
    x = np.linspace(0, ncols-1, ncols)*csz+xllc
    y = np.linspace(0, nrows-1, nrows)*csz+yllc
    X, Y = np.meshgrid(x, y)
    cos = math.cos(math.pi*planeinclinationAngleDeg/180)
    sin = math.sin(math.pi*planeinclinationAngleDeg/180)
    X1 = X/cos
    Y1 = Y
    relTh = Hini * (1 - X1*X1/(L_x*L_x) - Y*Y/(L_y*L_y))
    relTh = np.where(relTh < 0, 0, relTh)

    relDict = {'relTh': relTh, 'X1': X1, 'Y1': Y1, 'demOri': demOri, 'X': X, 'Y': Y,
               'cos': cos, 'sin': sin}

    return relDict


def prepareParticlesFieldscom1DFA(fields, particles, relDict, simiDict, axis):
    """ get fields and particles dictionaries for given time step, for com1DFA domain origin is set to 0,0
        for particles - so info on domain is required

        Parameters
        -----------

        Fields: dictionary
            fields dictionary
        Particles: dictionary
            particles dictionary
        relDict: dict
            dictionary with info on release area
        simiDict: dict
            dictionary with center location in x for similarity solution
        axis: str
            axis (x or y) for profile

        Returns
        --------
        com1DFASol: dict
            dictionary with location of particles, flow depth, flow velocity,
            fields, and index for x or y cut of domain at the required time step

    """

    dem = relDict['dem']
    demOri = relDict['demOri']
    cos = relDict['cos']
    sin = relDict['sin']
    xCenter = simiDict['xCenter']

    # get info on DEM extent
    nrows = dem['header']['nrows']
    xllc = demOri['header']['xllcenter']
    yllc = demOri['header']['yllcenter']
    csz = dem['header']['cellsize']

    if axis == 'xaxis':
        ind = np.where(((particles['y']+yllc > -csz/2) & (particles['y']+yllc < csz/2)))
        indFinal = int(nrows * 0.5) -1
    elif axis == 'yaxis':
        ind = np.where(((particles['x']+xllc > xCenter-csz/2) & (particles['x']+xllc < xCenter+csz/2)))
        indFinal = int(np.round((xCenter - xllc)/csz) + 1)

    x = particles['x'][ind]+xllc
    y = particles['y'][ind]+yllc
    h = particles['h'][ind]
    ux = particles['ux'][ind]
    uy = particles['uy'][ind]
    uz = particles['uz'][ind]
    v = DFAtls.norm(ux, uy, uz)

    com1DFASol = {'x': x, 'y': y, 'h': h, 'v': v, 'vx': ux, 'vy': uy, 'vz': uz, 'indFinal': indFinal, 'fields': fields}

    return com1DFASol
