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
    _, _, _, _, _, _, simDF = com1DFA.com1DFAMain(avalancheDir, cfgMain, cfgFile=simiSolCfg, relThField=relTh)

    # simDF = cfgUtils.createConfigurationInfo(avalancheDir, standardCfg='', writeCSV=False)

    # compute the similartiy solution (this corresponds to our reference)
    log.info('Computing similarity solution')
    solSimi = mainSimilaritySol(simiSolCfg)

    # now compare the simulations to the reference
    # first fetch info about all the simulations performed (and maybe order them)
    varParList = cfg['ANALYSIS']['varParList'].split('|')
    ascendingOrder = cfg['ANALYSIS']['ascendingOrder']
    # load info for all configurations and order them
    simDF = simDF.sort_values(by=varParList, ascending=ascendingOrder)
    simDF = postProcessSimiSol(avalancheDir, cfgMain, cfg['SIMISOL'], simDF, solSimi, outDirTest)
    outAna1Plots.plotError(simDF, outDirTest)
    outAna1Plots.plotErrorLog(simDF, outDirTest, cfg['SIMISOL'])


def mainSimilaritySol(simiSolCfg):
    """ Compute similarity solution
    Parameters
    -----------
    simiSolCfg: pathlib path
        path to simiSol configuration file
    Returns
    ---------
    solSimi: dictionary
        similarity solution:
            time: time array (without dimention)
            Time: time array (with dimention)
            g_sol: g array
            g_p_sol: first derivativ of g array
            f_sol: f array
            f_p_sol: first derivativ of f array

    """

    # Load configuration
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
    t_1 = 0.1         # start time for ode solvers, end time for early time sol (we need t_1<<1)
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
    t_early = np.append(t_early, t_1)
    solSimi = calcEarlySol(t_early, earthPressureCoefficients, x_0, zeta, delta, eps_x, eps_xy, eps_y)

    # Runge-Kutta integration away from the singularity
    # initial conditions
    t_start = t_1
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


def computeFCoeff(K_x, K_y, zeta, delta, eps_x, eps_xy, eps_y):
    """ Compute coefficients eq 3.2 for the function F
        Parameters
        -----------
        K_x: float
            Kx earth pressure coef
        K_y: float
            Ky earth pressure coef
        zeta: float
            slope angle
        delta: float
            friction angle
        eps_x: float
            scale in x dir
        eps_xy: float
            scale in x/y dir
        eps_y: float
            scale in y dir
        Returns
        ---------
        A, B, C, D, E: floats
            coefficients of eq 3.2
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
    """ Compute the early solution for 0<t<t_1
        to avoid singularity in the Runge-Kutta integration process
        Parameters
        -----------
        t: numpy array
            time array
        earthPressureCoefficients: numpy array
            earth Pressure Coefficients
        x_0: numpy array
            initial condition
        zeta: float
            slope angle
        delta: float
            friction angle
        eps_x: float
            scale in x dir
        eps_xy: float
            scale in x/y dir
        eps_y: float
            scale in y dir
        Returns
        ---------
        solSimi: dictionary
            similarity solution (for early times):
                time: time array (without dimention)
                g_sol: g array
                g_p_sol: first derivativ of g array
                f_sol: f array
                f_p_sol: first derivativ of f array

    """

    # early solution exists only if first derivative of f at t=0 is zero
    assert x_0[3] == 0, "f'(t=0)=f_p0 must be equal to 0"
    K_x, K_y = computeEarthPressCoeff(x_0, earthPressureCoefficients)
    A, B, C, D, E = computeFCoeff(K_x, K_y, zeta, delta, eps_x, eps_xy, eps_y)
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
            initial condition, column vector of size 4
        earthPressureCoefficients: numpy array
            earth Pressure Coefficients
        zeta: float
            slope angle
        delta: float
            friction angle
        eps_x: float
            scale in x dir
        eps_xy: float
            scale in x/y dir
        eps_y: float
            scale in y dir

        Returns:
        F: numpy array
            column vector of size 4
    """

    K_x, K_y = computeEarthPressCoeff(x, earthPressureCoefficients)
    A, B, C, D, E = computeFCoeff(K_x, K_y, zeta, delta, eps_x, eps_xy, eps_y)
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
        Parameters
        -----------
        solver: `ode` instance
            `ode` instance corresponding to out ODE
        dt: float
            time step
        t_end: float
            end time
        solSimi: dictionary
            similarity solution (for early times):
                time: time array (without dimention)
                g_sol: g array
                g_p_sol: first derivativ of g array
                f_sol: f array
                f_p_sol: first derivativ of f array
        Returns
        ---------
        solSimi: dictionary
            similarity solution (copleted with all time steps):
                time: time array (without dimention)
                g_sol: g array
                g_p_sol: first derivativ of g array
                f_sol: f array
                f_p_sol: first derivativ of f array
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


def computeH(solSimi, x1, y1, i, L_x, L_y, H, AminusC):
    """ get flow depth from f and g solutions
        Parameters
        -----------
        solSimi: dictionary
            similarity solution
        x1: numpy array
            x coordinate location desiered for the solution
        y1: numpy array
            y coordinate location desiered for the solution
        i: int
            time simHash
        L_x: float
            scale in x dir
        L_y: float
            scale in y dir
        H: float
            scale in z dir
        AminusC:
            A-C coefficient

        Returns
        --------
        h: numpy array
            h similarity solution at (x1, y1)

    """
    time = solSimi['time']
    g_sol = solSimi['g_sol']
    f_sol = solSimi['f_sol']
    y1 = -(y1/L_y)**2/(f_sol[i])**2
    x1 = -(x1/L_x-AminusC/2*(time[i])**2)**2/(g_sol[i])**2
    h = H*(1+x1+y1)/(f_sol[i]*g_sol[i])

    return h


def computeU(solSimi, x1, y1, i, L_x, U, AminusC):
    """ get flow velocity in x direction from f and g solutions
        Parameters
        -----------
        solSimi: dictionary
            similarity solution
        x1: numpy array
            x coordinate location desiered for the solution
        y1: numpy array
            y coordinate location desiered for the solution
        i: int
            time simHash
        L_x: float
            scale in x dir
        U: float
            x velocity component scale
        AminusC:
            A-C coefficient

        Returns
        --------
        u: numpy array
            u similarity solution at (x1, y1)
    """
    time = solSimi['time']
    g_sol = solSimi['g_sol']
    g_p_sol = solSimi['g_p_sol']
    u = U*(AminusC*time[i]+(x1/L_x-AminusC/2*(time[i])**2)*g_p_sol[i]/g_sol[i])

    return u


def computeV(solSimi, x1, y1, i, L_y, V):
    """ get flow velocity in y direction from f and g solutions
        Parameters
        -----------
        solSimi: dictionary
            similarity solution
        x1: numpy array
            x coordinate location desiered for the solution
        y1: numpy array
            y coordinate location desiered for the solution
        i: int
            time simHash
        L_y: float
            scale in y dir
        V: float
            y velocity component scale

        Returns
        --------
        v: numpy array
            v similarity solution at (x1, y1)
    """
    f_sol = solSimi['f_sol']
    f_p_sol = solSimi['f_p_sol']
    v = V*y1/L_y*f_p_sol[i]/f_sol[i]

    return v


def computeXC(solSimi, x1, y1, i, L_x, AminusC):
    """ get center of mass location
        Parameters
        -----------
        solSimi: dictionary
            similarity solution
        x1: numpy array
            x coordinate location desiered for the solution
        y1: numpy array
            y coordinate location desiered for the solution
        i: int
            time simHash
        L_x: float
            scale in x dir
        AminusC:
            A-C coefficient

        Returns
        --------
        xc: numpy array
            x position of the center of the similarity solution pile
    """
    time = solSimi['time']
    xc = L_x*AminusC/2*(time[i])**2

    return xc


def getSimiSolParameters(solSimi, header, ind_time, cfgSimi, Hini, gravAcc):
    """ get flow depth, flow velocity and center location of flow mass of similarity solution
        for required time step

        Parameters
        -----------
        solSimi: dict
            similarity solution
        header: dict
            header dictionary with info about the extend and cell size
        ind_time: int
            simHash for required time step in similarity solution
        cfg: dict
            configuration
        Hini: float
            initial release depth
        gravAcc: float
            gravity acceleration

        Returns
        --------
        simiDict: dict
            dictionary of similiarty solution with flow depth, flow velocity,
            and center location in x for required time step
        """

    L_x = cfgSimi.getfloat('L_x')
    L_y = cfgSimi.getfloat('L_y')
    bedFrictionAngleDeg = cfgSimi.getfloat('bedFrictionAngle')
    planeinclinationAngleDeg = cfgSimi.getfloat('planeinclinationAngle')

    # Set parameters
    Pi = math.pi
    zeta = planeinclinationAngleDeg * Pi /180       # plane inclination
    delta = bedFrictionAngleDeg * Pi /180           # basal angle of friction

    cos = math.cos(zeta)
    sin = math.sin(zeta)

    # get info on DEM extent
    ncols = header['ncols']
    nrows = header['nrows']
    xllc = header['xllcenter']
    yllc = header['yllcenter']
    csz = header['cellsize']
    x = np.linspace(0, ncols-1, ncols)*csz+xllc
    y = np.linspace(0, nrows-1, nrows)*csz+yllc
    X, Y = np.meshgrid(x, y)
    X1 = X/cos
    Y1 = Y

    # Dimensioning parameters
    U = np.sqrt(gravAcc*L_x)
    V = np.sqrt(gravAcc*L_y)

    # A-C
    A = np.sin(zeta)
    C = np.cos(zeta) * np.tan(delta)
    AminusC = A - C

    # get simi sol
    hSimi = computeH(solSimi, X1, Y1, ind_time, L_x, L_y, Hini, AminusC)
    hSimi = np.where(hSimi <= 0, 0, hSimi)
    uxSimi = computeU(solSimi, X1, Y1, ind_time, L_x, U, AminusC)
    uxSimi = np.where(hSimi <= 0, 0, uxSimi)
    uySimi = computeV(solSimi, X1, Y1, ind_time, L_y, V)
    uySimi = np.where(hSimi <= 0, 0, uySimi)
    vSimi = DFAtls.norm(uxSimi, uySimi, 0*uySimi)
    xCenter = computeXC(solSimi, X1, Y1, ind_time, L_x, AminusC)*cos

    simiDict = {'hSimi': hSimi, 'vSimi': vSimi, 'vxSimi': uxSimi*cos, 'vySimi': uySimi, 'vzSimi': -uxSimi*sin,
                'xCenter': xCenter, 'cos': cos, 'sin': sin}

    return simiDict

##########################
# Analyze and compare analytic to numerical solution
#########################


def postProcessSimiSol(avalancheDir, cfgMain, cfgSimi, simDF, solSimi, outDirTest):
    # loop on all the simulations and make the comparison to reference
    for simHash, simDFrow in simDF.iterrows():
        simName = simDFrow['simName']
        # fetch the simulation results
        particlesList, Tsave = com1DFA.readPartFromPickle(avalancheDir, simName=simName, flagAvaDir=True, comModule='com1DFA')
        fieldsList, fieldHeader = com1DFA.readFields(avalancheDir, ['FD', 'FV', 'Vx', 'Vy', 'Vz'], simName=simName, flagAvaDir=True, comModule='com1DFA')
        simDF.loc[simHash, 'Npart'] = particlesList[0]['Npart']
        # analyze and compare results
        hEL2Array, hELMaxArray, vEL2Array, vELMaxArray = analyzeResults(particlesList, fieldsList, solSimi, fieldHeader,
                                                                        cfgSimi, outDirTest, simHash, simDFrow)
        # add result of error analysis
        # save results in the simDF
        tSave = cfgSimi.getfloat('tSave')
        ind_t = min(np.searchsorted(Tsave, tSave), min(len(Tsave)-1, len(fieldsList)-1))
        simDF.loc[simHash, 'hErrorL2'] = hEL2Array[ind_t]
        simDF.loc[simHash, 'vErrorL2'] = vEL2Array[ind_t]
        simDF.loc[simHash, 'hErrorLMax'] = hELMaxArray[ind_t]
        simDF.loc[simHash, 'vErrorLMax'] = vELMaxArray[ind_t]
        # +++++++++POSTPROCESS++++++++++++++++++++++++
        # -------------------------------
        # if cfgMain['FLAGS'].getboolean('showPlot'):
        #     outAna1Plots.plotContoursSimiSol(particlesList, fieldsList, solSimi, relDict, cfgSimi, outDirTest)

        # outAna1Plots.showSaveTimeSteps(cfgMain, cfgSimi, particlesList, fieldsList, solSimi, Tsave, fieldHeader,
        #                                outDirTest, simHash, simDFrow)

    simDF.to_pickle(outDirTest / 'results.p')

    return simDF


def analyzeResults(particlesList, fieldsList, solSimi, fieldHeader, cfgSimi, outDirTest, simHash, simDFrow):
    """Compare analytical and com1DFA results
        Parameters
        -----------
        particlesList: list
            list of particles dictionaries
        fieldsList: list
            list of fields dictionaries
        solSimi: dictionary
            similarity solution:
                time: time array (without dimention)
                Time: time array (with dimention)
                g_sol: g array
                g_p_sol: first derivativ of g array
                f_sol: f array
                f_p_sol: first derivativ of f array
        fieldHeader: dict
            header dictionary with info about the extend and cell size
        cfgSimi: dict
            simisol confguration section
        outDirTest: pathlib path
            path output directory (where to save the figures)
        simHash: str
            com1DFA simulation id

        Returns
        --------
        hErrorL2Array: numpy array
            L2 error on flow depth for saved time steps
        hErrorLMaxArray: numpy array
            LMax error on flow depth for saved time steps
        vErrorL2Array: numpy array
            L2 error on flow velocity for saved time steps
        vErrorLMaxArray: numpy array
            LMax error on flow velocity for saved time steps

    """
    relTh = simDFrow['relTh']
    gravAcc = simDFrow['gravAcc']
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
        simiDict = getSimiSolParameters(solSimi, fieldHeader, ind_time, cfgSimi, relTh, gravAcc)
        cellSize = fieldHeader['cellsize']
        cosAngle = simiDict['cos']
        hSimi = simiDict['hSimi']
        hNumerical = field['FD']
        vSimi = {'fx': simiDict['vxSimi'], 'fy': simiDict['vySimi'], 'fz': simiDict['vzSimi']}
        vNumerical = {'fx': field['Vx'], 'fy': field['Vy'], 'fz': field['Vz']}
        hErrorL2, hErrorL2Rel, hErrorLmax, hErrorLmaxRel = anaTools.normL2Scal(hSimi, hNumerical, cellSize, cosAngle)
        hErrorL2Array[count] = hErrorL2Rel
        hErrorLMaxArray[count] = hErrorLmaxRel
        log.debug("L2 error on the Flow Depth at t=%.2f s is : %.4f" % (t, hErrorL2))
        vErrorL2, vErrorL2Rel, vErrorLmax, vErrorLmaxRel = anaTools.normL2Vect(vSimi, vNumerical, cellSize, cosAngle)
        vErrorL2Array[count] = vErrorL2Rel
        vErrorLMaxArray[count] = vErrorLmaxRel
        log.debug("L2 error on the Flow velocity at t=%.2f s is : %.4f" % (t, vErrorL2))
        count = count + 1
    # outAna1Plots.plotErrorTime(time, hErrorL2Array, hErrorLMaxArray, vErrorL2Array, vErrorLMaxArray, outDirTest,
    #                            simHash, simDFrow)

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
    csz = cfg.getfloat('GENERAL', 'meshCellSize')
    x, y, xNew, yNew, diffExtentX, diffExtentY = geoTrans.getMeshXY(demOri, cellSizeNew=csz)

    nrows = len(yNew)
    ncols = len(xNew)
    xllc = demOri['header']['xllcenter']
    yllc = demOri['header']['yllcenter']

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


def prepareParticlesFieldscom1DFA(fields, particles, header, simiDict, axis):
    """ get fields and particles dictionaries for given time step, for com1DFA domain origin is set to 0,0
        for particles - so info on domain is required

        Parameters
        -----------

        Fields: dictionary
            fields dictionary
        Particles: dictionary
            particles dictionary
        header: dict
            header dictionary with info about the extend and cell size
        simiDict: dict
            dictionary with center location in x for similarity solution
        axis: str
            axis (x or y) for profile

        Returns
        --------
        com1DFASol: dict
            dictionary with location of particles, flow depth, flow velocity,
            fields, and simHash for x or y cut of domain at the required time step

    """

    xCenter = simiDict['xCenter']
    # get info on DEM extent
    ncols = header['ncols']
    nrows = header['nrows']
    xllc = header['xllcenter']
    yllc = header['yllcenter']
    csz = header['cellsize']

    xArrayFields = np.linspace(xllc, xllc+(ncols-1)*csz, ncols)
    yArrayFields = np.linspace(yllc, yllc+(nrows-1)*csz, nrows)

    if axis == 'xaxis':
        ind = np.where(((particles['y']+yllc > -csz) & (particles['y']+yllc < csz)))
        indFinal = int(nrows * 0.5) -1
    elif axis == 'yaxis':
        ind = np.where(((particles['x']+xllc > xCenter-csz) & (particles['x']+xllc < xCenter+csz)))
        indFinal = int(np.round((xCenter - xllc)/csz) + 1)

    x = particles['x'][ind]+xllc
    y = particles['y'][ind]+yllc
    h = particles['h'][ind]
    ux = particles['ux'][ind]
    uy = particles['uy'][ind]
    uz = particles['uz'][ind]
    v = DFAtls.norm(ux, uy, uz)

    com1DFASol = {'x': x, 'y': y, 'h': h, 'v': v, 'vx': ux, 'vy': uy, 'vz': uz, 'indFinal': indFinal,
                  'xArrayFields': xArrayFields, 'yArrayFields': yArrayFields, 'fields': fields}

    return com1DFASol