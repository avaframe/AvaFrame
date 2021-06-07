""" This script calculates the similarity solution for a gliding avalanche on
a inclined plane according to similarity solution from :
Hutter, K., Siegel, M., Savage, S.B. et al.
Two-dimensional spreading of a granular avalanche down an inclined plane
Part I. theory. Acta Mechanica 100, 37–68 (1993).
https://doi.org/10.1007/BF01176861
"""
import numpy as np
from scipy.integrate import ode
import math


Pi = math.pi
gravAcc = 9.81
# Set Parameters

planeinclinationAngle = 35.0 * Pi /180  # plane inclination
bedFrictionAngle = 30 * Pi /180  # basal angle of friction. 16deg in samosAT
internalFrictionAngle = 30 * Pi /180  # internal angle of friction phi>delta. 35deg in samosAT
# Set Parameters
zeta = planeinclinationAngle  # plane inclination
phi = internalFrictionAngle  # basal angle of friction. 16deg in samosAT
delta = bedFrictionAngle  # internal angle of friction phi>delta. 35deg in samosAT
# Dimentioning parameters
H = 8
L_x = 80.0
L_y = 80.0
U = np.sqrt(gravAcc*L_x)
V = np.sqrt(gravAcc*L_y)
T = np.sqrt(L_x/gravAcc)
# calculate aspect ratios
eps_x = H/L_x
eps_y = H/L_y
eps_xy = L_y/L_x
# Full scale end time
T_end = 55
# Non dimentional time for similarity sim calculation
t_1 = 0.1  # start time for ode solvers
t_end = T_end/T  # end time
dt_early = 0.01  # time step for early sol
dt = 0.01  # time step for early sol

# Initial conditions [g0 g_p0 f0 f_p0]
x_0 = [1.0, 0.0, 1.0, 0.0]  # here a circle as start point

##############################################################################


def define_earth_press_coeff():
    """Define earth pressure coefficients
    returns a vector of size 6 containing:
    [Kxact Kxpass Kyact(Kxact) Kypass(Kxact) Kyact(Kxpass) Kypass(Kxpass)]
    """

    cos2phi = np.cos(phi)**2
    cos2delta = np.cos(delta)**2
    tan2delta = np.tan(delta)**2
    root1 = np.sqrt(1.0 - cos2phi / cos2delta)
    earthPressureCoefficients = np.empty((6, 1))
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


def compute_earth_press_coeff(x, earthPressureCoefficients):
    """Compute earth pressure coefficients function of sng of f and g
    i.e depending on if we are in the active or passive case
    """
    g_p = x[1]
    f_p = x[3]
    if g_p >= 0:
        K_x = earthPressureCoefficients[0]
        if f_p >= 0:
            K_y = earthPressureCoefficients[2]
        else:
            print('ky passive')
            K_y = earthPressureCoefficients[3]
    else:
        print('kx passive')
        K_x = earthPressureCoefficients[1]
        if f_p >= 0:
            K_y = earthPressureCoefficients[4]
        else:
            print('ky passive')
            K_y = earthPressureCoefficients[5]
    return [K_x, K_y]


def compute_F_coeff(x, K_x, K_y):
    """Compute coefficients for function F
    coefficients for the simplified mode, eq 3.1
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
    return [A, B, C, D, E]


def calc_early_sol(t, earthPressureCoefficients):
    """ Compute the early solution for 0<t<t_1 to avoid singularity in the
    Runje-Kutta integration process
    """
    # early solution exists only if first derivative of f at t=0 is zero
    assert x_0[3] == 0, "f'(t=0)=f_p0 must be equal to 0"
    [K_x, K_y] = compute_earth_press_coeff(x_0, earthPressureCoefficients)
    [A, B, C, D, E] = compute_F_coeff(x_0, K_x, K_y)
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


def F_function(t, x, earthPressureCoefficients):
    """Calculate right hand side of the differential equation :
    dx/dt = F(x,t) F is discribed in Hutter 1993.
    inputs:
    t -- curent time
    x -- column vector of size 4
    output:
    F -- column vector of size 4

    """
    global A, C
    [K_x, K_y] = compute_earth_press_coeff(x, earthPressureCoefficients)
    [A, B, C, D, E] = compute_F_coeff(x, K_x, K_y)
    # print(A,B,C,D)
    u_c = (A - C)*t
    # print(u_c)
    g = x[0]
    g_p = x[1]
    f = x[2]
    f_p = x[3]

    dx0 = g_p
    dx1 = 2*B/(g**2*f)
    dx2 = f_p
    if C == 0:
        dx3 = 2*D/(2*f**2)
    else:
        dx3 = 2*D/(2*f**2)-C*f_p/u_c
    F = [dx0, dx1, dx2, dx3]
    return F


def ode_solver(solver, t_end, solSimi):
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


def h(solSimi, x1, y1, i):
    time = solSimi['time']
    g_sol = solSimi['g_sol']
    f_sol = solSimi['f_sol']
    y1 = -(y1/L_y)**2/(f_sol[i])**2
    x1 = -(x1/L_x-(A-C)/2*(time[i])**2)**2/(g_sol[i])**2
    z = H*(1+x1+y1)/(f_sol[i]*g_sol[i])
    return z


def u(solSimi, x1, y1, i):
    time = solSimi['time']
    g_sol = solSimi['g_sol']
    g_p_sol = solSimi['g_p_sol']
    z = U*((A-C)*time[i]+(x1/L_x-(A-C)/2*(time[i])**2)*g_p_sol[i]/g_sol[i])
    return z


def v(solSimi, x1, y1, i):
    f_sol = solSimi['f_sol']
    f_p_sol = solSimi['f_p_sol']
    z = V*y1/L_y*f_p_sol[i]/f_sol[i]
    return z


def xc(solSimi, x1, y1, i):
    time = solSimi['time']
    z = L_x*(A-C)/2*(time[i])**2
    return z


def write_raster_file(z_mat, x1, y1,FileName_ext):
    # Define domain characteristics
    ncols = len(x1)
    nrows = len(y1)
    xllcorner = 0.0
    yllcorner = 0.0
    cellsize = 5
    noDATA = -9999

    with open(FileName_ext,'w') as f:
        f.write('ncols  %d\n' % (ncols))
        f.write('nrows  %d\n' % (nrows))
        f.write('xllcorner  %.02f\n' % (xllcorner))
        f.write('yllcorner %.02f\n' % (yllcorner))
        f.write('cellsize  %d\n' % (cellsize))
        f.write('NODATA_value %d\n' % (noDATA))
        for line in z_mat:
            np.savetxt(f, line, fmt='%.2f')


def create_raster_file(solSimi, x1, y1, i, FileName):
    # Define domain characteristics
    ncols = np.shape(x1)[0]
    nrows = np.shape(y1)[0]
    xllcorner = 0.0
    yllcorner = 0.0
    cellsize = 5
    noDATA = -9999.00

    zh = h(solSimi, x1, y1, i)
    zu = u(solSimi, x1, y1, i)
    zu = np.where(zh > 0.0, zu, 0)
    zv = v(solSimi, x1, y1, i)
    zv = np.where(zh > 0.0, zv, 0)
    zh = np.where(zh > 0.0, zh, 0)
    # Save elevation data to .asc file and add header lines

    FileName_ext = FileName + '_fd.asc'
    z_mat = np.matrix(zh)
    write_raster_file(z_mat, x1, y1, FileName_ext)

    FileName_ext = FileName + '_vx.asc'
    z_mat = np.matrix(zu)
    write_raster_file(z_mat, x1, y1, FileName_ext)

    FileName_ext = FileName + '_vy.asc'
    z_mat = np.matrix(zv)
    write_raster_file(z_mat, x1, y1, FileName_ext)


def runSimilarity():
    """ Run main model"""
    earthPressureCoefficients = define_earth_press_coeff()
    # Early time solution
    t_early = np.arange(0, t_1, dt_early)
    solSimi = calc_early_sol(t_early, earthPressureCoefficients)
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
    solver = ode(F_function)
    solver.set_integrator('dopri5')
    solver.set_f_params(earthPressureCoefficients)
    solver.set_initial_value(x_1, t_start)
    solSimi = ode_solver(solver, t_end, solSimi)

    # x1 = x0/np.cos(zeta)
    # y1 = y0

    # plt.grid(gridflag); plt.xlabel(xlabeltext); plt.ylabel(ylabeltext); plt.title(titletext)
    # ax.set_aspect('equal', adjustable='datalim')

    Time = solSimi['time']*T
    solSimi['Time'] = Time
    # t = 10
    #
    # ind_time = np.searchsorted(Time, t)
    #
    # FileName = './OutputAna_release01_t%012.4fs' %t
    # # FileName = 'SimilaritySolAnalytic/OutputAna_release01_t%012.4fs' %t
    # create_raster_file(solSimi, x1, y1, ind_time, FileName)

    return solSimi
