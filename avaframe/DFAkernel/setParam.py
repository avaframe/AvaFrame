import math
# acceleration of gravity [m/s²]
gravAcc = 9.81
# density of snow [kg/m³]
rho = 200
rhoEnt = 200
# mass per particle [kg]
massPerPart = 5000

velMagMin = 1.0e-6
depMin = 1.0e-6

# time step [s]
dt = 0.1
# save every [s]
dtSave = 1
# End time
Tend = 30

# SamosAt friction model
mu = 0.155
Rs0 = 0.222
kappa = 0.43
R = 0.05
B = 4.13

# Resistance force parameters
hRes = 1
# entrainment Erosion Energy
entEroEnergy = 0
entShearResistance = 0
entDefResistance = 0

# SPH kernel
# use "spiky" kernel: w = (h - r)**3 * 10/(pi*h**5)
rKernel = 5
facKernel = 10.0 / (math.pi * pow(rKernel, 5.0))
dfacKernel = -3.0 * facKernel
d2facKernel = -2.0 * dfacKernel
