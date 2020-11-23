import numpy as np


def SamosATfric(cfg, v, p, h):
    rho = cfg.getfloat('rho')
    Rs0 = cfg.getfloat('Rs0')
    mu = cfg.getfloat('mu')
    kappa = cfg.getfloat('kappa')
    B = cfg.getfloat('B')
    R = cfg.getfloat('R')
    Rs = rho * v * v / (p + 0.001)
    div = h / R
    div = np.where(div < 1.0, 1, div)
    # if(div < 1.0):
    #     div = 1.0
    div = np.log(div) / kappa + B
    tau = p * mu * (1.0 + Rs0 / (Rs0 + Rs)) + rho * v * v / (div * div)
    return tau
