import numpy as np
import math
from scipy.optimize import minimize
import matplotlib.pyplot as plt


# Local imports
import avaframe.com1DFA.DFAtools as DFAtls
import avaframe.out3Plot.plotUtils as pU


def runSplitPartTest():
    # splitting pattern
    # pattern 1: particle in the center and others arrond at distance alpha * R0
    # pattern 2: same as pattern 1 without the particle in the center
    splittingPattern = 2
    # numper of particles after splitting
    nPartSplit = 7
    # initialize C, b and Q
    Ch = 0
    bh = np.zeros((nPartSplit, 1))
    Qh = np.zeros((nPartSplit, nPartSplit))
    CgradH = 0
    bgradH = np.zeros((nPartSplit, 1))
    QgradH = np.zeros((nPartSplit, nPartSplit))

    # kernel radius
    R0 = 1
    minRKern = 0.01
    # number of steps f the integral approximation
    nSteps = 101
    dx = R0*(nSteps-1)/2
    area = dx*dx
    # create grid for integral approximation
    x = np.linspace(-2*R0, 2*R0, num=nSteps)
    xGrid, yGrid = np.meshgrid(x, x)
    zGrid = np.zeros(np.shape(xGrid))

    N = 50
    M = 50
    EPSILON = np.linspace(0.1, 0.8, M)
    R0NewList = np.linspace(0.3, 1, N)
    errorHEps = np.zeros((N, M))
    errorgradHEps = np.zeros((N, M))
    M0 = np.zeros((N, M))
    M1 = np.zeros((N, M))
    for R0New, indR0 in zip(R0NewList, range(N)):
        for eps, indEps in zip(EPSILON, range(M)):
            print(eps, R0New)
            hTot = np.zeros(np.shape(xGrid))
            if splittingPattern == 1:
                xPart = np.zeros(nPartSplit)
                yPart = np.zeros(nPartSplit)
                zPart = np.zeros(nPartSplit)
                alpha = 2*math.pi/(nPartSplit-1)*np.arange(nPartSplit-1)
                xEps = eps*np.cos(alpha)
                yEps = eps*np.sin(alpha)
                xPart[1:] = xEps
                yPart[1:] = yEps
            elif splittingPattern == 2:
                alpha = 2*math.pi/nPartSplit*np.arange(nPartSplit)
                xPart = eps*np.cos(alpha)
                yPart = eps*np.sin(alpha)
                zPart = np.zeros(nPartSplit)

            # first compute the grad(W_i) in x and y dir
            rx = xGrid
            ry = yGrid
            rz = zGrid
            hi, gradiX, gradiY, gradiZ = getGradW(R0, minRKern, rx, ry, rz)
            h2i = hi*hi
            normGradi = DFAtls.norm2(gradiX, gradiY, gradiZ)
            Ch = np.sum(h2i)
            CgradH = np.sum(normGradi)

            for k in range(nPartSplit):
                rkx = xGrid - xPart[k]
                rky = yGrid - yPart[k]
                rkz = zGrid - zPart[k]
                hk, gradkX, gradkY, gradkZ = getGradW(R0New, minRKern, rkx, rky, rkz)
                gradiGradk = DFAtls.scalProd(gradiX, gradiY, gradiZ, gradkX, gradkY, gradkZ)
                bh[k] = np.sum(hk*hi)
                bgradH[k] = np.sum(gradiGradk)
                for l in range(nPartSplit):
                    rlx = xGrid - xPart[l]
                    rly = yGrid - yPart[l]
                    rlz = zGrid - zPart[l]
                    hl, gradlX, gradlY, gradlZ = getGradW(R0New, minRKern, rlx, rly, rlz)
                    gradkGradl = DFAtls.scalProd(gradkX, gradkY, gradkZ, gradlX, gradlY, gradlZ)
                    sum = np.sum(gradkGradl)
                    Qh[k, l] = np.sum(hl*hk)
                    QgradH[k, l] = sum
                    # if l < k:
                    #     Q[l, k] = -sum
            # print(C, b, Q)
            if splittingPattern == 1:
                # solve the min problem to get the lambdas
                # constraint: sum of lambda should be 1
                cons = ({'type': 'eq', 'fun': lambda x:  x[0] + (nPartSplit-1)*x[1] - 1})
                x0 = np.ones((2, 1)) / nPartSplit
                bnds = ((0.001, None), (0.1369, None))
                res = minimize(error2, x0, args=(Ch, bh, Qh, nPartSplit), method='SLSQP', bounds=bnds, constraints=cons, tol=None)
                errH = error2(res.x, Ch, bh, Qh, nPartSplit)
                errgradH = error2(res.x, CgradH, bgradH, QgradH, nPartSplit)
                lam = np.ones(nPartSplit) * res.x[1]
                lam[0] = res.x[0]
                M0[indR0, indEps] = res.x[0]
                M1[indR0, indEps] = res.x[1]
            elif splittingPattern == 2:
                lam = np.ones((nPartSplit, 1)) / nPartSplit
                errH = error(lam, Ch, bh, Qh)
                errgradH = error(lam, CgradH, bgradH, QgradH)
            errorHEps[indR0, indEps] = errH/nSteps/nSteps
            errorgradHEps[indR0, indEps] = errgradH/nSteps/nSteps

            # print(lam)
            # for k in range(nPartSplit):
            #     rkx = xGrid - xPart[k]
            #     rky = yGrid - yPart[k]
            #     rkz = zGrid - zPart[k]
            #     hk, gradkX, gradkY, gradkZ = getGradW(R0New, minRKern, rkx, rky, rkz)
            #     hTot = hTot + lam[k]*hk
            #
            # rx = xGrid
            # ry = yGrid
            # rz = zGrid
            # hi, gradiX, gradiY, gradiZ = getGradW(R0, minRKern, rx, ry, rz)
            # print(np.nanmax(hTot), np.nanmax(hi))
            # cmap, _, ticks, norm = pU.makeColorMap(pU.cmapDepth, np.nanmin(hi), np.nanmax(hi), continuous=True)
            # fig, ax = plt.subplots(figsize=(3*pU.figW, 3*pU.figH))
            # Cp1 = ax.contour(xGrid, yGrid, hTot, levels=20,  cmap=cmap)
            # Cp2 = ax.contour(xGrid, yGrid, hi, levels=20,  cmap=cmap, linestyles='dashed')
            # ax.clabel(Cp1, inline=1, fontsize=8)
            # ax.clabel(Cp2, inline=1, fontsize=8)
            # plt.show()
            # plt.close(fig)

    cmap, _, ticks, norm = pU.makeColorMap(pU.cmapDepth, np.nanmin(M0), np.nanmax(M0), continuous=True)

    fig1, ax1 = plt.subplots(figsize=(pU.figW, pU.figH))
    ref0, im = pU.NonUnifIm(ax1, EPSILON, R0NewList , M0, 'epsilon', 'R0',
                            extent=[EPSILON.min(), EPSILON.max(), R0NewList.min(), R0NewList.max()],
                            cmap=cmap, norm=None)
    Cp1 = ax1.contour(EPSILON, R0NewList, M0, levels=20, colors='k')
    pU.addColorBar(im, ax1, ticks, '')
    cmap, _, ticks, norm = pU.makeColorMap(pU.cmapDepth, np.nanmin(M1), np.nanmax(M1), continuous=True)
    fig2, ax2 = plt.subplots(figsize=(pU.figW, pU.figH))
    ref0, im = pU.NonUnifIm(ax2, EPSILON, R0NewList, M1, 'epsilon', 'R0',
                            extent=[EPSILON.min(), EPSILON.max(), R0NewList.min(), R0NewList.max()],
                            cmap=cmap, norm=None)
    Cp1 = ax2.contour(EPSILON, R0NewList, M1, levels=20, colors='k')
    pU.addColorBar(im, ax2, ticks, '')

    cmap, _, ticks, norm = pU.makeColorMap(pU.cmapDepth, np.nanmin(np.log10(errorHEps)), np.nanmax(np.log10(errorHEps)), continuous=True)
    fig3, ax3 = plt.subplots(figsize=(pU.figW, pU.figH))
    ref0, im = pU.NonUnifIm(ax3, EPSILON, R0NewList, np.log10(errorHEps), 'epsilon', 'R0',
                            extent=[EPSILON.min(), EPSILON.max(), R0NewList.min(), R0NewList.max()],
                            cmap=cmap, norm=None)
    Cp1 = ax3.contour(EPSILON, R0NewList, np.log10(errorHEps), levels=20, colors='k')
    pU.addColorBar(im, ax3, ticks, '')

    cmap, _, ticks, norm = pU.makeColorMap(pU.cmapDepth, np.nanmin(np.log10(errorgradHEps)), np.nanmax(np.log10(errorgradHEps)), continuous=True)
    fig4, ax4 = plt.subplots(figsize=(pU.figW, pU.figH))
    ref0, im = pU.NonUnifIm(ax4, EPSILON, R0NewList, np.log10(errorgradHEps), 'epsilon', 'R0',
                            extent=[EPSILON.min(), EPSILON.max(), R0NewList.min(), R0NewList.max()],
                            cmap=cmap, norm=None)
    Cp1 = ax4.contour(EPSILON, R0NewList, np.log10(errorgradHEps), levels=20, colors='k')
    pU.addColorBar(im, ax4, ticks, '')

    plt.show()


def getGradW(R0, minRKern, rx, ry, rz):
    r = DFAtls.norm(rx, ry, rz)
    facKernel = 10.0 / (math.pi * R0**5)
    dfacKernel = - 3.0 * facKernel
    indSmall = np.where(r < minRKern * R0)
    # impose a minimum distance between particles
    rx[indSmall] = minRKern * R0 * rx[indSmall]
    ry[indSmall] = minRKern * R0 * ry[indSmall]
    rz[indSmall] = minRKern * R0 * rz[indSmall]
    r[indSmall] = minRKern * R0
    indBig = np.where(r >= R0)
    hr = R0 - r
    h = facKernel * hr * hr * hr
    dwdr = dfacKernel * hr * hr
    mdwdrr = dwdr / r
    mdwdrr[indBig] = 0
    h[indBig] = 0

    gradX = mdwdrr*rx
    gradY = mdwdrr*ry
    gradZ = mdwdrr*rz
    return h, gradX, gradY, gradZ


def error(x, C, b, Q):
    return C - 2*np.matmul(np.transpose(x), b) + np.matmul(np.matmul(np.transpose(x), Q), x)


def error2(x, C, b, Q, nPartSplit):
    xTemp = np.ones((nPartSplit, 1))*x[1]
    xTemp[0] = x[0]
    err = (C - 2*np.matmul(np.transpose(xTemp), b) + np.matmul(np.matmul(np.transpose(xTemp), Q), xTemp))[0, 0]
    return err


if __name__ == '__main__':
    runSplitPartTest()
