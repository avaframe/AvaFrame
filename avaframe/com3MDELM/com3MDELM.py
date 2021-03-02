"""
    Main logic for MDELM computational module

    This file is part of Avaframe.
"""

import logging
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.image import NonUniformImage
import math

# Local imports
import avaframe.in3Utils.geoTrans as geoTrans
import avaframe.in2Trans.ascUtils as IOf
from avaframe.out3Plot.plotUtils import *

# create local logger
log = logging.getLogger(__name__)
debugPlot = False


def prepareDEM(demSource, resampleResolution):
    # Read input data for MDELM
    dem = IOf.readRaster(demSource)

    # make grid and change resolution
    header = dem['header']
    ncols = header.ncols
    nrows = header.nrows
    cellsize = header.cellsize
    n = int(np.floor((ncols-1)/resampleResolution))
    m = int(np.floor((nrows-1)/resampleResolution))
    print(ncols, nrows, n, m)
    csz = cellsize*resampleResolution
    log.info("Using new celle size: %s m", csz)
    xgrid = np.linspace(header.xllcorner, header.xllcorner+(n-1)*csz, n)
    ygrid = np.linspace(header.yllcorner, header.yllcorner+(m-1)*csz, m)
    PointsX, PointsY = np.meshgrid(xgrid, ygrid)
    Points = {}
    Points['x'] = PointsX
    Points['y'] = PointsY

    Points, _ = geoTrans.projectOnRaster(dem, Points,
                                                      interp='nearest')
    fig = plt.figure(figsize=(figW, figH))
    ax1 = plt.subplot(111)
    cmap = cmapDEM
    im0 = NonUniformImage(ax1, extent=[xgrid.min(), xgrid.max(),
                                       ygrid.min(), ygrid.max()], cmap=cmap)
    # im.set_interpolation('bilinear')
    im0.set_data(xgrid, ygrid, Points['z'])
    ref1 = ax1.images.append(im0)
    cbar = ax1.figure.colorbar(im0, ax=ax1, use_gridspec=True)
    cbar.ax.set_ylabel('Altitude [m]')
    ax1.title.set_text('DEM')
    ax1.set_xlim([xgrid.min(), xgrid.max()])
    ax1.set_ylim([ygrid.min(), ygrid.max()])
    ax1.set_xlabel(r'$x\;[m]$')
    ax1.set_ylabel(r'$y\;[m]$')
    # plt.show()
    plt.close(fig)

    return Points, csz


def calcV2(x, y, z, v2d, d, g, mu):
    xd = x[1, 1]
    yd = y[1, 1]
    zd = z[1, 1]
    # estimate local distance
    DeltaLLoc = np.sqrt(np.square(x-xd)+np.square(y-yd)+np.square(z-zd))
    DeltaSLoc = np.sqrt(np.square(x-xd)+np.square(y-yd))

    # compute slope angle
    bx = (z[1, 2]-zd)/(x[1, 2]-xd)
    bx2 = bx*bx
    by = (z[2, 1]-zd)/(y[2, 1]-yd)
    by2 = by*by
    c = np.sqrt(1+bx2+by2)

    # JT's formulation
    # v2rest = v2d - 2*g*((z-zd) + mu*((sn-sd)))
    # Matthias formulation
    # calculate V2 in the reciver cells
    v2rest = v2d - 2*g*((z-zd) + mu*DeltaLLoc/c)
    # calculate GradV2 in the reciver cells
    DeltaSLoc[1, 1] = 1
    gradV2 = -2*g*((z-zd) + mu*DeltaLLoc/c)/DeltaSLoc
    # update GradV2 with forbiden cells
    # if cell is donor
    gradV2 = np.where(np.isnan(d), np.nan, gradV2)
    # if cell has nan V2
    gradV2 = np.where(np.isnan(v2rest), np.nan, gradV2)
    # if cell has negative velocity2
    gradV2 = np.where(v2rest < 0, np.nan, gradV2)
    # update V2 with forbiden cells
    # if cell is donor
    v2rest = np.where(np.isnan(d), 0.0, v2rest)
    # if cell has nan V2
    v2rest = np.where(np.isnan(v2rest), 0.0, v2rest)
    # if cell has negative velocity2
    v2rest = np.where(v2rest < 0, 0.0, v2rest)

    return DeltaLLoc, DeltaSLoc, gradV2, v2rest


def wFunk(x):
    fx = np.arctan(x) + math.pi/2
    fx = fx/(math.pi/2)
    fx = np.power(fx, 1)
    return fx


def com3MDELMMain(cfgPath, cfgSetup):
    """
    """
    log.info("Running com3MDELMMain model on DEM \n \t %s \n",
             cfgPath['demSource'])
    resampleResolution = float(cfgSetup['distance'])
    xx0 = float(cfgSetup['xx0'])
    yy0 = float(cfgSetup['yy0'])
    m0 = float(cfgSetup['m0'])
    mmin = float(cfgSetup['mmin'])
    v20 = float(cfgSetup['v20'])
    xx0 = float(cfgSetup['xx0'])
    mu = float(cfgSetup['mu'])
    g = float(cfgSetup['g'])
    log.info("Reading and preparing DEM")
    Points, csz = prepareDEM(cfgPath['demSource'], resampleResolution)

    # determine start indices
    X = Points['x']
    Y = Points['y']
    Z = Points['z']
    ny, nx = np.shape(X)
    indx0 = np.argmin(abs(X[1, :]-xx0))
    indy0 = np.argmin(abs(Y[:, 1]-yy0))
    z0 = Z[indy0, indx0]
    x0 = X[indy0, indx0]
    y0 = Y[indy0, indx0]
    xPath = np.array([x0])
    yPath = np.array([y0])
    zPath = np.array([z0])
    sPath = np.array([0])

    log.info("Running with standard start point X0=%3.2f m, Y0=%3.2f m at an altitude of Z0=%3.2f m \n", x0, y0, z0)
    # initialize list that can grow
    xyDonIndList = np.array([indy0, indx0])
    xyDonIndList = xyDonIndList.reshape(2, 1)

    # Donor cellsize
    D = np.zeros((ny, nx))

    # velocity squared
    V2 = np.zeros((ny, nx))
    # mass
    MFlow = np.zeros((ny, nx))
    Mstep = np.zeros((ny, nx))
    Matrest = np.zeros((ny, nx))
    peakMass = np.zeros((ny, nx))
    # distance measurement
    S = np.zeros((ny, nx))
    Sglob = np.zeros((ny, nx))

    Lloc = np.zeros((ny, nx))
    Lglob = np.zeros((ny, nx))

    Ekin = np.zeros((ny, nx))
    Epot = np.zeros((ny, nx))

    # INITIAL CONDITIONS
    MFlow[indy0, indx0] = m0
    V2[indy0, indx0] = v20
    D[indy0, indx0] = np.nan
    # get initial conditions
    ekin = 0.5*m0*v20
    epot = m0*g*z0
    Ekin[indy0, indx0] = ekin
    Epot[indy0, indx0] = epot
    V2Path = np.array([v20])
    EkinPath = np.array([ekin])
    EpotPath = np.array([epot])

    log.info("conservation properties iteration nb %03.0f : s= %3.2f m, v2= %3.2f (m/s)^2, m= %3.2f kg, ekin=  %3.2f J ",
             0, sPath[-1], v20, m0, ekin)

    # global estimate of shortest distances
    Sglob = np.sqrt(np.square(X-x0)+np.square(Y-y0))
    Lglob = np.sqrt(np.square(X-x0)+np.square(Y-y0)+np.square(Z-z0))

    # define how to check all neighbors
    # ---------------------------------------
    # |----NW-----||-----N-----||----NE-----|
    # | i_dy0 + 1 || i_dy0 + 1 || i_dy0 + 1 |
    # | i_dx0 - 1 || i_dx0 + 0 || i_dx0 + 1 |
    # ---------------------------------------
    # |-----W-----||----  -----||-----E-----|
    # | i_dy0 + 0 ||   i_dy0   || i_dy0 + 0 |
    # | i_dx0 - 1 ||   i_dx0   || i_dx0 + 1 |
    # ---------------------------------------
    # |----SW-----||-----S-----||----SE-----|
    # | i_dy0 - 1 || i_dy0 - 1 || i_dy0 - 1 |
    # | i_dx0 - 1 || i_dx0 + 0 || 1i_nxdx0 + 1 |
    # ---------------------------------------
    indN = np.array([-1, 0, 1])
    indgrid = np.ix_(indN+indy0, indN+indx0)
    iter = 0
    iterate = True
    while iterate:
        # initialize step Fields: S, V2, M and Ekin
        Sstep = np.ones((ny, nx))*np.nan
        Lstep = np.ones((ny, nx))*np.nan
        V2step = np.zeros((ny, nx))
        Mpstep = Mstep
        Mstep = np.zeros((ny, nx))
        Ekinstep = np.zeros((ny, nx))
        R = np.zeros((ny, nx))
        xyDonIndListP = xyDonIndList
        nd = np.shape(xyDonIndListP)[1]
        # delete list of future donors
        xyDonIndList = np.empty((2, 0), int)

        iter = iter + 1

        ######################################################
        # Calc S and L and Estimate possible velocity and mass
        ######################################################
        log.debug("Distance, velocity and mass estimation loop")
        for id in range(nd):
            indx = xyDonIndListP[1, id]
            indy = xyDonIndListP[0, id]
            if (indx == 0) or (indx == nx-1) or (indy == 0) or (indy == ny-1):
                log.warning("\n +++ Approaching border of computational domain - computation abborted +++ \n")
                iterate = False
                break
            indgrid = np.ix_(indN+indy, indN+indx)
            x = X[indgrid]
            y = Y[indgrid]
            z = Z[indgrid]
            d = D[indgrid]
            sn = Sstep[indgrid]
            # at step p (previous one, iter-1)
            md = MFlow[indy, indx]
            v2d = V2[indy, indx]
            sd = S[indy, indx]
            ld = Lloc[indy, indx]

            # estimate local distance, V2 and GradV2
            DLLoc, DSLoc, gradV2, v2rest = calcV2(x, y, z, v2d, d, g, mu)

            sloc = sd + DSLoc
            lloc = ld + DLLoc
            # find the shortest path to the cell
            Sstep[indgrid] = np.fmin(Sstep[indgrid], sloc)
            Lstep[indgrid] = np.fmin(Lstep[indgrid], lloc)
            # previous donors forbiden
            Sstep[indgrid] = np.where(np.isnan(d), np.nan, Sstep[indgrid])
            Lstep[indgrid] = np.where(np.isnan(d), np.nan, Lstep[indgrid])
            # Sstep[indgrid] = np.where(Sglob[indgrid] < Sglob[indy, indx],
            #                           np.nan, Sstep[indgrid])

            # update velocity
            V2step[indgrid] = v2rest

            # create weights for mass spreading calculation

            ind = np.where(np.isnan(gradV2))
            weight = wFunk(gradV2)
            weight[ind] = 0
            sumWeight = np.nansum(weight)
            # if there are some options to move
            # estimate the mass transfer
            if sumWeight > 0:
                mn = weight/sumWeight*md
                Mstep[indgrid] = Mstep[indgrid] + mn
                Mstep[indy, indx] = Mstep[indy, indx] - md

        if iterate == False:
            log.warning("\n +++ Approaching border of computational domain - computation abborted +++ \n")
            break
        ######################
        # Mass correction loop
        ######################
        log.debug("Mass correction loop")
        R = np.zeros((ny, nx))
        Smes = np.zeros((ny, nx))
        Msteppot = Mstep
        V2steppot = V2step
        V2step = np.zeros((ny, nx))
        MV2step = np.zeros((ny, nx))
        Mstep = np.zeros((ny, nx))
        # delete list because it will get shorter
        xyDonIndList = np.empty((2, 0), int)
        for id in range(nd):
            indx = xyDonIndListP[1, id]
            indy = xyDonIndListP[0, id]
            indgrid = np.ix_(indN+indy, indN+indx)
            x = X[indgrid]
            y = Y[indgrid]
            z = Z[indgrid]
            d = D[indgrid]
            sn = Sstep[indgrid]
            Rtmp = np.zeros((ny, nx))

            # at step p (previous one, iter-1)
            md = MFlow[indy, indx]
            v2d = V2[indy, indx]
            sd = S[indy, indx]

            # estimate local distance, V2 and GradV2
            DLLoc, DSLoc, gradV2, v2rest = calcV2(x, y, z, v2d, d, g, mu)

            sloc = sd + DSLoc

            gradV2Inter = gradV2
            v2Inter = v2rest
            gradV2 = np.where(Msteppot[indgrid] <= mmin, np.nan, gradV2)
            v2rest = np.where(Msteppot[indgrid] <= mmin, 0.0, v2rest)

            # create weights for mass spreading calculation
            ind = np.where(np.isnan(gradV2))
            weight = wFunk(gradV2)
            weight[ind] = 0
            sumWeight = np.nansum(weight)
            if sumWeight > 0.0:
                # spreading
                mn = weight/sumWeight*md
                Mstep[indgrid] = Mstep[indgrid] + mn
                Mstep[indy, indx] = Mstep[indy, indx] - md
                MV2step[indgrid] = MV2step[indgrid] + mn*v2rest
                Rtmp[indgrid] = np.where(((Mstep[indgrid] > 0) &
                                         (R[indgrid] == 0)), 1, Rtmp[indgrid])
                R[indgrid] = np.where(Mstep[indgrid] > 0, np.nan, R[indgrid])
                newD = np.where(Rtmp[indgrid] == 1)
                newDInd = np.array([indgrid[0].reshape(1, 3)[0][newD[0]], (indgrid[1])[0][newD[1]]])
                xyDonIndList = np.append(xyDonIndList, newDInd, axis=1)
                Smes[indgrid] = Smes[indgrid] + sloc*mn  # *v2rest
            else:
                # point motion
                # find index of max.. of quantity
                ind = np.where(np.isnan(gradV2Inter))
                toMaximize = wFunk(gradV2Inter)
                toMaximize[ind] = 0
                log.debug("no spreading for this donor:%d", id)
                if np.nanmax(toMaximize) > 0:
                    Ind = np.where(toMaximize == np.nanmax(toMaximize))
                    maxIndCol = Ind[0][0]
                    maxIndRow = Ind[1][0]
                    log.debug("point motion for this donor:%d", id)
                    if np.shape(Ind)[1] > 1:
                        test = Sglob[indgrid]
                        Indmax = np.where(test[Ind] == np.nanmax(test[Ind]))
                        maxIndCol = Ind[0][Indmax[0][0]]
                        maxIndRow = Ind[1][Indmax[0][0]]
                    indxn = indx + (maxIndRow - 1)
                    indyn = indy + (maxIndCol - 1)

                    newInd = np.array([indyn, indxn]).reshape(2, 1)
                    xyDonIndList = np.append(xyDonIndList, newInd, axis=1)
                    Mstep[indyn, indxn] = Mstep[indyn, indxn] + md
                    Mstep[indy, indx] = Mstep[indy, indx] - md
                    MV2step[indyn, indxn] = MV2step[indyn, indxn] + md*v2Inter[maxIndCol, maxIndRow]
                    R[indyn, indxn] = np.nan
                    Smes[indyn, indxn] = Smes[indyn, indxn] + sloc[maxIndCol, maxIndRow]*md  # *v2Inter[maxIndCol, maxIndRow]
                else:
                    log.debug("no motion for this donor:%d", id)
                    Matrest[indy, indx] = Matrest[indy, indx] + md
                    Mstep[indy, indx] = Mstep[indy, indx] - md
                    R[indy, indx] = np.nan

        if np.shape(xyDonIndList)[1]>0:
            # update donor cells
            D = D + R
            # update velocity and mass
            V2step = np.zeros((ny, nx))
            MFlow = MFlow + Mstep
            peakMass = np.maximum(peakMass, MFlow)
            MM = np.where(np.isnan(MFlow), 0.0, MFlow)
            IndNull = np.where(MM > 0.0)
            V2step[IndNull] = np.divide(MV2step[IndNull], MM[IndNull])
            V2step = np.where(np.isnan(V2step), 0.0, V2step)
            V2 = V2 + V2step

            # update distance
            # IndNull = np.where(MV2step > 0.0)
            # np.divide(Smes[IndNull], MV2step[IndNull])
            Smes[IndNull] = np.divide(Smes[IndNull], MM[IndNull])
            # S = np.where(MFlow > 0.0, Smes, S)
            S = np.where(MFlow > 0.0, Smes, S)  # + np.where(M<=0.0, S, S)

            # update energy
            Ekinstep = 0.5*MV2step
            Epotstep = MFlow*g*Z
            Etotstep = Ekinstep + Epotstep
            Ekin = Ekin + Ekinstep
            Epot = Epot + Epotstep
            Etot = Ekin + Epot
            EkinSumCoE = np.sum(Ekinstep)
            EpotSumCoE = np.sum(Epotstep)
            EtotSumCoE = EkinSumCoE + EpotSumCoE

            pond = Ekinstep
            pondSum = EkinSumCoE
            zcoE = np.nansum(pond*Z)/pondSum
            # calculate XY locations of computation step
            xcoE = np.sum(pond*X)/pondSum
            ycoE = np.sum(pond*Y)/pondSum
            # calculate global distance of coE coordinates
            scoE = np.sum(pond*Smes)/pondSum
            v2coE = np.sum(pond*V2)/pondSum

            # append x, y and update distance
            xPath = np.append(xPath, xcoE)
            yPath = np.append(yPath, ycoE)
            zPath = np.append(zPath, zcoE)
            sPath = np.append(sPath, scoE)
            V2Path = np.append(V2Path, v2coE)

            # update energy
            EkinPath = np.append(EkinPath, EkinSumCoE)
            EpotPath = np.append(EpotPath, EpotSumCoE)

            log.info("conservation properties iteration nb %03.0f : s= %3.2f m, v2= %3.2f (m/s)^2, m= %3.2f kg, ekin=  %3.2f J ",
            iter, scoE, v2coE, np.sum(MFlow + Matrest), EkinSumCoE)
        else:
            iterate = False

        # if iter>0:
        #     fig = plt.figure(figsize=(2*figW, figH))
        #     ax1 = plt.subplot(311)
        #     cmap = cmapPlasma
        #     cmap.set_under(color='w')
        #     x = Points['x'][0, :]
        #     y = Points['y'][:, 0]
        #     im0 = NonUniformImage(ax1, extent=[x.min(), x.max(), y.min(), y.max()], cmap=cmap)
        #     im0.set_clim(vmin=0.000000001)
        #     # im.set_interpolation('bilinear')
        #     im0.set_data(x, y, MFlow)
        #     ref1 = ax1.images.append(im0)
        #     cbar = ax1.figure.colorbar(im0, ax=ax1, use_gridspec=True)
        #     cbar.ax.set_ylabel('Kinetic Energy [J]')
        #     ax1.title.set_text('Energy')
        #
        #     ax1.plot(xPath, yPath, 'k--', label='avapath')
        #     ax1.set_xlim([x.min(), x.max()])
        #     ax1.set_ylim([y.min(), y.max()])
        #     ax1.set_xlabel(r'$x\;[m]$')
        #     ax1.set_ylabel(r'$y\;[m]$')
        #     ax1.legend()
        #
        #     ax2 = plt.subplot(312)
        #     im = NonUniformImage(ax2, extent=[x.min(), x.max(), y.min(), y.max()], cmap=cmap)
        #     im.set_clim(vmin=0.000000001)
        #     cmap.set_under(color='w')
        #     # im.set_interpolation('bilinear')
        #     im.set_data(x, y, Ekin)
        #     ref1 = ax2.images.append(im)
        #     cbar = ax2.figure.colorbar(im, ax=ax2, use_gridspec=True)
        #     cbar.ax.set_ylabel('distance [m]')
        #     ax2.title.set_text('S')
        #
        #     ax2.plot(xPath, yPath, 'k--', label='avapath')
        #     ax2.set_xlim([x.min(), x.max()])
        #     ax2.set_ylim([y.min(), y.max()])
        #     ax2.set_xlabel(r'$x\;[m]$')
        #     ax2.set_ylabel(r'$y\;[m]$')
        #     ax2.legend()
        #
        #     ax3 = plt.subplot(313)
        #     im = NonUniformImage(ax3, extent=[x.min(), x.max(), y.min(), y.max()], cmap=cmap)
        #     im.set_clim(vmin=0.000000001)
        #     # im.set_interpolation('bilinear')
        #     im.set_data(x, y, D)
        #     ref1 = ax3.images.append(im)
        #     cbar = ax3.figure.colorbar(im, ax=ax3, use_gridspec=True)
        #     cbar.ax.set_ylabel('distance [m]')
        #     ax3.title.set_text('S')
        #     ax3.set_xlim([x.min(), x.max()])
        #     ax3.set_ylim([y.min(), y.max()])
        #
        #     plt.show()
        #
        # input("Press Enter to continue...")

    # plotting
    avapath = {}
    avapath['x'] = np.array([xx0, 4000])
    avapath['y'] = np.array([yy0, 0])
    dem = IOf.readRaster(cfgPath['demSource'])
    AvaProfile, projPoint = geoTrans.prepareLine(dem, avapath, distance=10, Point=None)
    # plt.close(fig)
    # plt.close(fig1)
    fig = plt.figure(figsize=(2*figW, figH))
    ax1 = plt.subplot(211)
    cmap = cmapPlasma
    cmap.set_under(color='w')
    x = Points['x'][0, :]
    y = Points['y'][:, 0]
    im0 = NonUniformImage(ax1, extent=[x.min(), x.max(), y.min(), y.max()], cmap=cmap)
    im0.set_clim(vmin=0.000000001)
    # im.set_interpolation('bilinear')
    im0.set_data(x, y, Ekin)
    ref1 = ax1.images.append(im0)
    cbar = ax1.figure.colorbar(im0, ax=ax1, use_gridspec=True)
    cbar.ax.set_ylabel('Kinetic Energy [J]')
    ax1.title.set_text('Energy')
    Cp = ax1.contour(X, Y, Z, levels=10, colors='k')
    # ax1.clabel(Cp, colors='k', inline=1, fontsize=5)
    ax1.plot(xPath, yPath, 'k--', label='avapath')
    # ax1.axis('equal')
    ax1.set_xlim([x.min(), x.max()])
    ax1.set_ylim([y.min(), y.max()])
    ax1.set_xlabel(r'$x\;[m]$')
    ax1.set_ylabel(r'$y\;[m]$')
    ax1.legend()

    ax2 = plt.subplot(212)
    # ax2.title.set_text('Path on DEM')
    # cmap = cmapDEM
    # im = NonUniformImage(ax2, extent=[x.min(), x.max(), y.min(), y.max()], cmap=cmap)
    # # im.set_interpolation('bilinear')
    # im.set_clim(vmin=-0.000000001)
    # im.set_data(x, y, Z)
    # ref0 = ax2.images.append(im)
    # cbar = ax2.figure.colorbar(im, ax=ax2, use_gridspec=True)
    # cbar.ax.set_ylabel('Altitude [m]')
    #
    # ax2.plot(xPath, yPath, 'k', label='avapath')
    # ax2.set_xlim([x.min(), x.max()])
    # ax2.set_ylim([y.min(), y.max()])
    # ax2.set_xlabel(r'$x\;[m]$')
    # ax2.set_ylabel(r'$y\;[m]$')
    # ax2.legend()

    ax2.plot(sPath, zPath, 'k-', label='Avalanche profile')
    # ax2.plot(AvaProfile['s'], AvaProfile['z'], 'r-', label='Avalanche profile')
    f = zPath[0] - mu * sPath
    ax2.plot(sPath, f, '-', color='b', label='AlphaLine')
    Zene = zPath + V2Path/(2*g)
    scat = ax2.scatter(sPath, Zene, marker='s', cmap=cmap, s=2*ms, c= EkinPath, label='Total energy height')
    cbar2 = ax2.figure.colorbar(scat, ax=ax2, use_gridspec=True)
    cbar2.ax.set_ylabel('Kinetic Energy [J]')
    ax2.axvline(x=sPath[-1], color='k',
    linewidth=1, linestyle='-.', label='Run out point')
    ax2.set_xlabel('s [m]', fontsize=fs)
    ax2.set_ylabel('Altitude [m]', fontsize=fs)
    # ax1.set_title('title subplot 1')
    ax2.legend()
    fig.tight_layout()

    fig1 = plt.figure(figsize=(2*figW, figH))
    ax1 = plt.subplot(311)
    cmap = cmapPlasma
    cmap.set_under(color='w')
    x = Points['x'][0, :]
    y = Points['y'][:, 0]
    im0 = NonUniformImage(ax1, extent=[x.min(), x.max(), y.min(), y.max()], cmap=cmap)
    im0.set_clim(vmin=0.000000001)
    # im.set_interpolation('bilinear')
    im0.set_data(x, y, peakMass)
    ref1 = ax1.images.append(im0)
    cbar = ax1.figure.colorbar(im0, ax=ax1, use_gridspec=True)
    cbar.ax.set_ylabel('M [kg]')
    # ax1.title.set_text('M')
    Cp1 = ax1.contour(X, Y, Z, levels=10, colors='k')
    # ax1.clabel(Cp, colors='k', inline=1, fontsize=5)
    ax1.plot(xPath, yPath, 'k--', label='avapath')
    # ax1.axis('equal')
    ax1.set_xlim([x.min(), x.max()])
    ax1.set_ylim([y.min(), y.max()])
    ax1.set_xlabel(r'$x\;[m]$')
    ax1.set_ylabel(r'$y\;[m]$')
    ax1.legend()

    ax2 = plt.subplot(312)
    cmap = cmapPlasma
    cmap.set_under(color='w')
    x = Points['x'][0, :]
    y = Points['y'][:, 0]
    im2 = NonUniformImage(ax2, extent=[x.min(), x.max(), y.min(), y.max()], cmap=cmap)
    im2.set_clim(vmin=0.000000001)
    # im.set_interpolation('bilinear')
    im2.set_data(x, y, np.sqrt(V2))
    ref2 = ax2.images.append(im2)
    cbar2 = ax2.figure.colorbar(im2, ax=ax2, use_gridspec=True)
    cbar2.ax.set_ylabel('V [m/s^2]')
    # ax2.title.set_text('V')
    Cp2 = ax2.contour(X, Y, Z, levels=10, colors='k')
    # ax1.clabel(Cp, colors='k', inline=1, fontsize=5)
    ax2.plot(xPath, yPath, 'k--', label='avapath')
    # ax1.axis('equal')
    ax2.set_xlim([x.min(), x.max()])
    ax2.set_ylim([y.min(), y.max()])
    ax2.set_xlabel(r'$x\;[m]$')
    ax2.set_ylabel(r'$y\;[m]$')
    ax2.legend()

    ax3 = plt.subplot(313)
    cmap = cmapPlasma
    cmap.set_under(color='w')
    x = Points['x'][0, :]
    y = Points['y'][:, 0]
    im3 = NonUniformImage(ax3, extent=[x.min(), x.max(), y.min(), y.max()], cmap=cmap)
    im3.set_clim(vmin=0.000000001)
    # im.set_interpolation('bilinear')
    im3.set_data(x, y, np.sqrt(S))
    ref3 = ax3.images.append(im3)
    cbar3 = ax3.figure.colorbar(im3, ax=ax3, use_gridspec=True)
    cbar3.ax.set_ylabel('S [m]')
    # ax3.title.set_text('S')
    Cp3 = ax3.contour(X, Y, Z, levels=10, colors='k')
    # ax1.clabel(Cp, colors='k', inline=1, fontsize=5)
    ax3.plot(xPath, yPath, 'k--', label='avapath')
    # ax1.axis('equal')
    ax3.set_xlim([x.min(), x.max()])
    ax3.set_ylim([y.min(), y.max()])
    ax3.set_xlabel(r'$x\;[m]$')
    ax3.set_ylabel(r'$y\;[m]$')
    ax3.legend()

    fig1.tight_layout()
    plt.show()
