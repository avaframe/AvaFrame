"""
    Main logic for MDELM computational module

    This file is part of Avaframe.
"""

import os
import glob
import pickle
import logging
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.image import NonUniformImage
from matplotlib.collections import LineCollection

# Local imports
import avaframe.in2Trans.geoTrans as geoTrans
import avaframe.in3Utils.ascUtils as IOf
from avaframe.out3SimpPlot.plotSettings import *

# create local logger
log = logging.getLogger(__name__)
debugPlot = False


def com3MDELMMain(cfgPath, cfgSetup):
    """
    """
    log.info("Running com3MDELMMain model on DEM \n \t %s \n",
             cfgPath['demSource'])
    resampleResolution = float(cfgSetup['distance'])
    xx0 = float(cfgSetup['xx0'])
    yy0 = float(cfgSetup['yy0'])
    m0 = float(cfgSetup['m0'])
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
    xyIndList = np.array([indy0, indx0])
    xyIndList = xyIndList.reshape(1, 2)
    # get initial conditions
    M = np.zeros((1,))
    M[0] = m0
    Ekin = -np.ones((ny, nx))
    Epot = np.zeros((ny, nx))
    ekin = 0.5*m0*v20
    epot = m0*g*z0
    Ekin[indy0, indx0] = ekin
    Epot[indy0, indx0] = epot
    V2Path = np.array([v20])
    EkinPath = np.array([ekin])
    EpotPath = np.array([epot])
    # coE = [i_iteration,x_coE,y_coE,z_coE, s_coE, ,v2_coE ,m_coE, ekin_sum_coE]
    log.info("conservation properties iteration nb %03.0f : s= %3.2f m, v2= %3.2f (m/s)^2, m= %3.2f kg, ekin=  %3.2f J ",
             0, sPath[-1], v20, m0, ekin)

    # Donor cellsize
    D = np.zeros((ny, nx))

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
    v2p = v20
    iterate = True
    while iterate:
        indx = xyIndList[iter, 1]
        indy = xyIndList[iter, 0]
        D[indy, indx] = np.nan
        iter = iter + 1
        if (indx == 0) or (indx == nx-1) or (indy == 0) or (indy == ny-1):
            log.info("\n +++ Approaching border of computational domain - computation abborted +++ \n")
            iterate = False
            break

        indgrid = np.ix_(indN+indy, indN+indx)
        z = Z[indgrid]
        x = X[indgrid]
        y = Y[indgrid]
        d = D[indgrid]
        # at step n
        xd = x[1, 1]
        yd = y[1, 1]
        zd = z[1, 1]
        if iter<2:
            xdm2 = xPath[-1]
            ydm2 = yPath[-1]
            zdm2 = zPath[-1]
            v2pm2 = V2Path[-1]
        else:
            xdm2 = xPath[-2]
            ydm2 = yPath[-2]
            zdm2 = zPath[-2]
            v2pm2 = V2Path[-2]
        DeltaL = np.sqrt(np.square(x-xd)+np.square(y-yd)+np.square(z-zd))
        DeltaLm2 = np.sqrt(np.square(x-xdm2)+np.square(y-ydm2)+np.square(z-zdm2))
        DeltaS = np.sqrt(np.square(x-xd)+np.square(y-yd))
        bx = (z[1, 2]-zd)/(x[1, 2]-xd)
        bx2 = bx*bx
        by = (z[2, 1]-zd)/(y[2, 1]-yd)
        by2 = by*by
        c = np.sqrt(1+bx2+by2)
        print('deltaL')
        print(DeltaL)
        V2next = v2p - 2*g*((z-zd) + mu*DeltaL/c)
        V2nextm2 = v2pm2 - 2*g*((z-zdm2) + mu*DeltaLm2/c)
        V2nextJT = v2p - 2*g*((z-zd) + mu*DeltaS)

        ekin = 0.5*m0*V2next
        epot = m0*g*z
        etot = ekin + epot
        deltaEtot = (etot - etot[1][1])/DeltaL
        # find index of max.. of quantity
        toMaximize = V2next - d
        toMaximize = np.where(V2next<0, np.nan, toMaximize)

        print('V2next')
        print(V2next)
        print('etot')
        print(etot)
        print('deltaEtot')
        print(deltaEtot)
        print('V2nextm2')
        print(V2nextm2)
        Ind = np.where(toMaximize == np.nanmax(toMaximize))
        try:
            maxIndCol = Ind[0][0]
            maxIndRow = Ind[1][0]
        except IndexError:
            maxIndCol = 1
            maxIndRow = 1
        if np.shape(Ind)[1]>1:
            maxIndCol = Ind[0][1]
            maxIndRow = Ind[1][1]
        indxn = indx + (maxIndRow - 1)
        indyn = indy + (maxIndCol - 1)

        newInd = np.array([indyn, indxn]).reshape(1, 2)
        xyIndList = np.append(xyIndList, newInd, axis=0)
        if np.nanmax(V2next-d) <= 0:
            iterate = False
            v2n = 0
            v2p = v2n
        else:
            v2n = V2next[maxIndCol, maxIndRow]
            v2p = v2n
            # append x, y and update distance
            xPath = np.append(xPath, x[maxIndCol, maxIndRow])
            yPath = np.append(yPath, y[maxIndCol, maxIndRow])
            zPath = np.append(zPath, z[maxIndCol, maxIndRow])
            sPath = np.append(sPath, sPath[-1] + DeltaS[maxIndCol, maxIndRow])
            V2Path = np.append(V2Path, v2n)

            # update energy
            Ekin[indyn, indxn] = ekin[maxIndCol, maxIndRow]
            Epot[indyn, indxn] = epot[maxIndCol, maxIndRow]
            EkinPath = np.append(EkinPath, ekin[maxIndCol, maxIndRow])
            EpotPath = np.append(EpotPath, ekin[maxIndCol, maxIndRow])
            log.info("conservation properties iteration nb %03.0f : s= %3.2f m, v2= %3.2f (m/s)^2, m= %3.2f kg, ekin=  %3.2f J ",
                      iter, sPath[-1], v2n, m0, ekin[maxIndCol, maxIndRow])



    #------------------------------------------------------#
    # plotting
    avapath = {}
    avapath['x'] = np.array([xx0, 4000])
    avapath['y'] = np.array([yy0, 0])
    dem = IOf.readRaster(cfgPath['demSource'])
    AvaProfile, projPoint = geoTrans.prepareLine(dem, avapath, distance=10, Point=None)
    # plt.close(fig)
    # plt.close(fig1)
    fig = plt.figure(figsize=(2*figW, figH), dpi=figReso)
    ax1 = plt.subplot(121)
    cmap = cmapPlasma
    cmap.set_under(color='w')
    x = Points['x'][0, :]
    y = Points['y'][:, 0]
    im0 = NonUniformImage(ax1, extent=[x.min(), x.max(), y.min(), y.max()], cmap=cmap)
    im0.set_clim(vmin=-0.000000001)
    # im.set_interpolation('bilinear')
    im0.set_data(x, y, Ekin)
    ref1 = ax1.images.append(im0)
    cbar = ax1.figure.colorbar(im0, ax=ax1, use_gridspec=True)
    cbar.ax.set_ylabel('Kinetic Energy [J]')
    ax1.title.set_text('Energy')
    ax1.set_xlim([x.min(), x.max()])
    ax1.set_ylim([y.min(), y.max()])
    ax1.set_xlabel(r'$x\;[m]$')
    ax1.set_ylabel(r'$y\;[m]$')

    ax2 = plt.subplot(122)
    ax2.title.set_text('Path on DEM')
    cmap = cmapDEM
    im = NonUniformImage(ax2, extent=[x.min(), x.max(), y.min(), y.max()], cmap=cmap)
    # im.set_interpolation('bilinear')
    im.set_clim(vmin=-0.000000001)
    im.set_data(x, y, Z)
    ref0 = ax2.images.append(im)
    cbar = ax2.figure.colorbar(im, ax=ax2, use_gridspec=True)
    cbar.ax.set_ylabel('Altitude [m]')

    ax2.plot(xPath, yPath, 'k', label='avapath')
    ax2.set_xlim([x.min(), x.max()])
    ax2.set_ylim([y.min(), y.max()])
    ax2.set_xlabel(r'$x\;[m]$')
    ax2.set_ylabel(r'$y\;[m]$')
    ax2.legend()

    fig.tight_layout()

    fig1 = plt.figure(figsize=(figW, figH), dpi=figReso)
    ax1 = plt.subplot(111)
    cmap = cmapPlasma
    ax1.plot(sPath, zPath, 'k-', linewidth=lw/2, label='Avalanche profile')
    # ax1.plot(AvaProfile['s'], AvaProfile['z'], 'r-', linewidth=lw/2, label='Avalanche profile')
    f = zPath[0] - mu * sPath
    ax1.plot(sPath, f, '-', color='b', linewidth=lw/2, label='AlphaLine')
    Zene = zPath + V2Path/(2*g)
    scat = ax1.scatter(sPath, Zene, marker='s', cmap=cmap, s=2*ms, c= Zene, label='Total energy height')
    # points = np.array([sPath, Zene]).T.reshape(-1, 1, 2)
    # segments = np.concatenate([points[:-1], points[1:]], axis=1)
    # norm = plt.Normalize(Zene.min(), Zene.max())
    # lc = LineCollection(segments, cmap=cmap, norm=norm)
    # lc.set_array(Zene)
    # lc.set_linewidth(2)
    # line = ax1.add_collection(lc)
    cbar = ax1.figure.colorbar(scat, ax=ax1, use_gridspec=True)
    cbar.ax.set_ylabel('Kinetic Energy [J]')
    # ax1.colorbar()
    ax1.axvline(x=sPath[-1], color='k',
    linewidth=1, linestyle='-.', label='Run out point')
    ax1.set_xlabel('s [m]', fontsize=fs)
    ax1.set_ylabel('Altitude [m]', fontsize=fs)
    # ax1.set_title('title subplot 1')
    ax1.legend()

    plt.show()


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
    print(np.min(xgrid), np.max(xgrid), header.xllcorner, header.xllcorner+(ncols-1)*cellsize)
    ygrid = np.linspace(header.yllcorner, header.yllcorner+(m-1)*csz, m)
    print(np.min(ygrid), np.max(ygrid), header.yllcorner, header.yllcorner+(nrows-1)*cellsize)
    PointsX, PointsY = np.meshgrid(xgrid, ygrid)
    Points = {}
    Points['x'] = PointsX
    Points['y'] = PointsY

    Points, itot, ioob = geoTrans.projectOnRasterVect(dem, Points, interp='nearest')

    fig = plt.figure(figsize=(figW, figH), dpi=figReso)
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
    plt.close(fig)
    # plt.show()

    return Points, csz

def com3MDELMMain2(cfgPath, cfgSetup):
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
    R = np.zeros((ny, nx))

    # velocity squared
    V2_glob = np.zeros((ny, nx))
    V2 = np.zeros((ny, nx))
    # mass
    M = np.zeros((ny, nx))
    Mstep = np.zeros((ny, nx))
    Mpstep = np.zeros((ny, nx))
    Marrest = np.zeros((ny, nx))
    Mflow = np.zeros((ny, nx))
    mstepsum = 0;
    # distance measurement
    Sloc = np.zeros((ny, nx))
    Sglob = np.zeros((ny, nx))

    Lloc = np.zeros((ny, nx))
    Lglob = np.zeros((ny, nx))

    Ekin = np.zeros((ny, nx))
    #Ekin = -np.ones((ny, nx))
    Epot = np.zeros((ny, nx))

    # INITIAL CONDITIONS
    M[indy0, indx0] = m0
    MFlow = M;
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
    # coE = [i_iteration,x_coE,y_coE,z_coE, s_coE, ,v2_coE ,m_coE, ekin_sum_coE]
    log.info("conservation properties iteration nb %03.0f : s= %3.2f m, v2= %3.2f (m/s)^2, m= %3.2f kg, ekin=  %3.2f J ",
             0, sPath[-1], v20, m0, ekin)
    # WE CAN CALCULATE THEM HERE OR IN THE LOOPS -> EVERYTHING MOVED TO LOOPS NOW
    #global estimate of shortest distan0ceS
    Sglob = np.sqrt(np.square(X-x0)+np.square(Y-y0))
    Lglob = np.sqrt(np.square(X-x0)+np.square(Y-y0)+np.square(Z-z0))
    #global estimate of v2
    V2pot = 2*g*(z0 - Z)
    V2dis = 2*g*(mu * Lglob)
    V2glob = np.ones((ny, nx))*v20 + (V2pot-V2dis)
    #there cannot be negative velofigurecity
    V2glob = np.where(V2glob<0, 0, V2glob)

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
        #initialize step Fields: S, V2, M and Ekin
        Sstep = np.ones((ny, nx))*np.nan
        Lstep = np.ones((ny, nx))*np.nan
        V2step = np.zeros((ny, nx))
        Mpstep = Mstep;
        Mstep = np.zeros((ny, nx))
        Ekin_step = np.zeros((ny, nx))
        R = np.zeros((ny, nx))
        xyDonIndListP = xyDonIndList
        nd = np.shape(xyDonIndListP)[1]
        print(xyDonIndListP)
        print(np.shape(xyDonIndListP))
        #delete list of future donors
        xyDonIndList = np.empty((2, 0), int)

        iter = iter + 1
        ####################
        # Measuring distance
        ####################
        log.info("Distance measurement loop")
        for id in range(nd):
            indx = xyDonIndListP[1, id]
            indy = xyDonIndListP[0, id]
            if (indx == 0) or (indx == nx-1) or (indy == 0) or (indy == ny-1):
                log.info("\n +++ Approaching border of computational domain - computation abborted +++ \n")
                iterate = False
                break
            indgrid = np.ix_(indN+indy, indN+indx)
            x = X[indgrid]
            y = Y[indgrid]
            z = Z[indgrid]
            d = D[indgrid]
            lglob = Lglob[indgrid]
            sglob = Sglob[indgrid]

            # at step p (previous one, iter-1)
            xd = X[indy, indx]
            yd = Y[indy, indx]
            zd = Z[indy, indx]
            llocp = Lloc[indy, indx]
            slocp = Sloc[indy, indx]

            DeltaLLoc = np.sqrt(np.square(x-xd)+np.square(y-yd)+np.square(z-zd))
            DeltaSLoc = np.sqrt(np.square(x-xd)+np.square(y-yd))

            sloc = slocp + DeltaSLoc
            lloc = llocp + DeltaLLoc

            Sstep[indgrid] = np.fmin(Sstep[indgrid], sloc)
            Lstep[indgrid] = np.fmin(Lstep[indgrid], lloc)

        if iterate == False:
            log.info("\n +++ Approaching border of computational domain - computation abborted +++ \n")
            break

        #####################################
        # Estimate possible velocity and mass
        #####################################
        log.info("velocity and mass estimation loop")
        for id in range(nd):
            indx = xyDonIndListP[1, id]
            indy = xyDonIndListP[0, id]
            indgrid = np.ix_(indN+indy, indN+indx)
            x = X[indgrid]
            y = Y[indgrid]
            z = Z[indgrid]
            d = D[indgrid]
            # at step p (previous one, iter-1)
            md = M[indy, indx]
            v2d = V2[indy, indx]
            xd = X[indy, indx]
            yd = Y[indy, indx]
            zd = Z[indy, indx]

            DeltaLLoc = np.sqrt(np.square(x-xd)+np.square(y-yd)+np.square(z-zd))
            DeltaSLoc = np.sqrt(np.square(x-xd)+np.square(y-yd))

            bx = (z[1, 2]-zd)/(x[1, 2]-xd)
            bx2 = bx*bx
            by = (z[2, 1]-zd)/(y[2, 1]-yd)
            by2 = by*by
            c = np.sqrt(1+bx2+by2)

            v2rest = v2d - 2*g*((z-zd) + mu*DeltaLLoc/c)
            v2rest = np.where(np.isnan(d), 0.0, v2rest)
            v2rest = np.where(v2rest<0, 0.0, v2rest)
            V2nextJT = v2d - 2*g*((z-zd) + mu*DeltaSLoc)

            V2step[indgrid] = v2rest
            print('velocity, donor', id)
            print(V2step[indgrid])
            v2rsum = np.nansum(v2rest)
            # if there are some options to move
            if v2rsum>0:
                mn = V2step/v2rsum*md
            else:
                mn = np.zeros((ny, nx))

            Mstep[indgrid] = Mstep[indgrid] + mn[indgrid]
            Mstep[indy, indx] = Mstep[indy, indx] - md
            R[indgrid] = np.where(Mstep[indgrid]>0, 1, R[indgrid])
            newD = np.where(R[indgrid] == 1)
            newDInd = np.array([indgrid[0].reshape(1,3)[0][newD[0]],(indgrid[1])[0][newD[1]]])
            xyDonIndList = np.append(xyDonIndList, newDInd, axis=1)

        for id in range(nd):
            indx = xyDonIndListP[1, id]
            indy = xyDonIndListP[0, id]
            indgrid = np.ix_(indN+indy, indN+indx)
            print('mass recevers estimate for donor', id)
            print(Mstep[indgrid])

        ######################
        # Mass correction loop
        ######################
        log.info("Mass correction loop")
        Rpot = R
        R = np.zeros((ny, nx))
        Msteppot = Mstep
        V2steppot = V2step
        V2step = np.zeros((ny, nx))
        MV2step = np.zeros((ny, nx))
        Mstep = np.zeros((ny, nx))
        xyDonIndList = np.empty((2, 0), int) # delete list because it will get shorter
        for id in range(nd):
            indx = xyDonIndListP[1, id]
            indy = xyDonIndListP[0, id]
            indgrid = np.ix_(indN+indy, indN+indx)
            x = X[indgrid]
            y = Y[indgrid]
            z = Z[indgrid]
            d = D[indgrid]
            # at step p (previous one, iter-1)
            md = M[indy, indx]
            v2d = V2[indy, indx]
            xd = X[indy, indx]
            yd = Y[indy, indx]
            zd = Z[indy, indx]
            Rtmp = np.zeros((ny, nx))

            DeltaLLoc = np.sqrt(np.square(x-xd)+np.square(y-yd)+np.square(z-zd))
            DeltaSLoc = np.sqrt(np.square(x-xd)+np.square(y-yd))

            bx = (z[1, 2]-zd)/(x[1, 2]-xd)
            bx2 = bx*bx
            by = (z[2, 1]-zd)/(y[2, 1]-yd)
            by2 = by*by
            c = np.sqrt(1+bx2+by2)

            v2rest = v2d - 2*g*((z-zd) + mu*DeltaLLoc/c)
            v2rest = np.where(np.isnan(d), 0.0, v2rest)
            v2rest = np.where(v2rest<0, 0.0, v2rest)

            V2step[indgrid] = v2rest


            V2step[indgrid] = np.where(Msteppot[indgrid]<=mmin, 0.0, v2rest)
            v2rsum = np.sum(V2step[indgrid])
            # print('mass estimation, donor', id)
            # print(Msteppot[indgrid])
            # print('new velocity, donor', id)
            # print(V2step[indgrid])
            v2rsum = np.sum(V2step[indgrid])
            # Redistribute the mass
            if v2rsum>0:
                mn = V2step/v2rsum*md
                Mstep[indgrid] = Mstep[indgrid] + mn[indgrid]
                Mstep[indy, indx] = Mstep[indy, indx] - md
                MV2step[indgrid] = MV2step[indgrid] + Mstep[indgrid]*V2step[indgrid]
                Rtmp[indgrid] = np.where((Mstep[indgrid]>0) & (R[indgrid]==0), 1, Rtmp[indgrid])
                R[indgrid] = np.where(Mstep[indgrid]>0, 1, R[indgrid])
                newD = np.where(Rtmp[indgrid] == 1)
                newDInd = np.array([indgrid[0].reshape(1,3)[0][newD[0]],(indgrid[1])[0][newD[1]]])
                xyDonIndList = np.append(xyDonIndList, newDInd, axis=1)
            else:
                # find index of max.. of quantity
                v2rest = V2steppot[indgrid]
                toMaximize = v2rest

                print('no spreading for this donor:', id)
                print(v2rest)
                Ind = np.where(toMaximize == np.nanmax(toMaximize))
                try:
                    maxIndCol = Ind[0][0]
                    maxIndRow = Ind[1][0]
                except IndexError:
                    maxIndCol = 1
                    maxIndRow = 1
                if np.shape(Ind)[1]>1:
                    maxIndCol = Ind[0][1]
                    maxIndRow = Ind[1][1]
                indxn = indx + (maxIndRow - 1)
                indyn = indy + (maxIndCol - 1)

                newInd = np.array([indyn, indxn]).reshape(2,1)
                xyDonIndList = np.append(xyDonIndList, newInd, axis=1)
                D[indyn, indxn] = np.nan
                V2step[indyn, indxn] = V2steppot[indyn, indxn]
                Mstep[indyn, indxn] = Mstep[indyn, indxn] + md
                MV2step[indgrid] = MV2step[indgrid] + Mstep[indgrid]*V2step[indgrid]
                # print('mass estimation, donor', id)
                # print(Mstep[indgrid])
                # print('new velocity, donor', id)
                # print(V2step[indgrid])
                R[indyn, indxn] = 1;
                if np.nanmax(v2rest) <= 0:
                    iterate = False
                    v2n = 0
                else:
                    v2n = v2rest[maxIndCol, maxIndRow]
                    ekin = 0.5*m0*v2rest
                    epot = m0*g*z
                    etot = ekin + epot
        M = M + Mstep
        MM = np.where(np.isnan(M), 0.0, M)
        V2 = np.divide(MV2step,MM)
        V2 = np.where(np.isnan(V2), 0.0, V2)
        for id in range(nd):
            indx = xyDonIndListP[1, id]
            indy = xyDonIndListP[0, id]
            indgrid = np.ix_(indN+indy, indN+indx)
            print('mass recevers final for donor', id)
            print(M[indgrid])
            print('velocity recevers final for donor', id)
            print(MV2step[indgrid])

        # M = M + Mstep
        Sloc = np.where(M>0, Sstep, Sloc)
        Ekinstep = 0.5*M*V2step
        Ekin = Ekin + Ekinstep
        V2 = V2step
        D = D + R
        EkinSumCoE = np.sum(Ekinstep)
        zcoE = np.sum(Ekinstep*Z)/EkinSumCoE
        # calculate XY locations of computation step
        xcoE = np.sum(Ekinstep*X)/EkinSumCoE
        ycoE = np.sum(Ekinstep*Y)/EkinSumCoE
        # calculate global distance of coE coordinates
        # s_coE = sqrt((x_coE-x0).^2+(y_coE-y0).^2);
        # calculate coE v2 , m
        scoE   = np.nansum(Ekinstep*Sstep)/EkinSumCoE
        v2coE  = np.sum(Ekinstep*V2step)/EkinSumCoE

        # append x, y and update distance
        xPath = np.append(xPath, xcoE)
        yPath = np.append(yPath, ycoE)
        zPath = np.append(zPath, zcoE)
        sPath = np.append(sPath, scoE)
        V2Path = np.append(V2Path, v2coE)

        # update energy
        EkinPath = np.append(EkinPath, EkinSumCoE)

        log.info("conservation properties iteration nb %03.0f : s= %3.2f m, v2= %3.2f (m/s)^2, m= %3.2f kg, ekin=  %3.2f J ",
        iter, scoE, v2coE, np.sum(M), EkinSumCoE)




        input("Press Enter to continue...")

        #------------------------------------------------------#
        # plotting
        avapath = {}
        avapath['x'] = np.array([xx0, 4000])
        avapath['y'] = np.array([yy0, 0])
        dem = IOf.readRaster(cfgPath['demSource'])
        AvaProfile, projPoint = geoTrans.prepareLine(dem, avapath, distance=10, Point=None)
        # plt.close(fig)
        # plt.close(fig1)
        fig = plt.figure(figsize=(2*figW, figH), dpi=figReso)
        ax1 = plt.subplot(121)
        cmap = cmapPlasma
        cmap.set_under(color='w')
        x = Points['x'][0, :]
        y = Points['y'][:, 0]
        im0 = NonUniformImage(ax1, extent=[x.min(), x.max(), y.min(), y.max()], cmap=cmap)
        im0.set_clim(vmin=-0.000000001)
        # im.set_interpolation('bilinear')
        im0.set_data(x, y, Ekin)
        ref1 = ax1.images.append(im0)
        cbar = ax1.figure.colorbar(im0, ax=ax1, use_gridspec=True)
        cbar.ax.set_ylabel('Kinetic Energy [J]')
        ax1.title.set_text('Energy')
        ax1.set_xlim([x.min(), x.max()])
        ax1.set_ylim([y.min(), y.max()])
        ax1.set_xlabel(r'$x\;[m]$')
        ax1.set_ylabel(r'$y\;[m]$')

        ax2 = plt.subplot(122)
        ax2.title.set_text('Path on DEM')
        cmap = cmapDEM
        im = NonUniformImage(ax2, extent=[x.min(), x.max(), y.min(), y.max()], cmap=cmap)
        # im.set_interpolation('bilinear')
        im.set_clim(vmin=-0.000000001)
        im.set_data(x, y, Z)
        ref0 = ax2.images.append(im)
        cbar = ax2.figure.colorbar(im, ax=ax2, use_gridspec=True)
        cbar.ax.set_ylabel('Altitude [m]')

        ax2.plot(xPath, yPath, 'k', label='avapath')
        ax2.set_xlim([x.min(), x.max()])
        ax2.set_ylim([y.min(), y.max()])
        ax2.set_xlabel(r'$x\;[m]$')
        ax2.set_ylabel(r'$y\;[m]$')
        ax2.legend()

        fig.tight_layout()

        fig1 = plt.figure(figsize=(figW, figH), dpi=figReso)
        ax1 = plt.subplot(111)
        cmap = cmapPlasma
        ax1.plot(sPath, zPath, 'k-', linewidth=lw/2, label='Avalanche profile')
        # ax1.plot(AvaProfile['s'], AvaProfile['z'], 'r-', linewidth=lw/2, label='Avalanche profile')
        f = zPath[0] - mu * sPath
        ax1.plot(sPath, f, '-', color='b', linewidth=lw/2, label='AlphaLine')
        Zene = zPath + V2Path/(2*g)
        scat = ax1.scatter(sPath, Zene, marker='s', cmap=cmap, s=2*ms, c= EkinPath, label='Total energy height')
        # points = np.array([sPath, Zene]).T.reshape(-1, 1, 2)
        # segments = np.concatenate([points[:-1], points[1:]], axis=1)
        # norm = plt.Normalize(Zene.min(), Zene.max())
        # lc = LineCollection(segments, cmap=cmap, norm=norm)
        # lc.set_array(Zene)
        # lc.set_linewidth(2)
        # line = ax1.add_collection(lc)
        cbar = ax1.figure.colorbar(scat, ax=ax1, use_gridspec=True)
        cbar.ax.set_ylabel('Kinetic Energy [J]')
        # ax1.colorbar()
        ax1.axvline(x=sPath[-1], color='k',
        linewidth=1, linestyle='-.', label='Run out point')
        ax1.set_xlabel('s [m]', fontsize=fs)
        ax1.set_ylabel('Altitude [m]', fontsize=fs)
        # ax1.set_title('title subplot 1')
        ax1.legend()

        plt.show()
