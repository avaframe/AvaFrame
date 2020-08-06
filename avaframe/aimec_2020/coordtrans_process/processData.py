# -*- coding: utf-8 -*-

#==============================================================================
# packages
#==============================================================================
import os, sys
from time import sleep as sleep
import math
import numpy as np
import pandas as pd
import copy
import functools
from multiprocessing import Pool
from matplotlib import pyplot as plt
import matplotlib
import warnings
from PyQt5 import QtCore
from PyQt5 import QtGui
import coordtrans_process.IO_functionality as IOf
# coordtrans_process.
#==============================================================================
# bresenham
#==============================================================================
def bresenham(x0, y0, x1, y1, cs):
    """
    RASTERIZE - bresenham algorithmus - JT 2011

    input: x0, y0, x1, y1,cellsize
    output: array of x y coodinates of cells hit inbetween

    C IMPLEMENTIERUNG von http://de.wikipedia.org/wiki/Bresenham-Algorithmus
    void line(int x0, int y0, int x1, int y1)
     {
       int dx =  abs(x1-x0), sx = x0<x1 ? 1 : -1;
       int dy = -abs(y1-y0), sy = y0<y1 ? 1 : -1;
       int err = dx+dy, e2; /* error value e_xy */

       for(;;){  /* loop */
         setPixel(x0,y0);
         if (x0==x1 && y0==y1) break;
         e2 = 2*err;
         if (e2 > dy) { err += dy; x0 += sx; } /* e_xy+e_x > 0 */
         if (e2 < dx) { err += dx; y0 += sy; } /* e_xy+e_y < 0 */
       }
     }
    """
    # normalize Cellsize cs to 1
    x0 = round(x0/cs)
    x1 = round(x1/cs)
    y0 = round(y0/cs)
    y1 = round(y1/cs)

    dx = abs(x1-x0)
    dy = abs(y1-y0)
    sx = np.sign(x1-x0) # step in x direction
    sy = np.sign(y1-y0) # step in y direction
    err = dx-dy

    z = []
    while True:
        z.append([x0*cs, y0*cs])
        if x0 == x1 and y0 == y1: # if no step exists we are already there
            break
        e2 = 2*err
        if (e2 > -dy):
            err -= dy
            x0 += sx
        if (e2 < dx):
            err += dx
            y0 += sy

    return z

#==============================================================================
# path2domain
#==============================================================================
def path2domain(x, y, w, csz):
    """
    path2domain
    Usage:
        [xp, yp, ] = path2domain(x, y, w, csz)

       load xydata.txt
       width=200;
       [x,y]=thalweg2path(xydata(:,1),xydata(:,2),width);
       C=[1:length(x(:,1))]'; C=[C C];
       figure(1),clf,pcolor(x,y,C),shading flat, axis('equal'),colorbar

       Input:
           x, y:   Polyline Coordinates from file
           w:      Pfadbreite
           csz:    Zellgröße
       Output:
           xp, yp: Arrays determining a path of width w along a polyline

       Uwe Schlifkowitz/ BFW, June 2011
    """
#    Difference between x- bzw. y-Coordinates of Polyline
#    first and last  Vertex: Difference between this and the next
#    other vertices: Differenz zwischen previous and next
    dx = np.array((x[1]-x[0]))
    dy = np.array((y[1]-y[0]))
    for i in range(2, len(x)):
        dx = np.append(dx, (x[i]-x[i-2])/2.)
        dy = np.append(dy, (y[i]-y[i-2])/2.)
#        dy += [(y[i]-y[i-2])/2]
    dx = np.append(dx, x[len(x)-1]-x[len(x)-2])
    dy = np.append(dy, y[len(x)-1]-y[len(x)-2])

#    Direction of normal vector of difference,
#    a.k.a. bisecting line of angle
    d = np.arctan2(dy, dx) + math.pi/2

#    x- und y-Coordinates (left and right) of path edges,
#    total width w
#    x-KOO[left right]
    OX = np.array((x + w * np.cos(d),
                   x + w * np.cos(d + math.pi)))
#    y-KOO[left right]
    OY = np.array((y + w * np.sin(d),
                   y + w * np.sin(d + math.pi)))

#    AK 2013
#    x- und y-Coordinates (left and right) of path edges,
#    total width w + shift for area/rastersize
    # x-KOO[[left],[right]]
    OOX = np.array((x + (w+csz/2) * np.cos(d),
                    x + (w+csz/2) * np.cos(d + math.pi)))
    # y-KOO[[left],[right]]
    OOY = np.array((y + (w+csz/2) * np.sin(d),
                    y + (w+csz/2) * np.sin(d + math.pi)))

#    x-KOO[[left right], ... ,[left right]]
    OOXX = np.zeros((len(OOX[0])*2, len(OOX)))
#    y-KOO[[left right], ... ,[left right]]
    OOYY = np.zeros((len(OOY[0])*2, len(OOY)))
#        vorwärts
    OOXX[0:-1:2, 0] = OOX[0] + csz/2. * np.cos(d + 1./2*math.pi)
    OOXX[0:-1:2, 1] = OOX[1] + csz/2. * np.cos(d + 1./2*math.pi)
    OOYY[0:-1:2, 0] = OOY[0] + csz/2. * np.sin(d + 1./2*math.pi)
    OOYY[0:-1:2, 1] = OOY[1] + csz/2. * np.sin(d + 1./2*math.pi)
#        rückwärts
    OOXX[1::2, 0] = OOX[0] + csz/2. * np.cos(d + 3./2*math.pi)
    OOXX[1::2, 1] = OOX[1] + csz/2. * np.cos(d + 3./2*math.pi)
    OOYY[1::2, 0] = OOY[0] + csz/2. * np.sin(d + 3./2*math.pi)
    OOYY[1::2, 1] = OOY[1] + csz/2. * np.sin(d + 3./2*math.pi)

    return OX, OY, OOXX, OOYY

#==============================================================================
# process_data_ind
#==============================================================================
def processDataInd(fnames, polyfnames, w, with_doku, doku_name, with_damages, damages_name,
         visu, outpath, out, dpp_threshold):
    """
    process data ind
    this function is used to process the rasterdata such that it can be
    analysed with the methods for a regular grid
    data given in a regulare grid is projected on a nonuniform grid given by
    a polyline

    JT Fischer, Uwe Schlifkowitz BFW 2010-2012
    AK BFW 2014

    input: names of rasterfiles, poly names, path width
    ouput: structure{x coordinate along new raster, y coordinate, rasterdata}
    """

    statusoutput = 1
#    fnames:       full name of data file including path and extension
#    polyfnames:   file with polylinie

    aval_data = np.array(([None]*5))

    m = 0
    n = 0
    m_total = 0
    n_total = 0
    m_alt = 0

    fname = fnames[0]
    print('[PD] Data-file %s analysed' % fname)
    #read data
    header = IOf.readASCheader(fname)
    # print(header)
#    header = {}
#    with open(fname) as f:
#        for i in xrange(6):
#            (key, val) = (f.readline()).split()
#            header[key] = val
# keys: 'ncols', 'nrows', 'xllcenter', 'yllcenter','cellsize', 'NODATA_value'
#    ncols = header[1] unused
#    nrows = header[2] unused
#    xllcenter = float(header['xllcenter'])
#    yllcenter = float(header['yllcenter'])
#    cellsize = float(header['cellsize'])
    xllcenter = header.xllcenter
    yllcenter = header.yllcenter
    cellsize = header.cellsize

    rasterdata = np.flipud(np.loadtxt(fname, skiprows=6))
    rasterdata[rasterdata == header.noDataValue] = np.NaN

#    pdata = np.loadtxt(polyfnames, delimiter=',', skiprows=1) # Achtung: CoSiCa
    pdata = pd.read_csv(polyfnames, delimiter=',') #
    if 'X' in pdata.columns: x_path = np.array(pdata.X)
    if 'x' in pdata.columns: x_path = np.array(pdata.x)
    if 'Y' in pdata.columns: y_path = np.array(pdata.Y)
    if 'y' in pdata.columns: y_path = np.array(pdata.y)
    if 'Z' in pdata.columns: z_path = np.array(pdata.Z)
    if 'z' in pdata.columns: z_path = np.array(pdata.z)

#    output which file is analyzed
    if statusoutput:
        print('[PD] Creating new raster along polyline: %s' % polyfnames)

#     erzeugung der eckpuntke der segmente des unregelmaesigen gitters,
#     Domain Boundaries DB
#     input: mittlerer path
#     output: eckpunkte für punkt entlang der linie
    DB_x_rl, DB_y_rl, DB_x_csz, DB_y_csz = path2domain(x_path, y_path,
                                                       w/2., cellsize)

    # plt.figure(2)
    # plt.plot(DB_x_rl, DB_y_rl,'r')
    # plt.plot(DB_x_csz, DB_y_csz,'b')
    # plt.plot(x_path, y_path,'k')
    # plt.show()
#     Shift path to raster because raster is not in global grid
    DB_x_rl -= xllcenter
    DB_y_rl -= yllcenter
#    for calculation cellsize
    DB_x_csz -= xllcenter
    DB_y_csz -= yllcenter

    # plt.figure(3)
    # plt.plot(DB_x_rl, DB_y_rl,'r')
    # plt.plot(DB_x_csz, DB_y_csz,'b')
    # plt.plot(x_path-xllcenter, y_path-yllcenter,'k')
    # plt.show()
#    use bresemham algorithm to determine the new raster
    for i in range(len(DB_x_rl[0])):
#   for each segment check the number of CROSS-cells
        n = bresenham(DB_x_rl[0, i], DB_y_rl[0, i],
                      DB_x_rl[1, i], DB_y_rl[1, i], cellsize)
        n_total = max(n_total, len(n))
#   number of raster cells of edges parallel to polyline
    m = np.zeros(len(DB_x_rl[0])-1).astype('int')
    for i in range(len(DB_x_rl[0])-1):
        # left edge
        zl = bresenham(DB_x_rl[0, i], DB_y_rl[0, i],
                       DB_x_rl[0, i+1], DB_y_rl[0, i+1], cellsize)
        # right edge
        zr = bresenham(DB_x_rl[1, i], DB_y_rl[1, i],
                       DB_x_rl[1, i+1], DB_y_rl[1, i+1], cellsize)
        m[i] = max(len(zl), len(zr))
#    delete the lines that are double at ech segment connection
    m_total = sum(m) - (len(DB_x_rl[0])-2)

#    Calculation of segments
    new_raster = np.zeros((2, m_total, n_total)) + np.NaN
    # new raster filled with NaN

#    Each dataset needs its own s_coord, l_coord. This is saved to a cell array
#    and returned from this function for use in other parts of the
#    program package

#    l_coord is distance from polyline
    l_coord = np.linspace(-w/2, w/2, n_total)
    s_coord = np.zeros(m_total)
    ds1 = 0
    if statusoutput: # output which file is analyzed
        print('[PD] Transferring data from old to new raster')
    for i in range(len(DB_x_rl[0])-1): # for each segment
        # Division of side edges in n segments
        # x-/y-values of lines are linspace(x0,x1,m) etc.
        # DB_x_rl, DB_y_rl are values from polyline2path
        bxl = np.linspace(DB_x_rl[0][i], DB_x_rl[0][i+1], m[i]) # left
        byl = np.linspace(DB_y_rl[0][i], DB_y_rl[0][i+1], m[i])

        bxr = np.linspace(DB_x_rl[1][i], DB_x_rl[1][i+1], m[i]) # right
        byr = np.linspace(DB_y_rl[1][i], DB_y_rl[1][i+1], m[i])

        # => generation of grid points for each segment
        # grid points can be assigned to original raster data, using
        # rasterize function
        new_rastersegment = np.zeros((2, m[i], n_total))

        for j in range(m[i]):
            x = np.linspace(bxl[j], bxr[j], n_total) # line coordinates x
            y = np.linspace(byl[j], byr[j], n_total) # line coordinates y
            for k in range(n_total):
                # x,y-Koordinaten of cells on line
                # xy_coord = bresenham(x(k),y(k),x(k),y(k),cellsize);
                xy_coord = [round(x[k]/cellsize) * cellsize,
                            round(y[k]/cellsize) * cellsize]
                # cell coordinates of new raster
                xy_ind = [xy_coord[0]/cellsize +1, xy_coord[1]/cellsize +1]
                # translate coordinate of cell to cell index
                # THIS IS THE NEAREST NEIGHBOUR APPROXIMATION
                # Assign pressure data to loc
                new_rastersegment[:, j, k] = [xy_ind[0], xy_ind[1]]

#        For each segment following the first we must delete the first
#        line since it is identical to the last line of the previous
#        segment.

#        s_coord = x-Coordinate along Polyline.
#        % Start of Polylinie at x = 0
        m_neu = m[i]

        if (i == 0):
            # Distance from starting point of current segment
            # from beginning of polyline
            ds0 = 0
        else:
            ds0 += ds1

        ds1 = math.sqrt((x_path[i+1]-x_path[i])**2 +
                        (y_path[i+1]-y_path[i])**2)
        new_raster[:, m_alt:m_alt+m_neu, :] = [new_rastersegment[0],
                                               new_rastersegment[1]]

        s_coord[m_alt:m_alt+m_neu] = np.linspace(ds0, ds0+ds1, m[i])
        m_alt = m_alt+m_neu-1
        if (i==0):
            s_coordmin = s_coord[0]
        s_coord -= s_coordmin

#    calclation of cellsize (for area)
    new_raster_area = np.zeros((m_total, n_total)) + np.NaN
    sum_mi = 0

#    if m_total % 2 == 0:
#        seg_boundary_lines = len(DB_x_csz)-2
#    else:
    seg_boundary_lines = len(DB_x_csz)-1
    for i in range(seg_boundary_lines): # for offset segment boundary lines
        # DB_x_rl DB_y_rl for i --> Koordinaten
        x_DB_i = np.linspace(DB_x_csz[i][0], DB_x_csz[i][1], n_total+1)
        y_DB_i = np.linspace(DB_y_csz[i][0], DB_y_csz[i][1], n_total+1)
        # DB_x_rl DB_y_rl for i+1 --> Koordinaten
        x_DB_ii = np.linspace(DB_x_csz[i+1][0], DB_x_csz[i+1][1], n_total+1)
        y_DB_ii = np.linspace(DB_y_csz[i+1][0], DB_y_csz[i+1][1], n_total+1)

        if i % 2 == 0: # i gerade
            for j in range(n_total):
                x_seg_j = [x_DB_i[j], x_DB_ii[j]]
                y_seg_j = [y_DB_i[j], y_DB_ii[j]]

                x_seg_jj = [x_DB_i[j+1], x_DB_ii[j+1]]
                y_seg_jj = [y_DB_i[j+1], y_DB_ii[j+1]]

                k = 0
                a = sum_mi
                new_raster_area[a][j] = 1./2 * ((y_seg_j[k]-y_seg_jj[k+1]) *
                                                (x_seg_jj[k]-x_seg_j[k+1]) +
                                                (y_seg_j[k+1]-y_seg_jj[k]) *
                                                (x_seg_j[k]-x_seg_jj[k+1]))

            sum_mi = sum_mi+1
        else: # i ungerade
            m_i = round(i/2) # m for each segment
            for j in range(n_total):
                x_seg_j = np.linspace(x_DB_i[j], x_DB_ii[j], m[m_i]-1)
                y_seg_j = np.linspace(y_DB_i[j], y_DB_ii[j], m[m_i]-1)

                x_seg_jj = np.linspace(x_DB_i[j+1], x_DB_ii[j+1], m[m_i]-1)
                y_seg_jj = np.linspace(y_DB_i[j+1], y_DB_ii[j+1], m[m_i]-1)

                for k in range(m[m_i]-2):
                    a = sum_mi+k
                    # print(np.shape(new_raster_area))
                    # print(a,j)
                    new_raster_area[a, j] = 1./2*((y_seg_j[k]-y_seg_jj[k+1]) *
                                                  (x_seg_jj[k]-x_seg_j[k+1]) +
                                                  (y_seg_j[k+1]-y_seg_jj[k]) *
                                                  (x_seg_j[k]-x_seg_jj[k+1]))

            sum_mi = a + 1

    print('[PD] Size of rasterdata- old: %d x %d - new: %d x %d' % (
    np.size(rasterdata, 0), np.size(rasterdata,1),
    np.size(new_raster, 1), np.size(new_raster, 2)))

    aval_data[0] = header
    aval_data[1] = s_coord
    aval_data[2] = l_coord
    aval_data[3] = new_raster
    aval_data[4] = abs(new_raster_area)

    # visu
    figure_width = 2*5
    figure_height = 2*4
    lw = 1

    fig = plt.figure(figsize=(figure_width, figure_height), dpi=150)

#    for figure: referenz-simulation bei p_lim=1
    new_rasterdata = rasterdata
    masked_array = np.ma.masked_where(new_rasterdata==0, new_rasterdata)
    cmap = copy.copy(matplotlib.cm.jet)
    cmap.set_bad('w', 1.)

    n, m = np.shape(new_rasterdata)
    xx, yy = np.meshgrid(np.arange(m), np.arange(n))

    ref1 = plt.imshow(masked_array, vmin=new_rasterdata.min(),
               vmax=new_rasterdata.max(),
               origin='lower',
               cmap=cmap,
               label='pressure data',
               aspect='auto',
               extent=[xx.min()*cellsize+xllcenter, xx.max()*cellsize+xllcenter,
                       yy.min()*cellsize+yllcenter, yy.max()*cellsize+yllcenter])
    plt.autoscale(False)
    ref2 = plt.plot(x_path, y_path,
                    'b-', linewidth = lw, label = 'flow path')
    ref3 = plt.plot(DB_x_rl+xllcenter, DB_y_rl+yllcenter,
                    'g-', linewidth = lw, label = 'domain')
    ref3 = plt.plot(DB_x_rl.T+xllcenter, DB_y_rl.T+yllcenter,
                    'g-', linewidth = lw, label = 'domain')
    refs = [ref2[0], ref3[0]]

    labels = ['flow path', 'domain']
         # falls doku vorhanden
    if with_doku and isinstance(doku_name, basestring):
        depdata = np.loadtxt(doku_name)
        xdep = depdata[:, 0]
        ydep = depdata[:, 1]
        ref4 = plt.plot(xdep, ydep,
                        'r', linewidth = lw, label = 'doku data')
        refs.append(ref4[0])
        labels.append('doku data')

    if with_damages and isinstance(damages_name[0], basestring):
        damage_name = damages_name
#        for damage_name in damages_name:
        damdata = np.loadtxt(damage_name)
        if len(damdata.shape) == 1:
            xdam = damdata[0]
            ydam = damdata[1]
            ref5 = plt.plot(xdam, ydam, 'm*', linewidth = lw, label = 'damages')
        else:
            xdam = damdata[:, 0]
            ydam = damdata[:, 1]
            ref5 = plt.plot(xdam, ydam, 'm-*', linewidth = lw, label = 'damages')
        refs.append(ref5[0])
        labels.append('damages')

    plt.legend(refs, labels, loc=0)
    plt.xlim([xx.min()*cellsize+xllcenter, xx.max()*cellsize+xllcenter])
    plt.ylim([yy.min()*cellsize+yllcenter, yy.max()*cellsize+yllcenter])
    plt.xlabel('x [m]')
    plt.ylabel('y [m]')
    cbh = plt.colorbar()
    cbh.set_label('peak pressure [kPa]')
    # plt.show()
    if out:
        pro_name = fnames[0].split('/')[-3] # CoSiCa-samos-structure
#        pro_name = fnames[0].split('/')[-5] + '_' + fnames[0].split('/')[-2] # DAKUMO_structure
        outname_fin = ''.join([outpath, '/pics/', pro_name,
                               '_dptr', str(int(dpp_threshold)),
                               '_simulationxy','.pdf'])
        if not os.path.exists(os.path.dirname(outname_fin)):
            os.makedirs(os.path.dirname(outname_fin))
        fig.savefig(outname_fin, transparent=True)

    plt.close(fig)

    return aval_data


#==============================================================================
# rasterize
#==============================================================================
def rasterize(xcoords, ycoords, width, height):
    """
    Simple function to rasterize a polygon
    using the Qt toolkit (QPainter on QImage)
    Takes a list the x-coordiantes, a list of the y-coordinates,
    the width (integer) and the height (integer)

    Depends on PyQt4.QtCore and PyQt4.QtGui
    """

    polygon = QtGui.QPolygonF()

    for x,y in zip(xcoords, ycoords):
        polygon.append(QtCore.QPointF(x,y))

    return rasterizePolygon(polygon, QtCore.QSize(width, height))

#==============================================================================
# rasterizePolygon
#==============================================================================
def rasterizePolygon(polygon, size):
    """
    Simple function to rasterize a polygon
    Takes a QPolygonF and QSize as input
    """

    #use QPainter to print the polygon
    painterpath = QtGui.QPainterPath()
    painterpath.addPolygon(polygon)
    img = QtGui.QImage(size, QtGui.QImage.Format_ARGB32)
    img.fill(QtGui.QColor(0, 0, 0, 255))
    p = QtGui.QPainter(img)
    p.fillPath(painterpath, QtGui.QColor(1, 1, 1, 255))
    p.end()

    #extract pointer to QImage data structure
    ptr = img.bits()
    ptr.setsize(img.byteCount())

    #reshape and copy the  QImage data structures
    return np.copy(np.asarray(ptr).reshape(img.height(), img.width(), 4)[:,:,0])



def poly2mask_simple(ydep,xdep,ncols,nrows):
    mask = np.zeros((nrows,ncols))
    xyframe = bresenham(xdep[0],ydep[0],xdep[1],ydep[1],1)
    xyframe = np.delete(xyframe,-1,0)
    xyframe = np.transpose(xyframe)
    for i in range(1,len(xdep)-1):
        xyline = bresenham(xdep[i],ydep[i],xdep[i+1],ydep[i+1],1)
        xyline = np.delete(xyline,-1,0) # letzer punkt ist erster punkt der neun linie
        xyline = np.transpose(xyline)
        xyframe = np.hstack((xyframe,xyline));

    xyline = bresenham(xdep[-1],ydep[-1],xdep[0],ydep[0],1)
    xyline = np.delete(xyline,-1,0)
    xyline = np.transpose(xyline)
    xyframe = np.hstack((xyframe,xyline));
    for i in range(0,len(xyframe[0,:])) :
        mask[xyframe[0,i],xyframe[1,i]] = 1

    # inneres des polygons mit einsen füllen
    # [i, j] = find(mask);
    # i, j = np.where(mask == 1)
    i = xyframe[0]
    j = xyframe[1]
    mv, nv = np.meshgrid(np.linspace(0,nrows-1,nrows),np.linspace(0,ncols-1,ncols)) # create index space
    mask = inpolygon(mv,nv,i,j)
    mask = np.transpose(mask)
    return mask

def inpolygon(X, Y, xv, yv):

    npol = len(xv)
    lx = np.shape(X)[0]
    ly = np.shape(Y)[1]
    IN = np.zeros(np.shape(X))
    j = npol-1
    for i in range(npol-1):
        delta_xv = xv[j] - xv[i]
        delta_yv = yv[j] - yv[i]
        ## distance = [distance from (X,Y) to edge] * length(edge)
        distance = delta_xv*(Y-yv[i]) - (X-xv[i])*delta_yv
        ## is Y between the y-values of edge i,j
        ##        AND (X,Y) on the left of the edge ?
        for ii in range(lx):
            for jj in range(ly):
                if (((yv[i]<=Y[ii][jj] and Y[ii][jj]<yv[j]) or (yv[j]<=Y[ii][jj] and Y[ii][jj]<yv[i]) ) and 0< distance[ii][jj]*delta_yv):
                    if IN[ii][jj]==0:
                        IN[ii][jj] = 1
                    else:
                        IN[ii][jj] = 0
        j = i
    for i in range(npol-1):
        IN[yv[i]][xv[i]] = 1

    return IN


#==============================================================================
# create_mask
#==============================================================================
def createMask(rasterIndData, polydepfnames, p_lim):
    """
    creating rasterdata for deposit

    this function is used to create rasterdata for the deposit-polyline (or deposit-map)
    and transfer the rasterdata from samos-raster --> deskewed_raster
    for grid point in --> 1
    for grid point outside the polyline --> 0

    Andreas Kofler, 2013

    input: names of rasterfiles, name of depositfile, new_raster_ind,
        new_raster_pressure, new_raster_depth
    ouput: structure{x coordinate along new raster,
                     y coordinate, rasterdata_xind, rasterdata_yind}
    """

    header = rasterIndData[0]
#    create mask from deposit-polygon
    ncols = int(header.ncols)
    nrows = int(header.nrows)
    xllcenter = float(header.xllcenter)
    yllcenter = float(header.yllcenter)
    cellsize = float(header.cellsize)

    if (QtCore.QString(polydepfnames).endsWith('.asc') or QtCore.QString(polydepfnames).endsWith('.txt')):
        print('[CM] use mask from rasterfile: %s' % polydepfnames.split('/')[-1])
        mask = np.flipud(np.loadtxt(str(polydepfnames), skiprows=6))
        mask[np.where(mask < p_lim)] = 0
        mask[np.where(mask >= p_lim)] = 1
    elif QtCore.QString(polydepfnames).endsWith('.xyz'):
        print('[CM] create mask from polyline: %s' % polydepfnames.split('/')[-1])
        depdata =  np.loadtxt(polydepfnames)
        if len(depdata.shape) == 1:
            mask = np.zeros((ncols, nrows))
            mask[int((depdata[0]-xllcenter)/cellsize), int((depdata[1]-yllcenter)/cellsize)] = 1
            mask = mask.T
        elif len(depdata) == 2:
            mask = np.zeros((ncols, nrows))
            xdep = (depdata[:, 0]-xllcenter)/cellsize
            ydep = (depdata[:, 1]-yllcenter)/cellsize
            nn = bresenham(xdep[0], ydep[0], xdep[1], ydep[1], cellsize)
            xx = np.array(np.linspace(xdep[0], xdep[1], len(nn)), dtype=int)
            yy = np.array(np.linspace(ydep[0], ydep[1], len(nn)), dtype=int)
            mask[xx, yy] = 1
            mask = mask.T
        else:
            xdep = (depdata[:, 0]-xllcenter)/cellsize
            ydep = (depdata[:, 1]-yllcenter)/cellsize
            #Teste ob Polygon geschlossen und im Fall schließen
            if xdep[0] != xdep[-1] or ydep[0] != ydep[-1]:
                print ('[CM]: Polygonzug nicht geschlossen, 1. Punkt wird hinten angehängt')
                xdep = np.append(xdep, xdep[0])
                ydep = np.append(ydep, ydep[0])
            if len(xdep) != len(ydep):
                warnings.warn('[CM] Y und X Werte des Polygons sind nicht gleich lang')
            # mask = rasterize(xdep, ydep, ncols, nrows)
            mask = poly2mask_simple(xdep, ydep, ncols, nrows)
            print ('[CM]: we are in the else')
    else:
        warnings.warn('file ending not known. options are txt, asc, xyz.')
        return

#   transformation mask in new raster (like in assign_data)
    new_mask = np.zeros((len(rasterIndData[1]), len(rasterIndData[2])))
    xy_oldind = rasterIndData[3].astype('int')

    print('[CM] Transferring mask from old to new raster')

    i_oob = 0  # counter out of bounds
    i_ib = 0   # counter in bounds
    for x_ind in range(new_mask.shape[0]):
        for y_ind in range(new_mask.shape[1]):
            i_ib += 1
            try:
                new_mask[x_ind, y_ind] = mask[xy_oldind[1][x_ind, y_ind]][xy_oldind[0][x_ind, y_ind]]
            except:
                i_oob += 1
                new_mask[x_ind][y_ind] = np.NaN
    print('[AD] %d raster values transferred and %d are out of original raster bounds!' % (i_ib-i_oob, i_oob))

#    aval_data[0] = header
#    aval_data[1] = rasterIndData[1]
#    aval_data[2] = rasterIndData[2]
#    aval_data[3] = new_mask
    aval_data = new_mask

    return aval_data

#==============================================================================
# transform
#==============================================================================
def transform(fname, rasterIndData, statusoutput=0):
    name = fname.split('/')

    xy_oldind = rasterIndData[3].astype('int')

    header = IOf.readASCheader(fname)
    xllcenter = header.xllcenter
    yllcenter = header.yllcenter
    cellsize = header.cellsize

    rasterdata = np.flipud(np.loadtxt(fname, skiprows=6))
    rasterdata[rasterdata == header.noDataValue] = np.NaN

#        out of bounds counter
    i_oob = 0
    i_ib  = 0

    new_raster = np.zeros((len(rasterIndData[1]), len(rasterIndData[2])))

    for x_ind in range(new_raster.shape[0]):
        for y_ind in range(new_raster.shape[1]):
            i_ib += 1
            try:
                new_raster[x_ind, y_ind] = rasterdata[xy_oldind[1][x_ind, y_ind]][xy_oldind[0][x_ind, y_ind]]
            except:
                i_oob += 1
                new_raster[x_ind, y_ind] = np.NaN
    if statusoutput:
        print('[AD] Data-file: %s - %d raster values transferred - %d out of original raster bounds!' % (name[-1], i_ib-i_oob, i_oob))

#        aval_data[topo_num][0] = header
#        aval_data[topo_num][1] = rasterIndData[1]
#        aval_data[topo_num][2] = rasterIndData[2]
#        aval_data[topo_num][3] = new_raster
    return new_raster

#==============================================================================
# assign_data
#==============================================================================
def assignData(fnames, rasterIndData):
    """
    process data

    this function is used to process the rasterdata such that it can be
    analysed with the methods for a regular grid
    data given in a regulare grid is projected on a nonuniform grid given by
    a polyline

    JT Fischer, Uwe Schlifkowitz BFW 2010-2012
    AK BFW 2014

    input: names of rasterfiles, poly names, path width
    ouput: structure{x coordinate along new raster, y coordinate, rasterdata_xind, rasterdata_yind}

    fnames:       full name of data file including path and extension
    polyfnames:   file with polylinie
    rasterIndData: raster of new shape containing indices from corresponding points of old raster
    """
# if started without arguments
#    if (nargin == 0):
#        return None

    maxtopo = len(fnames)
#    aval_data = np.array([[None for m in xrange(4)] for n in xrange(maxtopo)])
    aval_data = np.array(([None] * maxtopo))

    print( '[AD] Transfer data of %d file(s) from old to new raster' % maxtopo)

    pool = Pool()
    aval_data = pool.map(functools.partial(transform, rasterIndData = rasterIndData, statusoutput = 1), fnames)
    pool.close()
    pool.join()

    return aval_data

#==============================================================================
# rasterize example
#==============================================================================
if __name__ == '__main__':
    import matplotlib.pyplot as plt
    a = [1.1, 5.3, 8, 8.5, 2]
    b = [1.1, 5, 2, 8.8, 8.0]
    arr = rasterize(a,b, 10, 15);
    arr1 = poly2mask_simple(a,b, 10, 15);
    pressurefileList = ['/home/matthiastonnel/Documents/gitea/AvaFrame_PSAM_UH/NameOfAvalanche/Outputs/dfa_pressure/000001.txt']
    depthfileList = ['/home/matthiastonnel/Documents/gitea/AvaFrame_PSAM_UH/NameOfAvalanche/Outputs/dfa_depth/000001.txt']
    ava_data = processDataInd(pressurefileList, '/home/matthiastonnel/Documents/gitea/AvaFrame_PSAM_UH/NameOfAvalanche/Inputs/avalanche_path.xyz', 300, False, 'doku_name', False, 'damages_name',
         'visu', './', False, 'dpp_threshold')
    deskewedRasterPressure = assignData(pressurefileList, ava_data)
    deskewedRasterDepth = assignData(depthfileList, ava_data)

    fig, ax = plt.subplots()
    ax.imshow(deskewedRasterPressure[0], cmap=plt.cm.jet, interpolation='nearest')
    plt.show()


    # a = np.append(a,a[0])
    # b = np.append(b,b[0])
    # fig, ax = plt.subplots()
    # ax.imshow(arr, cmap=plt.cm.jet, interpolation='nearest')
    # plt.plot(a,b,'k')
    # # ax.imshow(deskewedRasterPressure[0], cmap=plt.cm.jet, interpolation='nearest')
    # fig, ax = plt.subplots()
    # ax.imshow(arr1, cmap=plt.cm.jet, interpolation='nearest')
    # plt.plot(a,b,'k')
    # # ax.imshow(deskewedRasterPressure[0], cmap=plt.cm.jet, interpolation='nearest')
    # plt.show()
