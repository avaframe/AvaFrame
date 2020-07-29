import os
import logging
# create local logger
# change log level in calling module to DEBUG to see log messages
log = logging.getLogger(__name__)
import pickle
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import datetime

colors = ["#393955", "#8A8A9B", "#E9E940"]
mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=colors)


<<<<<<< HEAD
def readABresults(saveOutPath,name):
    savename = name + '_com2AB_param.pickle'
    save_file = os.path.join(saveOutPath, savename)
    with open(save_file, 'rb') as handle:
        eqParams = pickle.load(handle)
    savename = name + 'com2AB_out.pickle'
    save_file = os.path.join(saveOutPath, savename)
=======
def readABresults(saveOutPath):
    # TODO: the rest is plotting + analytical stuff
    # I (FSO) suggest we just write it out and move these functions to IO module?
    save_file = os.path.join(saveOutPath, 'com2AB_param.pickle')
    with open(save_file, 'rb') as handle:
        eqParams = pickle.load(handle)
    save_file = os.path.join(saveOutPath, 'com2AB_out.pickle')
>>>>>>> 6175d0dc6ad423c90515b6bb4a0afc37978790f6
    with open(save_file, 'rb') as handle:
        eqOut = pickle.load(handle)

    return eqParams, eqOut

<<<<<<< HEAD
def processABresults(eqParams, eqOut):
=======
def processABresults(eqParams, eqOut, saveOutPath):
>>>>>>> 6175d0dc6ad423c90515b6bb4a0afc37978790f6

    s = eqOut['s']
    x = eqOut['x']
    y = eqOut['y']
    z = eqOut['z']
    CuSplit = eqOut['CuSplit']
    ids_10Point = eqOut['ids_10Point']
    poly = eqOut['poly']
    beta = eqOut['beta']
    alpha = eqOut['alpha']
    SDs = eqOut['SDs']
    alphaSD = eqOut['alphaSD']
    indSplit = eqOut['indSplit']
    # Line down to alpha
    f = z[0] + np.tan(np.deg2rad(-alpha)) * s
    fplus1SD = z[0] + np.tan(np.deg2rad(-alphaSD[0])) * s
    fminus1SD = z[0] + np.tan(np.deg2rad(-alphaSD[1])) * s
    fminus2SD = z[0] + np.tan(np.deg2rad(-alphaSD[2])) * s

    # First it calculates f - g and the corresponding signs
    # using np.sign. Applying np.diff reveals all
    # the positions, where the sign changes (e.g. the lines cross).
    ids_alpha = np.argwhere(np.diff(np.sign(f - z))).flatten()
    ids_alphaP1SD = np.argwhere(np.diff(np.sign(fplus1SD - z))).flatten()
    ids_alphaM1SD = np.argwhere(np.diff(np.sign(fminus1SD - z))).flatten()
    ids_alphaM2SD = np.argwhere(np.diff(np.sign(fminus2SD - z))).flatten()

<<<<<<< HEAD

    # Only get the first index past the splitpoint
    try:
        ids_alpha = ids_alpha[s[ids_alpha] > CuSplit][0]
    except:
        log.info('Alpha out of profile')
=======
    # Only get the first index past the splitpoint
    try:
        ids_alpha = ids_alpha[s[ids_alpha] > CuSplit][0]
    except log.info('Alpha out of profile'):
>>>>>>> 6175d0dc6ad423c90515b6bb4a0afc37978790f6
        ids_alpha = None

    try:
        ids_alphaP1SD = ids_alphaP1SD[s[ids_alphaP1SD] > CuSplit][0]
<<<<<<< HEAD
    except:
        log.info('+1 SD above beta point')
=======
    except log.info('+1 SD above beta point'):
>>>>>>> 6175d0dc6ad423c90515b6bb4a0afc37978790f6
        ids_alphaP1SD = None

    try:
        ids_alphaM1SD = ids_alphaM1SD[s[ids_alphaM1SD] > CuSplit][0]
<<<<<<< HEAD
    except:
        log.info('-1 SD out of profile')
=======
    except log.info('-1 SD out of profile'):
>>>>>>> 6175d0dc6ad423c90515b6bb4a0afc37978790f6
        ids_alphaM1SD = None

    try:
        ids_alphaM2SD = ids_alphaM2SD[s[ids_alphaM2SD] > CuSplit][0]
<<<<<<< HEAD
    except:
        log.info('-2 SD out of profile')
=======
    except log.info('-2 SD out of profile'):
>>>>>>> 6175d0dc6ad423c90515b6bb4a0afc37978790f6
        ids_alphaM2SD = None


    eqOut['f'] = f
    eqOut['ids_alpha'] = ids_alpha
    eqOut['ids_alphaP1SD'] = ids_alphaP1SD
    eqOut['ids_alphaM1SD'] = ids_alphaM1SD
    eqOut['ids_alphaM2SD'] = ids_alphaM2SD

    return eqOut

<<<<<<< HEAD

def writeABpostOut(header, rasterdata, Avapath, SplitPoint, saveOutPath):
    """ Loops on the given Avapath and runs AlpahBeta Postprocessing
    """
    NameAva = Avapath['Name']
    StartAva = Avapath['Start']
    LengthAva = Avapath['Length']
    CoordAva = Avapath['Coord']

    CoordSplit = SplitPoint['Coord']

    for i in range(len(NameAva)):
        name = NameAva[i]
        OutPath = saveOutPath + 'Outputs/'
        start = StartAva[i]
        end = start + LengthAva[i] - 1
        avapath = CoordAva[:,int(start):int(end)]
        eqParams, eqOut = readABresults(OutPath,name)
        eqPost = processABresults(eqParams, eqOut)
        ParameterSet = eqParams['ParameterSet']
        # Plot the whole profile with beta, alpha ... points and lines
        plotSaveResults(header, rasterdata, avapath, CoordSplit, eqPost, ParameterSet, saveOutPath)

        WriteResults(eqPost, ParameterSet, saveOutPath)
=======
def writeABpostOut(header, rasterdata, avapath, splitPoint, saveOutPath):

    eqParams, eqOut = readABresults(saveOutPath)
    eqPost = processABresults(eqParams, eqOut, saveOutPath)
    ParameterSet = eqParams['ParameterSet']
    # Plot the whole profile with beta, alpha ... points and lines
    plotSaveResults(header, rasterdata, avapath, splitPoint, eqPost, ParameterSet, saveOutPath)

    WriteResults(eqPost, ParameterSet, saveOutPath)
>>>>>>> 6175d0dc6ad423c90515b6bb4a0afc37978790f6


def plotSaveResults(header, rasterdata,avapath, splitPoint, eqOutput, ParameterSet, saveOutPath):
    """TODO: - doc
    - split into save and plot"""

    s = eqOutput['s']
    x = eqOutput['x']
    y = eqOutput['y']
    z = eqOutput['z']
    CuSplit = eqOutput['CuSplit']
    ids_10Point = eqOutput['ids_10Point']
    poly = eqOutput['poly']
    beta = eqOutput['beta']
    alpha = eqOutput['alpha']
    SDs = eqOutput['SDs']
    alphaSD = eqOutput['alphaSD']
    f = eqOutput['f']
    ids_alpha = eqOutput['ids_alpha']
    ids_alphaP1SD = eqOutput['ids_alphaP1SD']
    ids_alphaM1SD = eqOutput['ids_alphaM1SD']
    ids_alphaM2SD = eqOutput['ids_alphaM2SD']
    indSplit = eqOutput['indSplit']

<<<<<<< HEAD

=======
>>>>>>> 6175d0dc6ad423c90515b6bb4a0afc37978790f6
    # Plot raster and path
    fig1, ax1 = plt.subplots()
    cmap = mpl.cm.Greys
    cmap.set_bad(color='white')
    im1 = plt.imshow(rasterdata, cmap, origin='lower')
    divider = make_axes_locatable(ax1)
    cax = divider.append_axes("right", size="5%", pad=0.1)
    cb1 = fig1.colorbar(im1, cax=cax)
    path1 = ax1.plot((x-header.xllcorner)/header.cellsize,
                     (y-header.yllcorner)/header.cellsize)
    ax1.plot((x-header.xllcorner)/header.cellsize,
             (y-header.yllcorner)/header.cellsize, 'k')
    ax1.plot((x[indSplit]-header.xllcorner)/header.cellsize,
             (y[indSplit]-header.yllcorner)/header.cellsize, '.',
             color='0.6', label='Projection header, rasterdata, avapath, splitPoint, of Split Point on ava path')
    ax1.plot((splitPoint[0]-header.xllcorner)/header.cellsize,
             (splitPoint[1]-header.yllcorner)/header.cellsize, '.',
             color='0.3', label='Split point''Projection of Split Point on ava path')
    plt.show()



    # Plot the whole profile with beta, alpha ... points and lines
    plt.close("all")
    fig = plt.figure(4, figsize=(10, 6))
    titleText = 'Profile'
    plt.title(titleText)

    xlabelText = 'Distance [m]\nBeta: '+str(round(beta, 1)) + '$^\circ$' + \
        '  Alpha: '+str(round(alpha, 1)) + '$^\circ$'
    plt.xlabel(xlabelText, multialignment='center')

    plt.ylabel('Height [m]')

    plt.plot(s, z, '-', label='DEM')
    plt.plot(s, poly(s), ':', label='QuadFit')
    plt.axvline(x=s[indSplit], color='0.7',
                linewidth=1, linestyle='--', label='Split point')
    plt.axvline(x=s[ids_10Point], color='0.8',
                linewidth=1, linestyle='-.', label='Beta')

    plt.plot(s, f, '-', label='AlphaLine')
<<<<<<< HEAD
    if ids_alphaM1SD:
=======
    if ids_alpha:
>>>>>>> 6175d0dc6ad423c90515b6bb4a0afc37978790f6
        plt.plot(s[ids_alphaM1SD], z[ids_alphaM1SD], 'x', markersize=8,
                 label='Alpha - 1SD')
    if ids_alphaM2SD:
        plt.plot(s[ids_alphaM2SD], z[ids_alphaM2SD], 'x', markersize=8,
                 label='Alpha - 2SD')

    ax = plt.gca()

    versionText = datetime.datetime.now().strftime("%d.%m.%y") + \
        '; ' + 'AlphaBeta ' + ParameterSet
    plt.text(0, 0, versionText, fontsize=8, verticalalignment='bottom',
             horizontalalignment='left', transform=ax.transAxes,
             color='0.5')
    # plt.text(-0.2, 0, 'matplotlib -2', \
    #          verticalalignment='center', transform=ax.transAxes)

    plt.gca().set_aspect('equal', adjustable='box')
    plt.grid(linestyle=':', color='0.9')
    plt.legend(frameon=False)
    plt.draw()
    plt.show()
    log.info('Saving profile figure to: %s', saveOutPath)
    save_file = os.path.join(saveOutPath, 'AlphaBeta_.pdf')
    plt.savefig(save_file)
    plt.close(fig)
    plt.close("all")


def WriteResults(eqOutput, ParameterSet, saveOutPath):
    """ Write com2AB results to file """
    s = eqOutput['s']
    x = eqOutput['x']
    y = eqOutput['y']
    z = eqOutput['z']
    CuSplit = eqOutput['CuSplit']
    ids_10Point = eqOutput['ids_10Point']
    poly = eqOutput['poly']
    beta = eqOutput['beta']
    alpha = eqOutput['alpha']
    SDs = eqOutput['SDs']
    alphaSD = eqOutput['alphaSD']
    f = eqOutput['f']
    ids_alpha = eqOutput['ids_alpha']
    ids_alphaP1SD = eqOutput['ids_alphaP1SD']
    ids_alphaM1SD = eqOutput['ids_alphaM1SD']
    ids_alphaM2SD = eqOutput['ids_alphaM2SD']
    indSplit = eqOutput['indSplit']


    log.info('Parameter Set %s\n' % ParameterSet)
    log.info('Alpha point (x,y,z,s) in [m]:(%.2f,%.2f,%.2f,%.2f) and angle in [°] : %.2f\n' % (
        x[ids_alpha], y[ids_alpha], z[ids_alpha], s[ids_alpha], alpha))
    log.info('Beta point (x,y,z,s) in [m]:(%.2f,%.2f,%.2f,%.2f) and angle in [°] : %.2f\n' % (
        x[ids_10Point], y[ids_10Point], z[ids_10Point], s[ids_10Point], beta))
<<<<<<< HEAD
    if ids_alphaM1SD:
        log.info('alphaM1SD point (x,y,z,s) in [m]:(%.2f,%.2f,%.2f,%.2f) and angle in [°] : %.2f\n' % (
            x[ids_alphaM1SD], y[ids_alphaM1SD], z[ids_alphaM1SD], s[ids_alphaM1SD], alphaSD[1]))
    else:
        log.info('alphaM1SD point out of profile\n')
    if ids_alphaM2SD:
        log.info('alphaM2SD point (x,y,z,s) in [m]:(%.2f,%.2f,%.2f,%.2f) and angle in [°] : %.2f\n' % (
            x[ids_alphaM2SD], y[ids_alphaM2SD], z[ids_alphaM2SD], s[ids_alphaM2SD], alphaSD[2]))
    else:
        log.info('alphaM2SD point out of profile\n')
=======
    log.info('alphaM1SD point (x,y,z,s) in [m]:(%.2f,%.2f,%.2f,%.2f) and angle in [°] : %.2f\n' % (
        x[ids_alphaM1SD], y[ids_alphaM1SD], z[ids_alphaM1SD], s[ids_alphaM1SD], alphaSD[1]))
    log.info('alphaM2SD point (x,y,z,s) in [m]:(%.2f,%.2f,%.2f,%.2f) and angle in [°] : %.2f\n' % (
        x[ids_alphaM2SD], y[ids_alphaM2SD], z[ids_alphaM2SD], s[ids_alphaM2SD], alphaSD[2]))
>>>>>>> 6175d0dc6ad423c90515b6bb4a0afc37978790f6
    log.info('alphaP1SD point (x,y,z,s) in [m]:(%.2f,%.2f,%.2f,%.2f) and angle  in [°] : %.2f\n' % (
        x[ids_alphaP1SD], y[ids_alphaP1SD], z[ids_alphaP1SD], s[ids_alphaP1SD], alphaSD[0]))

    FileName_ext = saveOutPath + 'results_python.txt'
    with open(FileName_ext, 'w') as outfile:
        outfile.write('Parameter Set %s\n' % ParameterSet)
        outfile.write('Alpha point (x,y,z,s) in [m]:(%.2f,%.2f,%.2f,%.2f) and angle in [°] : %.2f\n' % (
            x[ids_alpha], y[ids_alpha], z[ids_alpha], s[ids_alpha], alpha))
        outfile.write('Beta point (x,y,z,s) in [m]:(%.2f,%.2f,%.2f,%.2f) and angle in [°] : %.2f\n' % (
            x[ids_10Point], y[ids_10Point], z[ids_10Point], s[ids_10Point], beta))
<<<<<<< HEAD
        if ids_alphaM1SD:
            outfile.write('alphaM1SD point (x,y,z,s) in [m]:(%.2f,%.2f,%.2f,%.2f) and angle in [°] : %.2f\n' % (
                x[ids_alphaM1SD], y[ids_alphaM1SD], z[ids_alphaM1SD], s[ids_alphaM1SD], alphaSD[1]))
        else:
            log.info('alphaM1SD point out of profile\n')
        if ids_alphaM2SD:
            outfile.write('alphaM2SD point (x,y,z,s) in [m]:(%.2f,%.2f,%.2f,%.2f) and angle in [°] : %.2f\n' % (
                x[ids_alphaM2SD], y[ids_alphaM2SD], z[ids_alphaM2SD], s[ids_alphaM2SD], alphaSD[2]))
        else:
            log.info('alphaM2SD point out of profile\n')
=======
        outfile.write('alphaM1SD point (x,y,z,s) in [m]:(%.2f,%.2f,%.2f,%.2f) and angle in [°] : %.2f\n' % (
            x[ids_alphaM1SD], y[ids_alphaM1SD], z[ids_alphaM1SD], s[ids_alphaM1SD], alphaSD[1]))
        outfile.write('alphaM2SD point (x,y,z,s) in [m]:(%.2f,%.2f,%.2f,%.2f) and angle in [°] : %.2f\n' % (
            x[ids_alphaM2SD], y[ids_alphaM2SD], z[ids_alphaM2SD], s[ids_alphaM2SD], alphaSD[2]))
>>>>>>> 6175d0dc6ad423c90515b6bb4a0afc37978790f6
        outfile.write('alphaP1SD point (x,y,z,s) in [m]:(%.2f,%.2f,%.2f,%.2f) and angle in [°] : %.2f\n' % (
            x[ids_alphaP1SD], y[ids_alphaP1SD], z[ids_alphaP1SD], s[ids_alphaP1SD], alphaSD[0]))
