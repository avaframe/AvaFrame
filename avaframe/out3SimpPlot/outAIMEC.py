import sys
import os
import logging
import math
import numpy as np
import scipy as sp
import copy
from matplotlib import pyplot as plt
import matplotlib
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
from matplotlib import cm
from matplotlib.image import NonUniformImage


# Local imports
import avaframe.in2Trans.shpConversion as shpConv
import avaframe.in2Trans.geoTrans as geoTrans
import avaframe.in3Utils.ascUtils as IOf

# create local logger
log = logging.getLogger(__name__)

def visu_transfo(raster_transfo, input_data, cfgPath, cfgFlags):
    """
    Plot and save the domain transformation figure
    """
    # read paths
    pathResult = cfgPath['pathResult']
    project_name = cfgPath['dirName']
    # read rasterdata
    sourceData = input_data['sourceData']
    header = sourceData['header']
    xllcenter = header.xllcenter
    yllcenter = header.yllcenter
    cellsize = header.cellsize
    rasterdata = sourceData['rasterData']
    # read avaPath with scale
    Avapath = input_data['Avapath']
    x_path = Avapath['x']*cellsize+xllcenter
    y_path = Avapath['y']*cellsize+yllcenter
    # read domain boundarries with scale
    DB = input_data['DB']
    DB_x_l = DB['DB_x_l']*cellsize+xllcenter
    DB_x_r = DB['DB_x_r']*cellsize+xllcenter
    DB_y_l = DB['DB_y_l']*cellsize+yllcenter
    DB_y_r = DB['DB_y_r']*cellsize+yllcenter

    figure_width = 2*10
    figure_height = 2*5
    lw = 1

    fig = plt.figure(figsize=(figure_width, figure_height), dpi=150)

#    for figure: referenz-simulation bei p_lim=1
    ax1 = plt.subplot(121)
    indBeta = raster_transfo['indBeta']
    xx = raster_transfo['x'][indBeta]
    yy = raster_transfo['y'][indBeta]
    new_rasterdata = rasterdata
    masked_array = new_rasterdata #np.ma.masked_where(np.isnan(new_rasterdata), new_rasterdata)
    cmap = copy.copy(matplotlib.cm.jet)
    cmap.set_under(color='w')
    cmap.set_bad(color='k')

    n, m = np.shape(new_rasterdata)
    x = np.arange(m)*cellsize+xllcenter
    y = np.arange(n)*cellsize+yllcenter
    im = NonUniformImage(ax1, extent=[x.min(), x.max(), y.min(), y.max()], cmap=cmap)
    # im.set_interpolation('bilinear')
    im.set_clim(vmin=0.000000001)
    im.set_data(x, y, masked_array)
    ref1 = ax1.images.append(im)
    cbar = ax1.figure.colorbar(im, ax=ax1, use_gridspec=True)
    plt.autoscale(False)
    ref0 = plt.plot(xx, yy, 'ro', label='Beta point')
    ref2 = plt.plot(x_path, y_path,
                    'b-', linewidth=lw, label='flow path')
    ref3 = plt.plot(DB_x_l, DB_y_l,
                    'g-', linewidth=lw, label='domain')
    ref3 = plt.plot(DB_x_r, DB_y_r,
                    'g-', linewidth=lw, label='domain')
    ref3 = plt.plot([DB_x_l, DB_x_r], [DB_y_l, DB_y_r],
                    'g-', linewidth=lw, label='domain')
    refs = [ref0[0], ref2[0], ref3[0]]

    labels = ['Beta point', 'flow path', 'domain']
    ax1.title.set_text('XY Domain')
    ax1.legend(refs, labels, loc=0)
    ax1.set_xlim([x.min(), x.max()])
    ax1.set_ylim([y.min(), y.max()])
    ax1.set_xlabel('x [m]')
    ax1.set_ylabel('y [m]')

    ax2 = plt.subplot(122)
    ax2.title.set_text('sl Domain \n Black = out of raster')
    isosurf = copy.deepcopy(input_data['aval_data'])
    l_coord = raster_transfo['l']
    s_coord = raster_transfo['s']
    ref1 = ax2.axhline(y=s_coord[indBeta], color='r', linewidth=1,
                       linestyle='-', label='Beta point')
    masked_array = isosurf #np.ma.array(isosurf,mask=np.isnan(isosurf))
    im = NonUniformImage(ax2, extent=[l_coord.min(), l_coord.max(),
                                      s_coord.min(), s_coord.max()], cmap=cmap)
    # im.set_interpolation('bilinear')
    im.set_clim(vmin=0.000000001)
    im.set_data(l_coord, s_coord, masked_array)
    ref0 = ax2.images.append(im)
    cbar = ax2.figure.colorbar(im, ax=ax2, use_gridspec=True)
    cbar.ax.set_ylabel('peak pressure [kPa]')
    ax2.set_xlim([l_coord.min(), l_coord.max()])
    ax2.set_ylim([s_coord.min(), s_coord.max()])
    ax2.set_xlabel('l [m]')
    ax2.set_ylabel('s [m]')
    ax2.legend(loc=0)

    fig.tight_layout()
    if cfgFlags.getboolean('plotFigure'):
        plt.show()
    else:
        plt.ioff()
    if cfgFlags.getboolean('savePlot'):
        outname_fin = ''.join([pathResult, '/pics/', project_name,
                               '_domTransfo', '.pdf'])
        if not os.path.exists(os.path.dirname(outname_fin)):
            os.makedirs(os.path.dirname(outname_fin))
        fig.savefig(outname_fin, transparent=True)

    plt.close(fig)


def visu_runout(raster_transfo, inputPlot, cfgPath, cfgFlags):
    """
    Plot and save the Peak Pressure  distribution after coord transfo
    """
    # read paths
    pathResult = cfgPath['pathResult']
    project_name = cfgPath['dirName']
    # read data
    s_coord = raster_transfo['s']
    l_coord = raster_transfo['l']
    indBeta = raster_transfo['indBeta']
    rasterArea = raster_transfo['rasterArea']
    dataPressure = inputPlot['dataPressure']
    rasterdataPres = dataPressure[0]  # ana3AIMEC.makeRasterAverage(dataPressure)
    runout = inputPlot['runout']
    runout_mean = inputPlot['runout_mean']
    p_lim = inputPlot['pressureLimit']

    p_mean = inputPlot['p_mean']
    p_median = inputPlot['p_median']
    p_percentile = inputPlot['p_percentile']

    figure_width = 3*5
    figure_height = 3*3

    fig = plt.figure(figsize=(figure_width, figure_height), dpi=150)
    ax1 = plt.subplot(121)
    ax1.title.set_text('Peak Pressure 2D plot for the reference')
    ref1 = ax1.axhline(y=s_coord[indBeta], color='k', linewidth=1,
                       linestyle='-', label='Beta point')
    ref1 = ax1.axhline(y=np.max(runout), color='r', linewidth=1,
                       linestyle='-', label='runout max')
    ref2 = ax1.axhline(y=np.average(runout), color='y', linewidth=1,
                       linestyle='-', label='runout mean')
    ref3 = ax1.axhline(y=np.min(runout), color='g', linewidth=1,
                       linestyle='-', label='runout min')
    # ref3 = ax1.plot(np.zeros(np.shape(s_coord)), s_coord,'.r', linewidth=0.1)
    isosurf = copy.deepcopy(rasterdataPres)
    xx, yy = np.meshgrid(l_coord, s_coord)
    masked_array = np.ma.masked_where(isosurf == 0, isosurf)
    cmap = copy.copy(matplotlib.cm.jet)
    cmap.set_bad('w', 1.)
    im = NonUniformImage(ax1, extent=[xx.min(), xx.max(), yy.min(), yy.max()], cmap=cmap)
    # im.set_interpolation('bilinear')
    im.set_data(l_coord, s_coord, masked_array)
    ref0 = ax1.images.append(im)
    cbar = ax1.figure.colorbar(im, ax=ax1, use_gridspec=True)
    cbar.ax.set_ylabel('peak pressure [kPa]')
    plt.autoscale(False)
    ax1.set_xlim([xx.min(), xx.max()])
    ax1.set_ylim([yy.min(), yy.max()])
    ax1.set_xlabel('l [m]')
    ax1.set_ylabel('s [m]')
    ax1.legend(loc=0)

    ax2 = plt.subplot(122)
    ax2.title.set_text('Peak Pressure distribution along the path between runs')
    ax2.fill_betweenx(s_coord, p_percentile[2], p_percentile[0],
                      facecolor=[.8, .8, .8], alpha=0.5, label='quantiles')
    ref1 = mpatches.Patch(alpha=0.5, color=[.8, .8, .8])
    ref2 = ax2.plot(p_median, s_coord, color='r', linewidth=2, label='median')
    ref3 = ax2.plot(p_mean, s_coord, color='b', linewidth=1, label='mean')
    # ref3 = mlines.Line2D([], [], color='b', linewidth=2)
    ax2.set_ylabel('s [m]')
    ax2.set_ylim([yy.min(), yy.max()])
    ax2.set_xlim(auto=True)
    ax2.set_xlabel('Pmax(s) [kPa]')
    ax2.legend(loc=0)

    fig.tight_layout()

    if cfgFlags.getboolean('savePlot'):
        outname_fin = ''.join([pathResult, '/pics/', project_name, '_dptr',
                               str(int(p_lim)), '_slComparison', '.pdf'])
        if not os.path.exists(os.path.dirname(outname_fin)):
            os.makedirs(os.path.dirname(outname_fin))
        fig.savefig(outname_fin, transparent=True)

    if cfgFlags.getboolean('plotFigure'):
        plt.show()
    else:
        plt.ioff()

    plt.close(fig)


def result_write(cfgPath, cfgSetup, resAnalysis):
    """
    This function is used as standart output file

    example path: 'log/doublestep_5_5m_calibration_results.txt'
    example header: 'runout, maxvel, rho, mu, tau_0, R_s^0, kappa, R, B, q \n'

    INPUT: data, path, header
    """

    project_name = cfgPath['project_name']
    path_name = cfgPath['path_name']
    dem_name = os.path.basename(cfgPath['demSource'])
    data_name = [os.path.basename(name) for name in cfgPath['pressurefileList']]
    domainWidth = cfgSetup['domainWidth']
    pressureLimit = cfgSetup['pressureLimit']

    runout = resAnalysis['runout']
    AMPP = resAnalysis['AMPP']
    MMPP = resAnalysis['MMPP']
    AMD = resAnalysis['AMD']
    MMD = resAnalysis['MMD']
    deltaH = resAnalysis['deltaH']
    elevRel = resAnalysis['elevRel']
    relMass = resAnalysis['relMass']
    entMass = resAnalysis['entMass']
    growthIndex = resAnalysis['growthIndex']
    growthGrad = resAnalysis['growthGrad']

    legend = ['fileNr', 'runout', 'elevRel', 'deltaH', 'AMPP',
              'MMPP', 'entMass', 'growthIndex', 'AMD', 'MMD']
    resfile = [runout, elevRel, deltaH, AMPP, MMPP, entMass, growthIndex, AMD, MMD]

    header = ''.join(['project_name: ',  project_name, '\n',
                      'path: ', path_name, '\n',
                      'dhm: ', dem_name, '\n',
                      'domain_width: ', str(domainWidth), '\n',
                      'pressure_limit: ', str(pressureLimit), '\n',
                      'release_mass: ', str(relMass[0]), '\n'])

    outname = ''.join([cfgPath['pathResult'], os.path.sep,
                       'Results_pl', str(pressureLimit),
                       '_w', str(domainWidth), '.txt'])


#    chekc if folder exists / create
    if not os.path.exists(os.path.dirname(outname)):
        os.makedirs(os.path.dirname(outname))

    output = resfile
    log.info('write output file: %s' % outname)
    fid = open(outname, 'w')
    fid.write(header)
    # write table legend
    for j in range(len(legend)):
        fid.write('{:<12s}'.format(legend[j]))
    fid.write('\n')
    # write table values
    for i in range(len(output[0])):
        tmp = os.path.basename(data_name[i])
        name = os.path.splitext(tmp)[0]
        fid.write('{:<12s}'.format(name))
        for j in range(len(output)):
            try:
                fid.write('{:<12.3f}'.format(output[j][i]))
            except:
                fid.write('{:<12}'.format('NaN'))
        fid.write('\n')
    fid.close()

    log.info('File written: %s' % outname)


def colorvar(k, k_end, colorflag, disp=0):
    """
    jt colorvariation editor - JT 2012
    determine how color changes from runnumber = 1 to runnumber = runlength
    input: runnumber,runlength,colorflag,verbose
    output: [R G B]

    possible colorflags:
    'ry': red to yellow
    'bb': blue to light blue
    'pw': pink to white
    'kg' black to green
    'bp' blue to pink
    'pb' pink to blue
    'gy' green to yellow
    'cw' cyan to white
    'kr' black to red
    'gb' green to blue
    'rp' red to pink
    'yw' yellow to white
    'kb' black to blue
    'kc' black to cyan
    'kp' black to pink
    'kw' black to white
    """

    colors = {
        'ry': [1., k/k_end, 0.],  # rot zu gelb
        'bb': [0., k/k_end, 1.],  # blau zu hellbalu
        'pw': [0.8, k/k_end, 0.8],  # pink zu weiss
        'kg': [0., k/k_end, 0.],  # schwarz zu gruen
        'bp': [k/k_end, 0., 1.],  # blau zu pink
        'pb': [1.-k/k_end, 1., 0.],  # blau zu pink
        'gy': [k/k_end, 1., 0.],  # green zu yellow
        'cw': [k/k_end, 1., 1.],  # cyan zu weiss
        'kr': [k/k_end, 1., 1.],  # black to red
        'gb': [0., 1., k/k_end],  # gruen zu blau
        'rp': [1., 0., k/k_end],  # rot zu pink
        'yw': [1., 1., k/k_end],  # yellow to white
        'kb': [0., 0., k/k_end],  # black tp blue
        'kc': [0., k/k_end, k/k_end],  # black zu cyan
        'kp': [k/k_end, 0., k/k_end],  # black zu pink
        'kw': [1.-k/k_end, 1.-k/k_end, 1.-k/k_end]  # black to white
    }

    colornames = {
        'ry': 'red to yellow',
        'bb': 'blue to cyan',
        'pw': 'pink to white',
        'kg': 'black to green',
        'bp': 'blue to pink',
        'pb': 'blue to pink',
        'gy': 'green to yellow',
        'cw': 'cyan to white',
        'kr': 'black to red',
        'gb': 'green to blue',
        'rp': 'rot to pink',
        'yw': 'yellow to white',
        'kb': 'black to blue',
        'kc': 'black to cyan',
        'kp': 'black to pink',
        'kw': 'black to white'
    }

    if colorflag.lower() in colors:
        farbe = colors.get(colorflag.lower())
        if k == 0:
            log.info('Color is: %s' % colornames.get(colorflag.lower()))
    else:
        farbe = [0, 0, 0]
        if k == 0:
            log.info('Color is black')

    return farbe


def result_visu(cfgPath, resAnalysis, doku, GI, dpp_threshold):
    """
    Visualize results in a nice way
    Jan-Thomas Fischer BFW 2010-2012
    AK BFW 2014-2015
    """

    fnames = cfgPath['pressurefileList']
    rasterSource = cfgPath['demSource']
    ProfileLayer = cfgPath['profileLayer']
    outpath = cfgPath['pathResult']
    DefaultName = cfgPath['project_name']

    runout = resAnalysis['runout']
    mean_max_dpp = resAnalysis['AMPP']
    max_max_dpp = resAnalysis['MMPP']

    cvar = ['ry', 'bb', 'pw', 'gy']
    colorflag = cvar[0]

    figure_width = 7*2
    figure_height = 4*2
    fs = 20
    mks = 10
    lw = 2
    # includes flag for y axis -
    # 1 = rddp
    # 2 = frontal shape
    # 3 = groth index
    # 4 = runout for intrapraevent
    # 5 = pressure data
    flag = 3
    if (len(fnames) > 100):
        plot_density = 1
    else:
        plot_density = 0

    if flag == 1:
        log.info('Visualizing pressure data')
        tipo = 'rapp'
        data = mean_max_dpp / mean_max_dpp[0]
        yaxis_label = 'rAPP [-]'
        ytick_increment = 0.25
        ymax = 3
    elif flag == 2:
        log.info('Visualizing EGU growth index data')
        tipo = 'GI'
        data = GI
        yaxis_label = 'growth index [GI]'
        ytick_increment = 2
    elif flag == 3:
        log.info('Visualizing pressure data')
        tipo = 'rmpp'
        data = max_max_dpp / max_max_dpp[0]
        yaxis_label = 'rMPP [-]'
        ytick_increment = 0.1
        ymax = max(data[1:])+(max(data[1:])-min(data[1:]))*0.1
        ymin = min(data[1:])-(max(data[1:])-min(data[1:]))*0.1
    else:
        log.error('Wrong flag')
        return None

    # read data
    dem = IOf.readRaster(rasterSource)
    header = dem['header']
    xllcenter = header.xllcenter
    yllcenter = header.yllcenter
    cellsize = header.cellsize

    rasterdata = dem['rasterData']

    Avapath = shpConv.readLine(ProfileLayer, DefaultName, dem['header'])
    AvaProfile, SplitPoint = geoTrans.prepareLine(dem, Avapath, distance=10)
    x_path = AvaProfile['x']
    y_path = AvaProfile['y']
    z_path = AvaProfile['z']
    s_path = AvaProfile['s']

    xlim_prof_axis = max(s_path) + 50

    # Final result diagram - z_profile+data
    fig = plt.figure(figsize=(figure_width, figure_height), dpi=300)

    markers = ['+', 'o', 'x', '*', 's', 'd', '^', 'v', '>', '<', 'p', 'h', '.',
               '^', 'v', '>', '<', 'p', 'h', '.']
    mk = 0

#    show flow path
    ax1 = fig.add_subplot(111)
#    plt.xlim([0, xlim_prof_axis])
#    plt.ylim([0, math.ceil(max(data)+0.25)])
#    plt.ylim([0, ymax])
#    plt.yticks(np.arange([0, math.ceil(max(data)+0.25), ytick_increment]))
    ax1.set_ylabel(yaxis_label, color='b', fontsize=2*fs)
    ax1.set_xlabel(''.join(['s [m] - runout with ', str(dpp_threshold),
                            ' kPa threshold']), color='black', fontsize=2*fs)
    if plot_density:  # estimate 2D histogram --> create pcolormesh
        nbins = 100
        H, xedges, yedges = np.histogram2d(runout, data, bins=nbins)
        H = np.flipud(np.rot90(H))
        Hmasked = np.ma.masked_where(H == 0, H)
        data_density = plt.pcolormesh(xedges, yedges, Hmasked, cmap=cm.Blues)
#        data_density = plt.pcolormesh(xedges, yedges, Hmasked, cmap=cm.cool)
        cbar = plt.colorbar(data_density, orientation='horizontal')
        cbar.ax.set_ylabel('Counts')
    ax2 = ax1.twinx()
    ax2.set_ylabel('z [m]', color='g', fontsize=2*fs)
    ax2.plot(s_path, z_path, color='green', label='path', linestyle='--', linewidth=2*lw)
    plt.xlim([0, xlim_prof_axis])
    plt.ylim([math.floor(min(z_path)/10)*10, math.ceil(max(z_path)/10)*10])
    if not plot_density:
        for k in range(len(runout)):
            topo_name = fnames[k].split('/')[-1]
            pfarbe = colorvar(float(k), len(runout), colorflag)
            if k == 0:
                ax1.plot(runout[k], data[k], marker='+',
                         markersize=2*mks, color='g', label=topo_name)
    #            plt.yticks(np.arange([0,5000,250]))
                # Make the y-tick labels of first axes match the line color.
                for tl in ax1.get_yticklabels():
                    tl.set_color('b')
            else:
                ax1.plot(runout[k], data[k], label=topo_name, marker=markers[mk],
                         markersize=mks, color=pfarbe, linewidth=lw)
            mk = mk+1
            if mk == len(markers):
                mk = 1
    plt.grid('on')

    pro_name = fnames[0].split('/')[-3]
    outname_fin = ''.join([outpath, '/pics/', pro_name, '_dptr',
                           str(int(dpp_threshold)), '_', tipo, '.pdf'])

    if not os.path.exists(os.path.dirname(outname_fin)):
        os.makedirs(os.path.dirname(outname_fin))
    fig.savefig(outname_fin, transparent=True)

    plt.close(fig)

    # Final result diagram - roc-plots
    rTP = (np.array(doku[0]) / (float(doku[0][0]) + float(doku[1][0]))).astype(float)
    rFP = (np.array(doku[2]) / (float(doku[2][0]) + float(doku[3][0]))).astype(float)


#    rFP = (np.array(doku[2]) / (float(doku[0][0]) + float(doku[1][0]))).astype(float)

    fig = plt.figure(figsize=(figure_width, figure_height), dpi=300)

    mk = 0
    ax1 = fig.add_subplot(111)
    ax1.set_ylabel('True positive rate', fontsize=2*fs)
    ax1.set_xlabel('False positive rate', fontsize=2*fs)
    if plot_density:  # estimate 2D histogram --> create pcolormesh
        nbins = 100
        H, xedges, yedges = np.histogram2d(rFP, rTP, bins=nbins)
        H = np.flipud(np.rot90(H))
        Hmasked = np.ma.masked_where(H == 0, H)
#        data_density = plt.pcolormesh(xedges, yedges, Hmasked, cmap=cm.Blues)
        data_density = plt.pcolormesh(xedges, yedges, Hmasked, cmap=cm.cool)
        cbar = plt.colorbar(data_density, orientation='horizontal')
        cbar.ax.set_ylabel('hit rate density')
    if not plot_density:
        for k in range(len(rTP)):
            topo_name = fnames[k].split('/')[-1]
            pfarbe = colorvar(float(k), len(rTP), colorflag)
            ax1.plot(rFP[k], rTP[k], label=topo_name, marker=markers[mk],
                     markersize=mks, color=pfarbe, linewidth=lw)
            mk = mk+1
            if mk == len(markers):
                mk = 0
    plt.xlim([0, max(1, max(rFP))])
    plt.ylim([0, 1])
    plt.grid('on')

    outname_fin = ''.join([outpath, '/pics/', pro_name, '_dptr',
                           str(int(dpp_threshold)), '_ROC.pdf'])

    if not os.path.exists(os.path.dirname(outname_fin)):
        os.makedirs(os.path.dirname(outname_fin))
    fig.savefig(outname_fin, transparent=True)

    plt.close(fig)

    return
