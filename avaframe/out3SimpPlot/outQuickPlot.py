"""

This is a simple function for a quick plot of datasets and option for comparison
between two datasets of identical shape. Also plots cross and longprofiles.

set desired values in

This file is part of Avaframe.

"""

import matplotlib.pyplot as plt
from avaframe.in3Utils import fileHandlerUtils as fU
import numpy as np
import os
import seaborn as sns
import logging
import shutil
import glob


# create local logger
# change log level in calling module to DEBUG to see log messages
log = logging.getLogger(__name__)


def prepareData(avaDir, inputDir):
    """ Prepare data files to be plotted """

    # Input directory and load input datasets
    datafiles = glob.glob(inputDir+os.sep + '*.asc')

    # make directory of input file names and datasets
    data = {'files' : [], 'values' : [], 'names' : [], 'resType' : [], 'simType' : []}
    for m in range(len(datafiles)):
        data['files'].append(datafiles[m])
        data['values'].append( np.loadtxt(datafiles[m], skiprows=6))
        name = os.path.splitext(os.path.basename(datafiles[m]))[0]
        data['names'].append(name)
        data['simType'].append(name.split('_')[1])
        data['resType'].append(name.split('_')[3])

    return data


def getRefData(avaDir, outputDir, suffix):
    """ Grab reference data to be plotted with quick plot function """

    # Input directory and load input datasets
    ava = avaDir.split(os.sep)[1]
    refDir = os.path.join('..', 'benchmarks', ava)

    dataRefFiles = glob.glob(refDir+os.sep + '*%s.asc' % suffix)

    for files in dataRefFiles:
        shutil.copy2(files, outputDir)

    log.info('Reference files copied from directory: %s' % refDir)


def quickPlot(avaDir, suffix, com1DFAOutput, simName):
    """ Plot two raster datasets of same dimension """

    # Create required directories
    workDir = os.path.join(avaDir, 'Work', 'out3SimplPlot')
    fU.makeADir(workDir)

    # Setup input from com1DFA
    fU.getDFAData(avaDir, com1DFAOutput, workDir, suffix)

    # Get data from reference run
    getRefData(avaDir, workDir, suffix)

    # prepare data
    data = prepareData(avaDir, workDir)

    print('In here suffix is: ', suffix, ' and sim Name', simName)
    # sava file, name, and data to dictionary
    indSuffix = []
    for m in range(len(data['files'])):
        if data['resType'][m] == suffix and data['simType'][m] == simName:
            indSuffix.append(m)

    data1 = data['values'][indSuffix[0]]
    data2 = data['values'][indSuffix[1]]
    ny = data1.shape[0]
    nx = data1.shape[1]
    log.info('dataset1: %s' % data['files'][indSuffix[0]])
    log.info('dataset2: %s' % data['files'][indSuffix[1]])

    # Location of Profiles
    ny_loc = int(nx *0.5)
    nx_loc = int(ny *0.5)

    # Plot data
    # Figure 1 shows the result parameter data
    sns.set()
    fig = plt.figure()
    fig.set_size_inches(15, 5)
    ax1 = fig.add_subplot(131)
    ax2 = fig.add_subplot(132)
    ax3 = fig.add_subplot(133)
    cmap = sns.cubehelix_palette(8, start=.5, rot=-.75, as_cmap=True)
    sns.heatmap(data1, cmap=cmap, ax=ax1)
    sns.heatmap(data2, cmap=cmap, ax=ax2)
    cmapdiv = sns.color_palette("RdBu_r")
    sns.heatmap(data1-data2, cmap=cmapdiv, ax=ax3)
    ax1.set_title('%s' % data['names'][indSuffix[0]])
    ax2.set_title('%s' % data['names'][indSuffix[1]])
    ax3.set_title('Difference ref-sim')

    # Fgiure 2 cross and lonprofile
    fig, ax = plt.subplots(ncols=2, figsize=(15, 5))
    ax[0].plot(data1[:, ny_loc], 'k', linewidth=4, label='Reference')
    ax[0].plot(data2[:, ny_loc], 'b--', label='Simulation')
    ax[0].set_xlabel('Location across track [nrows]')
    ax[0].set_ylabel('Result parameter %s' % suffix, fontsize=12)
    ax[0].set_title('Cross profile at y =  %d' % ny_loc)
    plt.legend()
    ax[1].plot(data1[nx_loc, :], 'k', linewidth=4, label='Reference')
    ax[1].plot(data2[nx_loc, :], 'b--', label='Simulation')
    ax[1].set_xlabel('Location along track [ncols]')
    ax[1].set_ylabel('Result parameter %s' % suffix, fontsize=12 )
    ax[1].set_title('Long profile at x =  %d' % nx_loc)

    ax[0].legend()
    ax[1].legend()

    plt.show()
