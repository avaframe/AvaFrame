"""
    make plots for in1Data module

"""

# load python modules
import os
import numpy as np
import matplotlib.pyplot as plt
import glob

import avaframe.out3Plot.plotUtils as pU


# flagLim if x axis limited
flagLim = False


def plotDist(workingDir, CDF, a, b, c, cfg, flagShow):
    """ plot the CDF """

    # Generate plot of CDF
    halfLine = np.zeros(len(CDF)) + 0.5
    fig = plt.figure()
    plt.plot(halfLine, 'k--')
    plt.plot(CDF)
    plt.title('%s CDF, a=%.3f, b=%.3f, c=%.3f m' % (cfg['distType'], a, b, c))
    plt.xlabel('support')
    plt.ylabel('CDF')
    plt.grid()

    # save and or plot fig
    plotName = 'CDF_%s_%s' % (cfg['name'], cfg['distType'])
    plotPath = pU.saveAndOrPlot({'pathResult': workingDir}, plotName, fig)


def plotSample(workingDir, sample, cfg, flagShow):
    """ Generate bar plot of sample values """

    xSteps = np.arange(len(sample))
    fig = plt.figure()
    plt.title('Sampled values from distribution')
    plt.ylabel('Sampled values')
    plt.bar(xSteps, sample)

    # save and or plot fig
    plotName = 'samples_%s_%s' % (cfg['name'], cfg['distType'])
    plotPath = pU.saveAndOrPlot({'pathResult': workingDir}, plotName, fig)


def plotSamplePDF(workingDir, sampleVect, kdeDict, PDF, cfg, flagShow):
    """ make comparison plot of desired PDF and approximated sample PDF """
    x = np.linspace(0, 1, 10000)
    fig, ax1 = plt.subplots()
    fig.suptitle('Desired PDF vs. retrieved sample´s PDF')
    ax1.plot(x, PDF, 'k', linewidth=4, label='Desired PDF')
    ax1.legend(loc='upper left')
    ax2 = ax1.twiny()
    for key in kdeDict:
        ax2.plot(sampleVect, kdeDict[key](sampleVect), label=cfg['bwMethod'])
    ax2.legend(loc='upper right')

    # save and or plot fig
    plotName = 'PDFcompare%s_%s' % (cfg['name'], cfg['distType'])
    plotPath = pU.saveAndOrPlot({'pathResult': workingDir}, plotName, fig)


def plotEmpCDF(workingDir, CDF, CDFEmp, xSample, cfg, methodAbbr, flagShow, x=''):
    """ make a comparison plot of desired CDF and empirical CDF of sample """

    halfLine = np.zeros(len(CDF)) + 0.5
    if len(x) == 0:
        x = np.linspace(float(cfg['a']), float(cfg['c']), len(CDF))
    fig = plt.figure()
    plt.title('Desired CDF vs. retrieved sample´s CDF- %s' % methodAbbr)
    plt.plot(x, halfLine, 'k--')
    plt.plot(x, CDF, 'g', label='Desired CDF')
    plt.plot(xSample, CDFEmp, 'b*', label='Actual CDF')

    # save and or plot fig
    plotName = 'CDFcompare%s_%s_%s' % (methodAbbr, cfg['name'], cfg['distType'])
    plotPath = pU.saveAndOrPlot({'pathResult': workingDir}, plotName, fig)


def plotEmpPDF(workingDir, PDF, sampleVect, cfg, flagShow, x=''):
    """ make a comparison plot of desired CDF and empirical CDF of sample """

    if len(x) == 0:
        x = np.linspace(float(cfg['a']), float(cfg['c']), len(PDF))
    fig = plt.figure()
    plt.title('Desired PDF vs. sample histogram')
    bins = int(int(cfg['sampleSize'])*0.25)
    plt.hist(sampleVect, bins, density=True, label='sample')
    plt.plot(x, PDF, 'k--', label='Desired PDF')
    plt.xlabel('sample values')
    plt.legend()

    # save and or plot fig
    plotName = 'PDFcompare_%s_%s' % (cfg['name'], cfg['distType'])
    plotPath = pU.saveAndOrPlot({'pathResult': workingDir}, plotName, fig)


def plotECDF(workingDir, CDF, sample, cfg, methodAbbr, flagShow):
    """ make a comparison plot of desired CDF and empirical CDF of sample """

    halfLine = np.zeros(len(CDF)) + 0.5
    x = np.linspace(float(cfg['a']), float(cfg['c']), len(CDF))
    fig, ax = plt.subplots()
    plt.suptitle('Desired CDF vs. retrieved sample´s CDF- %s' % methodAbbr)
    ax.plot(x, halfLine, 'k--')
    ax.plot(x, CDF, 'g', label='Desired CDF')
    # plot the cumulative histogram
    n_bins = int(int(cfg['sampleSize']) * 0.25)
    n, bins, patches = ax.hist(sample, n_bins, density=True, histtype='step',
                           cumulative=True, label='Empirical')

    # save and or plot fig
    plotName = 'CDFcompare%s_%s_%s' % (methodAbbr, cfg['name'], cfg['distType'])
    plotPath = pU.saveAndOrPlot({'pathResult': workingDir}, plotName, fig)


def plotAreaShpError(xFeat, yFeat, nParts, pathDict):
    """ plot polygon parts of polygon read from shp file to check if holes

        Parameters
        ------------
        xFeat, yFeat: numpy array
            x, y coordinates of polygon
        nParts: list
            indices of parts of polygon (2 make only an outer polygon as entire number of points is added)
        pathDict: dict
            dictionary with info on outDir, outFileName and title of plot
    """

    # create plot for each lineFeature
    fig = plt.figure(figsize=(pU.figW, pU.figH))
    for indParts, valP in enumerate(nParts[:-1]):
        if indParts == 0:
            plt.plot(xFeat[0:nParts[indParts + 1]], yFeat[0:nParts[indParts + 1]],
                     label=('part %d' % indParts))
        else:
            plt.plot(xFeat[valP:nParts[indParts + 1]], yFeat[valP:nParts[indParts + 1]],
                     label=('part %d' % indParts))

    plt.title(pathDict['title'])
    plt.legend()
    # save and or plot
    pU.saveAndOrPlot(pathDict, pathDict['outFileName'], fig)
