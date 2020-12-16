"""
    make plots for in1Data module

    This file is part of Avaframe.

"""

# load python modules
import os
import numpy as np
import matplotlib.pyplot as plt
import glob


def plotDist(workingDir, CDF, a, b, c, cfg, flagShow):
    """ plot the CDF """

    # Generate plot of CDF
    halfLine = np.zeros(len(CDF)) + 0.5
    fig = plt.figure()
    plt.plot(halfLine, 'k--')
    plt.plot(CDF)
    plt.title('Pert CDF, a=%.2f, b=%.2f, c=%.2f m' % (a,b,c))
    plt.xlabel('support')
    plt.ylabel('CDF')
    plt.grid()

    if flagShow:
        plt.show()

    # save fig
    fig.savefig(os.path.join(workingDir, 'CDF_%s.png' % cfg['name']))
    plt.close('all')


def plotSample(workingDir, sample, cfg, flagShow):
    """ Generate bar plot of sample values """

    xSteps = np.arange(len(sample))
    fig = plt.figure()
    plt.title('Sampled values from distribution')
    plt.ylabel('Sampled values')
    plt.bar(xSteps,sample)

    if flagShow:
        plt.show()

    # save fig
    fig.savefig(os.path.join(workingDir, 'samples_%s.png' % cfg['name']))
    plt.close('all')


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

    if flagShow:
        plt.show()

    # save fig
    fig.savefig(os.path.join(workingDir, 'PDFcompare_%s.png' % cfg['name']))
    plt.close('all')


def plotEmpCDF(workingDir, CDF, CDFEmp, xSample, cfg, methodAbbr, flagShow):
    """ make a comparison plot of desired CDF and empirical CDF of sample """

    halfLine = np.zeros(len(CDF)) + 0.5
    x = np.linspace(float(cfg['a']), float(cfg['c']), len(CDF))
    fig = plt.figure()
    plt.title('Desired CDF vs. retrieved sample´s CDF- %s' % methodAbbr)
    plt.plot(x, halfLine, 'k--')
    plt.plot(x, CDF, 'g', label='Desired CDF')
    plt.plot(xSample, CDFEmp, 'b*', label='Actual CDF')

    if flagShow:
        plt.show()

    # save fig
    fig.savefig(os.path.join(workingDir, 'CDFcompare%s_%s.png' % (methodAbbr, cfg['name'])))
    plt.close('all')


def plotEmpPDF(workingDir, PDF, sampleVect, cfg, flagShow):
    """ make a comparison plot of desired CDF and empirical CDF of sample """

    x = np.linspace(float(cfg['a']), float(cfg['c']), len(PDF))
    fig = plt.figure()
    plt.title('Desired PDF vs. sample histogram')
    bins = int(int(cfg['sampleSize'])*0.25)
    plt.hist(sampleVect, bins, density=True, label='sample')
    plt.plot(x, PDF, 'k--', label='Desired PDF')
    plt.xlabel('sample values')
    plt.legend()

    if flagShow:
        plt.show()

    # save fig
    fig.savefig(os.path.join(workingDir, 'PDFcompare_%s.png' % (cfg['name'])))
    plt.close('all')


def plotECDF(workingDir, CDF, sample, cfg, methodAbbr, flagShow):
    """ make a comparison plot of desired CDF and empirical CDF of sample """

    halfLine = np.zeros(len(CDF)) + 0.5
    x = np.linspace(float(cfg['a']), float(cfg['c']), len(CDF))
    fig, ax = plt.subplots()
    plt.suptitle('Desired CDF vs. retrieved sample´s CDF- %s' % methodAbbr)
    ax.plot(x, halfLine, 'k--')
    ax.plot(x, CDF, 'g', label='Desired CDF')
    # plot the cumulative histogram
    n_bins = int(int(cfg['sampleSize'])* 0.25)
    n, bins, patches = ax.hist(sample, n_bins, density=True, histtype='step',
                           cumulative=True, label='Empirical')
    if flagShow:
        plt.show()

    # save fig
    fig.savefig(os.path.join(workingDir, 'CDFcompare%s_%s.png' % (methodAbbr, cfg['name'])))
    plt.close('all')


def plotHistRunout(outDir, Lrun, ECDF, xSample, cfg, pLim, distType, flagShow):
    """ plot a histogram of a sample and the empirical CDF, works for one or more samples """

    # determine x-axis extent
    start = cfg['start']
    end = cfg['end']
    cs = cfg['cs']
    # compute number of bins
    step = int((end - start) / cs)

    # make plot
    fig, ax = plt.subplots()
    plt.suptitle('Sample histogram')
    binsEd = np.linspace(start+cs, end, step)
    alphaval = 1.0
    if len(Lrun.shape) == 1:
        ax.hist(Lrun, binsEd, density=True, label='sample')
    else:
        for m in range(Lrun.shape[1]):
            alphaval = 1. / (m+1)
            ax.hist(Lrun[:,m], binsEd, density=True, alpha=alphaval, label='sample')
    ax.set_xlim(start, end)
    tickNames = ax.get_xticks().tolist()
    for m in range(len(tickNames)):
        tickNames[m] = tickNames[m] - cfg['start']
    ax.set_xticklabels(tickNames)
    ax.set_xlabel('sample values')
    ax.legend(loc='upper left')

    # add second y axis
    ax2 = ax.twinx()
    if len(xSample.shape) == 1:
        ax2.plot(xSample, ECDF, 'k--', label='ECDF')
    else:
        for m in range(xSample.shape[1]):
            ax2.plot(xSample[:,m], ECDF[:,m], 'k--', label='ECDF')
    ax2.legend(loc='upper right')

    # shpw plot
    if flagShow:
        plt.show()

    # save fig
    fig.savefig(os.path.join(outDir, 'sampleHistogram_pLim%s_dist%s.png' % (pLim, distType)))
    plt.close('all')
