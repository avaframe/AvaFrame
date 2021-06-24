"""

This is a simple function for computing a mask of all peak files of one parameter that exceed a particular threshold

"""

import os
import numpy as np
import logging
from matplotlib import pyplot as plt
import pathlib

import avaframe.out3Plot.plotUtils as pU
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import fileHandlerUtils as fU
import avaframe.in2Trans.ascUtils as IOf

# create local logger
# change log level in calling module to DEBUG to see log messages
log = logging.getLogger(__name__)


def probAnalysis(avaDir, cfg, module, parametersDict='', inputDir=''):
    """ Compute propability map of a given set of simulation result exceeding a particular threshold and save to outDir

        Parameters
        ----------
        avaDir: str
            path to avalanche directory
        cfg : dict
            configuration read from ini file of probAna function
        module
            computational module that was used to run the simulations
        parametersDict: dict
            dictionary with simulation parameters to filter simulations
        inputDir : str
            optional - path to directory where data that should be analysed can be found, required if not in module results
    """

    # get filename of module
    modName = pathlib.Path(module.__file__).stem

    # set output directory
    outDir = pathlib.Path(avaDir, 'Outputs', 'ana4Stats')
    fU.makeADir(outDir)

    # fetch all result files and filter simulations according to parametersDict
    simNameList = cfgUtils.filterSims(avaDir, parametersDict, specDir=inputDir)
    if inputDir == '':
        inputDir = pathlib.Path(avaDir, 'Outputs', modName, 'peakFiles')
        flagStandard = True
        peakFiles, _ = fU.makeSimDict(inputDir, simID='simHash', avaDir=avaDir)
    else:
        peakFiles, _ = fU.makeSimDict(inputDir, avaDir=avaDir)

    # get header info from peak files - this should be the same for all peakFiles
    header = IOf.readASCheader(peakFiles['files'][0])
    cellSize = header.cellsize
    nRows = header.nrows
    nCols = header.ncols
    xllcenter = header.xllcenter
    yllcenter = header.yllcenter
    noDataValue = header.noDataValue

    # Initialise array for computations
    probSum = np.zeros((nRows, nCols))
    count = 0

    # Loop through peakFiles and compute probability
    for m in range(len(peakFiles['names'])):

        # only take simulations that match filter criteria from parametersDict
        if peakFiles['simName'][m] in simNameList:
            # Load peak field for desired peak field parameter
            if peakFiles['resType'][m] == cfg['GENERAL']['peakVar']:

                # Load data
                fileName = peakFiles['files'][m]
                data = np.loadtxt(fileName, skiprows=6)
                dataLim = np.zeros((nRows, nCols))

                log.info('File Name: %s , simulation parameter %s ' % (fileName, cfg['GENERAL']['peakVar']))

                # Check if peak values exceed desired threshold
                dataLim[data > float(cfg['GENERAL']['peakLim'])] = 1.0
                probSum = probSum + dataLim
                count = count + 1

    # Create probability map ranging from 0-1
    probMap = probSum / count
    unit = pU.cfgPlotUtils['unit%s' % cfg['GENERAL']['peakVar']]
    log.info('probability analysis performed for peak parameter: %s and a peak value threshold of: %s %s' % (cfg['GENERAL']['peakVar'], cfg['GENERAL']['peakLim'], unit))
    log.info('%s peak fields added to analysis' % count)

    # # Save to .asc file
    avaName = os.path.basename(avaDir)
    outFileName = '%s_probMap%s.asc' % (avaName, cfg['GENERAL']['peakLim'])
    outFile = outDir / outFileName
    IOf.writeResultToAsc(header, probMap, outFile)
