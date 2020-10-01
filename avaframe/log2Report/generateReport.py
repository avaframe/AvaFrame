"""
    Generate a markdown report for data provided in dictionary

    This file is part of Avaframe.
"""

# Load modules
import os
import glob
import logging
import numpy as np
import shutil
from avaframe.in3Utils import fileHandlerUtils as fU

# create local logger
# change log level in calling module to DEBUG to see log messages
log = logging.getLogger(__name__)

def writeColumns(dict, key, pfile):
    """ Create block with columns for each key and value pair of dict """

    # Write header for block - key of the reporDict
    pfile.write('### %s \n' % key)
    for value in dict:
        pfile.write('| %s ' % value)
    pfile.write('| \n')
    for value in dict:
        pfile.write('| ----------')
    pfile.write('| \n')
    for value in dict:
        pfile.write('| %s ' % dict[value])
    pfile.write('| \n')
    pfile.write(' \n')


def writeReportFile(reportD, pfile):
    """ Write text to file """

    for key in reportD:

        # HEADER BLOCK
        if key == 'headerLine':
            pfile.write('# %s \n' % reportD['headerLine'])
        if key == 'simName':
            pfile.write('### Simulation name: *%s* \n'  % reportD['simName'])

        # SIMULATION BLOCK
        if key == 'simParameters':
            pfile.write('| Parameters | Values | \n')
            pfile.write('| ---------- | ------ | \n')
            for value in reportD[key]:
                pfile.write('| %s | %s | \n' % (value, reportD[key][value]))
            pfile.write(' \n')

        # INPUTS BLOCK
        if key == 'Release area' or key == 'Entrainment area' or key == 'Resistance area':
            writeColumns(reportD[key], key, pfile)

        # IMAGE BLOCK
        if key == 'images':
            pfile.write('### Images \n')
            for value in reportD['images']:
                pfile.write('##### Figure:   %s \n'  % value)
                pfile.write('![%s](%s) \n' % (value, reportD['images'][value]))

        # TEXT BLOCK
        if key == 'text':
            pfile.write('### Additional Info \n')
            for value in reportD['text']:
                pfile.write('##### Topic:   %s \n'  % value)
                pfile.write('%s \n' % (reportD['text'][value]))


def writeReport(avaDir, reportDictList, plotDict, reportOneFile=True):
    """ Write a report for simulation """

    outDir = os.path.join(avaDir, 'Outputs', 'com1DFA', 'reports')

    if reportOneFile:
        # Start writing markdown style report
        with open(os.path.join(outDir, 'fullSimulationReport.md'), 'w') as pfile:

            # Loop through all simulations
            for reportD in reportDictList:

                # extract additional info from log file
                reportD['simParameters'].update(fU.extractParameterInfo(avaDir, reportD['simName']))

                # add plot info to general report Dict
                reportD['images'] = plotDict[reportD['simName']]

                # Write report file
                writeReportFile(reportD, pfile)


    else:

        # Loop through all simulations
        for reportD in reportDictList:

            # extract additional info from log file
            reportD['simParameters'].update(fU.extractParameterInfo(avaDir, reportD['simName']))

            # add plot info to general report Dict
            reportD['images'] = plotDict[reportD['simName']]

            # Start writing markdown style report
            with open(os.path.join(outDir, '%s.md' % reportD['simName']), 'w') as pfile:

                # Write report file
                writeReportFile(reportD, pfile)
