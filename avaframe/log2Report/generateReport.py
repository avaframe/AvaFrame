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

    # Loop through keys and perform action according to value found in type
    for key in reportD:
        for subKey in reportD[key]:

            # HEADER BLOCK
            # Title
            if reportD[key][subKey] == 'title':
                for value in reportD[key]:
                    if value != 'type':
                        pfile.write('# %s \n' % reportD[key][value])
                        break
            # Simulation name
            if reportD[key][subKey] == 'simName':
                for value in reportD[key]:
                    if value != 'type':
                        pfile.write('### Simulation name: *%s* \n' % reportD[key][value])

            # PARAMETER BLOCK
            # Table listing all the key : value pairs in rows
            if reportD[key][subKey] == 'list':
                pfile.write('### %s \n' % key)
                pfile.write('| Parameters | Values | \n')
                pfile.write('| ---------- | ------ | \n')
                for value in reportD[key]:
                    if value != 'type':
                        pfile.write('| %s | %s | \n' % (value, reportD[key][value]))
                pfile.write(' \n')
            # Table listing the key : value pairs in columns
            if reportD[key][subKey] == 'columns':
                pfile.write('### %s \n' % key)
                for value in reportD[key]:
                    if value != 'type':
                        pfile.write('| %s ' % value)
                pfile.write('| \n')
                for value in reportD[key]:
                    if value != 'type':
                        pfile.write('| ----------')
                pfile.write('| \n')
                for value in reportD[key]:
                    if value != 'type':
                        pfile.write('| %s ' % reportD[key][value])
                pfile.write('| \n')
                pfile.write(' \n')

            # IMAGE BLOCK
            if reportD[key][subKey] == 'image':
                pfile.write('### %s \n' % key)
                for value in reportD[key]:
                    if value != 'type':
                        pfile.write('##### Figure:   %s \n' % value)
                        pfile.write('![%s](%s) \n' % (value, reportD[key][value]))

            # TEXT BLOCK
            if reportD[key][subKey] == 'text':
                pfile.write('### %s \n' % key)
                for value in reportD[key]:
                    if value != 'type':
                        pfile.write('##### Topic:  %s \n' % value)
                        pfile.write('%s \n' % (reportD[key][value]))


def writeReport(outDir, reportDictList, cfgFLAGS, plotDict=''):
    """ Write a report for simulation """

    if cfgFLAGS.getboolean('reportOneFile'):
        # Start writing markdown style report
        with open(os.path.join(outDir, 'fullSimulationReport.md'), 'w') as pfile:

            # Loop through all simulations
            for reportD in reportDictList:

                # extract additional info from log file
                #reportD['simParameters'].update(fU.extractParameterInfo(avaDir, reportD['simName']))

                if plotDict != '':
                    # add plot info to general report Dict
                    reportD['Simulation Results'] = plotDict[reportD['simName']['name']]
                    reportD['Simulation Results'].update({'type' : 'image'})
                # Write report file
                writeReportFile(reportD, pfile)

    else:

        # Loop through all simulations
        for reportD in reportDictList:

            # extract additional info from log file
            #reportD['simParameters'].update(fU.extractParameterInfo(avaDir, reportD['simName']))

            if plotDict != '':
                # add plot info to general report Dict
                reportD['Simulation Results'] = plotDict[reportD['simName']['name']]
                reportD['Simulation Results'].update({'type' : 'image'})

            # Start writing markdown style report
            with open(os.path.join(outDir, '%s.md' % reportD['simName']), 'w') as pfile:

                # Write report file
                writeReportFile(reportD, pfile)
