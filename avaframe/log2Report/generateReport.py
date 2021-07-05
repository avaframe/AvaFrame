"""
    Generate a markdown report for data provided in dictionary

"""

# Load modules
import os
import logging

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
    """ Create markdown report with blocks, tables, list according to type key

        Parameters
        ----------
        reportD : dict
            report dictionary
        pfile : file object
            file where markdown format report is written to

    """

    # Loop through keys and perform action according to value found in type
    for key in reportD:
        for subKey in reportD[key]:
            # HEADER BLOCK
            # Title
            if reportD[key][subKey] == 'title':
                for value in reportD[key]:
                    if value != 'type':
                        pfile.write('# %s \n' % reportD[key][value])
                        pfile.write(' \n')
                        break

            # Avalanche name
            if reportD[key][subKey] == 'avaName':
                for value in reportD[key]:
                    if value != 'type':
                        pfile.write('### Avalanche directory: *%s* \n' % reportD[key][value])
                        pfile.write(' \n')

            # Simulation name
            if reportD[key][subKey] == 'simName':
                for value in reportD[key]:
                    if value != 'type':
                        pfile.write('### Simulation name: *%s* \n' % reportD[key][value])
                        pfile.write(' \n')
            # Time info
            if reportD[key][subKey] == 'time':
                for value in reportD[key]:
                    if value != 'type':
                        pfile.write('### Date: %s \n' % reportD[key][value])
                        pfile.write(' \n')
                        break
            # PARAMETER BLOCK
            # Table listing all the key : value pairs in rows
            if reportD[key][subKey] == 'list':
                pfile.write('### %s \n' % key)
                pfile.write(' \n')
                pfile.write('| Parameters | Values | \n')
                pfile.write('| ---------- | ------ | \n')
                for value in reportD[key]:
                    if value != 'type':
                        pfile.write('| %s | %s | \n' % (value, reportD[key][value]))
                pfile.write(' \n')
            # Table listing the key : value pairs in columns
            if reportD[key][subKey] == 'columns':
                pfile.write('### %s \n' % key)
                pfile.write(' \n')
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
                pfile.write(' \n')
                for value in reportD[key]:
                    if value != 'type':
                        pfile.write('##### Figure:   %s \n' % value)
                        pfile.write(' \n')
                        pfile.write('![%s](%s) \n' % (value, reportD[key][value]))
                        pfile.write(' \n')
                        pfile.write(' \n')

            # TEXT BLOCK
            if reportD[key][subKey] == 'text':
                pfile.write('### %s \n' % key)
                pfile.write(' \n')
                for value in reportD[key]:
                    if value != 'type':
                        pfile.write('##### Topic:  %s \n' % value)
                        pfile.write(' \n')
                        pfile.write('%s \n' % (reportD[key][value]))
                        pfile.write(' \n')


def writeReport(outDir, reportDictList, cfgFLAGS, plotDict=''):
    """ Write a report in markdown format for simulations, saved to outDir

        Parameters
        ----------
        outDir : str
            path to output directory
        reportDictList : list
            list of report dictionaries from simulations
        cfgFLAGS : dict
            configuration dictionary
        plotDict : dict
            optional dictionary with info on plots that shall be included in report

    """

    if cfgFLAGS.getboolean('reportOneFile'):
        # Start writing markdown style report
        with open(os.path.join(outDir, 'fullSimulationReport.md'), 'w') as pfile:

            # Loop through all simulations
            for reportD in reportDictList:

                if plotDict != '' and ('simName' in reportD):
                    # add plot info to general report Dict
                    reportD['Simulation Results'] = plotDict[reportD['simName']['name']]
                    reportD['Simulation Results'].update({'type': 'image'})
                # Write report file
                writeReportFile(reportD, pfile)

    else:

        # Loop through all simulations
        for reportD in reportDictList:

            if plotDict != '':
                # add plot info to general report Dict
                reportD['Simulation Results'] = plotDict[reportD['simName']['name']]
                reportD['Simulation Results'].update({'type': 'image'})

            # Start writing markdown style report
            with open(os.path.join(outDir, '%s.md' % reportD['simName']), 'w') as pfile:

                # Write report file
                writeReportFile(reportD, pfile)
