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

# create local logger
# change log level in calling module to DEBUG to see log messages
log = logging.getLogger(__name__)



def writeReport(avaDir, reportDict):
    """ Write a report for simulation """

    outDir = os.path.join(avaDir, 'Outputs', 'com1DFA', 'reports')

    # Load dictionary with report data
    parameterD = reportDict['simParameters']
    simName = reportDict['simName']

    with open(os.path.join(outDir, '%s.md' % simName), 'w') as pfile:
        pfile.write('# %s \n' % reportDict['headerLine'])
        pfile.write('### Simulation name: *%s* \n'  % simName)
        pfile.write('| Parameters | Values | \n')
        pfile.write('| ---------- | ------ | \n')
        for value in parameterD:
            pfile.write('| %s | %s | \n' % (value, parameterD[value]))
        pfile.write(' \n')

        if parameterD['release area'] != '':
            releaseD = reportDict['releaseArea']
            pfile.write('### Release areas \n')
            for value in releaseD:
                pfile.write('| %s ' % value)
            pfile.write('| \n')
            for value in releaseD:
                pfile.write('| ----------')
            pfile.write('| \n')
            for value in releaseD:
                pfile.write('| %s ' % releaseD[value])
            pfile.write('| \n')
            pfile.write(' \n')

        if parameterD['entrainment area'] != '':
            entrainD = reportDict['entrainmentArea']
            pfile.write('### Entrainment areas \n')
            for value in entrainD:
                pfile.write('| %s ' % value)
            pfile.write('| \n')
            for value in entrainD:
                pfile.write('| ----------')
            pfile.write('| \n')
            for value in entrainD:
                pfile.write('| %s ' % entrainD[value])
            pfile.write('| \n')
            pfile.write(' \n')


        if parameterD['resistance area'] != '':
            resistD = reportDict['resistanceArea']
            pfile.write('### Resistance areas \n')
            lenRelVal = len(resistD)
            for value in resistD:
                pfile.write('| %s ' % value)
            pfile.write('| \n')
            for value in resistD:
                pfile.write('| ----------')
            pfile.write('| \n')
            for value in resistD:
                pfile.write('| %s ' % resistD[value])
            pfile.write('| \n')
            pfile.write(' \n')

        if 'images' in reportDict:
            imageD = reportDict['images']
            pfile.write('### Images \n')
            for value in imageD:
                pfile.write('##### Figure:   %s \n'  % value)
                pfile.write('![%s](%s) \n' % (value, imageD[value]))

        if 'text' in reportDict:
            textD = reportDict['text']
            pfile.write('### Additional Info \n')
            for value in textD:
                pfile.write('##### Topic:   %s \n'  % value)
                pfile.write('%s \n' % (textD[value]))
