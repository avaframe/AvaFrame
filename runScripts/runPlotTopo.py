"""
    Run script for plotting DEMs in module in3Utils
"""

import logging

# Local imports
from avaframe.in3Utils import cfgUtils, logUtils
from avaframe.out3Plot import outTopo as oT

# log file name; leave empty to use default runLog.log
logName = 'plotTopo'

# Load avalanche directory from general configuration file
cfgMain = cfgUtils.getGeneralConfig()
avalancheDir = cfgMain['MAIN']['avalancheDir']

plotList = ['data/avaPyramid',
            'data/avaParabola',
            'data/avaHelix',
            'data/avaHelixChannel',
            'data/avaHockeyChannel',
            'data/avaHockeySmall',
            'data/avaBowl',
            'data/avaInclinedPlane',
            'data/avaFlatPlane',
            'data/avaWog',
            'data/avaGar',
            'data/avaHit',
            'data/avaKot',
            'data/avaMal',
            'data/avaAlr'
            ]

# Start logging
log = logUtils.initiateLogger(avalancheDir, logName)
log.info('MAIN SCRIPT')

for topo in plotList:
    cfgMain['MAIN']['avalancheDir'] = topo
    # Plot topogrpahy
    oT.plotDEM3D(cfgMain)
