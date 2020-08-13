"""Run script for module com2AB
"""

# Local imports
from avaframe.com2AB import com2AB
from avaframe.out3SimpPlot import outAB
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils

# -----------Required settings-----------------

# TODO: move this to a main cfg file! (alongside flags for intermediate data/results/plots)
# and have function cfgUtils.readGlobalCfg that allows for a local_MainCfg

# Avalanche directory; see doc.avaframe.org for setup
# TODO: full path needed?
avalancheDir = 'data/avaSlide'


# log file name; leave empty to use default runLog.log
logName = 'runCom2AB'

# ---------------------------------------------

# Start logging
log = logUtils.initiateLogger(avalancheDir, logName)
log.info('MAIN SCRIPT')
log.info('Current avalanche: %s', avalancheDir)

# Load all input Parameters from config file
# get the configuration of an already imported module
cfg = cfgUtils.getModuleConfig(com2AB)

cfgSetup = cfg['ABSETUP']
cfgFlags = cfg['FLAGS']

# Dump config to log file
logUtils.writeCfg2Log(cfg,'com2AB')

# Extract input file locations
cfgPath = com2AB.readABinputs(avalancheDir)

log.info("Running com2ABMain model on test case DEM \n %s \n with profile \n %s ",
         cfgPath['demSource'], cfgPath['profileLayer'])

# Read input data for ALPHABETA
dem = com2AB.readRaster(cfgPath['demSource'])
avaPath = com2AB.readAvaPath(cfgPath['profileLayer'], cfgPath['defaultName'], dem['header'])
splitPoint = com2AB.readSplitPoint(cfgPath['splitPointSource'], dem['header'])

# Calculate ALPHABETA
com2AB.com2ABMain(dem, avaPath, splitPoint, cfgPath['saveOutPath'],
                  cfgSetup)


# Analyse/ plot/ write results #
plotFile, writeFile = outAB.writeABpostOut(dem,
                                           avaPath, splitPoint,
                                           cfgPath['saveOutPath'],
                                           cfgFlags)

log.info('Plotted to: %s', plotFile)
log.info('Data written: %s', writeFile)
