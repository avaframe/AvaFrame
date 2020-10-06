"""Run script for module com2AB
"""

# Local imports
from avaframe.com2AB import com2AB
from avaframe.out3Plot import outAB
import avaframe.in2Trans.shpConversion as shpConv
import avaframe.in3Utils.ascUtils as IOf
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils


# log file name; leave empty to use default runLog.log
logName = 'runCom2AB'

# Load avalanche directory from general configuration file
cfgMain = cfgUtils.getGeneralConfig()
avalancheDir = cfgMain['MAIN']['avalancheDir']

# Start logging
log = logUtils.initiateLogger(avalancheDir, logName)
log.info('MAIN SCRIPT')
log.info('Current avalanche: %s', avalancheDir)

# Load all input Parameters from config file
# get the configuration of an already imported module
# write config to log file
cfg = cfgUtils.getModuleConfig(com2AB)

cfgSetup = cfg['ABSETUP']
cfgFlags = cfg['FLAGS']

# Extract input file locations
cfgPath = com2AB.readABinputs(avalancheDir)

log.info("Running com2ABMain model on DEM \n \t %s \n \t with profile \n \t %s ",
         cfgPath['demSource'], cfgPath['profileLayer'])

# Read input data for ALPHABETA
dem = IOf.readRaster(cfgPath['demSource'])
avaPath = shpConv.readLine(cfgPath['profileLayer'], cfgPath['defaultName'], dem)
splitPoint = shpConv.readPoints(cfgPath['splitPointSource'], dem)

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
