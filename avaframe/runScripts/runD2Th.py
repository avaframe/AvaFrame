"""
Run conversion from depth to thickness using DEM
"""

import pathlib

# Local imports
import avaframe.in2Trans.ascUtils as IOf
from avaframe.in1Data import getInput as gI
import avaframe.in2Trans.transfromFields as tF
import avaframe.out3Plot.outTransformPlots as oT
from avaframe.in3Utils import fileHandlerUtils as fU
from avaframe.in3Utils import cfgUtils


def runD2Th(avaDir, comMod, resType, num, profileAxis, profileIndex):
    # fetch dem
    dem = gI.readDEM(avaDir)

    # directory with peak files
    inDir = pathlib.Path(avaDir, "Outputs", comMod, "peakFiles")
    resFiles = list(inDir.glob("*_%s.asc" % resType))

    # loop over resType files found
    for rF in resFiles:
        # read depth field to dict
        depthField = IOf.readRaster(rF)

        # convert depth to thickness using dem
        thicknessDict, depthRasterResized, slopeAngleField = tF.convertDepthToThickness(depthField, dem, num=8)

        # create output directory
        outDir = pathlib.Path(avaDir, "Outputs", comMod, "peakFiles", "transformed")
        fU.makeADir(outDir)
        pName = rF.stem.split("_pft")[0] + "transformed" + "_%s.asc" % resType

        # create plot
        oT.plotDepthToThickness(
            depthRasterResized, thicknessDict['rasterData'], slopeAngleField, profileAxis, profileIndex, outDir, pName
        )

        # write thickness to file
        outFile = outDir / pName
        IOf.writeResultToAsc(thicknessDict["header"], thicknessDict['rasterData'], outFile, flip=True)


if __name__ == "__main__":
    # +++++++++REQUIRED+++++++++++++
    # log file name; leave empty to use default runLog.log
    logName = "runDepthToThickness"
    comMod = "com1DFA"
    resType = "pft"
    num = 6
    profileAxis = "x"
    profileIndex = None
    # ++++++++++++++++++++++++++++++

    # fetch input directory
    cfgMain = cfgUtils.getGeneralConfig()
    avaDir = cfgMain["MAIN"]["avalancheDir"]

    # call conversion
    runD2Th(avaDir, comMod, resType, num, profileAxis, profileIndex)
