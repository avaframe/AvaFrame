import rasterio
import numpy as np
import math
import logging
import pathlib
from rasterio.features import rasterize
from shapely.geometry import shape, mapping
from avaframe.in2Trans.shpConversion import SHP2Array
from avaframe.in1Data import getInput as gI
from avaframe.in3Utils import fileHandlerUtils as fU
import avaframe.in2Trans.rasterUtils as IOf

# Configure logging
logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger(__name__)


def scarpAnalysisMain(cfg, baseDir):
    """Run the scarp analysis using parameters from an .ini cfguration file and input from a directory

    Parameters
    ----------
    cfg : cfgparser
        with all required fields in scarpCfg.ini
    baseDir: str
        path to directory of data to analyze

    Returns
    -------
    None
    """

    # Get input from baseDir and check if all files exist
    if cfg["INPUT"].getboolean("useShapefiles"):

        # Set directories for inputs, outputs and current work
        inputDir = pathlib.Path(baseDir, "Inputs")

        # get the path to the perimeter shapefile
        perimeterShapefilePath, periFile = gI.getAndCheckInputFiles(
            inputDir, "POLYGONS", "scarp perimeter", fileExt="shp", fileSuffix="_perimeter"
        )
        if periFile:
            log.info("Perimeterfile is: %s" % perimeterShapefilePath)
        else:
            log.error("Perimeter shapefile not found in %s", inputDir)

        # get the path to the coordinates shapefile
        coordinatesShapefilePath, coordFile = gI.getAndCheckInputFiles(
            inputDir, "POINTS", "scarp coordinates", fileExt="shp", fileSuffix="_coordinates"
        )
        if coordFile:
            log.info("Coordinate shapefile is: %s" % coordinatesShapefilePath)
        else:
            log.error("Coordinate shapefile not found in %s", inputDir)

    else:
        log.error(
            "Shapefile option not selected. Please set 'useShapefiles' to True in the configuration file. A "
            "raster version is currently not implemented."
        )

    # Initialize feature parameters
    planeFeatures = []
    ellipsoidFeatures = []

    # Read the DEM
    dem = gI.readDEM(baseDir)
    dem["rasterData"] = np.flipud(dem["rasterData"])

    n = dem["header"]["nrows"]
    m = dem["header"]["ncols"]

    periData = readPerimeterSHP(perimeterShapefilePath, dem["header"]["transform"], (n, m))

    SHPdata = SHP2Array(coordinatesShapefilePath)

    log.debug("Feature shapefile data loaded: %s", SHPdata)

    method = cfg["SETTINGS"]["method"].lower()

    if method == "plane":
        # Read required attributes directly from the shapefile's attribute table
        try:
            planesZseed = list(map(float, SHPdata['zseed']))
            planesDip = list(map(float, SHPdata['dip']))
            planesSlope = list(map(float, SHPdata['slopeangle']))
        except KeyError as e:
            raise ValueError(f"Required attribute '{e.args[0]}' not found in shapefile. Make sure 'zseed', 'dip', and 'slope' fields exist.")

        if not (len(planesZseed) == len(planesDip) == len(planesSlope) == SHPdata["nFeatures"]):
            raise ValueError("Mismatch between number of features and extracted plane attributes in the shapefile.")

        if not (len(planesZseed) == len(planesDip) == len(planesSlope) == SHPdata["nFeatures"]):
            raise ValueError(
                "Mismatch between number of shapefile features and plane parameters in the .ini file."
            )

        for i in range(SHPdata["nFeatures"]):
            xSeed = SHPdata["x"][int(SHPdata["Start"][i])]
            ySeed = SHPdata["y"][int(SHPdata["Start"][i])]
            zSeed = planesZseed[i]
            dip = planesDip[i]
            slopeAngle = planesSlope[i]
            planeFeatures.extend([xSeed, ySeed, zSeed, dip, slopeAngle])

        features = ",".join(map(str, planeFeatures))
        log.debug("Plane features extracted and combined: %s", features)

    elif method == "ellipsoid":
        try:
            ellipsoidsMaxDepth = list(map(float, SHPdata['maxdepth']))
            ellipsoidsSemiMajor = list(map(float, SHPdata['semimajor']))
            ellipsoidsSemiMinor = list(map(float, SHPdata['semiminor']))
        except KeyError as e:
            raise ValueError(f"Required attribute '{e.args[0]}' not found in shapefile. Ensure the fields 'maxdepth', 'semimajor', and 'semiminor' exist.")

        if not (
            len(ellipsoidsMaxDepth)
            == len(ellipsoidsSemiMajor)
            == len(ellipsoidsSemiMinor)
            == SHPdata["nFeatures"]
        ):
            raise ValueError(
                "Mismatch between number of shapefile features and ellipsoid parameters in the .ini file."
            )

        for i in range(SHPdata["nFeatures"]):
            xCenter = SHPdata["x"][int(SHPdata["Start"][i])]
            yCenter = SHPdata["y"][int(SHPdata["Start"][i])]
            maxDepth = ellipsoidsMaxDepth[i]
            semiMajor = ellipsoidsSemiMajor[i]
            semiMinor = ellipsoidsSemiMinor[i]
            ellipsoidFeatures.extend([xCenter, yCenter, maxDepth, semiMajor, semiMinor])

        features = ",".join(map(str, ellipsoidFeatures))
        log.debug("Ellipsoid features extracted and combined: %s", features)

    if method == 'plane':
        scarpData = calculateScarpWithPlanes(
            dem["rasterData"], periData, dem["header"]["transform"], features
        )
    elif method == 'ellipsoid':
        scarpData = calculateScarpWithEllipsoids(
            dem["rasterData"], periData, dem["header"]["transform"], features
        )
    else:
        raise ValueError("Unsupported method. Choose 'plane' or 'ellipsoid'.")

    hRelData = dem["rasterData"] - scarpData

    # create output directory and files
    outDir = pathlib.Path(baseDir)
    outDir = outDir / "Outputs" / "com6RockAvalanche" / "scarp"
    fU.makeADir(outDir)

    # Set file names and write
    elevationOut = outDir / "scarpElevation"
    hRelOut = outDir / "scarpHRel"
    IOf.writeResultToRaster(dem["header"], scarpData, elevationOut)
    IOf.writeResultToRaster(dem["header"], hRelData, hRelOut)

def readPerimeterSHP(perimeterShapefilePath, elevTransform, elevShape):
    """Read perimeter from shapefile and return as numpy array.

    Parameters
    ----------
    perimeterShapefilePath : str
        Path to the perimeter shapefile.
    elevTransform : Affine
        transformation of the elevation raster.
    elevShape : tuple
        Shape (height, width) of the elevation raster.

    Returns
    -------
    periData : np.ndarray
        2D array representing the perimeter mask.
    """

    log.debug("Processing perimeter shapefile: %s", perimeterShapefilePath)
    SHPdata = SHP2Array(perimeterShapefilePath)
    log.debug("Perimeter shapefile data loaded: %s", SHPdata)

    shapes = []
    for i in range(SHPdata["nFeatures"]):
        start = int(SHPdata["Start"][i])
        length = int(SHPdata["Length"][i])
        coords = [(SHPdata["x"][j], SHPdata["y"][j]) for j in range(start, start + length)]
        poly = shape({"type": "Polygon", "coordinates": [coords]})
        shapes.append(poly)

    periData = rasterize(
        [(mapping(poly), 1) for poly in shapes],
        out_shape=elevShape,
        transform=elevTransform,
        fill=0,
        all_touched=True,
        dtype=np.uint8,
    )
    log.debug("Perimeter shapefile rasterized.")

    return periData


def calculateScarpWithPlanes(elevData, periData, elevTransform, planes):
    """Calculate the scarp using sliding planes.

    Parameters
    ----------
    elevData : np.ndarray
        The elevation data as a 2D array.
    periData : np.ndarray
        The perimeter data as a 2D array, which defines the extent of the scarp.
    elevTransform : Affine
        The affine transformation matrix of the raster (used to convert pixel to geographic coordinates).
    planes : str
        Comma-separated string defining sliding planes (xseed, yseed, zseed, dip, slope).

    Returns
    -------
    np.ndarray
        A 2D array with the calculated scarp values.
    """
    n, m = elevData.shape
    scarpData = np.zeros_like(elevData, dtype=np.float32)

    planes = list(map(float, planes.split(',')))
    nPlanes = int(len(planes) / 5)

    xSeed = [planes[0]]
    ySeed = [planes[1]]
    zSeed = [planes[2]]
    dip = [planes[3]]
    slope = [planes[4]]

    betaX = [math.tan(math.radians(slope[0])) * math.cos(math.radians(dip[0]))]
    betaY = [math.tan(math.radians(slope[0])) * math.sin(math.radians(dip[0]))]

    for i in range(1, nPlanes):
        xSeed.append(planes[5 * i])
        ySeed.append(planes[5 * i + 1])
        zSeed.append(planes[5 * i + 2])
        dip.append(planes[5 * i + 3])
        slope.append(planes[5 * i + 4])
        betaX.append(math.tan(math.radians(slope[i])) * math.cos(math.radians(dip[i])))
        betaY.append(math.tan(math.radians(slope[i])) * math.sin(math.radians(dip[i])))

    for row in range(n):
        for col in range(m):
            west, north = rasterio.transform.xy(elevTransform, row, col, offset='center')

            scarpVal = zSeed[0] + (north - ySeed[0]) * betaY[0] - (west - xSeed[0]) * betaX[0]
            for k in range(1, nPlanes):
                scarpVal = max(scarpVal, zSeed[k] + (north - ySeed[k]) * betaY[k] - (west - xSeed[k]) * betaX[k])

            if periData[row, col] > 0:
                scarpData[row, col] = min(elevData[row, col], scarpVal)
            else:
                scarpData[row, col] = elevData[row, col]

    return scarpData

def calculateScarpWithEllipsoids(elevData, periData, elevTransform, ellipsoids):
    """Calculate the scarp using sliding circles (rotational ellipsoids).

    Parameters
    ----------
    elevData : np.ndarray
        The elevation data as a 2D array.
    periData : np.ndarray
        The perimeter data as a 2D array, which defines the extent of the scarp.
    elevTransform : Affine
        The affine transformation matrix of the raster (used to convert pixel to geographic coordinates).
    ellipsoids : str
        Comma-separated string defining ellipsoids (x_center, y_center, max_depth, semi_major_axis, semi_minor_axis).

    Returns
    -------
    np.ndarray
        A 2D array with the calculated scarp values.
    """
    n, m = elevData.shape
    scarpData = np.zeros_like(elevData, dtype=np.float32)

    ellipsoids = list(map(float, ellipsoids.split(',')))
    nEllipsoids = int(len(ellipsoids) / 5)

    xCenter = [ellipsoids[0]]
    yCenter = [ellipsoids[1]]
    maxDepth = [ellipsoids[2]]
    semiMajor = [ellipsoids[3]]
    semiMinor = [ellipsoids[4]]

    for i in range(1, nEllipsoids):
        xCenter.append(ellipsoids[5 * i])
        yCenter.append(ellipsoids[5 * i + 1])
        maxDepth.append(ellipsoids[5 * i + 2])
        semiMajor.append(ellipsoids[5 * i + 3])
        semiMinor.append(ellipsoids[5 * i + 4])

    for row in range(n):
        for col in range(m):
            west, north = rasterio.transform.xy(elevTransform, row, col, offset='center')
            scarpVal = elevData[row, col]
            for k in range(nEllipsoids):
                x2 = ((west - xCenter[k]) / semiMajor[k]) ** 2
                y2 = ((north - yCenter[k]) / semiMinor[k]) ** 2
                if x2 + y2 < 1:
                    scarpVal = min(scarpVal, elevData[row, col] - maxDepth[k] * (1 - (x2 + y2)))
            if periData[row, col] > 0:
                scarpData[row, col] = scarpVal
            else:
                scarpData[row, col] = elevData[row, col]

    return scarpData
