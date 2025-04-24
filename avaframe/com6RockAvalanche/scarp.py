import rasterio
import numpy as np
import math
import configparser
import logging
from avaframe.in2Trans.shpConversion import SHP2Array
from rasterio.features import rasterize
from shapely.geometry import shape, mapping

# Configure logging
logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger(__name__)

def runScarpAnalysis(configFile):
    # TODO: rename this to scarpMain or similar
    """Run the scarp analysis using parameters from an .ini configuration file.

    Parameters
    ----------
    configFile : str
        Path to the configuration file in .ini format.

    Returns
    -------
    None
    """
    config = configparser.ConfigParser()
    config.read(configFile)

    # Read input parameters from the configuration file
    elevation = config['INPUT']['elevation']
    perimeterInput = config['INPUT']['perimeter']
    shapefilePath = config['INPUT'].get('shapefile', '').strip()
    perimeterShapefilePath = config['INPUT'].get('perimeter_shapefile', '').strip()
    features = config['INPUT'].get('features', '')

    # TODO: outputs should go to Outputs/com6RockAvalanche, so these configs should be removed
    elevScarp = config['OUTPUT']['elevscarp']
    hRelease = config['OUTPUT']['hrelease']

    method = config['SETTINGS']['method'].lower()

    # Initialize feature parameters
    planeFeatures = []
    ellipsoidFeatures = []

    with rasterio.open(elevation) as elevSrc:
        elevData = elevSrc.read(1)
        elevTransform = elevSrc.transform
        elevCRS = elevSrc.crs
        n, m = elevData.shape

    periData = readPerimeter(perimeterInput, perimeterShapefilePath, elevTransform, (n, m), elevCRS)

    if shapefilePath:
        SHPdata = SHP2Array(shapefilePath)
        log.debug("Feature shapefile data loaded: %s", SHPdata)

        if method == 'plane':
            planesZseed = list(map(float, config['SETTINGS']['planes_zseed'].split(',')))
            planesDip = list(map(float, config['SETTINGS']['planes_dip'].split(',')))
            planesSlope = list(map(float, config['SETTINGS']['planes_slope'].split(',')))

            if not (len(planesZseed) == len(planesDip) == len(planesSlope) == SHPdata['nFeatures']):
                raise ValueError("Mismatch between number of shapefile features and plane parameters in the .ini file.")

            for i in range(SHPdata['nFeatures']):
                xSeed = SHPdata['x'][int(SHPdata['Start'][i])]
                ySeed = SHPdata['y'][int(SHPdata['Start'][i])]
                zSeed = planesZseed[i]
                dip = planesDip[i]
                slopeAngle = planesSlope[i]
                planeFeatures.extend([xSeed, ySeed, zSeed, dip, slopeAngle])

            features = ','.join(map(str, planeFeatures))
            log.debug("Plane features extracted and combined: %s", features)

        elif method == 'ellipsoid':
            ellipsoidsMaxDepth = list(map(float, config['SETTINGS']['ellipsoids_max_depth'].split(',')))
            ellipsoidsSemiMajor = list(map(float, config['SETTINGS']['ellipsoids_semi_major'].split(',')))
            ellipsoidsSemiMinor = list(map(float, config['SETTINGS']['ellipsoids_semi_minor'].split(',')))

            if not (len(ellipsoidsMaxDepth) == len(ellipsoidsSemiMajor) == len(ellipsoidsSemiMinor) == SHPdata['nFeatures']):
                raise ValueError("Mismatch between number of shapefile features and ellipsoid parameters in the .ini file.")

            for i in range(SHPdata['nFeatures']):
                xCenter = SHPdata['x'][int(SHPdata['Start'][i])]
                yCenter = SHPdata['y'][int(SHPdata['Start'][i])]
                maxDepth = ellipsoidsMaxDepth[i]
                semiMajor = ellipsoidsSemiMajor[i]
                semiMinor = ellipsoidsSemiMinor[i]
                ellipsoidFeatures.extend([xCenter, yCenter, maxDepth, semiMajor, semiMinor])

            features = ','.join(map(str, ellipsoidFeatures))
            log.debug("Ellipsoid features extracted and combined: %s", features)

    else:
        if not features:
            raise ValueError("No features provided. Please specify features in the .ini file or provide a shapefile.")

    scarpData = np.zeros_like(elevData, dtype=np.float32)
    hRelData = np.zeros_like(elevData, dtype=np.float32)

    if method == 'plane':
        scarpData = calculateScarpWithPlanes(elevData, periData, elevTransform, features)
    elif method == 'ellipsoid':
        scarpData = calculateScarpWithEllipsoids(elevData, periData, elevTransform, features)
    else:
        raise ValueError("Unsupported method. Choose 'plane' or 'ellipsoid'.")

    hRelData = elevData - scarpData
    saveRaster(elevScarp, scarpData, elevTransform, elevCRS)
    saveRaster(hRelease, hRelData, elevTransform, elevCRS)

def readPerimeter(perimeterInput, perimeterShapefilePath, elevTransform, elevShape, elevCRS):
    """Read perimeter from raster or shapefile and return as numpy array.

    Parameters
    ----------
    perimeterInput : str
        Path to the perimeter raster file.
    perimeterShapefilePath : str
        Path to the perimeter shapefile.
    elevTransform : Affine
        transformation of the elevation raster.
    elevShape : tuple
        Shape (height, width) of the elevation raster.
    elevCrs : rasterio.crs.CRS
        Coordinate Reference System of the elevation raster.

    Returns
    -------
    periData : np.ndarray
        2D array representing the perimeter mask.
    """
    if perimeterShapefilePath:
        log.debug("Processing perimeter shapefile: %s", perimeterShapefilePath)
        SHPdata = SHP2Array(perimeterShapefilePath)
        log.debug("Perimeter shapefile data loaded: %s", SHPdata)

        shapes = []
        for i in range(SHPdata['nFeatures']):
            start = int(SHPdata['Start'][i])
            length = int(SHPdata['Length'][i])
            coords = [(SHPdata['x'][j], SHPdata['y'][j]) for j in range(start, start + length)]
            poly = shape({'type': 'Polygon', 'coordinates': [coords]})
            shapes.append(poly)

        periData = rasterize(
            [(mapping(poly), 1) for poly in shapes],
            out_shape=elevShape,
            transform=elevTransform,
            fill=0,
            all_touched=True,
            dtype=np.uint8
        )
        log.debug("Perimeter shapefile rasterized.")

    elif perimeterInput.lower().endswith(('.tif', '.tiff', '.asc')):
        log.debug("Reading perimeter raster: %s", perimeterInput)
        with rasterio.open(perimeterInput) as periSrc:
            periData = periSrc.read(1)
        log.debug("Perimeter raster loaded.")
    else:
        raise ValueError("Unsupported perimeter input format. Provide a raster (.tif, .asc) or a shapefile (.shp).")

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

def saveRaster(filename, data, transform, crs):
    # TODO: use the functions we have in avaframe.in3Utils.rasterUtils.writeResultToRaster
    """Save raster data to a specified format.

    Parameters
    ----------
    output_path : str
        The file path where the output raster will be saved.
    data : np.ndarray
        The 2D array of data to be saved.
    elevTransform : Affine
        The affine transformation matrix of the elevation raster.
    elevCRS : rasterio.crs.CRS
        The Coordinate Reference System of the elevation raster.

    Returns
    -------
    None
    """
    with rasterio.open(
        filename,
        'w',
        driver='GTiff',
        height=data.shape[0],
        width=data.shape[1],
        count=1,
        dtype=data.dtype,
        crs=crs,
        transform=transform,
    ) as dst:
        dst.write(data, 1)
    log.debug("Raster saved: %s", filename)
