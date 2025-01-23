# -*- coding: utf-8 -*-
"""
Created on Tue Jan 21 11:43:18 2025

@author: Domi
"""

import rasterio
import numpy as np
import configparser
from rasterio.features import rasterize
from shapely.geometry import shape, mapping
from in2Trans.shpConversion import SHP2Array
import logging

# Configure logging
logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger(__name__)

def generateMuXsiRasters(configPath):
    """
    Generate raster files for μ and ξ based on input DEM and shapefiles.

    Parameters
    ----------
    configPath : str
        Path to the configuration file.

    Returns
    -------
    None
    """
    # Load configuration
    config = configparser.ConfigParser()
    config.read(configPath)

    demPath = config['INPUT']['dem']
    muShapefile = config['INPUT']['mu_shapefile']
    xsiShapefile = config['INPUT']['xsi_shapefile']
    defaultMu = float(config['DEFAULTS']['default_mu'])
    defaultXsi = float(config['DEFAULTS']['default_xsi'])
    muOutputPath = config['OUTPUT']['mu_raster']
    xsiOutputPath = config['OUTPUT']['xsi_raster']

    # Read DEM
    with rasterio.open(demPath) as demSrc:
        demData = demSrc.read(1)
        demTransform = demSrc.transform
        demCrs = demSrc.crs
        demShape = demData.shape

    # Helper function to rasterize shapefiles
    def rasterizeShapefile(shapefilePath, defaultValue, attributeName):
        if not shapefilePath:
            return np.full(demShape, defaultValue, dtype=np.float32)

        shpData = SHP2Array(shapefilePath)
        shapes = []
        for i in range(shpData['nFeatures']):
            start = int(shpData['Start'][i])
            length = int(shpData['Length'][i])
            coords = [(shpData['x'][j], shpData['y'][j]) for j in range(start, start + length)]
            poly = shape({'type': 'Polygon', 'coordinates': [coords]})
            value = shpData['attributes'][i][attributeName]  # Extract the attribute value
            shapes.append((mapping(poly), value))

        rasterData = rasterize(
            shapes,
            out_shape=demShape,
            transform=demTransform,
            fill=defaultValue,
            all_touched=True,
            dtype=np.float32
        )
        return rasterData

    # Generate μ and ξ rasters
    log.info("Rasterizing μ shapefile.")
    muRaster = rasterizeShapefile(muShapefile, defaultMu, "mu")

    log.info("Rasterizing ξ shapefile.")
    xsiRaster = rasterizeShapefile(xsiShapefile, defaultXsi, "xsi")

    # Save output rasters
    def saveRaster(outputPath, data):
        with rasterio.open(
            outputPath,
            'w',
            driver='GTiff',
            height=data.shape[0],
            width=data.shape[1],
            count=1,
            dtype=data.dtype,
            crs=demCrs,
            transform=demTransform,
        ) as dst:
            dst.write(data, 1)

    log.info("Saving μ raster to %s", muOutputPath)
    saveRaster(muOutputPath, muRaster)

    log.info("Saving ξ raster to %s", xsiOutputPath)
    saveRaster(xsiOutputPath, xsiRaster)

    log.info("Raster generation completed.")
