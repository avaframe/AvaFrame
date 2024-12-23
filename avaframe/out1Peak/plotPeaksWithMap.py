import rasterio
import rasterio.plot
from rasterio.enums import Resampling
import matplotlib.pyplot as plt
import contextily as ctx
import numpy as np
import geopandas as gpd
from rasterio.crs import CRS
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.offsetbox import TextArea, DrawingArea, OffsetImage, AnnotationBbox


def plotASCwithWmts(ASCFile, cropShape=None):
    # Open the GeoTIFF file
    with rasterio.open(ASCFile, 'r+') as src:
        src.crs = CRS.from_epsg(31287)
        # Read the first band
        data = src.read(1, resampling=Resampling.bilinear)

    # Replace no-data values with NaN for better visualization
    data = np.where(data == src.nodata, np.nan, data)
    data = np.where(data < 0.001, np.nan, data)

    # Plot the GeoTIFF
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.set_title("GeoTIFF/ASC with WMTS Background")

    # Plot the GeoTIFF
    modRes = rasterio.plot.show(data, transform=src.transform, ax=ax, alpha=1, cmap="viridis", zorder=5)

    if cropShape:
        focus = gpd.read_file(cropShape)
        focus.plot(ax=ax, zorder=12, edgecolor="red", linewidth=2, facecolor="none")
        extent = focus.total_bounds
        ax.set_xlim(extent[0], extent[2])
        ax.set_ylim(extent[1], extent[3])

    # Add a WMTS background (orthophoto from a WMTS server)
    # ctx.add_basemap(ax, crs=src.crs, source=ctx.providers.BasemapAT.grau, zorder=4)
    ctx.add_basemap(ax, crs=src.crs, source=ctx.providers.BasemapAT.orthofoto, zorder=4)

    # Plot second time, otherwise basemap is on top, in case zorder is not working... (FSO)
    # rasterio.plot.show(data, transform=src.transform, ax=ax, alpha=0.9, cmap="viridis")

    # North Arrow
    x, y, arrow_length = 0.95, 0.95, 0.1
    ax.annotate('N', xy=(x, y), xytext=(x, y - arrow_length),
                arrowprops=dict(facecolor='black', width=1, headwidth=5),
                ha='center', va='center', fontsize=10,
                xycoords=ax.transAxes, zorder=6)

    # Colorbar
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    im = modRes.get_images()[0]
    fig.colorbar(im, cax=cax)

    # logo = plt.imread('logo.png')

    # newax = fig.add_axes([0.21, 0.21, 0.22, 0.22], anchor='SE')
    # newax.imshow(logo, alpha=0.6, zorder=20)
    # newax.axis('off')

    # imagebox = OffsetImage(logo, zoom=0.8)
    # ab = AnnotationBbox(imagebox, (0, 0), zorder=15)
    # ax.add_artist(ab)

    # ax.figure.figimage(logo, 200, 200, alpha=.5, zorder=15, resize=True)

    # Show the plot
    plt.show()


baseAscFile = 'pft.asc'
focusFile = 'focus.shp'

# Replace with the path to your GeoTIFF file
plotASCwithWmts(baseAscFile, focusFile)
# plotASCwithWmts(baseAscFile)
