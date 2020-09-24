import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.colors as mcolors
from matplotlib.image import NonUniformImage

from colorspace.colorlib import HCL
from colorspace.CVD import CVD
from colorspace.CVD import desaturate
from colorspace import specplot, sequential_hcl

def hex_to_rgb(value):
    '''
    Converts hex to rgb colours
    value: string of 6 characters representing a hex colour.
    Returns: list length 3 of RGB values'''
    value = value.strip("#") # removes hash symbol if present
    lv = len(value)
    return tuple(int(value[i:i + lv // 3], 16) for i in range(0, lv, lv // 3))


def rgb_to_dec(value):
    '''
    Converts rgb to decimal colours (i.e. divides each value by 256)
    value: list (length 3) of RGB values
    Returns: list (length 3) of decimal values'''
    return [v/256 for v in value]

def get_continuous_cmap(hex_list, float_list=None):
    ''' creates and returns a listed color map.
        If float_list is not provided, colour map graduates linearly between each color in hex_list.
        If float_list is provided, each color in hex_list is mapped to the respective location in float_list.

        Parameters
        ----------
        hex_list: list of hex code strings
        float_list: list of floats between 0 and 1, same length as hex_list. Must start with 0 and end with 1.

        Returns
        ----------
        colour map'''
    rgb_list = [rgb_to_dec(hex_to_rgb(i)) for i in hex_list]
    if float_list:
        pass
    else:
        float_list = list(np.linspace(0,1,len(rgb_list)))

    cdict = dict()
    for num, col in enumerate(['red', 'green', 'blue']):
        col_list = [[float_list[i], rgb_list[i][num], rgb_list[i][num]] for i in range(len(float_list))]
        cdict[col] = col_list
    cmp = mcolors.LinearSegmentedColormap('my_cmp', segmentdata=cdict, N=256)
    cmp = mcolors.ListedColormap(rgb_list)
    return cmp

def createColorMap(lev, threshold, h, c=[10, 80], l=[10, 80], power=[1, 1], test=False):
    """
    Create a listed colormap from hlc info
    Inputs:
            -lev: levels associated to the colors
            -threshold: list of thresholds of len n
            -h: list of hues of len n+1
            -c=[10, 80]: max and min croma values
            -l=[10, 80]: max and min lightness values
            -power=[1, 1]: power for sampling of croma and lightness
            -test=False: True is test plot wanted
    Outputs: Returns the new colormap
            -cmap: listed color map
            -colors: hex color list
            -norm: norm for plotting
    """
    H = np.repeat(h[0], len(lev))
    for t, hh in zip(threshold, h[1:1+len(threshold)]):
        H[np.where(np.asarray(lev) >= t)] = hh
    C = np.power(np.linspace(0, 1, len(lev), dtype = float), power[0]) * (c[1]-c[0]) + c[0]
    L = l[1] - np.power(np.linspace(0, 1, len(lev), dtype = float), 1) * (l[1]-l[0])
    # Create a HCL color object
    cols = HCL(H, C, L)
    # Load colors
    colors  = cols.colors()
    cmap = get_continuous_cmap(colors)
    norm = mcolors.BoundaryNorm(lev, cmap.N)
    if test:
        testColormap(colors, lev, norm, type='hex')

    return colors, cmap, norm

def makeColorMap(colors, lev, levMin, levMax):
    """
    Extract part of a listed colormap
    Inputs:
            -colors: list of hex colors
            -lev: levels associated to the colors
            -levMin: min of the new levels
            -levMax: max of the new levels
    Outputs: Returns the new colormap
            -newCmap: new color map
            -newColors: new hex color list
            -newLev: new levels list
            -newNorm: new norm for plotting
    """
    indStart = np.where(np.asarray(lev) <= levMin)[0][-1]
    indEnd = np.where(np.asarray(lev) >= levMax)[0][0]
    newLev = lev[indStart:indEnd]
    newColors = colors[indStart:indEnd]
    newCmap = get_continuous_cmap(newColors)
    newNorm = mcolors.BoundaryNorm(newLev, newCmap.N)

    return newCmap, newColors, newLev, newNorm

def testColormap(colors, lev, norm, type='hex'):
    """ Function that plots a given color map for normal vision, greyscale
        and alterated vision
    """
    if type=='hex':
        cmap = get_continuous_cmap(colors)
    elif type=='cmap':
        cmap = colors
        colors = []
        for i in range(cmap.N):
            rgb = cmap(i)[:3] # will return rgba, we take only first 3 so we get rgb
            colors.append(str(mcolors.rgb2hex(rgb)))


    colorsDeut = CVD(colors, "deutan").colors()
    cmapDeut = get_continuous_cmap(colorsDeut)
    colorsProt = CVD(colors, "protan").colors()
    cmapProt = get_continuous_cmap(colorsProt)
    colorsTrit = CVD(colors, "tritan").colors()
    cmapTrit = get_continuous_cmap(colorsTrit)
    colorsDesat = desaturate(colors)
    cmapDesat = get_continuous_cmap(colorsDesat)

    # make some data to plot
    levMin = min(lev)
    levMax = max(lev)

    xx = np.linspace(levMin, levMax, 500, dtype = float)
    yy = np.linspace(levMin, levMax, 500, dtype = float)
    x ,y = np.meshgrid(xx, yy)
    z = np.sqrt(x**2+y**2)

    fig = plt.figure(figsize=(10, 10))

    ax1 = plt.subplot(111)
    im1 = NonUniformImage(ax1, extent=[x.min(), x.max(), y.min(), y.max()], cmap=cmap, norm=norm, origin='lower')
    im1.set_clim(vmin=levMin,vmax=levMax)
    im1.set_data(x[0,:], y[:,0], z)
    ref1 = ax1.images.append(im1)
    cbar = ax1.figure.colorbar(im1, extend='both',ax=ax1)
    ax1.set_xlim(levMin,levMax)
    ax1.set_ylim(levMin,levMax)
    ax1.set_title('Normal vision')

    fig2 = plt.figure(figsize=(10, 10))
    ax2 = plt.subplot(221)
    im2 = NonUniformImage(ax2, extent=[x.min(), x.max(), y.min(), y.max()], cmap=cmapDesat, norm=norm, origin='lower')
    im2.set_clim(vmin=levMin,vmax=levMax)
    im2.set_data(x[0,:], y[:,0], z)
    ref2 = ax2.images.append(im2)
    cbar = ax2.figure.colorbar(im2, ax=ax2)
    ax2.set_xlim(levMin,levMax)
    ax2.set_ylim(levMin,levMax)
    ax2.set_title('Desaturated version')

    ax3 = plt.subplot(222)
    im3 = NonUniformImage(ax3, extent=[x.min(), x.max(), y.min(), y.max()], cmap=cmapDeut, norm=norm, origin='lower')
    im3.set_clim(vmin=levMin,vmax=levMax)
    im3.set_data(x[0,:], y[:,0], z)
    ref3 = ax3.images.append(im3)
    cbar = ax3.figure.colorbar(im3, ax=ax3)
    ax3.set_xlim(levMin,levMax)
    ax3.set_ylim(levMin,levMax)
    ax3.set_title('Deuteranope vision')

    ax4 = plt.subplot(223)
    im4 = NonUniformImage(ax4, extent=[x.min(), x.max(), y.min(), y.max()], cmap=cmapProt, norm=norm, origin='lower')
    im4.set_clim(vmin=levMin,vmax=levMax)
    im4.set_data(x[0,:], y[:,0], z)
    ref4 = ax4.images.append(im4)
    cbar = ax4.figure.colorbar(im4, ax=ax4)
    ax4.set_xlim(levMin,levMax)
    ax4.set_ylim(levMin,levMax)
    ax4.set_title('Proteranope vision')

    ax5 = plt.subplot(224)
    im5 = NonUniformImage(ax5, extent=[x.min(), x.max(), y.min(), y.max()], cmap=cmapTrit, norm=norm, origin='lower')
    im5.set_clim(vmin=levMin,vmax=levMax)
    im5.set_data(x[0,:], y[:,0], z)
    ref5 = ax5.images.append(im5)
    cbar = ax5.figure.colorbar(im5, ax=ax5)
    ax5.set_xlim(levMin,levMax)
    ax5.set_ylim(levMin,levMax)
    ax5.set_title('Triteranope vision')
    plt.show()





if __name__ == "__main__":
    # cmap = sequential_hcl(h = [180, 310], c = [0, 90, 0], l = [95, 30], power = [1.5,1.5])
    # colors = cmap(len(lev))
    #testColormap(colors, type='hex')


    # multi sequential colormap for pressure
    levP  = [0., 0.25, 0.50, 0.75, 1.0, 1.25, 1.50, 1.75, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5,
            5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 14.0, 16.0, 20.0, 25.0, 30.0, 35.0,
            40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 125.0, 150.0, 175.0, 200.0]
    ticksP=[0, 1, 3, 5, 10, 20, 40, 60, 100, 150, 200]
    threshold = [1, 3, 5, 10]
    h = [140, 180, 250, 300, 350]
    colorsP, cmapP, normP = createColorMap(lev, threshold, h, c=[10, 80], l=[10, 80], power=[1, 1], test=True)

    cmap, newColors, newLev, norm = makeColorMap(colorsP, levP, 0.0, 15.0)



    # make up some randomly distributed data
    x, y = np.mgrid[0:10:0.01, 0:10:0.01]
    z = x
    fig = plt.figure(figsize=(10, 10))
    ax1 = plt.subplot(111)
    im1 = NonUniformImage(ax1, extent=[x.min(), x.max(), y.min(), y.max()], cmap=cmap, norm=norm, origin='lower')
    im1.set_clim(vmin=0,vmax=10)
    im1.set_data(x[:,0], y[0,:], z)
    ref0 = ax1.images.append(im1)
    cbar = ax1.figure.colorbar(im1, extend='both',ax=ax1, ticks=[0, 1, 3, 5, 10, 20, 30, 50, 70, 100])
    ax1.set_xlim(0,10)
    ax1.set_ylim(0,10)
    ax1.set_title('Normal vision')
    plt.show()
