import numpy as np
import matplotlib.colors as mcolors
from matplotlib import cm


def hex_to_rgb(value):
    '''
    Converts hex to rgb colours
    value: string of 6 characters representing a hex colour.
    Returns: list length 3 of RGB values'''
    value = value.strip("#")  # removes hash symbol if present
    lv = len(value)
    return tuple(int(value[i:i + lv // 3], 16) for i in range(0, lv, lv // 3))


def rgb_to_dec(value):
    '''
    Converts rgb to decimal colours (i.e. divides each value by 256)
    value: list (length 3) of RGB values
    Returns: list (length 3) of decimal values'''
    return [v/256 for v in value]


def get_continuous_cmap(hex_list, continuous=False, float_list=None):
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
        float_list = list(np.linspace(0, 1, len(rgb_list)))
    if continuous:
        cdict = dict()
        for num, col in enumerate(['red', 'green', 'blue']):
            col_list = [[float_list[i], rgb_list[i][num], rgb_list[i][num]]
                        for i in range(len(float_list))]
            cdict[col] = col_list

        cmp = mcolors.LinearSegmentedColormap('my_cmp', segmentdata=cdict, N=256)
    else:
        cmp = mcolors.ListedColormap(rgb_list)
    return cmp


def makeColorMap(colormap, levMax, continuous=False):
    """
    Extract part of a listed colormap
    Inputs:
            -colors: list of hex colors
            -lev: levels associated to the colors
            -levMax: max of the new levels
    Outputs: Returns the new colormap
            -newCmap: new color map
            -newColors: new hex color list
            -newLev: new levels list
            -newNorm: new norm for plotting
    """
    colors = colormap['colors']
    cmap = colormap['cmap']
    lev = colormap['lev']
    ticks = colormap['ticks']
    if continuous:
        newNorm = None
        newLev = None
        newColors = None
        newCmap = cmap
        ticks = None

    else:
        # indStart = np.where(np.asarray(lev) <= levMin)[0][-1]
        indEnd = np.where(np.asarray(lev) > levMax)[0][0]
        newLev = lev[:indEnd+1]
        newColors = colors[:indEnd]
        newCmap = get_continuous_cmap(newColors, continuous=False)
        newNorm = mcolors.BoundaryNorm(newLev, newCmap.N)

    return newCmap, newColors, newLev, newNorm, ticks
