import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.colors as mcolors

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
    ''' creates and returns a color map that can be used in heat map figures.
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
    return cmp


def makeColorMap():
    ## Define choosen color palette first
    lev  = [0., 0.01, 0.1, 0.25, 0.50, 0.75, 1.0, 1.25, 1.50, 1.75, 2.0,
                2.5, 3.0, 4.0, 5.0, 7.0, 10.0, 15.0, 20.0, 60.0]
    lev = 20*np.power(np.linspace(0, 1, 40, dtype = float), 1)


    H = np.repeat(180, len(lev))
    # H[np.where(np.asarray(lev) >= 1.)] = 180 # Blueish above 2.0 inches
    H[np.where(np.asarray(lev) >= 5.)] = 250   # Reddish above 5.0 inches
    H[np.where(np.asarray(lev) >= 10.)] = 310   # Reddish above 5.0 inches
    C = np.power(np.linspace(0, 1, len(lev), dtype = float), 1) * 90
    L = 90 - np.power(np.linspace(0, 1, len(lev), dtype = float), 1) * 75
    # Create a HCL color object
    cols = HCL(H, C, L)
    # Load colors
    colors  = cols.colors()
    colors = sequential_hcl(h = [180, 310], c = [0, 90, 0], l = [95, 30], power = [1.5,1.5])(len(lev))

    colorsDeut = CVD(colors, "deutan").colors()
    colorsProt = CVD(colors, "protan").colors()
    colorsTrit = CVD(colors, "tritan").colors()
    colorsDesat = desaturate(colors)

    # # make up some randomly distributed data
    # x, y = np.mgrid[0:20:0.05, 0:20:0.05]
    # z = 20-np.sqrt(x**2+y**2)
    # z[np.where(z < 0)] = 0
    #
    #
    # fig = plt.figure(figsize=(10, 10))
    #
    # ax1 = plt.subplot(111)
    # CS = ax1.contour(x,y,z,len(colors)-1,linewidths=0.5,colors='k')
    # im1 = ax1.contourf(x,y,z,len(colors)-1,colors=colors,
    # vmax=abs(z).max(), vmin=-abs(z).max())
    # cbar = ax1.figure.colorbar(im1, ax=ax1, use_gridspec=True)
    # # cbar.ax.set_ylabel('peak pressure [kPa]')
    # ax1.set_xlim(0,20)
    # ax1.set_ylim(0,20)
    # ax1.set_title('Normal vision')
    #
    # fig2 = plt.figure(figsize=(10, 10))
    # ax2 = plt.subplot(221)
    # CS = ax2.contour(x,y,z,len(colors)-1,linewidths=0.5,colors='k')
    # im2 = ax2.contourf(x,y,z,len(colors)-1,colors=colorsDesat,
    # vmax=abs(z).max(), vmin=-abs(z).max())
    # cbar = ax2.figure.colorbar(im2, ax=ax2, use_gridspec=True)
    # # cbar.ax.set_ylabel('peak pressure [kPa]')
    # ax2.set_xlim(0,20)
    # ax2.set_ylim(0,20)
    # ax2.set_title('Desaturated version')
    #
    # ax3 = plt.subplot(222)
    # CS = ax3.contour(x,y,z,len(colors)-1,linewidths=0.5,colors='k')
    # im3 = ax3.contourf(x,y,z,len(colors)-1,colors=colorsDeut,
    # vmax=abs(z).max(), vmin=-abs(z).max())
    # cbar = ax3.figure.colorbar(im3, ax=ax3, use_gridspec=True)
    # # cbar.ax.set_ylabel('peak pressure [kPa]')
    # ax3.set_xlim(0,20)
    # ax3.set_ylim(0,20)
    # ax3.set_title('Deuteranope vision')
    #
    # ax4 = plt.subplot(223)
    # CS = ax4.contour(x,y,z,len(colors)-1,linewidths=0.5,colors='k')
    # im4 = ax4.contourf(x,y,z,len(colors)-1,colors=colorsProt,
    # vmax=abs(z).max(), vmin=-abs(z).max())
    # cbar = ax4.figure.colorbar(im4, ax=ax4, use_gridspec=True)
    # # cbar.ax.set_ylabel('peak pressure [kPa]')
    # ax4.set_xlim(0,20)
    # ax4.set_ylim(0,20)
    # ax4.set_title('Proteranope vision')
    #
    # ax5 = plt.subplot(224)
    # CS = ax5.contour(x,y,z,len(colors)-1,linewidths=0.5,colors='k')
    # im5 = ax5.contourf(x,y,z,len(colors)-1,colors=colorsTrit,
    # vmax=abs(z).max(), vmin=-abs(z).max())
    # cbar = ax5.figure.colorbar(im5, ax=ax5, use_gridspec=True)
    # # cbar.ax.set_ylabel('peak pressure [kPa]')
    # ax5.set_xlim(0,20)
    # ax5.set_ylim(0,20)
    # ax5.set_title('Triteranope vision')
    # plt.show()
    return colors





if __name__ == "__main__":
    colors = makeColorMap()
    x, y = np.mgrid[0:20:0.05, 0:20:0.05]
    z = 20-np.sqrt(x**2+y**2)
    z[np.where(z < 0)] = 0
    hex_list = colors #['#E2E2E2', '#DADEDD', '#D1D9D8', '#C8D5D3', '#BFD0CE', '#B6CCC9', '#ADC8C4', '#A4C3BF', '#9BBFBA', '#92BAB5', '#A1AEC7', '#9BA9C4', '#94A3C1', '#8E9EBE', '#8799BC', '#8094B9', '#7A8FB6', '#738BB3', '#6C86B1', '#6581AE', '#A2669F', '#9F609C', '#9C5999', '#995396', '#964C93', '#934590', '#913D8D', '#8E358B', '#8C2B88', '#8A2086', '#880F84', '#860082', '#850081', '#840080', '#84007F', '#840080', '#870082', '#8B0086', '#950090', '#AA00A3']
    fig, ax = plt.subplots(1,1)
    im = ax.imshow(z, cmap=get_continuous_cmap(hex_list))
    fig.colorbar(im)
    ax.yaxis.set_major_locator(plt.NullLocator()) # remove y axis ticks
    ax.xaxis.set_major_locator(plt.NullLocator()) # remove x axis ticks
    plt.show()
