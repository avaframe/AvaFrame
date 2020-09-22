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


def makeColorMap(plotFlag=False):
    ## Define choosen color palette first
    lev  = [0., 0.01, 0.1, 0.25, 0.50, 0.75, 1.0, 1.25, 1.50, 1.75, 2.0,
                2.5, 3.0, 4.0, 5.0, 7.0, 10.0, 15.0, 20.0, 60.0]
    lev = np.power(np.linspace(0, 1, 20, dtype = float), 1)


    H = np.repeat(180, len(lev))
    H[np.where(np.asarray(lev) >= 0.2)] = 250   # Reddish above 5.0 inches
    H[np.where(np.asarray(lev) >= 0.5)] = 310   # Reddish above 5.0 inches
    C = np.power(np.linspace(0, 1, len(lev), dtype = float), 1) * 90
    L = 90 - np.power(np.linspace(0, 1, len(lev), dtype = float), 1) * 75
    # Create a HCL color object
    cols = HCL(H, C, L)
    # Load colors
    colors  = cols.colors()
    cmap = get_continuous_cmap(colors)

    return cmap, colors

def testColormap(colors, type='hex'):
    if type=='hex':
        pass
    elif type=='cmap':
        cmap = colors
        colors = []
        for i in range(cmap.N):
            rgb = cmap(i)[:3] # will return rgba, we take only first 3 so we get rgb
            colors.append(str(mcolors.rgb2hex(rgb)))


    colorsDeut = CVD(colors, "deutan").colors()
    colorsProt = CVD(colors, "protan").colors()
    colorsTrit = CVD(colors, "tritan").colors()
    colorsDesat = desaturate(colors)

    # make up some randomly distributed data
    x, y = np.mgrid[0:1:0.01, 0:1:0.01]
    z = x #np.sqrt(x**2+y**2)
    z[np.where(z < 0)] = 0


    fig = plt.figure(figsize=(10, 10))

    ax1 = plt.subplot(111)
    im1 = ax1.contourf(x,y,z,len(colors)-1,colors=colors,vmax=abs(z).max(), vmin=-abs(z).max())
    # im1 = NonUniformImage(ax1, extent=[x.min(), x.max(), y.min(), y.max()], cmap=cmap)
    # im.set_interpolation('bilinear')
    # im1.set_clim(vmin=0.000000001)
    # im1.set_data(x[:,0], y[0,:], z)
    # ref0 = ax1.images.append(im1)
    cbar = ax1.figure.colorbar(im1, ax=ax1, use_gridspec=True)
    # cbar.ax.set_ylabel('peak pressure [kPa]')
    ax1.set_xlim(0,1)
    ax1.set_ylim(0,1)
    ax1.set_title('Normal vision')

    fig2 = plt.figure(figsize=(10, 10))
    ax2 = plt.subplot(221)
    # CS = ax2.contour(x,y,z,len(colors)-1,linewidths=0.5,colors='k')
    im2 = ax2.contourf(x,y,z,len(colors)-1,colors=colorsDesat,
    vmax=abs(z).max(), vmin=-abs(z).max())
    cbar = ax2.figure.colorbar(im2, ax=ax2, use_gridspec=True)
    # cbar.ax.set_ylabel('peak pressure [kPa]')
    ax2.set_xlim(0,1)
    ax2.set_ylim(0,1)
    ax2.set_title('Desaturated version')

    ax3 = plt.subplot(222)
    # CS = ax3.contour(x,y,z,len(colors)-1,linewidths=0.5,colors='k')
    im3 = ax3.contourf(x,y,z,len(colors)-1,colors=colorsDeut,
    vmax=abs(z).max(), vmin=-abs(z).max())
    cbar = ax3.figure.colorbar(im3, ax=ax3, use_gridspec=True)
    # cbar.ax.set_ylabel('peak pressure [kPa]')
    ax3.set_xlim(0,1)
    ax3.set_ylim(0,1)
    ax3.set_title('Deuteranope vision')

    ax4 = plt.subplot(223)
    # CS = ax4.contour(x,y,z,len(colors)-1,linewidths=0.5,colors='k')
    im4 = ax4.contourf(x,y,z,len(colors)-1,colors=colorsProt,
    vmax=abs(z).max(), vmin=-abs(z).max())
    cbar = ax4.figure.colorbar(im4, ax=ax4, use_gridspec=True)
    # cbar.ax.set_ylabel('peak pressure [kPa]')
    ax4.set_xlim(0,1)
    ax4.set_ylim(0,1)
    ax4.set_title('Proteranope vision')

    ax5 = plt.subplot(224)
    # CS = ax5.contour(x,y,z,len(colors)-1,linewidths=0.5,colors='k')
    im5 = ax5.contourf(x,y,z,len(colors)-1,colors=colorsTrit,
    vmax=abs(z).max(), vmin=-abs(z).max())
    cbar = ax5.figure.colorbar(im5, ax=ax5, use_gridspec=True)
    # cbar.ax.set_ylabel('peak pressure [kPa]')
    ax5.set_xlim(0,1)
    ax5.set_ylim(0,1)
    ax5.set_title('Triteranope vision')
    plt.show()





if __name__ == "__main__":
    pow = 2
    lev = np.linspace(0, 1, 30, dtype = float)
    # lev[np.where(lev <= 1)] = np.power(lev[np.where(lev <= 1)], pow)
    # lev[np.where(lev > 1)] = pow*lev[np.where(lev > 1)]-1
    # lev = 10*lev
    # lev  = [0., 0.25, 0.50, 0.75, 1.0, 1.50, 2.0,
    #         2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0,
    #         12.0, 14.0, 16.0, 20.0, 25.0, 30.0, 35.0, 40.0, 50.0]
    print(lev)

    # fig = plt.figure(figsize=(10, 10))
    #
    # ax1 = plt.subplot(111)
    # ax1.plot(np.linspace(0, 1, len(lev), dtype = float), lev)
    # plt.show()


    # cmap = sequential_hcl(h = [180, 310], c = [0, 90, 0], l = [95, 30], power = [1.5,1.5])
    # colors = cmap(len(lev))
    #
    # testColormap(colors, type='hex')
    #
    # cmap = sequential_hcl(h = [180, 310], c = [0, 90, 0], l = [95, 30], power = [1.5,1.5])
    # colors = cmap(len(lev))
    #
    # testColormap(colors, type='hex')

    H = np.repeat(140, len(lev))
    H[np.where(np.asarray(lev) >= 0.2)] = 180   # Reddish above 5.0 inches
    H[np.where(np.asarray(lev) >= 0.4)] = 250   # Reddish above 5.0 inches
    H[np.where(np.asarray(lev) >= 0.6)] = 300
    H[np.where(np.asarray(lev) >= 0.8)] = 350
    C = np.power(np.linspace(0, 1, len(lev), dtype = float), 1) * 70 + 10
    L = 80 - np.power(np.linspace(0, 1, len(lev), dtype = float), 1) * 70
    # Create a HCL color object
    cols = HCL(H, C, L)
    # Load colors
    colors  = cols.colors()
    cmap = get_continuous_cmap(colors)

    # testColormap(colors, type='hex')


    # make up some randomly distributed data
    x, y = np.mgrid[0:10:0.01, 0:10:0.01]
    z = x #np.sqrt(x**2+y**2)
    z[np.where(z < 0)] = 0


    fig = plt.figure(figsize=(10, 10))

    ax1 = plt.subplot(111)
    im1 = ax1.contourf(x,y,z,len(colors)-1,colors=colors,vmax=abs(z).max(), vmin=-abs(z).max())
    # im1 = NonUniformImage(ax1, extent=[x.min(), x.max(), y.min(), y.max()], cmap=cmap)
    # im.set_interpolation('bilinear')
    # im1.set_clim(vmin=0.000000001)
    # im1.set_data(x[:,0], y[0,:], z)
    # ref0 = ax1.images.append(im1)
    cbar = ax1.figure.colorbar(im1, ax=ax1, use_gridspec=True)
    # cbar.ax.set_ylabel('peak pressure [kPa]')
    ax1.set_xlim(0,10)
    ax1.set_ylim(0,10)
    ax1.set_title('Normal vision')
    plt.show()
