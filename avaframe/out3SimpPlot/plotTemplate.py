import sys
import os
import math
import numpy as np
import scipy as sp
import copy
from matplotlib import pyplot as plt
import matplotlib
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
from matplotlib import cm
from matplotlib.image import NonUniformImage
import seaborn as sns





x1 = np.linspace(0,10,21)
x2 = np.linspace(0,10,21)
x2 = x2 + 0.5*np.random.rand(np.shape(x2)[0])

y1 = np.cos(x1)
y2 = np.cos(x1+1)

X, Y = np.meshgrid(x1,x1)

Z = np.cos(X)*np.cos(Y)
Z[0:4,0:4] = np.NaN


figW = 20
figH = 20
figReso = 150
lw = 1
ms = 2
fs = 12


sns.set()
sns.set_style("dark")
cmap = sns.cubehelix_palette(8, start=.5, rot=-.75, as_cmap=True)
cmap.set_under(color='b')
cmap.set_over(color='r')
cmap.set_bad(color='k')
cmap1 = sns.cubehelix_palette(8, as_cmap=True)
cmap1.set_bad(color='k')

cmapdiv = sns.color_palette("RdBu_r")
# cmapdiv.set_bad(color='k')
# sns.set_style("dark", {'xtick.bottom': True,
#  'ytick.left': True})

def linePlot():


    fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(figW, figH), dpi=figReso)
    fig.suptitle('Big Title', fontsize=fs)


    ax[0,0].plot(x1, y1, 'k-o', linewidth=lw, markersize=ms, label='Reference')
    ax[0,0].plot(x1, y2, 'b--', label='Simulation')
    ax[0,0].set_xlabel('x axis label', fontsize=fs)
    ax[0,0].set_ylabel('y axis label', fontsize=fs)
    ax[0,0].set_title('title subplot 1')
    ax[0,0].legend()


    im = ax[0,1].imshow(Z, cmap=cmap, origin='lower', vmin=-0.8, vmax=0.8, extent=[0,10,0,10])
    cbar = fig.colorbar(im, ax=ax[0,1], extend='both')
    cbar.set_label('Label name size 10',size=10)
    ax[0,1].set_xlabel('x axis label', fontsize=fs)
    ax[0,1].set_ylabel('y axis label', fontsize=fs)
    ax[0,1].set_title('title subplot 2\n Black = No data')

    im1 = NonUniformImage(ax[1,0], cmap=cmap)
    im1.set_data(x1, x2, Z)
    im1.set_clim(vmin=-0.8, vmax=0.8)
    ax[1,0].images.append(im1)
    cbar = fig.colorbar(im1, ax=ax[1,0], extend='both')
    cbar.set_label('Label name size 15',size=15)
    ax[1,0].set_xlabel('x axis label', fontsize=fs)
    ax[1,0].set_ylabel('y axis label', fontsize=fs)
    ax[1,0].set_xlim([x1.min(), x1.max()])
    ax[1,0].set_ylim([x2.min(), x2.max()])


    im = ax[1,1].imshow(Z, cmap=cmap1, origin='lower', extent=[0,10,0,10])
    cbar = fig.colorbar(im, ax=ax[1,1])
    cbar.set_label('Label name size 10',size=10)
    ax[1,1].set_xlabel('x axis label', fontsize=fs)
    ax[1,1].set_ylabel('y axis label', fontsize=fs)

    # fig.tight_layout()
    # plt.show()

    sns.set_style("dark")

    fig = plt.figure(figsize=(figW, figH), dpi=figReso)
    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(222)
    ax3 = fig.add_subplot(223)

    im0 = ax1.imshow(Z, cmap=cmap, origin='lower', vmin=-0.8, vmax=0.8, extent=[0,10,0,10])
    cbar = fig.colorbar(im0, ax=ax1, extend='both')
    cbar.set_label('Label name size 10',size=10)
    # ax1.set_xlabel('x axis label', fontsize=fs)
    # ax1.set_ylabel('y axis label', fontsize=fs)
    ax1.set_title('title subplot 2\n Black = No data')


    im1 = ax2.imshow(Z, cmap=cmap1, origin='lower', vmin=-0.8, vmax=0.8, extent=[0,10,0,10])
    cbar = fig.colorbar(im1, ax=ax2, extend='both')
    cbar.set_label('Label name size 10',size=10)
    # ax2.set_xlabel('x axis label', fontsize=fs)
    # ax2.set_ylabel('y axis label', fontsize=fs)
    ax2.set_title('title subplot 2\n Black = No data')


    im2 = NonUniformImage(ax3, cmap=cmap)
    im2.set_data(x1, x2, Z)
    im2.set_clim(vmin=-0.8, vmax=0.8)
    ax3.images.append(im2)
    cbar = fig.colorbar(im2, ax=ax3, extend='both')
    cbar.set_label('Label name size 10',size=10)
    # ax3.set_xlabel('x axis label', fontsize=fs)
    # ax3.set_ylabel('y axis label', fontsize=fs)
    ax3.set_xlim([x1.min(), x1.max()])
    ax3.set_ylim([x2.min(), x2.max()])


    plt.show()


if __name__ == '__main__':
    linePlot()
