"""
Auxiliary functions for plotting particle results

Created on Tue Jan 09 2023
@author: Höller JF
"""


# Load modules
import numpy as np                                  # Vektormathematik
import matplotlib as mpl                            # Get / use colormaps
import matplotlib.pyplot as plt                     # Plot created plots
from matplotlib.tri import Triangulation, TriAnalyzer, UniformTriRefiner, LinearTriInterpolator
from matplotlib.offsetbox import AnchoredText
import pathlib                                      # Work with smart paths
import math                                         # Round up to nearest 10
import pandas as pd
import scipy as sp                                  # Interpolation and datafiltering


# Temporary modules


# Local imports
import PlotParticles_auxiliary_MF_JFH as ppam       # Auxiliary functions for plotting


####   ----------------------------------------------------------------------------   ###
#                                  Auxiliary Plot Functions                             #
####   ----------------------------------------------------------------------------   ###


def GenPlotSettings(figcolor):
    params = {#'backend': 'TKAgg',
            'axes.edgecolor': figcolor,
            'axes.labelcolor': figcolor,
            'axes.grid' : True,
            'grid.color' : 'gray',
            'grid.linestyle': '--',
            'text.color' : figcolor,              # Color of title, subtitle
            'xtick.color' : figcolor,
            'ytick.color' : figcolor,
            #'font.size': 8,
            'legend.fontsize': 6,
            #'ytick.minor.visible': True,
            #'axes.prop_cycle' : plt.cycler("color", plt.cm.tab20.colors),
            }
    plt.rcParams.update(params)
    mpl.use('TKAgg')
    # Check rc.params: plt.rcParams.keys()



def PlotForLatex():
    params = {'figure.figsize': [5.34213, 2.3],     # [8,4],
            'font.size': 11.0,
            'font.family': 'serif',
            'font.serif': [],                       # 'Computer Modern Roman'
            'axes.titlesize': 'medium',
            'figure.titlesize': 'medium',
            'text.usetex': True,
            'pgf.texsystem': 'lualatex',
            'text.latex.preamble': "\n".join([      # Preamble for Latex used in Backend
                r'\usepackage{amsmath}',            # Extensions for the mathematical formula set -> Liefert Fehlermeldung
                r'\usepackage{siunitx}',            # 
                r'\sisetup{locale = DE,per-mode = symbol,group-separator=.,group-minimum-digits=4,group-digits=integer}',
                r'\usepackage{lmodern}',
                ]),
            'pgf.rcfonts': True,                    # don't setup fonts from rc parameters
            'pgf.preamble': "\n".join([             # Preamble for .pgf Export
                #r'\usepackage{fontspec}',
                r'\usepackage{amsmath}',            # Extensions for the mathematical formula set
                r'\usepackage{siunitx}',
                r'\sisetup{locale = DE,per-mode = symbol,group-separator=.,group-minimum-digits=4,group-digits=integer}',
                r'\usepackage{unicode-math}',       # unicode math setup
                r'\usepackage{lmodern}',
                #r'\setmainfont{sans}',
                ]),
                # r"\usepackage[utf8x]{inputenc}",    # use utf8 fonts 
                # r"\usepackage[T1]{fontenc}",        # plots will be generated
            }
    plt.rcParams.update(params)
    #mpl.use('pgf')
    # https://matplotlib.org/stable/tutorials/text/pgf.html


def PropertyLabel(property):
    """ Create Y-label for specific particle properties
    
    Parameters
    ----------
    property: str
        Property name

    Returns
    ----------
        Proper property name for labeling an axis
    """
    match property:
        case 't':
            property = r'$t$ [s]'
            #  property = r'Time $t$ [s]'
        case 'umag':
            property = r'$v$ [\si{\m\per\s}]'     # Velocity magnitude 
            #property = r'$v_{mag}$ [m/s]'     # Velocity magnitude
        case 'umag_max':
            property = r'$v_\text{pfv}$ [\si{\m\per\s}]'     # Velocity magnitude 
            # property = r'peak flow velocity $v_\text{pfv}$ [\si{\m\per\s}]'     # Velocity magnitude 
        case 's':
            property = r'Horizontal distance $s_\text{xy}$ [m]'
        case 'sAimec':
            property = r'$s_\text{XY}$ [m]'
            # property = r'Travel length on the avalanche path $s_\text{Aimec}$ [m]'
        case 'sAimec_max':
            property = r'max $s_\text{XY}$ [m]'
        case 'lAimec':
            property = r'$l_\text{XY}$ [m]'
        case 'travelAngle':
            property = r'Travel angle [°]'
        case 'h':
            property = r'Flow depth $h$ [m]'
        case 'dAimec':
            property = r'$d_\text{CoM}$ [m]'
            # property = r'Distance between particle and CoM $d_\text{CoM}$ [m]'
        case 'dsCoM':
            # property = r'Path distance between particle and CoM $d_\text{s,CoM}$ [m]'
            property = r'rel. distance $d_\text{s,CoM}$ [m]'
        case 'dsCoM_max':
            property = r'max rel. distance $d_\text{s,CoM}$ [m]'
        case 'sParticle':
            property = r'$s^\text{particle}$ [m]'
            # property = r'Particle travel length $s^\text{particle}$ [m]'
        case other:
            print('Lable for the particle property %s is not yet implemented!' % (property))
    return property



def LablePlot(ax,xlabel,ylabel,*args):
    """ Label matplotlib plot
    
    Parameters
    ----------
    ax: ax
        Axis Handle
    xlabel: str
    ylable: str
    *args:
        title: str
    """
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if args:
        ax.set_title(args[0])



def PlotSettings(ax,aspect='square',factor=1.1):
    """ Define standard plot settings
    
    Parameters
    ----------
    ax: ax
        Axis Handle
    aspect: str
        Choose axis layout ('square','centre','optimise')
    factor: float
        Choose distance to the edge 
    """
    # ax.axis('square')
    # ax.set_aspect('equal')
    # ax.set_aspect(1)
    match aspect:
        case 'square':
            ax.axis('square')

        case 'centre':
            ax.set_aspect('equal')

            # Check whick delta is larger
            YLim = ax.get_ylim()
            deltay = YLim[1]-YLim[0]
            XLim = ax.get_xlim()
            deltax = XLim[1]-XLim[0]

            if deltax>deltay:
                delta = deltax*factor
            else:
                delta = deltay*factor
            
            # Calculate centre point
            centrex = XLim[0]+deltax/2
            centrey = YLim[0]+deltay/2
            ax.set_xlim([centrex-delta/2,centrex+delta/2]) # add 5 %
            ax.set_ylim([centrey-delta/2,centrey+delta/2]) # add 5 %

        case 'optimise':
            ax.set_aspect('equal')

            YLim = ax.get_ylim()
            deltay = (YLim[1]-YLim[0])

            XLim = ax.get_xlim()
            deltax = (XLim[1]-XLim[0])

            # Calculate centre point
            centrex = XLim[0]+deltax/2
            centrey = YLim[0]+deltay/2

            ax.set_xlim([centrex-deltax*factor/2,centrex+deltax*factor/2]) # add 5 %
            ax.set_ylim([centrey-deltay*factor/2,centrey+deltay*factor/2]) # add 5 %
        case _:
            print('Aspect setting not programmed!')
            return



def ScatterPlot(ax,particleList,property,t,max=True,cmap='Spectral_r',triangulation=False,**kwargs):
    """ Create scatter plot for the particles in a given time step in the sl-coordinate system
    
    Parameters
    ----------
    ax: ax
        Axis Handle
    particleList: list of dict
        list of dictionaries
    property: str
        Property name
    t: int
        number of time step
    max: bool
        Activate if the particles should be colored by the max value in every time steps
    cmap: str
        name of a color map (matplotlib)
    **kwargs:
        vmax: float
            Enter vmax, for consistant colorsheme in multiple plots
        vmin: float
            Enter vmin, for consistant colorsheme in multiple plots
        ARcoords: list of tuples
            list of coordinates for calculated the corresponing z-value uon the new mesh
        ARprop: list of floats
    """
    x = np.array(particleList[t]['lAimec'])
    y = np.array(particleList[t]['sAimec'])

    if max:
        # Find max value of every particle for coloring the scatter plots
        z = ppam.MaxList(particleList,property)
    else:
        z = particleList[t][property]

    if not triangulation:     # Standard Scatter Plot (ggf. ergänzen: particleList[0]['nPart']<500:)
        if 'vmax' in kwargs:
            scat = ax.scatter(x,y, c=z, s=4, cmap=cmap,vmax=kwargs['vmax'])
        else:
            scat = ax.scatter(x,y, c=z, s=4, cmap=cmap)
    else:           # Use Contour Plot, when to many points in one plot
        if not kwargs is None:
            scat = plot_triangulation(ax,x,y,z,cmap,**kwargs)
        else:
            scat = plot_triangulation(ax,x,y,z,cmap)

    return scat



def plot_triangulation(ax,x,y,z,cmap,refine=True,**kwargs):
    """ Create Contourplot for a given scatterplot
    
    Parameters
    ----------
    ax: ax
        Axis Handle
    x,y,z: np.array
        array of point coordinates
    cmap: str
        name of a color map (matplotlib)
    refine: boolean
        Refine Mesh by interpolation
        Decide if the Data should be interpolated ~can lead to increased max/min results
    **kwargs:
        vmax: float
            Enter vmax, for consistant colorsheme in multiple plots
        vmin: float
            Enter vmin, for consistant colorsheme in multiple plots
        ARcoords: list of tuples
            list of coordinates for calculated the corresponing z-value uon the new mesh
        ARprop: list of floats
    """
    # ---------------------------------------------------------------------------- #
    # Settings
    # ---------------------------------------------------------------------------- #
    # Quantile limit for the distance or maximum (rel.) distance between the points
    # for calculating the concave hull ('quantile', 'relative, 'absolute')
    alpha = 0.98
    calc = 'quantile'

    # Number of recursive subdivisions of the initial mesh for smooth plots.
    # Values >3 might result in a very high number of triangles for the refine
    # mesh: new triangles numbering = (4**subdiv)*ntri
    subdiv = 3

    # Minimum circle ratio - border triangles with circle ratio below this will be
    # masked if they touch a border. Suggested value 0.01; use -1 to keep all
    # triangles.
    min_circle_ratio = .01

    # Round points to a grid (here only y-coordinates)
    roundy = True
    raster_size_y = 1

    # User options for plots
    plot_refi_tri = False       # plot of refined triangulation
    plot_tri = False            # plot of base triangulation
    plot_masked_tri = False     # plot of excessively flat excluded triangles
    plot_scatter = False        # plot original and new points as scatter plot

    # Graphical options for tricontouring
    nlevel = 41 #21 # 16

    # ---------------------------------------------------------------------------- #
    # Step 1: Load Dataframe
    # ---------------------------------------------------------------------------- #
    df_all = pd.DataFrame({'x':x, 'y':y, 'z':z})
    # define a consistent colormap
    if 'vmax' not in kwargs:
        vmax = z.max()
    else:
        vmax = kwargs['vmax']
    if 'vmin' not in kwargs:
        vmin = z.min()
    else:
        vmin = kwargs['vmin']
    levels = np.linspace(vmin, vmax, nlevel)


    # ---------------------------------------------------------------------------- #
    # Step 2: Removing quasi-double points before meshing by rounding the points to a grid
    # ---------------------------------------------------------------------------- #
    # Rounding not optimal -> Better interpolate
    if roundy:
        df_all['y'] = ppam.baseround(df_all['y'],raster_size_y)
    

    # ---------------------------------------------------------------------------- #
    # Step 3: Remove duplicate points and calculate mean z0
    # ---------------------------------------------------------------------------- #
    dfmask = df_all.duplicated(subset=['x', 'y'],keep=False)
    df_dupl = df_all[dfmask]

    index_list = df_dupl.groupby(list(df_dupl)[0:2]).apply(lambda i: tuple(i.index)).tolist()

    for i in index_list:
        i_list = list(i)
        middle = sum(df_all.iloc[i_list]['z'])/len(i_list)
        df_all.loc[i_list,'z'] = middle

    df_red = df_all.drop_duplicates()


    # ---------------------------------------------------------------------------- #
    # Step 4: Meshing with Delaunay triangulation
    # ---------------------------------------------------------------------------- #
    x_tri = np.array(df_red['x'])
    y_tri = np.array(df_red['y'])
    z_tri = np.array(df_red['z'])
    tri = Triangulation(x_tri, y_tri)


    # ---------------------------------------------------------------------------- #
    # Step 5: Create Concav hull, input max edge length
    # ---------------------------------------------------------------------------- #
    apply_mask(x_tri,y_tri,tri,alpha,calc)   # Löscht keine Dreiecke, stellt nur deren Sichtbarkeit aus (mask[i]=True)
    #, alpha=20

    # ---------------------------------------------------------------------------- #
    # Step 6: masking badly shaped triangles at the border of the triangular mesh
    # ---------------------------------------------------------------------------- #
    mask = TriAnalyzer(tri).get_flat_tri_mask(min_circle_ratio)     # Eliminate excessively flat border triangles from the triangulation
    tri.set_mask(mask)


    # ---------------------------------------------------------------------------- #
    # Step 7: refining the data
    # ---------------------------------------------------------------------------- #
    if refine:
        refiner = UniformTriRefiner(tri)
        tri_refi, z_refi = refiner.refine_field(z_tri, subdiv=subdiv)
    else:
        tri_refi = tri
        z_refi = z_tri


    # ---------------------------------------------------------------------------- #
    # Step 8: Plotting
    # ---------------------------------------------------------------------------- #
    # 1) plot of the refined (computed) data contours:
    plot = ax.tricontourf(tri_refi, z_refi, levels=levels, cmap=cmap)  #, vmax=vmax
    
    # 2) plot of the fine mesh on which interpolation was done:
    if plot_refi_tri:
        ax.triplot(tri_refi, color='0.97',linewidth=0.1)
    
    # 3) plot of the initial 'coarse' mesh:
    if plot_tri:
        ax.triplot(tri, color='0.7',linewidth=0.5)
    
    # 4) plot of the unvalidated triangles from naive Delaunay Triangulation:
    if plot_masked_tri:
        # for the demo: loading the 'flat' triangles for plot
        flat_tri = Triangulation(x_tri, y_tri)
        flat_tri.set_mask(~mask)
        ax.triplot(flat_tri, color='red',linewidth=0.5)

    # 5) plot the used points
    if plot_scatter:
        # Original points
        # ax.scatter(x, y, c=z, s=0.5, cmap=cmap, vmax=vmax)
        ax.scatter(x, y, s=0.5, color='k')
        # Masked points
        #ax.scatter(x_tri, y_tri, s=1, color='k')

    # Analyse Mesh -> Get z-value of point on mesh
    if 'ARcoords' in kwargs:
        newz=[]
        deltaprop=[]
        for i,AR in enumerate(kwargs['ARcoords']):
            newz.append(CalcZValueOfTri(AR,tri_refi,z_refi))
            deltaprop.append((newz[i]-kwargs['ARprop'][i])/kwargs['ARprop'][i])
        deltaprop


    return plot



def CalcZValueOfTri(query_point,tri,z_tri):
    """ Calculate the z-Value for a given point on a Triangulated Net
    
    Parameters:
    -----------
    x,y: float
        x and y coordinates of the point
    
    tri:
        triangulation net
    z_tri:
        corresponding z-values for the triangulates

    Returns:
        interpolated_z_value
    """
    #query_point = np.array([x,y])
    interpolator = LinearTriInterpolator(tri,z_tri)
    interpolated_z_value = interpolator(*query_point)
    print("Der Z-Wert für den Punkt ", query_point, "ist ", interpolated_z_value)

    return float(interpolated_z_value)



def CalcHist2d(DataXflat,DataYflat,nBins):
    """ Calculate 2D Histogram

    Parameters
    ----------
    DataXflat, DataYflat: lst of arrays
        Flattened list of arrays
    nBins: int or [int int]
        The bin specification

    Returns
    -------
    X,Y: numpy.ndarray
        Mesh
    H: numpy.ndarray
        2D field of counts per section
    """
    # Create 2d histogram with numpy and plot with matplotlib (faster calculation speed)
    H, xedges, yedges = np.histogram2d(DataXflat,DataYflat,bins=nBins)
    H = H.T
    X, Y = np.meshgrid(xedges, yedges)
    
    return X,Y,H



def PlotHist2D(ax,X,Y,H,fig,colorbar=True,**kwargs):
    """ Plot 2D Histogram

    Parameters
    ----------
    ax: ax
        Axis handle
    X,Y: numpy.ndarray
        Mesh
    H: numpy.ndarray
        2D field of counts per section
    fig: matplotlib fig
        Figure Handle
    colorbar: boolean
        Decide if a colorbar should be created
    **kwargs:
        vmax: float
    """
    #hist = axs.pcolormesh(X,Y,H,cmap = cmapSpectral,vmin=1)
    if colorbar:
        if 'vmax' in kwargs:
            hist = ax.pcolormesh(X,Y,H,cmap = Colormap(),norm = mpl.colors.LogNorm(vmax=kwargs['vmax']))       # vmax=30000
        else:
            hist = ax.pcolormesh(X,Y,H,cmap = Colormap(),norm = mpl.colors.LogNorm())
        # Plot Settings
        cbar = fig.colorbar(hist,ax=ax)
        cbar.set_label('Anzahl Partikel')       # Counts
    else:
        if 'vmax' in kwargs:
            ax.pcolormesh(X,Y,H,cmap = Colormap(),norm = mpl.colors.LogNorm(vmax=kwargs['vmax']))
        else:
            ax.pcolormesh(X,Y,H,cmap = Colormap(),norm = mpl.colors.LogNorm())



def Colormap(color='bone_r'):
    """ Define Colormap

    Parameters
    ----------
    color: str
        name of a colormap of Matplotlib
    """
    #cmapModified = copy.copy(plt.cm.bone_r)
    #cmapModified.set_under('w')

    cmapModified = mpl.colormaps[color]

    return cmapModified



def plotStats(ax,x,y,q=[0.8, 0.5, 0.2],label=True):
    """ Calculation of median and quantile for each time step
    
    Parameters
    ----------
    ax: ax
        Axis handle
    x,y:
        Data
    q: list
        List of quantiles
    label: bool
        de-/activate labels
    """
    sortedX,sortedY = ppam.SortLists(x,y)
    
    Quants = ppam.Quantiles(sortedY,q)
    # label = [(str(round(q[0]*100)) + r'\% Quantile'),
    #         r'Median',
    #         (str(round(q[-1]*100)) + r'\% Quantile')]

    # cmap = mpl.cm.get_cmap('tab20c')
    # ax.set_prop_cycle(color=[cmap(6/20),cmap(7/20)])

    # Calculate Max Y-Value and corresponding X-Value
    # median=ppam.Quantiles(sortedY,0.5)
    # max(median)
    # sortedX[median.index(max(median))]


    color='darkorange'
    style = [':','-','--']

    if label:
        label = [r'Minimalwerte', r'Median', r'Maximalwerte']
        ax.plot(sortedX,Quants,color=color,linewidth=0.8,label=label)
    else:
        ax.plot(sortedX,Quants,color=color,linewidth=0.8)
    for i,j in enumerate(ax.lines):
        j.set_linestyle(style[i])
    # ax.legend()



def plotAvaNodeStats(ax,AvaNodes,xresult,yresult,q=[0,0.5,1],stats=True):
    """ Plot and calculate AvaNode results (Min, Median, Max)

    Parameters
    ----------
    ax: ax
        Axis handle
    AvaNodes: list
        list of AvaNode dictionaries
    """
    color='darkorange'
    if not stats:
        styleAR = ['-','-.','--',':','-','-.','--',':']
    ARX = []
    ARY = []
    steps = 1000
    # get max length of AvaNode results
    max_x = max([AR[-1][xresult] for AR in AvaNodes])
    xi = np.linspace(0,max_x,steps)

    for i,AR in enumerate(AvaNodes):
        x = [i[xresult] for i in AR]
        y = [i[yresult] for i in AR]

        # Add 0 after the last time step up to the max time
        if not max(x) == max_x:
            x.extend([x[-1]+0.1,max_x])      #### ADD: Interpolate results to 0!
            y.extend([0,0])

        # Interpolate AR results to uniform time steps
        yi = np.interp(xi, x, y)              # linear interpolation

        # Add all results to one list for median calculation
        ARX.extend(xi)
        ARY.extend(yi)

        if not stats:
            ax.plot(xi,yi,linestyle=styleAR[i],linewidth=0.8,label='AvaNode '+str(i))
    

    if stats:
        # Median und Umhüllende Berechnen
        sortedX,sortedY = ppam.SortLists(ARX,ARY)
        Quants = ppam.Quantiles(sortedY,q)
        # ax.plot(sortedX,Quants,'b--',linewidth=0.5)

        # Line smoothing with a Savitzky-Golay filter
        label = [r'Minimalwerte', r'Median', r'Maximalwerte']
        style = [':','-','--']
        y_smooth = []
        for i in range(len(Quants[0])):
            temp_Quants = [j[i] for j in Quants]
            y_smooth.append(sp.signal.savgol_filter(temp_Quants, window_length=20, polyorder=3, mode="nearest"))
            
            ax.plot(sortedX,y_smooth[i],color=color,linestyle=style[i],linewidth=0.8,label=label[i])
        
        ax.fill_between(sortedX,y_smooth[0],y_smooth[-1],color=color,edgecolor="none",alpha=0.3,label='AvaNode-Ergebnisse')

    YLim = ax.get_ylim()
    ax.set_ylim([0,YLim[-1]])



def savePlot(fig,name,avaDir,comModule="com1DFA",LatexExport=False):
    """ Create consistent plot styling

    Parameters
    ----------
    fig: matplotlib fig
        Figure Handle
    name: str
        Name of the figure
    avaDir: str
        path to avalanche directory
    comModule: str
        module that computed the particles
    """
    if LatexExport:
        FigTypes = ['.pdf','.pgf','.eps']
        #FigTypes = ['.pdf','.eps']
    else:
        # FigTypes = ['.pdf','.png']
        FigTypes = ['.png']

    for FigType in FigTypes:
        figPath = pathlib.Path(avaDir,'Outputs',comModule,'particleTracing',name).with_suffix(FigType)
        pathlib.Path(figPath).parent.mkdir(parents=True, exist_ok=True)
        # log.info('Saved Figure to: %s' % figPath)
        if FigType == '.png':
            fig.savefig(figPath, dpi=300)
        else:
            fig.savefig(figPath)



def PlotText(ax,text,loc='upper left'):
    """ Add a text patch to the plot
    
    Parameters
    ----------
    ax: ax
        Axis handle
    text: str
        Text to plot
    loc: str
        Location of the patch
    """
    if not isinstance(ax,list) and not isinstance(ax,np.ndarray):
        ax = [ax]
    
    if text=='avaNordketteVU' or text=='avaNordkette':
        size=7 #7
    else:
        size=5

    for axis in ax:
        # text = text+'\phantom{test}'
        at = AnchoredText(text, prop=dict(size=size), frameon=True, loc=loc)
        at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
        # at.patch.set_boxstyle("round",pad=0.3,rounding_size=0.2)

        # bbox = at.get_bbox_patch()
        # bbox.set_boxstyle("round", pad=(0, 0.5, 0, 0), rounding_size=0.2)

        # at.set_bbox_to_anchor((1.1, 0.5))
        axis.add_artist(at)



def apply_mask(x,y,triang, alpha=0.9, calc='quantile'):
    """ Creates concav hull around the points
    Mask off unwanted triangles -> according to the used alpha

    Parameters
    ----------
    x: np.array
        array of x coordinates
    y: np.array
        array of y coordinates
    triang: 'Triangulation' object (matplotlib)
        triangulation object
    alpha: int
        max limit for the distance between the points for concave hull
        mask every triangle above alpha
    calc: str
        Interpretation methode for the alpha value
        'quantile', 'relative', 'absolute'
    """
    # Mask triangles with sidelength bigger some alpha
    # (nicht als Distanz, sonder relativ)
    triangles = triang.triangles
    
    # Mask off unwanted triangles.
    xtri = x[triangles] - np.roll(x[triangles], 1, axis=1)
    ytri = y[triangles] - np.roll(y[triangles], 1, axis=1)
    maxi = np.max(np.sqrt(xtri**2 + ytri**2), axis=1)
    match calc:
        case 'quantile':
            limit = np.quantile(maxi,alpha)     # Calculate max. absolute distance between two points
        case 'relative':
            maxi = maxi/maxi.max()
            limit = alpha
        case 'absolute':
            limit = alpha

    # Apply masking
    # -> If (rel.) distance between two points is larger than alpha will be masked (not visualised)
    # The larger alpha, the less masked triangles
    triang.set_mask(maxi > limit)

    return triang



def plotNodeDict(ax,PartDict,prop,i,**kwargs):
    """ Plot results of a list of particle dictionaries

    Parameters
    ----------
    ax: ax
        axis handle
    PartDict: list of dict
    prop: list of str
        ['propx','propy']
    i: int
        particle number for labeling
    **kwargs:
        color: matplotlib color
            linecolor
    """
    AvaNodeNames = ['1.1','2.1','2.3','3.1','4.10','5.7','5.9','5.10']
    styles = ['-',':',(0, (1, 10)),'--',(5, (10, 3)),'-.',(0, (3, 5, 1, 5, 1, 5)),(0, (3, 1, 1, 1, 1, 1))]

    x = [step[prop[0]] for step in PartDict]
    y = [step[prop[1]] for step in PartDict]

    # Color anpassen...
    # cmap = mpl.cm.get_cmap('tab20c')
    # ax.set_prop_cycle(color=[cmap(6/20),cmap(7/20)])

    if 'color' in kwargs:
        ax.plot(x,y,linestyle=styles[i],c=kwargs['color'],label='AvaNode '+AvaNodeNames[i],linewidth=1.2)
        ax.plot(x,y,linestyle=styles[i],c='k',linewidth=0.1)
    else:
        ax.plot(x,y,'b--',label='AvaNode '+AvaNodeNames[i],linewidth=0.5)
