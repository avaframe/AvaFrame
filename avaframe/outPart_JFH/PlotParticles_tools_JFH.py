"""
Script for plotting particle results

Created on Tue Nov 28 2022
@author: Höller JF
"""


# Load modules
import numpy as np                                  # Vektormathematik
import matplotlib as mpl                            # Get / use colormaps
import matplotlib.pyplot as plt                     # Plot created plots
import matplotlib.ticker as tick

# Temporary modules
import scipy as sp                                  # Interpolation and datafiltering


# Local imports
import PlotParticles_auxiliary_PF_JFH as ppap       # Auxiliary functions for plotting
import PlotParticles_auxiliary_MF_JFH as ppam       # Auxiliary functions for plotting
import avaframe.com1DFA.DFAtools as DFAtls          # e.g. calculate Norm
import avaframe.out3Plot.outCom1DFA as outCom1DFA   # Plot tracectory plot


####   ----------------------------------------------------------------------------   ###
#                                     Plot Functions                                    #
####   ----------------------------------------------------------------------------   ###



####   ----------------------------------------------------------------------------   ###
#                                           4                                           #
#      Geschwindigkeitsdichte (Unterschiede Simulation und AvaNode herausarbeiten)      #
#                                                                                       #
####   ----------------------------------------------------------------------------   ###


def CompareHistogram2d(particleList,property='umag',nBins=50,plotStats=False,LatexExport=False,**kwargs):
    """ Create Histogram in every timestep

    Parameters
    ----------
    particleList: list of dict
        list of dictionaries
    property: str or (list of str)
        property of the y-axis (e.g. "umag")
    nBins: int or [int int]
        The bin specification
    **kwargs:
        SimName: str
            Name of the simulation
        ProjectName: str
            Name of the project
        avaDir: str
            path to avalanche directory
        AvaNodes: list
            list of AvaNode dictionaries
    """
    if isinstance(particleList,dict):
        return
    
    if not isinstance(property,list):
        property = [property]

    for prop in property:
        if not prop in particleList[0]:
            continue

        fig,axs = plt.subplots(1, 2,figsize=[5.34213, 2.7])    #,figsize=[5.34213, 2.8]
        
        if 'SimName' in kwargs:
            fig.suptitle('Comparison of the '
            + ppap.PropertyLabel(prop) +
            ' in the time- and thalweg-coordinatesystem\nSimulation: ' +
            kwargs['SimName'])
        elif not LatexExport:
            fig.suptitle('Comparison of the ' + ppap.PropertyLabel(prop) + ' in the time- and thalweg-coordinatesystem')

        if 'ProjectName' in kwargs:
            ppap.PlotText(axs,kwargs['ProjectName'],loc='upper right')


        #######################     Time-Diagram     #######################

        DataX = [[i['t']]*i['nPart'] for i in particleList]
        DataY = [i[prop] for i in particleList]
        xBins = [i['t'] for i in particleList]
        xBins.append(2*xBins[-1]-xBins[-2])         # Add last bin edge for histogram

        # Flatten list of arrays:
        DataXflat1 = ppam.flattenList(DataX)
        DataYflat1 = ppam.flattenList(DataY)
        nyBins = int(len(particleList))               # Anzahl Kästchen in y-Richtung

        # Calculate histogram
        X1,Y1,H1 = ppap.CalcHist2d(DataXflat1,DataYflat1,[xBins,nyBins])


        #######################      Thalwegdiagram     #######################

        DataX2 = [i['sAimec'] for i in particleList]
        DataY2 = [i[prop] for i in particleList]

        # Flatten list of arrays:
        DataXflat2 = ppam.flattenList(DataX2)
        DataYflat2 = ppam.flattenList(DataY2)

        # Calculate histogram
        # nBins=200
        # nxBins2 = np.linspace(0,500,nBins)

        X2,Y2,H2 = ppap.CalcHist2d(DataXflat2,DataYflat2,nBins)    # nBins # [nxBins2,nyBins]
        # Calculate vmax
        vmax = max(H1.max(),H2.max())

        # Plot both histograms
        ppap.PlotHist2D(axs[0],X1,Y1,H1,fig,colorbar=False,vmax=vmax)
        ppap.PlotHist2D(axs[1],X2,Y2,H2,fig,colorbar=True,vmax=vmax)

        # Plot Quantiles of AvaFrame particles
        if plotStats:
            ppap.plotStats(axs[0],DataXflat1,DataYflat1,q=[0, 0.5, 1])
            ppap.plotStats(axs[1],DataXflat2,DataYflat2,q=[0, 0.5, 1],label=False)
            
            lines_labels = [ax.get_legend_handles_labels() for ax in fig.axes]
            lines, labels = [sum(lol, []) for lol in zip(*lines_labels)]
            plt.figlegend(lines, labels, loc = 'lower center', ncol=3, labelspacing=0.)


        # # Calculate standard deviation
        # sd = [np.std(i) for i in DataY]
        # fig2,ax2 = plt.subplots()
        # ax2.plot(sd)
        # # bzw. Thalweg:
        # sX,sY = ppam.SortLists(DataXflat2,DataYflat2)
        # sd2 = [np.std(i) for i in sY]
        # fig3,ax3 = plt.subplots()
        # ax3.plot(sX,sd2)


        # Calculate Max Y-Value and corresponding X-Value
        # max(DataYflat1)
        # DataXflat1[(DataYflat1.index(max(DataYflat1)))]   # t_max
        # DataXflat2[(DataYflat2.index(max(DataYflat2)))]   # s_xy,max


        #######################      AvaNode-results     #######################

        if 'AvaNodes' in kwargs:
            stats=True
            for i,xresult in enumerate(['t','sAimec']):
                ppap.plotAvaNodeStats(axs[i],kwargs['AvaNodes'],xresult,prop,stats=stats)
            # Centre legend
            lines_labels = [ax.get_legend_handles_labels() for ax in fig.axes]
            lines, labels = [sum(lol, []) for lol in zip(*lines_labels)]
            if stats:
                # Reorder plot labels
                order = [3, 0, 1, 2]
                # pass handle & labels lists along with order as below
                plt.figlegend([lines[0:4][i] for i in order], [labels[0:4][i] for i in order], loc = 'lower center', ncol=4, labelspacing=0.)
            else:
                #plt.figlegend(lines[0:8], labels[0:8], loc = 'lower center', ncol=4, labelspacing=0.,borderaxespad=0.1)
                #plt.subplots_adjust(bottom=0.85)
                # Shrink current axis's height by 10% on the bottom
                box = axs[0].get_position()
                axs[0].set_position([box.x0, box.y0, box.width, box.height * 0.7])
                box = axs[1].get_position()
                axs[1].set_position([box.x0, box.y0, box.width, box.height * 0.7])

                # Put a legend below current axis
                # axs[0].legend(loc='upper center',
                #               bbox_to_anchor=(0.5, -0.05),
                #               fancybox=True, ncol=4)
                
                plt.figlegend(lines[0:8], labels[0:8], loc = 'lower center', ncol=4, labelspacing=0.,borderaxespad=0.1)



        #######################      Plot settings     #######################

        ppap.LablePlot(axs[0],ppap.PropertyLabel('t'),ppap.PropertyLabel(prop),'Zeitliche Entwicklung')        # temporal evolution / Zeitliche Entwicklung
        ppap.LablePlot(axs[1],ppap.PropertyLabel('sAimec'),ppap.PropertyLabel(prop),'Räumliche Entwicklung')    # spacial evolution / Räumliche Entwicklung
        

        if 'AvaDir' in kwargs:
            if 'AvaNodes' in kwargs:
                ppap.savePlot(fig,kwargs['ProjectName']+'_2DHist_AR_'+prop,kwargs['AvaDir'],comModule="Aimec",LatexExport=LatexExport)
            else:
                ppap.savePlot(fig,kwargs['ProjectName']+'_2DHist_'+prop,kwargs['AvaDir'],comModule="Aimec",LatexExport=LatexExport)




####   ----------------------------------------------------------------------------   ###
#                                           5                                           #
#     Max Geschwindigkeit im Auslösegebiet gegenüber der absoluten Fallhöhe stellen     #
#                Linker Plot Geschw. <-> Rechter Plot absolute Fallhöhe                 #
#                                                                                       #
####   ----------------------------------------------------------------------------   ###


def CompareTrajectory(particleList,property,LatexExport=False,**kwargs):
    """ Compare relativ distance of the particles to the centre of mass in the release area
    to the relativ distance of the particles to the centre of mass in the final area

    Parameters
    ----------
    particleList: list of dict
        list of dictionaries
    property: str
        property for the result visualisation
    LatexExport: bool
        Decide whether the graphics should be exported for use in Latex
    **kwargs:
        SimName: str
            Name of the simulation
        ProjectName: str
            Name of the project
        avaDir: str
            path to avalanche directory
    """
    if not property in particleList[0]:
        print('Particle property %s is not defined!' % (property))
        exit()
    
    fig, axs = plt.subplots(1,2,layout="constrained")           # figsize=[8.5,4] -> moved to ppap.PlotForLatex()

    if 'SimName' in kwargs:
        fig.suptitle('Max. particle results\nSimulation: ' + kwargs['SimName'])
    elif not LatexExport:
        fig.suptitle('Max. particle results')


    Tri = True
    statAusw = False
    # Linien entsprechend der Ergebnisse einfärben
    match property:
        case 'umag':
            maxValues = True
            colormap = 'Spectral_r'
            match kwargs['ProjectName']:
                case 'avaParabola':
                    min_z = 0
                    max_z = 40
                case 'avaHelix':
                    min_z = 0
                    max_z = 40
                case 'avaNordkette':
                    min_z = 0
                    max_z = 25
                case 'avaNordketteVU':      # max = 19.146
                    Tri = False
                    min_z = 0
                    max_z = 20
        case 'sAimec':
            maxValues = True
            colormap = 'plasma'
            # colormap = 'PuBu'
            min_z = 200
            max_z = 1200
            match kwargs['ProjectName']:
                case 'avaParabola':
                    min_z = 200
                    max_z = 500
                case 'avaHelix':
                    min_z = 300
                    max_z = 1100
                case 'avaNordkette':
                    min_z = 250
                    max_z = 550
                case 'avaNordketteVU':      # max = 432.32 m
                    min_z = 0
                    max_z = 500
        case 'dsCoM':
            maxValues = False
            colormap = 'RdYlBu_r'
            match kwargs['ProjectName']:
                case 'avaParabola':
                    min_z = -200
                    max_z = 200
                case 'avaHelix':
                    min_z = -500
                    max_z = 500
                case 'avaNordkette':
                    min_z = -250
                    max_z = 250
                case 'avaNordketteVU':
                    min_z = -320
                    max_z = 320


    if 'ProjectName' in kwargs:
        ppap.PlotText([axs[0],axs[1]],kwargs['ProjectName'])


    vmaxbool = True
    if vmaxbool:
        scat1 = ppap.ScatterPlot(axs[0],particleList,property,t=0,max=maxValues,cmap=colormap,triangulation=True,vmin=min_z,vmax=max_z)
        scat2 = ppap.ScatterPlot(axs[1],particleList,property,t=-1,max=maxValues,cmap=colormap,triangulation=Tri,vmin=min_z,vmax=max_z)

    else:
        scat1 = ppap.ScatterPlot(axs[0],particleList,property,t=0,max=maxValues,cmap=colormap,triangulation=True)
        scat2 = ppap.ScatterPlot(axs[1],particleList,property,t=-1,max=maxValues,cmap=colormap,triangulation=True)


    if property=='dsCoM':
        axs[0].scatter(particleList[0]['lCoM'],particleList[0]['sCoM'],c='r',label='CoM')
        axs[1].scatter(particleList[-1]['lCoM'],particleList[-1]['sCoM'],c='r',label='CoM')




    # Add AvaRange Results
    if 'AvaNodes' in kwargs:
        statAusw = True
        # Settings:
        AvaNodeNames = ['1.1','2.1','2.3','3.1','4.10','5.7','5.9','5.10']
        alpha = 0.7
        size = 2
        color = 'k'

        if maxValues:
            z = np.array([ppam.MaxList(AR,property) for AR in kwargs['AvaNodes']])

        for i,AR in enumerate(kwargs['AvaNodes']):
            # Release Area:
            x1 = AR[0]['lAimec']
            y1 = AR[0]['sAimec']
            # Run-Out Area:
            x2 = AR[-1]['lAimec']
            y2 = AR[-1]['sAimec']

            if maxValues:
                j = (z[i]-min_z)/(max_z-min_z)

                circle = plt.Circle((x1, y1),2.5,alpha=alpha,color=scat1.cmap(j))
                axs[0].add_patch(circle)
                axs[0].scatter(x1,y1,c=z[i],cmap=colormap,s=3,vmin=min_z,vmax=max_z)
                axs[0].scatter(x1,y1,color=color,s=1)

                if property=='sAimec':
                    #axs[1].scatter(x2,y2,color=color,s=size)
                    print('test')
                else:
                    circle = plt.Circle((x2,y2),10,alpha=alpha,color=scat1.cmap(j))
                    axs[1].add_patch(circle)
                    axs[1].scatter(x2,y2,c=z[i],cmap=colormap,s=3,vmin=min_z,vmax=max_z)
                axs[1].scatter(x2,y2,color=color,s=1)
            else:
                axs[0].scatter(x1,y1,color=color,marker='^',s=size)
                axs[1].scatter(x2,y2,color=color,marker='^',s=size)

            # Label AvaNodes:
            axs[0].annotate(AvaNodeNames[i],(x1+1.5,y1-1),fontsize=4,color=color)
            lst = [2,5,6,1]
            if i in lst:
                axs[1].annotate(AvaNodeNames[i],(x2,y2),color=color,xytext=(x2+90, y2-10),arrowprops=dict(linewidth=0.5,color=color,mutation_scale=4,
                    arrowstyle="->"),fontsize=4)
            else:
                axs[1].annotate(AvaNodeNames[i],(x2,y2),color=color,xytext=(x2-60, y2-10),arrowprops=dict(linewidth=0.5,color=color,mutation_scale=4,
                    arrowstyle="->"),fontsize=4)


    # Maximalwerte der Avanodes:
    # max([max([step['umag'] for step in p]) for p in kwargs['AvaNodes']])
    # Postion im Auslösegebiet:
    # ARposition = [(p[0]['sAimec'],p[0]['lAimec']) for p in kwargs['AvaNodes']]


    if vmaxbool:
        ticks = np.linspace(min_z,max_z,5)  #6
        cbar = fig.colorbar(scat2,ax=axs[-1],ticks=ticks)           # If colorbar is not properly formated try ax=axs
    else:
        cbar = fig.colorbar(scat2,ax=axs[-1])
    # cbar.ax.yaxis.set_major_formatter(tick.FormatStrFormatter('%.2f'))

    if maxValues:
        cbar.set_label(ppap.PropertyLabel(property+'_max'))
    else:
        cbar.set_label(ppap.PropertyLabel(property))


    # Plot settings
    language = 'deutsch'
    match language:
        case 'english':
            axs[0].set_title(r'particles in release area $t = \SI{0}{\s}$')
            axs[1].set_title(r'particles in run-out area $t = \SI{'
                            + str(round(particleList[-1]['t'],1))
                            + r'}{\s}$')
        case 'deutsch':
            axs[0].set_title(r'Auslösegebiet $t = \SI{0}{\s}$')
            axs[1].set_title(r'Ablagerungsgebiet $t = \SI{'
                            + str(round(particleList[-1]['t'],1))
                            + r'}{\s}$')


    if property == 'sAimec':
        cbar.ax.invert_yaxis()

    if property == 'dsCoM':
        axs[0].legend()
        axs[1].legend()
        cbar.ax.invert_yaxis()

    for a in axs:
        ppap.PlotSettings(a,aspect='centre')
        a.invert_yaxis()
        if kwargs['ProjectName']=='avaNordketteVU':
            a.invert_xaxis()
        ppap.LablePlot(a,ppap.PropertyLabel('lAimec'),ppap.PropertyLabel('sAimec'))

    if statAusw:
        # Stauchung und Dehnung der Lawine berechnen:
        x0 = np.array(particleList[0]['lAimec'])
        y0 = np.array(particleList[0]['sAimec'])
        x1 = np.array(particleList[-1]['lAimec'])
        y1 = np.array(particleList[-1]['sAimec'])
        kx = (x1.max()-x1.min())/(x0.max()-x0.min())
        ky = (y1.max()-y1.min())/(y0.max()-y0.min())
        print('kx = '+str(kx))
        print('ky = '+str(ky))
        
        # Auswertung des relativen Abstands zum CoM:
        ds_AR=np.array([p[0]['dsCoM'] for p in kwargs['AvaNodes']])
        de_AR=np.array([p[-1]['dsCoM'] for p in kwargs['AvaNodes']])
        d_AR = abs(de_AR)-abs(ds_AR)
        median_AR = np.quantile(d_AR,0.5)

        # Numerische Partikel
        ds_part = particleList[0]['dsCoM']
        de_part = particleList[-1]['dsCoM']
        d_part = abs(de_part)-abs(ds_part)
        # Median
        median_part = np.quantile(d_part,0.5)

        nKleiner = zaehle_negative_elemente(d_part)

        # Berechnung des Korrelationskoeffizienten
        korrelation_AR = np.corrcoef(ds_AR, de_AR)[0, 1]
        korrelation_part = np.corrcoef(ds_part, de_part)[0, 1]

        print("Korrelationskoeffizient der AvaNodes:", korrelation_AR)
        print("Korrelationskoeffizient der numerischen Partikel:", korrelation_part)


    # plt.show()
    if 'AvaDir' in kwargs:
        if 'AvaNodes' in kwargs:
            ppap.savePlot(fig,kwargs['ProjectName']+'_ARandSim_'+property,kwargs['AvaDir'],comModule="Aimec",LatexExport=LatexExport)
        else:
            ppap.savePlot(fig,kwargs['ProjectName']+'_'+property,kwargs['AvaDir'],comModule="Aimec",LatexExport=LatexExport)



# Anzahl der Elemente, bei denen der Abstand kleiner wird:
def zaehle_negative_elemente(liste):
    zaehler = 0
    for element in liste:
        if element < 0:
            zaehler += 1
    return zaehler



####   ----------------------------------------------------------------------------   ###
#                                           6                                           #
#                                     Sankey Chart                                      #
#                                                                                       #
####   ----------------------------------------------------------------------------   ###



def SankeyChart(particleList,property='umag',colorscheme='max',LatexExport=False,**kwargs):
    """ Create Sankey flowchart for visualising particle flow
    
    Parameter
    ---------
    particleList: list of dict
        list of particle dictionaries
    property: str
        property for the result visualisation
    colorscheme: str
        'max',
    LatexExport: bool
        Decide whether the graphics should be exported for use in Latex
    **kwargs:
        ProjectName: str
            Name of the project
        AvaDir: str
    """
    fig,ax = plt.subplots(figsize=[5.34213, 3])

    alpha = 0.1

    if 'ProjectName' in kwargs:
        ppap.PlotText(ax,kwargs['ProjectName'],loc='upper right')
        match kwargs['ProjectName']:
            case 'avaParabola':
                alpha = 0.05
            case 'avaHelix':
                alpha = 0.02
            case 'avaWog':
                alpha = 0.005
            case 'avaNordkette':
                alpha = 0.7
            case 'avaNordketteVU':
                alpha = 0.15
    
    #######################      AvaFrame-results     #######################
    
    # Linien entsprechend der Ergebnisse einfärben
    match property:
        case 'umag':
            colormap = 'Spectral_r'
            match kwargs['ProjectName']:
                case 'avaParabola':
                    min_z = 0
                    max_z = 40
                case 'avaHelix':
                    min_z = 0
                    max_z = 40
                case 'avaNordkette':
                    min_z = 0
                    max_z = 25
                case 'avaNordketteVU':
                    min_z = 0
                    max_z = 20
        case 'sAimec':
            colormap = 'plasma'
            #colormap = 'PuBu'
            match kwargs['ProjectName']:
                case 'avaParabola':
                    min_z = 200
                    max_z = 500
                case 'avaHelix':
                    min_z = 300
                    max_z = 1100
                case 'avaWog':
                    min_z = 0
                    max_z = 500
                case 'avaNordkette':
                    min_z = 250
                    max_z = 550
                case 'avaNordketteVU':
                    min_z = 0
                    max_z = 500
                    alpha = 0.3
        case 'dsCoM':
            colormap = 'RdYlBu_r'            
            match kwargs['ProjectName']:
                case 'avaParabola':
                    min_z = -200
                    max_z = 200
                case 'avaHelix':
                    min_z = -500
                    max_z = 500
                case 'avaWog':
                    min_z = 0
                    max_z = 500
                case 'avaNordkette':
                    min_z = -250
                    max_z = 250
                case 'avaNordketteVU':
                    min_z = -320
                    max_z = 320


    match colorscheme:
        case 'max':
            z = ppam.MaxList(particleList,property)
            cmap = mpl.cm.get_cmap(colormap)
            znorm = (z-z.min())/(z.max()-z.min())
            norm = mpl.colors.Normalize(vmin=z.min(),vmax=z.max())
            sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
            sm.set_array([])

        case 'explicit':
            z = ppam.MaxList(particleList,property)
            cmap = mpl.cm.get_cmap(colormap)
            znorm = (z-min_z)/(max_z-min_z)
            norm = mpl.colors.Normalize(vmin=min_z,vmax=max_z)
            sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
            sm.set_array([])
        
    if kwargs['ProjectName'] == 'avaNordketteVU':
        flawelem = [37,38,88,89,90,91,92,93,152,153,154,
                    157,158,160,165,166,222,223,224,225,226,
                    227,228,229,291]

    t = [step['t'] for step in particleList]
    for i in range(particleList[0]['nPart']):
        y = [step['sAimec'][i] for step in particleList]
        if kwargs['ProjectName'] == 'avaNordketteVU':
            if any(element == i for element in flawelem):
                continue

        ax.plot(t,y,c=cmap(znorm[i]),linewidth=0.7,alpha=alpha)       # Parabola: 0.01, Helix: 0.02, avaWog:0.005 ,Nordkette: 0.15

    if colorscheme=='explicit':
        ticks = np.linspace(min_z,max_z,5)
        cbar = fig.colorbar(sm,ax=ax,ticks=ticks)
    else:
        cbar = fig.colorbar(sm,ax=ax)
    cbar.set_label(ppap.PropertyLabel(property+'_max'))

    cbar.solids.set(alpha=1)    # alpha


    #######################      AvaNode-results     #######################
    
    if 'AvaNodes' in kwargs:        
        for i,AR in enumerate(kwargs['AvaNodes']):
            # Get max result
            zmax = max([a[property] for a in AR])
            # zmax = (zmax-z.min())/(z.max()-z.min())
            zmax = (zmax-min_z)/(max_z-min_z)
            ppap.plotNodeDict(ax,AR,['t','sAimec'],i,color=cmap(zmax),alpha=alpha)
        ax.legend()

    
    #######################      Centre of Mass     #######################
    if property == 'dsCoM':
        CoM = [i['sCoM'] for i in particleList]
        ax.plot(t,CoM,'--r',linewidth=1,label='CoM')
        ax.legend(loc='lower left')

    
    #######################      Plot settings     #######################
    if property=='dsCoM' or property=='sAimec':
        cbar.ax.invert_yaxis()

    ax.invert_yaxis()
    ppap.LablePlot(ax,'t [s]',ppap.PropertyLabel('sAimec'))
    # fig.set_rasterized(True)
    # plt.show()
    if 'AvaDir' in kwargs:
            if 'AvaNodes' in kwargs:
                ppap.savePlot(fig,kwargs['ProjectName']+'_TT_AR_'+property,kwargs['AvaDir'],comModule="Aimec",LatexExport=LatexExport)
            else:
                ppap.savePlot(fig,kwargs['ProjectName']+'_TT_'+property,kwargs['AvaDir'],comModule="Aimec",LatexExport=LatexExport)




def ThalwegTime(particleList,property='umag',colorscheme='max',LatexExport=False,**kwargs):
    """ Create Thalweg-Time diagramm for visualising particle flow
    
    Parameter
    ---------
    particleList: list of dict
        list of particle dictionaries
    property: str
        property for the result visualisation
    colorscheme: str
        'max',
    LatexExport: bool
        Decide whether the graphics should be exported for use in Latex
    **kwargs:
        ProjectName: str
            Name of the project
        AvaDir: str
    """
    fig,ax = plt.subplots(figsize=[5.34213, 3])
    alpha = 0.1
    s = 12

    if 'ProjectName' in kwargs:
        ppap.PlotText(ax,kwargs['ProjectName'],loc='upper right')
        match kwargs['ProjectName']:
            case 'avaParabola':
                alpha = 0.05
            case 'avaHelix':
                alpha = 0.02
            case 'avaWog':
                alpha = 0.005
            case 'avaNordkette':
                alpha = 1
            case 'avaNordketteVU':
                alpha = 0.15

        
    #######################      AvaFrame-results     #######################
    
    # Color lines according to the results
    match property:
        case 'umag':
            colormap = 'Spectral_r'
            min_z = 5
            max_z = 35
        case 'sAimec':
            colormap = 'PuBu'
            #colormap = 'YlGnBu'
            min_z = 200
            max_z = 1200
        case 'dsCoM':
            colormap = 'RdYlBu_r'
            match kwargs['ProjectName']:
                case 'avaParabola':
                    min_z = -200
                    max_z = 200
                case 'avaHelix':
                    min_z = -500
                    max_z = 500
                case 'avaWog':
                    min_z = 0
                    max_z = 500
                case 'avaNordkette':
                    min_z = -250
                    max_z = 250
                case 'avaNordketteVU':
                    min_z = -320    # min -321.76
                    max_z = 320

    # Ignore flawed elements
    if kwargs['ProjectName'] == 'avaNordketteVU':
        flawelem = [37,38,88,89,90,91,92,93,152,153,154,
                    157,158,160,165,166,222,223,224,225,226,
                    227,228,229,291]
        t = [step['t'] for step in particleList]*(particleList[0]['nPart']-len(flawelem))
    else:
        t = [step['t'] for step in particleList]*particleList[0]['nPart']
        
    y = []
    c = []
    for i in range(particleList[0]['nPart']):
        if kwargs['ProjectName'] == 'avaNordketteVU':
            if any(element == i for element in flawelem):
                continue
            
        y.extend([step['sAimec'][i] for step in particleList])
        c.extend([step[property][i] for step in particleList])

    
    vmax = True
    if vmax:
        scat=ax.scatter(t,y,c=c,s=s,cmap=colormap,vmin=min_z,vmax=max_z,alpha=alpha,edgecolors='none')
        ticks = np.linspace(min_z,max_z,5)  #6
        cbar = fig.colorbar(scat,ax=ax,ticks=ticks)           # If colorbar is not properly formated try ax=axs
    else:
        scat=ax.scatter(t,y,c=c,s=s,cmap=colormap,alpha=alpha,edgecolors='none')
        cbar = fig.colorbar(scat,ax=ax)
    
    
    if property == 'dsCoM':
        CoM = [i['sCoM'] for i in particleList]
        ax.plot(t[0:len(particleList)],CoM,'--r',linewidth=1,label='CoM')
        ax.legend(loc='upper right')



    #######################      AvaNode-results     #######################
    
    if 'AvaNodes' in kwargs:
        AvaNodeNames = ['1.1','2.1','2.3','3.1','4.10','5.7','5.9','5.10']
        styles = ['-',':',(0, (1, 10)),'--',(5, (10, 3)),'-.',(0, (3, 5, 1, 5, 1, 5)),(0, (3, 1, 1, 1, 1, 1))]
        for i,AR in enumerate(kwargs['AvaNodes']):
            tAR = []
            yAR = []
            cAR = []
            for j,tstep in enumerate(AR):
                tAR.append(tstep['t'])
                yAR.append(tstep['sAimec'])
                cAR.append(tstep[property])
            scatAR=ax.scatter(tAR,yAR,c=cAR,s=s,cmap=colormap,vmin=min_z,vmax=max_z,alpha=1,edgecolors='none')
            ax.plot(tAR,yAR,linestyle=styles[i],label='AvaNode '+AvaNodeNames[i],color='grey',linewidth=0.2)


    #######################      Plot settings     #######################


    if property=='dsCoM' or property=='sAimec':
        cbar.ax.invert_yaxis()

    #cbar.ax.yaxis.set_major_formatter(tick.FormatStrFormatter('%.2f'))
    
    cbar.solids.set(alpha=1)
    cbar.set_label(ppap.PropertyLabel(property))
    # or: # cbar.set_alpha(1)
    # cbar.draw_all()


    ax.invert_yaxis()
    ax.legend(loc='lower left')
    ppap.LablePlot(ax,'t [s]',ppap.PropertyLabel('sAimec'))

    # plt.show()
    if 'AvaDir' in kwargs:
            if 'AvaNodes' in kwargs:
                ppap.savePlot(fig,kwargs['ProjectName']+'_TT_AR_'+property,kwargs['AvaDir'],comModule="Aimec",LatexExport=LatexExport)
            else:
                ppap.savePlot(fig,kwargs['ProjectName']+'_TT_'+property,kwargs['AvaDir'],comModule="Aimec",LatexExport=LatexExport)





####   ----------------------------------------------------------------------------   ###
#                                           7                                           #
#                                   Plot Particles xy                                   #
#                                                                                       #
####   ----------------------------------------------------------------------------   ###



def plotParticles(particleList,AvaNodes,**kwargs):
    """ Plot all particles in the release and runout area and add the AvaNodes

    Parameters
    ----------
    particleList: list of dict
        list of dictionaries
    AvaNodes: list
    dem: dict
        dem dictionary          -> nur ergenzen, wenn Datstellung im xy-Koordinatensystem
    **kwargs:
        SimName: str
            Name of the simulation
        ProjectName: str
            Name of the project
        radius: float
            radius of the particle area
        dem: dict
            dem dictionary (normal information will be added in this function)
        AvaDir: str
            path to avalanche directory
    """    
    fig, ax = plt.subplots(layout="constrained",figsize = [5.34213,5.34213])           # figsize=[8.5,4] -> moved to ppap.PlotForLatex()

    if 'SimName' in kwargs:
        fig.suptitle('Max. particle results\nSimulation: ' + kwargs['SimName'])

    if 'ProjectName' in kwargs:
        ppap.PlotText(ax,kwargs['ProjectName'],loc='upper right')


    # Add dem
    if 'dem' in kwargs:
        # Transform DEM, 0,0 in the left bottom corner
        kwargs['dem']['header']['xllcenter'] = 0
        kwargs['dem']['header']['yllcenter'] = 0
        ax = outCom1DFA.addDem2Plot(ax,kwargs['dem'],what='slope')

    ax.scatter(particleList[0]['x'],particleList[0]['y'],s=1.5,label='Partikel der Simulation')     #,color='k'
    
    # Add AvaNodes as circles
    AvaNodeNames = ['1.1','2.1','2.3','3.1','4.10','5.7','5.9','5.10']
    # AvaColor = [0,4,5,8,12,16,17,18]
    AvaColor = [16,4,6,8,12,0,1,2]
    cmap = mpl.cm.get_cmap('tab20c')
    for i,AR in enumerate(AvaNodes):
        circle = plt.Circle((AR[0]['x'], AR[0]['y']), kwargs['radius'],alpha=0.3,color=cmap(AvaColor[i]))     # color='r'
        ax.add_patch(circle)
        ax.scatter(AR[0]['x'], AR[0]['y'],s=3,label='AvaNode '+AvaNodeNames[i],color=cmap(AvaColor[i]))


    # Plot settings
    language = 'deutsch'
    match language:
        case 'english':
            ax.set_title(r'particles in release area $t = \SI{0}{\s}$')
        case 'deutsch':
            ax.set_title(r'Auslösegebiet $t = \SI{0}{\s}$')

    ax.set_aspect('equal')

    ppap.LablePlot(ax,'x [m]','y [m]')
    ax.legend()

    if 'AvaDir' in kwargs:
        ppap.savePlot(fig,kwargs['ProjectName']+'_particles',kwargs['AvaDir'],comModule="Aimec",LatexExport=True)





####   ----------------------------------------------------------------------------   ###
#                                           8                                           #
#                                  Plot distance to CoM                                 #
#                                                                                       #
####   ----------------------------------------------------------------------------   ###



def plotdsCoM(AvaNodes,**kwargs):
    """ Plot the distance between the particles and the CoM

    Parameters
    ----------
    AvaNodes: list of dict
        list of particle dictionaries
    particleList: list of dict
        list of particle dictionaries
    **kwargs:
        ProjectName: str
            Name of the project
        AvaDir: str
    """
    fig,ax = plt.subplots(figsize = [5.34213,3])
    
    # get max time of AvaNode results
    tmax = max([AR[-1]['t'] for AR in AvaNodes])

    AvaNodeNames = ['1.1','2.1','2.3','3.1','4.10','5.7','5.9','5.10']
    styles = ['-',':',(0, (1, 10)),'--',(5, (10, 3)),'-.',(0, (3, 5, 1, 5, 1, 5)),(0, (3, 1, 1, 1, 1, 1))]
    AvaColor = [16,4,6,8,12,0,1,2]
    cmap = mpl.cm.get_cmap('tab20c')
    
    for i,particles in enumerate(AvaNodes):
        x = [tstep['t'] for tstep in particles]
        y = [tstep['dsCoM'] for tstep in particles]

        # Interpolate AR results to uniform time steps
        y_smooth = sp.signal.savgol_filter(y, window_length=30, polyorder=3, mode="nearest")

        ax.plot(x,y_smooth,color=cmap(AvaColor[i]),label='AvaNode '+AvaNodeNames[i])      #*zip(*res_AR)
        # ax.plot(x,y_smooth,color=cmap(AvaColor[i]),linestyle=styles[i],label='AvaNode '+AvaNodeNames[i])      #*zip(*res_AR)
    
    ax.plot([0,tmax],[0,0],'r--',linewidth=1,label='CoM')


    #######################      Plot settings     #######################
    ax.invert_yaxis()
    ax.legend(loc='upper left')
    ppap.LablePlot(ax,'t [s]',ppap.PropertyLabel('dsCoM'))

    # plt.show()
    if 'AvaDir' in kwargs:
        ppap.savePlot(fig,kwargs['ProjectName']+'_t_dsCoM',kwargs['AvaDir'],comModule="Aimec",LatexExport=True)



def plotdsCoMnumeric(particleList,**kwargs):
    """
    
    """
    fig,ax = plt.subplots(figsize = [5.34213,3])
    ds_part = particleList[0]['dsCoM']
    de_part = particleList[-1]['dsCoM']
    # Calculate Correlation
    correlation_part = np.corrcoef(ds_part, de_part)[0, 1]

    ax.scatter(ds_part,de_part,s=1,label='numerische Partikel: r = '+str(round(correlation_part,2)))
    ax.scatter(0,0,c='r',label='CoM')

    if 'AvaNodes' in kwargs:
        ds_AR=np.array([p[0]['dsCoM'] for p in kwargs['AvaNodes']])
        de_AR=np.array([p[-1]['dsCoM'] for p in kwargs['AvaNodes']])
        correlation_AR = np.corrcoef(ds_AR, de_AR)[0, 1]
        # np.corrcoef(np.concatenate((ds_AR[:2], ds_AR[3:])), np.concatenate((de_AR[:2], de_AR[3:])))

        ax.scatter(ds_AR,de_AR,s=4,label='AvaNodes: r ='+str(round(correlation_AR,2)))

    
    # ,'Korrelationskoeffizient r = '+str(round(correlation_part,2))

    #######################      Plot settings     #######################
    #ax.invert_yaxis()
    ax.legend(loc='lower right')
    ppap.LablePlot(ax,r'rel. distance $d_\text{s,CoM,start}$ [m]',r'rel. distance $d_\text{s,CoM,end}$ [m]')

    # plt.show()
    if 'AvaDir' in kwargs:
        ppap.savePlot(fig,kwargs['ProjectName']+'_cor_dsCoM',kwargs['AvaDir'],comModule="Aimec",LatexExport=True)


