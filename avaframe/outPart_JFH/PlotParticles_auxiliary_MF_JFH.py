"""
Auxiliary functions for calculating the particle results

Created on Tue Jan 09 2023
@author: Höller JF
"""


# Load modules
import numpy as np                                  # Vektormathematik
import pathlib                                      # Work with smart paths
import pandas as pd                                 # Read CSV files
from scipy.interpolate import CubicSpline           # Cubische Spline Interpolation 1D

# Temporary modules

# Local imports
from avaframe.in3Utils import cfgUtils              # Load GeneralConfiguration
from avaframe.ana3AIMEC import dfa2Aimec, ana3AIMEC, aimecTools
import avaframe.in2Trans.ascUtils as IOf            # Read raster of a dem
import avaframe.ana3AIMEC.aimecTools as aT
import avaframe.in3Utils.geoTrans as geoTrans       # Transform coordinates


####   ----------------------------------------------------------------------------   ###
#                                  Auxiliary Math Functions                             #
####   ----------------------------------------------------------------------------   ###


def velocityMagnitude(particlesList):
    """ Calculation of the maximum velocity in each time step

    Parameters
    ----------
    particleList: list of dict
        list of dictionaries

    Returns
    ----------
    particlesList: list of dict
        updated list of dictionaries
    """
    UX = np.array([d['ux'] for d in particlesList],dtype='object')
    UY = np.array([d['uy'] for d in particlesList],dtype='object')
    UZ = np.array([d['uz'] for d in particlesList],dtype='object')
    UXYZ = (UX*UX + UY*UY + UZ*UZ)
    UXYZ = UXYZ.tolist()                                            # Change type, because np only supports nested arrays with equal lengths
    UMAG = [np.sqrt(d) for d in UXYZ]

    # UX = [d['ux'] for d in particlesList]
    # UY = [d['uy'] for d in particlesList]
    # UZ = [d['uz'] for d in particlesList]

    for i in range(len(UMAG)):
        particlesList[i]['umag'] = UMAG[i]

    return particlesList



def flattenList(List):
    """ Function for flattening lists

    Parameters
    ----------
    List: list
    """
    out = []
    for sublist in List:
        out.extend(sublist)
    return out



def Quantiles(Data,quantiles):
    """ Calculate the median and the q-th quantiles for the data set

    Parameters
    ----------
    Data: list of lists
        list of lists [timesteps][ParticleResult]
    quantiles: array_like or int
        Quantile or sequence of quantiles to compute, which must be between 1 and 100 inclusive.

    Returns
    -------
    Quant: list of lists
        Data Set every list stands for a Quantile
    """
    Quant = [np.quantile(Data[i],quantiles) for i in range(len(Data))]

    # Evtl. Polynom erstellen, dass durch die Punkte verläuft bzw. sich diesen annähert?

    return Quant



def ReducedList(a):
    """ Calculated the middle value of every list element
    
    Parameters
    ----------
    a: list
    """
    b = np.array(a[1:])
    a = np.array(a)

    c = abs(a[:-1]-b)/2
    newa = a[:-1] + c
    return newa



def SortLists(x,y):
    """ Sort list x according to y  
    
    Parameters
    ----------
    x,y: list

    Returns
    -------
    newX: List
    newY: List of lists
    """
    temp = sorted(zip(x,y))

    newX = [temp[0][0]]
    newY = [[]]

    n = 0
    for i in temp:
        if i[0] > newX[-1]:
            n=n+1
            newX.extend([i[0]])
            newY.append([])
        newY[n].extend([i[1]])

    return newX,newY



def MinMax(Data):
    """ Calculate Max and Min of a dictionarry in the first and last row (time step)

    Parameters
    ----------

    Data: array of arrays

    Returns
    -------
    MaxMin: list
        List with the max and min value in the first and last row (time step)
    """
    MinMax = []
    MinMax.append(min(Data[0]))
    MinMax.append(max(Data[0]))
    MinMax.append(min(Data[-1]))
    MinMax.append(max(Data[-1]))

    return MinMax



def CentreOfMass(particleList,coords=['x','y']):
    """ Calculate the centre of mass (weighted average) of a particle cluster
    projected on the xy-plane
    
    Parameters
    ----------
    particleList: list of dict
        Particle dictionary for a specific time step
    coords: list of str
        ['x','y']

    Returns
    -------
    CoM: list
        Centre of mass for multiple particles in spaces
    """
    if isinstance(particleList,dict):
        particleList = [particleList]
    
    for particleDict in particleList:
        m = particleDict['m']

        x = particleDict[coords[0]]
        y = particleDict[coords[1]]
        #z = particleDict[coords[2]]
        xs = np.average(x, weights=m)
        ys = np.average(y, weights=m)
        #zs = np.average(z, weights=m)

        particleDict[coords[0] + 'CoM'] = xs
        particleDict[coords[1] + 'CoM'] = ys

    return particleList



def ARCentreOfMass(particles,Simulation):
    """ Calculate Centre Of Mass in every time step of the measurement
    
    
    """
    t=[]
    lCoM=[]
    sCoM=[]

    for i,tstep in enumerate(Simulation):
        t.append(tstep['t'])
        lCoM.append(tstep['lCoM'])
        sCoM.append(tstep['sCoM'])
    # t=np.array(t)
    # lCoM=np.array(lCoM)
    # sCoM=np.array(sCoM)

    # Cubic Spline for t and lCoM
    csl = CubicSpline(t,lCoM)
    # Cubic Spline for t and sCoM
    css = CubicSpline(t,sCoM)


    for i,tstep in enumerate(particles):
        # interpolate lAimec of CoM
        tstep['lCoM']=float(csl(tstep['t']))
        # interpolate sAimec of CoM
        tstep['sCoM']=float(css(tstep['t']))



def transformSL(particleList,avalancheDir,dem,AvaNode=False):
    """ Transform xy coordinates into sl coordinates

    Parameter
    ---------
    particleList: list of dict
        Particles dictionary for a specific time step
    avalancheDir: str
        path to avalanche directory
    dem: dict
        dem dictionary with header and raster data

    Returns
    -------
    particleList: list of dicts
        updated particle dictionary for a specific time step
    """
    if AvaNode:
        xname = 'x'
        yname = 'y'
        lname = 'lAimec'
        sname = 'sAimec'
    else:
        xname = 'xCoM'
        yname = 'yCoM'
        lname = 'lCoM'
        sname = 'sCoM'

    cfg = cfgUtils.getModuleConfig(ana3AIMEC)
    cfgSetup = cfg['AIMECSETUP']
    anaMod = cfgSetup['anaMod']

    #Setup input from computational module
    inputsDF, resTypeList = dfa2Aimec.mainDfa2Aimec(avalancheDir, anaMod, cfg)

    #define reference simulation
    refSimRowHash, refSimName, inputsDF, colorParameter, valRef = aimecTools.fetchReferenceSimNo(avalancheDir, inputsDF, anaMod,cfg)
    pathDict = {'refSimRowHash': refSimRowHash, 'refSimName': refSimName, 'compType': ['singleModule', anaMod],
        'colorParameter': colorParameter, 'resTypeList': resTypeList, 'valRef': valRef}
    pathDict = aimecTools.readAIMECinputs(avalancheDir, pathDict, dirName=anaMod)
    pathDict = aimecTools.checkAIMECinputs(cfgSetup, inputsDF, pathDict)

    # read reference file and raster and config
    refResultSource = inputsDF.loc[refSimRowHash, cfgSetup['runoutResType']]
    refRaster = IOf.readRaster(refResultSource)
    refHeader = refRaster['header']

    rasterTransfo = aT.makeDomainTransfo(pathDict, dem, refHeader['cellsize'], cfgSetup)

    # Transform xy Schwerpunkt in Thalweg (sl)
    x = np.array([p[xname] for p in particleList])+particleList[0]['xllcenter']
    y = np.array([p[yname] for p in particleList])+particleList[0]['yllcenter']


    for i in range(len(particleList)):
        distance = np.sqrt((x[i]-rasterTransfo['gridx'])**2 + (y[i]-rasterTransfo['gridy'])**2)
        (sIndex, lIndex) = np.unravel_index(np.argmin(distance, axis=None), distance.shape)
        particleList[i][lname] = rasterTransfo['l'][lIndex]
        particleList[i][sname] = rasterTransfo['s'][sIndex]

    return particleList



def CalcDistCoM(particleList,Euclidean=True):
    """ Calculate the direct distance between the particle and the CoM in every time step
    Coordinate system: projected sl (Thalweg)
    
    Parameters
    ----------
    particleList: list of dict
        Particles dictionary for a specific time step
    Euclidean: boolean
        Calculate the euclidean distance between the CoM and the particle?

    Returns
    -------
    particleList: list of dicts
        updated particle dictionary for a specific time step
    """

    for tstep in particleList:
        x = tstep['lAimec']
        y = tstep['sAimec']
        tstep['dlCoM']= tstep['lAimec'] - tstep['lCoM']
        tstep['dsCoM']= tstep['sAimec'] - tstep['sCoM']

        if Euclidean:
            tstep['dAimec'] = EuclideanDist(x,y,[tstep['lCoM'],tstep['sCoM']])

    return particleList



def EuclideanDist(x,y,point,relative=False):
    """ Calculate euclidean distance between a list of points and a centre point
    
    Parameters
    ----------
    x,y: list
        List of coordinates
    point: list
        [x,y]
    
    Returns
    -------
    distances: list
        List of the distances between the points and the centre point
    """
    sP = list(zip(x,y))

    sP= np.array(sP)
    point = np.array(point)

    distances = np.linalg.norm(sP - point, ord=2, axis=1.)

    if relative:
        maxDist = max(distances)
        distances = distances/maxDist

    return distances



def CalcTravelLength(particleList):
    """ Calculate the travel length of every particle in every time step
    in the sl-coordinate system

    Parameters
    ----------
    particleList: list of dict
        Particles dictionary for a specific time step

    Returns
    -------
    particleList: list of dicts
        updated particle dictionary for a specific time step
    """
    for tstep in particleList:
        tstep['sParticle'] = np.array(tstep['sAimec']) - np.array(particleList[0]['sAimec'])

    return particleList



def readDEM(avaDir,demName):
    """ read the ascii DEM file from a provided avalanche directory

    Parameters
    ----------
    avaDir : str
        path to avalanche directory

    Returns
    -------
    dem : dict
        dict with header and raster data
    """
    # get dem file name
    demSource = pathlib.Path(avaDir, 'Inputs', demName)

    dem = IOf.readRaster(demSource)
    return(dem)



def GetPointOnDEM(Point,dem):
    """ Project point on DEM

    Parameters
    ----------
    Point: str
        centerTrackPartPoint of the location of the particles to track
        (x|y coordinates)
    dem: dict
        dem dictionary

    Returns
    ----------
    PointOnDEM: dict
        point dictionary:
            x : x coordinate
            y : y coordinate
            z : z coordinate
    """
    Point = Point.split('|')
    PointOnDEM = {'x': np.array([float(Point[0])]),
                            'y': np.array([float(Point[1])])}

    PointOnDEM, _ = geoTrans.projectOnRaster(dem, PointOnDEM, interp='bilinear')
    PointOnDEM['x'] = (PointOnDEM['x']
                        - dem['header']['xllcenter'])
    PointOnDEM['y'] = (PointOnDEM['y']
                        - dem['header']['yllcenter'])
    return PointOnDEM



def baseround(x, base=5):
    return base * round(x/base)



def readAvaRange(AvaRangeDir,FileNames,avalancheDir,dem,xllcenter,yllcenter):
    """ Read AvaRange result files
    
    Parameters
    ----------
    AvaRangeDir:
        path to the folder where the files are saved
    FileNames: list of str
        list with the file names that should be imported
    avalancheDir: str
        path to avalanche directory
    dem: dict
        dem dictionary (normal information will be added in this function)
    xllcenter,yllcenter: float
        Center of local coordinatesystem

    Return
    ------
    ARList: list of dict
        List of dictionaries with the AvaNode results
    """
    ARList = []
    # tsyncs = [384,401,715,491,1395,460,290,613]         # Index of release time
    # start_time_list = [37.5, 39.16, 71.05, 49-0.25, 137.55, 46-.585, 28+0.1152, 61-0.461]
    tsyncs = [375,391,705,487,1376,454,281,605]         # Index of release time
    # tsynce = [410,410,410,410,410,410,410,410]         # Index of runout time
    # end_time_list = [80+2.55, 82-1.071, 112-5.11, 98-3.94, 189-4.373, 85-4.95, 70-4.82, 102-4.56]
    tsynce = [451,410,360,454,472,347,371,370]         # Index of runout time
    factor = [1000,1000,1000,1000,1000000,1000000,1000000,1000000]

    for i in range(len(FileNames)):
        df = pd.read_csv(pathlib.Path(AvaRangeDir,FileNames[i]))
        df['Time'] = df['Time']/factor[i]                                   # change from micro seconds to seconds
        
        # remove first row before the release of the avalanche
        df.drop(index=df.index[0:tsyncs[i]],inplace = True)
        # Daten auf t=0 synkronisieren
        t0 = df['Time'].min()
        df['Time'] = df['Time']-t0
        # remove last timesteps
        df.drop(labels=df.index[tsynce[i]:],inplace = True)

        # Rename columns
        df.rename(columns = {'Time':'t','East':'x','North':'y','v':'umag'}, inplace = True)

        # add xllcenter and yllcenter to the dataframe
        df['xllcenter'] = xllcenter
        df['yllcenter'] = yllcenter
        # Transform in local coordinate system
        df['x'] = df['x'] - xllcenter
        df['y'] = df['y'] - yllcenter

        data = df.to_dict('records')                                    # list of dict
        ARList.append(data)

        # Transform coordinates into the s/l-coordinate system
        ARList[i] = transformSL(ARList[i],avalancheDir,dem,AvaNode=True)
        
    return ARList



def MaxList(ListOfDict,key):
    """ Calculate list of maximums of a list of dictionary
    
    Parameters
    ----------
    ListOfDict: list of dict
        list of dictionaries
    key: str
        key
    """
    ml = np.array([i[key] for i in ListOfDict])
    
    if isinstance(ml[0],np.ndarray):
        # Finde das absolute Maximum jeder Spalte in `ml`
        abs_max = np.amax(np.abs(ml), axis=0)
        # Setze das Vorzeichen der Werte in jeder Spalte auf das Vorzeichen des absoluten Maximums
        signs = np.sign(ml[np.argmax(np.abs(ml), axis=0), np.arange(ml.shape[1])])
        # Finde das Maximum jedes Spaltenvektors in `signed_max`
        result = abs_max * signs
    else:
        result = max(ml, key=abs)

    return result





