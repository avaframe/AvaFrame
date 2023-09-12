"""
Run script for tracing and plotting particles

@author: Höller JF
"""

# Clear terminal
import os
clear = lambda: os.system('cls')
clear()


# Load modules
import pathlib                                      # Work with smart paths
import matplotlib                                   # Log Scale, access colormap
import matplotlib.pyplot as plt                     # Plot created plots
import csv
#mpl.use('pgf')

# Temporary modules
# import time                                         # Stop Time

# Local imports
from avaframe.in3Utils import cfgUtils              # Load GeneralConfiguration
from avaframe.in3Utils import logUtils              # Tools for writing log-files
from avaframe.com1DFA import particleTools          # Tools for reading particles from a Python Pikle
import PlotParticles_tools_JFH as ppt               # Further Ploting tools
import PlotParticles_auxiliary_PF_JFH as ppap       # Auxiliary functions for plotting
import PlotParticles_auxiliary_MF_JFH as ppam       # Auxiliary functions for plotting

# Create PDFs from PGFs
from matplotlib.backends.backend_pgf import FigureCanvasPgf
matplotlib.backend_bases.register_backend('pdf', FigureCanvasPgf)





####   ----------------------------------------------------------------------------   ###
#                                                                                       #
#                                     Configurations                                    #
#                                                                                       #
####   ----------------------------------------------------------------------------   ###



Project = 6
# Load AvaRange reuslts
AvaRange = False
match Project:
    case 1:     # avaParabola
        avalancheDir = r'C:\Users\frede\Documents\Programmdateien\AvaFrame\Simulationen_Masterarbeit\avaParabola'
        centerTrackPartPoint = '0|0'
    case 2:     # avaHelix
        avalancheDir = r'C:\Users\frede\Documents\Programmdateien\AvaFrame\Simulationen_Masterarbeit\avaHelix'
        centerTrackPartPoint = '2333|-4872'
    case 3:     # avaHockeyChannel      # evtl. streichen
        avalancheDir = r'C:\Users\frede\Documents\Programmdateien\AvaFrame\Simulationen_Masterarbeit\avaHockeyChannel'
        centerTrackPartPoint = '1350|-4000'
    case 4:     # avaWog
        avalancheDir = r'C:\Users\frede\Documents\Programmdateien\AvaFrame\Simulationen_Masterarbeit\avaWog'
        centerTrackPartPoint = '169300|362450'
    case 5:     # avaNordkette
        avalancheDir = r'C:\Users\frede\Documents\Programmdateien\AvaFrame\Simulationen_Masterarbeit\Vergleich_Nordkette\avaNordkette'
        centerTrackPartPoint = '79400|242000'
        # Load AvaRange reuslts
        AvaRange = True
    case 6:     # avaNordkette
        avalancheDir = r'C:\Users\frede\Documents\Programmdateien\AvaFrame\Simulationen_Masterarbeit\Vergleich_Nordkette\avaNordketteVU'
        # avalancheDir = r'C:\Users\frede\Desktop\avaNordketteVU'
        centerTrackPartPoint = '79400|242000'
        # Load AvaRange reuslts
        AvaRange = True

# Path of the AvaRange results
AvaRangeDir = r'C:\Users\frede\Documents\Programmdateien\AvaFrame\Simulationen_Masterarbeit\Vergleich_Nordkette\AvaNode_data\Export'


## Choose relevant particles
TrackMethod = 1
# Method 1: Koordinates + Radius
radius = 5
# Method 2: Explicit particle choice
PartID = [1,2,3,12718, 3831, 12717, 12999]


# Do the particle files contain the transformed coordinates?
Aimec = True
# Select tracked properties
properties = ['x','y','z','umag','lAimec','sAimec','dAimec']


# Plot Color
figcolor = 'black'
# Decide whether the graphics should be exported for use in Latex
LatexExport = True


# Decide which plots should be generated
plots = [5]

# Select the property to plot
# prop = 'umag'
prop = 'sAimec'
# prop = 'dsCoM'





####   ----------------------------------------------------------------------------   ###
#                                                                                       #
#                                      Run Script                                       #
#                                                                                       #
####   ----------------------------------------------------------------------------   ###



# ---------------------------------------------------------------------------- #
# Step 1: Define the general settings
# ---------------------------------------------------------------------------- #

# Load avalanche directory from general configuration file
# cfgMain = cfgUtils.getGeneralConfig()
# avalancheDir = cfgMain['MAIN']['avalancheDir']
FileName = (pathlib.PurePath(avalancheDir)).name

# log file name; leave empty to use default runLog.log
logName = 'runPlotParticles'
# Start logging
log = logUtils.initiateLogger(avalancheDir, logName)
log.info('Plot Particles')
log.info('Current avalanche: %s', avalancheDir)


# General plot settings
ppap.GenPlotSettings(figcolor)
ppap.PlotForLatex()

# Path to the saved input files (Python Dictionary)
if Aimec:
    inDirPart = pathlib.Path(avalancheDir, 'Outputs', 'Aimec', 'ParticlesProjected', 'particles')
else:
    inDirPart = pathlib.Path(avalancheDir, 'Outputs', 'com1DFA', 'particles')



# ---------------------------------------------------------------------------- #
# Step 2: Scan configuration file and DEM
# ---------------------------------------------------------------------------- #

# Which module is used for the calculation of the avalanche simulation?
comModule='com1DFA'
# Read configuration file
simDF = cfgUtils.createConfigurationInfo(avalancheDir, comModule)
# Get simulation names
SimNames = simDF['simName'].tolist()

# Select simulation number (if multiple Simulations are saved in the folder)
SimNumber = 0
# Read DEM
demName=simDF['DEM']
dem = ppam.readDEM(avalancheDir,demName[SimNumber])

# Einlesen aller Python Dictionaries
fullparticlesList, timeStepInfo = particleTools.readPartFromPickle(inDirPart,SimNames[SimNumber])


# Export List of dict to csv
def export_to_csv(data_list, filename_template, columns=None, path=''):
    for i, data_dict in enumerate(data_list):
        #filename = filename_template.format(i)
        filename = f"{filename_template}_{i}.csv"
        with open(f"{path}/{filename}", 'w', newline='') as file:
            writer = csv.writer(file)

            # Schreiben Sie die ausgewählten Spaltenüberschriften
            if columns:
                writer.writerow(columns)
            else:
                writer.writerow(data_dict.keys())

            # Schreiben Sie die Datenzeilen für die ausgewählten Spalten
            if columns:
                selected_values = [data_dict[col] for col in columns]
                for row in zip(*selected_values):
                    writer.writerow(row)
            else:
                for row in zip(*data_dict.values()):
                    writer.writerow(row)

exprotcsv=False
if exprotcsv:
    # Pfad und Dateiname für die CSV-Datei
    file_path = r'C:\Users\frede\Desktop\avaNordketteVU\Outputs\exported_fullparticlesList'

    # Extrahieren der Feldnamen aus dem ersten Dictionary
    keys = ['x','y','z','lAimec','sAimec']

    export_to_csv(fullparticlesList, 'mydata', columns=keys,path=file_path)



# ---------------------------------------------------------------------------- #
# Step 3: Calculation of additional result variables
# ---------------------------------------------------------------------------- #

# Calculate velocityMagnitude
if 'umag' in properties:
    fullparticlesList = ppam.velocityMagnitude(fullparticlesList)

# Calculate centre of mass in every time step
fullparticlesList = ppam.CentreOfMass(fullparticlesList)
# Transform xy-coordinates into sl-coordinates
fullparticlesList = ppam.transformSL(fullparticlesList,avalancheDir,dem)
# Calculate Distance to CoM in every time step
fullparticlesList = ppam.CalcDistCoM(fullparticlesList)
# Calculate travel length of every particle
fullparticlesList = ppam.CalcTravelLength(fullparticlesList)



# ---------------------------------------------------------------------------- #
# Step 4: Load AvaRange Results
# ---------------------------------------------------------------------------- #

if AvaRange:
    # Read all available filenames
    AvaRangeFiles = os.listdir(AvaRangeDir)
    propsAR = ['Index','Timestep','North','East','z']       # Wenn möglich entfernen


    AvaNodesList = ppam.readAvaRange(AvaRangeDir,AvaRangeFiles,avalancheDir,dem,
                                fullparticlesList[0]['xllcenter'],fullparticlesList[0]['yllcenter'])


    for particles in AvaNodesList:
        # Add CoM
        ppam.ARCentreOfMass(particles,fullparticlesList)
        # Calculate distance to CoM
        particles = ppam.CalcDistCoM(particles,Euclidean=False)

    plotdsCoM = False
    if plotdsCoM:
        ppt.plotdsCoM(AvaNodesList,ProjectName=FileName,AvaDir=avalancheDir)




# ---------------------------------------------------------------------------- #
# Step 5: Get relevant Particle IDs for particle tracing
# ---------------------------------------------------------------------------- #

match TrackMethod:
    case 1:
        if AvaRange:
            particles2TrackList = []
            for i in range(len(AvaNodesList)):
                # Search for particles around special coordinates
                PointOnDEM = {'x': AvaNodesList[i][0]['x'],
                                'y': AvaNodesList[i][0]['y']}
                particles2Track = particleTools.findParticles2Track(fullparticlesList[0], PointOnDEM, radius)[0]
                particles2TrackList.append(particles2Track)
        else:
            PointOnDEM = ppam.GetPointOnDEM(centerTrackPartPoint,dem)
            particles2Track = particleTools.findParticles2Track(fullparticlesList[0], PointOnDEM, radius)[0]
    case 2:
        particles2Track = PartID
    case other:
        log.error('Please select a valid particle tracing method!')
        exit()


    

# ---------------------------------------------------------------------------- #
# Step 6: Create a dictionary of the tracked particles and their properties
# ---------------------------------------------------------------------------- #

# Create index list of tracked particles (for each time step in the Dictonary):
particlesList, nPartTracked = particleTools.getTrackedParticles(fullparticlesList, particles2Track)

# Reduce the dictionary to the tracked particles and the desired properties
trackedPartProp = particleTools.getTrackedParticlesProperties(particlesList, nPartTracked, properties)

# Add centre of mass
trackedPartProp['lCoM'] = [None]*len(particlesList)
trackedPartProp['sCoM'] = [None]*len(particlesList)
for p in particlesList:
    trackedPartProp['lCoM'].append(p['lCoM'])
    trackedPartProp['sCoM'].append(p['sCoM'])






# ---------------------------------------------------------------------------- #
# Step 7: Create plots
# ---------------------------------------------------------------------------- #
# if AvaRange:
#     ppt.plotdsCoMnumeric(fullparticlesList,ProjectName=FileName,AvaDir=avalancheDir,AvaNodes=AvaNodesList)
# else:
#     ppt.plotdsCoMnumeric(fullparticlesList,ProjectName=FileName,AvaDir=avalancheDir)


for p in plots:
    match p:     
        ##########     Auswertung der Geschwindigkeitsgrößen     ##########
        case 4:         # 4. Geschwindigkeitsdichte (Vergleich mit AvaNodes -> Statistische Auswertung, wie groß ist die Abweichung der Nodes zu der Verteilung der Simulation)
            # Relativer Abstand zum Schwerpunkt (Trajektorale Entfernung) -> Dichteplot im s/l-Koordinatensystem
            #ppt.CompareHistogram2d(fullparticlesList,'umag',LatexExport=LatexExport,ProjectName=FileName,AvaDir=avalancheDir)
            # ppt.CompareHistogram2d(fullparticlesList,'umag',SimName=SimNames[SimNumber],ProjectName=FileName)
            if AvaRange:
                ppt.CompareHistogram2d(fullparticlesList,'umag',LatexExport=LatexExport,ProjectName=FileName,AvaDir=avalancheDir,AvaNodes=AvaNodesList)
            else:
                ppt.CompareHistogram2d(fullparticlesList,'umag',LatexExport=LatexExport,plotStats=True,ProjectName=FileName,AvaDir=avalancheDir)

        case 5:         # 5. Ergebnisse im Anbruch- und Ablagerungsgebiet gegenüber stellen
            if AvaRange:
                ppt.CompareTrajectory(fullparticlesList,property=prop,ProjectName=FileName,AvaDir=avalancheDir,LatexExport=LatexExport,AvaNodes=AvaNodesList,radius=radius)
            else:
                ppt.CompareTrajectory(fullparticlesList,property=prop,ProjectName=FileName,AvaDir=avalancheDir,LatexExport = LatexExport)

        case 6:         # 6. Sankey Diagramm erstellen, zur Visualisierung der Partikelströme
            if AvaRange:
                if prop=='dsCoM':
                    ppt.ThalwegTime(fullparticlesList,property=prop,colorscheme='explicit',LatexExport=False,AvaDir=avalancheDir,ProjectName=FileName,AvaNodes=AvaNodesList)
                else:
                    ppt.SankeyChart(fullparticlesList,property=prop,colorscheme='explicit',AvaNodes=AvaNodesList,LatexExport=False,AvaDir=avalancheDir,ProjectName=FileName)
            else:
                match prop:
                    case 'umag':
                        ppt.SankeyChart(fullparticlesList,property=prop,colorscheme='explicit',LatexExport=False,AvaDir=avalancheDir,ProjectName=FileName)
                    case 'sAimec':
                        ppt.SankeyChart(fullparticlesList,property=prop,colorscheme='explicit',LatexExport=False,AvaDir=avalancheDir,ProjectName=FileName)
                    case 'dsCoM':
                        # ppt.SankeyChart(fullparticlesList,property=prop,colorscheme='explicit',LatexExport=False,AvaDir=avalancheDir,ProjectName=FileName)
                        ppt.ThalwegTime(fullparticlesList,property=prop,colorscheme='explicit',LatexExport=False,AvaDir=avalancheDir,ProjectName=FileName)

        case 7:         # Stelle AvaNodes im Anburchgebiet dar
            if AvaRange:
                ppt.plotParticles(fullparticlesList,AvaNodesList,ProjectName=FileName,radius=radius,AvaDir=avalancheDir)        #dem=dem,




# Show plots
plt.show(block=True)
