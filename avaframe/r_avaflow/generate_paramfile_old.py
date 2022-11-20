"""
    Generating parameter file (param.txt) for r.avaflow simulations
"""

# Importing libraries
import os
import pathlib

#local imports
import avaframe.r_avaflow as r_avaflow

def generate_param(cfgAvaflow,cfgCom1DFA,avalancheDir):
    """
    Parameters
    ----------
    cfgAvaflow : configparser object
        main configuration of a r.avaflow simulation
    avalancheDir : str or pathlib Path
        path to avalanche data

    Returns
    -------
    Returns the parameter file (param.txt) which is needed for r.avaflow
    """
    
    #Load .ini-files
    cfgGen = cfgAvaflow["GENERAL"]
    cfgCom1DFA = cfgCom1DFA["GENERAL"]
      
    #Setting correct path to the parameter file
    paramdir2 = os.path.dirname(os.getcwd())

    paramdir = paramdir2 + "/avaframe/" + avalancheDir + "/Avaflow_Input/" + "param1.txt"
    print(paramdir)
    
    # start writing parameter file and initializing parameter file 
    p1file = open(paramdir, "w")
   
    # mainmapset = cfgGen["mainmapset"]
    # print(mainmapset, file=p1file)  # name of main mapset
    print("0", file = p1file)
    
    #Multiple model runs  (0 or 1, default: 0)
    mflag = cfgGen["mflag"]
    if mflag == "0":
        print("0", file = p1file) # identifier for single model run
    else:
        print("1", file= p1file)
    print("1", file = p1file) # mesh spacing (multiple of ascii cell size)
    
    #Prefix for output, files and folders
    prefix = cfgGen["prefix"]
    print(prefix, file = p1file)
    
    directory = pathlib.Path(r_avaflow.__file__).resolve().parent
    directory = os.path.join(paramdir2, "avaframe", avalancheDir,"Avaflow_Input/")
    print(directory, file = p1file)
    
    # Contorl: in orig "None"
    #print("%s_results/%s_ascii/" % (prefix, prefix), file=p1file)  # path and prefix for storing output maps
    #print("%s_results/%s_files/" % (prefix, prefix), file=p1file)  # path and prefix for storing output files
    #print("%s_results/%s_paraview/" % (prefix, prefix), file=p1file)  # path and prefix for storing output paraview files
    #print("%s_results/%s_plots/" % (prefix, prefix), file=p1file)  # path and prefix for storing output R plots
    
    print("None", file = p1file)
    print("None", file = p1file)
    print("None", file = p1file)
    print("None", file = p1file)
    print("None", file = p1file)
    
    #Phases - A maximum of three phases can be defined through shortcuts  (default: m)
       ## x = Mixture (Voellmy-type model)
       ## s = Solid (plastic behaviour, frictional, non-viscous)
       ## fs = Fine solid (plasticity-dominated viscoplastic behaviour, frictional and viscous)
       ## f = Fluid (viscosity-dominated viscoplastic, non-frictional, viscous)
       ## m = Multi-phase (P1: solid, P2: fine solid, P3: fluid)
       ## s,s,f = Multi-phase (P1: solid, P2: solid, P3: fluid)
    
    phases = cfgGen["phases"]
    if phases == "m":
        phases = "s,fs,f"
    phases = list(map(str, phases.split(",")))
        
    if len(phases) > 3:
            print("Error: number of phases")
    else:
        if len(phases) == 1:
            model = 1
        elif len(phases) == 2:
            print("Error: number of phases. Please use the multi-phase model (three phases) and provide no release data for the phase which is not needed.")
        elif len(phases) == 3:
            model = 7
    
        phasesnum = []
    #Phasetext - Optionally, comma-separated strings can be provided, describing the materials associated to each phase. e.g. of a multi-phase model would be (phasetext=Rock,Ice,Water)
    phasetext = []

    for i in range(0, len(phases)):
        if not (len(phases) == 1 and phases[i] == "x") and not phases[i] == "s" and not phases[i] == "fs" and not phases[i] == "f":
            print("Error: types of phases")
        elif phases[i] == "x":
            phasesnum.append(1)
            model = 0
            if not phasetext:
                phasetext.append("Mixture")
        elif phases[i] == "s":
            phasesnum.append(1)
            if not phasetext:
                phasetext.append("Solid")
        elif phases[i] == "fs":
            phasesnum.append(2)
            if not phasetext:
                phasetext.append("Fine")
        elif phases[i] == "f":
            phasesnum.append(3)
            if not phasetext:
                phasetext.append("Fluid")

    if model <= 3:
        phasesnum.append(0)
        phasesnum.append(0)
        phasetext.append("NA")
        phasetext.append("NA")
        
    print(str(model), file = p1file)
        
    for l in range(0,3):                           #phasesnum
         print(str(phasesnum[l]), file = p1file)


    #Map plots of pressure and kinetic energy (0 or 1, default: 0)
    aflag = cfgGen["aflag"]      
    print(aflag, file = p1file)

    #Paraview output  (0 or 1, default: 1)
    pflag = cfgGen["pflag"]      
    print(pflag, file = p1file)
     
    #map plots of impact wave or tsunami height  (0 or 1, default: 0)
    tflag = cfgGen["tflag"]      # control for tsunami output 
    print(tflag, file = p1file)
    
    #Gravitational force
    gravity = cfgGen["gravity"]
    print(gravity, file = p1file)
    
    #Numerical limiter (1=Minmod, 2=Superbee, 3=Woodward, 4=van Leer)
    limiter = cfgGen["limiter"]
    print(limiter, file = p1file)
    
    #Layer mode (0=no, 1=weak, 2=strong)
    layers = cfgGen["layers"]
    print(layers, file= p1file)

    controls = cfgGen["controls"]
    
    if controls == "0":
        controls = "0,0,1,0,0,0,0,0,0,2,0"
    controls = controls.split(",")
    print(controls[0], file=p1file)  # control for correction of flow height
    print(controls[1], file=p1file)  # control for diffusion control
    print(controls[2], file=p1file)  # control for curvature
    print(controls[3], file=p1file)  # control for surface control
    print(controls[4], file=p1file)  # control for entrainment and deposition
    print(controls[5], file=p1file)  # control for stopping
    print(controls[6], file=p1file)  # control for dynamic variation of friction
    print(controls[7], file=p1file)  # control for non-hydrostatic effects
    print(controls[8], file=p1file)  # control for phase separation
    print(controls[9], file=p1file)  # control for input hydrograph management      
    print(controls[10], file=p1file)  # control for deceleration management
    
    #Name of elevation ascii-file
    elevation = cfgGen["elevation"]
    print(elevation, file=p1file)  # name of elevation map
    
    #Release mass MIXTURE or PHASE 1 - ascii-file
    hrelease = cfgGen["hrelease"]
    hrelease1 = cfgGen["hrelease1"]
    if hrelease:
        print(hrelease, file=p1file)  # name of MIXTURE release height map
    elif hrelease1 != "0":
        print(hrelease1, file=p1file)  # name of PHASE 1 release height map
    else:
        print("None", file=p1file)
    
    #Release mass PHASE 2
    hrelease2 = cfgGen["hrelease2"]
    if hrelease2 != "0":
        print(hrelease2, file=p1file)  # name of PHASE 2 release height map
    else:
        print("None", file=p1file)
    
    #Release mass PHASE 3
    hrelease3 = cfgGen["hrelease3"]
    if hrelease3 != "0":
        print(hrelease3, file=p1file)  # name of PHASE 3 release height map
    else:
        print("None", file=p1file)
        
    #Fraction of PHASE 1 release heigh
    rhrelease1 = cfgGen["rhrelease1"]
    if rhrelease1 != "0":
        print(rhrelease1, file=p1file)  # fraction of PHASE 1 release height
    else:
        print("1.00", file=p1file)
        
    #Variation of release height
    vhrelease = cfgGen["vhrelease"]
    if vhrelease != "0":
        print(vhrelease, file=p1file) # variation of release height
    else:
        print("1.00", file=p1file)

    #Release velocity x MIXTURE or PHASE 1
    vinx1 = cfgGen["vinx1"]
    if vinx1 != "0":
        print(vinx1, file=p1file)  # name of MIXTURE or PHASE 1 release velocity in x direction map
    else:
        print("None", file=p1file)
       
    #Release velocity x PHASE 2
    vinx2 = cfgGen["vinx2"]   
    if vinx2 != "0":
        print(vinx2, file=p1file)  # name of PHASE 2 release velocity in x direction map
    else:
        print("None", file=p1file)
        
    #Release velocity x PHASE 3
    vinx3 = cfgGen["vinx3"]    
    if vinx3 != "0":
        print(vinx3, file=p1file)  # name of PHASE 3 release velocity in x direction map
    else:
        print("None", file=p1file)
        
    #Release velocity y MIXTURE or PHASE 1
    viny1 = cfgGen["viny1"]     
    if viny1 != "0":
        print(viny1, file=p1file)  # name of MIXTURE or PHASE 1 release velocity in y direction map
    else:
        print("None", file=p1file)
      
    #Release velocity y PHASE 2
    viny2 = cfgGen["viny2"]    
    if viny2 != "0":
        print(viny2, file=p1file)  # name of PHASE 2 release velocity in y direction map
    else:
        print("None", file=p1file)
        
    #Release velocity y PHASE 3
    viny3 = cfgGen["viny3"]
    if viny3 != "0":
        print(viny3, file=p1file)  # name of PHASE 3 release velocity in y direction map
    else:
        print("None", file=p1file)
    
    #Maximum height of entrainment MIXTURE or PHASE 1
    hentrmax = cfgGen["hentrmax"]   
    hentrmax1 = cfgGen["hentrmax1"]
    if hentrmax == 1:
        print(hentrmax, file=p1file)  # name of maximum height of MIXTURE entrainment map
    elif hentrmax1 != "0":
        print(hentrmax1, file=p1file)  # name of maximum height of PHASE 1 entrainment map
    else:
        print("None", file=p1file)
      
    #Maximum height of entrainment PHASE 2
    hentrmax2 = cfgGen["hentrmax2"]
    if hentrmax2 != "0":
        print(hentrmax2, file=p1file) # name of maximum height of PHASE 2 entrainment map
    else:
        print("None", file=p1file)
        
    #Maximum height of entrainment PHASE 3
    hentrmax3 = cfgGen["hentrmax3"]   
    if hentrmax3 != "0":
        print(hentrmax3, file=p1file)  # name of maximum height of PHASE 3 entrainment map
    else:
        print("None", file=p1file)
       
    #Fraction of PHASE 1 maximum height of entrainment
    rhentrmax1 = cfgGen["rhentrmax1"]
    if rhentrmax1 != "0":
        print(rhentrmax1, file=p1file)  # fraction of PHASE 1 maximum height of entrainment
    else:
        print("1.00", file=p1file)
        
    #Variation of maximum height of entrainment
    vhentrmax = cfgGen["vhentrmax"]
    if vhentrmax != "0":
        print("1.00", file=p1file) # variation of maximum height of entrainment
    else:
        print("1.00", file=p1file)
      
    #Zones
    zones = cfgGen["zones"]
    if zones != "0":
       print(zones, file=p1file)  # name of zones map
    else:
        print("None", file=p1file)  

    #Entrainment coefficient
    centr = cfgGen["centr"]                 
    if centr != "0":
        print(centr, file=p1file)  # name of entrainment coefficient map
    else:
        print("None", file=p1file)

    #Shear velocity coefficient
    cvshear = cfgGen["cvshear"]               
    if cvshear != "0":
        print(cvshear, file=p1file)  # name of shear velocity coefficient map
    else:
        print("None", file=p1file) 

    #Internal friction angle of MIXTURE or PHASE 1
    phi1= cfgGen["phi1"]              
    if phi1 != "0":
        print(phi1, file=p1file)  # name of internal friction angle of MIXTURE or PHASE 1 map
    else:
        print("None", file=p1file)  

    #Internal friction angle of PHASE 2
    phi2= cfgGen["phi2"]              
    if phi2 != "0":
        print(phi2, file=p1file)  # name of internal friction angle of PHASE 2 map
    else:
        print("None", file=p1file)    

    #Internal friction angle of PHASE 2
    phi3= cfgGen["phi3"]                
    if phi3 != "0":
        print(phi3, file=p1file)  # name of internal friction angle of PHASE 3 map
    else:
        print("None", file=p1file)  

    #Basal friction difference
    deltab= cfgGen["deltab"]              
    if deltab != "0":
        print(deltab, file=p1file)  # name of basal friction difference map
    else:
        print("None", file=p1file) 
       
    #Turbulent friction coefficient
    tufri= cfgGen["tufri"]            
    if tufri != "0":
        print(tufri, file=p1file)  # name of turbulent friction coefficient map
    else:
        print("None", file=p1file)
    
    #Basal friction angle of MIXTURE or PHASE 1
    delta1= cfgGen["delta1"]
    if delta1 != "0":
        print(delta1, file=p1file)  # name of basal friction angle of MIXTURE or PHASE 1 map
    else:
        print("None", file=p1file)
     
    #Basal friction angle of PHASE 2
    delta2= cfgGen["delta2"]
    if delta2 != "0":
        print(delta2, file=p1file)  # name of basal friction angle of PHASE 2 map
    else:
        print("None", file=p1file)  

    #Basal friction angle of PHASE 3
    delta3= cfgGen["delta3"]                  
    if delta3 != "0":
        print(delta3, file=p1file)  # name of basal friction angle of PHASE 3 map
    else:
        print("None", file=p1file)   

    #Kinematic viscosity of MIXTURE or PHASE 1
    ny1= cfgGen["ny1"]                 
    if ny1 != "0":
        print(ny1, file=p1file)  # name of viscosity of MIXTURE or PHASE 1 map
    else:
        print("None", file=p1file)  
     
    #Kinematic viscosity of PHASE 2
    ny2= cfgGen["ny2"]              
    if ny2 != "0":
        print(ny2, file=p1file)  # name of viscosity of PHASE 2 map
    else:
        print("None", file=p1file)  
       
    #Kinematic viscosity of PHASE 3
    ny3= cfgGen["ny3"]         
    if ny3 != "0":
        print(ny3, file=p1file)  # name of viscosity of PHASE 3 map
    else:
        print("None", file=p1file)  

    #Ambient drag
    ambdrag= cfgGen["ambdrag"]                   
    if ambdrag != "0":
        print(ambdrag, file=p1file)  # name of ambient drag coefficient map
    else:
        print("None", file=p1file)  

    #Fluid friction number
    flufri= cfgGen["flufri"]                   
    if flufri != "0":
        print(flufri, file=p1file)  # name of fluid friction coefficient map
    else:
        print("None", file=p1file)   
    
    #PHASE 1 - PHASE 2 transformation coefficient        
    ctrans12= cfgGen["ctrans12"] 
    if ctrans12 != "0":
        print(ctrans12, file=p1file)  # name of PHASE 1-PHASE 2 transformation coefficient map
    else:
        print("None", file=p1file)  

    #PHASE 1 - PHASE 3 transformation coefficient
    ctrans13= cfgGen["ctrans13"]                   
    if ctrans13 != "0":
        print(ctrans13, file=p1file)  # name of PHASE 1-PHASE 3 transformation coefficient map
    else:
        print("None", file=p1file)  
    
    #PHASE 2 - PHASE 3 transformation coefficient
    ctrans23= cfgGen["ctrans23"]         
    if ctrans23 != "0":
        print(ctrans23, file=p1file)  # name of PHASE 2-PHASE 3 transformation coefficient map
    else:
        print("None", file=p1file)  

    #Release time
    trelease= cfgGen["trelease"]             
    if trelease != "0":
        print(trelease, file=p1file)  # name of release time map
    else:
        print("None", file=p1file)  

    #Release stop time
    trelstop= cfgGen["trelstop"]                   
    if trelstop != "0":
        print(trelstop, file=p1file)  # name of release hydrograph stop time map
    else:
        print("None", file=p1file) 

    #Stopping time
    stoptime= cfgGen["stoptime"]                    
    if stoptime != "0":
        print(stoptime, file=p1file)  # name of stopping time map
    else:
        print("None", file=p1file)  

    #Time of initial sliding
    tslide = cfgGen["tslide"]                   
    if tslide != "0":
        print(tslide, file=p1file)  # name of time of initial sliding map
    else:
        print("None", file=p1file) 

    #Observed impact area
    impactarea = cfgGen["impactarea"]                    
    if impactarea != "0":
        print(impactarea, file=p1file)  # name of observed impact area map
    else:
        print("None", file=p1file)  

    #Observed height of deposit
    hdeposit = cfgGen["hdeposit"]                   
    if hdeposit != "0":
        print(hdeposit, file=p1file)  # name of height of observed deposit map
    else:
        print("None", file=p1file)
     
    #Orthophoto not yet implemented
    #orthophoto = cfgGen["orthophoto"]
    #if orthophoto:
        
        # grass.run_command("r.out.gdal", input=ortho_c1, output=ascpath + ortho_c1 + ".asc", format="AAIGrid", overwrite=True)
        # corrasc(ascpath + ortho_c1)     
        
        # grass.run_command("r.out.gdal", input=ortho_c2, output=ascpath + ortho_c2 + ".asc", format="AAIGrid", overwrite=True)
        # corrasc(ascpath + ortho_c2)
        
        # grass.run_command("r.out.gdal", input=ortho_c3, output=ascpath + ortho_c3 + ".asc", format="AAIGrid", overwrite=True)
        # corrasc(ascpath + ortho_c3)
        
        # print(ortho_c1, file=p1file)  # name of orthophoto channel 1 map
        # print(ortho_c2, file=p1file)  # name of orthophoto channel 2 map
        # print(ortho_c3, file=p1file)  # name of orthophoto channel 3 map
          
    #else:
    print("None", file=p1file)
    print("None", file=p1file)
    print("None", file=p1file)
    
    #Parameters of hydrograph profiles
    hydrocoords = cfgGen["hydrocoords"]
    if hydrocoords != "0":
        print(hydrocoords, file=p1file)  # hydrograph profile parameters
    else:
        print("None", file=p1file)
        
    #Pathes to input hydrograph files
    hydrograph = cfgGen["hydrograph"]
    if hydrograph != "0":
        print(hydrograph, file=p1file)  #path to hydrograph file
    else:
        print("None", file=p1file)
    
    #Path to adaptograph file
    adaptograph = cfgGen["adaptograph"]
    if adaptograph != "0":
        print(adaptograph, file=p1file)  # path to adaptograph file
    else:
        print("None", file=p1file)
    
    #Path to frictiograph file
    frictiograph = cfgGen["frictiograph"]
    if frictiograph != "0":
        print(frictiograph, file=p1file)  # path to frictiograph file
    else:
        print("None", file=p1file)
    
    #Path to transformograph file
    transformograph = cfgGen["transformograph"]
    if transformograph != "0":
        print(transformograph, file=p1file)  # path to transformograph file
    else:
        print("None", file=p1file)
     
     # Flow parameters
    mflag = cfgGen["mflag"] #Multiple model runs  (0 or 1, default: 0)
    density = cfgGen["density"] #Density for PHASE 1,2,3 (kg/m³)
    friction = cfgGen["friction"]#Internal & Basal friction angle of PHASE 1,2,3 (degrees) + Fluid friction number
    basal = cfgGen["basal"] # Log 10 of entrainment coefficient & Stopping criterion
    viscosity = cfgGen["viscosity"]#Log10 of kinematic viscosity of PHASE 1 (m²/s) & Yield strength of PHASE 1 (Pa) and same for PHASE 2 & 3
    transformation = cfgGen["transformation"] #Transformation coefficient PHASE 1 - PHASE 2, PHASE 1 - PHASE 3, PHASE 2 - PHASE 3
    dynfric = cfgGen["dynfric"]
    #Dynfric
       ## Minimum value of internal and basal friction
       ## Kinetic energy coefficient
       ## Phase fraction scaling exponent
       
    special = cfgGen["special"]
    #Special:
       ## Shear velocity coefficient for entrainment and deposition
       ## Basal friction difference for entrainment and deposition
       ## Maximum water content of deposit
       ## Ambient drag coefficient
       ## Virtual mass number
       ## l parameter related to the virtual mass coefficients
       ## n parameter related to the virtual mass coefficients
       ## Reduction factor for virtual mass coefficients
       ## Mass flux parameter for drag (m/s)
       ## Exponent for scaling of the fluid-like drag contributions to flow resistance
       ## Exponent for scaling of drag parameter with solid fraction
       ## Terminal velocity (m/s)
       ## Particle Reynolds number
       ## Exponent for drag
       ## Vertical PHASE 1 velocity distribution (0=no shearing, 3=parabolic profile)
       ## Vertical PHASE 2 velocity distribution (0=no shearing, 3=parabolic profile)
       ## Vertical PHASE 3 velocity distribution (0=no shearing, 3=parabolic profile)
       ## Vertical distribution of PHASE 1 (shape factor: 0=uniform, 3=parabolic)
       ## Vertical distribution of PHASE 2 (shape factor: 0=uniform, 3=parabolic)
       ## Vertical distribution of PHASE 3 (shape factor: 0=uniform, 3=parabolic)
       ## Exponent for mobility of PHASE 2 at interface with PHASE 1
       ## Exponent for mobility of PHASE 3 at interface with PHASE 1
       ## Exponent for mobility of PHASE 3 at interface with PHASE 2
       ## Suitably chosen numerical parameter for regularization
       ## Exponent for scaling of viscosity with fraction of phase (0=no scaling, 1=linear scaling)
       ## Maximum value of PHASE 1 passive earth pressure coefficient
       ## Maximum value of PHASE 2 passive earth pressure coefficient
       ## Maximum value of PHASE 3 passive earth pressure coefficient
       ## Shear factor for phase separation
       ## Multiplication factor for dispersion
       ## Coefficient for constraining dispersion and phase separation
       ## Maximum slope angle (degrees) considered as plane surface
       ## Criterion for maximum flow velocity
    
    
    if mflag == "0":

        if model == 0:

            if density == "0":
                density = "2000"
            if friction == "0":
                friction = "35,20,3.0"
            if basal == "0":
                basal = "-7.0,0.0"
            if special == "0":
                special = "0.05,0.0,1.0,4.0,1.0,200.0"
            if dynfric == "0":
                dynfric = "0.0,-6.0"
            lmax = 14

        elif model == 1:

            if density == "0":
                density = "2700"
            if friction == "0":
                friction = "35,20,0.05"
            if viscosity == "0":
                viscosity = "-9999,-9999"
            if basal == "0":
                basal = "-7.0,0.0"
            if special == "0":
                special = "0.05,0.0,0.0,1,10,1,1.0,4.0,1.0,200.0"
            if dynfric == "0":
                dynfric = "0.0,-6.0"
            lmax = 20

        if model == 7:

            if density == "0":
                density = "2700,1800,1000"
            if friction == "0":
                friction = "35,20,0,0,0,0,0.05"
            if viscosity == "0":
                viscosity = "-9999,-9999,-3.0,-9999,-3.0,0.0"
            if basal == "0":
                basal = "-7.0,0.0"
            if transformation == "0":
                transformation = "0.0,0.0,0.0"
            if special == "0":
                special = "0.05,0.0,0.333,0.0,10,0.12,1,1,1,3,1,0.1,1,1,1,1,1,0,0,0,1,1,1,10,0,1,1,1,0.0,1.0,4.0,1.0,200.0"
            if dynfric == "0":
                dynfric = "0.0,-6.0,0.0"
            lmax = 57 
       
    if model == 0:
        flowparam = density + "," + friction + "," + basal + "," + special + "," + dynfric
    elif model == 1:
        flowparam = density + "," + friction + "," + viscosity + "," + basal + "," + special + "," + dynfric
    else:
        flowparam = density + "," + friction + "," + viscosity + "," + basal + "," + transformation + "," + special + "," + dynfric

    flowparam = list(map(str, flowparam.split(",")))
    
    gt = [float(flowparam[0])]
    for l in range(1, lmax):
        gt.append(float(flowparam[l]))
          
    gt = []
    for l in range(0, lmax):
        gt.append(float(flowparam[l]))
      
    print(lmax, file = p1file)  # number of flow parameters
    for l in range(0, lmax):
        print(round(gt[l], 10), file=p1file)  # flow parameters
    
    #Thresholds:
       ## Threshold for display of flow height (m)
       ## Threshold for display of flow kinetic energy (J)
       ## Threshold for display of flow pressure (Pa)
       ## Threshold flow depth for simulation (m)
    thresholds = cfgGen["thresholds"]
    
    if thresholds == "0":
        thresholds = "0.1,10000,10000,0.001"
    thresholds = thresholds.split(",")
    thresholdsc = thresholds[3]
    
    print(thresholdsc, file=p1file)  # threshold of flow height (for computation)
    print(thresholds[0], file=p1file)  # threshold of flow height (for display)
    print(thresholds[1], file=p1file)  # threshold of flow kinetic energy
    print(thresholds[2], file=p1file)  # threshold of flow pressure
    
    #Time interval for writing output (s) & Time at which to stop simulation (s)
    times = cfgGen["times"]
    if times == "0":
            times = "10,300"
    times = times.split(",")
    tint = times[0]
    tstop = times[1]
    
    print(tint, file=p1file)  # time for writing output
    print(tstop, file=p1file)  # process duration at which to stop
    
    
    # Factor for slow motion
    slomo = cfgGen["slomo"]
    if slomo == "0":
        slomo = "1"

    elif slomo == "min":
        slomo = "60"

    elif slomo == "h":
        slomo = "3600"
        
    elif slomo == "d":
        slomo = "86400"
        
    elif slomo == "w":
        slomo = "604800"

    elif slomo == "mon":
        slomo = "2592000"

    elif slomo == "yr":
        slomo = "31536000"

    elif slomo == "g":
        slomo = "-31536000"
    
    print(slomo, file=p1file)  # factor for slow motion
    
    #Slidepar
       ## Search radius for initial sliding (m)
       ## Weighting exponent for initial sliding
       ## Coefficient for deformation during initial sliding
    slidepar = cfgGen["slidepar"]
    if slidepar == "0":
        slidepar = "0,0,0"
    slidepar = slidepar.split(",")
    print(slidepar[0], file=p1file)  # search radius for initial sliding
    print(slidepar[1], file=p1file)  # exponent for weighting for initial sliding
    print(slidepar[2], file=p1file)  # coefficient for deformation
    
    #CFL criterion & Length of first time step (s)
    cfl = cfgGen["cfl"]
    if cfl == "0":
        cfl = "0.40,0.001"
    cfl = cfl.split(",")
    print(cfl[0], file=p1file)  # cfl criterion
    print(cfl[1], file=p1file)  # maximum length of time step
    
    #Profile vertices (x and y coordinates)
    profile = cfgGen["profile"]
    if profile != "0":
        print(str(profile), file=p1file)  # profile vertices (x and y coordinates)
    else:
        print("None", file=p1file)
    
    #Profile vertices (x and y coordinates)
    ctrlpoints = cfgGen["ctrlpoints"]
    if ctrlpoints != "0":
        print(str(ctrlpoints), file=p1file) # control points (x and y coordinates)
    else:
        print("None", file=p1file)
    
    #Paraview:
       ## Minimum flow height for Paraview visualization (m)
       ## Reference flow height for Paraview visualization (m)
       ## Reference tsunami height for Paraview visualization (m)
       ## Minimum level for flow height contours in Paraview (m)
       ## Maximum level for flow height contours in Paraview (m)
       ## Interval for flow height contours in Paraview (m)
       ## Minimum level for elevation contours in Paraview (m)
       ## Maximum level for elevation contours in Paraview (m)
       ## Interval for elevation contours in Paraview (m)
       ## Weight of red colour for flow visualization in Paraview (neglected for multi-phase model)
       ## Weight of green colour for flow visualization in Paraview (neglected for multi-phase model)
       ## Weight of blue colour for flow visualization in Paraview (neglected for multi-phase model)
       ## Exponent of transparency curve for flow visualization in Paraview
    paraview = cfgGen["paraview"]
    if paraview == "0":
            #paraview = "0.1,5.0,5.0,1,100,2,-11000,9000,100,0.60,0.25,0.15,0.2,1.0,None,None,None"
            paraview = "0.1,5.0,5.0,1,100,2,-11000,9000,100,0.60,0.25,0.15,0.2" #adapted !!!
    paraview = paraview.split(",")
    for i in range(0, len(paraview)):
        print(paraview[i], file=p1file) #paraview parameter
        
    # Exaggeration factor for flow heights in profile plots: (0,1,undefined)
    exaggeration = cfgGen["exaggeration"]
    print(exaggeration, file = p1file)

    print(str("/usr/bin/paraview"), file = p1file)
    print(str("/usr/bin/R"), file=p1file)
    print("None", file=p1file)
     

    p1file.close() # closing parameter file
    
    print("Paramfile generated")
