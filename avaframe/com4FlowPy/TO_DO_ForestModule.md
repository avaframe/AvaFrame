#### Implementation of functionality from avaframe/FlowPy 'forest_detrainmente' branch

###### Implementation steps

* __Anpassen com4FlowPyCfg.ini__
    * ~~Add option to turn forest-interaction on/off~~
    * ~~Add Forest-Friction and Forest-Detrainment Parameters~~
    * Add forest input raster path
        * ~~custom paths (user-defined)~~
        * inside AvaFrame data-structure ()

* __Modifikation in flow_class.py__
    * ~~Übergabe der Forest-Interaction Parameter im Konstruktor~~
        * ~~ein einfachsten wsl. als dict {} oder tuple ()~~
    * Adding the forest parameters to Cell class (constructor)
    * ~~add functions / check implementation:~~
        * ~~def forest_detrainment(self)~~
            * 'no_detrainmnet_effect_zdelta' --> only needs to be computed once (i.e.
               simply add to constructor and option tuple)
        * ~~def  calc_z_delta(self)~~            * 'no_friciton_effect_zdelta' --> only calculate once/not at every funciton call!!!
            * only use forest info for calculation of alpha_calc if forest option is 'True'

* __Modifikation in flowCore.py__
    * ~~Anpassen der 'calculation(args)' Funktion (analog zu Infrastruktur)~~
    * ~~Anpassen der 'run(optTuple)' Funktion (analog zu Infrastruktur)~~
        * ~~Übergabe von forest parameters~~
        * check if this needs modifications in 'initializeProject' or similar???

* __Modifications in com4FlowPy.py__
    * ~~tbd.~~

* __Modifications in runCom4FlowPy.py__
    * ~~tbd.~~

###### Testing results against foreste_detrainment branch in avaframe/FlowPy