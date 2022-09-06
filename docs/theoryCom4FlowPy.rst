com4FlowPy Theory
=================

![Image](img/Motivation_2d.png)

*Fig. 1: Definition of angles and geometric measures for the calculation of Z<sup>&delta;</sup>, where s is the projected distance along the path and z(s) the corresponding altitude.*

The model equations that determine the run out in three dimensional terrain are mainly motivated with respect to simple, geometric, two dimensional concepts [0,3,4] in conjunction with ideas existing algorithms for flow routing in three dimensional terrain [1,2], controlling the main routing and final stopping of the flow.

Figure 1 summarizes the basic concept of a constant run out angle (alpha) with the corresponding geometric relations in two dimensions along a possible process path.

![tan_alpha](img/tan_alpha.png)

The local travel angle gamma is defined by the altitude difference and projected distance along the path, from the release point to the current location.

![tan_gamma](img/tan_gamma.png)

The angle delta is the difference between the local travel angle gamma and the runout angle alpha and is related to Z<sup>&delta;</sup>, so when delta equals zero or gamma equals alpha, the maximum runout distance is reached.

![z_alpha](img/z_alpha.png)

Z<sup>&alpha;</sup> can be interpreted as dissipation energy.

![z_gamma](img/z_gamma.png)

Z<sup>&gamma;</sup> is the altitude difference between the starting point and the current calculation step at the projected distance s.

Z<sup>&delta;</sup> is the difference between Z<sup>&gamma;</sup> and Z<sup>&alpha;</sup>, so when Z<sup>&delta;</sup> is lower or equal zero the stopping criterion is met and the flow stops. Z<sup>&delta;</sup> is associated to the process magnitude and can be interpreted as the kinetic energy or velocity of the process.

![z_delta](img/z_delta.png)

The major drawback of implementing the geometric runout angle concepts is that they require a predefined flow path in two dimensional terrain. To allow for an enhanced routing in three dimensional terrain without prior knowledge of the flow path we combine these concepts [4] with extensions of existing algorithms [1,2,3] that are described in the following sections.

## Spatial Input and Iterative Calculation Steps on the Path

In nature a GMF has one or more release areas that span over single or multiple release cells. Flow-Py computes the so called path, which is defined as the spatial extent of the routing from each release cell. Each release area (single raster cell in release area  layer) has it's own unique path (collection of raster cells), and a location on the terrain (a single raster cell) can belong to many paths. Flow-Py identifies the path with spatial iterations starting with a release area raster cell and only iterating over cells which receive routing flux.  The corresponding functions are implemented in the code in the flow_class.calc_distribution() function.

To route on the surface of the three dimensional terrain, operating on a quadrilateral grid, we implement the geometric concepts that have been sketched in the model motivation utilizing the following cell definitions:

​	![grid_overview](img/Neighbours.png)

*Fig. 2: Definition of parent, base, child and neighbors, as well as the indexing around the base.*

Each path calculation starts with a release cell and operates on the raster, requiring the definition of parent, base, child and neighbor cells (see Fig. 2).
The base cell is the cell being calculated on the current spatial iteration step. The 8 raster cells surrounding the base cell are called neighbor cells (n, i) which have the potential to be parents (supplying flux to base cell), or a child (receive flux from the base cell). 
In 2d the base cell corresponds to the cell/location at the distance s along the path in Fig. 1.

Every base has at least one parent cell, except in the first calculation step from the release cell, where we start our calculation, this would be at s = s<sub>0</sub> in Fig. 1.

During an iteration step a raster cell from the iteration list is identified as the current base cell. The routing flux is calculated across the base cell from the parent cell to possible child cells. The goal is to keep the spatial iteration steps to a minimum, which is achieved by only adding neighbor cells to the iteration list that have flux routed to them from the base cell and do not meet either of the stopping conditions. These cells are called child cells. Child cells that are not already on the iteration list are added to the list and flow_class python object is created for the raster cell. The child cells flow_class has the parent added to it as a source for routing flux. By being added to the iteration list the cell has been recognized as being part of the GMF path and will be the base cell for a future iteration step.

When the iteration list is empty and all potential children fulfill one of the stopping criteria:

- Z<sup>&delta;</sup> has to be smaller than zero: Z<sup>&delta;</sup> < 0,
- Routing Flux has to be smaller than the flux cut off: R<sub>i</sub> < R<sub>Stop</sub>,

the path calculation is finished. The required information is saved from the cell class to the summarizing output raster files. Then the calculation starts again for the next release cell and respective flow path. The spatial extent and magnitude for all release cells are summarized in the output raster files, which represent the overlay of all paths.

Every path is independent from the other, but depending on the information we want to extract, we save the highest values (e.g. Z<sup>&delta;</sup>) or sums (e.g.Cell Counts) of different paths to the output raster file.

### Z<sup>&delta;</sup>

For each base cell in a path we solve the equations (6,7 and 8) for every neighbor n, if Z<sub>bn</sub><sup>&delta;</sup> is higher than zero, this neighbor is defined as a potential child of this base, and routing  in this direction is possible.

![z_delta_i](img/z_delta_array.png)

Here S<sub>bn</sub> is the projected distance between the base and the neighbor.

As Z<sub>bn</sub><sup>&delta;</sup> can be interpreted as process magnitude (and kinetic energy or velocity respectively) it is possible to limit this value to a maximum. In comparison to process based modeling approaches this would correspond to maximum velocity induced by a velocity dependent turbulent friction term.

![z_delta_max](img/z_delta_max.png)

The projected path lengths, or total travel distance to one of the neighbors (S<sub>n</sub>) equals the path length to the base (S<sub>b</sub>) plus the path from base to the neighbor (S<sub>bn</sub>), which reads:

![S_bn](img/S_bn.png)

As there are many possibilities for the path from the starting point to the actual cell or base, the shortest path is taken into account, corresponding to the highest Z<sup>&delta;</sup> in the base. If Z<sup>&delta;</sup><sub>max</sub> is set to infinity, or as in the code to 8848 m (= Mount Everest), we can calculate the shortest path from the starting point to the base and yields the total projected travel distance:

![S_n_eq1](img/S_n_eq1.png)

This equations determine the routing and corresponding run out distance for the process, the next steps demonstrate how spreading is handled on the surface of the three dimensional terrain.

### Persistence based routing

The persistence contribution P<sub>i</sub> aims to reproduce the behavior of inertia, and takes the change in flow direction into account [3].
The direction contribution is scaled with the process magnitude Z<sup>&delta;</sup><sub>parent</sub>, such that the direction from a parent cell with higher process magnitude has more effect on the path routing and direction.

![](img/persistence.png)

The direction contributions D<sub>n</sub> are defined by the cosine of the angle between parent, base and child/neighbor minus pi:  

![](img/persistence_cos_function.png)

Therefore the direction contribution limits the maximum number of potential children to three, getting input via the persistence function from one parent.

In the first calculation step, at the release or start cell no parent cells are defined and the persistence is set to one. So the first calculation step is solely determined by the terrain contribution.

### Terrain based routing

The terrain based routing is solely dependent on the slope angle phi. The exponent exp allows to control the divergence of the spreading. 
The Holmgren (1994) algorithm [1] is used in different kind of models and works well for avalanches but also rockfall or soil slides. For avalanches an exponent of 8 shows good results. To reach a single flow in step terrain (rockfall, soil slides, steepest descend), an exponent of 75 is considered.

![Holmgrem](img/flow_direction.png)
*Holmgrem Algorithm from 1994 [1]*

To overcome the challenge of routing in flat or uphill terrain, we adapted the slope angle phi for the normalized terrain contribution to:

![Phi_Formula](img/Phi.png)

### Routing Flux 

The routing flux summarizes the persistence and terrain contributions according to Eq.(16):

![](img/flux.png)

where i is the direction and n are the neighbors from 1 to 8. R<sub>i</sub> is then the routing flux in direction i.
R<sub>b</sub> is the flux in the base, for a release cell or starting cell the flux of the base equals one. The result of Eq. (16) is a 3 x 3 array with assigned flux values. A normalization stage is then required to bring the sum of the R<sub>i</sub>'s to the value of R<sub>b</sub>. This aims at avoiding loss of flux [2].

### Flow Chart / Overview

In Fig. 3 the algorithm of the computational implementation is sketched, including function and files names with respect to the code in the repository.

The file main.py handles the input for the computation and splits the release layer in tiles and saves them in a release list. Then the main.py starts one process per tile, which calls the flow_core.py and starts the calculation for one release cell and the corresponding path. The number of processes is depending on the hardware setting (CPU and RAM).  Whenever a new cell is created flow_core.py calls flow_class.py and makes a new instance of this class, which is saved in the path. When the calculation in flow_core.py is finished it returns the path to main.py which saves the result to the output rasters. 

![Flow_Chart](img/Flow-Py_chart.png)

*Fig.3: Flow chart of the Flow-Py computational process and an overview of the files and what they manage.*

### References

[0] [Heim, A. (1932).]: Bergstürze und Menschenleben

[1] [Holmgren, P. (1994).](https://www.researchgate.net/publication/229484151_Multiple_flow_direction_algorithms_for_runoff_modelling_in_grid_based_elevation_models_An_empirical_evaluation) 
Multiple flow direction algorithms for runoff modelling in
grid based elevation models: an empirical evaluation. Hydrological Processes, 8:327–334.

[2] [Horton, P., Jaboyedoff, M.,
Rudaz, B., and Zimmermann, M. (2013).](https://nhess.copernicus.org/articles/13/869/2013/nhess-13-869-2013.pdf) 
Flow-R, a model for susceptibility mapping of debris
flows and other gravitational hazards at a regional scale. Natural Hazards and Earth System
Science, 13:869–885.

[3] [Gamma, P. (1999).](https://www.researchgate.net/publication/34432465_dfwalk-Ein_Murgang-Simulationsprogramm_zur_Gefahrenzonierung) dfwalk - Ein
Murgang-Simulationsprogramm zur Gefahrenzonierung. PhD thesis, Universität Bern.

[4] [Huber, A., Fischer, J. T., Kofler, A., and Kleemayr, K. (2016).] Using spatially
distributed statistical models for avalanche runout estimation. In International Snow Science Workshop, Breckenridge, Colorado, USA - 2016.  

## Contact and acknowledgment

For Questions contact:  
Michael Neuhauser, Austrian Research Centre for Forest: Michael.Neuhauser@bfw.gv.at  
Christopher D'Amboise, Austrian Research Centre for Forest: Christopher.DAmboise@bfw.gv.at

This study was carried out in the framework of the GreenRisk4Alps project
ASP635, funded by the European Regional Development Fund through the Interreg Alpine Space programme. Additional financial support from the AvaRange (www.AvaRange.org, international cooperation project “AvaRange - Particle Tracking in Snow Avalanches” supported by
the German Research Foundation (DFG) and the Austrian Science Fund (FWF, project number I 4274-N29) and the AvaFrame (www.AvaFrame.org, AvaFrame - The open Avalanche Framework is a cooperation between the Austrian Research Centre for Forests (Bundesforschungszentrum für Wald; BFW) and Austrian Avalanche and Torrent Service (Wildbach- und Lawinenverbauung; WLV) in conjunction with the Federal Ministry Republic of Austria: Agriculture, Regions and Tourism (BMLRT)) projects are greatly acknowledged.

## Citation:

Michael Neuhauser, Christopher D'Amboise, Michaela Teich, Andreas Kofler, Andreas Huber, Reinhard Fromm, & Jan Thomas Fischer. (2021, June 24). Flow-Py: routing and stopping of gravitational mass flows (Version 1.0). Zenodo. http://doi.org/10.5281/zenodo.5027275
