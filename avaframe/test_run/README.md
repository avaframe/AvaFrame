## This is an attempt to create the folder structure for an avaframe project and run different modules
---
11.08.2020

#### Initialization
The initialization sequence enables the creation of the avaframe folder structure.
python path/to/Main_avaframe.py --help
for details about what is available.

To initialize a directory:
python path/to/Main_avaframe.py -init pathAvalancheName=/path/to/NameofAvalanche ParameterSet
for example:
python avaframe/test_run/Main_avaframe.py -init pathAvalancheName=/home/matthiastonnel/Documents/github/AvaFrame/avaframe/data/TestAva1 locCfgAB=1

This creates the following folder structure:

    NameofAvalanche/
       Inputs/
          REL/    - release areas
          RES/    - resistance areas
          ENT/    - entrainment areas
          LINES/  - avalanche path line
          POINTS/ - split points
          .asc    - DEM
          avalanche_path.xyz
       Outputs/

       Work/


#### Running ALPHABETA
Move all the required input files to Inputs and run:
python path/to/Main_avaframe.py -alphabeta pathAvalancheName=/path/to/NameofAvalanche
for example:
python avaframe/test_run/Main_avaframe.py -alphabeta pathAvalancheName=/home/matthiastonnel/Documents/github/AvaFrame/avaframe/data/TestAva/
