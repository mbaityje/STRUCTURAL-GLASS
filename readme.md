# STRUCTURAL-GLASS  
# Molecular Dynamics of Supercooled liquids

A repository for MD simulations of supercooled liquids based on [HOOMD](http://glotzerlab.engin.umich.edu/hoomd-blue/).

This repository contains code aimed at the development of the following projects:

- **Generic Tutorials** I spent some time making some tutorials for hoomd, since I couldn't find any useful examples on the web.

- **Noise Correlation Functions**
In the Nishimori-Swanzig formalism of the projection operator technique, calculate the correlation functions of the dynamics in the space orthogonal to the relevant variables.

- **Metabasin Dynamics**
Follow the succession of the inherent structures corresponding to every single configuration along the dynamics, identifying metabasins and their mutual relation.



---
## Table Of Contents

- Installation
- Organization of the repository
- Projects
    - Noise Correlation
    - Metabasin Dynamics

---

## Installation
No installation is required.

---

## Organization of the repository
The repository contains the following directories:

#### ./TUTORIALS/
I made a series of tutorials. They are self-consistent well-commented pieces of code that face several tasks that are implemented in the production code.
They rely on two directories:  
-``./TUTORIALS/sample-states/``: contains some initial configurations.   
-``./TUTORIALS/test-output/``: contains the output of the tutorials.  

#### ./OUTPUT/  
This directory is not in the repository for storage reasons. However, it is needed, as the data are assumed to be in there.

#### ./THERMALIZE/  
This is where currently all the code is stored. There are two subfolders

- `./THERMALIZE/progs/`:  
All the programs.


- `./THERMALIZE/script/`:  
All the scripts.

#### ./LIB/  
In this directory I store the few modules that I made for the MD simulations:  

- `./LIB/module_measurements.py/`:  
Functions for measurements on the packing.

- `./LIB/module_potentials.py/`  
Functions for setting the potential.

- `./LIB/module_signals.py/`  
Functions for dealing with external kill signals.

- `./LIB/module_timelists.py/`  
Functions for dealing with times.

- `./LIB/module_debug.py/`:  
Functions for debugging.

- Other:  
Currently, there is some code taken from [Ian Dun](https://github.com/iansdunn) for fitting functions as a sum of exponentials. I haven't decided whether it is useful for me or not.


#### ./PLAYGROUND/  
This directory is where I test new ideas of code. It is meant to be the disordered place where things are born.

#### ./PLOTS/  
Directory devoted to making plots. Not all the plot-making programs are here (in some cases I like to end the simulation with a graph.png), but if a program is purely devoted to plots, then it will be here.

#### ./PAPERS/  
A directory containing some research papers.

#### ./UTILITIES/
Some scripts for doing useful straightforward things, such as extracting instant information from gsd files. Program names are self-explanatory.

#### ./FRANCOIS/  
Some code written by François Landes. Better versions are available in his github repository [his git nickname and repository].

---

# Projects

## Noise Correlations
System size is *N* = 1080.

Temperatures are *T* = **10.0**, **2.0**, (1.0), (0.8), **0.6**, (0.53), (0.49), (0.466), (0.45), (0.44), (0.43).

10 samples per temperature.

### Generating Thermalized Configurations

### Making sure the configurations are well-thermalized

### Running new NVE trajectories for correlation functions

### Calculating Average Trivial Correlation Functions

### Calculating Noise Correlation Functions

---


## Metabasin Dynamics
System size is *N* = 65.

Temperatures are *T* = **10.0**, **2.0**, (1.0), (0.8), **0.6**, (0.53), (0.49), (0.466), (0.45), (0.44), (0.43).

10 samples per temperature.

### Generating Thermalized Configurations

### Making sure the configurations are well-thermalized

### Running new NVE trajectories with inherent-structure sampling, and calculating barriers with the Ridge Method

### Calculating barriers with the Nudged Elastic Band method


---

