# STRUCTURAL-GLASS  
# Molecular Dynamics of Supercooled Liquids

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

---

# Projects

The code is written to be run either on my personal computer, or on a computing cluster (with the SLURM scheduler). This is determined through the `SYSTEM` flag in each script, which can be either `PennPuter` (no scheduler, run online) or `talapas` (with scheduler, run with SLURM-managed queues, option of using GPUs).

Here, I will describe how to manage the main projects, referring to some parameter choices that I made.

---
## Noise Correlations
System size is *N* = 1080.

Temperatures are *T* = 10.0, <s>5.0</s>, 2.0, 1.0, <s>0.8</s>, 0.6, <s>0.52, 0.49, 0.466, 0.45, 0.44, 0.43</s>.

<s>10</s> 5 samples per temperature.

### Generating Thermalized Configurations
```
# Access script directory
cd ./THERMALIZE/script/
# Edit script to change simulation parameters, such as temperature (TLIST) and number of samples (nsam)
emacs ThermalizeN1080.sh
# Launch thermalization script
bash ThermalizeN1080.sh
cd -
```

The script `ThermalizeN1080.sh` reads the file `./THERMALIZE/data/thermalizationtimes.txt`, which is a table of estimated thermalization times &tau;<sub>est</sub> (from previous runs and from literature arXiv:1203.3392 and arXiv:0805.3104), translates it into number of MD steps, and multiplies it by 10.
Then, it invokes the program `ReadAndThermalize.py` which takes care of running the remaining amount of steps.

### Making sure the configurations are well-thermalized
```
# Access script directory
cd ./THERMALIZE/script/
# Edit script to change simulation parameters, such as temperature (TLIST), number of samples (nsam), system size (N), and more
emacs CheckThermalizationAll.sh
# Calculate trajectories and self-intermediate scattering functions for each sample
bash CheckThermalizationAll.sh
# Make sure that the script MediasFkt.sh has the right parameters
emacs MediasFkt.sh
# Calculate average self-intermediate scattering functions
bash MediasFkt.sh
cd -
# Plot self-intermediate scattering functions
cd ./PLOTS/
emacs Fkt.gp # Make sure that the parameters are the correct ones
gnuplot Fkt.gp
```
The script `CheckThermalizationAll.sh` simply loops across parameter choices.
The script that checks the thermalization is `SelfIntermediateScatteringFunction.sh`, which 

- runs an NVE trajectory (labeled _ifr0) long &tau;<sub>est</sub>;

- calculates the self intermediate scattering function on the _ifr0 trajectory;

- runs 20 &tau;<sub>est</sub> time steps;

- runs another NVE trajectory (labeled _aftergap) long &tau;<sub>est</sub>;

- calculates the self intermediate scattering function on the _aftergap trajectory;

After this, averages are calculated with `MediasFkt.sh`, which puts the averages just outside the samples' directory. The functions can be plotted with `Fkt.gp`, which in the current state does not save the figures (because I don't want them), so it should be used with the `gnuplot` interactive shell.


### Running new NVE (or NVT) trajectories for correlation functions
Now that we have some well-thermalized configurations, we want to create the trajectories on which to perform measurements. In this case, the parameters can be both hard-coded, or inserted through command line.
The syntax is

`bash CreateTrajectories.sh [#-of-trajectories-per-sample] ["list of temperatures separated 
by spaces"] ["sizes separated by spaces"] ["samples separated by spaces"]`

so, for example, if I want to create 15 trajectories at *T*=1.0,0.49 for *N*=1080, for samples 0 through 4, I need to do

`bash CreateTrajectories.sh 15 "1.0 0.49" "1080" "0 1 2 3 4"`

The trajectories of each sample can be found in a subdirectory called `trajectories/`.

### Calculating Average Trivial Correlation Functions

#### Calculate Diagonal Correlations

#### Calculate force and momentum correlations

#### Calculate Self-intermediate scattering function

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

