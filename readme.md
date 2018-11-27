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
I made a series of tutorials. They are self-consistent thoroughly commented pieces of code that face several tasks that are implemented in the production code.
They rely on two directories:  
-``./TUTORIALS/sample-states/``: contains some initial configurations.   
-``./TUTORIALS/test-output/``: contains the output of the tutorials.  

#### ./OUTPUT/  
This directory is not in the repository for storage reasons. However, it is needed, as the data are assumed to be in there. The directories are organized by temperature, system size, and potential type (e.g. `./OUTPUT/T2.0/N65/xplor/`).

#### ./THERMALIZE/  
This is where the main code is stored, and also some data that the code can use. There are three subfolders

- `./THERMALIZE/progs/`:  
All the programs.


- `./THERMALIZE/script/`:  
All the scripts.

- `./THERMALIZE/data/`:  
Some very summarized data.

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
This directory is where I test new ideas of code. It is meant to be the disordered place where things are born. Apeiron.

#### ./PLOTS/  
Directory devoted to making plots. Not all the plot-making programs are here (in some cases I like to end the simulation with a graph.png), but if a program is purely devoted to plots, then it will be here.

#### ./PAPERS/  
A directory containing some random research papers.

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

Temperatures are *T* = 10.0, 5.0, 2.0, 1.0, 0.8, 0.7 0.6, 0.55, <s>0.52, 0.49, 0.466, 0.45, 0.44, 0.43</s>.

When doing a new temperature, remember to make the file params.in inside the appropriate output directory.

<s>10</s> 5 samples per temperature.

### Generating Thermalized Configurations
```
# Access script directory
cd ./THERMALIZE/script/
# Launch thermalization script for T=5.0,2.0 1.0, with 10 samples (default is 10 samples)
nsam=10 bash ThermalizeN1080.sh "5.0 2.0 1.0"
cd -
```

The script `ThermalizeN1080.sh` reads the file `./THERMALIZE/data/thermalizationtimes.txt`, which is a table of estimated thermalization times &tau;<sub>est</sub> (from previous runs and from literature [arXiv:1203.3392](https://arxiv.org/pdf/1203.3392.pdf) and [arXiv:0805.3104](https://arxiv.org/pdf/0805.3104.pdf)), translates it into number of MD steps, and multiplies it by 10.
Then, it invokes the program `ReadAndThermalize.py` which takes care of running the remaining amount of steps.

### Making sure the configurations are well-thermalized
```
# Access script directory
cd ./THERMALIZE/script/

# Edit script to change simulation parameters, such as temperature (TLIST), number of samples (nsam), system size (N), and more
emacs CheckThermalizationAll.sh

# Calculate trajectories and self-intermediate scattering functions for each sample
bash CheckThermalizationAll.sh

# Calculate average self-intermediate scattering functions
# Syntax: bash MediasFkt.sh "T1 T2 ... Tm" "N1 N2 ... Nn"
# Can also set the potential type by setting the environment variable pot_mode=xplor(default),shift, no_shift
pot_mode=xplor bash MediasFkt.sh "5.0 2.0 1.0" "1080"
cd -

# Plot self-intermediate scattering functions
cd ./PLOTS/
emacs Fkt.gp # Make sure that the parameters are the correct ones
gnuplot Fkt.gp
cd -
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

```
cd ./THERMALIZE/script/
bash CreateTrajectories.sh 15 "1.0 0.49" "1080" "0 1 2 3 4"

# Default thermostat is NVT Nose-Hoover. To change thermostat use the following syntax.
thermostat='NVE' bash CreateTrajectories.sh 15 "1.0 0.49" "1080" "0 1 2 3 4"

# Same for the potential mode (xplor is default, other options are shift, no_shift)
pot_mode='xplor' bash CreateTrajectories.sh 15 "1.0 0.49" "1080" "0 1 2 3 4"
cd -
```

The trajectories of each sample can be found in a subdirectory called `trajectories/`.

### Calculating Average Trivial Correlation Functions
From the previous trajectories, which output .npy files with positions, velocities and accelerations, we can calculate the following correlators:

- Calculate short-time correlations: force-force C<sup>FF</sup>(t), force-momentum C<sup>FP</sup>(t), momentum-momentum C<sup>PP</sup>(t). Note: I am actually measuring accelerations and velocities instead of forces and momenta, but in this case there is no difference because the mass is equal to one.

- Calculate Diagonal Correlations C<sub>d</sub>(t): this is the diagonal part of the self correlation as defined in [arXiv:1705.00036](https://arxiv.org/abs/1705.00036).

- Calculate Self-intermediate scattering function: it is not critical to calculate it again, but it can be useful to make sure everything is in order.

Of the ones just presented, only C<sup>FP</sup>(t) and C<sup>FF</sup>(t) are necessary to calculate the noise correlations.


The same program calculates all of them, but to save calculation time one can restrict the calculation to only few of them ( C<sub>d</sub>(t) is slow to calculate).

```
cd ./THERMALIZE/script
#arguments: <observables> <T-list> <N-list> <thermostat-list>

# To calculate everything at T=5.0; N=1080 with NVE thermostat
bash CalculateCorrelations.sh "--msd --Fkt --CFF --CFP --CPP --Cd" "5.0" "1080" "NVE"

# Only Diagonal correlations, limiting the input data to only two trajectories; at T=5.0,1.0; N=1080 with NVT thermostat
limit_input=2 bash CalculateCorrelations.sh "--Cd" "5.0 1.0" "1080" "NVT"
cd -
```
The output correlation functions are saved both in binary and text format, in the directory `./OUTPUT/T$T/N$N/`, with self-explanatory names.

### Calculating Noise Correlation Functions
At this point, the final step is reading the previously calculated correlation functions, and use them as kernels for calculating the noise correlation function. 

```
cd ./THERMALIZE/script
# arguments: <T-list> <N-list> <thermostat-list>
bash CalculateNoiseCorrelations.sh "5.0 1.0" "1080" "NVT"
cd -
```
Some options can be given to the script:

- `maxtime`: reduces the total integration time to maxtime

- `shiftCFP`: if not set, nothing happens. If set to anything, sets to zero the first point of CFP, since C<sup>FP</sup>(t=0)=0

- `softening`: if not set, nothing happens. If set to anything, introduces a damping term in CFP and CFF so that they are smaller at high times where the signal-to-noise ratio is low. The implementation I did is sloppy because I don't need this option in the long term.

- `tstar`: the non-selfconsistent integration is truncated at tstar (choose tstar so that it is the time after which the kernels are only noise). The implementation I did is sloppy because I don't need this option in the long term.

- `normalsc`: If set to anything, also calculates the noise correlation function *selfconsistently*.

- `lin`: If set to anything, also calculates the noise correlation function using a *linear* grid.

- `linsc`: If set to anything, also calculates the noise correlation function *selfconsistently* using a *linear* grid.

For example one can launch in the following way:

`maxtime=0.9 normalsc=1 lin=1 linsc=1 shiftCFP=1 softening=1 bash CalculateNoiseCorrelations.sh "5.0 1.0" "1080" "NVE"`

### Consistency checks on Correlations

The noise correlation is the memory kernel of the velocity correlation. Therefore, it must satisfy

C&#x307;<sup>P</sup>(t) = - &int;<sub>0</sub><sup>t</sup> K(t-s) C<sup>P</sup>(s) ds  .

Since 

C&#x307;<sup>P</sup>(t) = -C<sup>FP</sup>(t) = C<sup>PF</sup>(t),

the expression can be rewritten as

C<sup>FP</sup>(t) = &int;<sub>0</sub><sup>t</sup> K(t-s) C<sup>P</sup>(s) ds .

These last two relations can and should be verified numerically.

```
cd ./THERMALIZE/script
# bash CorrelationConsistency.sh "temperatures" "sizes" "thermostats"
bash CorrelationConsistency.sh "5.0 2.0" "1080" "NVE"
cd -
```

### Calculating Friction Coefficients

The friction coefficient is the integral of the autocorrelation function. We can calculate it both on *C*<sub>d</sub>(*t*) and on *K*(*t*).


```cd ./THERMALIZE/script
emacs CalculateFriction.sh #Put the right temperatures
bash CalculateFriction.sh
cd -
```

The frictions as a function of temperature can then be found in `./THERMALIZE/data/frictions.txt` and plotted through `./PLOTS/NoiseCorr.gp`.


### Yet not implemented
These are likely the next steps in the code development:

- **Jack-Knife** computation of the errors on the noise correlation function.



---


## Metabasin Dynamics
System size is *N* = 65.

Temperatures are *T* = **10.0**, **5.0**, **2.0**, **1.0**, **0.8**, **0.6**, (0.52), (0.49), (0.466), (0.45), (0.44), (0.43).

10 samples per temperature. When doing a new temperature, remember to make the file params.in inside the appropriate output directory.

The potential now is a `shift` Lennard-Jones. For this kind of runs, I chose shift because at some point I might have to calculate a bunch of Hessians.

### Generating Thermalized Configurations

First of all, we need to thermalize the system.

```
cd ./THERMALIZE/script
nsam=10 potential='shift' dt=0.0025 backupFreq=10000 bash ThermalizeForMetabasins.sh "5.0 2.0" "65"
cd -
```

### Making sure the configurations are well-thermalized

Then, we make sure that the configurations are well-thermalized.

```
# Access script directory
cd ./THERMALIZE/script/

# Edit script to change simulation parameters, such as temperature (TLIST), number of samples (nsam), system size (N), and more
emacs CheckThermalizationForMetabasins.sh

# Calculate trajectories and self-intermediate scattering functions for each sample
bash CheckThermalizationForMetabasins.sh

# Calculate average self-intermediate scattering functions
# Syntax: bash MediasFkt.sh "T1 T2 ... Tm" "N1 N2 ... Nn"
# Can also set the potential type by setting the environment variable pot_mode=xplor,shift(default), no_shift
pot_mode=shift bash MediasFkt.sh "5.0 2.0 1.0" "1080"
cd -

# Plot self-intermediate scattering functions
cd ./PLOTS/
emacs Fkt.gp # Make sure that the parameters are the correct ones
gnuplot Fkt.gp
cd -
```



### Running new NVE trajectories with inherent-structure sampling, and calculating barriers with the Ridge Method

### Calculating barriers with the Nudged Elastic Band method

### Abandoning a metabasin

Starting from a single deep metabasin, I run several trajectories.


---

