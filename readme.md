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
Currently, there is also some code taken from [Ian Dun](https://github.com/iansdunn) for fitting functions as a sum of exponentials.


#### ./PLAYGROUND/  
This directory is where I test new ideas of code. It is meant to be the disordered place where things are born. Apeiron.

#### ./PLOTS/  
Directory devoted to making plots. Not all the plot-making programs are here (in some cases I like to end the simulation with a graph.png), but if a program is purely devoted to plots, then it will be here.

#### ./PAPERS/  
A directory containing some random research papers.

#### ./UTILITIES/
Some scripts for doing useful straightforward things, such as extracting instant information from gsd files. Program names are self-explanatory.

#### ./EXPERIMENTS/  
Some runs I did to check something specific.

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

Temperatures are *T* = 10.0, 5.0, 2.0, 1.0, 0.8, 0.7 0.6, 0.55, 0.52, 0.5, 0.49, 0.48 0.47, 0.46, 0.45.

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

When launching from the cluster (**talapas**), it can be useful to set the queue and the total time for the simulation (backups are done regularly, so it is not a problem if the runtime is shorter than the thermalization time):
```
queue=longgpu simTime="1-23:00:00" nsam=10 bash ThermalizeN1080.sh "0.48"
```
Default is queue=gpu, simTime="0-02:00:00". The format is days-hours:minutes:seconds.

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
Even though the thermostat can be chosen, the trajectories used for the paper are all run with NVT.


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

To plot some of these correlations, do 
```
cd ./PLOTS
gnuplot tricial-correlations.gp
okular FIGURES/Fkt.eps
okular FIGURES/CFF.eps
okular FIGURES/CFP.eps
okular FIGURES/CPP.eps
cd -
```

Alternatively, one can use `CalculateCorrelationsJK.py`. This other program is very similar. It has the advantage of computing the jackknife blocks (to allow statistical analysis) and
of constructing blocks while reading, which allows for minimal memory usage. The down side is that it is less elastic, and all the input data (the trajectories) have to be very consistent, otherwise it will throw an error. The syntax is the same (there is the extra argument `lblo` for the length of the blocks), though the calculation of C<sub>d</sub>(t) is not implemented because not necessary.
The script to launch it is `CalculateCorrelationsJK.sh`, which has the same exact syntax as `CalculateCorrelations.sh`, and by default makes 10 JK blocks, which is not a lot, but it saves a lot of memory and computation time.

```
cd ./THERMALIZE/script
bash CalculateCorrelationsJK.sh "--CFF --CFP --CPP" "5.0" "1080" "NVT"
cd -
```

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

- `fits`: If set to anything, instead of making interpolations does fits.

For example one can launch in the following way:

`maxtime=0.9 fits=1 normalsc=1 lin=1 linsc=1 shiftCFP=1 softening=1 bash CalculateNoiseCorrelations.sh "5.0 1.0" "1080" "NVE"`

To plot the noise correlations, you can use the following gnuplot script:
```
cd ./PLOTS
gnuplot NoiseCorr.gp
cd -
```


If the JackKnife blocks of the trivial correlations were `CalculateCorrelationsJK.sh` produced, then the noise correlations can also be calculated with JK. Use the following script, with analogous syntax of its non-JK counterpart:
```
bash CalculateNoiseCorrelationsJK.sh "5.0 2.0" "1080" "NVT"
```
(new versions of the code do the jackknife blocks directly in `CalculateCorrelations.sh` --  the user must specify where they want to do with or without JK, through the flag `noJK`)


The jackknife blocks can be plotted individually with the following script
```
cd ./PLOTS
python PlotJKblocks.py -oK -T2.0 -N1080 -tNVT -M3 --nomean
cd -
```

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
bash CorrelationConsistency.sh "5.0 2.0" "1080" "NVT"
cd -
```
The program spits a figure (in the directory of the data) of the velocity autocorrelation calculated directly from the simulation, and through integration of the memory kernel. 



Consistenct checks can also be done with JackKnife:
```
bash CorrelationConsistencyJK.sh "5.0 2.0" "1080" "NVT"
#If you want to visualize the plots while checking:
showplots=1 bash CorrelationConsistencyJK.sh "5.0 2.0" "1080" "NVT"
```


### Calculate Diffusion constants

The diffusion constants are calculated in four ways 

- from the long-time behavior of the mean square displacement, *msd/(6t)*=*D*+*k*/*t*, where *k*/*t* are subleading corrections.

- by integrating the velocity correlations (this should give the same result as the previous point)

- by integrating *C*<sub>d</sub>(*t*): D = T/(integral)

- by integrating *K*(*t*): D = T/(integral)

```
cd ./PLOTS/
gnuplot msd.gp
cd -
```

The first generates figures (`Dmsd.eps`) in `./PLOTS/FIGURES` and data on the diffusion constants in `./THERMALIZE/data/D.txt`.


### Calculating Friction Coefficients

The friction coefficient is the integral of the autocorrelation function, divided by the temperature. We can calculate it both on *C*<sub>d</sub>(*t*) and on *K*(*t*).

The following script calculates the friction coefficient on the noise and diagonal correlation functions.
Further, the short-time behavior is fitted through a form *f*<sub>short</sub>(*x*)=*a*<sub>1</sub>/ cosh(*a*<sub>2</sub> *x*), as described e.g. in Eq.(3) of [arXiv:cond-mat/0109285](https://arxiv.org/abs/cond-mat/0109285).
The friction coefficient is calculated also on the short-time fitted function, and on the long-time one (i.e. total minus short).
All the quantities are output in the same directory of the correlation function.


```
cd ./THERMALIZE/script
bash CalculateFriction.sh "5.0 2.0 1.0 0.8 0.7 0.6 0.55 0.52 0.49 0.47 0.46" "1080" "NVT"
cd -
```

To plot the 2x3=6 friction coefficients and their fits, do
```
cd ./PLOTS
gnuplot friction.gp
```


(In the old version the frictions as a function of temperature could be found in `./THERMALIZE/data/frictions.txt` and plotted through `./PLOTS/NoiseCorr.gp`.)



### Long and short-time behavior of the correlation functions

To plot the curvatures of the short-time behavior of both correlation functions, just do:
```
cd ./PLOTS
gnuplot short-times.gp
```

For the long-time behavior use the following gnuplot script.
```
gnuplot long-times.gp
cd -
```

### Standard mode-coupling theory memory function

Calculate the memory function by passing to Fourier space, and factorizing the correlator, and assuming that the normal dynamics is similar to the orthogonal dynamics.

```
cd ./THERMALIZE/script
#Calculate the contribution for each k (the list of n vectors is contained in ../data/n.txt)
ntCd=25 ntw=20 bash CalculateVertexCorrelations.sh "5.0" "1080" "NVT"
#Integrate in dk
showplots=1 bash CalculateKvertex.sh "5.0" "1080" "NVT"
cd -
```

To plot the function,
```
cd ./PLOTS
gnuplot NoiseCorrCompare3.gp
cd -
```



### Yet to do

- Some stuff needs to be pushed, and the readme updated.

--- 

--- 

## Metabasin Dynamics
System size is *N* = 65.

Temperatures are *T* = **10.0**, **5.0**, **2.0**, **1.0**, **0.8**, **0.6**, **0.49**, **0.46**.

10 samples per temperature. When doing a new temperature, remember to make the file `params.in` inside the appropriate output directory.

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
cd ./THERMALIZE/script
#Launch for T=5.0, N=65, samples 0,4,6,7,8, shift potential, NVT thermostat
pot_mode='shift' thermostat='NVT' bash CheckThermalizationForMetabasins.sh "5.0" "65" "0 4 6 7 8"
pot_mode=shift bash MediasFkt.sh "2.0" "65" "NVT"
cd -
```

The first script produces, in the output directory of each sample, two files with the self-intermediate scattering function Fk(t). The sample average is obtained through the second script, that outputs the average curve. 
The curves can the be plotted with the gnuplot script `../PLOTS/Fkt.gp`.


### A rough estimate of tau alpha

We need a rough estimate of tau_alpha. This is easily done by doing

```
cd ./THERMALIZE/script
thermostat='NVT' bash CalculateTauAlpha.sh
cd -
```

The output goes to `./THERMALIZE/data/tau_alpha.txt`.

### Generating a long trajectory
We now want to simulate a very long trajectory (1000 tau_alpha) on sample 0, in order to make the metabasin analysis on it. The way we do it is in chunks, so that from each chunk we can extract the IS trajectory independently.

The following script sequentially creates a trajectory chunk, and then minimizes it, before passing to the following chunk. Of course, it is easily modified in order to run all chunks first, and then minimize them all. It is mainly a matter of resource availability (memory or C(G)PU).

```
cd ./THERMALIZE/script
bash ChunkSectAll.sh "0.6" "65"
cd -
```

This program saves in the output directory (`......./chunkIS/`) the list of inherent structures, the list of ridges (this can be disactivated for more speed), the thermal trajectory, and a **dictionary containing the inherent structures and the time at which they occur (to do)**.



### Metabasin trajectory
At this point we have the inherent structure trajectory. Now we want to identify the metabasins.

```
python IdentifyMetabasins.py ~/STRUCTURAL-GLASS/OUTPUT/T0.7/N65/shift/S0/chunksIS/elistIS.txt
```

---

### Running new NVE trajectories with inherent-structure sampling, and calculating barriers with the Ridge Method

### Calculating barriers with the Nudged Elastic Band method

### Abandoning a metabasin

Starting from a single deep metabasin, I run several trajectories.


---

