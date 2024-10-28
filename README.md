# ABOUT

This repository hosts _MCQsim_, a multi-cycle earthquake rupture simulator for generation the of earthquake sequences along predefined faults or fault systems of arbitrary geometric complexity. _MCQsim_ uses graphical user interfaces (GUIs) for most of data input and output. These GUIs were developed in _MATLAB_. Stand-alone versions of these GUIs, along with the _MATLAB_ GUIs can be found here, e.g., under _"01_Mstuff"_. _MCQsim_ itself is written in **_C_** and parallelized via **_MPI_**. It works on large-scale computational infrastructure as well as personal computers. The code can be found under _"00_Cstuff"_. 

The most recent version of _MCQsim_, provided here, employs _QuadTree_ and _H-matrix_ to a) decrease the size of the stiffness kernel, and b) decrease computation time. This update dramatically increased code efficiency and scalability, especially for larger-scale simulations with 100k+ fault elements. Future releases are planned to include pore-pressure changes as well as surface topography, further guided by the specific needs by our fellow users (i.e., you).

**Please check the following publication for details on _MCQsim_ model formulation and cite when using the code**

> Zielke, O., and P.M. Mai (2023). MCQsim: A Multicycle Earthquake Simulator. Bull. Seismol. Soc. Am. 113(3), 889-908, doi:10.1785/0120220248. [*.pdf (open access)](https://pubs.geoscienceworld.org/ssa/bssa/article/113/3/889/622487/MCQsim-A-Multicycle-Earthquake-Simulator)

**Contact:** olaf.zielke@kaust.edu.sa

## SAMPLE OUTPUTS
### Single event and earthquake catalogs

Following are a few example outputs, generated with MCQsim for a "toy example" listric normal fault that features geometric roughness and heterogeneous strength distribution. We present animated GIFs of the single-event quasi-dynamic evolution of a **M**7.54 event (slip and slip-rate). We further present the earthquake catalog and catalog derivates for the same fault.

|   M7.54 normal faulting EQ (slip)      | M7.54 normal faulting EQ (rate)  |
| -------------------------------------  | ------------- |
| ![animated GIF of toy example rupture](https://github.com/OlafZielke-EQ/MCQsim/blob/main/pagematerial/M754_slip.gif) | ![screenshot of website for easy download of GitHub folders](https://github.com/OlafZielke-EQ/MCQsim/blob/main/pagematerial/M754_rate.gif) |


![screenshot of toy example catalog and derivatives](https://github.com/OlafZielke-EQ/MCQsim/blob/main/pagematerial/ExampleCatalogOutput.png)

 ### Earthquake catalogs for PSHA/OpenQuake

Earthquake catalogs, created with _MCQsim_, can serve as an input for PSHA, providing an earthquake rupture forecast for a modeled fault (system). As such, they provide internally consistent, physics-based branches of the PSHA logic tree. We recently streamlined the post-processing workflow to feed _MCQsim_ earthquake catalogs into the [_OpenQuake_](https://www.globalquakemodel.org) PSHA engine. Below is an example seismic hazard map for toy model listric fault, using [_QGIS_](https://www.qgis.org) and [_OpenQuake's_](https://www.globalquakemodel.org) [_IRM toolkit_](https://docs.openquake.org/oq-irmt-qgis/v3.1.0/).

![screenshot of toy example PSHA analysis](https://github.com/OlafZielke-EQ/MCQsim/blob/main/pagematerial/MCQsim2PSHA.png)

## SCALING TESTS

Following are the results of strong scaling tests that were done on [Shaheen3@KAUST](https://www.hpc.kaust.edu.sa). We show them here a) to showcase _MCQsim_ performance on our HPC facilities, and b) to provide users with first-order guidence on simulation preparation (e.g., CPU number as a function of fault element number i.e., fault system size). Below, we show scaling in terms of speed up factor (relative to single-CPU simulation) and effciency for a two-fault simulation using 10k, 50k, 100k, and 200k fault elements. The simulations indicate that the speed up factor closely follows the ideal case, before scaling begins to come less efficient (at different CPU numbers for the different fault element numbers). The "odd" shapes in effciency curves for the 10k and 50k element simulations reflect the internal CPU organization within a single node of Shheen3 (each node consisting of 192 CPUs, single-node indicated by gray area). Given these scaling tests we suggest to use: 6-12 CPUs for 10k elements, 24-48 CPUs for 50k elements, 384-768 CPUs for 100k elements, and 768 to 1536 CPUs for 200k elements.

***IMPORTANT:*** Keep in mind that the results of this scaling test are specific to the model setup and _Shaheen3_. It is best-practice to perform these tests yourself on the HPC facilities you are going to use.

![Strong scaling plot for 4 different grid resolutions](https://github.com/OlafZielke-EQ/MCQsim/blob/main/pagematerial/ScalingPlot.png)

|   Element number      |  Resolved Elements | Stiffness matrix (H-mat)  |  Stiffness matrix (classical, N^2)  |
| --------------------- | ------------------ | ------------------------- | ----------------------------------- |
|  10,000               |    2,770           | 1.1 Gb                    | 3.6 Gb                              |
|  50,000               |    3,656           | 14.1 Gb                   | 90 Gb                               |
|  100,000              |    3,248           | 45.5 Gb                   | 360 Gb                              |
|  200,000              |    3,775           | 182.1 Gb                  | 1,440 Gb                            |

Our H-matrix implementation drastically reduces the number of "resolved elements" and therefore the number of individual fault interactions that need to be computed.

## DOWNLOAD
  Here is a simple way to download specific folders from GitHub, as opposed to downloading the entire project:

> open a browser tab and go to
> 
>  https://download-directory.github.io
> 
> then copy and paste the GitHub folder you want to download, e.g.,
> 
> https://github.com/OlafZielke-EQ/MCQsim/tree/main/TREAD/
> 
>![screenshot of website for easy download of GitHub folders](https://github.com/OlafZielke-EQ/MCQsim/blob/main/pagematerial/GitHubDownload.png)
>
> into the corresponding field and the download will start automatically.

# MCQSIM SIMULATION
The _MCQsim_ worklflow can be broken down into three components, namely pre-simulation, simulation, and post-simulation. The first and third component rely heavily on MATLAB-based GUIs that we developed for this purpose. The actual simulation is "terminal" based.

## MCQSIM INPUT FILES
_MCQsim_ simulations require a number of files to be located in the same folder as the executable. They are created with the _MATLAB_ input GUIs.

|   File name   (input)                  | File content  |
| -------------------------------------  | ------------- |
| FaultModel.txt                         | Parameter file. Contains details about how to run the simulation e.g., catalog length. ASCII format.  |
| FaultModel_FLT.txt                     | Summary file. Contains information about fault geometry. ASCII format.  |
| FaultModel_FLTtrig.dat                 | Data file. Contains the 3-D fault geometry, including the list of all triangular fault elements. Binary format.  |
| FaultModel_Summary_RoughStrength.txt   | Summary file. Contains information about used fault roughness, strength, and friction configuration. ASCII format.  |
| FaultModel_1_Rough.dat                 | Data file. Contains the 3-D fault geometry. This file is equivalent to _FLTtrig.dat except that fault roughness might be added. Binary format.  | 
| FaultModel_1_Strgth.dat                | Data file. Contains friction and strength parameter for each fault element. Binary format.  |
| FaultModel_BND.txt      _(optional)_   | Summary file. Contains information about of boundary surface geometry. ASCII format.  | 
| FaultModel_BNDtrig.dat  _(optional)_   | Data file. Contains the 3-D geometry of boundary surface triangular elements. Binary format. | 

### PARAMETER FILE  
The parameter file has a specific structure that **_must not_** be changed. **_Do not include or remove lines and do not enter spaces in the descriptions_**. Following is an example parameter file _"FaultModel.txt"_.

```
Run_Parameter_Information
--------------------------
InputName:                     FaultModel
Realization_Number:            1
CatalogTypeNr(1/2/3):          3
PlotCatalog2Screen(0/1):       1
RandSeedValue:                 1001
MinElementNum4Catalog:         1
ContinueCatalog(0/1):          0
LoadPrevious_Kh_mat(0/1/2):    1
Kh_mat_file_2_load:            FaultModel_Khmat.dat
--------------------------
StoreSTF4LargeEQs(0/1):        1
MinMagnitude4STF:              7.0f
UseRuptPropag(0/1):            0
MinMagnitude4RupProp:          6.5f
--------------------------
IntSeisLoadStep(days):         1.0f
pow2forInterSeisSteps:         8
--------------------------
Visco_AfterSlip(years):        3.0f
Visco_DeepRelax(years):        20.0f
--------------------------
HealingDuringCoseisFraction:   0.05f
OvershootFraction:             1.2f
PreStressFraction:             0.90f
MinCoSeisSlipRate(m/s):        5.0E-3
--------------------------
EQrecordLength(years):         5000.0f
--------------------------
```

A brief explanation of the parameter file entries
|   Parameter information              | Description  |
| ------------------------------------ | ------------- |
| InputName:                    |   This is the "root" part of file name e.g., _"FaultModel"_ of file _"FaultModel_FLT.txt"_. All input files of the simulation need to share this input (see input files example above). |
| Realization_Number:           |   _MCQsim_ was built to facilitate parameter space investigations, where multiple realizations can be generated with input GUI "_B_MCQsim_RoughStrgth_" (while only changing selected parameter). The realization number refers to an index in the file name of the _"*Roughn.dat"_ and _"*Strgth.dat"_ file. For example. For example, for realization 1 would, these files would be named _"FaultModel_1_Roughn.dat"_ and _"FaultModel_1_Strgth.dat"_; for realization 42, this would be _"FaultModel_45_Roughn.dat"_ and _"FaultModel_45_Strgth.dat"_|
| CatalogTypeNr(1/2/3):         |    Catalog type refers to the amount of information that is provided in the catalogs. **_Type 1_** contains only general event characteristics such as time, magnitude, hypo-center, moment-centroid, rupture area, mean slip, and mean stress drop. **_Type 2_** further contains the event's moment rate function. **_Type 3_** further contains information for each fault element that was active (slipping) during the event. This includes the element's rupture start time, slip in strike- and dip-direction, coseismic (static) stress change, and stability type. _MCQsim_ first generates a "RAW" catalog that corrresponds to _Type 3_. This "RAW"catalog can be processed as part of the post-simulation workflow to generate other catalog types that are smaller in size (not every investigation may have the need to consider element-based information).  |
| PlotCatalog2Screen(0/1):      |    This switch (0 = no, 1 = yes) indicates if catalog information should be printed to screen during the simulation.  It is typically left at 1. |
| RandSeedValue:                |    If desired, some parameters can be randomized during simulation. For example, friction parameters of activated fault elements may be permitted to vary (within specified range) at the end of the coseismic phase. Reusing the same seed value allows to create the exact earthquake sequence. **_IMPORTANT_** each process/CPU gets its own seed value (the process' MPIrank is added to _"RandSeedValue"_). Hence, in order to create the exact same catalog, it is important to not only use the same seed value, but also the same number of processes. |
| MinElementNum4Catalog:        |    Minimum number of fault elements that need to be activated in an event to add the event to the catalog. If set to 1, all events are included in the catalog; if set to 10, only events with 10+ active elements are included in the catalog.  **_IMPORTANT_** The element number refers **_only_** to the writing part (the smaller events still exist in the simulation, they are just not written to file).  | 
| ContinueCatalog(0/1):         |   This switch (0 = no, 1 = yes) indicates if an existing catalog should be extended.  |
| LoadPrevious_Kh_mat(0/1/2):   |   This switch indicates if a new stiffness matrix _"*Khmat.dat"_ should be created or not, and whether it should be saved to file (e.g., to be used when continuing an existing catalog). **_switch = 0_** means to create a new stiffness matrix but **_NOT_** storing it to file. **_switch = 2_** means to create a new stiffness matrix and storing it to file. **_switch = 2_** means to load an existing stiffness matrix. A typical setup for this and the previous entry would be _ContinueCatalog 0_ and LoadPrevious_Khmat 1_ or _ContinueCatalog 1_ and LoadPrevious_Khmat 2_ |
| Kh_mat_file_2_load:           |    The file name used to either read (if _LoadPrevious_Kh_mat  2_) or write (if _LoadPrevious_Kh_mat  1_) the stiffness matrix _"*Khmat.dat"_  |
| StoreSTF4LargeEQs(0/1):       |   This switch (0 = no, 1 = yes) indicates whether the quasi-dynamic rupture evolution of individual events should be saved to file.  |
| MinMagnitude4STF:             |   The minimum magnitude an event must have in order to write its quasi-dynamic rupture evolution to file.  |
| UseRuptPropag(0/1):           |   Not used in current code version.   |
| MinMagnitude4RupProp:         |   Not used in current code version.   | 
| IntSeisLoadStep(days):        |   The minimum time step (in days) used for post-seismic i.e., interseismic loading and stress redistribution.   |
| pow2forInterSeisSteps:        |   _MCQsim_ uses the minimum time step (previous row) and a binary search to determine the loading step to the next event. The value of _pow2forInterSeisSteps_ refers to just that. For example with a value of 8, _MCQsim_ first jumps 2^8 days ahead in time to identify if the combination of inter-seismic loading and post-seismic stress redistribuiton caused at least one element to fail. If not, loading proceeds. If yes, then the loading step is halfed, etc.    |
| Visco_AfterSlip(years):       |   Postseismic response ALONG THE FAULT, employing a simple Maxwell spring-dashpot model. This is the time it takes to decrease a postseismic signal to 1/e (approx. 35%) of its original value. Providing this measure in years (rather than using viscosity) seems more intuitive. For example the postseismic signal/effect (here after slip) decreased by ~65% within the first 3 years. Higher value corresponds to higher viscosity. | 
| Visco_DeepRelax(years):       |   Postseismic response ALONG THE BOUNDARY (if used), employing a simple Maxwell spring-dashpot model. This is the time it takes to decrease a postseismic signal to 1/e (approx. 35%) of its original value. Providing this measure in years (rather than using viscosity) seems more intuitive. For example the postseismic signal/effect (here visco-elastic relaxation) decreased by ~65% within the first 20 years. Higher value corresponds to higher viscosity.  | 
| HealingDuringCoseisFraction:  |   Fraction by which the friction coefficient of activated (unstable or conditionally stable) fault elements increases coseismically (during each coseismic iteration step) once these elements stopped to slip in the event (while the rupture itself is still ongoing).   | 
| OvershootFraction:            |   Describes the ratio of dyanamic and arrest friction coefficient.   |
| PreStressFraction:            |   Faults are pre-stressed at the beginning of a simulation to this fraction of their frictional strength (unless an existing catalog is continued). This is done to speed up the intial "burn in" phase i.e., decreasing the computational time for the modeled fault system to reach its dynamic equilibrium.  |
| MinCoSeisSlipRate(m/s):       |   Minimum slip rate to be considered coseismic slip.   |
| EQrecordLength(years):        |   Record length in years. If a catalog is continued, this value refers to the added length and not to the total length.  |


## RUNNING MCQSIM
Assuming all files are present, a simulation can be started via:
>
>  _mpirun -n 6 ./MCQsim FaultModel.txt_

This command will start MCQsim as an MPI run, using n = 6 CPUs. The parameter file _"FaultModel.txt"_ contains relevant information to start/control the simulation. _MCQsim_ then runs and produces a number of output files. Most importantly the _" *Catalog.dat"_ file and a number of _" *.srfb"_ files.

## MCQSIM OUTPUT FILES
_MCQsim_ creates a number of output files, as listed below. They can be visualized using our MATLAB-based GUIs. Catalog files can be opened with GUIs starting with a _"Z"_ such as _"Z_Catalog_MCQsim_TREAD"_ or _"Z_MCQsim_LoadCatalog"_. They also allow you to export the catalog to either _*.mat_ or _*.txt_ file. The PreRun file as well as the quasi-dynamic rupture models of individual events can be opened with GUIs starting with _"Y"_ such as _"Y"_PreRunAndSRF_MCQsim_Tread_.

|   File name   (output)               | File content  |
| ------------------------------------ | ------------- |
| FaultModel_BrchInfo.dat              | Data file. Contains information about the QuadTree structure. Useful for _OpenQuake_ pipeline. Binary format.  |
| FaultModel_Khmat.dat _(optional)_ | Data file. Contains stiffness matrix for elastic interaction. Allows continuing an existing simulation/catalog without re-calculating Kh. Binary format. !**_large file size_**!  |
| FaultModel_PreRunData.dat            | Data file. Friction properties, fault strength, long-term slip-rate of fault elements at start of simulation. Binary format.  |
| FaultModel_PostRunState.dat          | Data file. Stress state after simulation is over. Allows continuing/extending an existing catalog. Binary format.  |
| FaultModel_RAWCatalog.dat            | Data file. Earthquake catalog. Binary format.  |
| FaultModel_M7.38472_t3845.8374.srfb  | Data file. Quasi-dynamic earthquake rupture for selected large event. Name provides event magnitude and time. Tthis is a binary version of a *.srf file (standard ruputre format; Graves, 2002). Binary format.  | 



