# ABOUT

This repository hosts _MCQsim_, a multi-cycle earthquake rupture simulator for generation the of earthquake sequences along predefined faults or fault systems of arbitrary geometric complexity. _MCQsim_ uses graphical user interfaces (GUIs) for most of data input and output. These GUIs were developed in _MATLAB_. Stand-alone versions of these GUIs, along with the _MATLAB_ GUIs can be found here, e.g., under _"01_Mstuff"_. _MCQsim_ itself is written in **_C_** and parallelized via **_MPI_**. It works on large-scale computational infrastructure as well as personal computers. The code can be found under _"00_Cstuff"_. 

The most recent version of _MCQsim_, provided here, employs _QuadTree_ and _H-matrix_ to decrease the size of the stiffness kernel. This update dramatically increased code efficiency and scalability, especially for larger-scale simulations with 100k+ fault elements. Future releases are planned to include pore-pressure changes as well as surface topography, further guided by the specific needs by our fellow users (i.e., you).

**Please check the following publication for details on _MCQsim_ model formulation and cite when using the code**

> Zielke, O., and P.M. Mai (2023). MCQsim: A Multicycle Earthquake Simulator. Bull. Seismol. Soc. Am. 113(3), 889-908. doi: 10.1785/0120220248

**Contact:** olaf.zielke@kaust.edu.sa

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

|   File name   (input)                 | File content  |
| ------------------------------------- | ------------- |
| FaultModel.txt                        | Parameter file. Contains details about how to run the simulation e.g., catalog length. ASCII format.  |
| FaultModel_FLT.txt                    | Summary file. Contains information about fault geometry. ASCII format.  |
| FaultModel_FLTtrig.dat                | Data file. Contains the 3-D fault geometry, including the list of all triangular fault elements. Binary format.  |
| FaultModel_Rough.dat                  | Data file. Contains the 3-D fault geometry. This file is equivalent to _FLTtrig.dat except that fault roughness might be added. Binary format.  | 
| FaultModel_Strgth.dat                 | Data file. Contains friction and strength parameter for each fault element. Binary format.  |
| FaultModel_BND.txt      _(optional)_  | Summary file. Contains information about of boundary surface geometry. ASCII format.  | 
| FaultModel_BNDtrig.dat  _(optional)_  | Data file. Contains the 3-D geometry of boundary surface triangular elements. Binary format. | 

## RUNNING MCQSIM
Assuming all files are present, a simulation can be started via:
>
>  _mpirun -n 6 ./MCQsim FaultModel.txt_

This command will start MCQsim as an MPI run, using n = 6 CPUs. The parameter file _"FaultModel.txt"_ contains relevant information to start/control the simulation. _MCQsim_ then runs and produces a number of output files. Most importantly the _" *Catalog.dat"_ file and a number of _" *.srfb"_ files.

## MCQSIM OUTPUT FILES
|   File name   (output)               | File content  |
| ------------------------------------ | ------------- |
| FaultModel_BrchInfo.dat              | Data file. Contains information about the QuadTree structure. Useful for _OpenQuake_ pipeline. Binary format.  |
| FaultModel_Kh.dat _(optional)_ | Data file. Contains stiffness matrix for elastic interaction. Allows continuing an existing simulation/catalog without re-calculating. Binary format. !**_large file size_**!  |
| FaultModel_PreRunData.dat            | Data file. Friction properties, fault strength, long-term slip-rate of fault elements at start of simulation. Binary format.  |
| FaultModel_PostRunState.dat          | Data file. Stress state after simulation is over. Allows continuing/extending an existing catalog. Binary format.  |
| FaultModel_RAWCatalog.dat            | Data file. Earthquake catalog. Binary format.  |
| FaultModel_M7.38472_t3845.8374.srfb  | Data file. Quasi-dynamic earthquake rupture for selected large event. Name provides event magnitude and time. Tthis is a binary version of a *.srf file (standard ruputre format; Graves, 2002). Binary format.  | 

### PARAMETER FILE  
  The parameter file has a specific structure that must not be changed (e.g., by adding/shifting lines etc.). Values can be changed however. Here an example of what "ParameterFile.txt" contains. Line numbers to the left are added for reference only and not part of the file.

![this is an example parameter file](https://github.com/OlafZielke-EQ/MCQsim/blob/main/pagematerial/ParameterFileScreenShot.png)

here a brief explanation for selected entries:

_08ContinueCatalog_      - a switch ( 0 / 1) telling if an existing catalog is meant to be extended (1 == continue catalog)

_09LoadPrevious_Kh_mat_  - a switch ( 0 / 1 / 2) telling if a new Kh_matrix is computed or not and whether it is stored to file or not; 0 == new Kh_mat without saving; 1 == new Kh_mat with saving; 2 = use existing Kh_mat. A typical combination of entries in line 8 and 9 woudl be (a) line 8 is 0 and line 9 is 1; or (b) line 8 is 1 and line 9 is 2.

_10Kh_mat_file_2_load_   - name of the _Kh_mat file_ that is either written or re-used (depending on choice in line 9)

_12StoreSTF4LargeEQs_    - a switch (0 / 1) telling if the detailed quasi-dynamic rupture model of selected earthquqakes should be written to file or not (1 == write to file)

_13MinMagnitude4ST_      - threshold magnitude to store the quasi-dynamic rupture model; this is included to save disk space and because smaller events might be less interesting in that regard. IMPORTANT: "smaller" in this context mainly relates to number of elements that failed in the event. The connection between number of failed elements and this minimum magnitude is therefore the mesh (i.e. fault element) size.

_14/15_                 - currently not implemented
