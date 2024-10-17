This repository hosts MCQsim, a multi-cycle earthquake rupture simulator for generation of earthquake sequences along predefined faults or fault systems of arbitrary geometric complexity.
The updated code version includes QuadTree and H-matrix approximation of the stiffness kernel to improve code efficiency and scalability.



Please check out the following publication for further details:

Zielke, O., and P.M. Mai (2023). MCQsim: A Multicycle Earthquake Simulator. Bull. Seismol. Soc. Am. 113(3), 889-908. doi: 10.1785/0120220248
Contact: olaf.zielke@kaust.edu.sa



DOWNLOADING MCQsim
  To download individual folders:

  Open a browser tab and go to
  https://download-directory.github.io

  Copy and paste the GitHub site you want to download, such as:
  https://github.com/OlafZielke-EQ/MCQsim/tree/main/forTREAD/IO_standalone_MAC
  for the OSX/MAC version of the standalone I/O GUIs 

  or:
  https://github.com/OlafZielke-EQ/MCQsim/tree/main/00_Cstuff
  for the MCQsim source code 

  This will download the folders.



INSTALLATION
MCQSIM INPUT/OUTPUT
  I/O includes model and parameter setup (input) and various ways to visualize simulations results (output, e.g., individual earthquakes and earthquake catalogs).
  Input/output for MCQsim is done through MATLAB-based GUIs i.e., standalone versions thereof (i.e., no MATLAB licence required). The latter are currently for OSX/MAC and Ubuntu.
 
  for MATLAB: start the *.p files by typing for example:
  A_Geometry_MCQsim_TREAD
  in the MATLAB command window and then press enter.

  for standalone GUIs:
  Install by double-clicking the *.app file (e.g., Ainput_MCQsim.app) and then follow the instructions. 
  IMPORTANT: The installation process includes downloading MATLAB resources that are needed to run the GUIs (automatically done during installation, accessing the mathworks.com page).
  IMPORTANT: It appears that sudo rights are necessary to install the standalones on Ubuntu.



MCQSIM CODE COMPILATION
  MCQsim is written in C, parallelized with MPI. It further uses BLAS (provided via GSL or else). These libraries must be accessible during compilation/run-time.
  MCQsim compilation/installation is done in the command window i.e., terminal.

Compilation on OSX/MAC
  mpicc   MCQsim_Main.c   MCQsim_StrainHS.c   MCQsim_StrainFS.c   -lm  -lgslcblas  -lmpich  -Wall  -O3  -o  MCQsim

  this compiliation assumes use of BLAS provided by GSL (hence, GSL needs to be accessible as well).

Compilation on OSX/MAC
  mpicc   MCQsim_Main.c   MCQsim_StrainHS.c   MCQsim_StrainFS.c   -lm  -lgslcblas  -lmpich  -Wall  -O3  -o  MCQsim

  It might be necessary to check if the required modules are available/accessible. Check via:

  module avail

  and then load the required modules (such as GSL and/or MPI) for example

  module load mpich
  module load gsl



MCQSIM SIMULATION
  Assuming you have all required input and parameter files defined, you can start the simulation via:

  mpirun -n 6 ./MCQsim ParameterFile.txt

  This command will start MCQsim as an MPI run, using n = 6 CPUs. The parameter file contains relevant information to start/control the simulation.

PARAMETER FILE  
  The parameter file has a specific structure that must not be changed (e.g., by adding/shifting lines etc.). Values can be changed however. Here an example of what "ParameterFile.txt" contains. Line numbers to the left are added for reference only and not part of the file.

01--------------------------
02InputName:                     MyFault
03Realization_Number:            1
04CatalogTypeNr:                 3
05PlotCatalog2Screen:            1
06RandSeedValue:                 1001
07MinElementNum4Catalog:         1
08ContinueCatalog:               0
09LoadPrevious_Kh_mat:           1
10Kh_mat_file_2_load:            MyFault_Khmat.dat
11--------------------------
12StoreSTF4LargeEQs(0/1):        1
13MinMagnitude4STF:              7.3f
14UseRuptPropag(0/1):            0
15MinMagnitude4RupProp:          6.0f
16--------------------------
17IntSeisLoadStep(days):         1.0f
18pow2forInterSeisSteps:         9
19--------------------------
20Visco_AfterSlip(years):        3.0f
21Visco_DeepRelax(years):        20.0f
22--------------------------
23HealingDuringCoseis:           0.025f
24OvershootFraction:             1.2f
25PreStressFraction:             0.85f
26MinCoSeisSlipRate(m/s):        5.0E-3
27--------------------------
28EQrecordLength(years):         30000.0f
29--------------------------

here a brief explanation for selected entries:

08ContinueCatalog      - a switch ( 0 / 1) telling if an existing catalog is meant to be extended (1 == continue catalog)
09LoadPrevious_Kh_mat  - a switch ( 0 / 1 / 2) telling if a new Kh_matrix is computed or not and whether it is stored to file or not; 0 == new Kh_mat without saving; 1 == new Kh_mat with saving; 2 = use existing Kh_mat
A typical combination of entries in line 8 and 9 woudl be (a) line 8 is 0 and line 9 is 1; or (b) line 8 is 1 and line 9 is 2.
10Kh_mat_file_2_load   - name of the Kh_mat file that is either written or re-used (depending on choice in line 9)
12StoreSTF4LargeEQs    - a switch (0 / 1) telling if the detailed quasi-dynamic rupture model of selected earthquqakes should be written to file or not (1 == write to file)
13MinMagnitude4ST      - threshold magnitude to store the quasi-dynamic rupture model; this is included to save disk space and because smaller events might be less interesting in that regard. IMPORTANT: "smaller" in this context mainly relates to number of elements that failed in the event. The connection between number of failed elements and this minimum magnitude is therefore the mesh (i.e. fault element) size.
14/15                  - currently not implemented
