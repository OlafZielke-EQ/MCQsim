# MCQsim

This repository hosts MCQsim, a multi-cycle earthquake rupture simulator for generation the of earthquake sequences along predefined faults or fault systems of arbitrary geometric complexity. The updated version of this code includes QuadTree and H-matrix approximation of the stiffness kernel to improve code efficiency and scalability.

**Please check out the following publication for further details:**

> Zielke, O., and P.M. Mai (2023). MCQsim: A Multicycle Earthquake Simulator. Bull. Seismol. Soc. Am. 113(3), 889-908. doi: 10.1785/0120220248

**Contact:** olaf.zielke@kaust.edu.sa



## DOWNLOAD
  MCQsim main code and graphical user interfaces for input and ouput I/O can be found here. A simple way to download entire folders from GitHub is as follows:

  open a browser tab and go to
  
 > https://download-directory.github.io

  copy and paste the GitHub folder you want to download, such as
  
 > https://github.com/OlafZielke-EQ/MCQsim/tree/main/TREAD/
  
  the download will then start automatically







  

## MCQSIM SIMULATION
  Assuming you have all required input and parameter files defined, you can start the simulation via:

 > _mpirun -n 6 ./MCQsim ParameterFile.txt_

  This command will start MCQsim as an MPI run, using n = 6 CPUs. The parameter file contains relevant information to start/control the simulation.

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
