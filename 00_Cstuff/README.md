# MCQSIM COMPILATION
  MCQsim is written in C and parallelized with MPI. It further uses BLAS (provided via GSL or else). These libraries must be accessible during compilation/run-time. MCQsim compilation/installation is done in the command window i.e., terminal.

**Compilation on OSX/MAC**

> _mpicc   MCQsim_Main.c   MCQsim_StrainHS.c   MCQsim_StrainFS.c   -lm  -lgslcblas  -lmpich  -Wall  -O3  -o  MCQsim_

  this compiliation assumes use of BLAS provided by GSL (hence, GSL needs to be accessible as well).

**Compilation on Ubuntu**

 > _mpicc   MCQsim_Main.c   MCQsim_StrainHS.c   MCQsim_StrainFS.c   -lm  -lgslcblas  -lmpich  -Wall  -O3  -o  MCQsim_

  It might be necessary to check if the required modules are available/accessible. Check via:

> _module avail_

  and then load the required modules (such as GSL and/or MPI) for example

> _module load mpich_
  
 > _module load gsl_
