# MCQSIM COMPILATION
_MCQsim_ is written in _C_ and parallelized with _MPI_. It further uses _CBLAS_ (provided via _GSL_ or otherwise). These libraries must be accessible during compilation/run-time. _MCQsim_ compilation/installation is done in the command window i.e., terminal.

***IMPORTANT:*** You need to make a small change the file _MCQsim_Main.c_ depending on whether you access _CBLAS_ trough _GSL_ or not. Change the _#include_ libary accordingly to:
```
//#include <gsl/gsl_cblas.h>
#include <cblas.h>
```
and compile with:
> _mpicc   MCQsim_Main.c   MCQsim_StrainHS.c   MCQsim_StrainFS.c   -lm  ***-lcblas***  ~~-lgslcblas~~ -lmpich  -Wall  -O3  -o  MCQsim_

if _CBLAS_ is accessed trough machine installation or:
```
#include <gsl/gsl_cblas.h>
//#include <cblas.h>
```
and compile with:
> _mpicc   MCQsim_Main.c   MCQsim_StrainHS.c   MCQsim_StrainFS.c   -lm  ~~-lcblas~~ ***-lgslcblas***  -lmpich  -Wall  -O3  -o  MCQsim_

if _CBLAS_ is accessed trough _GSL_.


***IMPORTANT:*** For LINUX, it might be necessary to check if the required modules/libraries are available and accessible. You can check via:

> _module avail_

and then load the required modules (such as GSL and/or MPI) for example with:

> _module load mpich_
>  
> _module load gsl_
