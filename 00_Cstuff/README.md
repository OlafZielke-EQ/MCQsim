> [!IMPORTANT]
> We are currently in the process of "renovating" the code architecture to improve efficiency/scalability etc. This also includes the input/out GUIs. To ensure that code versions are consistent, we decided to take the current version down for the time being (i.e., the time we are busy with "renovation". Please don't hesitate to contact me for further questions/comments. 
# MCQSIM COMPILATION
_MCQsim_ is written in _C_ and parallelized with _MPI_. It further uses _CBLAS_ (provided via _GSL_ or otherwise). These libraries must be accessible during compilation/run-time. _MCQsim_ compilation/installation is done in the command window i.e., terminal.


> [!IMPORTANT]
> You need to make a small change the file _MCQsim_Main.c_ depending on whether you access _CBLAS_ trough _GSL_ or ortherwise. Change the _#include_ libary accordingly
>
> _CBLAS_ is accessed trough machine installation:
> ```
> //#include <gsl/gsl_cblas.h>
> #include <cblas.h>
> ```
> compile with:
> 
> _cc   MCQsim_Main.c   MCQsim_StrainHS.c   MCQsim_StrainFS.c   -lm ~~-lgslcblas~~ ~~-lmpich~~  -Wall  -O3  -o  MCQsim_
>
> _CBLAS_ is accessed trough _GSL_:
> ```
> #include <gsl/gsl_cblas.h>
> //#include <cblas.h>
> ```
> compile with:
> 
> _mpicc   MCQsim_Main.c   MCQsim_StrainHS.c   MCQsim_StrainFS.c   -lm  ***-lgslcblas***  -lmpich  -Wall  -O3  -o  MCQsim_

For example, I compile using the first option (without _GSL_) when compiling on our HPC facilities, _Shaheen3/KAUST_. I compile using the second option (with _GSL_) when compiling on OSX or Ubuntu. Compilation creates an executable file called _"MCQsim"_.



Whether simulations are done within half-space of full space is controlled by the following statement in _MCQsim_Main.c_.

use half-space:
```
#define USEHALFSPACE            1u
```
use full space:
```
#define USEHALFSPACE            0u
```


> [!IMPORTANT]
> For LINUX, it might be necessary to check if the required modules/libraries are available and accessible. You might be able to check via:
> 
> _module avail_
>
> and then load the required modules (such as GSL and/or MPI) for example with:
> 
> _module load mpich_
>  
> _module load gsl_
> 
> Additionally: Some people reported issues/conficts between _MPICH_ and _OPENMPI_ installations. Either of these MPI libaries can be used, but compilation might give you errors if both are installed and incorrectly linked.
