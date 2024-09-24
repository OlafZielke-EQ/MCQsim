#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <mpi.h>
#include <time.h>
#include <limits.h>

#include <gsl/gsl_cblas.h>

#define USEVPVSVELOCITY         1u
#define USEHALFSPACE            1u 
#define MINBRANCHSZEFACT        1.0f 
#define CUTDISTANCE             0.66f  
#define KH_TOL1                 0.0f 
#define KH_TOL2                 0.5f 
#define KH_FULL_DIST            15.0f
#define KH_TOL_DIST             30.0f
#define MININTSEISSTRESSCHANGE  50000.0f 

#define INTERNALFRICTION        0.75f 
#define INTERNALCOHESION        5.0E+6f 

#define MAXITERATION4BOUNDARY   200u 
#define MAXMOMRATEFUNCLENGTH    5000u 
#define NeighborFactor          2.0f 

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define MASK 123459876

#define MAX(A,B) ( ((A) >   (B))  *(A)   +  ((A) <= (B))  *(B)    )
#define MIN(A,B) ( ((A) <   (B))  *(A)   +  ((A) >= (B))  *(B)    )


extern void  StrainHS_Nikkhoo(float Stress[6], float Strain[6], float X, float Y, float Z, float P1[3], float P2[3], float P3[3], float SS, float Ds, float Ts, const float mu, const float lambda);
extern void  StrainFS_Nikkhoo(float Stress[6], float Strain[6], float X, float Y, float Z, float P1[3], float P2[3], float P3[3], float SS, float Ds, float Ts, const float mu, const float lambda);

float            ran0(long *idum);
float        FloatPow(float Base, unsigned int Exp);
void     ScaleStress6(float fStressIn[6], float fScaleFact);
void     RotStressG2L(float fOutTens[3], float RotMat[9], float fGlbStrss[6]);
void    SubtractVect3(float fDiff[3], float fMinu[3], float fSubt[3]);
void     CrossProduct(float fXProd[3], float fVect1[3], float fVect2[3]);
float    VectorLength(float fTempVect[3]);
void       NrmlzeVect(float fTempVect[3]);
void FindMainStressDir(float fSMaxVect[3], float fStrIn[6]);
void GetStkAndDipVect(float fNrm[3], float fStk[3], float fDip[3]);
void GetGlobVertsForRectangle(float fP1[3], float fP2[3], float fP3[3], float fP4[3], float *fKh_Brch, unsigned int uStart);


int main(int argc, char **argv)
{   if ( (argc  > 3 ) || (argc < 2) ) {   fprintf(stdout,"Input Error\n Please start the code in the following way:\n\n mpirun -np 4 ./MCQsimV2p5 RunParaFile.txt");   exit(10);    }
    
    unsigned int i,   j,   k;
    int iRANK,   iSIZE;
    double time_taken,   time_InterSeis = 0.0,   time_CoSeis = 0.0;
    clock_t     timer0,   timerI,   timerC;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &iRANK);
    MPI_Comm_size(MPI_COMM_WORLD, &iSIZE);
    MPI_Status  STATUS;
    MPI_Offset  OFFSETall;
    MPI_File    fp_KHMAT;
    MPI_File    fp_PREVRUN;
    MPI_File    fp_STF;
    
    double dAddedTime = 0.0,   dRecordLength = 0.0,   dMeanRecurTime = 0.0; 
    long lSeed;
    unsigned int uFPNum = 0u,   uFVNum = 0u,   uBPNum = 0u,   uBVNum = 0u,   uEQcntr = 0u,   uEQtickr = 0u, uFSegN = 0u,   uBSegN = 0u,   uMaxPNum,   uMaxSNum;
    unsigned int uRunNum = 0u,   uCatType = 0u,   uUseTimeProp = 0u,   uChgBtwEQs = 0u,   uContPrevRun = 0u,   uLoadPrev_Khmat = 0u,   uPlotCatalog2Screen = 0u,   uMinElemNum4Cat = 0u,   uLoadStep_POW2 = 0u, uStoreSTF4LargeEQs = 0u;
    float fPreStressFract = 0.0f,   fOvershootFract = 0.0f,   fMinCoSeisSlipRate = 0.0f,   fIntSeisLoadStep = 0.0f,   fMinMag4STF = 0.0f,   fMinMag4Prop = 0.0f,   fAfterSlipTime = 0.0f,   fDeepRelaxTime = 0.0f,   fHealFact = 0.0f;
    float fPSeisTime = 0.0f;
    float fFltSide,   fBndSide,   fElemArea,   fBoundArea,  fUnitSlipF,   fUnitSlipB;
    
    float fModPara[11];
    char cInputName[512],   cOutputName[512],   cKhmatFileName[512],   cPrevStateName[512],   cNameAppend[512],   cSTFName[512],   cOutputNameCUT[512];
    
    {   FILE        *fp1,                   *fp2,                       *fp3; 
        char        cFileName1[512],        cFileName2[512],            cFileName3[512];
        char        cTempVals[512];
        strcpy(cFileName1,  argv[1]);
        
        if ((fp1 = fopen(cFileName1,"r")) == NULL)      {   fprintf(stdout,"Error -cant open %s file in initializeVariables function \n",cFileName1);      exit(10);     }
        if (fgets(cTempVals, 512, fp1)    != NULL)      {                                                           }
        if (fgets(cTempVals, 512, fp1)    != NULL)      {                                                           }
        if (fgets(cTempVals, 512, fp1)    != NULL)      {   sscanf(cTempVals,"%*s %s", cInputName);                 }
        if (fgets(cTempVals, 512, fp1)    != NULL)      {   sscanf(cTempVals,"%*s %u", &uRunNum);                   }
        if (fgets(cTempVals, 512, fp1)    != NULL)      {   sscanf(cTempVals,"%*s %u", &uCatType);                  }
        if (fgets(cTempVals, 512, fp1)    != NULL)      {   sscanf(cTempVals,"%*s %u", &uPlotCatalog2Screen);       }
        if (fgets(cTempVals, 512, fp1)    != NULL)      {   sscanf(cTempVals,"%*s %ld",&lSeed);                     }
        if (fgets(cTempVals, 512, fp1)    != NULL)      {   sscanf(cTempVals,"%*s %u", &uMinElemNum4Cat);           }
        if (fgets(cTempVals, 512, fp1)    != NULL)      {   sscanf(cTempVals,"%*s %u", &uContPrevRun);              }
        if (fgets(cTempVals, 512, fp1)    != NULL)      {   sscanf(cTempVals,"%*s %u", &uLoadPrev_Khmat);           }
        if (fgets(cTempVals, 512, fp1)    != NULL)      {   sscanf(cTempVals,"%*s %s", cKhmatFileName);             }
        if (fgets(cTempVals, 512, fp1)    != NULL)      {                                                           }
        if (fgets(cTempVals, 512, fp1)    != NULL)      {   sscanf(cTempVals,"%*s %u", &uStoreSTF4LargeEQs);        }
        if (fgets(cTempVals, 512, fp1)    != NULL)      {   sscanf(cTempVals,"%*s %f", &fMinMag4STF);               }
        if (fgets(cTempVals, 512, fp1)    != NULL)      {   sscanf(cTempVals,"%*s %u", &uUseTimeProp);              }
        if (fgets(cTempVals, 512, fp1)    != NULL)      {   sscanf(cTempVals,"%*s %f", &fMinMag4Prop);              }
        if (fgets(cTempVals, 512, fp1)    != NULL)      {                                                           }
        if (fgets(cTempVals, 512, fp1)    != NULL)      {   sscanf(cTempVals,"%*s %f", &fIntSeisLoadStep);          }
        if (fgets(cTempVals, 512, fp1)    != NULL)      {   sscanf(cTempVals,"%*s %u", &uLoadStep_POW2);            }
        if (fgets(cTempVals, 512, fp1)    != NULL)      {                                                           }
        if (fgets(cTempVals, 512, fp1)    != NULL)      {   sscanf(cTempVals,"%*s %f", &fAfterSlipTime);            }
        if (fgets(cTempVals, 512, fp1)    != NULL)      {   sscanf(cTempVals,"%*s %f", &fDeepRelaxTime);            }
        if (fgets(cTempVals, 512, fp1)    != NULL)      {                                                           }
        if (fgets(cTempVals, 512, fp1)    != NULL)      {   sscanf(cTempVals,"%*s %f", &fHealFact);                 }
        if (fgets(cTempVals, 512, fp1)    != NULL)      {   sscanf(cTempVals,"%*s %f", &fOvershootFract);           }
        if (fgets(cTempVals, 512, fp1)    != NULL)      {   sscanf(cTempVals,"%*s %f", &fPreStressFract);           }
        if (fgets(cTempVals, 512, fp1)    != NULL)      {   sscanf(cTempVals,"%*s %f", &fMinCoSeisSlipRate);        }
        if (fgets(cTempVals, 512, fp1)    != NULL)      {                                                           }
        if (fgets(cTempVals, 512, fp1)    != NULL)      {   sscanf(cTempVals,"%*s %lf",&dRecordLength);             }

        uContPrevRun    = (uContPrevRun == 1u) ? 1u : 0u;
        uLoadPrev_Khmat = (uContPrevRun == 1u) ? 2u : uLoadPrev_Khmat; 
        
        fclose(fp1);
        strcpy(cFileName2, cInputName);     sprintf(cNameAppend, "_%u_Roughn.dat",uRunNum);   strcat(cFileName2,cNameAppend);
        strcpy(cFileName3, cInputName);     strcat(cFileName3,"_BNDtrig.dat"); 
        
        if ((fp2 = fopen(cFileName2,"rb"))     == NULL)     {    printf("Error -cant open *_Roughn.dat file.  in Initialize variables...   %s\n", cFileName2);      exit(10);     }
        if (fread(&uFPNum,     sizeof(unsigned int),   1, fp2) != 1)       {   exit(10);   }
        if (fread(&uFVNum,     sizeof(unsigned int),   1, fp2) != 1)       {   exit(10);   }
        if (fread(fModPara,    sizeof(float),          4, fp2) != 4)       {   exit(10);   }
        if (fread(&uChgBtwEQs, sizeof(unsigned int),   1, fp2) != 1)       {   exit(10);   }
        
        fclose(fp2);
        
        uBPNum   = 0u;               uBVNum   = 0u;
        if ((fp3 = fopen(cFileName3,"rb"))     != NULL)
        {   if (fread(&uBPNum,   sizeof(unsigned int),   1, fp3) != 1)       {   exit(10);   }
            if (fread(&uBVNum,   sizeof(unsigned int),   1, fp3) != 1)       {   exit(10);   }
            fclose(fp3);
        }
    }
    
    uMaxPNum = MAX(uFPNum, uBPNum);
    int  uFBASEelem = (int)(uFPNum/iSIZE); 
    int  uFADDelem  = (int)(uFPNum%iSIZE);
    int  uBBASEelem = (int)(uBPNum/iSIZE); 
    int  uBADDelem  = (int)(uBPNum%iSIZE);
    
    int   *iFSTART  = calloc(iSIZE, sizeof *iFSTART ); 
    int   *iBSTART  = calloc(iSIZE, sizeof *iBSTART );
    int   *iFOFFSET = calloc(iSIZE, sizeof *iFOFFSET );
    int   *iBOFFSET = calloc(iSIZE, sizeof *iBOFFSET );
    
    for (i = 0u; i < iSIZE;     i++)      {   iFOFFSET[i]     = uFBASEelem;                            }
    for (i = 0u; i < uFADDelem; i++)      {   iFOFFSET[i]    += 1u;                                    }
    for (i = 1u; i < iSIZE;     i++)      {   iFSTART[i]      = iFSTART[i-1] + iFOFFSET[i-1];          }
    
    for (i = 0u; i < iSIZE;     i++)      {   iBOFFSET[i]     = uBBASEelem;                            }
    for (i = 0u; i < uBADDelem; i++)      {   iBOFFSET[i]    += 1u;                                    }
    for (i = 1u; i < iSIZE;     i++)      {   iBSTART[i]      = iBSTART[i-1] + iBOFFSET[i-1];          }
    
    
    unsigned int *uF_temp   = calloc( 7*uFPNum,   sizeof *uF_temp );
    unsigned int *uB_temp   = calloc( 4*uBPNum,   sizeof *uB_temp );
    
    float        *fF_temp   = calloc(16*uFPNum,   sizeof *fF_temp );
    float        *fB_temp   = calloc(16*uBPNum,   sizeof *fB_temp );
    float        *fFV_temp  = calloc( 4*uFVNum,   sizeof *fFV_temp );
    float        *fBV_temp  = calloc( 3*uBVNum,   sizeof *fBV_temp );
    
    float        *fFRef     = calloc(14*iFOFFSET[iRANK],  sizeof *fFRef   );
    float        *fFFric    = calloc( 6*iFOFFSET[iRANK],  sizeof *fFFric  );
    float        *fFEvent   = calloc(17*iFOFFSET[iRANK],  sizeof *fFEvent );
    unsigned int *uFEvent   = calloc( 5*iFOFFSET[iRANK],  sizeof *uFEvent );
    unsigned int *uFActivEl = calloc( 1*iFOFFSET[iRANK],  sizeof *uFActivEl ); 
    
    float        *fBEvent   = calloc( 9*iBOFFSET[iRANK],  sizeof *fBEvent ); 
    float        *fFTempVal = calloc( 5*iFOFFSET[iRANK],  sizeof *fFTempVal );
    float        *fBTempVal = calloc( 3*iBOFFSET[iRANK],  sizeof *fBTempVal );
    float        *fMRFvals   = calloc( MAXMOMRATEFUNCLENGTH, sizeof *fMRFvals);
    
    {   unsigned int *uTempF = calloc(uFPNum, sizeof *uTempF );
        unsigned int *uTempB = calloc(uBPNum, sizeof *uTempB );
        float *fTempF        = calloc(uFPNum, sizeof *fTempF );
        float *fTempB        = calloc(uBPNum, sizeof *fTempB );
        float *fTempFv       = calloc(uFVNum, sizeof *fTempFv );
        float *fTempBv       = calloc(uBVNum, sizeof *fTempBv );
        
        unsigned int uTemp0,   uTemp1;
        float  fTemp,   fMeanHeight = 0.0f;
        float  fP1[3],   fP2[3],   fP3[3],   fP1P2[3],   fP1P3[3],   fP2P3[3],   fArea,   fP1Pc[3],   fP2Pc[3],   fTempVect[3];
        float  fNrm[3],   fStk[3],   fDip[3];
        char   cAppend[512];
        char   cFileName1[512],   cFileName2[512],   cFileName3[512];
        FILE   *fp1,   *fp2,   *fp3;
        
        
        strcpy(cFileName1,cInputName);      sprintf(cAppend, "_%u_Roughn.dat",uRunNum);     strcat(cFileName1,cAppend);
        strcpy(cFileName2,cInputName);      sprintf(cAppend, "_%u_Strgth.dat",uRunNum);     strcat(cFileName2,cAppend);
        strcpy(cFileName3,cInputName);      strcat(cFileName3,"_BNDtrig.dat");
        
        if ((fp1 = fopen(cFileName1,"rb"))     == NULL)     {   printf("Error -cant open %s LoadInputParameter function...\n",cFileName1);      exit(10);     }
        
        fseek(fp1, (4*sizeof(float)+3*sizeof(unsigned int)), SEEK_SET); 
        
        if (fread(uTempF, sizeof(unsigned int), uFPNum, fp1) != uFPNum)       {   exit(10); }
        for (i = 0u; i < uFPNum; i++)       {   uF_temp[i*7 +0] = uTempF[i] -1u;            }
        if (fread(uTempF, sizeof(unsigned int), uFPNum, fp1) != uFPNum)       {   exit(10); }
        for (i = 0u; i < uFPNum; i++)       {   uF_temp[i*7 +1] = uTempF[i] -1u;            }
        if (fread(uTempF, sizeof(unsigned int), uFPNum, fp1) != uFPNum)       {   exit(10); }
        for (i = 0u; i < uFPNum; i++)       {   uF_temp[i*7 +2] = uTempF[i] -1u;            }
        if (fread(uTempF, sizeof(unsigned int), uFPNum, fp1) != uFPNum)       {   exit(10); }
        for (i = 0u; i < uFPNum; i++)       {   uF_temp[i*7 +3] = uTempF[i] -1u;            }
        if (fread(uTempF, sizeof(unsigned int), uFPNum, fp1) != uFPNum)       {   exit(10); }
        for (i = 0u; i < uFPNum; i++)       {   uF_temp[i*7 +4] = uTempF[i] -1u;            }
        if (fread(uTempF, sizeof(unsigned int), uFPNum,fp1) != uFPNum)        {   exit(10); }
        for (i = 0u; i < uFPNum; i++)       {   uF_temp[i*7 +5] = uTempF[i];                }
        
        if (fread(fTempF, sizeof(float), uFPNum,fp1) != uFPNum)      {   exit(10);          }
        for (i = 0u; i < uFPNum; i++)       {   fF_temp[i*16 +4] = fTempF[i];               }
        if (fread(fTempF, sizeof(float), uFPNum,fp1) != uFPNum)      {   exit(10);          }
        for (i = 0u; i < uFPNum; i++)       {   fF_temp[i*16 +5] = fTempF[i];               }
        
        if (fread(fTempF, sizeof(float), uFPNum,fp1) != uFPNum)     {   exit(10);       }
        for (i = 0u; i < uFPNum; i++)       {   fF_temp[i*16 +0] = fTempF[i]*1000.0f;   }
        if (fread(fTempF, sizeof(float), uFPNum,fp1) != uFPNum)     {   exit(10);       }
        for (i = 0u; i < uFPNum; i++)       {   fF_temp[i*16 +1] = fTempF[i]*1000.0f;   }
        if (fread(fTempF, sizeof(float), uFPNum,fp1) != uFPNum)     {   exit(10);       }
        for (i = 0u; i < uFPNum; i++)       {   fF_temp[i*16 +2] = fTempF[i]*1000.0f;   }
        
        if (fread(fTempFv, sizeof(float), uFVNum,fp1) != uFVNum)      {   exit(10);     }
        for (i = 0u; i < uFVNum; i++)       {   fFV_temp[i*4 +0] = fTempFv[i]*1000.0f;  }
        if (fread(fTempFv, sizeof(float), uFVNum,fp1) != uFVNum)      {   exit(10);     }
        for (i = 0u; i < uFVNum; i++)       {   fFV_temp[i*4 +1] = fTempFv[i]*1000.0f;  }
        if (fread(fTempFv, sizeof(float), uFVNum,fp1) != uFVNum)      {   exit(10);     }
        for (i = 0u; i < uFVNum; i++)       {   fFV_temp[i*4 +2] = fTempFv[i]*1000.0f;  }
        if (fread(fTempFv, sizeof(float), uFVNum,fp1) != uFVNum)      {   exit(10);     }
        for (i = 0u; i < uFVNum; i++)       {   fFV_temp[i*4 +3] = fTempFv[i]*1000.0f;  }
        
        fclose(fp1);
        
        if ((fp2 = fopen(cFileName2,"rb"))     == NULL)     {   printf("Error -cant open %s LoadInputParameter function...\n",cFileName2);      exit(10);     }
        
        fseek(fp2, sizeof(unsigned int), SEEK_SET); 
    
        if (fread(fTempF, sizeof(float), uFPNum,fp2) != uFPNum) {   exit(10);           }
        for (i = iFSTART[iRANK]; i < (iFSTART[iRANK]+iFOFFSET[iRANK]); i++)       {   fFRef[(i-iFSTART[iRANK])*14 +3] = fTempF[i];          }
        if (fread(fTempF, sizeof(float), uFPNum,fp2) != uFPNum) {   exit(10);           }
        for (i = iFSTART[iRANK]; i < (iFSTART[iRANK]+iFOFFSET[iRANK]); i++)       {   fFRef[(i-iFSTART[iRANK])*14 +5] = fTempF[i];          }
        if (fread(fTempF, sizeof(float), uFPNum,fp2) != uFPNum) {   exit(10);           }
        for (i = iFSTART[iRANK]; i < (iFSTART[iRANK]+iFOFFSET[iRANK]); i++)       {   fFRef[(i-iFSTART[iRANK])*14 +1] = fTempF[i];          }
        if (fread(fTempF, sizeof(float), uFPNum,fp2) != uFPNum) {   exit(10);           }
        for (i = iFSTART[iRANK]; i < (iFSTART[iRANK]+iFOFFSET[iRANK]); i++)       {   fFRef[(i-iFSTART[iRANK])*14 +7] = fTempF[i];          }
        if (fread(fTempF, sizeof(float), uFPNum,fp2) != uFPNum) {   exit(10);           }
        for (i = iFSTART[iRANK]; i < (iFSTART[iRANK]+iFOFFSET[iRANK]); i++)       {   fFRef[(i-iFSTART[iRANK])*14 +4] = fTempF[i]*0.01f*fFRef[(i-iFSTART[iRANK])*14 +3];          }
        if (fread(fTempF, sizeof(float), uFPNum,fp2) != uFPNum) {   exit(10);           }
        for (i = iFSTART[iRANK]; i < (iFSTART[iRANK]+iFOFFSET[iRANK]); i++)       {   fFRef[(i-iFSTART[iRANK])*14 +6] = fTempF[i]*0.01f*fFRef[(i-iFSTART[iRANK])*14 +5];          }
        if (fread(fTempF, sizeof(float), uFPNum,fp2) != uFPNum) {   exit(10);           }
        for (i = iFSTART[iRANK]; i < (iFSTART[iRANK]+iFOFFSET[iRANK]); i++)       {   fFRef[(i-iFSTART[iRANK])*14 +8] = fTempF[i]*0.01f*fFRef[(i-iFSTART[iRANK])*14 +7];          }
        if (fread(fTempF, sizeof(float), uFPNum,fp2) != uFPNum) {   exit(10);           }
        for (i = iFSTART[iRANK]; i < (iFSTART[iRANK]+iFOFFSET[iRANK]); i++)       {   fFRef[(i-iFSTART[iRANK])*14 +3] += fTempF[i];          } 
        if (fread(fTempF, sizeof(float), uFPNum,fp2) != uFPNum) {   exit(10);           }
        for (i = iFSTART[iRANK]; i < (iFSTART[iRANK]+iFOFFSET[iRANK]); i++)       {   fFRef[(i-iFSTART[iRANK])*14 +5] += fTempF[i];          } 
        if (fread(fTempF, sizeof(float), uFPNum,fp2) != uFPNum) {   exit(10);           }
        for (i = iFSTART[iRANK]; i < (iFSTART[iRANK]+iFOFFSET[iRANK]); i++)       {   fFRef[(i-iFSTART[iRANK])*14 +1] += fTempF[i];          } 
        if (fread(fTempF, sizeof(float), uFPNum,fp2) != uFPNum) {   exit(10);           }
        for (i = iFSTART[iRANK]; i < (iFSTART[iRANK]+iFOFFSET[iRANK]); i++)       {   fFRef[(i-iFSTART[iRANK])*14 +7] += fTempF[i];          } 
        
        fclose(fp2);
        
        if ((fp3 = fopen(cFileName3,"rb")) != NULL)
        {   
            fseek(fp3, (2*sizeof(unsigned int)), SEEK_SET);
            
            if (fread(uTempB,    sizeof(unsigned int), uBPNum,fp3) != uBPNum)  {   exit(10); }
            for (i = 0u; i < uBPNum; i++)       {   uB_temp[i*4 +0] = uTempB[i] -1u;         }
            if (fread(uTempB,    sizeof(unsigned int), uBPNum,fp3) != uBPNum)  {   exit(10); }
            for (i = 0u; i < uBPNum; i++)       {   uB_temp[i*4 +1] = uTempB[i] -1u;         }
            if (fread(uTempB,    sizeof(unsigned int), uBPNum,fp3) != uBPNum)  {   exit(10); }
            for (i = 0u; i < uBPNum; i++)       {   uB_temp[i*4 +2] = uTempB[i] -1u;         }
            if (fread(uTempB,    sizeof(unsigned int), uBPNum,fp3) != uBPNum)  {   exit(10); }
            for (i = 0u; i < uBPNum; i++)       {   uB_temp[i*4 +3] = uTempB[i] -1u;         }
            
            if (fread(fTempB,    sizeof(float), uBPNum,fp3) != uBPNum)  {   exit(10);       }
            for (i = 0u; i < uBPNum; i++)       {   fB_temp[i*16 +0] = fTempB[i]*1000.0f;   }
            if (fread(fTempB,    sizeof(float), uBPNum,fp3) != uBPNum)  {   exit(10);       }
            for (i = 0u; i < uBPNum; i++)       {   fB_temp[i*16 +1] = fTempB[i]*1000.0f;   }
            if (fread(fTempB,    sizeof(float), uBPNum,fp3) != uBPNum)  {   exit(10);       }
            for (i = 0u; i < uBPNum; i++)       {   fB_temp[i*16 +2] = fTempB[i]*1000.0f;   }
            
            if (fread(fTempBv,    sizeof(float), uBVNum,fp3) != uBVNum)  {   exit(10);      }
            for (i = 0u; i < uBVNum; i++)       {   fBV_temp[i*3 +0] = fTempBv[i]*1000.0f;  }
            if (fread(fTempBv,    sizeof(float), uBVNum,fp3) != uBVNum)  {   exit(10);      }
            for (i = 0u; i < uBVNum; i++)       {   fBV_temp[i*3 +1] = fTempBv[i]*1000.0f;  }
            if (fread(fTempBv,    sizeof(float), uBVNum,fp3) != uBVNum)  {   exit(10);      }
            for (i = 0u; i < uBVNum; i++)       {   fBV_temp[i*3 +2] = fTempBv[i]*1000.0f;  }
            
            if (fread(fTempB,    sizeof(float), uBPNum,fp3) != uBPNum)  {   exit(10);       }
            for (i = 0u; i < uBPNum; i++)       {   fB_temp[i*16 +3] = fTempB[i]*0.001f;   }
            if (fread(fTempB,    sizeof(float), uBPNum,fp3) != uBPNum)  {   exit(10);       }
            for (i = 0u; i < uBPNum; i++)       {   fB_temp[i*16 +4] = fTempB[i]*0.001f;   }
            if (fread(fTempB,    sizeof(float), uBPNum,fp3) != uBPNum)  {   exit(10);       }
            for (i = 0u; i < uBPNum; i++)       {   fB_temp[i*16 +5] = fTempB[i]*0.001f;   }
            
            fclose(fp3);
        }
        else    {   fprintf(stdout,"\nWarning -cant open  %s; Continue without boundary surface\n",cFileName3);     }
        
        fFltSide = 0.0f;   fBndSide = 0.0f;   fElemArea = 0.0f,   fBoundArea = 0.0f;
        
        for (i = 0u; i < uFPNum; i++)
        {   
            memcpy(fP1, &fFV_temp[uF_temp[i*7 +0]*4 +0],  3u*sizeof(float));
            memcpy(fP2, &fFV_temp[uF_temp[i*7 +1]*4 +0],  3u*sizeof(float));
            memcpy(fP3, &fFV_temp[uF_temp[i*7 +2]*4 +0],  3u*sizeof(float));
            uFSegN = MAX(uFSegN, uF_temp[i*7 +3]);
            
            SubtractVect3(fP1P2, fP2, fP1);                     SubtractVect3(fP1P3, fP3, fP1);                         SubtractVect3(fP2P3, fP3, fP2);
            SubtractVect3(fP1Pc, &fF_temp[i*16 +0], fP1);       SubtractVect3(fP2Pc, &fF_temp[i*16 +0], fP2);
            
            CrossProduct(fTempVect,fP1P2, fP1Pc);               fArea = 0.5f*VectorLength(fTempVect);                   fMeanHeight += ((2.0f*fArea)/VectorLength(fP1P2));
            CrossProduct(fTempVect,fP1P3, fP1Pc);               fArea = 0.5f*VectorLength(fTempVect);                   fMeanHeight += ((2.0f*fArea)/VectorLength(fP1P3));
            CrossProduct(fTempVect,fP2P3, fP2Pc);               fArea = 0.5f*VectorLength(fTempVect);                   fMeanHeight += ((2.0f*fArea)/VectorLength(fP2P3));
            
            CrossProduct(fNrm,fP1P2, fP1P3);                    fF_temp[i*16 +15] = 0.5f*VectorLength(fNrm); 
            
            fElemArea += fF_temp[i*16 +15];
            fFltSide  += (VectorLength(fP1P2) + VectorLength(fP1P3) +VectorLength(fP2P3))/3.0f;
            
            NrmlzeVect(fNrm);
            
            uTemp0 = (fNrm[2] < 0.0f)*1u + (fNrm[2] >= 0.0f)*0u;
            uTemp1 = uF_temp[i*7 +1];
            uF_temp[i*7 +1] = (uTemp0 == 0u)*uF_temp[i*7 +1] + (uTemp0 != 0u)*uF_temp[i*7 +2];
            uF_temp[i*7 +2] = (uTemp0 == 0u)*uF_temp[i*7 +2] + (uTemp0 != 0u)*uTemp1;
            
            GetStkAndDipVect(fNrm, fStk, fDip);
            
            memcpy(&fF_temp[i*16 +6], fNrm,  3u*sizeof(float));
            memcpy(&fF_temp[i*16 +9], fStk,  3u*sizeof(float));
            memcpy(&fF_temp[i*16+12], fDip,  3u*sizeof(float));
            
        }
        uFSegN   += 1u;
        
        fFltSide    = (uFPNum > 0u) ? fFltSide/(float)uFPNum           : fFltSide;
        fElemArea   = (uFPNum > 0u) ? fElemArea/(float)uFPNum          : fElemArea;
        fMeanHeight = (uFPNum > 0u) ?(fMeanHeight/((float)uFPNum*3.0f)): fMeanHeight; 
        
        for (i = 0u; i < uBPNum; i++)
        {   memcpy(fP1, &fBV_temp[uB_temp[i*4 +0]*3 +0],  3u*sizeof(float));
            memcpy(fP2, &fBV_temp[uB_temp[i*4 +1]*3 +0],  3u*sizeof(float));
            memcpy(fP3, &fBV_temp[uB_temp[i*4 +2]*3 +0],  3u*sizeof(float));
            uBSegN = MAX(uBSegN, uB_temp[i*4 +3]);
            
            SubtractVect3(fP1P2, fP2, fP1);                     SubtractVect3(fP1P3, fP3, fP1);                         SubtractVect3(fP2P3, fP3, fP2);
            
            CrossProduct(fNrm,fP1P2, fP1P3);                    fB_temp[i*16 +15] = 0.5f*VectorLength(fNrm); 
            
            fBoundArea  += fB_temp[i*16 +15];
            fBndSide    += (VectorLength(fP1P2) + VectorLength(fP1P3) +VectorLength(fP2P3))/3.0f;
            
            NrmlzeVect(fNrm);
            
            uTemp0 = (fNrm[2] < 0.0f)*1u + (fNrm[2] >= 0.0f)*0u;
            uTemp1 = uB_temp[i*4 +1];
            uB_temp[i*4 +1] = (uTemp0 == 0u)*uB_temp[i*4 +1] + (uTemp0 != 0u)*uF_temp[i*4 +2];
            uB_temp[i*4 +2] = (uTemp0 == 0u)*uB_temp[i*4 +2] + (uTemp0 != 0u)*uTemp1;
            
            GetStkAndDipVect(fNrm, fStk, fDip);
            
            memcpy(&fB_temp[i*16 +6], fNrm,  3u*sizeof(float));
            memcpy(&fB_temp[i*16 +9], fStk,  3u*sizeof(float));
            memcpy(&fB_temp[i*16+12], fDip,  3u*sizeof(float));
            
        }
        uBSegN   += 1u;
        
        fBndSide    = (uBPNum > 0u) ? fBndSide/(float)uBPNum   : fBndSide;
        fBoundArea  = (uBPNum > 0u) ? fBoundArea/(float)uBPNum : fBoundArea;
        
        uMaxSNum   = MAX(uFSegN, uBSegN);
        
        fUnitSlipF = fFltSide *1.0E-5f; 
        fUnitSlipB = (fBoundArea/fElemArea)*fUnitSlipF;
        
        fModPara[2] = fModPara[2]*1.0E+9f; 
        fModPara[4] =(2.0f*fModPara[2]*fModPara[3])/(1.0f-2.0f*fModPara[3]); 
        fTemp       = (fModPara[0] > 0.0f)*fModPara[0] + (fModPara[0] <= 0.0f)*2700.0f; 
        fModPara[5] = sqrtf((fModPara[4] +2.0f*fModPara[2])/fTemp); 
        fModPara[6] = (USEVPVSVELOCITY == 1u)*(fModPara[5]/sqrtf(fModPara[2]/fTemp)) + (USEVPVSVELOCITY != 1u)*1.0f; 
        fModPara[7] = (2.0f*fMeanHeight)/fModPara[5];
        fModPara[8] = (uUseTimeProp == 1u)*fModPara[7] + (uUseTimeProp != 1u)*FLT_MAX;
        fModPara[9] = sqrtf(fModPara[2]/fTemp); 
        fModPara[10]= fMinCoSeisSlipRate*fModPara[7];
        
        if (iRANK == 0) 
        {   fprintf(stdout,"----------------------\n");
            fprintf(stdout,"Number of RANKS: %u\n",iSIZE);
            fprintf(stdout,"----------------------\n");
            fprintf(stdout,"System info: Byte Size for FLOAT: %lu     unsigned INT: %lu    \n", sizeof(float), sizeof(unsigned int));
            fprintf(stdout,"----------------------\n");
            fprintf(stdout,"FileName:               %s\n",cInputName);                      fprintf(stdout,"RunNumber:              %u\n",uRunNum);
            fprintf(stdout,"CatalogType:            %u\n",uCatType);                        fprintf(stdout,"PlotCat2Screen:         %u\n",uPlotCatalog2Screen);
            fprintf(stdout,"SeedLocation:           %ld\n",lSeed);                          fprintf(stdout,"MinElem4Catalog:        %u\n",uMinElemNum4Cat);
            fprintf(stdout,"ContPreviousRun:        %u\n",uContPrevRun);                    fprintf(stdout,"LoadPrev_Kh_mat:        %u\n",uLoadPrev_Khmat);
            fprintf(stdout,"Kh_mat file name:       %s\n",cKhmatFileName);                  
            fprintf(stdout,"----------------------\n");  
            fprintf(stdout,"SaveSTF4LargeEQ:        %u\n",uStoreSTF4LargeEQs);              fprintf(stdout,"MinMag2SaveSTF:         %f\n",fMinMag4STF); 
            fprintf(stdout,"UseTime4RuptProp:       THIS IS NOT USED/INCLUDED IN V2p5\n");  fprintf(stdout,"MinMag2UseRuptProp:     THIS IS NOT USED/INCLUDED IN V2p5\n");
            fprintf(stdout,"UseHalfSpace:           %u\n",USEHALFSPACE);
            fprintf(stdout,"----------------------\n"); 
            fprintf(stdout,"IntSeisTimeStep:        %fdays\n",fIntSeisLoadStep);            fprintf(stdout,"LoadSteps_POW2:         %u\n",uLoadStep_POW2);
            fprintf(stdout,"----------------------\n"); 
            fprintf(stdout,"ViscAftSlipTime:        %fyears\n",fAfterSlipTime);             fprintf(stdout,"ViscDeepRelaxTime:      %fyears\n",fDeepRelaxTime);
            fprintf(stdout,"----------------------\n"); 
            fprintf(stdout,"CoSeisHealFract:        %f\n",fHealFact);                       fprintf(stdout,"OvershootFract:         %f\n",fOvershootFract);
            fprintf(stdout,"PreStressFract:         %f\n",fPreStressFract);                 fprintf(stdout,"MinCoSeisSlipRate(m/s): %f\n",fMinCoSeisSlipRate);
            fprintf(stdout,"----------------------\n"); 
            fprintf(stdout,"RecLength:              %lf\n",dRecordLength);               
            fprintf(stdout,"----------------------\n"); 
            fprintf(stdout,"Elastic properties\n");
            fprintf(stdout,"Medium density (kg/m^3):   %e\n",fModPara[0]);              fprintf(stdout,"AddedNormalStress (MPa):   %e\n",fModPara[1]);
            fprintf(stdout,"ShearModulus (Pa):         %e\n",fModPara[2]);              fprintf(stdout,"PoissonRatio:              %e\n",fModPara[3]);
            fprintf(stdout,"ChangeFricBtwEQs:          %u\n",uChgBtwEQs);               fprintf(stdout,"Lambda (Pa):               %e\n",fModPara[4]);
            fprintf(stdout,"P-waveVelocity (m/s):      %e\n",fModPara[5]);              fprintf(stdout,"vP/vS ratio:               %e\n",fModPara[6]);
            fprintf(stdout,"unit slip (m, fault):      %e\n",fUnitSlipF);               fprintf(stdout,"unit slip (m, bound):      %e\n",fUnitSlipB);
            fprintf(stdout,"min. slip 2 start EQ (m):  %e\n",fModPara[10]); 
            fprintf(stdout,"coseis. timestep (s):      %e\n",fModPara[7]);              fprintf(stdout,"mean height (m):           %e\n",fMeanHeight);  
            
            float fTemp;
            fTemp = 100.0f - expf(-1.0f/fAfterSlipTime)*100.0f; 
            fprintf(stdout,"\nFractional post-seismic change during first year of after-slip: %2.4f percent released \n",fTemp);
            fTemp = 100.0f - expf(-1.0f/fDeepRelaxTime)*100.0f; 
            fprintf(stdout,"Fractional post-seismic change during first year of 'deep relaxation': %2.4f percent released \n",fTemp);
            
            fprintf(stdout,"FaultPatchNumber %d     FaultVertexNumber %d\n", uFPNum, uFVNum);
            fprintf(stdout,"BoundPatchNumber %d     BoundVertexNumber %d\n", uBPNum, uBVNum);
            fprintf(stdout,"FaultSideLenth and ElementArea: %5.2fm  and %4.3fkm^2;   BoundarySideLength and ElementArea: %5.2fm  and %4.3fkm^2\n\n",   fFltSide, fElemArea*1.0E-6f, fBndSide,fBoundArea*1.0E-6f);
        } 
    } 
    
    lSeed += (long)iRANK;
    
    {   unsigned int uGlobPos;
        float fGlbStrss[6], fOutTens[6],   fNrm[3];
        
        for (i = 0u; i < iFOFFSET[iRANK]; i++) 
        {   uGlobPos   = i + iFSTART[iRANK];
            
            fFRef[i*14 +0]   = fF_temp[uGlobPos*16 +15];
            fFRef[i*14 +1]   = fModPara[0] * 9.81f * fF_temp[uGlobPos*16 +2] +(fModPara[1]*1.0E+6f); 
            fFEvent[i*17 +2] = fFRef[i*14 +1];
            fFRef[i*14 +9]   = fF_temp[uGlobPos*16 +0];
            fFRef[i*14 +10]  = fF_temp[uGlobPos*16 +1];
            fFRef[i*14 +11]  = fF_temp[uGlobPos*16 +2];
            
            if (uF_temp[uGlobPos*7 +5] == 1u)
            {   fFEvent[i*17 +8]  = 0.0f;
                fFEvent[i*17 +9]  = 0.0f;
                fFEvent[i*17 +13] = 1.0E-3f*fF_temp[uGlobPos*16 +4] *cosf(fF_temp[uGlobPos*16 +5]*(M_PI/180.0f));
                fFEvent[i*17 +14] = 1.0E-3f*fF_temp[uGlobPos*16 +4] *sinf(fF_temp[uGlobPos*16 +5]*(M_PI/180.0f));
                
            } 
            else
            {   fNrm[0]      = sinf(fF_temp[uGlobPos*16 +5]*(M_PI/180.0f));
                fNrm[1]      = cosf(fF_temp[uGlobPos*16 +5]*(M_PI/180.0f));
                fNrm[2]      = 0.0f;
                fGlbStrss[0] = 1.0E+6f*fF_temp[uGlobPos*16 +4] *fNrm[0]*fNrm[0]; 
                fGlbStrss[1] = 1.0E+6f*fF_temp[uGlobPos*16 +4] *fNrm[0]*fNrm[1]; 
                fGlbStrss[2] = 1.0E+6f*fF_temp[uGlobPos*16 +4] *fNrm[0]*fNrm[2]; 
                fGlbStrss[3] = 1.0E+6f*fF_temp[uGlobPos*16 +4] *fNrm[1]*fNrm[1]; 
                fGlbStrss[4] = 1.0E+6f*fF_temp[uGlobPos*16 +4] *fNrm[1]*fNrm[2]; 
                fGlbStrss[5] = 1.0E+6f*fF_temp[uGlobPos*16 +4] *fNrm[2]*fNrm[2]; 
                
                RotStressG2L(fOutTens, &fF_temp[uGlobPos*16 +6], fGlbStrss);
                
                fFEvent[i*17 +8]  = fOutTens[0]; 
                fFEvent[i*17 +9]  = fOutTens[1]; 
                fFEvent[i*17 +13] = 0.0f;
                fFEvent[i*17 +14] = 0.0f;
    }   }   }
    
    
    unsigned int *uKh_FFcnt  = calloc(iFOFFSET[iRANK], sizeof *uKh_FFcnt );
    unsigned int **uKh_FFps2 = NULL;
    float **fKh_FFvalstk = NULL;
    float **fKh_FFvaldip = NULL;
    float **fKh_FFvalnrm = NULL;
    uKh_FFps2    = calloc(iFOFFSET[iRANK], sizeof *uKh_FFps2 );
    fKh_FFvalstk = calloc(iFOFFSET[iRANK], sizeof *fKh_FFvalstk );
    fKh_FFvaldip = calloc(iFOFFSET[iRANK], sizeof *fKh_FFvaldip );
    fKh_FFvalnrm = calloc(iFOFFSET[iRANK], sizeof *fKh_FFvalnrm );
    
    unsigned int *uKh_FBcnt = calloc(iFOFFSET[iRANK], sizeof *uKh_FBcnt );
    unsigned int **uKh_FBps2 = NULL;
    float **fKh_FBvalstk = NULL;
    float **fKh_FBvaldip = NULL;
    uKh_FBps2    = calloc(iFOFFSET[iRANK], sizeof *uKh_FBps2 );
    fKh_FBvalstk = calloc(iFOFFSET[iRANK], sizeof *fKh_FBvalstk );
    fKh_FBvaldip = calloc(iFOFFSET[iRANK], sizeof *fKh_FBvaldip );
    
    unsigned int *uKh_BBcnt = calloc(iBOFFSET[iRANK], sizeof *uKh_BBcnt );
    unsigned int **uKh_BBps2 = NULL;
    float **fKh_BBvalstk = NULL;
    float **fKh_BBvaldip = NULL;
    float **fKh_BBvalnrm = NULL;
    uKh_BBps2    = calloc(iBOFFSET[iRANK], sizeof *uKh_BBps2 );
    fKh_BBvalstk = calloc(iBOFFSET[iRANK], sizeof *fKh_BBvalstk );
    fKh_BBvaldip = calloc(iBOFFSET[iRANK], sizeof *fKh_BBvaldip );
    fKh_BBvalnrm = calloc(iBOFFSET[iRANK], sizeof *fKh_BBvalnrm );
    
    unsigned int *uKh_BFcnt = calloc(iBOFFSET[iRANK], sizeof *uKh_BFcnt );
    unsigned int **uKh_BFps2 = NULL;
    float **fKh_BFvalstk = NULL;
    float **fKh_BFvaldip = NULL;
    float **fKh_BFvalnrm = NULL;
    uKh_BFps2    = calloc(iBOFFSET[iRANK], sizeof *uKh_BFps2 );
    fKh_BFvalstk = calloc(iBOFFSET[iRANK], sizeof *fKh_BFvalstk );
    fKh_BFvaldip = calloc(iBOFFSET[iRANK], sizeof *fKh_BFvaldip );
    fKh_BFvalnrm = calloc(iBOFFSET[iRANK], sizeof *fKh_BFvalnrm );
    
    
    if ((uLoadPrev_Khmat == 0) || (uLoadPrev_Khmat == 1)) 
    {   
        unsigned int u_MaxBrLvl    = 20u; 
        unsigned int uF_TotalOffs   = 0u; 
        unsigned int uB_TotalOffs   = 0u; 
        unsigned int *uF_OTElem    = calloc((uFPNum*u_MaxBrLvl), sizeof *uF_OTElem );
        unsigned int *uB_OTElem    = calloc((uBPNum*u_MaxBrLvl), sizeof *uB_OTElem );
        unsigned int *uF_OTPrtChld = calloc((4*uFPNum*5),       sizeof *uF_OTPrtChld ); 
        unsigned int *uB_OTPrtChld = calloc((4*uBPNum*5),       sizeof *uB_OTPrtChld );
        
        {   timer0 = clock();
            unsigned int *u_OTElem_temp = calloc((uMaxPNum*u_MaxBrLvl), sizeof *u_OTElem_temp ); 
            unsigned int *u_BrIds        = calloc(uMaxPNum,             sizeof *u_BrIds  ); 
            unsigned int *u_nBrIds       = calloc(uMaxPNum,             sizeof *u_nBrIds ); 
            float *f_BrLims              = calloc(uMaxPNum*4,           sizeof *f_BrLims ); 
            float *f_nBrLims             = calloc(uMaxPNum*4,           sizeof *f_nBrLims ); 
            
            unsigned int *u_sElem        = calloc(uMaxPNum,             sizeof *u_sElem ); 
            unsigned int *u_newBr        = calloc(4*2,                  sizeof *u_newBr ); 
            
            unsigned int *u_ElemCntr     = calloc(       1*uMaxSNum,    sizeof *u_ElemCntr);
            unsigned int *u_ElSegIDs     = calloc(uMaxPNum*uMaxSNum,    sizeof *u_ElSegIDs);
            float *f_SegVals             = calloc(       6*uMaxSNum,    sizeof *f_SegVals);
            float *fCentRot              = calloc(uMaxPNum*3,           sizeof *fCentRot);
            float fNrm[3],   fStk[3],   fDip[3],   fShftCent[3],   fMaxMinVals[6];
            
            unsigned int uBrNum,   uBrLevel,   unBrNum,   uElem,   uNxtSegID,   uTemp,   uBrCount;
            float fTemp,   fBrSize;
            
            for (i = 0u; i < uFPNum; i++)
            {   uNxtSegID = uF_temp[i*7 +3];
                u_ElSegIDs[uNxtSegID*uMaxPNum + u_ElemCntr[uNxtSegID]] = i;
                u_ElemCntr[uNxtSegID] += 1u;
                
                f_SegVals[uNxtSegID*6 +0] += fF_temp[i*16 +0];
                f_SegVals[uNxtSegID*6 +1] += fF_temp[i*16 +1];
                f_SegVals[uNxtSegID*6 +2] += fF_temp[i*16 +2];
                
                fTemp = (fF_temp[i*16 +8] >= 0.0f)*1.0f + (fF_temp[i*16 +8] < 0.0f)*-1.0f;
                f_SegVals[uNxtSegID*6 +3] += (fF_temp[i*16 +6]*fTemp);
                f_SegVals[uNxtSegID*6 +4] += (fF_temp[i*16 +7]*fTemp);
                f_SegVals[uNxtSegID*6 +5] += (fF_temp[i*16 +8]*fTemp);
            }
            
            for (i = 0u; i < uFSegN; i++)
            {   u_OTElem_temp = realloc(u_OTElem_temp,    u_ElemCntr[i]*u_MaxBrLvl* sizeof *u_OTElem_temp);
                memset(u_OTElem_temp, 0, u_ElemCntr[i]*u_MaxBrLvl* sizeof(unsigned int) );
                memset(u_BrIds,       0, u_ElemCntr[i]*            sizeof(unsigned int) );
                
                f_SegVals[i*6 +0] /= (float)u_ElemCntr[i];        f_SegVals[i*6 +1] /= (float)u_ElemCntr[i];        f_SegVals[i*6 +2] /= (float)u_ElemCntr[i];
                f_SegVals[i*6 +3] /= (float)u_ElemCntr[i];        f_SegVals[i*6 +4] /= (float)u_ElemCntr[i];        f_SegVals[i*6 +5] /= (float)u_ElemCntr[i];
                memcpy(fNrm, &f_SegVals[i*6 +3], 3*sizeof(float));
                GetStkAndDipVect(fNrm, fStk, fDip);
                
                fTemp = -FLT_MAX;
                
                for (j = 0; j < u_ElemCntr[i]; j++)
                {   fShftCent[0] = fF_temp[ u_ElSegIDs[i*uMaxPNum +j]*16 +0] - f_SegVals[i*6 +0];
                    fShftCent[1] = fF_temp[ u_ElSegIDs[i*uMaxPNum +j]*16 +1] - f_SegVals[i*6 +1];
                    fShftCent[2] = fF_temp[ u_ElSegIDs[i*uMaxPNum +j]*16 +2] - f_SegVals[i*6 +2];
                    
                    fCentRot[j*3 +0] = fNrm[0]*fShftCent[0] + fNrm[1]*fShftCent[1] + fNrm[2]*fShftCent[2];
                    fCentRot[j*3 +1] = fStk[0]*fShftCent[0] + fStk[1]*fShftCent[1] + fStk[2]*fShftCent[2];
                    fCentRot[j*3 +2] = fDip[0]*fShftCent[0] + fDip[1]*fShftCent[1] + fDip[2]*fShftCent[2];
                    fTemp = MAX(fTemp, fCentRot[j*3 +2]);
                }
                for (j = 0; j < u_ElemCntr[i]; j++)     {   fCentRot[j*3 +2] -= fTemp;          }
                fMaxMinVals[0] = fCentRot[0*3 +1];                                fMaxMinVals[2] = fCentRot[0*3 +1];
                fMaxMinVals[3] = fCentRot[0*3 +2];                                fMaxMinVals[5] = fCentRot[0*3 +2];
                
                for (j = 1; j < u_ElemCntr[i]; j++)
                {   fMaxMinVals[0] = MIN(fMaxMinVals[0], fCentRot[j*3 +1]);       fMaxMinVals[2] = MAX(fMaxMinVals[2], fCentRot[j*3 +1]);
                    fMaxMinVals[3] = MIN(fMaxMinVals[3], fCentRot[j*3 +2]);       fMaxMinVals[5] = MAX(fMaxMinVals[5], fCentRot[j*3 +2]);
                }
                
                fMaxMinVals[0]  -= fFltSide;                    fMaxMinVals[2] += fFltSide;
                fMaxMinVals[3]  -= fFltSide;                    fMaxMinVals[5] += fFltSide;
                fBrSize          = MAX((fMaxMinVals[2] -fMaxMinVals[0]), (fMaxMinVals[5] -fMaxMinVals[3]));
                f_BrLims[0*4 +0] = -0.5f*fBrSize;              f_BrLims[0*4 +1] = 0.5f*fBrSize;
                f_BrLims[0*4 +2] = -1.0f*fBrSize;              f_BrLims[0*4 +3] = 0.0f*fBrSize;
                
                uBrLevel = 0u;
                uBrCount = 0u;
                uBrNum   = 1u;
                
                while (0.5f*fBrSize >= MINBRANCHSZEFACT*fFltSide) 
                {   uBrLevel += 1u; 
                    fBrSize  *= 0.5f; 
                    unBrNum   = 0u; 
                    
                    for (j = 0u; j < uBrNum; j++)
                    {   fMaxMinVals[0] = f_BrLims[j*4 +0];      fMaxMinVals[2] = f_BrLims[j*4 +1];       fMaxMinVals[1] = 0.5f*(fMaxMinVals[0] +fMaxMinVals[2]);
                        fMaxMinVals[3] = f_BrLims[j*4 +2];      fMaxMinVals[5] = f_BrLims[j*4 +3];       fMaxMinVals[4] = 0.5f*(fMaxMinVals[3] +fMaxMinVals[5]);
                        
                        uElem = 0u; 
                        memset(u_sElem, 0, u_ElemCntr[i]* sizeof(unsigned int) ); 
                        memset(u_newBr, 0, 4*2* sizeof(unsigned int) ); 
                        
                        for (k = 0u; k < u_ElemCntr[i]; k++)
                        {   uTemp = (u_OTElem_temp[u_MaxBrLvl*k +uBrLevel-1] == u_BrIds[j])*1u + (u_OTElem_temp[u_MaxBrLvl*k +uBrLevel-1] != u_BrIds[j])*0u;
                            u_sElem[uElem] = (uTemp == 1u)*k          + (uTemp != 1u)*u_sElem[uElem];
                            uElem          = (uTemp == 1u)*(uElem+1u) + (uTemp != 1u)*uElem;
                        }
                        
                        for (k = 0u; k < uElem; k++)
                        {   if ( (fCentRot[u_sElem[k]*3 +1] >= fMaxMinVals[1])       &&  (fCentRot[u_sElem[k]*3 +2] >= fMaxMinVals[4]) )
                            {   if (u_newBr[0*4 +0] == 0u)
                                {   uBrCount += 1u;      u_newBr[0*4 +0] = 1u;            u_newBr[1*4 +0] = uBrCount;
                                    f_nBrLims[unBrNum*4 +0] = fMaxMinVals[1];             f_nBrLims[unBrNum*4 +1] = fMaxMinVals[2];
                                    f_nBrLims[unBrNum*4 +2] = fMaxMinVals[4];             f_nBrLims[unBrNum*4 +3] = fMaxMinVals[5];
                                    u_nBrIds[unBrNum] = uBrCount;                         unBrNum += 1u;
                                    uF_OTPrtChld[(u_BrIds[j]+uF_TotalOffs)*5 +0] +=1u;    uF_OTPrtChld[(u_BrIds[j]+uF_TotalOffs)*5 +uF_OTPrtChld[(u_BrIds[j]+uF_TotalOffs)*5 +0]] = (uBrCount +uF_TotalOffs);
                                }
                                u_OTElem_temp[u_MaxBrLvl*u_sElem[k] +uBrLevel] = u_newBr[1*4 +0];
                            }
                            else if ( (fCentRot[u_sElem[k]*3 +1] < fMaxMinVals[1])  &&  (fCentRot[u_sElem[k]*3 +2] >= fMaxMinVals[4]) )
                            {   if (u_newBr[0*4 +1] == 0u)
                                {   uBrCount += 1u;      u_newBr[0*4 +1] = 1u;            u_newBr[1*4 +1] = uBrCount;
                                    f_nBrLims[unBrNum*4 +0] = fMaxMinVals[0];             f_nBrLims[unBrNum*4 +1] = fMaxMinVals[1];
                                    f_nBrLims[unBrNum*4 +2] = fMaxMinVals[4];             f_nBrLims[unBrNum*4 +3] = fMaxMinVals[5];
                                    u_nBrIds[unBrNum] = uBrCount;                         unBrNum += 1u;
                                    uF_OTPrtChld[(u_BrIds[j]+uF_TotalOffs)*5 +0] +=1u;    uF_OTPrtChld[(u_BrIds[j]+uF_TotalOffs)*5 +uF_OTPrtChld[(u_BrIds[j]+uF_TotalOffs)*5 +0]] = (uBrCount +uF_TotalOffs);
                                }
                                u_OTElem_temp[u_MaxBrLvl*u_sElem[k] +uBrLevel] = u_newBr[1*4 +1];
                            }
                            else if ( (fCentRot[u_sElem[k]*3 +1] < fMaxMinVals[1])  &&  (fCentRot[u_sElem[k]*3 +2] < fMaxMinVals[4]) )
                            {   if (u_newBr[0*4 +2] == 0u)
                                {   uBrCount += 1u;      u_newBr[0*4 +2] = 1u;            u_newBr[1*4 +2] = uBrCount;
                                    f_nBrLims[unBrNum*4 +0] = fMaxMinVals[0];             f_nBrLims[unBrNum*4 +1] = fMaxMinVals[1];
                                    f_nBrLims[unBrNum*4 +2] = fMaxMinVals[3];             f_nBrLims[unBrNum*4 +3] = fMaxMinVals[4];
                                    u_nBrIds[unBrNum] = uBrCount;                         unBrNum += 1u;
                                    uF_OTPrtChld[(u_BrIds[j]+uF_TotalOffs)*5 +0] +=1u;    uF_OTPrtChld[(u_BrIds[j]+uF_TotalOffs)*5 +uF_OTPrtChld[(u_BrIds[j]+uF_TotalOffs)*5 +0]] = (uBrCount +uF_TotalOffs);
                                }
                                u_OTElem_temp[u_MaxBrLvl*u_sElem[k] +uBrLevel] = u_newBr[1*4 +2];
                            }
                            else if ( (fCentRot[u_sElem[k]*3 +1] >=  fMaxMinVals[1])  &&  (fCentRot[u_sElem[k]*3 +2] < fMaxMinVals[4]) )
                            {   if (u_newBr[0*4 +3] == 0u)
                                {   uBrCount += 1u;      u_newBr[0*4 +3] = 1u;            u_newBr[1*4 +3] = uBrCount;
                                    f_nBrLims[unBrNum*4 +0] = fMaxMinVals[1];             f_nBrLims[unBrNum*4 +1] = fMaxMinVals[2];
                                    f_nBrLims[unBrNum*4 +2] = fMaxMinVals[3];             f_nBrLims[unBrNum*4 +3] = fMaxMinVals[4];
                                    u_nBrIds[unBrNum] = uBrCount;                         unBrNum += 1u;
                                    uF_OTPrtChld[(u_BrIds[j]+uF_TotalOffs)*5 +0] +=1u;    uF_OTPrtChld[(u_BrIds[j]+uF_TotalOffs)*5 +uF_OTPrtChld[(u_BrIds[j]+uF_TotalOffs)*5 +0]] = (uBrCount +uF_TotalOffs);
                                }
                                u_OTElem_temp[u_MaxBrLvl*u_sElem[k] +uBrLevel] = u_newBr[1*4 +3];
                    }   }   }
                    memcpy(u_BrIds,  u_nBrIds,    uMaxPNum* sizeof(unsigned int) );
                    memcpy(f_BrLims, f_nBrLims, 4*uMaxPNum* sizeof(float) );
                    uBrNum = unBrNum;
                }
                
                uBrCount   += 1u;
                uBrLevel   += 2u;
                for (j = 0; j < u_ElemCntr[i]; j++ )
                {   
                    u_OTElem_temp[u_MaxBrLvl*j + uBrLevel -1] = (uBrCount+ j);
                    
                    uTemp = u_ElSegIDs[i*uMaxPNum +j]*u_MaxBrLvl;
                    uF_OTElem[uTemp +0] = uBrLevel; 
                    
                    for (k = 0; k < uBrLevel; k++)
                    {   uF_OTElem[uTemp +(u_MaxBrLvl -uBrLevel +k)] = u_OTElem_temp[u_MaxBrLvl*j +k] +uF_TotalOffs;
                }   }
                
                uF_TotalOffs += (uBrCount + u_ElemCntr[i] ); 
                
            }
            
            if (uBPNum > 0u)
            {   memset(f_SegVals,  0,        6*uMaxSNum* sizeof(float));
                memset(u_ElemCntr, 0,        1*uMaxSNum* sizeof(unsigned int));
                memset(u_ElSegIDs, 0, uMaxPNum*uMaxSNum* sizeof(unsigned int));
                
                for (i = 0u; i < uBPNum; i++)
                {   uNxtSegID = uB_temp[i*4 +3];
                    u_ElSegIDs[uNxtSegID*uMaxPNum + u_ElemCntr[uNxtSegID]] = i;
                    u_ElemCntr[uNxtSegID] += 1u;
                    
                    f_SegVals[uNxtSegID*6 +0] += fB_temp[i*16 +0];
                    f_SegVals[uNxtSegID*6 +1] += fB_temp[i*16 +1];
                    f_SegVals[uNxtSegID*6 +2] += fB_temp[i*16 +2];
                    
                    fTemp = (fB_temp[i*16 +8] >= 0.0f)*1.0f + (fB_temp[i*16 +8] < 0.0f)*-1.0f;
                    f_SegVals[uNxtSegID*6 +3] += (fB_temp[i*16 +6]*fTemp);
                    f_SegVals[uNxtSegID*6 +4] += (fB_temp[i*16 +7]*fTemp);
                    f_SegVals[uNxtSegID*6 +5] += (fB_temp[i*16 +8]*fTemp);
                }
                
                for (i = 0u; i < uBSegN; i++)
                {   
                    u_OTElem_temp = realloc(u_OTElem_temp,    u_ElemCntr[i]*u_MaxBrLvl* sizeof *u_OTElem_temp);
                    memset(u_OTElem_temp,    0, u_ElemCntr[i]*u_MaxBrLvl* sizeof(unsigned int) );
                    memset(u_BrIds,          0, u_ElemCntr[i]*            sizeof(unsigned int) );
                    
                    f_SegVals[i*6 +0] /= (float)u_ElemCntr[i];        f_SegVals[i*6 +1] /= (float)u_ElemCntr[i];        f_SegVals[i*6 +2] /= (float)u_ElemCntr[i];
                    f_SegVals[i*6 +3] /= (float)u_ElemCntr[i];        f_SegVals[i*6 +4] /= (float)u_ElemCntr[i];        f_SegVals[i*6 +5] /= (float)u_ElemCntr[i];
                    memcpy(fNrm, &f_SegVals[i*6 +3], 3*sizeof(float));
                    GetStkAndDipVect(fNrm, fStk, fDip);
                    
                    fTemp = -FLT_MAX;
                    
                    for (j = 0; j < u_ElemCntr[i]; j++)
                    {   fShftCent[0] = fB_temp[ u_ElSegIDs[i*uMaxPNum +j]*16 +0] - f_SegVals[i*6 +0];
                        fShftCent[1] = fB_temp[ u_ElSegIDs[i*uMaxPNum +j]*16 +1] - f_SegVals[i*6 +1];
                        fShftCent[2] = fB_temp[ u_ElSegIDs[i*uMaxPNum +j]*16 +2] - f_SegVals[i*6 +2];
                        
                        fCentRot[j*3 +0] = fNrm[0]*fShftCent[0] + fNrm[1]*fShftCent[1] + fNrm[2]*fShftCent[2];
                        fCentRot[j*3 +1] = fStk[0]*fShftCent[0] + fStk[1]*fShftCent[1] + fStk[2]*fShftCent[2];
                        fCentRot[j*3 +2] = fDip[0]*fShftCent[0] + fDip[1]*fShftCent[1] + fDip[2]*fShftCent[2];
                        
                        fTemp = MAX(fTemp, fCentRot[j*3 +2]);
                    }
                    for (j = 0; j < u_ElemCntr[i]; j++)     {   fCentRot[j*3 +2] -= fTemp;      }
                    fMaxMinVals[0] = fCentRot[0*3 +1];                                fMaxMinVals[2] = fCentRot[0*3 +1];
                    fMaxMinVals[3] = fCentRot[0*3 +2];                                fMaxMinVals[5] = fCentRot[0*3 +2];
                    
                    for (j = 1; j < u_ElemCntr[i]; j++)
                    {   fMaxMinVals[0] = MIN(fMaxMinVals[0], fCentRot[j*3 +1]);       fMaxMinVals[2] = MAX(fMaxMinVals[2], fCentRot[j*3 +1]);
                        fMaxMinVals[3] = MIN(fMaxMinVals[3], fCentRot[j*3 +2]);       fMaxMinVals[5] = MAX(fMaxMinVals[5], fCentRot[j*3 +2]);
                    }
                    
                    fMaxMinVals[0]  -= fBndSide;                    fMaxMinVals[2] += fBndSide;
                    fMaxMinVals[3]  -= fBndSide;                    fMaxMinVals[5] += fBndSide;
                    fBrSize          = MAX((fMaxMinVals[2] -fMaxMinVals[0]), (fMaxMinVals[5] -fMaxMinVals[3]));
                    f_BrLims[0*4 +0] = -0.5f*fBrSize;              f_BrLims[0*4 +1] = 0.5f*fBrSize;
                    f_BrLims[0*4 +2] = -1.0f*fBrSize;              f_BrLims[0*4 +3] = 0.0f*fBrSize;
                    
                    uBrLevel = 0u; 
                    uBrCount = 0u;
                    uBrNum   = 1u;
                    
                    while (0.5f*fBrSize >= MINBRANCHSZEFACT*fBndSide) 
                    {   uBrLevel += 1u; 
                        fBrSize  *= 0.5f; 
                        unBrNum   = 0u; 
                        
                        for (j = 0u; j < uBrNum; j++)
                        {   fMaxMinVals[0] = f_BrLims[j*4 +0];      fMaxMinVals[2] = f_BrLims[j*4 +1];       fMaxMinVals[1] = 0.5f*(fMaxMinVals[0] +fMaxMinVals[2]);
                            fMaxMinVals[3] = f_BrLims[j*4 +2];      fMaxMinVals[5] = f_BrLims[j*4 +3];       fMaxMinVals[4] = 0.5f*(fMaxMinVals[3] +fMaxMinVals[5]);
                            
                            uElem = 0u; 
                            memset(u_sElem, 0, u_ElemCntr[i]* sizeof(unsigned int) ); 
                            memset(u_newBr, 0, 4*2* sizeof(unsigned int) ); 
                            
                            for (k = 0u; k < u_ElemCntr[i]; k++)
                            {   uTemp = (u_OTElem_temp[u_MaxBrLvl*k +uBrLevel-1] == u_BrIds[j])*1u + (u_OTElem_temp[u_MaxBrLvl*k +uBrLevel-1] != u_BrIds[j])*0u;
                                u_sElem[uElem] = (uTemp == 1u)*k          + (uTemp != 1u)*u_sElem[uElem];
                                uElem          = (uTemp == 1u)*(uElem+1u) + (uTemp != 1u)*uElem;
                            }
                            
                            for (k = 0u; k < uElem; k++)
                            {   if ( (fCentRot[u_sElem[k]*3 +1] >= fMaxMinVals[1])       &&  (fCentRot[u_sElem[k]*3 +2] >= fMaxMinVals[4]) )
                                {   if (u_newBr[0*4 +0] == 0u)
                                    {   uBrCount += 1u;      u_newBr[0*4 +0] = 1u;            u_newBr[1*4 +0] = uBrCount;
                                        f_nBrLims[unBrNum*4 +0] = fMaxMinVals[1];             f_nBrLims[unBrNum*4 +1] = fMaxMinVals[2];
                                        f_nBrLims[unBrNum*4 +2] = fMaxMinVals[4];             f_nBrLims[unBrNum*4 +3] = fMaxMinVals[5];
                                        u_nBrIds[unBrNum] = uBrCount;                         unBrNum += 1u;
                                        uB_OTPrtChld[(u_BrIds[j]+uB_TotalOffs)*5 +0] +=1u;    uB_OTPrtChld[(u_BrIds[j]+uB_TotalOffs)*5 +uB_OTPrtChld[(u_BrIds[j]+uB_TotalOffs)*5 +0]] = (uBrCount +uB_TotalOffs);
                                    }
                                    u_OTElem_temp[u_MaxBrLvl*u_sElem[k] +uBrLevel] = u_newBr[1*4 +0];
                                }
                                else if ( (fCentRot[u_sElem[k]*3 +1] < fMaxMinVals[1])  &&  (fCentRot[u_sElem[k]*3 +2] >= fMaxMinVals[4]) )
                                {   if (u_newBr[0*4 +1] == 0u)
                                    {   uBrCount += 1u;      u_newBr[0*4 +1] = 1u;            u_newBr[1*4 +1] = uBrCount;
                                        f_nBrLims[unBrNum*4 +0] = fMaxMinVals[0];             f_nBrLims[unBrNum*4 +1] = fMaxMinVals[1];
                                        f_nBrLims[unBrNum*4 +2] = fMaxMinVals[4];             f_nBrLims[unBrNum*4 +3] = fMaxMinVals[5];
                                        u_nBrIds[unBrNum] = uBrCount;                         unBrNum += 1u;
                                        uB_OTPrtChld[(u_BrIds[j]+uB_TotalOffs)*5 +0] +=1u;    uB_OTPrtChld[(u_BrIds[j]+uB_TotalOffs)*5 +uB_OTPrtChld[(u_BrIds[j]+uB_TotalOffs)*5 +0]] = (uBrCount +uB_TotalOffs);
                                    }
                                    u_OTElem_temp[u_MaxBrLvl*u_sElem[k] +uBrLevel] = u_newBr[1*4 +1];
                                }
                                else if ( (fCentRot[u_sElem[k]*3 +1] < fMaxMinVals[1])  &&  (fCentRot[u_sElem[k]*3 +2] < fMaxMinVals[4]) )
                                {   if (u_newBr[0*4 +2] == 0u)
                                    {   uBrCount += 1u;      u_newBr[0*4 +2] = 1u;            u_newBr[1*4 +2] = uBrCount;
                                        f_nBrLims[unBrNum*4 +0] = fMaxMinVals[0];             f_nBrLims[unBrNum*4 +1] = fMaxMinVals[1];
                                        f_nBrLims[unBrNum*4 +2] = fMaxMinVals[3];             f_nBrLims[unBrNum*4 +3] = fMaxMinVals[4];
                                        u_nBrIds[unBrNum] = uBrCount;                         unBrNum += 1u;
                                        uB_OTPrtChld[(u_BrIds[j]+uB_TotalOffs)*5 +0] +=1u;    uB_OTPrtChld[(u_BrIds[j]+uB_TotalOffs)*5 +uB_OTPrtChld[(u_BrIds[j]+uB_TotalOffs)*5 +0]] = (uBrCount +uB_TotalOffs);
                                    }
                                    u_OTElem_temp[u_MaxBrLvl*u_sElem[k] +uBrLevel] = u_newBr[1*4 +2];
                                }
                                else if ( (fCentRot[u_sElem[k]*3 +1] >=  fMaxMinVals[1])  &&  (fCentRot[u_sElem[k]*3 +2] < fMaxMinVals[4]) )
                                {   if (u_newBr[0*4 +3] == 0u)
                                    {   uBrCount += 1u;      u_newBr[0*4 +3] = 1u;            u_newBr[1*4 +3] = uBrCount;
                                        f_nBrLims[unBrNum*4 +0] = fMaxMinVals[1];             f_nBrLims[unBrNum*4 +1] = fMaxMinVals[2];
                                        f_nBrLims[unBrNum*4 +2] = fMaxMinVals[3];             f_nBrLims[unBrNum*4 +3] = fMaxMinVals[4];
                                        u_nBrIds[unBrNum] = uBrCount;                         unBrNum += 1u;
                                        uB_OTPrtChld[(u_BrIds[j]+uB_TotalOffs)*5 +0] +=1u;    uB_OTPrtChld[(u_BrIds[j]+uB_TotalOffs)*5 +uB_OTPrtChld[(u_BrIds[j]+uB_TotalOffs)*5 +0]] = (uBrCount +uB_TotalOffs);
                                    }
                                    u_OTElem_temp[u_MaxBrLvl*u_sElem[k] +uBrLevel] = u_newBr[1*4 +3];
                        }   }   }
                        memcpy(u_BrIds, u_nBrIds,     uMaxPNum* sizeof(unsigned int) );
                        memcpy(f_BrLims, f_nBrLims, 4*uMaxPNum* sizeof(float) );
                        uBrNum = unBrNum;
                    }
                    
                    uBrCount   += 1u;
                    uBrLevel   += 2u;
                    for (j = 0; j < u_ElemCntr[i]; j++ )
                    {   
                        u_OTElem_temp[u_MaxBrLvl*j + uBrLevel -1] = (uBrCount+ j);
                        
                        uTemp = u_ElSegIDs[i*uMaxPNum +j]*u_MaxBrLvl;
                        uB_OTElem[uTemp +0] = uBrLevel; 
                        
                        for (k = 0; k < uBrLevel; k++)
                        {   uB_OTElem[uTemp +(u_MaxBrLvl -uBrLevel +k)] = u_OTElem_temp[u_MaxBrLvl*j +k] +uB_TotalOffs;
                    }   }
                    
                    uB_TotalOffs += (uBrCount + u_ElemCntr[i] ); 
                    
            }   }
            
            timer0 = clock() - timer0;        time_taken  = ((double)timer0)/CLOCKS_PER_SEC;
            MPI_Allreduce(MPI_IN_PLACE, &time_taken, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); 
            time_taken /=(double)iSIZE;
            if (iRANK == 0) 
            {   fprintf(stdout,"Total RunTime for OctTree (seconds): %6.4f   (min. branch size factor: %2.2f)\nBranchCountF: %7u\nBranchCountB:%7u\n",time_taken, MINBRANCHSZEFACT, uF_TotalOffs, uB_TotalOffs);
            }
        } 
        
        
        {   timer0 = clock();
            unsigned int uGlobPos,   uUsdNum,   uSrcElem,  uTemp0,   uTemp1,   uTemp2,   uPrvBrLv,   uTstBrId,   uPrvBrId;
            float fTemp0,   fTemp1,   fTemp2,   fTemp3,   fTemp4,   fTemp5,   fDist,  fDistF;
            float fRcv[3],   fSrc[3],   fP1[3],   fP2[3],   fP3[3],   fP4[3],   fStressStk[6],   fStressDip[6],   fStressNrm[6],   fStressStkT[6],   fStressDipT[6],   fStressNrmT[6],   fStrain[6];
            float fRotMat[9], fTempVect[3],   fPrvStress[9],   fTstStress[9];
            float fNrm[3],   fStk[3],   fDip[3];

            unsigned int *uFUsdEl = calloc(  uFPNum, sizeof *uFUsdEl );
            unsigned int *uFPrvEl = calloc(  uFPNum, sizeof *uFPrvEl );
            unsigned int *uFTstEl = calloc(  uFPNum, sizeof *uFTstEl );
            float *fKh_FFBrch     = calloc(21*uF_TotalOffs, sizeof *fKh_FFBrch );
            unsigned int *uFUsdBrch = calloc( uF_TotalOffs, sizeof *uFUsdBrch );
            unsigned int *uBUsdEl = calloc(  uBPNum, sizeof *uBUsdEl );
            unsigned int *uBPrvEl = calloc(  uBPNum, sizeof *uBPrvEl );
            unsigned int *uBTstEl = calloc(  uBPNum, sizeof *uBTstEl );
            float *fKh_BBBrch     = calloc(21*uB_TotalOffs, sizeof *fKh_BBBrch );
            unsigned int *uBUsdBrch = calloc( uB_TotalOffs, sizeof *uBUsdBrch );
            
            
            for (j = 0u; j < uFPNum; j++)
            {   uTemp0 = uF_OTElem[u_MaxBrLvl*j +0];
                for (k = 0u; k < uTemp0; k++)
                {   uTemp1 = uF_OTElem[j*u_MaxBrLvl +(u_MaxBrLvl -uTemp0 +k)]*21u;
                    fKh_FFBrch[uTemp1 +15] = FLT_MAX;           fKh_FFBrch[uTemp1 +16] = -FLT_MAX;
                    fKh_FFBrch[uTemp1 +17] = FLT_MAX;           fKh_FFBrch[uTemp1 +18] = -FLT_MAX;
                    fKh_FFBrch[uTemp1 +19] = FLT_MAX;           fKh_FFBrch[uTemp1 +20] = -FLT_MAX;
            }   }
            
            for (j = 0u; j < uFPNum; j++)
            {   uTemp0 = uF_OTElem[u_MaxBrLvl*j +0];
                memcpy(fP1, &fFV_temp[uF_temp[j*7 +0]*4 +0],  3u*sizeof(float));
                memcpy(fP2, &fFV_temp[uF_temp[j*7 +1]*4 +0],  3u*sizeof(float));
                memcpy(fP3, &fFV_temp[uF_temp[j*7 +2]*4 +0],  3u*sizeof(float));
                fTemp0 = MIN(fP1[0], fP2[0]);       fTemp0 = MIN(fTemp0, fP3[0]);            fTemp1 = MAX(fP1[0], fP2[0]);       fTemp1 = MAX(fTemp1, fP3[0]);
                fTemp2 = MIN(fP1[1], fP2[1]);       fTemp2 = MIN(fTemp2, fP3[1]);            fTemp3 = MAX(fP1[1], fP2[1]);       fTemp3 = MAX(fTemp3, fP3[1]);
                fTemp4 = MIN(fP1[2], fP2[2]);       fTemp4 = MIN(fTemp4, fP3[2]);            fTemp5 = MAX(fP1[2], fP2[2]);       fTemp5 = MAX(fTemp5, fP3[2]);
                for (k = 0u; k < uTemp0; k++)
                {   uTemp1 = uF_OTElem[j*u_MaxBrLvl +(u_MaxBrLvl -uTemp0 +k)]*21u;
                    
                    
                    
                    fKh_FFBrch[uTemp1 +0] += fF_temp[j*16 +0];       fKh_FFBrch[uTemp1 +1] += fF_temp[j*16 +1];    fKh_FFBrch[uTemp1 +2] += fF_temp[j*16 +2];
                    fKh_FFBrch[uTemp1 +3] += fF_temp[j*16 +15];
                    fKh_FFBrch[uTemp1 +4] += fF_temp[j*16 +6];       fKh_FFBrch[uTemp1 +5] += fF_temp[j*16 +7];    fKh_FFBrch[uTemp1 +6] += fF_temp[j*16 +8];
                    fKh_FFBrch[uTemp1 +13]+= 1.0f;
                    fKh_FFBrch[uTemp1 +15] = MIN(fKh_FFBrch[uTemp1 +15], fTemp0);           fKh_FFBrch[uTemp1 +16] = MAX(fKh_FFBrch[uTemp1 +16], fTemp1);
                    fKh_FFBrch[uTemp1 +17] = MIN(fKh_FFBrch[uTemp1 +17], fTemp2);           fKh_FFBrch[uTemp1 +18] = MAX(fKh_FFBrch[uTemp1 +18], fTemp3);
                    fKh_FFBrch[uTemp1 +19] = MIN(fKh_FFBrch[uTemp1 +19], fTemp4);           fKh_FFBrch[uTemp1 +20] = MAX(fKh_FFBrch[uTemp1 +20], fTemp5);
            }   }
            
            for (j = 0u; j < uF_TotalOffs; j++)
            {   fKh_FFBrch[j*21 +0] /= fKh_FFBrch[j*21 +13];    fKh_FFBrch[j*21 +1] /= fKh_FFBrch[j*21 +13];        fKh_FFBrch[j*21 +2] /= fKh_FFBrch[j*21 +13];
                memcpy(fNrm, &fKh_FFBrch[j*21 +4], 3u*sizeof(float));
                NrmlzeVect(fNrm);
                memcpy(&fKh_FFBrch[j*21 +4], fNrm, 3u*sizeof(float));
                
                GetStkAndDipVect(fNrm, fStk, fDip);
                
                memcpy(&fKh_FFBrch[j*21 +7], fStk, 3u*sizeof(float));
                memcpy(&fKh_FFBrch[j*21+10], fDip, 3u*sizeof(float));
                
                fTemp0 = MAX(fabs(fKh_FFBrch[j*21 +15] - fKh_FFBrch[j*21 +0]),  fabs(fKh_FFBrch[j*21 +16] - fKh_FFBrch[j*21 +0]));
                fTemp1 = MAX(fabs(fKh_FFBrch[j*21 +17] - fKh_FFBrch[j*21 +1]),  fabs(fKh_FFBrch[j*21 +18] - fKh_FFBrch[j*21 +1]));
                fTemp2 = MAX(fabs(fKh_FFBrch[j*21 +19] - fKh_FFBrch[j*21 +2]),  fabs(fKh_FFBrch[j*21 +20] - fKh_FFBrch[j*21 +2]));
                
                fKh_FFBrch[j*21 +14] = sqrtf(fTemp0*fTemp0 +fTemp1*fTemp1 +fTemp2*fTemp2);
                
            }
            
            for (j = 0u; j < uBPNum; j++)
            {   uTemp0 = uB_OTElem[u_MaxBrLvl*j +0];
                for (k = 0u; k < uTemp0; k++)
                {   uTemp1 = uB_OTElem[u_MaxBrLvl*j +(u_MaxBrLvl -uTemp0 +k)]*21u;
                    fKh_BBBrch[uTemp1 +15] = FLT_MAX;           fKh_BBBrch[uTemp1 +16] = -FLT_MAX;
                    fKh_BBBrch[uTemp1 +17] = FLT_MAX;           fKh_BBBrch[uTemp1 +18] = -FLT_MAX;
                    fKh_BBBrch[uTemp1 +19] = FLT_MAX;           fKh_BBBrch[uTemp1 +20] = -FLT_MAX;
            }   }
            
            for (j = 0u; j < uBPNum; j++)
            {   uTemp0 = uB_OTElem[u_MaxBrLvl*j +0];
                memcpy(fP1, &fBV_temp[uB_temp[j*4 +0]*3 +0],  3u*sizeof(float));
                memcpy(fP2, &fBV_temp[uB_temp[j*4 +1]*3 +0],  3u*sizeof(float));
                memcpy(fP3, &fBV_temp[uB_temp[j*4 +2]*3 +0],  3u*sizeof(float));
                fTemp0 = MIN(fP1[0], fP2[0]);       fTemp0 = MIN(fTemp0, fP3[0]);            fTemp1 = MAX(fP1[0], fP2[0]);       fTemp1 = MAX(fTemp1, fP3[0]);
                fTemp2 = MIN(fP1[1], fP2[1]);       fTemp2 = MIN(fTemp2, fP3[1]);            fTemp3 = MAX(fP1[1], fP2[1]);       fTemp3 = MAX(fTemp3, fP3[1]);
                fTemp4 = MIN(fP1[2], fP2[2]);       fTemp4 = MIN(fTemp4, fP3[2]);            fTemp5 = MAX(fP1[2], fP2[2]);       fTemp5 = MAX(fTemp5, fP3[2]);
                for (k = 0u; k < uTemp0; k++)
                {   uTemp1 = uB_OTElem[u_MaxBrLvl*j +(u_MaxBrLvl -uTemp0 +k)]*21u;
                    
                    
                    
                    fKh_BBBrch[uTemp1 +0] += fB_temp[j*16 +0];       fKh_BBBrch[uTemp1 +1] += fB_temp[j*16 +1];    fKh_BBBrch[uTemp1 +2] += fB_temp[j*16 +2];
                    fKh_BBBrch[uTemp1 +3] += fB_temp[j*16 +15];
                    fKh_BBBrch[uTemp1 +4] += fB_temp[j*16 +6];       fKh_BBBrch[uTemp1 +5] += fB_temp[j*16 +7];    fKh_BBBrch[uTemp1 +6] += fB_temp[j*16 +8];
                    fKh_BBBrch[uTemp1 +13]+= 1.0f;
                    fKh_BBBrch[uTemp1 +15] = MIN(fKh_BBBrch[uTemp1 +15], fTemp0);           fKh_BBBrch[uTemp1 +16] = MAX(fKh_BBBrch[uTemp1 +16], fTemp1);
                    fKh_BBBrch[uTemp1 +17] = MIN(fKh_BBBrch[uTemp1 +17], fTemp2);           fKh_BBBrch[uTemp1 +18] = MAX(fKh_BBBrch[uTemp1 +18], fTemp3);
                    fKh_BBBrch[uTemp1 +19] = MIN(fKh_BBBrch[uTemp1 +19], fTemp4);           fKh_BBBrch[uTemp1 +20] = MAX(fKh_BBBrch[uTemp1 +20], fTemp5);
            }   }
             
            for (j = 0u; j < uB_TotalOffs; j++)
            {   fKh_BBBrch[j*21 +0] /= fKh_BBBrch[j*21 +13];    fKh_BBBrch[j*21 +1] /= fKh_BBBrch[j*21 +13];        fKh_BBBrch[j*21 +2] /= fKh_BBBrch[j*21 +13];
                memcpy(fNrm, &fKh_BBBrch[j*21 +4], 3u*sizeof(float));
                NrmlzeVect(fNrm);
                memcpy(&fKh_BBBrch[j*21 +4], fNrm, 3u*sizeof(float));
                
                GetStkAndDipVect(fNrm, fStk, fDip);
                
                memcpy(&fKh_BBBrch[j*21 +7], fStk, 3u*sizeof(float));
                memcpy(&fKh_BBBrch[j*21+10], fDip, 3u*sizeof(float));
                
                fTemp0 = MAX(fabs(fKh_BBBrch[j*21 +15] - fKh_BBBrch[j*21 +0]),  fabs(fKh_BBBrch[j*21 +16] - fKh_BBBrch[j*21 +0]));
                fTemp1 = MAX(fabs(fKh_BBBrch[j*21 +17] - fKh_BBBrch[j*21 +1]),  fabs(fKh_BBBrch[j*21 +18] - fKh_BBBrch[j*21 +1]));
                fTemp2 = MAX(fabs(fKh_BBBrch[j*21 +19] - fKh_BBBrch[j*21 +2]),  fabs(fKh_BBBrch[j*21 +20] - fKh_BBBrch[j*21 +2]));

                fKh_BBBrch[j*21 +14] = sqrtf(fTemp0*fTemp0 +fTemp1*fTemp1 +fTemp2*fTemp2);
                
            }
            
            
            for (i = 0u; i < iFOFFSET[iRANK]; i++) 
            {   
                uKh_FFps2[i]    = calloc(     uFPNum, sizeof *uKh_FFps2[i] );
                fKh_FFvalstk[i] = calloc(   2*uFPNum, sizeof *fKh_FFvalstk[i] );
                fKh_FFvaldip[i] = calloc(   2*uFPNum, sizeof *fKh_FFvaldip[i] );
                fKh_FFvalnrm[i] = calloc(   2*uFPNum, sizeof *fKh_FFvalnrm[i] );
                
                uGlobPos = i+ iFSTART[iRANK];
                memcpy(fRcv, &fF_temp[uGlobPos*16 +0], 3u*sizeof(float));
                memcpy(fRotMat, &fF_temp[uGlobPos*16 +6], 9u*sizeof(float)); 
                
                uKh_FFcnt[i]  = 0u;
                uUsdNum       = 0u;
                memset(uFUsdEl, 0,         uFPNum*sizeof(unsigned int));
                memset(uFUsdBrch, 0, uF_TotalOffs*sizeof(unsigned int));
                uTstBrId = 0u;      uSrcElem = 0u;
                
                while (uUsdNum < uFPNum) 
                {   
                    for (j = 0u; j < uFPNum; j++)               {       if (uFUsdEl[j] == 0u)    {      uSrcElem = j;       break;          }           } 
                    

                    memcpy(fSrc, &fF_temp[uSrcElem*16 +0], 3u*sizeof(float));
                    
                    SubtractVect3(fTempVect, fRcv, fSrc);               fDist  = VectorLength(fTempVect);
                    
                    fDistF = (fDist >= CUTDISTANCE*fFltSide)*1.0f                    + (fDist < CUTDISTANCE*fFltSide)*0.0f; 
                    fDistF = (uF_temp[uGlobPos*7 +4] == uF_temp[uSrcElem*7 +4])*1.0f + (uF_temp[uGlobPos*7 +4] != uF_temp[uSrcElem*7 +4])*fDistF;
                    uF_temp[uGlobPos*7 +6] = (fDistF == 1.0f)*uF_temp[uGlobPos*7 +6] + (fDistF != 1.0f)*1u;
                    
                    uTemp0 = uF_OTElem[uSrcElem*u_MaxBrLvl +0];
                    uTemp1 = 0u;
                    
                    for (j = 1u; j < uTemp0; j++)
                    {   uPrvBrLv = j;
                        uPrvBrId = uF_OTElem[uSrcElem*u_MaxBrLvl +(u_MaxBrLvl -uTemp0 +uPrvBrLv)];
                        
                        memcpy(fSrc, &fKh_FFBrch[uPrvBrId*21 +0], 3u*sizeof(float));
                        SubtractVect3(fTempVect, fRcv, fSrc);               fDist  = VectorLength(fTempVect);
                        
                        uTemp1  = (uFUsdBrch[uPrvBrId]                                   == 0u)*1u     + (uFUsdBrch[uPrvBrId]                                      != 0u)*0u;
                        uTemp1  = (uF_OTPrtChld[uPrvBrId*5 +0]                            > 1u)*uTemp1 + (uF_OTPrtChld[uPrvBrId*5 +0]                              <= 1u)*0u;
                        uTemp1  = (fDist > (fKh_FFBrch[uPrvBrId*21 +14]+KH_FULL_DIST*fFltSide))*uTemp1 + (fDist <= (fKh_FFBrch[uPrvBrId*21 +14]+KH_FULL_DIST*fFltSide))*0u;
                        
                        if (uTemp1 == 1u)         {       break;      }
                    } 
                    
                    if (uTemp1 == 0u) 
                    {   memcpy(fP1, &fFV_temp[uF_temp[uSrcElem*7 +0]*4 +0], 3u*sizeof(float));
                        memcpy(fP2, &fFV_temp[uF_temp[uSrcElem*7 +1]*4 +0], 3u*sizeof(float));
                        memcpy(fP3, &fFV_temp[uF_temp[uSrcElem*7 +2]*4 +0], 3u*sizeof(float));
                        
                        if (USEHALFSPACE == 1u)
                        {   StrainHS_Nikkhoo(fStressStk, fStrain, fRcv[0], fRcv[1], fRcv[2], fP1, fP2, fP3, fUnitSlipF, 0.0f, 0.0f, fModPara[2], fModPara[4]);
                            StrainHS_Nikkhoo(fStressDip, fStrain, fRcv[0], fRcv[1], fRcv[2], fP1, fP2, fP3, 0.0f, fUnitSlipF, 0.0f, fModPara[2], fModPara[4]);
                            StrainHS_Nikkhoo(fStressNrm, fStrain, fRcv[0], fRcv[1], fRcv[2], fP1, fP2, fP3, 0.0f, 0.0f, fUnitSlipF, fModPara[2], fModPara[4]);
                        }
                        else
                        {   StrainFS_Nikkhoo(fStressStk, fStrain, fRcv[0], fRcv[1], fRcv[2], fP1, fP2, fP3, fUnitSlipF, 0.0f, 0.0f, fModPara[2], fModPara[4]);
                            StrainFS_Nikkhoo(fStressDip, fStrain, fRcv[0], fRcv[1], fRcv[2], fP1, fP2, fP3, 0.0f, fUnitSlipF, 0.0f, fModPara[2], fModPara[4]);
                            StrainFS_Nikkhoo(fStressNrm, fStrain, fRcv[0], fRcv[1], fRcv[2], fP1, fP2, fP3, 0.0f, 0.0f, fUnitSlipF, fModPara[2], fModPara[4]);
                        }
                        
                        ScaleStress6(fStressStk, (fDistF/fUnitSlipF));      RotStressG2L(&fPrvStress[0], fRotMat, fStressStk);
                        ScaleStress6(fStressDip, (fDistF/fUnitSlipF));      RotStressG2L(&fPrvStress[3], fRotMat, fStressDip);
                        ScaleStress6(fStressNrm, (fDistF/fUnitSlipF));      RotStressG2L(&fPrvStress[6], fRotMat, fStressNrm);
                        
                        if (uGlobPos == uSrcElem)
                        {   
                            fFEvent[i*17 +4] = fPrvStress[0];
                            fFEvent[i*17 +5] = fPrvStress[4];
                            fFEvent[i*17 +6] = fPrvStress[8];
                            fFRef[i*14 +2]   = MAX(fabs(fFEvent[i*17 +4]),fabs(fFEvent[i*17 +5])); 
                            fFRef[i*14 +2]   = MAX(fFRef[i*14 +2],        fabs(fFEvent[i*17 +6])); 
                            
                            fTemp0 = ran0(&lSeed);      fTemp0 = fTemp0*2.0f -1.0f;
                            fFEvent[i*17 +7] = fFRef[i*14 +7] + (fFRef[i*14 +8] * fTemp0); 
                            
                            fTemp0 = ran0(&lSeed);      fTemp0 = fTemp0*2.0f -1.0f;
                            fFFric[i*6 +0]   = fFRef[i*14 +3] + (fFRef[i*14 +4] * fTemp0); 
                            fFEvent[i*17 +3] = fFFric[i*6 +0]; 
                            
                            fTemp0 = ran0(&lSeed);      fTemp0 = fTemp0*2.0f -1.0f;
                            fFFric[i*6 +1]   = fFRef[i*14 +5] + (fFRef[i*14 +6] * fTemp0); 
                                
                            fTemp0 =(fFFric[i*6 +0] - fFFric[i*6 +1]) *-1.0f*fFRef[i*14 +1];; 
                            fTemp1 = fTemp0/fFRef[i*14 +2];
                            uFEvent[i*5 +0]  = (fTemp0 >= 0.0f)*2u + (fTemp0 < 0.0f)*3u;
                            uFEvent[i*5 +0]  = (fTemp1 > fFEvent[i*17 +7])*1u + (fTemp1 <= fFEvent[i*17 +7])*uFEvent[i*5 +0];
                            
                            fTemp0 = fFFric[i*6 +1]*-1.0f*fFRef[i*14 +1] + fFRef[i*14 +2]*fFEvent[i*17 +7]; 
                            fTemp1 = fFFric[i*6 +0]*-1.0f*fFRef[i*14 +1]; 
                            fTemp2 = MAX(fTemp0, fTemp1)/(-1.0f*fFRef[i*14 +1]); 
                            fFFric[i*6 +2]   = fTemp2 - fOvershootFract*(fTemp2-fFFric[i*6 +1]); 
                            
                            fFEvent[i*17 +10]= (fFFric[i*6 +1] - fFFric[i*6 +2])*-1.0f*fFRef[i*14 +1];
                            
                        }
                        
                        fKh_FFvalstk[i][uKh_FFcnt[i]*2 +0] = fPrvStress[0];         fKh_FFvalstk[i][uKh_FFcnt[i]*2 +1] = fPrvStress[3];
                        fKh_FFvaldip[i][uKh_FFcnt[i]*2 +0] = fPrvStress[1];         fKh_FFvaldip[i][uKh_FFcnt[i]*2 +1] = fPrvStress[4];
                        fKh_FFvalnrm[i][uKh_FFcnt[i]*2 +0] = fPrvStress[2];         fKh_FFvalnrm[i][uKh_FFcnt[i]*2 +1] = fPrvStress[5];
                        
                        uFUsdEl[uSrcElem]       = 1u;
                        uKh_FFps2[i][uSrcElem]  = uKh_FFcnt[i];
                        uKh_FFcnt[i]           += 1u;
                        
                        uUsdNum = 0u;
                        for (j = 0u; j < uFPNum;  j++)      {       uUsdNum +=  uFUsdEl[j];                                     }
                        
                        for (j = 1u; j < uTemp0; j++)       {       uFUsdBrch[ uF_OTElem[uSrcElem*u_MaxBrLvl +(u_MaxBrLvl -uTemp0 +j)] ] = 1u;       }
                    }
                    else
                    {   uPrvBrLv -= 1u;
                        uTemp2    = 1u;
                        while (uTemp2 == 1u)
                        {   
                            uPrvBrId = uF_OTElem[uSrcElem*u_MaxBrLvl +(u_MaxBrLvl -uTemp0 +uPrvBrLv)];
                            uTemp1   = 0u;
                            
                            while (uPrvBrLv < (uTemp0 -2))
                            {   uPrvBrLv += 1u;
                                uPrvBrId  = uF_OTElem[uSrcElem*u_MaxBrLvl +(u_MaxBrLvl -uTemp0 +uPrvBrLv)];
                                
                                memcpy(fSrc, &fKh_FFBrch[uPrvBrId*21 +0], 3u*sizeof(float));
                                SubtractVect3(fTempVect, fRcv, fSrc);               fDist  = VectorLength(fTempVect);
                                
                                uTemp1  = (uFUsdBrch[uPrvBrId]                                   == 0u)*1u     + (uFUsdBrch[uPrvBrId]                                      != 0u)*0u;
                                uTemp1  = (uF_OTPrtChld[uPrvBrId*5 +0]                            > 1u)*uTemp1 + (uF_OTPrtChld[uPrvBrId*5 +0]                              <= 1u)*0u;
                                uTemp1  = (fDist > (fKh_FFBrch[uPrvBrId*21 +14]+KH_FULL_DIST*fFltSide))*uTemp1 + (fDist <= (fKh_FFBrch[uPrvBrId*21 +14]+KH_FULL_DIST*fFltSide))*0u;
                                
                                if (uTemp1 == 1u)         {       break;      }
                            }
                            
                            if (uTemp1 == 0u)
                            {   memcpy(fP1, &fFV_temp[uF_temp[uSrcElem*7 +0]*4 +0], 3u*sizeof(float));
                                memcpy(fP2, &fFV_temp[uF_temp[uSrcElem*7 +1]*4 +0], 3u*sizeof(float));
                                memcpy(fP3, &fFV_temp[uF_temp[uSrcElem*7 +2]*4 +0], 3u*sizeof(float));
                                
                                if (USEHALFSPACE == 1u)
                                {   StrainHS_Nikkhoo(fStressStk, fStrain, fRcv[0], fRcv[1], fRcv[2], fP1, fP2, fP3, fUnitSlipF, 0.0f, 0.0f, fModPara[2], fModPara[4]);
                                    StrainHS_Nikkhoo(fStressDip, fStrain, fRcv[0], fRcv[1], fRcv[2], fP1, fP2, fP3, 0.0f, fUnitSlipF, 0.0f, fModPara[2], fModPara[4]);
                                    StrainHS_Nikkhoo(fStressNrm, fStrain, fRcv[0], fRcv[1], fRcv[2], fP1, fP2, fP3, 0.0f, 0.0f, fUnitSlipF, fModPara[2], fModPara[4]);
                                }
                                else
                                {   StrainFS_Nikkhoo(fStressStk, fStrain, fRcv[0], fRcv[1], fRcv[2], fP1, fP2, fP3, fUnitSlipF, 0.0f, 0.0f, fModPara[2], fModPara[4]);
                                    StrainFS_Nikkhoo(fStressDip, fStrain, fRcv[0], fRcv[1], fRcv[2], fP1, fP2, fP3, 0.0f, fUnitSlipF, 0.0f, fModPara[2], fModPara[4]);
                                    StrainFS_Nikkhoo(fStressNrm, fStrain, fRcv[0], fRcv[1], fRcv[2], fP1, fP2, fP3, 0.0f, 0.0f, fUnitSlipF, fModPara[2], fModPara[4]);
                                }
                                
                                ScaleStress6(fStressStk, (fDistF/fUnitSlipF));          RotStressG2L(&fPrvStress[0], fRotMat, fStressStk);
                                ScaleStress6(fStressDip, (fDistF/fUnitSlipF));          RotStressG2L(&fPrvStress[3], fRotMat, fStressDip);
                                ScaleStress6(fStressNrm, (fDistF/fUnitSlipF));          RotStressG2L(&fPrvStress[6], fRotMat, fStressNrm);
                                
                                fKh_FFvalstk[i][uKh_FFcnt[i]*2 +0] = fPrvStress[0];         fKh_FFvalstk[i][uKh_FFcnt[i]*2 +1] = fPrvStress[3];
                                fKh_FFvaldip[i][uKh_FFcnt[i]*2 +0] = fPrvStress[1];         fKh_FFvaldip[i][uKh_FFcnt[i]*2 +1] = fPrvStress[4];
                                fKh_FFvalnrm[i][uKh_FFcnt[i]*2 +0] = fPrvStress[2];         fKh_FFvalnrm[i][uKh_FFcnt[i]*2 +1] = fPrvStress[5];
                                
                                uFUsdEl[uSrcElem]       = 1u;
                                uKh_FFps2[i][uSrcElem]  = uKh_FFcnt[i];
                                uKh_FFcnt[i]           += 1u;
                                
                                uUsdNum = 0u;
                                for (j = 0u; j < uFPNum;  j++)     {       uUsdNum +=  uFUsdEl[j];         }
                                
                                for (j = 1u; j < uTemp0; j++)      {       uFUsdBrch[ uF_OTElem[uSrcElem*u_MaxBrLvl +(u_MaxBrLvl -uTemp0 +j)] ] = 1u;       }
                                uTemp2 = 0u;
                            }
                            else
                            {   
                                GetGlobVertsForRectangle(fP1, fP2, fP3, fP4, fKh_FFBrch, uPrvBrId*21u);
                                
                                if (USEHALFSPACE == 1u)
                                {   StrainHS_Nikkhoo(fStressStkT, fStrain, fRcv[0], fRcv[1], fRcv[2], fP1, fP2, fP4, fUnitSlipF, 0.0f, 0.0f, fModPara[2], fModPara[4]);
                                    StrainHS_Nikkhoo(fStressDipT, fStrain, fRcv[0], fRcv[1], fRcv[2], fP1, fP2, fP4, 0.0f, fUnitSlipF, 0.0f, fModPara[2], fModPara[4]);
                                    StrainHS_Nikkhoo(fStressNrmT, fStrain, fRcv[0], fRcv[1], fRcv[2], fP1, fP2, fP4, 0.0f, 0.0f, fUnitSlipF, fModPara[2], fModPara[4]);
                                    
                                    StrainHS_Nikkhoo(fStressStk, fStrain, fRcv[0], fRcv[1], fRcv[2], fP2, fP3, fP4, fUnitSlipF, 0.0f, 0.0f, fModPara[2], fModPara[4]);
                                    StrainHS_Nikkhoo(fStressDip, fStrain, fRcv[0], fRcv[1], fRcv[2], fP2, fP3, fP4, 0.0f, fUnitSlipF, 0.0f, fModPara[2], fModPara[4]);
                                    StrainHS_Nikkhoo(fStressNrm, fStrain, fRcv[0], fRcv[1], fRcv[2], fP2, fP3, fP4, 0.0f, 0.0f, fUnitSlipF, fModPara[2], fModPara[4]);
                                    
                                    fStressStk[0] += fStressStkT[0];        fStressStk[1] += fStressStkT[1];        fStressStk[2] += fStressStkT[2];
                                    fStressStk[3] += fStressStkT[3];        fStressStk[4] += fStressStkT[4];        fStressStk[5] += fStressStkT[5];
                                    fStressDip[0] += fStressDipT[0];        fStressDip[1] += fStressDipT[1];        fStressDip[2] += fStressDipT[2];
                                    fStressDip[3] += fStressDipT[3];        fStressDip[4] += fStressDipT[4];        fStressDip[5] += fStressDipT[5];
                                    fStressNrm[0] += fStressNrmT[0];        fStressNrm[1] += fStressNrmT[1];        fStressNrm[2] += fStressNrmT[2];
                                    fStressNrm[3] += fStressNrmT[3];        fStressNrm[4] += fStressNrmT[4];        fStressNrm[5] += fStressNrmT[5];
                                }
                                else
                                {   StrainFS_Nikkhoo(fStressStkT, fStrain, fRcv[0], fRcv[1], fRcv[2], fP1, fP2, fP4, fUnitSlipF, 0.0f, 0.0f, fModPara[2], fModPara[4]);
                                    StrainFS_Nikkhoo(fStressDipT, fStrain, fRcv[0], fRcv[1], fRcv[2], fP1, fP2, fP4, 0.0f, fUnitSlipF, 0.0f, fModPara[2], fModPara[4]);
                                    StrainFS_Nikkhoo(fStressNrmT, fStrain, fRcv[0], fRcv[1], fRcv[2], fP1, fP2, fP4, 0.0f, 0.0f, fUnitSlipF, fModPara[2], fModPara[4]);
                                    
                                    StrainFS_Nikkhoo(fStressStk, fStrain, fRcv[0], fRcv[1], fRcv[2], fP2, fP3, fP4, fUnitSlipF, 0.0f, 0.0f, fModPara[2], fModPara[4]);
                                    StrainFS_Nikkhoo(fStressDip, fStrain, fRcv[0], fRcv[1], fRcv[2], fP2, fP3, fP4, 0.0f, fUnitSlipF, 0.0f, fModPara[2], fModPara[4]);
                                    StrainFS_Nikkhoo(fStressNrm, fStrain, fRcv[0], fRcv[1], fRcv[2], fP2, fP3, fP4, 0.0f, 0.0f, fUnitSlipF, fModPara[2], fModPara[4]);
                                    
                                    fStressStk[0] += fStressStkT[0];        fStressStk[1] += fStressStkT[1];        fStressStk[2] += fStressStkT[2];
                                    fStressStk[3] += fStressStkT[3];        fStressStk[4] += fStressStkT[4];        fStressStk[5] += fStressStkT[5];
                                    fStressDip[0] += fStressDipT[0];        fStressDip[1] += fStressDipT[1];        fStressDip[2] += fStressDipT[2];
                                    fStressDip[3] += fStressDipT[3];        fStressDip[4] += fStressDipT[4];        fStressDip[5] += fStressDipT[5];
                                    fStressNrm[0] += fStressNrmT[0];        fStressNrm[1] += fStressNrmT[1];        fStressNrm[2] += fStressNrmT[2];
                                    fStressNrm[3] += fStressNrmT[3];        fStressNrm[4] += fStressNrmT[4];        fStressNrm[5] += fStressNrmT[5];
                                }
                                
                                ScaleStress6(fStressStk, (fDistF/fUnitSlipF));          ScaleStress6(fStressStk, (1.0f/fKh_FFBrch[uPrvBrId*21 +13]));       RotStressG2L(&fPrvStress[0], fRotMat, fStressStk);
                                ScaleStress6(fStressDip, (fDistF/fUnitSlipF));          ScaleStress6(fStressDip, (1.0f/fKh_FFBrch[uPrvBrId*21 +13]));       RotStressG2L(&fPrvStress[3], fRotMat, fStressDip);
                                ScaleStress6(fStressNrm, (fDistF/fUnitSlipF));          ScaleStress6(fStressNrm, (1.0f/fKh_FFBrch[uPrvBrId*21 +13]));       RotStressG2L(&fPrvStress[6], fRotMat, fStressNrm);
                                
                                uTemp1 = 1u;
                                for (j = 0; j < uF_OTPrtChld[uPrvBrId*5 +0]; j++)
                                {   uTstBrId = uF_OTPrtChld[uPrvBrId*5 +(j+1)];
                                    
                                    GetGlobVertsForRectangle(fP1, fP2, fP3, fP4, fKh_FFBrch, uTstBrId*21u);
                                    
                                    if (USEHALFSPACE == 1u)
                                    {   StrainHS_Nikkhoo(fStressStkT, fStrain, fRcv[0], fRcv[1], fRcv[2], fP1, fP2, fP4, fUnitSlipF, 0.0f, 0.0f, fModPara[2], fModPara[4]);
                                        StrainHS_Nikkhoo(fStressDipT, fStrain, fRcv[0], fRcv[1], fRcv[2], fP1, fP2, fP4, 0.0f, fUnitSlipF, 0.0f, fModPara[2], fModPara[4]);
                                        StrainHS_Nikkhoo(fStressNrmT, fStrain, fRcv[0], fRcv[1], fRcv[2], fP1, fP2, fP4, 0.0f, 0.0f, fUnitSlipF, fModPara[2], fModPara[4]);
                                        
                                        StrainHS_Nikkhoo(fStressStk, fStrain, fRcv[0], fRcv[1], fRcv[2], fP2, fP3, fP4, fUnitSlipF, 0.0f, 0.0f, fModPara[2], fModPara[4]);
                                        StrainHS_Nikkhoo(fStressDip, fStrain, fRcv[0], fRcv[1], fRcv[2], fP2, fP3, fP4, 0.0f, fUnitSlipF, 0.0f, fModPara[2], fModPara[4]);
                                        StrainHS_Nikkhoo(fStressNrm, fStrain, fRcv[0], fRcv[1], fRcv[2], fP2, fP3, fP4, 0.0f, 0.0f, fUnitSlipF, fModPara[2], fModPara[4]);
                                        
                                        fStressStk[0] += fStressStkT[0];        fStressStk[1] += fStressStkT[1];        fStressStk[2] += fStressStkT[2];
                                        fStressStk[3] += fStressStkT[3];        fStressStk[4] += fStressStkT[4];        fStressStk[5] += fStressStkT[5];
                                        fStressDip[0] += fStressDipT[0];        fStressDip[1] += fStressDipT[1];        fStressDip[2] += fStressDipT[2];
                                        fStressDip[3] += fStressDipT[3];        fStressDip[4] += fStressDipT[4];        fStressDip[5] += fStressDipT[5];
                                        fStressNrm[0] += fStressNrmT[0];        fStressNrm[1] += fStressNrmT[1];        fStressNrm[2] += fStressNrmT[2];
                                        fStressNrm[3] += fStressNrmT[3];        fStressNrm[4] += fStressNrmT[4];        fStressNrm[5] += fStressNrmT[5];
                                    }
                                    else
                                    {   StrainFS_Nikkhoo(fStressStkT, fStrain, fRcv[0], fRcv[1], fRcv[2], fP1, fP2, fP4, fUnitSlipF, 0.0f, 0.0f, fModPara[2], fModPara[4]);
                                        StrainFS_Nikkhoo(fStressDipT, fStrain, fRcv[0], fRcv[1], fRcv[2], fP1, fP2, fP4, 0.0f, fUnitSlipF, 0.0f, fModPara[2], fModPara[4]);
                                        StrainFS_Nikkhoo(fStressNrmT, fStrain, fRcv[0], fRcv[1], fRcv[2], fP1, fP2, fP4, 0.0f, 0.0f, fUnitSlipF, fModPara[2], fModPara[4]);
                                        
                                        StrainFS_Nikkhoo(fStressStk, fStrain, fRcv[0], fRcv[1], fRcv[2], fP2, fP3, fP4, fUnitSlipF, 0.0f, 0.0f, fModPara[2], fModPara[4]);
                                        StrainFS_Nikkhoo(fStressDip, fStrain, fRcv[0], fRcv[1], fRcv[2], fP2, fP3, fP4, 0.0f, fUnitSlipF, 0.0f, fModPara[2], fModPara[4]);
                                        StrainFS_Nikkhoo(fStressNrm, fStrain, fRcv[0], fRcv[1], fRcv[2], fP2, fP3, fP4, 0.0f, 0.0f, fUnitSlipF, fModPara[2], fModPara[4]);
                                        
                                        fStressStk[0] += fStressStkT[0];        fStressStk[1] += fStressStkT[1];        fStressStk[2] += fStressStkT[2];
                                        fStressStk[3] += fStressStkT[3];        fStressStk[4] += fStressStkT[4];        fStressStk[5] += fStressStkT[5];
                                        fStressDip[0] += fStressDipT[0];        fStressDip[1] += fStressDipT[1];        fStressDip[2] += fStressDipT[2];
                                        fStressDip[3] += fStressDipT[3];        fStressDip[4] += fStressDipT[4];        fStressDip[5] += fStressDipT[5];
                                        fStressNrm[0] += fStressNrmT[0];        fStressNrm[1] += fStressNrmT[1];        fStressNrm[2] += fStressNrmT[2];
                                        fStressNrm[3] += fStressNrmT[3];        fStressNrm[4] += fStressNrmT[4];        fStressNrm[5] += fStressNrmT[5];
                                    }
                                    
                                    ScaleStress6(fStressStk, (fDistF/fUnitSlipF));          ScaleStress6(fStressStk, (1.0f/fKh_FFBrch[uTstBrId*21 +13]));       RotStressG2L(&fTstStress[0], fRotMat, fStressStk);
                                    ScaleStress6(fStressDip, (fDistF/fUnitSlipF));          ScaleStress6(fStressDip, (1.0f/fKh_FFBrch[uTstBrId*21 +13]));       RotStressG2L(&fTstStress[3], fRotMat, fStressDip);
                                    ScaleStress6(fStressNrm, (fDistF/fUnitSlipF));          ScaleStress6(fStressNrm, (1.0f/fKh_FFBrch[uTstBrId*21 +13]));       RotStressG2L(&fTstStress[6], fRotMat, fStressNrm);
                                    
                                    memcpy(fSrc, &fKh_FFBrch[uTstBrId*21 +0], 3u*sizeof(float));
                                    SubtractVect3(fTempVect, fRcv, fSrc);
                                    fDist  = VectorLength(fTempVect);
                                    fTemp2 = (fDist/fFltSide - KH_FULL_DIST)/(KH_TOL_DIST - KH_FULL_DIST);
                                    fTemp2 = (fTemp2 >= 0.0f)*fTemp2 + (fTemp2 < 0.0f)*0.0f;
                                    fTemp2 = (fTemp2 <= 1.0f)*fTemp2 + (fTemp2 > 1.0f)*1.0f;
                                    fTemp2 = KH_TOL1 + (KH_TOL2-KH_TOL1)*fTemp2;
                                    
                                    fTemp0 = sqrtf( (fTstStress[0]-fPrvStress[0])*(fTstStress[0]-fPrvStress[0]) + (fTstStress[1]-fPrvStress[1])*(fTstStress[1]-fPrvStress[1]) + (fTstStress[2]-fPrvStress[2])*(fTstStress[2]-fPrvStress[2])  + (fTstStress[3]-fPrvStress[3])*(fTstStress[3]-fPrvStress[3]) + (fTstStress[4]-fPrvStress[4])*(fTstStress[4]-fPrvStress[4]) + (fTstStress[5]-fPrvStress[5])*(fTstStress[5]-fPrvStress[5])   + (fTstStress[6]-fPrvStress[6])*(fTstStress[6]-fPrvStress[6]) + (fTstStress[7]-fPrvStress[7])*(fTstStress[7]-fPrvStress[7]) + (fTstStress[8]-fPrvStress[8])*(fTstStress[8]-fPrvStress[8]));
                                    fTemp1 = sqrtf(                fPrvStress[0] *fPrvStress[0]                 +                fPrvStress[1] *fPrvStress[1]                 +                fPrvStress[2] *fPrvStress[2]                  +                fPrvStress[3] *fPrvStress[3]                 +                fPrvStress[4] *fPrvStress[4]                 +                fPrvStress[5] *fPrvStress[5]                   +                fPrvStress[6] *fPrvStress[6]                 +                fPrvStress[7] *fPrvStress[7]                 +                fPrvStress[8] *fPrvStress[8]);
                                    
                                    uTemp1 = ((fTemp0/fTemp1) <= fTemp2)*uTemp1 + ((fTemp0/fTemp1) > fTemp2)*0u; 
                                }
                                
                                if (uTemp1 == 1u)
                                {   fKh_FFvalstk[i][uKh_FFcnt[i]*2 +0] = fPrvStress[0];         fKh_FFvalstk[i][uKh_FFcnt[i]*2 +1] = fPrvStress[3];
                                    fKh_FFvaldip[i][uKh_FFcnt[i]*2 +0] = fPrvStress[1];         fKh_FFvaldip[i][uKh_FFcnt[i]*2 +1] = fPrvStress[4];
                                    fKh_FFvalnrm[i][uKh_FFcnt[i]*2 +0] = fPrvStress[2];         fKh_FFvalnrm[i][uKh_FFcnt[i]*2 +1] = fPrvStress[5];
                                    
                                    uUsdNum = 0u;
                                    for (j = 0u; j < uFPNum;  j++)
                                    {   uFUsdEl[j]      = (uF_OTElem[j*u_MaxBrLvl +(u_MaxBrLvl -uTemp0 +uPrvBrLv)] == uPrvBrId)*1u           + (uF_OTElem[j*u_MaxBrLvl +(u_MaxBrLvl -uTemp0 +uPrvBrLv)] != uPrvBrId)*uFUsdEl[j];
                                        uKh_FFps2[i][j] = (uF_OTElem[j*u_MaxBrLvl +(u_MaxBrLvl -uTemp0 +uPrvBrLv)] == uPrvBrId)*uKh_FFcnt[i] + (uF_OTElem[j*u_MaxBrLvl +(u_MaxBrLvl -uTemp0 +uPrvBrLv)] != uPrvBrId)*uKh_FFps2[i][j];
                                        uUsdNum        +=  uFUsdEl[j];
                                    }
                                    uKh_FFcnt[i] += 1u;
                                    
                                    for (j = 1u; j < uTemp0; j++)       {       uFUsdBrch[ uF_OTElem[uSrcElem*u_MaxBrLvl +(u_MaxBrLvl -uTemp0 +j)] ] = 1u;       }
                                    uTemp2 = 0u;
                }   }   }   }   }
                
                fKh_FFvalstk[i] = realloc(fKh_FFvalstk[i], 2*uKh_FFcnt[i]* sizeof *fKh_FFvalstk[i] );
                fKh_FFvaldip[i] = realloc(fKh_FFvaldip[i], 2*uKh_FFcnt[i]* sizeof *fKh_FFvaldip[i] );
                fKh_FFvalnrm[i] = realloc(fKh_FFvalnrm[i], 2*uKh_FFcnt[i]* sizeof *fKh_FFvalnrm[i] );
                
                if ((uPlotCatalog2Screen == 1)&&(iRANK == 0)) {    fprintf(stdout,"\nfault receiver:  %u/%u   with Flt/Flt:  %u of %u",i, iFOFFSET[iRANK], uKh_FFcnt[i], uFPNum);                  }
                
                
                uKh_FBps2[i]    = calloc(     uBPNum, sizeof *uKh_FBps2[i] );
                fKh_FBvalstk[i] = calloc(   3*uBPNum, sizeof *fKh_FBvalstk[i] );
                fKh_FBvaldip[i] = calloc(   3*uBPNum, sizeof *fKh_FBvaldip[i] );
                
                uKh_FBcnt[i]  = 0u;
                uUsdNum       = 0u;
                memset(uBUsdEl, 0,        uBPNum*sizeof(unsigned int));
                memset(uBUsdBrch, 0,uB_TotalOffs*sizeof(unsigned int));
                uTstBrId = 0u;      uSrcElem = 0u;
                
                while (uUsdNum < uBPNum) 
                {   
                    for (j = 0u; j < uBPNum; j++)               {       if (uBUsdEl[j] == 0u)    {      uSrcElem = j;       break;          }           } 
                    
                    memcpy(fSrc, &fB_temp[uSrcElem*16 +0], 3u*sizeof(float));
                    
                    SubtractVect3(fTempVect, fRcv, fSrc);               fDist  = VectorLength(fTempVect);
                    
                    fDistF = (fDist >= CUTDISTANCE*fBndSide)*1.0f + (fDist < CUTDISTANCE*fBndSide)*0.0f; 
                    uF_temp[uGlobPos*7 +6] = (fDistF == 1.0f)*uF_temp[uGlobPos*7 +6] + (fDistF != 1.0f)*1u;
                    
                    uTemp0 = uB_OTElem[uSrcElem*u_MaxBrLvl +0];
                    uTemp1 = 0u;
                    
                    for (j = 1u; j < uTemp0; j++)
                    {   uPrvBrLv = j;
                        uPrvBrId = uB_OTElem[uSrcElem*u_MaxBrLvl +(u_MaxBrLvl -uTemp0 +uPrvBrLv)];
                        
                        memcpy(fSrc, &fKh_BBBrch[uPrvBrId*21 +0], 3u*sizeof(float));
                        SubtractVect3(fTempVect, fRcv, fSrc);               fDist  = VectorLength(fTempVect);
                        
                        uTemp1  = (uBUsdBrch[uPrvBrId]                                   == 0u)*1u     + (uBUsdBrch[uPrvBrId]                                      != 0u)*0u;
                        uTemp1  = (uB_OTPrtChld[uPrvBrId*5 +0]                            > 1u)*uTemp1 + (uB_OTPrtChld[uPrvBrId*5 +0]                              <= 1u)*0u;
                        uTemp1  = (fDist > (fKh_BBBrch[uPrvBrId*21 +14]+KH_FULL_DIST*fBndSide))*uTemp1 + (fDist <= (fKh_BBBrch[uPrvBrId*21 +14]+KH_FULL_DIST*fBndSide))*0u;
                        
                        if (uTemp1 == 1u)         {       break;      }
                    } 
                    
                    if (uTemp1 == 0u) 
                    {   memcpy(fP1, &fBV_temp[uB_temp[uSrcElem*4 +0]*3 +0], 3u*sizeof(float));
                        memcpy(fP2, &fBV_temp[uB_temp[uSrcElem*4 +1]*3 +0], 3u*sizeof(float));
                        memcpy(fP3, &fBV_temp[uB_temp[uSrcElem*4 +2]*3 +0], 3u*sizeof(float));
                        
                        if (USEHALFSPACE == 1u)
                        {   StrainHS_Nikkhoo(fStressStk, fStrain, fRcv[0], fRcv[1], fRcv[2], fP1, fP2, fP3, fUnitSlipB, 0.0f, 0.0f, fModPara[2], fModPara[4]);
                            StrainHS_Nikkhoo(fStressDip, fStrain, fRcv[0], fRcv[1], fRcv[2], fP1, fP2, fP3, 0.0f, fUnitSlipB, 0.0f, fModPara[2], fModPara[4]);
                            StrainHS_Nikkhoo(fStressNrm, fStrain, fRcv[0], fRcv[1], fRcv[2], fP1, fP2, fP3, 0.0f, 0.0f, fUnitSlipB, fModPara[2], fModPara[4]);
                        }
                        else
                        {   StrainFS_Nikkhoo(fStressStk, fStrain, fRcv[0], fRcv[1], fRcv[2], fP1, fP2, fP3, fUnitSlipB, 0.0f, 0.0f, fModPara[2], fModPara[4]);
                            StrainFS_Nikkhoo(fStressDip, fStrain, fRcv[0], fRcv[1], fRcv[2], fP1, fP2, fP3, 0.0f, fUnitSlipB, 0.0f, fModPara[2], fModPara[4]);
                            StrainFS_Nikkhoo(fStressNrm, fStrain, fRcv[0], fRcv[1], fRcv[2], fP1, fP2, fP3, 0.0f, 0.0f, fUnitSlipB, fModPara[2], fModPara[4]);
                        }
                        
                        ScaleStress6(fStressStk, (fDistF/fUnitSlipB));          RotStressG2L(&fPrvStress[0], fRotMat, fStressStk);
                        ScaleStress6(fStressDip, (fDistF/fUnitSlipB));          RotStressG2L(&fPrvStress[3], fRotMat, fStressDip);
                        ScaleStress6(fStressNrm, (fDistF/fUnitSlipB));          RotStressG2L(&fPrvStress[6], fRotMat, fStressNrm);
                        
                        fKh_FBvalstk[i][uKh_FBcnt[i]*3 +0] = fPrvStress[0];         fKh_FBvalstk[i][uKh_FBcnt[i]*3 +1] = fPrvStress[3];       fKh_FBvalstk[i][uKh_FBcnt[i]*3 +2] = fPrvStress[6];
                        fKh_FBvaldip[i][uKh_FBcnt[i]*3 +0] = fPrvStress[1];         fKh_FBvaldip[i][uKh_FBcnt[i]*3 +1] = fPrvStress[4];       fKh_FBvaldip[i][uKh_FBcnt[i]*3 +2] = fPrvStress[7];
                        
                        uBUsdEl[uSrcElem]       = 1u;
                        uKh_FBps2[i][uSrcElem]  = uKh_FBcnt[i];
                        uKh_FBcnt[i]           += 1u;
                        
                        uUsdNum = 0u;
                        for (j = 0u; j < uBPNum;  j++)     {       uUsdNum +=  uBUsdEl[j];                                     }
                        
                        for (j = 1u; j < uTemp0; j++)      {       uBUsdBrch[ uB_OTElem[uSrcElem*u_MaxBrLvl +(u_MaxBrLvl -uTemp0 +j)] ] = 1u;       }
                    }
                    else
                    {   uPrvBrLv -= 1u;
                        uTemp2    = 1u;
                        while (uTemp2 == 1u)
                        {   
                            uPrvBrId = uB_OTElem[uSrcElem*u_MaxBrLvl +(u_MaxBrLvl -uTemp0 +uPrvBrLv)];
                            uTemp1   = 0u;
                            
                            while (uPrvBrLv < (uTemp0 -2))
                            {   uPrvBrLv += 1u;
                                uPrvBrId  = uB_OTElem[uSrcElem*u_MaxBrLvl +(u_MaxBrLvl -uTemp0 +uPrvBrLv)];
                                
                                memcpy(fSrc, &fKh_BBBrch[uPrvBrId*21 +0], 3u*sizeof(float));
                                SubtractVect3(fTempVect, fRcv, fSrc);               fDist  = VectorLength(fTempVect);
                                
                                uTemp1  = (uBUsdBrch[uPrvBrId]                                   == 0u)*1u     + (uBUsdBrch[uPrvBrId]                                      != 0u)*0u;
                                uTemp1  = (uB_OTPrtChld[uPrvBrId*5 +0]                            > 1u)*uTemp1 + (uB_OTPrtChld[uPrvBrId*5 +0]                              <= 1u)*0u;
                                uTemp1  = (fDist > (fKh_BBBrch[uPrvBrId*21 +14]+KH_FULL_DIST*fBndSide))*uTemp1 + (fDist <= (fKh_BBBrch[uPrvBrId*21 +14]+KH_FULL_DIST*fBndSide))*0u;
                                
                                if (uTemp1 == 1u)         {       break;      }
                            }
                            
                            if (uTemp1 == 0u)
                            {   memcpy(fP1, &fBV_temp[uB_temp[uSrcElem*4 +0]*3 +0], 3u*sizeof(float));
                                memcpy(fP2, &fBV_temp[uB_temp[uSrcElem*4 +1]*3 +0], 3u*sizeof(float));
                                memcpy(fP3, &fBV_temp[uB_temp[uSrcElem*4 +2]*3 +0], 3u*sizeof(float));
                                
                                if (USEHALFSPACE == 1u)
                                {   StrainHS_Nikkhoo(fStressStk, fStrain, fRcv[0], fRcv[1], fRcv[2], fP1, fP2, fP3, fUnitSlipB, 0.0f, 0.0f, fModPara[2], fModPara[4]);
                                    StrainHS_Nikkhoo(fStressDip, fStrain, fRcv[0], fRcv[1], fRcv[2], fP1, fP2, fP3, 0.0f, fUnitSlipB, 0.0f, fModPara[2], fModPara[4]);
                                    StrainHS_Nikkhoo(fStressNrm, fStrain, fRcv[0], fRcv[1], fRcv[2], fP1, fP2, fP3, 0.0f, 0.0f, fUnitSlipB, fModPara[2], fModPara[4]);
                                }
                                else
                                {   StrainFS_Nikkhoo(fStressStk, fStrain, fRcv[0], fRcv[1], fRcv[2], fP1, fP2, fP3, fUnitSlipB, 0.0f, 0.0f, fModPara[2], fModPara[4]);
                                    StrainFS_Nikkhoo(fStressDip, fStrain, fRcv[0], fRcv[1], fRcv[2], fP1, fP2, fP3, 0.0f, fUnitSlipB, 0.0f, fModPara[2], fModPara[4]);
                                    StrainFS_Nikkhoo(fStressNrm, fStrain, fRcv[0], fRcv[1], fRcv[2], fP1, fP2, fP3, 0.0f, 0.0f, fUnitSlipB, fModPara[2], fModPara[4]);
                                }
                                
                                ScaleStress6(fStressStk, (fDistF/fUnitSlipB));          RotStressG2L(&fPrvStress[0], fRotMat, fStressStk);
                                ScaleStress6(fStressDip, (fDistF/fUnitSlipB));          RotStressG2L(&fPrvStress[3], fRotMat, fStressDip);
                                ScaleStress6(fStressNrm, (fDistF/fUnitSlipB));          RotStressG2L(&fPrvStress[6], fRotMat, fStressNrm);
                                
                                fKh_FBvalstk[i][uKh_FBcnt[i]*3 +0] = fPrvStress[0];         fKh_FBvalstk[i][uKh_FBcnt[i]*3 +1] = fPrvStress[3];       fKh_FBvalstk[i][uKh_FBcnt[i]*3 +2] = fPrvStress[6];
                                fKh_FBvaldip[i][uKh_FBcnt[i]*3 +0] = fPrvStress[1];         fKh_FBvaldip[i][uKh_FBcnt[i]*3 +1] = fPrvStress[4];       fKh_FBvaldip[i][uKh_FBcnt[i]*3 +2] = fPrvStress[7];
                               
                                uBUsdEl[uSrcElem]       = 1u;
                                uKh_FBps2[i][uSrcElem]  = uKh_FBcnt[i];
                                uKh_FBcnt[i]           += 1u;
                                
                                uUsdNum = 0u;
                                for (j = 0u; j < uBPNum;  j++)      {       uUsdNum +=  uBUsdEl[j];         }
                                
                                for (j = 1u; j < uTemp0; j++)       {       uBUsdBrch[ uB_OTElem[uSrcElem*u_MaxBrLvl +(u_MaxBrLvl -uTemp0 +j)] ] = 1u;       }
                                uTemp2 = 0u;
                            }
                            else
                            {  
                                
                                GetGlobVertsForRectangle(fP1, fP2, fP3, fP4, fKh_BBBrch, uPrvBrId*21u);
                                
                                if (USEHALFSPACE == 1u)
                                {   StrainHS_Nikkhoo(fStressStkT, fStrain, fRcv[0], fRcv[1], fRcv[2], fP1, fP2, fP4, fUnitSlipB, 0.0f, 0.0f, fModPara[2], fModPara[4]);
                                    StrainHS_Nikkhoo(fStressDipT, fStrain, fRcv[0], fRcv[1], fRcv[2], fP1, fP2, fP4, 0.0f, fUnitSlipB, 0.0f, fModPara[2], fModPara[4]);
                                    StrainHS_Nikkhoo(fStressNrmT, fStrain, fRcv[0], fRcv[1], fRcv[2], fP1, fP2, fP4, 0.0f, 0.0f, fUnitSlipB, fModPara[2], fModPara[4]);
                                    
                                    StrainHS_Nikkhoo(fStressStk, fStrain, fRcv[0], fRcv[1], fRcv[2], fP2, fP3, fP4, fUnitSlipB, 0.0f, 0.0f, fModPara[2], fModPara[4]);
                                    StrainHS_Nikkhoo(fStressDip, fStrain, fRcv[0], fRcv[1], fRcv[2], fP2, fP3, fP4, 0.0f, fUnitSlipB, 0.0f, fModPara[2], fModPara[4]);
                                    StrainHS_Nikkhoo(fStressNrm, fStrain, fRcv[0], fRcv[1], fRcv[2], fP2, fP3, fP4, 0.0f, 0.0f, fUnitSlipB, fModPara[2], fModPara[4]);
                                    
                                    fStressStk[0] += fStressStkT[0];        fStressStk[1] += fStressStkT[1];        fStressStk[2] += fStressStkT[2];
                                    fStressStk[3] += fStressStkT[3];        fStressStk[4] += fStressStkT[4];        fStressStk[5] += fStressStkT[5];
                                    fStressDip[0] += fStressDipT[0];        fStressDip[1] += fStressDipT[1];        fStressDip[2] += fStressDipT[2];
                                    fStressDip[3] += fStressDipT[3];        fStressDip[4] += fStressDipT[4];        fStressDip[5] += fStressDipT[5];
                                    fStressNrm[0] += fStressNrmT[0];        fStressNrm[1] += fStressNrmT[1];        fStressNrm[2] += fStressNrmT[2];
                                    fStressNrm[3] += fStressNrmT[3];        fStressNrm[4] += fStressNrmT[4];        fStressNrm[5] += fStressNrmT[5];
                                }
                                else
                                {   StrainFS_Nikkhoo(fStressStkT, fStrain, fRcv[0], fRcv[1], fRcv[2], fP1, fP2, fP4, fUnitSlipB, 0.0f, 0.0f, fModPara[2], fModPara[4]);
                                    StrainFS_Nikkhoo(fStressDipT, fStrain, fRcv[0], fRcv[1], fRcv[2], fP1, fP2, fP4, 0.0f, fUnitSlipB, 0.0f, fModPara[2], fModPara[4]);
                                    StrainFS_Nikkhoo(fStressNrmT, fStrain, fRcv[0], fRcv[1], fRcv[2], fP1, fP2, fP4, 0.0f, 0.0f, fUnitSlipB, fModPara[2], fModPara[4]);
                                    
                                    StrainFS_Nikkhoo(fStressStk, fStrain, fRcv[0], fRcv[1], fRcv[2], fP2, fP3, fP4, fUnitSlipB, 0.0f, 0.0f, fModPara[2], fModPara[4]);
                                    StrainFS_Nikkhoo(fStressDip, fStrain, fRcv[0], fRcv[1], fRcv[2], fP2, fP3, fP4, 0.0f, fUnitSlipB, 0.0f, fModPara[2], fModPara[4]);
                                    StrainFS_Nikkhoo(fStressNrm, fStrain, fRcv[0], fRcv[1], fRcv[2], fP2, fP3, fP4, 0.0f, 0.0f, fUnitSlipB, fModPara[2], fModPara[4]);
                                    
                                    fStressStk[0] += fStressStkT[0];        fStressStk[1] += fStressStkT[1];        fStressStk[2] += fStressStkT[2];
                                    fStressStk[3] += fStressStkT[3];        fStressStk[4] += fStressStkT[4];        fStressStk[5] += fStressStkT[5];
                                    fStressDip[0] += fStressDipT[0];        fStressDip[1] += fStressDipT[1];        fStressDip[2] += fStressDipT[2];
                                    fStressDip[3] += fStressDipT[3];        fStressDip[4] += fStressDipT[4];        fStressDip[5] += fStressDipT[5];
                                    fStressNrm[0] += fStressNrmT[0];        fStressNrm[1] += fStressNrmT[1];        fStressNrm[2] += fStressNrmT[2];
                                    fStressNrm[3] += fStressNrmT[3];        fStressNrm[4] += fStressNrmT[4];        fStressNrm[5] += fStressNrmT[5];
                                }
                                
                                ScaleStress6(fStressStk, (fDistF/fUnitSlipB));          ScaleStress6(fStressStk, (1.0f/fKh_BBBrch[uPrvBrId*21 +13]));           RotStressG2L(&fPrvStress[0], fRotMat, fStressStk);
                                ScaleStress6(fStressDip, (fDistF/fUnitSlipB));          ScaleStress6(fStressDip, (1.0f/fKh_BBBrch[uPrvBrId*21 +13]));           RotStressG2L(&fPrvStress[3], fRotMat, fStressDip);
                                ScaleStress6(fStressNrm, (fDistF/fUnitSlipB));          ScaleStress6(fStressNrm, (1.0f/fKh_BBBrch[uPrvBrId*21 +13]));           RotStressG2L(&fPrvStress[6], fRotMat, fStressNrm);
                                
                                uTemp1 = 1u;
                                for (j = 0; j < uB_OTPrtChld[uPrvBrId*5 +0]; j++)
                                {   uTstBrId = uB_OTPrtChld[uPrvBrId*5 +(j+1)];
                                    
                                    GetGlobVertsForRectangle(fP1, fP2, fP3, fP4, fKh_BBBrch, uTstBrId*21u);
                                    
                                    if (USEHALFSPACE == 1u)
                                    {   StrainHS_Nikkhoo(fStressStkT, fStrain, fRcv[0], fRcv[1], fRcv[2], fP1, fP2, fP4, fUnitSlipB, 0.0f, 0.0f, fModPara[2], fModPara[4]);
                                        StrainHS_Nikkhoo(fStressDipT, fStrain, fRcv[0], fRcv[1], fRcv[2], fP1, fP2, fP4, 0.0f, fUnitSlipB, 0.0f, fModPara[2], fModPara[4]);
                                        StrainHS_Nikkhoo(fStressNrmT, fStrain, fRcv[0], fRcv[1], fRcv[2], fP1, fP2, fP4, 0.0f, 0.0f, fUnitSlipB, fModPara[2], fModPara[4]);
                                        
                                        StrainHS_Nikkhoo(fStressStk, fStrain, fRcv[0], fRcv[1], fRcv[2], fP2, fP3, fP4, fUnitSlipB, 0.0f, 0.0f, fModPara[2], fModPara[4]);
                                        StrainHS_Nikkhoo(fStressDip, fStrain, fRcv[0], fRcv[1], fRcv[2], fP2, fP3, fP4, 0.0f, fUnitSlipB, 0.0f, fModPara[2], fModPara[4]);
                                        StrainHS_Nikkhoo(fStressNrm, fStrain, fRcv[0], fRcv[1], fRcv[2], fP2, fP3, fP4, 0.0f, 0.0f, fUnitSlipB, fModPara[2], fModPara[4]);
                                        
                                        fStressStk[0] += fStressStkT[0];        fStressStk[1] += fStressStkT[1];        fStressStk[2] += fStressStkT[2];
                                        fStressStk[3] += fStressStkT[3];        fStressStk[4] += fStressStkT[4];        fStressStk[5] += fStressStkT[5];
                                        fStressDip[0] += fStressDipT[0];        fStressDip[1] += fStressDipT[1];        fStressDip[2] += fStressDipT[2];
                                        fStressDip[3] += fStressDipT[3];        fStressDip[4] += fStressDipT[4];        fStressDip[5] += fStressDipT[5];
                                        fStressNrm[0] += fStressNrmT[0];        fStressNrm[1] += fStressNrmT[1];        fStressNrm[2] += fStressNrmT[2];
                                        fStressNrm[3] += fStressNrmT[3];        fStressNrm[4] += fStressNrmT[4];        fStressNrm[5] += fStressNrmT[5];
                                    }
                                    else
                                    {   StrainFS_Nikkhoo(fStressStkT, fStrain, fRcv[0], fRcv[1], fRcv[2], fP1, fP2, fP4, fUnitSlipB, 0.0f, 0.0f, fModPara[2], fModPara[4]);
                                        StrainFS_Nikkhoo(fStressDipT, fStrain, fRcv[0], fRcv[1], fRcv[2], fP1, fP2, fP4, 0.0f, fUnitSlipB, 0.0f, fModPara[2], fModPara[4]);
                                        StrainFS_Nikkhoo(fStressNrmT, fStrain, fRcv[0], fRcv[1], fRcv[2], fP1, fP2, fP4, 0.0f, 0.0f, fUnitSlipB, fModPara[2], fModPara[4]);
                                        
                                        StrainFS_Nikkhoo(fStressStk, fStrain, fRcv[0], fRcv[1], fRcv[2], fP2, fP3, fP4, fUnitSlipB, 0.0f, 0.0f, fModPara[2], fModPara[4]);
                                        StrainFS_Nikkhoo(fStressDip, fStrain, fRcv[0], fRcv[1], fRcv[2], fP2, fP3, fP4, 0.0f, fUnitSlipB, 0.0f, fModPara[2], fModPara[4]);
                                        StrainFS_Nikkhoo(fStressNrm, fStrain, fRcv[0], fRcv[1], fRcv[2], fP2, fP3, fP4, 0.0f, 0.0f, fUnitSlipB, fModPara[2], fModPara[4]);
                                        
                                        fStressStk[0] += fStressStkT[0];        fStressStk[1] += fStressStkT[1];        fStressStk[2] += fStressStkT[2];
                                        fStressStk[3] += fStressStkT[3];        fStressStk[4] += fStressStkT[4];        fStressStk[5] += fStressStkT[5];
                                        fStressDip[0] += fStressDipT[0];        fStressDip[1] += fStressDipT[1];        fStressDip[2] += fStressDipT[2];
                                        fStressDip[3] += fStressDipT[3];        fStressDip[4] += fStressDipT[4];        fStressDip[5] += fStressDipT[5];
                                        fStressNrm[0] += fStressNrmT[0];        fStressNrm[1] += fStressNrmT[1];        fStressNrm[2] += fStressNrmT[2];
                                        fStressNrm[3] += fStressNrmT[3];        fStressNrm[4] += fStressNrmT[4];        fStressNrm[5] += fStressNrmT[5];
                                    }
                                    
                                    ScaleStress6(fStressStk, (fDistF/fUnitSlipB));          ScaleStress6(fStressStk, (1.0f/fKh_BBBrch[uTstBrId*21 +13]));           RotStressG2L(&fTstStress[0], fRotMat, fStressStk);
                                    ScaleStress6(fStressDip, (fDistF/fUnitSlipB));          ScaleStress6(fStressDip, (1.0f/fKh_BBBrch[uTstBrId*21 +13]));           RotStressG2L(&fTstStress[3], fRotMat, fStressDip);
                                    ScaleStress6(fStressNrm, (fDistF/fUnitSlipB));          ScaleStress6(fStressNrm, (1.0f/fKh_BBBrch[uTstBrId*21 +13]));           RotStressG2L(&fTstStress[6], fRotMat, fStressNrm);
                                    
                                    memcpy(fSrc, &fKh_BBBrch[uTstBrId*21 +0], 3u*sizeof(float));
                                    SubtractVect3(fTempVect, fRcv, fSrc);
                                    fDist  = VectorLength(fTempVect);
                                    fTemp2 = (fDist/fBndSide - KH_FULL_DIST)/(KH_TOL_DIST - KH_FULL_DIST);
                                    fTemp2 = (fTemp2 >= 0.0f)*fTemp2 + (fTemp2 < 0.0f)*0.0f;
                                    fTemp2 = (fTemp2 <= 1.0f)*fTemp2 + (fTemp2 > 1.0f)*1.0f;
                                    fTemp2 = KH_TOL1 + (KH_TOL2-KH_TOL1)*fTemp2;
                                    
                                    fTemp0 = sqrtf( (fTstStress[0]-fPrvStress[0])*(fTstStress[0]-fPrvStress[0]) + (fTstStress[1]-fPrvStress[1])*(fTstStress[1]-fPrvStress[1]) + (fTstStress[2]-fPrvStress[2])*(fTstStress[2]-fPrvStress[2])  + (fTstStress[3]-fPrvStress[3])*(fTstStress[3]-fPrvStress[3]) + (fTstStress[4]-fPrvStress[4])*(fTstStress[4]-fPrvStress[4]) + (fTstStress[5]-fPrvStress[5])*(fTstStress[5]-fPrvStress[5])   + (fTstStress[6]-fPrvStress[6])*(fTstStress[6]-fPrvStress[6]) + (fTstStress[7]-fPrvStress[7])*(fTstStress[7]-fPrvStress[7]) + (fTstStress[8]-fPrvStress[8])*(fTstStress[8]-fPrvStress[8]));
                                    fTemp1 = sqrtf(                fPrvStress[0] *fPrvStress[0]                 +                fPrvStress[1] *fPrvStress[1]                 +                fPrvStress[2] *fPrvStress[2]                  +                fPrvStress[3] *fPrvStress[3]                 +                fPrvStress[4] *fPrvStress[4]                 +                fPrvStress[5] *fPrvStress[5]                   +                fPrvStress[6] *fPrvStress[6]                 +                fPrvStress[7] *fPrvStress[7]                 +                fPrvStress[8] *fPrvStress[8]);
                                    
                                    uTemp1 = ((fTemp0/fTemp1) <= fTemp2)*uTemp1 + ((fTemp0/fTemp1) > fTemp2)*0u; 
                                }
                                
                                if (uTemp1 == 1u)
                                {   fKh_FBvalstk[i][uKh_FBcnt[i]*3 +0] = fPrvStress[0];         fKh_FBvalstk[i][uKh_FBcnt[i]*3 +1] = fPrvStress[3];       fKh_FBvalstk[i][uKh_FBcnt[i]*3 +2] = fPrvStress[6];
                                    fKh_FBvaldip[i][uKh_FBcnt[i]*3 +0] = fPrvStress[1];         fKh_FBvaldip[i][uKh_FBcnt[i]*3 +1] = fPrvStress[4];       fKh_FBvaldip[i][uKh_FBcnt[i]*3 +2] = fPrvStress[7];
                                    
                                    uUsdNum = 0u;
                                    for (j = 0u; j < uBPNum;  j++)
                                    {   uBUsdEl[j]      = (uB_OTElem[j*u_MaxBrLvl +(u_MaxBrLvl -uTemp0 +uPrvBrLv)] == uPrvBrId)*1u           + (uB_OTElem[j*u_MaxBrLvl +(u_MaxBrLvl -uTemp0 +uPrvBrLv)] != uPrvBrId)*uBUsdEl[j];
                                        uKh_FBps2[i][j] = (uB_OTElem[j*u_MaxBrLvl +(u_MaxBrLvl -uTemp0 +uPrvBrLv)] == uPrvBrId)*uKh_FBcnt[i] + (uB_OTElem[j*u_MaxBrLvl +(u_MaxBrLvl -uTemp0 +uPrvBrLv)] != uPrvBrId)*uKh_FBps2[i][j];
                                        uUsdNum        +=  uBUsdEl[j];
                                    }
                                    uKh_FBcnt[i] += 1u;
                                    
                                    for (j = 1u; j < uTemp0; j++)       {       uBUsdBrch[ uB_OTElem[uSrcElem*u_MaxBrLvl +(u_MaxBrLvl -uTemp0 +j)] ] = 1u;       }
                                    uTemp2 = 0u;
                }   }   }   }   }
                
                fKh_FBvalstk[i] = realloc(fKh_FBvalstk[i], 3*uKh_FBcnt[i]* sizeof *fKh_FBvalstk[i] );
                fKh_FBvaldip[i] = realloc(fKh_FBvaldip[i], 3*uKh_FBcnt[i]* sizeof *fKh_FBvaldip[i] );
                if ((uPlotCatalog2Screen == 1)&&(iRANK == 0))             {    fprintf(stdout,"   and Flt/Bnd:  %u/%u ",uKh_FBcnt[i], uBPNum);               }
            }
            
            
            for (i = 0u; i < iBOFFSET[iRANK]; i++) 
            {   
                uKh_BBps2[i]    = calloc(     uBPNum, sizeof *uKh_BBps2[i] );
                fKh_BBvalstk[i] = calloc(   3*uBPNum, sizeof *fKh_BBvalstk[i] );
                fKh_BBvaldip[i] = calloc(   3*uBPNum, sizeof *fKh_BBvaldip[i] );
                fKh_BBvalnrm[i] = calloc(   3*uBPNum, sizeof *fKh_BBvalnrm[i] );
                
                uGlobPos = i+ iBSTART[iRANK];
                memcpy(fRcv, &fB_temp[uGlobPos*16 +0], 3u*sizeof(float));      
                memcpy(fRotMat, &fB_temp[uGlobPos*16 +6], 9u*sizeof(float)); 
                
                uKh_BBcnt[i]  = 0u;
                uUsdNum       = 0u;
                memset(uBUsdEl, 0,         uBPNum*sizeof(unsigned int));
                memset(uBUsdBrch, 0, uB_TotalOffs*sizeof(unsigned int));
                uTstBrId = 0u;      uSrcElem = 0u;
                
                while (uUsdNum < uBPNum) 
                {   
                    for (j = 0u; j < uBPNum; j++)               {       if (uBUsdEl[j] == 0u)    {      uSrcElem = j;       break;          }           } 
                    
                    memcpy(fSrc, &fB_temp[uSrcElem*16 +0], 3u*sizeof(float));
                    
                    SubtractVect3(fTempVect, fRcv, fSrc);               fDist  = VectorLength(fTempVect);
                    
                    fDistF = (fDist >= CUTDISTANCE*fBndSide)*1.0f + (fDist < CUTDISTANCE*fBndSide)*0.0f; 
                    fDistF = (uB_temp[uGlobPos*4 +3] == uB_temp[uSrcElem*4 +3])*1.0f + (uB_temp[uGlobPos*4 +3] != uB_temp[uSrcElem*4 +3])*fDistF;
                    
                    uTemp0 = uB_OTElem[uSrcElem*u_MaxBrLvl +0];
                    uTemp1 = 0u;
                    
                    for (j = 1u; j < uTemp0; j++)
                    {   uPrvBrLv = j;
                        uPrvBrId = uB_OTElem[uSrcElem*u_MaxBrLvl +(u_MaxBrLvl -uTemp0 +uPrvBrLv)];
                        
                        memcpy(fSrc, &fKh_BBBrch[uPrvBrId*21 +0], 3u*sizeof(float));
                        SubtractVect3(fTempVect, fRcv, fSrc);               fDist  = VectorLength(fTempVect);
                        
                        uTemp1  = (uBUsdBrch[uPrvBrId]                                   == 0u)*1u     + (uBUsdBrch[uPrvBrId]                                      != 0u)*0u;
                        uTemp1  = (uB_OTPrtChld[uPrvBrId*5 +0]                            > 1u)*uTemp1 + (uB_OTPrtChld[uPrvBrId*5 +0]                              <= 1u)*0u;
                        uTemp1  = (fDist > (fKh_BBBrch[uPrvBrId*21 +14]+KH_FULL_DIST*fBndSide))*uTemp1 + (fDist <= (fKh_BBBrch[uPrvBrId*21 +14]+KH_FULL_DIST*fBndSide))*0u;
                        
                        if (uTemp1 == 1u)         {       break;      }
                    } 
                    
                    if (uTemp1 == 0u) 
                    {   memcpy(fP1, &fBV_temp[uB_temp[uSrcElem*4 +0]*3 +0], 3u*sizeof(float));
                        memcpy(fP2, &fBV_temp[uB_temp[uSrcElem*4 +1]*3 +0], 3u*sizeof(float));
                        memcpy(fP3, &fBV_temp[uB_temp[uSrcElem*4 +2]*3 +0], 3u*sizeof(float));
                        
                        if (USEHALFSPACE == 1u)
                        {   StrainHS_Nikkhoo(fStressStk, fStrain, fRcv[0], fRcv[1], fRcv[2], fP1, fP2, fP3, fUnitSlipB, 0.0f, 0.0f, fModPara[2], fModPara[4]);
                            StrainHS_Nikkhoo(fStressDip, fStrain, fRcv[0], fRcv[1], fRcv[2], fP1, fP2, fP3, 0.0f, fUnitSlipB, 0.0f, fModPara[2], fModPara[4]);
                            StrainHS_Nikkhoo(fStressNrm, fStrain, fRcv[0], fRcv[1], fRcv[2], fP1, fP2, fP3, 0.0f, 0.0f, fUnitSlipB, fModPara[2], fModPara[4]);
                        }
                        else
                        {   StrainFS_Nikkhoo(fStressStk, fStrain, fRcv[0], fRcv[1], fRcv[2], fP1, fP2, fP3, fUnitSlipB, 0.0f, 0.0f, fModPara[2], fModPara[4]);
                            StrainFS_Nikkhoo(fStressDip, fStrain, fRcv[0], fRcv[1], fRcv[2], fP1, fP2, fP3, 0.0f, fUnitSlipB, 0.0f, fModPara[2], fModPara[4]);
                            StrainFS_Nikkhoo(fStressNrm, fStrain, fRcv[0], fRcv[1], fRcv[2], fP1, fP2, fP3, 0.0f, 0.0f, fUnitSlipB, fModPara[2], fModPara[4]);
                        }
                        
                        ScaleStress6(fStressStk, (fDistF/fUnitSlipB));      RotStressG2L(&fPrvStress[0], fRotMat, fStressStk);
                        ScaleStress6(fStressDip, (fDistF/fUnitSlipB));      RotStressG2L(&fPrvStress[3], fRotMat, fStressDip);
                        ScaleStress6(fStressNrm, (fDistF/fUnitSlipB));      RotStressG2L(&fPrvStress[6], fRotMat, fStressNrm);
                        
                        if (uGlobPos == uSrcElem)
                        {   
                            fBEvent[i*9 +4] = fPrvStress[0];
                            fBEvent[i*9 +5] = fPrvStress[4];
                            fBEvent[i*9 +6] = fPrvStress[8];
                        }
                        
                        fKh_BBvalstk[i][uKh_BBcnt[i]*3 +0] = fPrvStress[0];         fKh_BBvalstk[i][uKh_BBcnt[i]*3 +1] = fPrvStress[3];         fKh_BBvalstk[i][uKh_BBcnt[i]*3 +2] = fPrvStress[6];
                        fKh_BBvaldip[i][uKh_BBcnt[i]*3 +0] = fPrvStress[1];         fKh_BBvaldip[i][uKh_BBcnt[i]*3 +1] = fPrvStress[4];         fKh_BBvaldip[i][uKh_BBcnt[i]*3 +2] = fPrvStress[7];
                        fKh_BBvalnrm[i][uKh_BBcnt[i]*3 +0] = fPrvStress[2];         fKh_BBvalnrm[i][uKh_BBcnt[i]*3 +1] = fPrvStress[5];         fKh_BBvalnrm[i][uKh_BBcnt[i]*3 +2] = fPrvStress[8];
                        
                        uBUsdEl[uSrcElem]       = 1u;
                        uKh_BBps2[i][uSrcElem]  = uKh_BBcnt[i];
                        uKh_BBcnt[i]           += 1u;
                        
                        uUsdNum = 0u;
                        for (j = 0u; j < uBPNum;  j++)      {       uUsdNum +=  uBUsdEl[j];                                     }
                        
                        for (j = 1u; j < uTemp0; j++)       {       uBUsdBrch[ uB_OTElem[uSrcElem*u_MaxBrLvl +(u_MaxBrLvl -uTemp0 +j)] ] = 1u;       }
                    }
                    else
                    {   uPrvBrLv -= 1u;
                        uTemp2    = 1u;
                        while (uTemp2 == 1u)
                        {   
                            uPrvBrId = uB_OTElem[uSrcElem*u_MaxBrLvl +(u_MaxBrLvl -uTemp0 +uPrvBrLv)];
                            uTemp1   = 0;
                            
                            while (uPrvBrLv < (uTemp0 -2))
                            {   uPrvBrLv += 1u;
                                uPrvBrId = uB_OTElem[uSrcElem*u_MaxBrLvl +(u_MaxBrLvl -uTemp0 +uPrvBrLv)];
                                
                                memcpy(fSrc, &fKh_BBBrch[uPrvBrId*21 +0], 3u*sizeof(float));
                                SubtractVect3(fTempVect, fRcv, fSrc);               fDist  = VectorLength(fTempVect);
                                
                                uTemp1  = (uBUsdBrch[uPrvBrId]                                   == 0u)*1u     + (uBUsdBrch[uPrvBrId]                                      != 0u)*0u;
                                uTemp1  = (uB_OTPrtChld[uPrvBrId*5 +0]                            > 1u)*uTemp1 + (uB_OTPrtChld[uPrvBrId*5 +0]                              <= 1u)*0u;
                                uTemp1  = (fDist > (fKh_BBBrch[uPrvBrId*21 +14]+KH_FULL_DIST*fBndSide))*uTemp1 + (fDist <= (fKh_BBBrch[uPrvBrId*21 +14]+KH_FULL_DIST*fBndSide))*0u;
                                
                                if (uTemp1 == 1u)         {       break;      }
                            }
                            
                            if (uTemp1 == 0u)
                            {   memcpy(fP1, &fBV_temp[uB_temp[uSrcElem*4 +0]*3 +0], 3u*sizeof(float));
                                memcpy(fP2, &fBV_temp[uB_temp[uSrcElem*4 +1]*3 +0], 3u*sizeof(float));
                                memcpy(fP3, &fBV_temp[uB_temp[uSrcElem*4 +2]*3 +0], 3u*sizeof(float));
                                
                                if (USEHALFSPACE == 1u)
                                {   StrainHS_Nikkhoo(fStressStk, fStrain, fRcv[0], fRcv[1], fRcv[2], fP1, fP2, fP3, fUnitSlipB, 0.0f, 0.0f, fModPara[2], fModPara[4]);
                                    StrainHS_Nikkhoo(fStressDip, fStrain, fRcv[0], fRcv[1], fRcv[2], fP1, fP2, fP3, 0.0f, fUnitSlipB, 0.0f, fModPara[2], fModPara[4]);
                                    StrainHS_Nikkhoo(fStressNrm, fStrain, fRcv[0], fRcv[1], fRcv[2], fP1, fP2, fP3, 0.0f, 0.0f, fUnitSlipB, fModPara[2], fModPara[4]);
                                }
                                else
                                {   StrainFS_Nikkhoo(fStressStk, fStrain, fRcv[0], fRcv[1], fRcv[2], fP1, fP2, fP3, fUnitSlipB, 0.0f, 0.0f, fModPara[2], fModPara[4]);
                                    StrainFS_Nikkhoo(fStressDip, fStrain, fRcv[0], fRcv[1], fRcv[2], fP1, fP2, fP3, 0.0f, fUnitSlipB, 0.0f, fModPara[2], fModPara[4]);
                                    StrainFS_Nikkhoo(fStressNrm, fStrain, fRcv[0], fRcv[1], fRcv[2], fP1, fP2, fP3, 0.0f, 0.0f, fUnitSlipB, fModPara[2], fModPara[4]);
                                }
                                
                                ScaleStress6(fStressStk, (fDistF/fUnitSlipB));          RotStressG2L(&fPrvStress[0], fRotMat, fStressStk);
                                ScaleStress6(fStressDip, (fDistF/fUnitSlipB));          RotStressG2L(&fPrvStress[3], fRotMat, fStressDip);
                                ScaleStress6(fStressNrm, (fDistF/fUnitSlipB));          RotStressG2L(&fPrvStress[6], fRotMat, fStressNrm);
                                
                                fKh_BBvalstk[i][uKh_BBcnt[i]*3 +0] = fPrvStress[0];         fKh_BBvalstk[i][uKh_BBcnt[i]*3 +1] = fPrvStress[3];         fKh_BBvalstk[i][uKh_BBcnt[i]*3 +2] = fPrvStress[6];
                                fKh_BBvaldip[i][uKh_BBcnt[i]*3 +0] = fPrvStress[1];         fKh_BBvaldip[i][uKh_BBcnt[i]*3 +1] = fPrvStress[4];         fKh_BBvaldip[i][uKh_BBcnt[i]*3 +2] = fPrvStress[7];
                                fKh_BBvalnrm[i][uKh_BBcnt[i]*3 +0] = fPrvStress[2];         fKh_BBvalnrm[i][uKh_BBcnt[i]*3 +1] = fPrvStress[5];         fKh_BBvalnrm[i][uKh_BBcnt[i]*3 +2] = fPrvStress[8];
                                
                                uBUsdEl[uSrcElem]       = 1u;
                                uKh_BBps2[i][uSrcElem]  = uKh_BBcnt[i];
                                uKh_BBcnt[i]           += 1u;
                                
                                uUsdNum = 0u;
                                for (j = 0u; j < uBPNum;  j++)      {       uUsdNum +=  uBUsdEl[j];         }
                                
                                for (j = 1u; j < uTemp0; j++)       {       uBUsdBrch[ uB_OTElem[uSrcElem*u_MaxBrLvl +(u_MaxBrLvl -uTemp0 +j)] ] = 1u;       }
                                uTemp2 = 0u;
                            }
                            else
                            {   
                                GetGlobVertsForRectangle(fP1, fP2, fP3, fP4, fKh_BBBrch, uPrvBrId*21u);
                                
                                if (USEHALFSPACE == 1u)
                                {   StrainHS_Nikkhoo(fStressStkT, fStrain, fRcv[0], fRcv[1], fRcv[2], fP1, fP2, fP4, fUnitSlipB, 0.0f, 0.0f, fModPara[2], fModPara[4]);
                                    StrainHS_Nikkhoo(fStressDipT, fStrain, fRcv[0], fRcv[1], fRcv[2], fP1, fP2, fP4, 0.0f, fUnitSlipB, 0.0f, fModPara[2], fModPara[4]);
                                    StrainHS_Nikkhoo(fStressNrmT, fStrain, fRcv[0], fRcv[1], fRcv[2], fP1, fP2, fP4, 0.0f, 0.0f, fUnitSlipB, fModPara[2], fModPara[4]);
                                    
                                    StrainHS_Nikkhoo(fStressStk, fStrain, fRcv[0], fRcv[1], fRcv[2], fP2, fP3, fP4, fUnitSlipB, 0.0f, 0.0f, fModPara[2], fModPara[4]);
                                    StrainHS_Nikkhoo(fStressDip, fStrain, fRcv[0], fRcv[1], fRcv[2], fP2, fP3, fP4, 0.0f, fUnitSlipB, 0.0f, fModPara[2], fModPara[4]);
                                    StrainHS_Nikkhoo(fStressNrm, fStrain, fRcv[0], fRcv[1], fRcv[2], fP2, fP3, fP4, 0.0f, 0.0f, fUnitSlipB, fModPara[2], fModPara[4]);
                                    
                                    fStressStk[0] += fStressStkT[0];        fStressStk[1] += fStressStkT[1];        fStressStk[2] += fStressStkT[2];
                                    fStressStk[3] += fStressStkT[3];        fStressStk[4] += fStressStkT[4];        fStressStk[5] += fStressStkT[5];
                                    fStressDip[0] += fStressDipT[0];        fStressDip[1] += fStressDipT[1];        fStressDip[2] += fStressDipT[2];
                                    fStressDip[3] += fStressDipT[3];        fStressDip[4] += fStressDipT[4];        fStressDip[5] += fStressDipT[5];
                                    fStressNrm[0] += fStressNrmT[0];        fStressNrm[1] += fStressNrmT[1];        fStressNrm[2] += fStressNrmT[2];
                                    fStressNrm[3] += fStressNrmT[3];        fStressNrm[4] += fStressNrmT[4];        fStressNrm[5] += fStressNrmT[5];
                                }
                                else
                                {   StrainFS_Nikkhoo(fStressStkT, fStrain, fRcv[0], fRcv[1], fRcv[2], fP1, fP2, fP4, fUnitSlipB, 0.0f, 0.0f, fModPara[2], fModPara[4]);
                                    StrainFS_Nikkhoo(fStressDipT, fStrain, fRcv[0], fRcv[1], fRcv[2], fP1, fP2, fP4, 0.0f, fUnitSlipB, 0.0f, fModPara[2], fModPara[4]);
                                    StrainFS_Nikkhoo(fStressNrmT, fStrain, fRcv[0], fRcv[1], fRcv[2], fP1, fP2, fP4, 0.0f, 0.0f, fUnitSlipB, fModPara[2], fModPara[4]);
                                    
                                    StrainFS_Nikkhoo(fStressStk, fStrain, fRcv[0], fRcv[1], fRcv[2], fP2, fP3, fP4, fUnitSlipB, 0.0f, 0.0f, fModPara[2], fModPara[4]);
                                    StrainFS_Nikkhoo(fStressDip, fStrain, fRcv[0], fRcv[1], fRcv[2], fP2, fP3, fP4, 0.0f, fUnitSlipB, 0.0f, fModPara[2], fModPara[4]);
                                    StrainFS_Nikkhoo(fStressNrm, fStrain, fRcv[0], fRcv[1], fRcv[2], fP2, fP3, fP4, 0.0f, 0.0f, fUnitSlipB, fModPara[2], fModPara[4]);
                                    
                                    fStressStk[0] += fStressStkT[0];        fStressStk[1] += fStressStkT[1];        fStressStk[2] += fStressStkT[2];
                                    fStressStk[3] += fStressStkT[3];        fStressStk[4] += fStressStkT[4];        fStressStk[5] += fStressStkT[5];
                                    fStressDip[0] += fStressDipT[0];        fStressDip[1] += fStressDipT[1];        fStressDip[2] += fStressDipT[2];
                                    fStressDip[3] += fStressDipT[3];        fStressDip[4] += fStressDipT[4];        fStressDip[5] += fStressDipT[5];
                                    fStressNrm[0] += fStressNrmT[0];        fStressNrm[1] += fStressNrmT[1];        fStressNrm[2] += fStressNrmT[2];
                                    fStressNrm[3] += fStressNrmT[3];        fStressNrm[4] += fStressNrmT[4];        fStressNrm[5] += fStressNrmT[5];
                                }
                                
                                ScaleStress6(fStressStk, (fDistF/fUnitSlipB));          ScaleStress6(fStressStk, (1.0f/fKh_BBBrch[uPrvBrId*21 +13]));           RotStressG2L(&fPrvStress[0], fRotMat, fStressStk);
                                ScaleStress6(fStressDip, (fDistF/fUnitSlipB));          ScaleStress6(fStressDip, (1.0f/fKh_BBBrch[uPrvBrId*21 +13]));           RotStressG2L(&fPrvStress[3], fRotMat, fStressDip);
                                ScaleStress6(fStressNrm, (fDistF/fUnitSlipB));          ScaleStress6(fStressNrm, (1.0f/fKh_BBBrch[uPrvBrId*21 +13]));           RotStressG2L(&fPrvStress[6], fRotMat, fStressNrm);
                                
                                uTemp1 = 1u;
                                for (j = 0; j < uB_OTPrtChld[uPrvBrId*5 +0]; j++)
                                {   uTstBrId = uB_OTPrtChld[uPrvBrId*5 +(j+1)];
                                    
                                    GetGlobVertsForRectangle(fP1, fP2, fP3, fP4, fKh_BBBrch, uTstBrId*21u);
                                    
                                    if (USEHALFSPACE == 1u)
                                    {   StrainHS_Nikkhoo(fStressStkT, fStrain, fRcv[0], fRcv[1], fRcv[2], fP1, fP2, fP4, fUnitSlipB, 0.0f, 0.0f, fModPara[2], fModPara[4]);
                                        StrainHS_Nikkhoo(fStressDipT, fStrain, fRcv[0], fRcv[1], fRcv[2], fP1, fP2, fP4, 0.0f, fUnitSlipB, 0.0f, fModPara[2], fModPara[4]);
                                        StrainHS_Nikkhoo(fStressNrmT, fStrain, fRcv[0], fRcv[1], fRcv[2], fP1, fP2, fP4, 0.0f, 0.0f, fUnitSlipB, fModPara[2], fModPara[4]);
                                        
                                        StrainHS_Nikkhoo(fStressStk, fStrain, fRcv[0], fRcv[1], fRcv[2], fP2, fP3, fP4, fUnitSlipB, 0.0f, 0.0f, fModPara[2], fModPara[4]);
                                        StrainHS_Nikkhoo(fStressDip, fStrain, fRcv[0], fRcv[1], fRcv[2], fP2, fP3, fP4, 0.0f, fUnitSlipB, 0.0f, fModPara[2], fModPara[4]);
                                        StrainHS_Nikkhoo(fStressNrm, fStrain, fRcv[0], fRcv[1], fRcv[2], fP2, fP3, fP4, 0.0f, 0.0f, fUnitSlipB, fModPara[2], fModPara[4]);
                                        
                                        fStressStk[0] += fStressStkT[0];        fStressStk[1] += fStressStkT[1];        fStressStk[2] += fStressStkT[2];
                                        fStressStk[3] += fStressStkT[3];        fStressStk[4] += fStressStkT[4];        fStressStk[5] += fStressStkT[5];
                                        fStressDip[0] += fStressDipT[0];        fStressDip[1] += fStressDipT[1];        fStressDip[2] += fStressDipT[2];
                                        fStressDip[3] += fStressDipT[3];        fStressDip[4] += fStressDipT[4];        fStressDip[5] += fStressDipT[5];
                                        fStressNrm[0] += fStressNrmT[0];        fStressNrm[1] += fStressNrmT[1];        fStressNrm[2] += fStressNrmT[2];
                                        fStressNrm[3] += fStressNrmT[3];        fStressNrm[4] += fStressNrmT[4];        fStressNrm[5] += fStressNrmT[5];
                                    }
                                    else
                                    {   StrainFS_Nikkhoo(fStressStkT, fStrain, fRcv[0], fRcv[1], fRcv[2], fP1, fP2, fP4, fUnitSlipB, 0.0f, 0.0f, fModPara[2], fModPara[4]);
                                        StrainFS_Nikkhoo(fStressDipT, fStrain, fRcv[0], fRcv[1], fRcv[2], fP1, fP2, fP4, 0.0f, fUnitSlipB, 0.0f, fModPara[2], fModPara[4]);
                                        StrainFS_Nikkhoo(fStressNrmT, fStrain, fRcv[0], fRcv[1], fRcv[2], fP1, fP2, fP4, 0.0f, 0.0f, fUnitSlipB, fModPara[2], fModPara[4]);
                                        
                                        StrainFS_Nikkhoo(fStressStk, fStrain, fRcv[0], fRcv[1], fRcv[2], fP2, fP3, fP4, fUnitSlipB, 0.0f, 0.0f, fModPara[2], fModPara[4]);
                                        StrainFS_Nikkhoo(fStressDip, fStrain, fRcv[0], fRcv[1], fRcv[2], fP2, fP3, fP4, 0.0f, fUnitSlipB, 0.0f, fModPara[2], fModPara[4]);
                                        StrainFS_Nikkhoo(fStressNrm, fStrain, fRcv[0], fRcv[1], fRcv[2], fP2, fP3, fP4, 0.0f, 0.0f, fUnitSlipB, fModPara[2], fModPara[4]);
                                        
                                        fStressStk[0] += fStressStkT[0];        fStressStk[1] += fStressStkT[1];        fStressStk[2] += fStressStkT[2];
                                        fStressStk[3] += fStressStkT[3];        fStressStk[4] += fStressStkT[4];        fStressStk[5] += fStressStkT[5];
                                        fStressDip[0] += fStressDipT[0];        fStressDip[1] += fStressDipT[1];        fStressDip[2] += fStressDipT[2];
                                        fStressDip[3] += fStressDipT[3];        fStressDip[4] += fStressDipT[4];        fStressDip[5] += fStressDipT[5];
                                        fStressNrm[0] += fStressNrmT[0];        fStressNrm[1] += fStressNrmT[1];        fStressNrm[2] += fStressNrmT[2];
                                        fStressNrm[3] += fStressNrmT[3];        fStressNrm[4] += fStressNrmT[4];        fStressNrm[5] += fStressNrmT[5];
                                    }
                                    
                                    ScaleStress6(fStressStk, (fDistF/fUnitSlipB));          ScaleStress6(fStressStk, (1.0f/fKh_BBBrch[uTstBrId*21 +13]));           RotStressG2L(&fTstStress[0], fRotMat, fStressStk);
                                    ScaleStress6(fStressDip, (fDistF/fUnitSlipB));          ScaleStress6(fStressDip, (1.0f/fKh_BBBrch[uTstBrId*21 +13]));           RotStressG2L(&fTstStress[3], fRotMat, fStressDip);
                                    ScaleStress6(fStressNrm, (fDistF/fUnitSlipB));          ScaleStress6(fStressNrm, (1.0f/fKh_BBBrch[uTstBrId*21 +13]));           RotStressG2L(&fTstStress[6], fRotMat, fStressNrm);
                                    
                                    memcpy(fSrc, &fKh_BBBrch[uTstBrId*21 +0], 3u*sizeof(float));
                                    SubtractVect3(fTempVect, fRcv, fSrc);
                                    fDist  = VectorLength(fTempVect);
                                    fTemp2 = (fDist/fBndSide - KH_FULL_DIST)/(KH_TOL_DIST - KH_FULL_DIST);
                                    fTemp2 = (fTemp2 >= 0.0f)*fTemp2 + (fTemp2 < 0.0f)*0.0f;
                                    fTemp2 = (fTemp2 <= 1.0f)*fTemp2 + (fTemp2 > 1.0f)*1.0f;
                                    fTemp2 = KH_TOL1 + (KH_TOL2-KH_TOL1)*fTemp2;
                                    
                                    fTemp0 = sqrtf( (fTstStress[0]-fPrvStress[0])*(fTstStress[0]-fPrvStress[0]) + (fTstStress[1]-fPrvStress[1])*(fTstStress[1]-fPrvStress[1]) + (fTstStress[2]-fPrvStress[2])*(fTstStress[2]-fPrvStress[2])  + (fTstStress[3]-fPrvStress[3])*(fTstStress[3]-fPrvStress[3]) + (fTstStress[4]-fPrvStress[4])*(fTstStress[4]-fPrvStress[4]) + (fTstStress[5]-fPrvStress[5])*(fTstStress[5]-fPrvStress[5])   + (fTstStress[6]-fPrvStress[6])*(fTstStress[6]-fPrvStress[6]) + (fTstStress[7]-fPrvStress[7])*(fTstStress[7]-fPrvStress[7]) + (fTstStress[8]-fPrvStress[8])*(fTstStress[8]-fPrvStress[8]));
                                    fTemp1 = sqrtf(                fPrvStress[0] *fPrvStress[0]                 +                fPrvStress[1] *fPrvStress[1]                 +                fPrvStress[2] *fPrvStress[2]                  +                fPrvStress[3] *fPrvStress[3]                 +                fPrvStress[4] *fPrvStress[4]                 +                fPrvStress[5] *fPrvStress[5]                   +                fPrvStress[6] *fPrvStress[6]                 +                fPrvStress[7] *fPrvStress[7]                 +                fPrvStress[8] *fPrvStress[8]);
                                    
                                    uTemp1 = ((fTemp0/fTemp1) <= fTemp2)*uTemp1 + ((fTemp0/fTemp1) > fTemp2)*0u; 
                                }
                                
                                if (uTemp1 == 1u)
                                {   fKh_BBvalstk[i][uKh_BBcnt[i]*3 +0] = fPrvStress[0];         fKh_BBvalstk[i][uKh_BBcnt[i]*3 +1] = fPrvStress[3];         fKh_BBvalstk[i][uKh_BBcnt[i]*3 +2] = fPrvStress[6];
                                    fKh_BBvaldip[i][uKh_BBcnt[i]*3 +0] = fPrvStress[1];         fKh_BBvaldip[i][uKh_BBcnt[i]*3 +1] = fPrvStress[4];         fKh_BBvaldip[i][uKh_BBcnt[i]*3 +2] = fPrvStress[7];
                                    fKh_BBvalnrm[i][uKh_BBcnt[i]*3 +0] = fPrvStress[2];         fKh_BBvalnrm[i][uKh_BBcnt[i]*3 +1] = fPrvStress[5];         fKh_BBvalnrm[i][uKh_BBcnt[i]*3 +2] = fPrvStress[8];
                                    
                                    uUsdNum = 0u;
                                    for (j = 0u; j < uBPNum;  j++)
                                    {   uBUsdEl[j]      = (uB_OTElem[j*u_MaxBrLvl +(u_MaxBrLvl -uTemp0 +uPrvBrLv)] == uPrvBrId)*1u           + (uB_OTElem[j*u_MaxBrLvl +(u_MaxBrLvl -uTemp0 +uPrvBrLv)] != uPrvBrId)*uBUsdEl[j];
                                        uKh_BBps2[i][j] = (uB_OTElem[j*u_MaxBrLvl +(u_MaxBrLvl -uTemp0 +uPrvBrLv)] == uPrvBrId)*uKh_BBcnt[i] + (uB_OTElem[j*u_MaxBrLvl +(u_MaxBrLvl -uTemp0 +uPrvBrLv)] != uPrvBrId)*uKh_BBps2[i][j];
                                        uUsdNum        +=  uBUsdEl[j];
                                    }
                                    uKh_BBcnt[i] += 1u;
                                    
                                    for (j = 1u; j < uTemp0; j++)       {       uBUsdBrch[ uB_OTElem[uSrcElem*u_MaxBrLvl +(u_MaxBrLvl -uTemp0 +j)] ] = 1u;       }
                                    uTemp2 = 0u;
                }   }   }   }   }
                
                fKh_BBvalstk[i] = realloc(fKh_BBvalstk[i], 3*uKh_BBcnt[i]* sizeof *fKh_BBvalstk[i] );
                fKh_BBvaldip[i] = realloc(fKh_BBvaldip[i], 3*uKh_BBcnt[i]* sizeof *fKh_BBvaldip[i] );
                fKh_BBvalnrm[i] = realloc(fKh_BBvalnrm[i], 3*uKh_BBcnt[i]* sizeof *fKh_BBvalnrm[i] );
                
                if ((uPlotCatalog2Screen == 1)&&(iRANK == 0)) {    fprintf(stdout,"\nboundary receiver:  %u/%u   with Bnd/Bnd:  %u of %u",i, iBOFFSET[iRANK], uKh_BBcnt[i], uBPNum);                  }
                
                
                uKh_BFps2[i]    = calloc(     uFPNum, sizeof *uKh_BFps2[i] );
                fKh_BFvalstk[i] = calloc(   2*uFPNum, sizeof *fKh_BFvalstk[i] );
                fKh_BFvaldip[i] = calloc(   2*uFPNum, sizeof *fKh_BFvaldip[i] );
                fKh_BFvalnrm[i] = calloc(   2*uFPNum, sizeof *fKh_BFvalnrm[i] );
                
                uKh_BFcnt[i]  = 0u;
                uUsdNum       = 0u;
                memset(uFUsdEl, 0,         uFPNum*sizeof(unsigned int));
                memset(uFUsdBrch, 0, uF_TotalOffs*sizeof(unsigned int));
                uTstBrId = 0u;      uSrcElem = 0u;
                
                while (uUsdNum < uFPNum) 
                {   
                    for (j = 0u; j < uFPNum; j++)               {       if (uFUsdEl[j] == 0u)    {      uSrcElem = j;       break;          }           } 
                    
                    memcpy(fSrc, &fF_temp[uSrcElem*16 +0], 3u*sizeof(float));
                    
                    SubtractVect3(fTempVect, fRcv, fSrc);               fDist  = VectorLength(fTempVect);
                    
                    fDistF = (fDist >= CUTDISTANCE*fFltSide)*1.0f + (fDist < CUTDISTANCE*fFltSide)*0.0f; 
                    
                    uTemp0 = uF_OTElem[uSrcElem*u_MaxBrLvl +0];
                    uTemp1 = 0u;
                    
                    for (j = 1u; j < uTemp0; j++)
                    {   uPrvBrLv = j;
                        uPrvBrId = uF_OTElem[uSrcElem*u_MaxBrLvl +(u_MaxBrLvl -uTemp0 +uPrvBrLv)];
                        
                        memcpy(fSrc, &fKh_FFBrch[uPrvBrId*21 +0], 3u*sizeof(float));
                        SubtractVect3(fTempVect, fRcv, fSrc);               fDist  = VectorLength(fTempVect);
                        
                        uTemp1  = (uFUsdBrch[uPrvBrId]                                   == 0u)*1u     + (uFUsdBrch[uPrvBrId]                                      != 0u)*0u;
                        uTemp1  = (uF_OTPrtChld[uPrvBrId*5 +0]                            > 1u)*uTemp1 + (uF_OTPrtChld[uPrvBrId*5 +0]                              <= 1u)*0u;
                        uTemp1  = (fDist > (fKh_FFBrch[uPrvBrId*21 +14]+KH_FULL_DIST*fFltSide))*uTemp1 + (fDist <= (fKh_FFBrch[uPrvBrId*21 +14]+KH_FULL_DIST*fFltSide))*0u;
                        
                        if (uTemp1 == 1u)         {       break;      }
                    } 
                    
                    if (uTemp1 == 0u) 
                    {   memcpy(fP1, &fFV_temp[uF_temp[uSrcElem*7 +0]*4 +0], 3u*sizeof(float));
                        memcpy(fP2, &fFV_temp[uF_temp[uSrcElem*7 +1]*4 +0], 3u*sizeof(float));
                        memcpy(fP3, &fFV_temp[uF_temp[uSrcElem*7 +2]*4 +0], 3u*sizeof(float));
                        
                        if (USEHALFSPACE == 1u)
                        {   StrainHS_Nikkhoo(fStressStk, fStrain, fRcv[0], fRcv[1], fRcv[2], fP1, fP2, fP3, fUnitSlipF, 0.0f, 0.0f, fModPara[2], fModPara[4]);
                            StrainHS_Nikkhoo(fStressDip, fStrain, fRcv[0], fRcv[1], fRcv[2], fP1, fP2, fP3, 0.0f, fUnitSlipF, 0.0f, fModPara[2], fModPara[4]);
                            StrainHS_Nikkhoo(fStressNrm, fStrain, fRcv[0], fRcv[1], fRcv[2], fP1, fP2, fP3, 0.0f, 0.0f, fUnitSlipF, fModPara[2], fModPara[4]);
                        }
                        else
                        {   StrainFS_Nikkhoo(fStressStk, fStrain, fRcv[0], fRcv[1], fRcv[2], fP1, fP2, fP3, fUnitSlipF, 0.0f, 0.0f, fModPara[2], fModPara[4]);
                            StrainFS_Nikkhoo(fStressDip, fStrain, fRcv[0], fRcv[1], fRcv[2], fP1, fP2, fP3, 0.0f, fUnitSlipF, 0.0f, fModPara[2], fModPara[4]);
                            StrainFS_Nikkhoo(fStressNrm, fStrain, fRcv[0], fRcv[1], fRcv[2], fP1, fP2, fP3, 0.0f, 0.0f, fUnitSlipF, fModPara[2], fModPara[4]);
                        }
                        
                        ScaleStress6(fStressStk, (fDistF/fUnitSlipF));              RotStressG2L(&fPrvStress[0], fRotMat, fStressStk);
                        ScaleStress6(fStressDip, (fDistF/fUnitSlipF));              RotStressG2L(&fPrvStress[3], fRotMat, fStressDip);
                        ScaleStress6(fStressNrm, (fDistF/fUnitSlipF));              RotStressG2L(&fPrvStress[6], fRotMat, fStressNrm);
                        
                        fKh_BFvalstk[i][uKh_BFcnt[i]*2 +0] = fPrvStress[0];         fKh_BFvalstk[i][uKh_BFcnt[i]*2 +1] = fPrvStress[3];
                        fKh_BFvaldip[i][uKh_BFcnt[i]*2 +0] = fPrvStress[1];         fKh_BFvaldip[i][uKh_BFcnt[i]*2 +1] = fPrvStress[4];
                        fKh_BFvalnrm[i][uKh_BFcnt[i]*2 +0] = fPrvStress[2];         fKh_BFvalnrm[i][uKh_BFcnt[i]*2 +1] = fPrvStress[5];
                        
                        uFUsdEl[uSrcElem]       = 1u;
                        uKh_BFps2[i][uSrcElem]  = uKh_BFcnt[i];
                        uKh_BFcnt[i]           += 1u;
                        
                        uUsdNum = 0u;
                        for (j = 0u; j < uFPNum;  j++)      {       uUsdNum +=  uFUsdEl[j];                                     }
                        
                        for (j = 1u; j < uTemp0; j++)       {       uFUsdBrch[ uF_OTElem[uSrcElem*u_MaxBrLvl +(u_MaxBrLvl -uTemp0 +j)] ] = 1u;       }
                    }
                    else
                    {   uPrvBrLv -= 1u;
                        uTemp2    = 1u;
                        while (uTemp2 == 1u)
                        {   
                            uPrvBrId = uF_OTElem[uSrcElem*u_MaxBrLvl +(u_MaxBrLvl -uTemp0 +uPrvBrLv)];
                            uTemp1   = 0u;
                            
                            while (uPrvBrLv < (uTemp0 -2))
                            {   uPrvBrLv += 1u;
                                uPrvBrId = uF_OTElem[uSrcElem*u_MaxBrLvl +(u_MaxBrLvl -uTemp0 +uPrvBrLv)];
                                
                                memcpy(fSrc, &fKh_FFBrch[uPrvBrId*21 +0], 3u*sizeof(float));
                                SubtractVect3(fTempVect, fRcv, fSrc);               fDist  = VectorLength(fTempVect);
                                
                                uTemp1  = (uFUsdBrch[uPrvBrId]                                   == 0u)*1u     + (uFUsdBrch[uPrvBrId]                                      != 0u)*0u;
                                uTemp1  = (uF_OTPrtChld[uPrvBrId*5 +0]                            > 1u)*uTemp1 + (uF_OTPrtChld[uPrvBrId*5 +0]                              <= 1u)*0u;
                                uTemp1  = (fDist > (fKh_FFBrch[uPrvBrId*21 +14]+KH_FULL_DIST*fFltSide))*uTemp1 + (fDist <= (fKh_FFBrch[uPrvBrId*21 +14]+KH_FULL_DIST*fFltSide))*0u;
                                
                                if (uTemp1 == 1u)         {       break;      }
                            }
                            
                            if (uTemp1 == 0u)
                            {   memcpy(fP1, &fFV_temp[uF_temp[uSrcElem*7 +0]*4 +0], 3u*sizeof(float));
                                memcpy(fP2, &fFV_temp[uF_temp[uSrcElem*7 +1]*4 +0], 3u*sizeof(float));
                                memcpy(fP3, &fFV_temp[uF_temp[uSrcElem*7 +2]*4 +0], 3u*sizeof(float));
                                
                                if (USEHALFSPACE == 1u)
                                {   StrainHS_Nikkhoo(fStressStk, fStrain, fRcv[0], fRcv[1], fRcv[2], fP1, fP2, fP3, fUnitSlipF, 0.0f, 0.0f, fModPara[2], fModPara[4]);
                                    StrainHS_Nikkhoo(fStressDip, fStrain, fRcv[0], fRcv[1], fRcv[2], fP1, fP2, fP3, 0.0f, fUnitSlipF, 0.0f, fModPara[2], fModPara[4]);
                                    StrainHS_Nikkhoo(fStressNrm, fStrain, fRcv[0], fRcv[1], fRcv[2], fP1, fP2, fP3, 0.0f, 0.0f, fUnitSlipF, fModPara[2], fModPara[4]);
                                }
                                else
                                {   StrainFS_Nikkhoo(fStressStk, fStrain, fRcv[0], fRcv[1], fRcv[2], fP1, fP2, fP3, fUnitSlipF, 0.0f, 0.0f, fModPara[2], fModPara[4]);
                                    StrainFS_Nikkhoo(fStressDip, fStrain, fRcv[0], fRcv[1], fRcv[2], fP1, fP2, fP3, 0.0f, fUnitSlipF, 0.0f, fModPara[2], fModPara[4]);
                                    StrainFS_Nikkhoo(fStressNrm, fStrain, fRcv[0], fRcv[1], fRcv[2], fP1, fP2, fP3, 0.0f, 0.0f, fUnitSlipF, fModPara[2], fModPara[4]);
                                }
                                
                                ScaleStress6(fStressStk, (fDistF/fUnitSlipF));          RotStressG2L(&fPrvStress[0], fRotMat, fStressStk);
                                ScaleStress6(fStressDip, (fDistF/fUnitSlipF));          RotStressG2L(&fPrvStress[3], fRotMat, fStressDip);
                                ScaleStress6(fStressNrm, (fDistF/fUnitSlipF));          RotStressG2L(&fPrvStress[6], fRotMat, fStressNrm);
                                
                                fKh_BFvalstk[i][uKh_BFcnt[i]*2 +0] = fPrvStress[0];         fKh_BFvalstk[i][uKh_BFcnt[i]*2 +1] = fPrvStress[3];
                                fKh_BFvaldip[i][uKh_BFcnt[i]*2 +0] = fPrvStress[1];         fKh_BFvaldip[i][uKh_BFcnt[i]*2 +1] = fPrvStress[4];
                                fKh_BFvalnrm[i][uKh_BFcnt[i]*2 +0] = fPrvStress[2];         fKh_BFvalnrm[i][uKh_BFcnt[i]*2 +1] = fPrvStress[5];
                                
                                uFUsdEl[uSrcElem]       = 1u;
                                uKh_BFps2[i][uSrcElem]  = uKh_BFcnt[i];
                                uKh_BFcnt[i]           += 1u;
                                
                                uUsdNum = 0u;
                                for (j = 0u; j < uFPNum;  j++)      {       uUsdNum +=  uFUsdEl[j];         }
                                
                                for (j = 1u; j < uTemp0; j++)       {       uFUsdBrch[ uF_OTElem[uSrcElem*u_MaxBrLvl +(u_MaxBrLvl -uTemp0 +j)] ] = 1u;       }
                                uTemp2 = 0u;
                            }
                            else
                            {   
                                GetGlobVertsForRectangle(fP1, fP2, fP3, fP4, fKh_FFBrch, uPrvBrId*21u);
                                
                                if (USEHALFSPACE == 1u)
                                {   StrainHS_Nikkhoo(fStressStkT, fStrain, fRcv[0], fRcv[1], fRcv[2], fP1, fP2, fP4, fUnitSlipF, 0.0f, 0.0f, fModPara[2], fModPara[4]);
                                    StrainHS_Nikkhoo(fStressDipT, fStrain, fRcv[0], fRcv[1], fRcv[2], fP1, fP2, fP4, 0.0f, fUnitSlipF, 0.0f, fModPara[2], fModPara[4]);
                                    StrainHS_Nikkhoo(fStressNrmT, fStrain, fRcv[0], fRcv[1], fRcv[2], fP1, fP2, fP4, 0.0f, 0.0f, fUnitSlipF, fModPara[2], fModPara[4]);
                                    
                                    StrainHS_Nikkhoo(fStressStk, fStrain, fRcv[0], fRcv[1], fRcv[2], fP2, fP3, fP4, fUnitSlipF, 0.0f, 0.0f, fModPara[2], fModPara[4]);
                                    StrainHS_Nikkhoo(fStressDip, fStrain, fRcv[0], fRcv[1], fRcv[2], fP2, fP3, fP4, 0.0f, fUnitSlipF, 0.0f, fModPara[2], fModPara[4]);
                                    StrainHS_Nikkhoo(fStressNrm, fStrain, fRcv[0], fRcv[1], fRcv[2], fP2, fP3, fP4, 0.0f, 0.0f, fUnitSlipF, fModPara[2], fModPara[4]);
                                    
                                    fStressStk[0] += fStressStkT[0];        fStressStk[1] += fStressStkT[1];        fStressStk[2] += fStressStkT[2];
                                    fStressStk[3] += fStressStkT[3];        fStressStk[4] += fStressStkT[4];        fStressStk[5] += fStressStkT[5];
                                    fStressDip[0] += fStressDipT[0];        fStressDip[1] += fStressDipT[1];        fStressDip[2] += fStressDipT[2];
                                    fStressDip[3] += fStressDipT[3];        fStressDip[4] += fStressDipT[4];        fStressDip[5] += fStressDipT[5];
                                    fStressNrm[0] += fStressNrmT[0];        fStressNrm[1] += fStressNrmT[1];        fStressNrm[2] += fStressNrmT[2];
                                    fStressNrm[3] += fStressNrmT[3];        fStressNrm[4] += fStressNrmT[4];        fStressNrm[5] += fStressNrmT[5];
                                }
                                else
                                {   StrainFS_Nikkhoo(fStressStkT, fStrain, fRcv[0], fRcv[1], fRcv[2], fP1, fP2, fP4, fUnitSlipF, 0.0f, 0.0f, fModPara[2], fModPara[4]);
                                    StrainFS_Nikkhoo(fStressDipT, fStrain, fRcv[0], fRcv[1], fRcv[2], fP1, fP2, fP4, 0.0f, fUnitSlipF, 0.0f, fModPara[2], fModPara[4]);
                                    StrainFS_Nikkhoo(fStressNrmT, fStrain, fRcv[0], fRcv[1], fRcv[2], fP1, fP2, fP4, 0.0f, 0.0f, fUnitSlipF, fModPara[2], fModPara[4]);
                                    
                                    StrainFS_Nikkhoo(fStressStk, fStrain, fRcv[0], fRcv[1], fRcv[2], fP2, fP3, fP4, fUnitSlipF, 0.0f, 0.0f, fModPara[2], fModPara[4]);
                                    StrainFS_Nikkhoo(fStressDip, fStrain, fRcv[0], fRcv[1], fRcv[2], fP2, fP3, fP4, 0.0f, fUnitSlipF, 0.0f, fModPara[2], fModPara[4]);
                                    StrainFS_Nikkhoo(fStressNrm, fStrain, fRcv[0], fRcv[1], fRcv[2], fP2, fP3, fP4, 0.0f, 0.0f, fUnitSlipF, fModPara[2], fModPara[4]);
                                    
                                    fStressStk[0] += fStressStkT[0];        fStressStk[1] += fStressStkT[1];        fStressStk[2] += fStressStkT[2];
                                    fStressStk[3] += fStressStkT[3];        fStressStk[4] += fStressStkT[4];        fStressStk[5] += fStressStkT[5];
                                    fStressDip[0] += fStressDipT[0];        fStressDip[1] += fStressDipT[1];        fStressDip[2] += fStressDipT[2];
                                    fStressDip[3] += fStressDipT[3];        fStressDip[4] += fStressDipT[4];        fStressDip[5] += fStressDipT[5];
                                    fStressNrm[0] += fStressNrmT[0];        fStressNrm[1] += fStressNrmT[1];        fStressNrm[2] += fStressNrmT[2];
                                    fStressNrm[3] += fStressNrmT[3];        fStressNrm[4] += fStressNrmT[4];        fStressNrm[5] += fStressNrmT[5];
                                }
                                
                                ScaleStress6(fStressStk, (fDistF/fUnitSlipF));          ScaleStress6(fStressStk, (1.0f/fKh_FFBrch[uPrvBrId*21 +13]));           RotStressG2L(&fPrvStress[0], fRotMat, fStressStk);
                                ScaleStress6(fStressDip, (fDistF/fUnitSlipF));          ScaleStress6(fStressDip, (1.0f/fKh_FFBrch[uPrvBrId*21 +13]));           RotStressG2L(&fPrvStress[3], fRotMat, fStressDip);
                                ScaleStress6(fStressNrm, (fDistF/fUnitSlipF));          ScaleStress6(fStressNrm, (1.0f/fKh_FFBrch[uPrvBrId*21 +13]));           RotStressG2L(&fPrvStress[6], fRotMat, fStressNrm);
                                
                                uTemp1 = 1u;
                                for (j = 0; j < uF_OTPrtChld[uPrvBrId*5 +0]; j++)
                                {   uTstBrId = uF_OTPrtChld[uPrvBrId*5 +(j+1)];
                                    
                                    GetGlobVertsForRectangle(fP1, fP2, fP3, fP4, fKh_FFBrch, uTstBrId*21u);
                                    
                                    if (USEHALFSPACE == 1u)
                                    {   StrainHS_Nikkhoo(fStressStkT, fStrain, fRcv[0], fRcv[1], fRcv[2], fP1, fP2, fP4, fUnitSlipF, 0.0f, 0.0f, fModPara[2], fModPara[4]);
                                        StrainHS_Nikkhoo(fStressDipT, fStrain, fRcv[0], fRcv[1], fRcv[2], fP1, fP2, fP4, 0.0f, fUnitSlipF, 0.0f, fModPara[2], fModPara[4]);
                                        StrainHS_Nikkhoo(fStressNrmT, fStrain, fRcv[0], fRcv[1], fRcv[2], fP1, fP2, fP4, 0.0f, 0.0f, fUnitSlipF, fModPara[2], fModPara[4]);
                                        
                                        StrainHS_Nikkhoo(fStressStk, fStrain, fRcv[0], fRcv[1], fRcv[2], fP2, fP3, fP4, fUnitSlipF, 0.0f, 0.0f, fModPara[2], fModPara[4]);
                                        StrainHS_Nikkhoo(fStressDip, fStrain, fRcv[0], fRcv[1], fRcv[2], fP2, fP3, fP4, 0.0f, fUnitSlipF, 0.0f, fModPara[2], fModPara[4]);
                                        StrainHS_Nikkhoo(fStressNrm, fStrain, fRcv[0], fRcv[1], fRcv[2], fP2, fP3, fP4, 0.0f, 0.0f, fUnitSlipF, fModPara[2], fModPara[4]);
                                        
                                        fStressStk[0] += fStressStkT[0];        fStressStk[1] += fStressStkT[1];        fStressStk[2] += fStressStkT[2];
                                        fStressStk[3] += fStressStkT[3];        fStressStk[4] += fStressStkT[4];        fStressStk[5] += fStressStkT[5];
                                        fStressDip[0] += fStressDipT[0];        fStressDip[1] += fStressDipT[1];        fStressDip[2] += fStressDipT[2];
                                        fStressDip[3] += fStressDipT[3];        fStressDip[4] += fStressDipT[4];        fStressDip[5] += fStressDipT[5];
                                        fStressNrm[0] += fStressNrmT[0];        fStressNrm[1] += fStressNrmT[1];        fStressNrm[2] += fStressNrmT[2];
                                        fStressNrm[3] += fStressNrmT[3];        fStressNrm[4] += fStressNrmT[4];        fStressNrm[5] += fStressNrmT[5];
                                    }
                                    else
                                    {   StrainFS_Nikkhoo(fStressStkT, fStrain, fRcv[0], fRcv[1], fRcv[2], fP1, fP2, fP4, fUnitSlipF, 0.0f, 0.0f, fModPara[2], fModPara[4]);
                                        StrainFS_Nikkhoo(fStressDipT, fStrain, fRcv[0], fRcv[1], fRcv[2], fP1, fP2, fP4, 0.0f, fUnitSlipF, 0.0f, fModPara[2], fModPara[4]);
                                        StrainFS_Nikkhoo(fStressNrmT, fStrain, fRcv[0], fRcv[1], fRcv[2], fP1, fP2, fP4, 0.0f, 0.0f, fUnitSlipF, fModPara[2], fModPara[4]);
                                        
                                        StrainFS_Nikkhoo(fStressStk, fStrain, fRcv[0], fRcv[1], fRcv[2], fP2, fP3, fP4, fUnitSlipF, 0.0f, 0.0f, fModPara[2], fModPara[4]);
                                        StrainFS_Nikkhoo(fStressDip, fStrain, fRcv[0], fRcv[1], fRcv[2], fP2, fP3, fP4, 0.0f, fUnitSlipF, 0.0f, fModPara[2], fModPara[4]);
                                        StrainFS_Nikkhoo(fStressNrm, fStrain, fRcv[0], fRcv[1], fRcv[2], fP2, fP3, fP4, 0.0f, 0.0f, fUnitSlipF, fModPara[2], fModPara[4]);
                                        
                                        fStressStk[0] += fStressStkT[0];        fStressStk[1] += fStressStkT[1];        fStressStk[2] += fStressStkT[2];
                                        fStressStk[3] += fStressStkT[3];        fStressStk[4] += fStressStkT[4];        fStressStk[5] += fStressStkT[5];
                                        fStressDip[0] += fStressDipT[0];        fStressDip[1] += fStressDipT[1];        fStressDip[2] += fStressDipT[2];
                                        fStressDip[3] += fStressDipT[3];        fStressDip[4] += fStressDipT[4];        fStressDip[5] += fStressDipT[5];
                                        fStressNrm[0] += fStressNrmT[0];        fStressNrm[1] += fStressNrmT[1];        fStressNrm[2] += fStressNrmT[2];
                                        fStressNrm[3] += fStressNrmT[3];        fStressNrm[4] += fStressNrmT[4];        fStressNrm[5] += fStressNrmT[5];
                                    }
                                    
                                    ScaleStress6(fStressStk, (fDistF/fUnitSlipF));          ScaleStress6(fStressStk, (1.0f/fKh_FFBrch[uTstBrId*21 +13]));           RotStressG2L(&fTstStress[0], fRotMat, fStressStk);
                                    ScaleStress6(fStressDip, (fDistF/fUnitSlipF));          ScaleStress6(fStressDip, (1.0f/fKh_FFBrch[uTstBrId*21 +13]));           RotStressG2L(&fTstStress[3], fRotMat, fStressDip);
                                    ScaleStress6(fStressNrm, (fDistF/fUnitSlipF));          ScaleStress6(fStressNrm, (1.0f/fKh_FFBrch[uTstBrId*21 +13]));           RotStressG2L(&fTstStress[6], fRotMat, fStressNrm);
                                    
                                    memcpy(fSrc, &fKh_FFBrch[uTstBrId*21 +0], 3u*sizeof(float));
                                    SubtractVect3(fTempVect, fRcv, fSrc);
                                    fDist  = VectorLength(fTempVect);
                                    fTemp2 = (fDist/fFltSide - KH_FULL_DIST)/(KH_TOL_DIST - KH_FULL_DIST);
                                    fTemp2 = (fTemp2 >= 0.0f)*fTemp2 + (fTemp2 < 0.0f)*0.0f;
                                    fTemp2 = (fTemp2 <= 1.0f)*fTemp2 + (fTemp2 > 1.0f)*1.0f;
                                    fTemp2 = KH_TOL1 + (KH_TOL2-KH_TOL1)*fTemp2;
                                    
                                    fTemp0 = sqrtf( (fTstStress[0]-fPrvStress[0])*(fTstStress[0]-fPrvStress[0]) + (fTstStress[1]-fPrvStress[1])*(fTstStress[1]-fPrvStress[1]) + (fTstStress[2]-fPrvStress[2])*(fTstStress[2]-fPrvStress[2])  + (fTstStress[3]-fPrvStress[3])*(fTstStress[3]-fPrvStress[3]) + (fTstStress[4]-fPrvStress[4])*(fTstStress[4]-fPrvStress[4]) + (fTstStress[5]-fPrvStress[5])*(fTstStress[5]-fPrvStress[5])   + (fTstStress[6]-fPrvStress[6])*(fTstStress[6]-fPrvStress[6]) + (fTstStress[7]-fPrvStress[7])*(fTstStress[7]-fPrvStress[7]) + (fTstStress[8]-fPrvStress[8])*(fTstStress[8]-fPrvStress[8]));
                                    fTemp1 = sqrtf(                fPrvStress[0] *fPrvStress[0]                 +                fPrvStress[1] *fPrvStress[1]                 +                fPrvStress[2] *fPrvStress[2]                  +                fPrvStress[3] *fPrvStress[3]                 +                fPrvStress[4] *fPrvStress[4]                 +                fPrvStress[5] *fPrvStress[5]                   +                fPrvStress[6] *fPrvStress[6]                 +                fPrvStress[7] *fPrvStress[7]                 +                fPrvStress[8] *fPrvStress[8]);
                                    
                                    uTemp1 = ((fTemp0/fTemp1) <= fTemp2)*uTemp1 + ((fTemp0/fTemp1) > fTemp2)*0u; 
                                }
                                
                                if (uTemp1 == 1u)
                                {   fKh_BFvalstk[i][uKh_BFcnt[i]*2 +0] = fPrvStress[0];         fKh_BFvalstk[i][uKh_BFcnt[i]*2 +1] = fPrvStress[3];
                                    fKh_BFvaldip[i][uKh_BFcnt[i]*2 +0] = fPrvStress[1];         fKh_BFvaldip[i][uKh_BFcnt[i]*2 +1] = fPrvStress[4];
                                    fKh_BFvalnrm[i][uKh_BFcnt[i]*2 +0] = fPrvStress[2];         fKh_BFvalnrm[i][uKh_BFcnt[i]*2 +1] = fPrvStress[5];
                                    
                                    uUsdNum = 0u;
                                    for (j = 0u; j < uFPNum;  j++)
                                    {   uFUsdEl[j]      = (uF_OTElem[j*u_MaxBrLvl +(u_MaxBrLvl -uTemp0 +uPrvBrLv)] == uPrvBrId)*1u           + (uF_OTElem[j*u_MaxBrLvl +(u_MaxBrLvl -uTemp0 +uPrvBrLv)] != uPrvBrId)*uFUsdEl[j];
                                        uKh_BFps2[i][j] = (uF_OTElem[j*u_MaxBrLvl +(u_MaxBrLvl -uTemp0 +uPrvBrLv)] == uPrvBrId)*uKh_BFcnt[i] + (uF_OTElem[j*u_MaxBrLvl +(u_MaxBrLvl -uTemp0 +uPrvBrLv)] != uPrvBrId)*uKh_BFps2[i][j];
                                        uUsdNum        +=  uFUsdEl[j];
                                    }
                                    uKh_BFcnt[i] += 1u;
                                    
                                    for (j = 1u; j < uTemp0; j++)       {       uFUsdBrch[ uF_OTElem[uSrcElem*u_MaxBrLvl +(u_MaxBrLvl -uTemp0 +j)] ] = 1u;       }
                                    uTemp2 = 0u;
                }   }   }   }   }
                
                fKh_BFvalstk[i] = realloc(fKh_BFvalstk[i], 2*uKh_BFcnt[i]* sizeof *fKh_BFvalstk[i] );
                fKh_BFvaldip[i] = realloc(fKh_BFvaldip[i], 2*uKh_BFcnt[i]* sizeof *fKh_BFvaldip[i] );
                fKh_BFvalnrm[i] = realloc(fKh_BFvalnrm[i], 2*uKh_BFcnt[i]* sizeof *fKh_BFvalnrm[i] );
                
                if ((uPlotCatalog2Screen == 1)&&(iRANK == 0))             {    fprintf(stdout,"   and Bnd/Flt:  %u/%u ",uKh_BFcnt[i], uFPNum);               }
                
                
            }
            
            
            MPI_Allreduce(MPI_IN_PLACE, uF_temp, 7*uFPNum, MPI_UNSIGNED, MPI_MAX, MPI_COMM_WORLD);
            
            timer0 = clock() - timer0;        time_taken  = ((double)timer0)/CLOCKS_PER_SEC;
            MPI_Allreduce(MPI_IN_PLACE, &time_taken, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); 
            time_taken /= (double)iSIZE;
            if (iRANK == 0)    {   fprintf(stdout,"\nTotal RunTime in seconds:    %6.4f\n",time_taken);     }
        } 
        MPI_Barrier( MPI_COMM_WORLD );
        
        
        if (uLoadPrev_Khmat == 1)
        {   if (iRANK  == 0)    {   fprintf(stdout,"now storing the Kh_ matrix\n");       }
            float *FTempVect = calloc(iFOFFSET[iRANK], sizeof *FTempVect);
            float *BTempVect = calloc(iBOFFSET[iRANK], sizeof *BTempVect);
            unsigned long long *lRankOffset = calloc(iSIZE, sizeof *lRankOffset );
            unsigned long long *lRankStart  = calloc(iSIZE, sizeof *lRankStart );
            unsigned long long  lLoclOffset, lTemp0;
            
            MPI_File_delete(cKhmatFileName, MPI_INFO_NULL);
            MPI_File_open(MPI_COMM_WORLD, cKhmatFileName, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fp_KHMAT);  
            if (iRANK == 0)
            {   MPI_File_write(fp_KHMAT,  &uFPNum,     1, MPI_UNSIGNED,  &STATUS);      MPI_File_write(fp_KHMAT,  &uFVNum,      1, MPI_UNSIGNED,  &STATUS);
                MPI_File_write(fp_KHMAT,  &uBPNum,     1, MPI_UNSIGNED,  &STATUS);      MPI_File_write(fp_KHMAT,  &uBVNum,      1, MPI_UNSIGNED,  &STATUS);
            } 
            OFFSETall = (4 * sizeof(unsigned int));
            
            MPI_File_write_at(fp_KHMAT,  (OFFSETall + iFSTART[iRANK]*sizeof(unsigned int)),  uKh_FFcnt,  iFOFFSET[iRANK], MPI_UNSIGNED,   &STATUS);         OFFSETall += (uFPNum * sizeof(unsigned int));
            MPI_File_write_at(fp_KHMAT,  (OFFSETall + iFSTART[iRANK]*sizeof(unsigned int)),  uKh_FBcnt,  iFOFFSET[iRANK], MPI_UNSIGNED,   &STATUS);         OFFSETall += (uFPNum * sizeof(unsigned int));
            MPI_File_write_at(fp_KHMAT,  (OFFSETall + iBSTART[iRANK]*sizeof(unsigned int)),  uKh_BBcnt,  iBOFFSET[iRANK], MPI_UNSIGNED,   &STATUS);         OFFSETall += (uBPNum * sizeof(unsigned int));
            MPI_File_write_at(fp_KHMAT,  (OFFSETall + iBSTART[iRANK]*sizeof(unsigned int)),  uKh_BFcnt,  iBOFFSET[iRANK], MPI_UNSIGNED,   &STATUS);         OFFSETall += (uBPNum * sizeof(unsigned int));
      
            for (i = 0u; i < iFOFFSET[iRANK]; i++)  {       FTempVect[i] = fFEvent[i*17 +4];      }
            MPI_File_write_at(fp_KHMAT,  (OFFSETall + iFSTART[iRANK]*sizeof(float)),  FTempVect,  iFOFFSET[iRANK], MPI_FLOAT,   &STATUS);                   OFFSETall += (uFPNum * sizeof(float));
            for (i = 0u; i < iFOFFSET[iRANK]; i++)  {       FTempVect[i] = fFEvent[i*17 +5];      }
            MPI_File_write_at(fp_KHMAT,  (OFFSETall + iFSTART[iRANK]*sizeof(float)),  FTempVect,  iFOFFSET[iRANK], MPI_FLOAT,   &STATUS);                   OFFSETall += (uFPNum * sizeof(float));
            for (i = 0u; i < iFOFFSET[iRANK]; i++)  {       FTempVect[i] = fFEvent[i*17 +6];      }
            MPI_File_write_at(fp_KHMAT,  (OFFSETall + iFSTART[iRANK]*sizeof(float)),  FTempVect,  iFOFFSET[iRANK], MPI_FLOAT,   &STATUS);                   OFFSETall += (uFPNum * sizeof(float));
            for (i = 0u; i < iBOFFSET[iRANK]; i++)  {       BTempVect[i] = fBEvent[i*9 +4];      }
            MPI_File_write_at(fp_KHMAT,  (OFFSETall + iBSTART[iRANK]*sizeof(float)),  BTempVect,  iBOFFSET[iRANK], MPI_FLOAT,   &STATUS);                   OFFSETall += (uBPNum * sizeof(float));
            for (i = 0u; i < iBOFFSET[iRANK]; i++)  {       BTempVect[i] = fBEvent[i*9 +5];      }
            MPI_File_write_at(fp_KHMAT,  (OFFSETall + iBSTART[iRANK]*sizeof(float)),  BTempVect,  iBOFFSET[iRANK], MPI_FLOAT,   &STATUS);                   OFFSETall += (uBPNum * sizeof(float));
            for (i = 0u; i < iBOFFSET[iRANK]; i++)  {       BTempVect[i] = fBEvent[i*9 +6];      }
            MPI_File_write_at(fp_KHMAT,  (OFFSETall + iBSTART[iRANK]*sizeof(float)),  BTempVect,  iBOFFSET[iRANK], MPI_FLOAT,   &STATUS);                   OFFSETall += (uBPNum * sizeof(float));
            
            lTemp0 = 0;         memset(lRankOffset, 0, iSIZE*sizeof(unsigned long long ));
            for (i = 0u; i < iFOFFSET[iRANK]; i++)  {   lTemp0 += ( uKh_FFcnt[i]*(                          6*sizeof(float) ) + uFPNum*sizeof(unsigned int));                   }
            lRankOffset[iRANK] = lTemp0;
            MPI_Allreduce(MPI_IN_PLACE, &lTemp0,       1,   MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, lRankOffset, iSIZE, MPI_UNSIGNED_LONG_LONG, MPI_MAX, MPI_COMM_WORLD);
            for (i = 1u; i < iSIZE; i++)            {   lRankStart[i]       = lRankStart[i-1] + lRankOffset[i-1];    }
            
            lLoclOffset = 0;
            for (i = 0u; i < iFOFFSET[iRANK]; i++)  
            {   MPI_File_write_at(fp_KHMAT,  (OFFSETall + lRankStart[iRANK] + lLoclOffset),     uKh_FFps2[i],   uFPNum,       MPI_UNSIGNED,   &STATUS);         lLoclOffset += (uFPNum        *sizeof(unsigned int));
                MPI_File_write_at(fp_KHMAT,  (OFFSETall + lRankStart[iRANK] + lLoclOffset),  fKh_FFvalstk[i], 2*uKh_FFcnt[i], MPI_FLOAT,      &STATUS);         lLoclOffset += (uKh_FFcnt[i]*2*sizeof(float));
                MPI_File_write_at(fp_KHMAT,  (OFFSETall + lRankStart[iRANK] + lLoclOffset),  fKh_FFvaldip[i], 2*uKh_FFcnt[i], MPI_FLOAT,      &STATUS);         lLoclOffset += (uKh_FFcnt[i]*2*sizeof(float));
                MPI_File_write_at(fp_KHMAT,  (OFFSETall + lRankStart[iRANK] + lLoclOffset),  fKh_FFvalnrm[i], 2*uKh_FFcnt[i], MPI_FLOAT,      &STATUS);         lLoclOffset += (uKh_FFcnt[i]*2*sizeof(float));
            }
            OFFSETall += lTemp0;
            
            lTemp0 = 0;         memset(lRankOffset, 0, iSIZE*sizeof(unsigned long long ));
            for (i = 0u; i < iFOFFSET[iRANK]; i++)  {   lTemp0 += ( uKh_FBcnt[i]*( 6*sizeof(float)) + uBPNum*sizeof(unsigned int));     }
            lRankOffset[iRANK] = lTemp0; 
            MPI_Allreduce(MPI_IN_PLACE, &lTemp0,       1,   MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, lRankOffset, iSIZE, MPI_UNSIGNED_LONG_LONG, MPI_MAX, MPI_COMM_WORLD);
            for (i = 1u; i < iSIZE; i++)            {   lRankStart[i]       = lRankStart[i-1] + lRankOffset[i-1];    }
            
            lLoclOffset = 0;
            for (i = 0u; i < iFOFFSET[iRANK]; i++)  
            {   MPI_File_write_at(fp_KHMAT,  (OFFSETall + lRankStart[iRANK] + lLoclOffset),  uKh_FBps2[i],      uBPNum,       MPI_UNSIGNED,   &STATUS);         lLoclOffset += (uBPNum        *sizeof(unsigned int));
                MPI_File_write_at(fp_KHMAT,  (OFFSETall + lRankStart[iRANK] + lLoclOffset),  fKh_FBvalstk[i], 3*uKh_FBcnt[i], MPI_FLOAT,      &STATUS);         lLoclOffset += (uKh_FBcnt[i]*3*sizeof(float));
                MPI_File_write_at(fp_KHMAT,  (OFFSETall + lRankStart[iRANK] + lLoclOffset),  fKh_FBvaldip[i], 3*uKh_FBcnt[i], MPI_FLOAT,      &STATUS);         lLoclOffset += (uKh_FBcnt[i]*3*sizeof(float));
            }
            OFFSETall += lTemp0;
            
            lTemp0 = 0;         memset(lRankOffset, 0, iSIZE*sizeof(unsigned long long ));
            for (i = 0u; i < iBOFFSET[iRANK]; i++)  {   lTemp0 += ( uKh_BBcnt[i]*( 9*sizeof(float)) + uBPNum*sizeof(unsigned int));     }
            lRankOffset[iRANK] = lTemp0; 
            MPI_Allreduce(MPI_IN_PLACE, &lTemp0,       1,   MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, lRankOffset, iSIZE, MPI_UNSIGNED_LONG_LONG, MPI_MAX, MPI_COMM_WORLD);
            for (i = 1u; i < iSIZE; i++)            {   lRankStart[i]       = lRankStart[i-1] + lRankOffset[i-1];    }
            
            lLoclOffset = 0;
            for (i = 0u; i < iBOFFSET[iRANK]; i++)  
            {   MPI_File_write_at(fp_KHMAT,  (OFFSETall + lRankStart[iRANK] + lLoclOffset),  uKh_BBps2[i],      uBPNum,       MPI_UNSIGNED,   &STATUS);         lLoclOffset += (uBPNum        *sizeof(unsigned int));
                MPI_File_write_at(fp_KHMAT,  (OFFSETall + lRankStart[iRANK] + lLoclOffset),  fKh_BBvalstk[i], 3*uKh_BBcnt[i], MPI_FLOAT,      &STATUS);         lLoclOffset += (uKh_BBcnt[i]*3*sizeof(float));
                MPI_File_write_at(fp_KHMAT,  (OFFSETall + lRankStart[iRANK] + lLoclOffset),  fKh_BBvaldip[i], 3*uKh_BBcnt[i], MPI_FLOAT,      &STATUS);         lLoclOffset += (uKh_BBcnt[i]*3*sizeof(float));
                MPI_File_write_at(fp_KHMAT,  (OFFSETall + lRankStart[iRANK] + lLoclOffset),  fKh_BBvalnrm[i], 3*uKh_BBcnt[i], MPI_FLOAT,      &STATUS);         lLoclOffset += (uKh_BBcnt[i]*3*sizeof(float));
            }
            OFFSETall += lTemp0;
            
            lTemp0 = 0;         memset(lRankOffset, 0, iSIZE*sizeof(unsigned long long ));
            for (i = 0u; i < iBOFFSET[iRANK]; i++)  {   lTemp0 += ( uKh_BFcnt[i]*( 6*sizeof(float)) + uFPNum*sizeof(unsigned int));     }
            lRankOffset[iRANK] = lTemp0; 
            MPI_Allreduce(MPI_IN_PLACE, &lTemp0,       1,   MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, lRankOffset, iSIZE, MPI_UNSIGNED_LONG_LONG, MPI_MAX, MPI_COMM_WORLD);
            for (i = 1u; i < iSIZE; i++)            {   lRankStart[i]       = lRankStart[i-1] + lRankOffset[i-1];    }
            
            lLoclOffset = 0;
            for (i = 0u; i < iBOFFSET[iRANK]; i++)  
            {   MPI_File_write_at(fp_KHMAT,  (OFFSETall + lRankStart[iRANK] + lLoclOffset),  uKh_BFps2[i],      uFPNum,       MPI_UNSIGNED,   &STATUS);         lLoclOffset += (uFPNum        *sizeof(unsigned int));
                MPI_File_write_at(fp_KHMAT,  (OFFSETall + lRankStart[iRANK] + lLoclOffset),  fKh_BFvalstk[i], 2*uKh_BFcnt[i], MPI_FLOAT,      &STATUS);         lLoclOffset += (uKh_BFcnt[i]*2*sizeof(float));
                MPI_File_write_at(fp_KHMAT,  (OFFSETall + lRankStart[iRANK] + lLoclOffset),  fKh_BFvaldip[i], 2*uKh_BFcnt[i], MPI_FLOAT,      &STATUS);         lLoclOffset += (uKh_BFcnt[i]*2*sizeof(float));
                MPI_File_write_at(fp_KHMAT,  (OFFSETall + lRankStart[iRANK] + lLoclOffset),  fKh_BFvalnrm[i], 2*uKh_BFcnt[i], MPI_FLOAT,      &STATUS);         lLoclOffset += (uKh_BFcnt[i]*2*sizeof(float));
            }
            
            MPI_Barrier( MPI_COMM_WORLD );
            MPI_File_close(&fp_KHMAT);
        } 
        
        
    } 
    else
    {   if (iRANK  == 0)    {   fprintf(stdout,"now loading the Kh_ matrix\n");       }
        float *FTempVect = calloc(iFOFFSET[iRANK], sizeof *FTempVect);
        float *BTempVect = calloc(iBOFFSET[iRANK], sizeof *BTempVect);
        unsigned int uFPNum_t,   uBPNum_t, uFVNum_t, uBVNum_t;
        unsigned long long *lRankOffset = calloc(iSIZE, sizeof *lRankOffset );
        unsigned long long *lRankStart  = calloc(iSIZE, sizeof *lRankStart );
        unsigned long long  lLoclOffset, lTemp0;
        float fTemp0,   fTemp1,   fTemp2;
                    
        MPI_File_open(MPI_COMM_WORLD, cKhmatFileName, MPI_MODE_RDONLY, MPI_INFO_NULL, &fp_KHMAT);  
        
        MPI_File_read(fp_KHMAT, &uFPNum_t,   1, MPI_UNSIGNED, &STATUS);       MPI_File_read(fp_KHMAT, &uFVNum_t,    1, MPI_UNSIGNED, &STATUS);
        MPI_File_read(fp_KHMAT, &uBPNum_t,   1, MPI_UNSIGNED, &STATUS);       MPI_File_read(fp_KHMAT, &uBVNum_t,    1, MPI_UNSIGNED, &STATUS);
        if ((uFPNum_t != uFPNum) || (uBPNum_t != uBPNum) || (uFVNum_t != uFVNum) || (uBVNum_t != uBVNum))
        {   fprintf(stdout,"Kh_matrix file was openend, but fault/boundary element numbers not matching => exit\n");        exit(10);    }
        
        OFFSETall = (4 * sizeof(unsigned int)) ;
        
        MPI_File_read_at(fp_KHMAT, (OFFSETall + iFSTART[iRANK]*sizeof(unsigned int)), uKh_FFcnt, iFOFFSET[iRANK], MPI_UNSIGNED, &STATUS);       OFFSETall += (uFPNum * sizeof(unsigned int));
        MPI_File_read_at(fp_KHMAT, (OFFSETall + iFSTART[iRANK]*sizeof(unsigned int)), uKh_FBcnt, iFOFFSET[iRANK], MPI_UNSIGNED, &STATUS);       OFFSETall += (uFPNum * sizeof(unsigned int));
        MPI_File_read_at(fp_KHMAT, (OFFSETall + iBSTART[iRANK]*sizeof(unsigned int)), uKh_BBcnt, iBOFFSET[iRANK], MPI_UNSIGNED, &STATUS);       OFFSETall += (uBPNum * sizeof(unsigned int));
        MPI_File_read_at(fp_KHMAT, (OFFSETall + iBSTART[iRANK]*sizeof(unsigned int)), uKh_BFcnt, iBOFFSET[iRANK], MPI_UNSIGNED, &STATUS);       OFFSETall += (uBPNum * sizeof(unsigned int));
  
        MPI_File_read_at(fp_KHMAT, (OFFSETall + iFSTART[iRANK]*sizeof(float)), FTempVect, iFOFFSET[iRANK], MPI_FLOAT, &STATUS);       OFFSETall += (uFPNum * sizeof(float));
        for (i = 0u; i < iFOFFSET[iRANK]; i++)  {   fFEvent[i*17 +4] = FTempVect[i];     }
        MPI_File_read_at(fp_KHMAT, (OFFSETall + iFSTART[iRANK]*sizeof(float)), FTempVect, iFOFFSET[iRANK], MPI_FLOAT, &STATUS);       OFFSETall += (uFPNum * sizeof(float));
        for (i = 0u; i < iFOFFSET[iRANK]; i++)  {   fFEvent[i*17 +5] = FTempVect[i];     }
        MPI_File_read_at(fp_KHMAT, (OFFSETall + iFSTART[iRANK]*sizeof(float)), FTempVect, iFOFFSET[iRANK], MPI_FLOAT, &STATUS);       OFFSETall += (uFPNum * sizeof(float));
        for (i = 0u; i < iFOFFSET[iRANK]; i++)  {   fFEvent[i*17 +6] = FTempVect[i];     }
        MPI_File_read_at(fp_KHMAT, (OFFSETall + iBSTART[iRANK]*sizeof(float)), BTempVect, iBOFFSET[iRANK], MPI_FLOAT, &STATUS);       OFFSETall += (uBPNum * sizeof(float));
        for (i = 0u; i < iBOFFSET[iRANK]; i++)  {   fBEvent[i*9 +4] = BTempVect[i];     }
        MPI_File_read_at(fp_KHMAT, (OFFSETall + iBSTART[iRANK]*sizeof(float)), BTempVect, iBOFFSET[iRANK], MPI_FLOAT, &STATUS);       OFFSETall += (uBPNum * sizeof(float));
        for (i = 0u; i < iBOFFSET[iRANK]; i++)  {   fBEvent[i*9 +5] = BTempVect[i];     }
        MPI_File_read_at(fp_KHMAT, (OFFSETall + iBSTART[iRANK]*sizeof(float)), BTempVect, iBOFFSET[iRANK], MPI_FLOAT, &STATUS);       OFFSETall += (uBPNum * sizeof(float));
        for (i = 0u; i < iBOFFSET[iRANK]; i++)  {   fBEvent[i*9 +6] = BTempVect[i];     }
        
        lTemp0 = 0;         memset(lRankOffset, 0, iSIZE*sizeof(unsigned long long ));
        
        for (i = 0u; i < iFOFFSET[iRANK]; i++) 
        {   uKh_FFps2[i]    = calloc(     uFPNum,       sizeof *uKh_FFps2[i] );
            fKh_FFvalstk[i] = calloc(   2*uKh_FFcnt[i], sizeof *fKh_FFvalstk[i] );
            fKh_FFvaldip[i] = calloc(   2*uKh_FFcnt[i], sizeof *fKh_FFvaldip[i] );
            fKh_FFvalnrm[i] = calloc(   2*uKh_FFcnt[i], sizeof *fKh_FFvalnrm[i] );
            lTemp0 += (uKh_FFcnt[i]*(                           6*sizeof(float) ) + uFPNum*sizeof(unsigned int));
        }
        lRankOffset[iRANK] = lTemp0;
        MPI_Allreduce(MPI_IN_PLACE, &lTemp0,       1,   MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, lRankOffset, iSIZE, MPI_UNSIGNED_LONG_LONG, MPI_MAX, MPI_COMM_WORLD);
        for (i = 1u; i < iSIZE; i++)            {   lRankStart[i]       = lRankStart[i-1] + lRankOffset[i-1];    }
        
        lLoclOffset = 0;
        for (i = 0u; i < iFOFFSET[iRANK]; i++)  
        {   MPI_File_read_at(fp_KHMAT,  (OFFSETall + lRankStart[iRANK] + lLoclOffset),     uKh_FFps2[i],   uFPNum,       MPI_UNSIGNED,   &STATUS);         lLoclOffset += (uFPNum  *sizeof(unsigned int));
            MPI_File_read_at(fp_KHMAT,  (OFFSETall + lRankStart[iRANK] + lLoclOffset),  fKh_FFvalstk[i], 2*uKh_FFcnt[i], MPI_FLOAT,      &STATUS);         lLoclOffset += (uKh_FFcnt[i]*2*sizeof(float));
            MPI_File_read_at(fp_KHMAT,  (OFFSETall + lRankStart[iRANK] + lLoclOffset),  fKh_FFvaldip[i], 2*uKh_FFcnt[i], MPI_FLOAT,      &STATUS);         lLoclOffset += (uKh_FFcnt[i]*2*sizeof(float));
            MPI_File_read_at(fp_KHMAT,  (OFFSETall + lRankStart[iRANK] + lLoclOffset),  fKh_FFvalnrm[i], 2*uKh_FFcnt[i], MPI_FLOAT,      &STATUS);         lLoclOffset += (uKh_FFcnt[i]*2*sizeof(float));
        }
        OFFSETall += lTemp0;
        
        lTemp0 = 0;         memset(lRankOffset, 0, iSIZE*sizeof(unsigned long long ));
        for (i = 0u; i < iFOFFSET[iRANK]; i++) 
        {   uKh_FBps2[i]    = calloc(     uBPNum,       sizeof *uKh_FBps2[i] );
            fKh_FBvalstk[i] = calloc(   3*uKh_FBcnt[i], sizeof *fKh_FBvalstk[i] );
            fKh_FBvaldip[i] = calloc(   3*uKh_FBcnt[i], sizeof *fKh_FBvaldip[i] );
            lTemp0 += (uKh_FBcnt[i]*( 6*sizeof(float)) + uBPNum*sizeof(unsigned int));
        }
        lRankOffset[iRANK] = lTemp0;
        MPI_Allreduce(MPI_IN_PLACE, &lTemp0,       1,   MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, lRankOffset, iSIZE, MPI_UNSIGNED_LONG_LONG, MPI_MAX, MPI_COMM_WORLD);
        for (i = 1u; i < iSIZE; i++)            {   lRankStart[i]       = lRankStart[i-1] + lRankOffset[i-1];    }
        
        lLoclOffset = 0;
        for (i = 0u; i < iFOFFSET[iRANK]; i++)  
        {   MPI_File_read_at(fp_KHMAT,  (OFFSETall + lRankStart[iRANK] + lLoclOffset),     uKh_FBps2[i],   uBPNum,       MPI_UNSIGNED,   &STATUS);         lLoclOffset += (uBPNum  *sizeof(unsigned int));
            MPI_File_read_at(fp_KHMAT,  (OFFSETall + lRankStart[iRANK] + lLoclOffset),  fKh_FBvalstk[i], 3*uKh_FBcnt[i], MPI_FLOAT,      &STATUS);         lLoclOffset += (uKh_FBcnt[i]*3*sizeof(float));
            MPI_File_read_at(fp_KHMAT,  (OFFSETall + lRankStart[iRANK] + lLoclOffset),  fKh_FBvaldip[i], 3*uKh_FBcnt[i], MPI_FLOAT,      &STATUS);         lLoclOffset += (uKh_FBcnt[i]*3*sizeof(float));
        }
        OFFSETall += lTemp0;
        
        lTemp0 = 0;         memset(lRankOffset, 0, iSIZE*sizeof(unsigned long long ));
        for (i = 0u; i < iBOFFSET[iRANK]; i++) 
        {   uKh_BBps2[i]    = calloc(     uBPNum,       sizeof *uKh_BBps2[i] );
            fKh_BBvalstk[i] = calloc(   3*uKh_BBcnt[i], sizeof *fKh_BBvalstk[i] );
            fKh_BBvaldip[i] = calloc(   3*uKh_BBcnt[i], sizeof *fKh_BBvaldip[i] );
            fKh_BBvalnrm[i] = calloc(   3*uKh_BBcnt[i], sizeof *fKh_BBvalnrm[i] );
            lTemp0 += (uKh_BBcnt[i]*( 9*sizeof(float)) + uBPNum*sizeof(unsigned int));
        }
        lRankOffset[iRANK] = lTemp0;
        MPI_Allreduce(MPI_IN_PLACE, &lTemp0,       1,   MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, lRankOffset, iSIZE, MPI_UNSIGNED_LONG_LONG, MPI_MAX, MPI_COMM_WORLD);
        for (i = 1u; i < iSIZE; i++)            {   lRankStart[i]       = lRankStart[i-1] + lRankOffset[i-1];    }
        
        lLoclOffset = 0;
        for (i = 0u; i < iBOFFSET[iRANK]; i++)  
        {   MPI_File_read_at(fp_KHMAT,  (OFFSETall + lRankStart[iRANK] + lLoclOffset),     uKh_BBps2[i],   uBPNum,       MPI_UNSIGNED,   &STATUS);         lLoclOffset += (uBPNum  *sizeof(unsigned int));
            MPI_File_read_at(fp_KHMAT,  (OFFSETall + lRankStart[iRANK] + lLoclOffset),  fKh_BBvalstk[i], 3*uKh_BBcnt[i], MPI_FLOAT,      &STATUS);         lLoclOffset += (uKh_BBcnt[i]*3*sizeof(float));
            MPI_File_read_at(fp_KHMAT,  (OFFSETall + lRankStart[iRANK] + lLoclOffset),  fKh_BBvaldip[i], 3*uKh_BBcnt[i], MPI_FLOAT,      &STATUS);         lLoclOffset += (uKh_BBcnt[i]*3*sizeof(float));
            MPI_File_read_at(fp_KHMAT,  (OFFSETall + lRankStart[iRANK] + lLoclOffset),  fKh_BBvalnrm[i], 3*uKh_BBcnt[i], MPI_FLOAT,      &STATUS);         lLoclOffset += (uKh_BBcnt[i]*3*sizeof(float));
        }
        OFFSETall += lTemp0;
        
        lTemp0 = 0;         memset(lRankOffset, 0, iSIZE*sizeof(unsigned long long ));
        for (i = 0u; i < iBOFFSET[iRANK]; i++) 
        {   uKh_BFps2[i]    = calloc(     uFPNum,       sizeof *uKh_BFps2[i] );
            fKh_BFvalstk[i] = calloc(   2*uKh_BFcnt[i], sizeof *fKh_BFvalstk[i] );
            fKh_BFvaldip[i] = calloc(   2*uKh_BFcnt[i], sizeof *fKh_BFvaldip[i] );
            fKh_BFvalnrm[i] = calloc(   2*uKh_BFcnt[i], sizeof *fKh_BFvalnrm[i] );
            lTemp0 += (uKh_BFcnt[i]*( 6*sizeof(float)) + uFPNum*sizeof(unsigned int));
        }
        lRankOffset[iRANK] = lTemp0;
        MPI_Allreduce(MPI_IN_PLACE, &lTemp0,       1,   MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, lRankOffset, iSIZE, MPI_UNSIGNED_LONG_LONG, MPI_MAX, MPI_COMM_WORLD);
        for (i = 1u; i < iSIZE; i++)            {   lRankStart[i]       = lRankStart[i-1] + lRankOffset[i-1];    }
        
        lLoclOffset = 0;
        for (i = 0u; i < iBOFFSET[iRANK]; i++)  
        {   MPI_File_read_at(fp_KHMAT,  (OFFSETall + lRankStart[iRANK] + lLoclOffset),     uKh_BFps2[i],   uFPNum,       MPI_UNSIGNED,   &STATUS);         lLoclOffset += (uFPNum  *sizeof(unsigned int));
            MPI_File_read_at(fp_KHMAT,  (OFFSETall + lRankStart[iRANK] + lLoclOffset),  fKh_BFvalstk[i], 2*uKh_BFcnt[i], MPI_FLOAT,      &STATUS);         lLoclOffset += (uKh_BFcnt[i]*2*sizeof(float));
            MPI_File_read_at(fp_KHMAT,  (OFFSETall + lRankStart[iRANK] + lLoclOffset),  fKh_BFvaldip[i], 2*uKh_BFcnt[i], MPI_FLOAT,      &STATUS);         lLoclOffset += (uKh_BFcnt[i]*2*sizeof(float));
            MPI_File_read_at(fp_KHMAT,  (OFFSETall + lRankStart[iRANK] + lLoclOffset),  fKh_BFvalnrm[i], 2*uKh_BFcnt[i], MPI_FLOAT,      &STATUS);         lLoclOffset += (uKh_BFcnt[i]*2*sizeof(float));
        }
        OFFSETall += lTemp0;
        
        MPI_Barrier( MPI_COMM_WORLD );
        MPI_File_close(&fp_KHMAT);
        
        for (i = 0u; i < iFOFFSET[iRANK]; i++)  
        {   
            fFRef[i*14 +2]    = MAX(fabs(fFEvent[i*17 +4]),fabs(fFEvent[i*17 +5])); 
            fFRef[i*14 +2]    = MAX(fFRef[i*14 +2],      fabs(fFEvent[i*17 +6])); 
            
            fTemp0 = ran0(&lSeed);      fTemp0 = fTemp0*2.0f -1.0f;
            fFEvent[i*17 +7] = fFRef[i*14 +7] + (fFRef[i*14 +8] * fTemp0); 
            
            fTemp0 = ran0(&lSeed);      fTemp0 = fTemp0*2.0f -1.0f;
            fFFric[i*6 +0] = fFRef[i*14 +3] + (fFRef[i*14 +4] * fTemp0); 
            fFEvent[i*17 +3] = fFFric[i*6 +0]; 
            
            fTemp0 = ran0(&lSeed);      fTemp0 = fTemp0*2.0f -1.0f;
            fFFric[i*6 +1] = fFRef[i*14 +5] + (fFRef[i*14 +6] * fTemp0); 
            
            fTemp0 =(fFFric[i*6 +0] - fFFric[i*6 +1]) *-1.0f*fFRef[i*14 +1]; 
            fTemp1 = fTemp0/fFRef[i*14 +2];
            uFEvent[i*5 +0] = (fTemp0 >= 0.0f)*2u + (fTemp0 < 0.0f)*3u;
            uFEvent[i*5 +0] = (fTemp1 > fFEvent[i*17 +7])*1u + (fTemp1 <= fFEvent[i*17 +7])*uFEvent[i*5 +0];
            
            fTemp0 = fFFric[i*6 +1]*-1.0f*fFRef[i*14 +1] + fFRef[i*14 +2]*fFEvent[i*17 +7]; 
            fTemp1 = fFFric[i*6 +0]*-1.0f*fFRef[i*14 +1]; 
            fTemp2 = MAX(fTemp0, fTemp1)/(-1.0f*fFRef[i*14 +1]); 
            fFFric[i*6 +2]  = fTemp2 - fOvershootFract*(fTemp2-fFFric[i*6 +1]); 
            
            fFEvent[i*17 +10] = (fFFric[i*6 +1] - fFFric[i*6 +2])*-1.0f*fFRef[i*14 +1];
            
        }
    } 
    
    unsigned int uMeanFFKlgth = 0u,          uMeanFBKlgth = 0u,          uMeanBBKlgth = 0u,           uMeanBFKlgth = 0u;
    unsigned int uMeanKLgths[4] = {0u,       0u,       0u,       0u};
    unsigned int uMin_KLgths[4] = {UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX};
    unsigned int uMax_KLgths[4] = {0u,       0u,       0u,       0u};
    for (i = 0u; i < iFOFFSET[iRANK]; i++)      
    {   uMeanKLgths[0] += uKh_FFcnt[i];      uMin_KLgths[0] = (uMin_KLgths[0] < uKh_FFcnt[i])*uMin_KLgths[0]  +  (uMin_KLgths[0] >= uKh_FFcnt[i])*uKh_FFcnt[i];      uMax_KLgths[0] = (uMax_KLgths[0] > uKh_FFcnt[i])*uMax_KLgths[0]  +  (uMax_KLgths[0] < uKh_FFcnt[i])*uKh_FFcnt[i];
        uMeanKLgths[1] += uKh_FBcnt[i];      uMin_KLgths[1] = (uMin_KLgths[1] < uKh_FBcnt[i])*uMin_KLgths[1]  +  (uMin_KLgths[1] >= uKh_FBcnt[i])*uKh_FBcnt[i];      uMax_KLgths[1] = (uMax_KLgths[1] > uKh_FBcnt[i])*uMax_KLgths[1]  +  (uMax_KLgths[1] < uKh_FBcnt[i])*uKh_FBcnt[i];
    }
    for (i = 0u; i < iBOFFSET[iRANK]; i++)      
    {   uMeanKLgths[2] += uKh_BBcnt[i];      uMin_KLgths[2] = (uMin_KLgths[2] < uKh_BBcnt[i])*uMin_KLgths[2]  +  (uMin_KLgths[2] >= uKh_BBcnt[i])*uKh_BBcnt[i];      uMax_KLgths[2] = (uMax_KLgths[2] > uKh_BBcnt[i])*uMax_KLgths[2]  +  (uMax_KLgths[2] < uKh_BBcnt[i])*uKh_BBcnt[i];
        uMeanKLgths[3] += uKh_BFcnt[i];      uMin_KLgths[3] = (uMin_KLgths[3] < uKh_BFcnt[i])*uMin_KLgths[3]  +  (uMin_KLgths[3] >= uKh_BFcnt[i])*uKh_BFcnt[i];      uMax_KLgths[3] = (uMax_KLgths[3] > uKh_BFcnt[i])*uMax_KLgths[3]  +  (uMax_KLgths[3] < uKh_BFcnt[i])*uKh_BFcnt[i];
    }
    uMeanKLgths[0] =                             uMeanKLgths[0]/iFOFFSET[iRANK];
    uMeanKLgths[1] =                             uMeanKLgths[1]/iFOFFSET[iRANK];
    uMeanKLgths[2] = (iBOFFSET[iRANK] > 0u)  ?  (uMeanKLgths[2]/iBOFFSET[iRANK])  :  0u;
    uMeanKLgths[3] = (iBOFFSET[iRANK] > 0u)  ?  (uMeanKLgths[3]/iBOFFSET[iRANK])  :  0u;
    MPI_Allreduce(MPI_IN_PLACE, &uMeanKLgths, 4, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &uMin_KLgths, 4, MPI_UNSIGNED, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &uMax_KLgths, 4, MPI_UNSIGNED, MPI_MAX, MPI_COMM_WORLD);
    uMeanKLgths[0] /= iSIZE;              uMeanKLgths[1] /= iSIZE;              uMeanKLgths[2] /= iSIZE;              uMeanKLgths[3] /= iSIZE;
    uMeanFFKlgth    = uMeanKLgths[0];     uMeanFBKlgth    = uMeanKLgths[1];     uMeanBBKlgth    = uMeanKLgths[2];     uMeanBFKlgth    = uMeanKLgths[3];
    if (iRANK == 0)
    {   fprintf(stdout,"\nAverage Kh_lengths:  FF %u    FB %u    BB %u    BF %u\n", uMeanFFKlgth, uMeanFBKlgth, uMeanBBKlgth, uMeanBFKlgth);
        fprintf(stdout,  "Minimum Kh_lengths:  FF %u    FB %u    BB %u    BF %u\n", uMin_KLgths[0], uMin_KLgths[1], uMin_KLgths[2], uMin_KLgths[3]);
        fprintf(stdout,  "Maximum Kh_lengths:  FF %u    FB %u    BB %u    BF %u\n\n", uMax_KLgths[0], uMax_KLgths[1], uMax_KLgths[2], uMax_KLgths[3]);
    }
    
    
    FILE *fpPrePost,   *fp_CATALOG,   *fp_CATALOGCUT;
    float *fFslip    = calloc(2*uFPNum, sizeof *fFslip); 
    float *fBslip    = calloc(3*uBPNum, sizeof *fBslip); 
    
    unsigned int *uSlipElCnt  = calloc(2, sizeof *uSlipElCnt );
    int *iStartPosF  = calloc(iSIZE, sizeof *iStartPosF);
    int *iOffstPosF  = calloc(iSIZE, sizeof *iOffstPosF);
    int *iStartPosB  = calloc(iSIZE, sizeof *iStartPosB);
    int *iOffstPosB  = calloc(iSIZE, sizeof *iOffstPosB);
    float *fSTFslip  = calloc(2*iFOFFSET[iRANK]*MAXMOMRATEFUNCLENGTH, sizeof *fSTFslip);
    float *fFslipL   = calloc(3*iFOFFSET[iRANK], sizeof *fFslipL);
    float *fFslipG   = calloc(3*uFPNum, sizeof *fFslipG); 
    float *fBslipL   = calloc(4*iBOFFSET[iRANK], sizeof *fBslipL);
    float *fBslipG   = calloc(4*uBPNum, sizeof *fBslipG);
    
    
    {   if (iRANK  == 0)    {   fprintf(stdout,"now determining tectonic loading function\n");       }
        unsigned int uGlobPos,   uTemp0,   uTemp1;
        float fTemp0,   fTemp1,   fTemp2,   fCmbSlp1,   fCmbSlp2;
        uTemp0   = 0u;
        fCmbSlp1 = 0.0f;
        for (i = 0u; i < iFOFFSET[iRANK]; i++)
        {   fTemp0   = fFEvent[i*17 +13];       fTemp1   = fFEvent[i*17 +14];
            fTemp2   = sqrtf(fTemp0*fTemp0 +fTemp1*fTemp1);
            fCmbSlp1+= fTemp2;
            uTemp0   = (fTemp2 <= 0.0f)*uTemp0 + (fTemp2 > 0.0f)*(uTemp0 +1u);
            
            fFEvent[i*17 +8] += (-1.0f*fTemp0 *fFEvent[i*17 +4]);
            fFEvent[i*17 +9] += (-1.0f*fTemp1 *fFEvent[i*17 +5]);
        } 
        
        MPI_Allreduce(MPI_IN_PLACE, &fCmbSlp1, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, &uTemp0,   1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
        fCmbSlp1 = (uTemp0 > 0u) ? (fCmbSlp1/(float)uTemp0)  : 0.0f;
        
        
        if (uBPNum > 0u)
        {   memset(iStartPosB, 0, iSIZE*sizeof(int) );           memset(iOffstPosB, 0, iSIZE*sizeof(int) );
            for (i = 0u; i < iBOFFSET[iRANK]; i++)
            {   uGlobPos = i +iBSTART[iRANK];
                fTemp0   =  -fB_temp[uGlobPos*16 +3];                      fTemp1 =  -fB_temp[uGlobPos*16 +4];          fTemp2 =  -fB_temp[uGlobPos*16 +5];
                if ((fabs(fTemp0) + fabs(fTemp1) + fabs(fTemp2)) > 0.0f)
                {   fBslipL[uSlipElCnt[1]*4 +0] = (float)uGlobPos;      fBslipL[uSlipElCnt[1]*4 +1] = fTemp0;       fBslipL[uSlipElCnt[1]*4 +2] = fTemp1;           fBslipL[uSlipElCnt[1]*4 +3] = fTemp2;
                    uSlipElCnt[1] += 1u;
            }   }
            iOffstPosB[iRANK] =uSlipElCnt[1]*4;
            
            MPI_Allreduce(MPI_IN_PLACE, uSlipElCnt, 2, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD); 
            MPI_Allreduce(MPI_IN_PLACE, iOffstPosB, iSIZE, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
            for (i = 1; i < iSIZE; i++)        {    iStartPosB[i]  = iStartPosB[i-1] +iOffstPosB[i-1];    } 
            MPI_Allgatherv(fBslipL, iOffstPosB[iRANK], MPI_FLOAT, fBslipG, iOffstPosB, iStartPosB, MPI_FLOAT, MPI_COMM_WORLD);
            
            for (i = 0u; i < iFOFFSET[iRANK]; i++)
            {   memset(fBslip, 0, 3*uBPNum*sizeof(float) );
                for (j = 0u; j < uSlipElCnt[1]; j++) 
                {   uTemp0 = (unsigned int)fBslipG[j*4 +0]; 
                    uTemp1 = uKh_FBps2[i][uTemp0]; 
                    fBslip[uTemp1*3 +0] += -fBslipG[j*4 +1];
                    fBslip[uTemp1*3 +1] += -fBslipG[j*4 +2];
                    fBslip[uTemp1*3 +2] += -fBslipG[j*4 +3];
                }
                fFEvent[i*17 +8] += cblas_sdot(3*uKh_FBcnt[i], fBslip, 1, fKh_FBvalstk[i], 1);
                fFEvent[i*17 +9] += cblas_sdot(3*uKh_FBcnt[i], fBslip, 1, fKh_FBvaldip[i], 1);
        }   }
        
        for (i = 0u; i < iFOFFSET[iRANK]; i++)
        {   fFTempVal[i*5 +0] = fFEvent[i*17 +8];          fFTempVal[i*5 +1] = fFEvent[i*17 +9];
            fFTempVal[i*5 +2] = 0.0f;                      fFTempVal[i*5 +3] = 0.0f;
        }
        
        for (k = 0u; k < MAXITERATION4BOUNDARY; k++)
 
        {   memset(uSlipElCnt, 0,  2*sizeof(unsigned int));
            memset(iStartPosF, 0, iSIZE*sizeof(int) );           memset(iOffstPosF, 0, iSIZE*sizeof(int) ); 
            for (i = 0u; i < iFOFFSET[iRANK]; i++)
            {   fTemp0 = -1.0f*fFTempVal[i*5 +0]/fFEvent[i*17 +4];       fFTempVal[i*5 +2] += fTemp0;
                fTemp1 = -1.0f*fFTempVal[i*5 +1]/fFEvent[i*17 +5];       fFTempVal[i*5 +3] += fTemp1;
                
                uGlobPos       = i +iFSTART[iRANK];
                fFslipL[uSlipElCnt[0]*3 +0] = (float)uGlobPos;        fFslipL[uSlipElCnt[0]*3 +1] = fTemp0;           fFslipL[uSlipElCnt[0]*3 +2] = fTemp1;
                uSlipElCnt[0] += 1u;
            }
            iOffstPosF[iRANK] = uSlipElCnt[0]*3;
            
            MPI_Allreduce(MPI_IN_PLACE, uSlipElCnt, 2, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, iOffstPosF, iSIZE, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
            for (i = 1; i < iSIZE; i++)         {       iStartPosF[i]  = iStartPosF[i-1] +iOffstPosF[i-1];         }
            MPI_Allgatherv(fFslipL, iOffstPosF[iRANK], MPI_FLOAT, fFslipG, iOffstPosF, iStartPosF, MPI_FLOAT, MPI_COMM_WORLD);
            
            for (i = 0u; i < iFOFFSET[iRANK]; i++)
            {   memset(fFslip, 0, 2*uFPNum*sizeof(float) );
                for (j = 0u; j < uSlipElCnt[0]; j++) 
                {   uTemp0 = (unsigned int)fFslipG[j*3 +0]; 
                    uTemp1 = uKh_FFps2[i][uTemp0]; 
                    fFslip[uTemp1*2 +0] += fFslipG[j*3 +1];
                    fFslip[uTemp1*2 +1] += fFslipG[j*3 +2];
                }
                fFTempVal[i*5 +0] += cblas_sdot(2*uKh_FFcnt[i], fFslip, 1, fKh_FFvalstk[i], 1);
                fFTempVal[i*5 +1] += cblas_sdot(2*uKh_FFcnt[i], fFslip, 1, fKh_FFvaldip[i], 1);
        }   }
        
        uTemp0   = 0u;
        fCmbSlp2 = 0.0f;
        for (i = 0u; i < iFOFFSET[iRANK]; i++) 
        {   fTemp0   = fFTempVal[i*5 +2];           fTemp1   = fFTempVal[i*5 +3];
            fTemp2   = sqrtf(fFEvent[i*17 +13]*fFEvent[i*17 +13] +fFEvent[i*17 +14]*fFEvent[i*17 +14]);
            
            uTemp0   = (fTemp2 <= 0.0f)*uTemp0   + (fTemp2 > 0.0f)*(uTemp0 +1u);
            fCmbSlp2 = (fTemp2 <= 0.0f)*fCmbSlp2 + (fTemp2 > 0.0f)*(fCmbSlp2 + sqrtf(fTemp0*fTemp0 +fTemp1*fTemp1));
            
            fFRef[i*14 +12]   = sqrtf(fTemp0*fTemp0 +fTemp1*fTemp1);
            fFEvent[i*17 +13] = 0.0f;
            fFEvent[i*17 +14] = 0.0f;
        } 
        
        MPI_Allreduce(MPI_IN_PLACE, &fCmbSlp2, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, &uTemp0,   1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
        fCmbSlp2 = (uTemp0 > 0u) ? (fCmbSlp2/(float)uTemp0) : 0.0f;
        
        if ((fCmbSlp1 > FLT_EPSILON) && (fCmbSlp2 > FLT_EPSILON))     {   fTemp0 = fCmbSlp1/fCmbSlp2;     }
        else                                                          {   fTemp0 = 1.0f;                  }
        
        for (i = 0u; i < iFOFFSET[iRANK]; i++) 
        {   fFEvent[i*17 +8]  *= fTemp0;        fFEvent[i*17 +9]  *= fTemp0;
            fFTempVal[i*5 +2] *= fTemp0;        fFTempVal[i*5 +3] *= fTemp0;
            fFRef[i*14 +12]   *= fTemp0;
        }
        if (iRANK == 0)     {   fprintf(stdout,"\nResulting average stressing rate on faults (mm/yr): %5.5f   and  %5.5f,   scalefactor: %5.5f\n",fCmbSlp1*1000.0f, fCmbSlp2*1000.0f, fTemp0);      }
        
        unsigned int *uTempF0= calloc(uFPNum, sizeof *uTempF0);
        unsigned int *uTempF1= calloc(uFPNum, sizeof *uTempF1);
        float *fTempF0  = calloc(uFPNum, sizeof *fTempF0 );
        float *fTempF1  = calloc(uFPNum, sizeof *fTempF1 );
        float *fTempF2  = calloc(uFPNum, sizeof *fTempF2 );
        float *fTempF3  = calloc(uFPNum, sizeof *fTempF3 );
        float *fTempF4  = calloc(uFPNum, sizeof *fTempF4 );
        float *fTempF5  = calloc(uFPNum, sizeof *fTempF5 );
        float *fTempB1  = calloc(uBPNum, sizeof *fTempB1 );
        float *fTempB2  = calloc(uBPNum, sizeof *fTempB2 );
        float *fTempB3  = calloc(uBPNum, sizeof *fTempB3 );
        unsigned int *uTempFL1 = calloc(iFOFFSET[iRANK], sizeof *uTempFL1 );
        float *fTempFL1 = calloc(iFOFFSET[iRANK], sizeof *fTempFL1 );
        float *fTempFL2 = calloc(iFOFFSET[iRANK], sizeof *fTempFL2 );
        float *fTempFL3 = calloc(iFOFFSET[iRANK], sizeof *fTempFL3 );
        float *fTempFL4 = calloc(iFOFFSET[iRANK], sizeof *fTempFL4 );
        float *fTempFL5 = calloc(iFOFFSET[iRANK], sizeof *fTempFL5 );
        float *fTempBL1 = calloc(iBOFFSET[iRANK], sizeof *fTempBL1 );
        float *fTempBL2 = calloc(iBOFFSET[iRANK], sizeof *fTempBL2 );
        float *fTempBL3 = calloc(iBOFFSET[iRANK], sizeof *fTempBL3 );
        strcpy(cPrevStateName, cInputName);      strcat(cPrevStateName,"_");      sprintf(cNameAppend, "%u",uRunNum);     strcat(cPrevStateName,cNameAppend);     strcat(cPrevStateName,"_PreRunData.dat");
        
        
        for (i = 0u; i < iFOFFSET[iRANK]; i++)      
        {   fTempFL1[i] = fFTempVal[i*5 +2];                         fTempFL2[i] = fFTempVal[i*5 +3];
            fTempFL3[i] = fFFric[i*6 +0] * -1.0*fFRef[i*14 +1];      fTempFL4[i] =(fFFric[i*6 +0] - fFFric[i*6 +1])* -1.0*fFRef[i*14 +1];
            fTempFL5[i] = fFEvent[i*17 +7];
            uTempFL1[i] = uFEvent[i*5 +0];
        }
        for (i = 0u; i < iBOFFSET[iRANK]; i++)
        {   fTempBL1[i] = fBTempVal[i*3 +0];                        fTempBL2[i] = fBTempVal[i*3 +1];                fTempBL3[i] = fBTempVal[i*3 +2];
        }
        MPI_Allgatherv(fTempFL1, iFOFFSET[iRANK], MPI_FLOAT, fTempF1, iFOFFSET, iFSTART, MPI_FLOAT, MPI_COMM_WORLD);
        MPI_Allgatherv(fTempFL2, iFOFFSET[iRANK], MPI_FLOAT, fTempF2, iFOFFSET, iFSTART, MPI_FLOAT, MPI_COMM_WORLD);
        MPI_Allgatherv(fTempFL3, iFOFFSET[iRANK], MPI_FLOAT, fTempF3, iFOFFSET, iFSTART, MPI_FLOAT, MPI_COMM_WORLD);
        MPI_Allgatherv(fTempFL4, iFOFFSET[iRANK], MPI_FLOAT, fTempF4, iFOFFSET, iFSTART, MPI_FLOAT, MPI_COMM_WORLD);
        MPI_Allgatherv(fTempFL5, iFOFFSET[iRANK], MPI_FLOAT, fTempF5, iFOFFSET, iFSTART, MPI_FLOAT, MPI_COMM_WORLD);
        MPI_Allgatherv(uTempFL1, iFOFFSET[iRANK], MPI_UNSIGNED,uTempF1, iFOFFSET, iFSTART, MPI_UNSIGNED, MPI_COMM_WORLD);
        MPI_Allgatherv(fTempBL1, iBOFFSET[iRANK], MPI_FLOAT, fTempB1, iBOFFSET, iBSTART, MPI_FLOAT, MPI_COMM_WORLD);
        MPI_Allgatherv(fTempBL2, iBOFFSET[iRANK], MPI_FLOAT, fTempB2, iBOFFSET, iBSTART, MPI_FLOAT, MPI_COMM_WORLD);
        MPI_Allgatherv(fTempBL3, iBOFFSET[iRANK], MPI_FLOAT, fTempB3, iBOFFSET, iBSTART, MPI_FLOAT, MPI_COMM_WORLD);       
        
        if (iRANK == 0)       
        {   
            if ((fpPrePost = fopen(cPrevStateName,"wb")) == NULL)        {   printf("Error -cant open %s  PrePostFile...\n",cPrevStateName);      exit(10);     }
            
            unsigned int *uTempB = calloc(uBPNum, sizeof *uTempB );
            float *fTempFV = calloc(uFVNum, sizeof *fTempFV );
            float *fTempBV = calloc(uBVNum, sizeof *fTempBV );
            
            fwrite( &uFPNum,  sizeof(unsigned int), 1, fpPrePost);
            fwrite( &uFVNum,  sizeof(unsigned int), 1, fpPrePost);
            fwrite( &uBPNum,  sizeof(unsigned int), 1, fpPrePost);
            fwrite( &uBVNum,  sizeof(unsigned int), 1, fpPrePost);
            for (i = 0u; i < uFPNum; i++)   {   uTempF0[i] = uF_temp[i*7 +0];   }     fwrite( uTempF0,     sizeof(unsigned int),  uFPNum, fpPrePost);
            for (i = 0u; i < uFPNum; i++)   {   uTempF0[i] = uF_temp[i*7 +1];   }     fwrite( uTempF0,     sizeof(unsigned int),  uFPNum, fpPrePost);
            for (i = 0u; i < uFPNum; i++)   {   uTempF0[i] = uF_temp[i*7 +2];   }     fwrite( uTempF0,     sizeof(unsigned int),  uFPNum, fpPrePost);
            for (i = 0u; i < uFPNum; i++)   {   uTempF0[i] = uF_temp[i*7 +3];   }     fwrite( uTempF0,     sizeof(unsigned int),  uFPNum, fpPrePost);
            for (i = 0u; i < uFPNum; i++)   {   uTempF0[i] = uF_temp[i*7 +4];   }     fwrite( uTempF0,     sizeof(unsigned int),  uFPNum, fpPrePost);
            
            for (i = 0u; i < uFVNum; i++)   {   fTempFV[i]= fFV_temp[i*4 +0];  }      fwrite( fTempFV,     sizeof(float),         uFVNum, fpPrePost);
            for (i = 0u; i < uFVNum; i++)   {   fTempFV[i]= fFV_temp[i*4 +1];  }      fwrite( fTempFV,     sizeof(float),         uFVNum, fpPrePost);
            for (i = 0u; i < uFVNum; i++)   {   fTempFV[i]= fFV_temp[i*4 +2];  }      fwrite( fTempFV,     sizeof(float),         uFVNum, fpPrePost);
            for (i = 0u; i < uFVNum; i++)   {   fTempFV[i]= fFV_temp[i*4 +3];  }      fwrite( fTempFV,     sizeof(float),         uFVNum, fpPrePost);
            
            for (i = 0u; i < uFPNum; i++)   {   fTempF0[i] = fF_temp[i*16 +0]; }      fwrite( fTempF0,     sizeof(float),         uFPNum, fpPrePost);
            for (i = 0u; i < uFPNum; i++)   {   fTempF0[i] = fF_temp[i*16 +1]; }      fwrite( fTempF0,     sizeof(float),         uFPNum, fpPrePost);
            for (i = 0u; i < uFPNum; i++)   {   fTempF0[i] = fF_temp[i*16 +2]; }      fwrite( fTempF0,     sizeof(float),         uFPNum, fpPrePost);
            
            for (i = 0u; i < uBPNum; i++)   {   uTempB[i] = uB_temp[i*4 +0];   }      fwrite( uTempB,      sizeof(unsigned int),  uBPNum, fpPrePost);
            for (i = 0u; i < uBPNum; i++)   {   uTempB[i] = uB_temp[i*4 +1];   }      fwrite( uTempB,      sizeof(unsigned int),  uBPNum, fpPrePost);
            for (i = 0u; i < uBPNum; i++)   {   uTempB[i] = uB_temp[i*4 +2];   }      fwrite( uTempB,      sizeof(unsigned int),  uBPNum, fpPrePost);
           
            for (i = 0u; i < uBVNum; i++)   {   fTempBV[i]= fBV_temp[i*3 +0];  }      fwrite( fTempBV,     sizeof(float),         uBVNum, fpPrePost);
            for (i = 0u; i < uBVNum; i++)   {   fTempBV[i]= fBV_temp[i*3 +1];  }      fwrite( fTempBV,     sizeof(float),         uBVNum, fpPrePost);
            for (i = 0u; i < uBVNum; i++)   {   fTempBV[i]= fBV_temp[i*3 +2];  }      fwrite( fTempBV,     sizeof(float),         uBVNum, fpPrePost);
              
            for (i = 0u; i < uFPNum; i++)   {   uTempF0[i] = uF_temp[i*7 +6];   }     fwrite( uTempF0,     sizeof(unsigned int),  uFPNum, fpPrePost);
            
            fwrite(fTempF1, sizeof(float), uFPNum, fpPrePost);                     
            fwrite(fTempF2, sizeof(float), uFPNum, fpPrePost);
            fwrite(fTempF3, sizeof(float), uFPNum, fpPrePost); 
            fwrite(fTempF4, sizeof(float), uFPNum, fpPrePost); 
            fwrite(fTempF5, sizeof(float), uFPNum, fpPrePost); 
            fwrite(uTempF1, sizeof(unsigned int), uFPNum, fpPrePost);
            fwrite(fTempB1, sizeof(float), uBPNum, fpPrePost); 
            fwrite(fTempB2, sizeof(float), uBPNum, fpPrePost); 
            fwrite(fTempB3, sizeof(float), uBPNum, fpPrePost); 
            
            fclose(fpPrePost); 
        }
        MPI_Barrier( MPI_COMM_WORLD );
    } 
    
    
    MPI_Barrier( MPI_COMM_WORLD );
    strcpy(cOutputName,    cInputName);      sprintf(cNameAppend, "_%u_RAWL3Catalog.dat",uRunNum);        strcat(cOutputName,    cNameAppend);
    strcpy(cPrevStateName, cInputName);      sprintf(cNameAppend, "_%u_PostRunState.dat",uRunNum);        strcat(cPrevStateName, cNameAppend);
    if      (uCatType == 1)     {   sprintf(cNameAppend, "_%u_L1Catalog.dat",uRunNum);                      }
    else if (uCatType == 2)     {   sprintf(cNameAppend, "_%u_L2Catalog.dat",uRunNum);                      }
    else if (uCatType == 3)     {   sprintf(cNameAppend, "_%u_L3Catalog.dat",uRunNum);                      }
    else                        {   fprintf(stdout,"Catalog level incorrect -make it 1/2/3");   exit(111);  }
    strcpy(cOutputNameCUT, cInputName); strcat(cOutputNameCUT, cNameAppend);
    
    if (uContPrevRun == 1u)
    {   
        unsigned int uDummy;
        unsigned int *uFvTemp   = calloc( iFOFFSET[iRANK],  sizeof *uFvTemp ); 
        float        *fFvTemp   = calloc( iFOFFSET[iRANK],  sizeof *fFvTemp ); 
        float        *fBvTemp   = calloc( iBOFFSET[iRANK],  sizeof *fBvTemp ); 
        
        MPI_File_open(MPI_COMM_WORLD, cPrevStateName, MPI_MODE_RDONLY, MPI_INFO_NULL, &fp_PREVRUN);
        
        MPI_File_read(fp_PREVRUN, &dAddedTime, 1, MPI_DOUBLE,             &STATUS);
        MPI_File_read(fp_PREVRUN, &uEQcntr,    1, MPI_UNSIGNED,           &STATUS);
        MPI_File_read(fp_PREVRUN, &uDummy,     1, MPI_UNSIGNED,           &STATUS);
        MPI_File_read(fp_PREVRUN, &fPSeisTime, 1, MPI_FLOAT,              &STATUS);
        
        OFFSETall = sizeof(double) +2*sizeof(unsigned int) +sizeof(float);
        
        MPI_File_read_at(fp_PREVRUN,   (OFFSETall + iFSTART[iRANK]*sizeof(unsigned int)), uFvTemp, iFOFFSET[iRANK], MPI_UNSIGNED, &STATUS); 
        OFFSETall += (uFPNum*sizeof(unsigned int));    for (i = 0u; i < iFOFFSET[iRANK]; i++)      {   uFEvent[i*5 +0] = uFvTemp[i];        }
        
        MPI_File_read_at(fp_PREVRUN,   (OFFSETall + iFSTART[iRANK]*sizeof(float)), fFvTemp, iFOFFSET[iRANK], MPI_FLOAT, &STATUS); 
        OFFSETall += (uFPNum*sizeof(float));           for (i = 0u; i < iFOFFSET[iRANK]; i++)      {   fFEvent[i*17 +0] = fFvTemp[i];        }
        
        MPI_File_read_at(fp_PREVRUN,   (OFFSETall + iFSTART[iRANK]*sizeof(float)), fFvTemp, iFOFFSET[iRANK], MPI_FLOAT, &STATUS); 
        OFFSETall += (uFPNum*sizeof(float));           for (i = 0u; i < iFOFFSET[iRANK]; i++)      {   fFEvent[i*17 +1] = fFvTemp[i];        }
        
        MPI_File_read_at(fp_PREVRUN,   (OFFSETall + iFSTART[iRANK]*sizeof(float)), fFvTemp, iFOFFSET[iRANK], MPI_FLOAT, &STATUS); 
        OFFSETall += (uFPNum*sizeof(float));           for (i = 0u; i < iFOFFSET[iRANK]; i++)      {   fFEvent[i*17 +2] = fFvTemp[i];        }
        
        MPI_File_read_at(fp_PREVRUN,   (OFFSETall + iFSTART[iRANK]*sizeof(float)), fFvTemp, iFOFFSET[iRANK], MPI_FLOAT, &STATUS); 
        OFFSETall += (uFPNum*sizeof(float));           for (i = 0u; i < iFOFFSET[iRANK]; i++)      {   fFEvent[i*17 +3] = fFvTemp[i];        }
        
        MPI_File_read_at(fp_PREVRUN,   (OFFSETall + iFSTART[iRANK]*sizeof(float)), fFvTemp, iFOFFSET[iRANK], MPI_FLOAT, &STATUS); 
        OFFSETall += (uFPNum*sizeof(float));           for (i = 0u; i < iFOFFSET[iRANK]; i++)      {   fFEvent[i*17 +7] = fFvTemp[i];        }
        
        MPI_File_read_at(fp_PREVRUN,   (OFFSETall + iFSTART[iRANK]*sizeof(float)), fFvTemp, iFOFFSET[iRANK], MPI_FLOAT, &STATUS); 
        OFFSETall += (uFPNum*sizeof(float));           for (i = 0u; i < iFOFFSET[iRANK]; i++)      {   fFEvent[i*17 +10] = fFvTemp[i];       }
        
        MPI_File_read_at(fp_PREVRUN,   (OFFSETall + iFSTART[iRANK]*sizeof(float)), fFvTemp, iFOFFSET[iRANK], MPI_FLOAT, &STATUS); 
        OFFSETall += (uFPNum*sizeof(float));           for (i = 0u; i < iFOFFSET[iRANK]; i++)      {   fFFric[i*6 +0] = fFvTemp[i];          }
        
        MPI_File_read_at(fp_PREVRUN,   (OFFSETall + iFSTART[iRANK]*sizeof(float)), fFvTemp, iFOFFSET[iRANK], MPI_FLOAT, &STATUS); 
        OFFSETall += (uFPNum*sizeof(float));           for (i = 0u; i < iFOFFSET[iRANK]; i++)      {   fFFric[i*6 +1] = fFvTemp[i];          }
        
        MPI_File_read_at(fp_PREVRUN,   (OFFSETall + iFSTART[iRANK]*sizeof(float)), fFvTemp, iFOFFSET[iRANK], MPI_FLOAT, &STATUS); 
        OFFSETall += (uFPNum*sizeof(float));           for (i = 0u; i < iFOFFSET[iRANK]; i++)      {   fFFric[i*6 +2] = fFvTemp[i];          }
        
        MPI_File_read_at(fp_PREVRUN,   (OFFSETall + iFSTART[iRANK]*sizeof(float)), fFvTemp, iFOFFSET[iRANK], MPI_FLOAT, &STATUS); 
        OFFSETall += (uFPNum*sizeof(float));           for (i = 0u; i < iFOFFSET[iRANK]; i++)      {   fFFric[i*6 +4] = fFvTemp[i];          }
        
        MPI_File_read_at(fp_PREVRUN,   (OFFSETall + iFSTART[iRANK]*sizeof(float)), fFvTemp, iFOFFSET[iRANK], MPI_FLOAT, &STATUS); 
        OFFSETall += (uFPNum*sizeof(float));           for (i = 0u; i < iFOFFSET[iRANK]; i++)      {   fFRef[i*14 +12] = fFvTemp[i];          }
        
        MPI_File_read_at(fp_PREVRUN,   (OFFSETall + iFSTART[iRANK]*sizeof(float)), fFvTemp, iFOFFSET[iRANK], MPI_FLOAT, &STATUS); 
        OFFSETall += (uFPNum*sizeof(float));           for (i = 0u; i < iFOFFSET[iRANK]; i++)      {   fFRef[i*14 +13] = fFvTemp[i];          }
        
        MPI_File_read_at(fp_PREVRUN,   (OFFSETall + iBSTART[iRANK]*sizeof(float)), fBvTemp, iBOFFSET[iRANK], MPI_FLOAT, &STATUS); 
        OFFSETall += (uBPNum*sizeof(float));           for (i = 0u; i < iBOFFSET[iRANK]; i++)      {   fBEvent[i*9 +0] = fBvTemp[i];         }
        
        MPI_File_read_at(fp_PREVRUN,   (OFFSETall + iBSTART[iRANK]*sizeof(float)), fBvTemp, iBOFFSET[iRANK], MPI_FLOAT, &STATUS); 
        OFFSETall += (uBPNum*sizeof(float));           for (i = 0u; i < iBOFFSET[iRANK]; i++)      {   fBEvent[i*9 +1] = fBvTemp[i];         }
        
        MPI_File_read_at(fp_PREVRUN,   (OFFSETall + iBSTART[iRANK]*sizeof(float)), fBvTemp, iBOFFSET[iRANK], MPI_FLOAT, &STATUS); 
        OFFSETall += (uBPNum*sizeof(float));           for (i = 0u; i < iBOFFSET[iRANK]; i++)      {   fBEvent[i*9 +2] = fBvTemp[i];         }
        
        MPI_File_read_at(fp_PREVRUN,   (OFFSETall + iBSTART[iRANK]*sizeof(float)), fBvTemp, iBOFFSET[iRANK], MPI_FLOAT, &STATUS); 
        OFFSETall += (uBPNum*sizeof(float));           for (i = 0u; i < iBOFFSET[iRANK]; i++)      {   fBEvent[i*9 +7] = fBvTemp[i];         }
        
        MPI_File_read_at(fp_PREVRUN,   (OFFSETall + iBSTART[iRANK]*sizeof(float)), fBvTemp, iBOFFSET[iRANK], MPI_FLOAT, &STATUS); 
        OFFSETall += (uBPNum*sizeof(float));           for (i = 0u; i < iBOFFSET[iRANK]; i++)      {   fBEvent[i*9 +8] = fBvTemp[i];         }
        
        MPI_File_close(&fp_PREVRUN);
        MPI_Barrier( MPI_COMM_WORLD );
        
        if (iRANK == 0)     {   if ((fp_CATALOG = fopen(cOutputName,"rb+"))     == NULL)     {   printf("Error -cant open %s RawCatalogFile...\n",cOutputName);      exit(10);     }           }
       
    }
    else
    {   
        if (iRANK == 0)
        {   if ((fp_CATALOG = fopen(cOutputName,"wb+"))     == NULL)     {   printf("Error -cant open %s  RawCatalogFile...\n",cOutputName);      exit(10);     }
            unsigned int uDummy = 3u;
            unsigned int *uTempF = calloc(uFPNum, sizeof *uTempF );
            float *fTempF  = calloc(uFPNum, sizeof *fTempF );
            float *fTempFV = calloc(uFVNum, sizeof *fTempFV );

            fwrite(&uEQcntr,  sizeof(unsigned int), 1, fp_CATALOG);
            fwrite(&uDummy,   sizeof(unsigned int), 1, fp_CATALOG);
            fwrite(&fModPara[2],     sizeof(float), 1, fp_CATALOG);
            fwrite(&fModPara[7],     sizeof(float), 1, fp_CATALOG);
            fwrite(&uFPNum,   sizeof(unsigned int), 1, fp_CATALOG);
            fwrite(&uFVNum,   sizeof(unsigned int), 1, fp_CATALOG);
        
            for (i = 0u; i < uFPNum; i++)   {   uTempF[i] = uF_temp[i*7 +0];   }     fwrite(uTempF, sizeof(unsigned int), uFPNum, fp_CATALOG);
            for (i = 0u; i < uFPNum; i++)   {   uTempF[i] = uF_temp[i*7 +1];   }     fwrite(uTempF, sizeof(unsigned int), uFPNum, fp_CATALOG);
            for (i = 0u; i < uFPNum; i++)   {   uTempF[i] = uF_temp[i*7 +2];   }     fwrite(uTempF, sizeof(unsigned int), uFPNum, fp_CATALOG);
            for (i = 0u; i < uFPNum; i++)   {   uTempF[i] = uF_temp[i*7 +3];   }     fwrite(uTempF, sizeof(unsigned int), uFPNum, fp_CATALOG);
            for (i = 0u; i < uFPNum; i++)   {   uTempF[i] = uF_temp[i*7 +4];   }     fwrite(uTempF, sizeof(unsigned int), uFPNum, fp_CATALOG);
            
            for (i = 0u; i < uFVNum; i++)   {   fTempFV[i]= fFV_temp[i*4 +0];  }     fwrite(fTempFV, sizeof(float), uFVNum, fp_CATALOG); 
            for (i = 0u; i < uFVNum; i++)   {   fTempFV[i]= fFV_temp[i*4 +1];  }     fwrite(fTempFV, sizeof(float), uFVNum, fp_CATALOG);
            for (i = 0u; i < uFVNum; i++)   {   fTempFV[i]= fFV_temp[i*4 +2];  }     fwrite(fTempFV, sizeof(float), uFVNum, fp_CATALOG);
            for (i = 0u; i < uFVNum; i++)   {   fTempFV[i]= fFV_temp[i*4 +3];  }     fwrite(fTempFV, sizeof(float), uFVNum, fp_CATALOG);
            
            for (i = 0u; i < uFPNum; i++)   {   fTempF[i] = fF_temp[i*16 +0]; }      fwrite(fTempF,  sizeof(float), uFPNum, fp_CATALOG);
            for (i = 0u; i < uFPNum; i++)   {   fTempF[i] = fF_temp[i*16 +1]; }      fwrite(fTempF,  sizeof(float), uFPNum, fp_CATALOG);
            for (i = 0u; i < uFPNum; i++)   {   fTempF[i] = fF_temp[i*16 +2]; }      fwrite(fTempF,  sizeof(float), uFPNum, fp_CATALOG);
    }   }
    
    free(uF_temp);          free(fFV_temp);         free(fF_temp);
    free(uB_temp);          free(fBV_temp);         free(fB_temp);
    
    if (uContPrevRun == 0u)
    {   float fTemp0,   fTemp1,   fTemp2;
        for (i = 0u; i < iFOFFSET[iRANK]; i++)
        {   fTemp0   = sqrtf(fFEvent[i*17 +8]*fFEvent[i*17 +8] + fFEvent[i*17 +9]*fFEvent[i*17 +9]); 
            fTemp1   = fFEvent[i*17 +3]*-1.0f*fFEvent[i*17 +2]; 
            fTemp2   = 0.5f*(1.0f - fPreStressFract)*(ran0(&lSeed)*2.0f -1.0f);
            fFEvent[i*17 +0] = ((fTemp1*(fPreStressFract+fTemp2))/fTemp0) * fFEvent[i*17 +8];
            fFEvent[i*17 +1] = ((fTemp1*(fPreStressFract+fTemp2))/fTemp0) * fFEvent[i*17 +9];
    }   } 
    
    
    
    
    unsigned int uTemp0,   uTemp1,   uGlobPos,   uEQstillOn,   uTotlRuptT,   uActElmG,   uActElmL,   uMRFlgth;
    float fNxtLoadStep, fNxtHalfStep,   fTemp0,   fTemp1,   fTemp2,   fTemp3,    fTemp4,   fTemp5,   fTemp6,   fTemp7;
    float fMaxLoadStep      = FloatPow(2.0f, uLoadStep_POW2);
    float fLoadingStepInYrs = fIntSeisLoadStep/365.25f;
    float fEQMagn,   fEQMomRel,   fHypoSlip,   fHypoLoc[3],   fEQMeanVals[7];
    double dTemp0;
    
    fMaxLoadStep *= fLoadingStepInYrs;
    
    long long int *lSTF_offset= calloc(iSIZE, sizeof *lSTF_offset);
    long long int *lSTF_start = calloc(iSIZE, sizeof *lSTF_start);
    unsigned int *uPatchID    = calloc(iFOFFSET[iRANK], sizeof *uPatchID);
    unsigned int *ut0ofPtch   = calloc(iFOFFSET[iRANK], sizeof *ut0ofPtch);
    unsigned int *uStabType   = calloc(iFOFFSET[iRANK], sizeof *uStabType);
    float        *fPtchDTau   = calloc(iFOFFSET[iRANK], sizeof *fPtchDTau);
    float        *fPtchSlpH   = calloc(iFOFFSET[iRANK], sizeof *fPtchSlpH);
    float        *fPtchSlpV   = calloc(iFOFFSET[iRANK], sizeof *fPtchSlpV);
    unsigned int *uTemp2Wrt1  = calloc(uFPNum,   sizeof *uTemp2Wrt1);
    unsigned int *uTemp2Wrt2  = calloc(uFPNum,   sizeof *uTemp2Wrt2);
    unsigned int *uTemp2Wrt3  = calloc(uFPNum,   sizeof *uTemp2Wrt3);
    float        *fTemp2Wrt1  = calloc(uFPNum,   sizeof *fTemp2Wrt1);
    float        *fTemp2Wrt2  = calloc(uFPNum,   sizeof *fTemp2Wrt2);
    float        *fTemp2Wrt3  = calloc(uFPNum,   sizeof *fTemp2Wrt3);
    
    dRecordLength      += dAddedTime;
    double dCurrTime    = dAddedTime;
    double dLastEQtime  = dAddedTime;
    
    unsigned int uPSeisSteps = 0u;
    float fMaxMagnitude      = 0.0f;
    
    for (j = 0u; j < iBOFFSET[iRANK]; j++)      {   fBEvent[j*9 +0] = 0.0f;        fBEvent[j*9 +1] = 0.0f;         fBEvent[j*9 +2] = 0.0f;          }
    
    
    if (iRANK == 0)     {   fprintf(stdout,"Starting the catalog...max load step (years) = %f\n",fMaxLoadStep);           }
    timer0 = clock();
    timerI = clock();
    
    while (dCurrTime <= dRecordLength)
    {   
        fNxtLoadStep = fMaxLoadStep;
        fNxtHalfStep = fNxtLoadStep;
        
        for (k = 0u; k <= (uLoadStep_POW2+1); k++)
        {   uPSeisSteps += (uEQcntr > 0u)*1u + (uEQcntr <= 0u)*0u;
            fHypoSlip  = 0.0f;      fHypoLoc[0] = 0.0f;       fHypoLoc[1] = 0.0f;      fHypoLoc[2] = 0.0f;
            
            memset(uSlipElCnt, 0,  2*sizeof(unsigned int));
            memset(iStartPosF, 0, iSIZE*sizeof(int) );           memset(iOffstPosF, 0, iSIZE*sizeof(int) );
            memset(iStartPosB, 0, iSIZE*sizeof(int) );           memset(iOffstPosB, 0, iSIZE*sizeof(int) );
            
            fTemp0 = MAX(0.0f, expf(-1.0f*(fPSeisTime + fNxtLoadStep)/fAfterSlipTime));
            
            for (i = 0u; i < iFOFFSET[iRANK]; i++)
            {   dTemp0 = (dCurrTime + (double)fNxtLoadStep)*(double)fFRef[i*14 +12]; 
                fFTempVal[i*5 +0] = (fFRef[i*14 +13] <= dTemp0)*(fFEvent[i*17 +0] +fNxtLoadStep*fFEvent[i*17 +8]) + (fFRef[i*14 +13] > dTemp0)*fFEvent[i*17 +0];
                fFTempVal[i*5 +1] = (fFRef[i*14 +13] <= dTemp0)*(fFEvent[i*17 +1] +fNxtLoadStep*fFEvent[i*17 +9]) + (fFRef[i*14 +13] > dTemp0)*fFEvent[i*17 +1];
                fFTempVal[i*5 +2] = fFEvent[i*17 +2];
                fFTempVal[i*5 +3] = fFFric[i*6 +0] + fTemp0*fFFric[i*6 +4];
                fFTempVal[i*5 +4] = 0.0f;
                
                fTemp1 = sqrtf( fFTempVal[i*5 +0]*fFTempVal[i*5 +0] + fFTempVal[i*5 +1]*fFTempVal[i*5 +1] );
                fTemp1 = (uFEvent[i*5 +0] == 1u)*0.0f + (uFEvent[i*5 +0] != 1u)*fTemp1; 
                fTemp3 = fFTempVal[i*5 +3] *-1.0f*fFTempVal[i*5 +2]; 
                
                fTemp5 = MAX(0.0f, (fTemp1 -fTemp3));
                
                if (fTemp5 > MININTSEISSTRESSCHANGE)
                {   fTemp3 = -1.0f*(((fTemp5/fTemp1)*fFTempVal[i*5 +0])/fFEvent[i*17 +4]);
                    fTemp4 = -1.0f*(((fTemp5/fTemp1)*fFTempVal[i*5 +1])/fFEvent[i*17 +5]);
                    fFTempVal[i*5 +4] = sqrtf(fTemp3*fTemp3 +fTemp4*fTemp4);
                    
                    uGlobPos       = i +iFSTART[iRANK];
                    fFslipL[uSlipElCnt[0]*3 +0] = (float)uGlobPos;        fFslipL[uSlipElCnt[0]*3 +1] = fTemp3;           fFslipL[uSlipElCnt[0]*3 +2] = fTemp4;
                    uSlipElCnt[0] += 1u;
                    
            }   }
            
            iOffstPosF[iRANK] = uSlipElCnt[0]*3;
            
            fTemp0 = MAX(0.0f, expf(-1.0f*(fPSeisTime +fNxtLoadStep)/fDeepRelaxTime));
            for (i = 0u; i < iBOFFSET[iRANK]; i++)
            {   fBTempVal[i*3 +0] = fBEvent[i*9 +0];         fBTempVal[i*3 +1] = fBEvent[i*9 +1];         fBTempVal[i*3 +2] = fBEvent[i*9 +2];
                fTemp1 = sqrtf( fBTempVal[i*3 +0]*fBTempVal[i*3 +0] + fBTempVal[i*3 +1]*fBTempVal[i*3 +1] ); 
                fTemp2 =   fabs(fBTempVal[i*3 +2]); 
                fTemp3 = fTemp0*fBEvent[i*9 +7]; 
                fTemp4 = fTemp0*fBEvent[i*9 +8]; 
                
                fTemp5 = MAX(0.0f, (fTemp1 -fTemp3));               fTemp6 = MAX(0.0f, (fTemp2 -fTemp4));
                fTemp7 = sqrtf(fTemp5*fTemp5 + fTemp6*fTemp6);
                
                if (fTemp7 >  MININTSEISSTRESSCHANGE)
                {   
                    fTemp1 = (fTemp1 <= 0.0f)*-1.0f + (fTemp1 > 0.0f)*fTemp1;
                    fTemp2 = (fTemp2 <= 0.0f)*-1.0f + (fTemp2 > 0.0f)*fTemp2;
                    fTemp3 = (fTemp1 <= 0.0f)*0.0f  + (fTemp1 > 0.0f)*(-1.0f*((fTemp5/fTemp1)*fBTempVal[i*3 +0])/fBEvent[i*9 +4]);
                    fTemp4 = (fTemp1 <= 0.0f)*0.0f  + (fTemp1 > 0.0f)*(-1.0f*((fTemp5/fTemp1)*fBTempVal[i*3 +1])/fBEvent[i*9 +5]);
                    fTemp5 = (fTemp2 <= 0.0f)*0.0f  + (fTemp2 > 0.0f)*(-1.0f*((fTemp6/fTemp2)*fBTempVal[i*3 +2])/fBEvent[i*9 +6]);
                    
                    uGlobPos      = i +iBSTART[iRANK];
                    fBslipL[uSlipElCnt[1]*4 +0] = (float)uGlobPos;        fBslipL[uSlipElCnt[1]*4 +1] = fTemp3;           fBslipL[uSlipElCnt[1]*4 +2] = fTemp4;         fBslipL[uSlipElCnt[1]*4 +3] = fTemp5;
                    uSlipElCnt[1] += 1u;
                    
            }   }
            
            iOffstPosB[iRANK] = uSlipElCnt[1]*4;
            
            MPI_Allreduce(MPI_IN_PLACE, uSlipElCnt, 2, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, iOffstPosF, iSIZE, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, iOffstPosB, iSIZE, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
            
            for (i = 1; i < iSIZE; i++)         {    iStartPosF[i]  = iStartPosF[i-1] +iOffstPosF[i-1];         iStartPosB[i]  = iStartPosB[i-1] +iOffstPosB[i-1];              }
            
            MPI_Allgatherv(fFslipL, iOffstPosF[iRANK], MPI_FLOAT, fFslipG, iOffstPosF, iStartPosF, MPI_FLOAT, MPI_COMM_WORLD);
            MPI_Allgatherv(fBslipL, iOffstPosB[iRANK], MPI_FLOAT, fBslipG, iOffstPosB, iStartPosB, MPI_FLOAT, MPI_COMM_WORLD);
            
            uEQstillOn = 0u;
            for (i = 0u; i < iFOFFSET[iRANK]; i++)
            {   
                if ( uSlipElCnt[1] < uKh_FBcnt[i])
                {   for (j = 0u; j < uSlipElCnt[1]; j++) 
                    {   uTemp0 = (unsigned int)fBslipG[j*4 +0]; 
                        uTemp1 = uKh_FBps2[i][uTemp0]; 
                        fTemp0 = fBslipG[j*4 +1];
                        fTemp1 = fBslipG[j*4 +2];
                        fTemp2 = fBslipG[j*4 +3];
                        
                        fFTempVal[i*5 +0] += (fTemp0*fKh_FBvalstk[i][uTemp1*3 +0] + fTemp1*fKh_FBvalstk[i][uTemp1*3 +1] + fTemp2*fKh_FBvalstk[i][uTemp1*3 +2]);
                        fFTempVal[i*5 +1] += (fTemp0*fKh_FBvaldip[i][uTemp1*3 +0] + fTemp1*fKh_FBvaldip[i][uTemp1*3 +1] + fTemp2*fKh_FBvaldip[i][uTemp1*3 +2]);
                }   }
                else
                {   memset(fBslip, 0, 3*uBPNum*sizeof(float) );
                    for (j = 0u; j < uSlipElCnt[1]; j++) 
                    {   uTemp0 = (unsigned int)fBslipG[j*4 +0]; 
                        uTemp1 = uKh_FBps2[i][uTemp0]; 
                        fBslip[uTemp1*3 +0] += fBslipG[j*4 +1];
                        fBslip[uTemp1*3 +1] += fBslipG[j*4 +2];
                        fBslip[uTemp1*3 +2] += fBslipG[j*4 +3];
                    }
                    fFTempVal[i*5 +0] += cblas_sdot(3*uKh_FBcnt[i], fBslip, 1, fKh_FBvalstk[i], 1);
                    fFTempVal[i*5 +1] += cblas_sdot(3*uKh_FBcnt[i], fBslip, 1, fKh_FBvaldip[i], 1);
                }
                
                if  (uSlipElCnt[0] < uKh_FFcnt[i])
                {   for (j = 0u; j < uSlipElCnt[0]; j++) 
                    {   uTemp0 = (unsigned int)fFslipG[j*3 +0]; 
                        uTemp1 = uKh_FFps2[i][uTemp0]; 
                        fTemp0 = fFslipG[j*3 +1];
                        fTemp1 = fFslipG[j*3 +2];
                        
                        fFTempVal[i*5 +0] += (fTemp0*fKh_FFvalstk[i][uTemp1*2 +0] + fTemp1*fKh_FFvalstk[i][uTemp1*2 +1]);
                        fFTempVal[i*5 +1] += (fTemp0*fKh_FFvaldip[i][uTemp1*2 +0] + fTemp1*fKh_FFvaldip[i][uTemp1*2 +1]);
                }   }
                else
                {   memset(fFslip, 0, 2*uFPNum*sizeof(float) );
                    for (j = 0u; j < uSlipElCnt[0]; j++) 
                    {   uTemp0 = (unsigned int)fFslipG[j*3 +0]; 
                        uTemp1 = uKh_FFps2[i][uTemp0]; 
                        fFslip[uTemp1*2 +0] += fFslipG[j*3 +1];
                        fFslip[uTemp1*2 +1] += fFslipG[j*3 +2];
                    }
                    fFTempVal[i*5 +0] += cblas_sdot(2*uKh_FFcnt[i], fFslip, 1, fKh_FFvalstk[i], 1);
                    fFTempVal[i*5 +1] += cblas_sdot(2*uKh_FFcnt[i], fFslip, 1, fKh_FFvaldip[i], 1);
                }
                
                fTemp1 = sqrtf( fFTempVal[i*5 +0]*fFTempVal[i*5 +0] + fFTempVal[i*5 +1]*fFTempVal[i*5 +1] );
                fTemp1 = (uFEvent[i*5 +0] == 1u)*fTemp1 + (uFEvent[i*5 +0] != 1u)*0.0f; 
                fTemp2 =    fFFric[i*6 +2] *-1.0f*fFTempVal[i*5 +2]; 
                fTemp3 = fFTempVal[i*5 +3] *-1.0f*fFTempVal[i*5 +2]; 
                fTemp5 = MAX(0.0f,(fTemp1 -fTemp3)) /fFRef[i*14 +2]; 
                
                uTemp0 = (fTemp5 > fModPara[10])*1u + (fTemp5 <= fModPara[10])*0u;
                uTemp0 = ((fTemp1-fTemp2) > fFEvent[i*17 +10])*uTemp0 + ((fTemp1-fTemp2) <= fFEvent[i*17 +10])*0u;
                
                uEQstillOn += uTemp0;
                fTemp5 = (uTemp0 == 1u)*fTemp5 + (uTemp0 != 1u)*0.0f; 
                
                if (fTemp5 > fHypoSlip)         {   fHypoSlip  = fTemp5;        fHypoLoc[0] = fFRef[i*14 + 9];       fHypoLoc[1] = fFRef[i*14 +10];      fHypoLoc[2] = fFRef[i*14 +11];         }
                
            }
            
            MPI_Allreduce(MPI_IN_PLACE, &uEQstillOn, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
            
            fNxtHalfStep /= 2.0f;
            if (uEQstillOn == 0u)
            {   if (k == 0u)                                    {                                               break;   }
                else if ((k > 0u) && (k < uLoadStep_POW2))      {   fNxtLoadStep += fNxtHalfStep;                        }
                else if (k == uLoadStep_POW2)                   {   fNxtLoadStep +=(2.0f*fNxtHalfStep);                  }
            }
            else
            {   if      (k < uLoadStep_POW2)                    {   fNxtLoadStep -= fNxtHalfStep;                        }
                else if (k == uLoadStep_POW2)                   {                                               break;   }
            }
            
        } 
        
        for (i = 0u; i < iBOFFSET[iRANK]; i++)
        {   
            if (uSlipElCnt[1] < uKh_BBcnt[i])
            {   for (j = 0u; j < uSlipElCnt[1]; j++) 
                {   uTemp0 = (unsigned int)fBslipG[j*4 +0]; 
                    uTemp1 = uKh_BBps2[i][uTemp0]; 
                    fTemp0 = fBslipG[j*4 +1];
                    fTemp1 = fBslipG[j*4 +2];
                    fTemp2 = fBslipG[j*4 +3];
                    
                    fBTempVal[i*3 +0] += (fTemp0*fKh_BBvalstk[i][uTemp1*3 +0] + fTemp1*fKh_BBvalstk[i][uTemp1*3 +1] + fTemp2*fKh_BBvalstk[i][uTemp1*3 +2]);
                    fBTempVal[i*3 +1] += (fTemp0*fKh_BBvaldip[i][uTemp1*3 +0] + fTemp1*fKh_BBvaldip[i][uTemp1*3 +1] + fTemp2*fKh_BBvaldip[i][uTemp1*3 +2]);
                    fBTempVal[i*3 +2] += (fTemp0*fKh_BBvalnrm[i][uTemp1*3 +0] + fTemp1*fKh_BBvalnrm[i][uTemp1*3 +1] + fTemp2*fKh_BBvalnrm[i][uTemp1*3 +2]);
            }   }
            else
            {   memset(fBslip, 0, 3*uBPNum*sizeof(float) );
                for (j = 0u; j < uSlipElCnt[1]; j++) 
                {   uTemp0 = (unsigned int)fBslipG[j*4 +0]; 
                    uTemp1 = uKh_BBps2[i][uTemp0]; 
                    fBslip[uTemp1*3 +0] += fBslipG[j*4 +1];
                    fBslip[uTemp1*3 +1] += fBslipG[j*4 +2];
                    fBslip[uTemp1*3 +2] += fBslipG[j*4 +3];
                }
                fBTempVal[i*3 +0] += cblas_sdot(3*uKh_BBcnt[i], fBslip, 1, fKh_BBvalstk[i], 1);
                fBTempVal[i*3 +1] += cblas_sdot(3*uKh_BBcnt[i], fBslip, 1, fKh_BBvaldip[i], 1);
                fBTempVal[i*3 +2] += cblas_sdot(3*uKh_BBcnt[i], fBslip, 1, fKh_BBvalnrm[i], 1);
            }
            
            if (uSlipElCnt[0] < uKh_BFcnt[i])
            {   for (j = 0u; j < uSlipElCnt[0]; j++) 
                {   uTemp0 = (unsigned int)fFslipG[j*3 +0]; 
                    uTemp1 = uKh_BFps2[i][uTemp0]; 
                    fTemp0 = fFslipG[j*3 +1];
                    fTemp1 = fFslipG[j*3 +2];
                    
                    fBTempVal[i*3 +0] += (fTemp0*fKh_BFvalstk[i][uTemp1*2 +0] + fTemp1*fKh_BFvalstk[i][uTemp1*2 +1]);
                    fBTempVal[i*3 +1] += (fTemp0*fKh_BFvaldip[i][uTemp1*2 +0] + fTemp1*fKh_BFvaldip[i][uTemp1*2 +1]);
                    fBTempVal[i*3 +2] += (fTemp0*fKh_BFvalnrm[i][uTemp1*2 +0] + fTemp1*fKh_BFvalnrm[i][uTemp1*2 +1]);
            }   }
            else
            {   memset(fFslip, 0, 2*uFPNum*sizeof(float) );
                for (j = 0u; j < uSlipElCnt[0]; j++) 
                {   uTemp0 = (unsigned int)fFslipG[j*3 +0]; 
                    uTemp1 = uKh_BFps2[i][uTemp0]; 
                    fFslip[uTemp1*2 +0] += fFslipG[j*3 +1];
                    fFslip[uTemp1*2 +1] += fFslipG[j*3 +2];
                }
               fBTempVal[i*3 +0] += cblas_sdot(2*uKh_BFcnt[i], fFslip, 1, fKh_BFvalstk[i], 1);
               fBTempVal[i*3 +1] += cblas_sdot(2*uKh_BFcnt[i], fFslip, 1, fKh_BFvaldip[i], 1);
               fBTempVal[i*3 +2] += cblas_sdot(2*uKh_BFcnt[i], fFslip, 1, fKh_BFvalnrm[i], 1);
        }   }
        
        fPSeisTime += fNxtLoadStep;
        dCurrTime  += (double)fNxtLoadStep;
        
        for (i = 0u; i < iFOFFSET[iRANK]; i++)
        {   fFEvent[i*17 +0] = fFTempVal[i*5 +0];       fFEvent[i*17 +1] = fFTempVal[i*5 +1];       fFEvent[i*17 +2] = fFTempVal[i*5 +2];
            fFEvent[i*17 +3] = fFTempVal[i*5 +3];       fFRef[i*14 +13] += fFTempVal[i*5 +4]; 
        }
        for (i = 0u; i < iBOFFSET[iRANK]; i++)
        {   fBEvent[i*9 +0] = fBTempVal[i*3 +0];        fBEvent[i*9 +1] = fBTempVal[i*3 +1];        fBEvent[i*9 +2] = fBTempVal[i*3 +2];
        }
        
        
        
        if (uEQstillOn > 0u)
        {   timerI = clock() - timerI;        time_InterSeis += (((double)timerI)/CLOCKS_PER_SEC);
            timerC = clock();
            uEQtickr += 1u;
            
            uMRFlgth  = 0u;             uActElmL  = 0u;         uTotlRuptT = 0u;            fEQMomRel = 0.0f;
            memset(uFActivEl, 0, iFOFFSET[iRANK]*sizeof(unsigned int));
            memset(fMRFvals,  0, MAXMOMRATEFUNCLENGTH*sizeof(float) );
            memset(fSTFslip,  0, 2*iFOFFSET[iRANK]*MAXMOMRATEFUNCLENGTH*sizeof(float));
            
            {   struct { float val;  int rank;  } in[1], out[1];
                in[0].val  = fHypoSlip;         in[0].rank = iRANK;
                
                MPI_Allreduce(in, out, 1, MPI_FLOAT_INT, MPI_MAXLOC, MPI_COMM_WORLD);
                
                if (out[0].rank != 0)
                {   if (iRANK == out[0].rank)   {   MPI_Send(fHypoLoc, 3, MPI_FLOAT,    0,        123, MPI_COMM_WORLD);               }
                    else if (iRANK == 0)        {   MPI_Recv(fHypoLoc, 3, MPI_FLOAT, out[0].rank, 123, MPI_COMM_WORLD, &STATUS);      }
                }
            } 
            
            for (i = 0u; i < iFOFFSET[iRANK]; i++)
            {   
                fFTempVal[i*5 +0] = fFEvent[i*17 +0];       fFTempVal[i*5 +1] = fFEvent[i*17 +1];       fFTempVal[i*5 +2] = fFEvent[i*17 +2];       fFTempVal[i*5 +3] = fFEvent[i*17 +3];
                
                fTemp0 = (fFFric[i*6 +1] *-1.0f *fFEvent[i*17 +2]) + fFRef[i*14 +2]*fFEvent[i*17 +7]; 
                fTemp0/= (-1.0f *fFEvent[i*17 +2]);
                fFEvent[i*17 +3] = MAX(fFEvent[i*17 +3], fTemp0);
                fFFric[i*6 +3]   = fFEvent[i*17 +3]; 
                fFFric[i*6 +5]   = fFEvent[i*17 +3]; 
                
            }
            
            
            while (uEQstillOn > 0u)
            {   uEQstillOn = 0u;
                memset(iStartPosF, 0, iSIZE*sizeof(int) );           memset(iOffstPosF, 0, iSIZE*sizeof(int) );
                
                for (i = 0u; i < iFOFFSET[iRANK]; i++)
                {   fTemp2 = fFFric[i*6 +3] - (fFFric[i*6 +3] - fFFric[i*6 +2]) *(fFEvent[i*17 +12]/fFEvent[i*17 +7]);
                    fTemp2 = MAX(fTemp2, fFFric[i*6 +2]);
                    fTemp3 = fFEvent[i*17 +3] + fHealFact*(fFFric[i*6 +5] - fFEvent[i*17 +3]); 
                    
                    fFEvent[i*17 +3] = (uFEvent[i*5 +1] != 1u)*fFEvent[i*17 +3]       + (uFEvent[i*5 +1] == 1u)*((fFEvent[i*17 +11] == 0.0f)*fTemp3  + (fFEvent[i*17 +11] != 0.0f)*fTemp2           );
                    fFFric[i*6 +3]   = (uFEvent[i*5 +1] != 1u)*fFFric[i*6 +3]         + (uFEvent[i*5 +1] == 1u)*((fFEvent[i*17 +11] == 0.0f)*fTemp3  + (fFEvent[i*17 +11] != 0.0f)*fFFric[i*6 +3]   );
                    fFEvent[i*17 +12]= (fFEvent[i*17 +11] != 0.0f)*fFEvent[i*17 +12]  + (fFEvent[i*17 +11] == 0.0f)*0.0f;
                    fFEvent[i*17 +11]= 0.0f;
                    
                    fTemp0           = sqrtf( fFEvent[i*17 +0]*fFEvent[i*17 +0] + fFEvent[i*17 +1]*fFEvent[i*17 +1] ); 
                    fFEvent[i*17 +2] = MIN(0.0f, fFEvent[i*17 +2]);
                    fTemp3           = INTERNALCOHESION + INTERNALFRICTION*-1.0f*fFEvent[i*17 +2]; 
                    uFEvent[i*5 +2]  = (fTemp0 > fTemp3)*1u + (fTemp0 <= fTemp3)*uFEvent[i*5 +2]; 
                    
                    uTemp0           = (uFEvent[i*5 +1] == 0u)*uFEvent[i*5 +2] + (uFEvent[i*5 +1] != 0u)*0u; 
                    
                    uFActivEl[uActElmL] = (uTemp0 == 1u)*i              + (uTemp0 != 1u)*uFActivEl[uActElmL]; 
                    uActElmL            = (uTemp0 == 1u)*(uActElmL+1u)  + (uTemp0 != 1u)*uActElmL;
                    uFEvent[i*5 +1]     = (uTemp0 == 1u)*2u             + (uTemp0 != 1u)*uFEvent[i*5 +1];
                    uFEvent[i*5 +3]     = (uTemp0 == 1u)*(uTotlRuptT)   + (uTemp0 != 1u)*uFEvent[i*5 +3];
                    uFEvent[i*5 +4]     = (uTemp0 == 1u)*1u             + (uTemp0 != 1u)*uFEvent[i*5 +4];
                    
                    fTemp1   = fFEvent[i*17 +3] *-1.0f*fFEvent[i*17 +2]; 
                    fTemp2   = fFFric[i*6 +2]   *-1.0f*fFEvent[i*17 +2]; 
                    fTemp3   = fTemp0 - fTemp1; 
                    fTemp4   = fTemp0 - fTemp2; 
                    fTemp5   = MAX(fTemp3, 0.0f) / fFRef[i*14 +2]; 
                    
                    uTemp0   = (uFEvent[i*5 +2] == 0u)*1u           + (uFEvent[i*5 +2] != 0u)*0u; 
                    uTemp0   = (fTemp5 >= fModPara[10])*uTemp0      + (fTemp5 < fModPara[10])*0u;
                    uTemp0   = (fTemp4 >= fFEvent[i*17 +10])*uTemp0 + (fTemp4 < fFEvent[i*17 +10])*0u;
                    
                    if (uTemp0 == 1u) 
                    {   
                        uFActivEl[uActElmL] = (uFEvent[i*5 +1] == 0u)*i            + (uFEvent[i*5 +1] != 0u)*uFActivEl[uActElmL];
                        uActElmL            = (uFEvent[i*5 +1] == 0u)*(uActElmL+1u)+ (uFEvent[i*5 +1] != 0u)*uActElmL;
                        uFEvent[i*5 +3]     = (uFEvent[i*5 +1] == 0u)*(uTotlRuptT) + (uFEvent[i*5 +1] != 0u)*uFEvent[i*5 +3];
                        uFEvent[i*5 +1]     = (uFEvent[i*5 +1] == 0u)*1u           + (uFEvent[i*5 +1] != 0u)*uFEvent[i*5 +1];
                        
                        uMRFlgth          = uTotlRuptT +1u;
                        uFEvent[i*5 +4]   = uTotlRuptT +1u - uFEvent[i*5 +3];
                        fFEvent[i*17 +15] = fabs( fModPara[9] *(fTemp3/fModPara[2]) *fModPara[7] );
                        
                        fTemp5 = -1.0f*((fTemp3/fTemp0) *fFEvent[i*17 +0]) / fFEvent[i*17 +4];
                        fTemp6 = -1.0f*((fTemp3/fTemp0) *fFEvent[i*17 +1]) / fFEvent[i*17 +5];
                        fTemp7 = fFEvent[i*17 +15]/sqrtf(fTemp5*fTemp5 +fTemp6*fTemp6);
                        fTemp7 = MIN(fTemp7, 1.0f);
                        fTemp5*= fTemp7;
                        fTemp6*= fTemp7;
                        
                        uGlobPos    = i +iFSTART[iRANK];
                        fFslipL[uEQstillOn*3 +0] = (float)uGlobPos;        fFslipL[uEQstillOn*3 +1] = fTemp5;               fFslipL[uEQstillOn*3 +2] = fTemp6;
                        uEQstillOn += 1u;
                        
                        uTemp0 = 2u*i*MAXMOMRATEFUNCLENGTH + 2u*(uTotlRuptT - uFEvent[i*5 +3]);
                        fSTFslip[uTemp0 +0] = fTemp5;                          fSTFslip[uTemp0 +1] = fTemp6;
                        
                        fFEvent[i*17 +11]   = sqrtf(fTemp5*fTemp5 +fTemp6*fTemp6);
                        fFEvent[i*17 +12]  += fFEvent[i*17 +11];
                        fFEvent[i*17 +13]  += fTemp5;
                        fFEvent[i*17 +14]  += fTemp6;
                        fEQMomRel          += (fFEvent[i*17 +11]*fFRef[i*14 +0]);
                        fMRFvals[uMRFlgth] += (fFEvent[i*17 +11]*fFRef[i*14 +0]);
                        
                }   }
                iOffstPosF[iRANK] = uEQstillOn*3;
                
                MPI_Allreduce(MPI_IN_PLACE, &uEQstillOn, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD); 
                MPI_Allreduce(MPI_IN_PLACE, iOffstPosF, iSIZE, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
                for (i = 1; i < iSIZE; i++)        {    iStartPosF[i]  = iStartPosF[i-1] +iOffstPosF[i-1];    } 
                MPI_Allgatherv(fFslipL, iOffstPosF[iRANK], MPI_FLOAT, fFslipG, iOffstPosF, iStartPosF, MPI_FLOAT, MPI_COMM_WORLD);
                
                for (i = 0u; i < iFOFFSET[iRANK]; i++)
                {   if (uEQstillOn < uKh_FFcnt[i])
                    {   for (j = 0u; j < uEQstillOn; j++)
                        {   uTemp0 = (unsigned int)fFslipG[j*3 +0]; 
                            uTemp1 = uKh_FFps2[i][uTemp0]; 
                            fTemp0 = fFslipG[j*3 +1]; 
                            fTemp1 = fFslipG[j*3 +2]; 
                            
                            fFEvent[i*17 +0] += (fTemp0*fKh_FFvalstk[i][uTemp1*2 +0] + fTemp1*fKh_FFvalstk[i][uTemp1*2 +1]);
                            fFEvent[i*17 +1] += (fTemp0*fKh_FFvaldip[i][uTemp1*2 +0] + fTemp1*fKh_FFvaldip[i][uTemp1*2 +1]);
                            fFEvent[i*17 +2] += (fTemp0*fKh_FFvalnrm[i][uTemp1*2 +0] + fTemp1*fKh_FFvalnrm[i][uTemp1*2 +1]);
                    }   }
                    else
                    {   memset(fFslip, 0, 2*uFPNum*sizeof(float) );
                        for (j = 0u; j < uEQstillOn; j++) 
                        {   uTemp0 = (unsigned int)fFslipG[j*3 +0]; 
                            uTemp1 = uKh_FFps2[i][uTemp0]; 
                            fFslip[uTemp1*2 +0] += fFslipG[j*3 +1];
                            fFslip[uTemp1*2 +1] += fFslipG[j*3 +2];
                        }
                        fFEvent[i*17 +0] += cblas_sdot(2*uKh_FFcnt[i], fFslip, 1, fKh_FFvalstk[i], 1);
                        fFEvent[i*17 +1] += cblas_sdot(2*uKh_FFcnt[i], fFslip, 1, fKh_FFvaldip[i], 1);
                        fFEvent[i*17 +2] += cblas_sdot(2*uKh_FFcnt[i], fFslip, 1, fKh_FFvalnrm[i], 1);
                }   }
                uEQstillOn = (uTotlRuptT < MAXMOMRATEFUNCLENGTH-1)*uEQstillOn  + (uTotlRuptT >= MAXMOMRATEFUNCLENGTH-1)*0u; 
                uTotlRuptT += 1u;
            }
            
            
            uActElmG = uActElmL;
            MPI_Allreduce(MPI_IN_PLACE, &uActElmG,   1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, &fEQMomRel,  1, MPI_FLOAT,    MPI_SUM, MPI_COMM_WORLD);
            
            fEQMagn = (log10f(fEQMomRel * fModPara[2]) -9.1f)/1.5f;
            
            
            if (uBPNum > 0)
            {   memset(uSlipElCnt, 0,  2*sizeof(unsigned int));
                memset(iStartPosF, 0, iSIZE*sizeof(int) );           memset(iOffstPosF, 0, iSIZE*sizeof(int) );
                for (i = 0u; i < iFOFFSET[iRANK]; i++)
                {   fTemp0 = fFEvent[i*17 +13];
                    fTemp1 = fFEvent[i*17 +14];
                    if ((fabs(fTemp0) + fabs(fTemp1)) > 0.0f)
                    {   uGlobPos       = i +iFSTART[iRANK];
                        fFslipL[uSlipElCnt[0]*3 +0] = (float)uGlobPos;        fFslipL[uSlipElCnt[0]*3 +1] = fTemp0;               fFslipL[uSlipElCnt[0]*3 +2] = fTemp1;
                        uSlipElCnt[0] += 1u;
                }   }
                iOffstPosF[iRANK] = uSlipElCnt[0]*3;
                
                MPI_Allreduce(MPI_IN_PLACE,uSlipElCnt, 2, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD); 
                MPI_Allreduce(MPI_IN_PLACE, iOffstPosF, iSIZE, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
                for (i = 1; i < iSIZE; i++)        {    iStartPosF[i]  = iStartPosF[i-1] +iOffstPosF[i-1];    } 
                MPI_Allgatherv(fFslipL, iOffstPosF[iRANK], MPI_FLOAT, fFslipG, iOffstPosF, iStartPosF, MPI_FLOAT, MPI_COMM_WORLD);
                
                for (i = 0u; i < iBOFFSET[iRANK]; i++)
                {   if (uSlipElCnt[0] < uKh_BFcnt[i])
                    {   for (j = 0u; j < uSlipElCnt[0]; j++)
                        {   uTemp0 = (unsigned int)fFslipG[j*3 +0]; 
                            uTemp1 = uKh_BFps2[i][uTemp0]; 
                            fTemp0 = fFslipG[j*3 +1]; 
                            fTemp1 = fFslipG[j*3 +2]; 
                            
                            fBEvent[i*9 +0] += (fTemp0*fKh_BFvalstk[i][uTemp1*2 +0] + fTemp1*fKh_BFvalstk[i][uTemp1*2 +1]);
                            fBEvent[i*9 +1] += (fTemp0*fKh_BFvaldip[i][uTemp1*2 +0] + fTemp1*fKh_BFvaldip[i][uTemp1*2 +1]);
                            fBEvent[i*9 +2] += (fTemp0*fKh_BFvalnrm[i][uTemp1*2 +0] + fTemp1*fKh_BFvalnrm[i][uTemp1*2 +1]);
                    }   }
                    else
                    {   memset(fFslip, 0, 2*uFPNum*sizeof(float) );
                        for (j = 0u; j < uSlipElCnt[0]; j++) 
                        {   uTemp0 = (unsigned int)fFslipG[j*3 +0]; 
                            uTemp1 = uKh_BFps2[i][uTemp0]; 
                            fFslip[uTemp1*2 +0] += fFslipG[j*3 +1];
                            fFslip[uTemp1*2 +1] += fFslipG[j*3 +2];
                        }
                        fBEvent[i*9 +0] +=  cblas_sdot(2*uKh_BFcnt[i], fFslip, 1, fKh_BFvalstk[i], 1);
                        fBEvent[i*9 +1] +=  cblas_sdot(2*uKh_BFcnt[i], fFslip, 1, fKh_BFvaldip[i], 1);
                        fBEvent[i*9 +2] +=  cblas_sdot(2*uKh_BFcnt[i], fFslip, 1, fKh_BFvalnrm[i], 1);
                    }
                    fBEvent[i*9 +7] = sqrtf(fBEvent[i*9 +0]*fBEvent[i*9 +0] + fBEvent[i*9 +1]*fBEvent[i*9 +1]);
                    fBEvent[i*9 +8] = fabs(fBEvent[i*9 +2]);
            }   }
            
            
            if (uActElmG >= uMinElemNum4Cat)
            {   uEQcntr += 1u;
                memset(iStartPosF,   0, iSIZE*sizeof(int));           memset(iOffstPosF, 0, iSIZE*sizeof(int));
                memset(fEQMeanVals, 0,     7*sizeof(float));
                
                for (i = 0u; i < uActElmL; i++)
                {   uTemp1   = uFActivEl[i];
                    
                    fTemp0          = sqrtf(fFEvent[uTemp1*17 +13]*fFEvent[uTemp1*17 +13] +fFEvent[uTemp1*17 +14]*fFEvent[uTemp1*17 +14]);
                    fEQMeanVals[0] += (fFRef[uTemp1*14 + 9]*fTemp0*fFRef[uTemp1*14 +0]); 
                    fEQMeanVals[1] += (fFRef[uTemp1*14 +10]*fTemp0*fFRef[uTemp1*14 +0]);
                    fEQMeanVals[2] += (fFRef[uTemp1*14 +11]*fTemp0*fFRef[uTemp1*14 +0]);
                    fEQMeanVals[3] +=  fFRef[uTemp1*14 +0]; 
                    fEQMeanVals[4] +=  fTemp0;
                    fEQMeanVals[5] +=  sqrtf(fFTempVal[uTemp1*5 +0]*fFTempVal[uTemp1*5 +0] + fFTempVal[uTemp1*5 +1]*fFTempVal[uTemp1*5 +1]) - sqrtf(fFEvent[uTemp1*17 +0]*fFEvent[uTemp1*17 +0] + fFEvent[uTemp1*17 +1]*fFEvent[uTemp1*17 +1]); 
                    fEQMeanVals[6] += (float)(uFEvent[uTemp1*5 +2]);
                    
                    uPatchID[i]  = uFActivEl[i] + iFSTART[iRANK];
                    ut0ofPtch[i] = uFEvent[uTemp1*5 +3];
                    uStabType[i] = uFEvent[uTemp1*5 +0];
                    
                    fPtchDTau[i] = sqrtf(fFTempVal[uTemp1*5 +0]*fFTempVal[uTemp1*5 +0] + fFTempVal[uTemp1*5 +1]*fFTempVal[uTemp1*5 +1]) - sqrtf(fFEvent[uTemp1*17 +0]*fFEvent[uTemp1*17 +0] + fFEvent[uTemp1*17 +1]*fFEvent[uTemp1*17 +1]); 
                    fPtchSlpH[i] = fFEvent[uTemp1*17 +13];
                    fPtchSlpV[i] = fFEvent[uTemp1*17 +14];
                }
                iOffstPosF[iRANK] = (int)uActElmL;
                
                MPI_Allreduce( MPI_IN_PLACE, &uMRFlgth,        1,            MPI_UNSIGNED, MPI_MAX, MPI_COMM_WORLD);
                MPI_Allreduce( MPI_IN_PLACE, iOffstPosF,     iSIZE,          MPI_INT, MPI_MAX, MPI_COMM_WORLD);
                MPI_Allreduce( MPI_IN_PLACE, fMRFvals, MAXMOMRATEFUNCLENGTH, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
                MPI_Allreduce( MPI_IN_PLACE, fEQMeanVals,      7,            MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
                
                fEQMeanVals[0] /= fEQMomRel;            fEQMeanVals[1] /= fEQMomRel;            fEQMeanVals[2] /= fEQMomRel;
                fEQMeanVals[4] /= (float) uActElmG;     fEQMeanVals[5] /= (float) uActElmG;
                
                for (i = 1; i < iSIZE; i++)        {    iStartPosF[i]  = iStartPosF[i-1] +iOffstPosF[i-1];    } 
                
                MPI_Gatherv(uPatchID,  iOffstPosF[iRANK], MPI_UNSIGNED, uTemp2Wrt1, iOffstPosF, iStartPosF, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
                MPI_Gatherv(ut0ofPtch, iOffstPosF[iRANK], MPI_UNSIGNED, uTemp2Wrt2, iOffstPosF, iStartPosF, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
                MPI_Gatherv(uStabType, iOffstPosF[iRANK], MPI_UNSIGNED, uTemp2Wrt3, iOffstPosF, iStartPosF, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
                MPI_Gatherv(fPtchDTau, iOffstPosF[iRANK], MPI_FLOAT,    fTemp2Wrt1, iOffstPosF, iStartPosF, MPI_FLOAT,    0, MPI_COMM_WORLD);
                MPI_Gatherv(fPtchSlpH, iOffstPosF[iRANK], MPI_FLOAT,    fTemp2Wrt2, iOffstPosF, iStartPosF, MPI_FLOAT,    0, MPI_COMM_WORLD);
                MPI_Gatherv(fPtchSlpV, iOffstPosF[iRANK], MPI_FLOAT,    fTemp2Wrt3, iOffstPosF, iStartPosF, MPI_FLOAT,    0, MPI_COMM_WORLD);
                
                if (iRANK == 0)
                {   fMaxMagnitude = (fEQMagn > fMaxMagnitude)*fEQMagn + (fEQMagn <= fMaxMagnitude)*fMaxMagnitude;
                    if (uPlotCatalog2Screen == 1u) { fprintf(stdout,"%5u  Time: %6.4lf  (%6.4lf since last)   RA: %8.2f   mSlip: %2.3f   mDtau: %3.3f   MRF %4u   Elem#: %6u   M: %3.2f; #broken: %6u   M_max: %3.2f\n", uEQcntr, dCurrTime, (dCurrTime-dLastEQtime), (fEQMeanVals[3]*1.0E-6f), fEQMeanVals[4], (fEQMeanVals[5]*1.0E-6f), uMRFlgth, uActElmG, fEQMagn,(unsigned int)(fEQMeanVals[6]), fMaxMagnitude);       }
                    
                    fseek(fp_CATALOG, 0, SEEK_END);
                    fwrite(&dCurrTime,  sizeof(double),       1, fp_CATALOG);
                    fwrite(&fEQMagn,    sizeof(float),        1, fp_CATALOG);
                    fwrite(fHypoLoc,    sizeof(float),        3, fp_CATALOG);
                    fwrite(fEQMeanVals, sizeof(float),        6, fp_CATALOG);
                    
                    fwrite(&uMRFlgth,   sizeof(unsigned int), 1, fp_CATALOG);
                    fwrite(fMRFvals,    sizeof(float), uMRFlgth, fp_CATALOG);
                    fwrite(&uActElmG,   sizeof(unsigned),     1, fp_CATALOG);
                    
                    fwrite( uTemp2Wrt1, sizeof(unsigned int), uActElmG, fp_CATALOG);
                    fwrite( uTemp2Wrt2, sizeof(unsigned int), uActElmG, fp_CATALOG);
                    fwrite( fTemp2Wrt1, sizeof(float),        uActElmG, fp_CATALOG);
                    fwrite( fTemp2Wrt2, sizeof(float),        uActElmG, fp_CATALOG);
                    fwrite( fTemp2Wrt3, sizeof(float),        uActElmG, fp_CATALOG);
                    fwrite( uTemp2Wrt3, sizeof(unsigned int), uActElmG, fp_CATALOG);
                    fseek(fp_CATALOG, 0, SEEK_SET);
                    fwrite(&uEQcntr,    sizeof(unsigned int), 1, fp_CATALOG);
                }
                
                dMeanRecurTime += (uEQcntr > 1u)*(dCurrTime - dLastEQtime) + (uEQcntr <= 1u)*0.0;
                dLastEQtime     = dCurrTime; 
            } 
            MPI_Barrier( MPI_COMM_WORLD );
            
            
            if ((uStoreSTF4LargeEQs == 1u) && (fEQMagn >= fMinMag4STF))
            {   memset(lSTF_offset, 0, iSIZE*sizeof(long long int) );           memset(lSTF_start, 0, iSIZE*sizeof(long long int) );
                long long int lTemp0 = 0;
                float *fSTFvals0 = calloc( MAXMOMRATEFUNCLENGTH, sizeof *fSTFvals0 );
                float *fSTFvals1 = calloc( MAXMOMRATEFUNCLENGTH, sizeof *fSTFvals1 );
                
                for (i = 0; i < iFOFFSET[iRANK]; i++)
                {   lSTF_offset[iRANK] = (uFEvent[i*5 +1] == 0u)*lSTF_offset[iRANK]  +  (uFEvent[i*5 +1] != 0u)*(lSTF_offset[iRANK] + 2*sizeof(float) +4*sizeof(unsigned int) +2*uFEvent[i*5 +4]*sizeof(float) );
                }
                
                MPI_Allreduce( MPI_IN_PLACE, lSTF_offset, iSIZE, MPI_LONG_LONG_INT, MPI_MAX, MPI_COMM_WORLD);
                for (i = 1; i < iSIZE; i++)        {    lSTF_start[i]  = lSTF_start[i-1] + lSTF_offset[i-1];    }
                
                strcpy(cSTFName,cInputName);      strcat(cSTFName,"_M");         sprintf(cNameAppend, "%f",fEQMagn);      strcat(cSTFName,cNameAppend);      strcat(cSTFName,"_t");         sprintf(cNameAppend, "%lf",dCurrTime);      strcat(cSTFName,cNameAppend);      strcat(cSTFName,".srfb");
                
                MPI_File_delete(cSTFName, MPI_INFO_NULL);
                MPI_File_open(MPI_COMM_WORLD, cSTFName, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fp_STF);
                
                if (iRANK == 0)
                {   MPI_File_write_at(fp_STF, lTemp0,  &uActElmG,     1, MPI_UNSIGNED,  &STATUS);          lTemp0 += 1*sizeof(unsigned int);
                    MPI_File_write_at(fp_STF, lTemp0,  &fModPara[7],  1, MPI_FLOAT,     &STATUS);          lTemp0 += 1*sizeof(float);
                    fTemp0 = -1.0f;
                    MPI_File_write_at(fp_STF, lTemp0,  &fTemp0,       1, MPI_FLOAT,     &STATUS);          lTemp0 += 1*sizeof(float);
                    MPI_File_write_at(fp_STF, lTemp0,  &fModPara[0],  1, MPI_FLOAT,     &STATUS);          lTemp0 += 1*sizeof(float);
                }
                
                lTemp0 = 1*sizeof(unsigned int) +3*sizeof(float); 
                
                for (i = 0; i < iFOFFSET[iRANK]; i++)
                {   if (uFEvent[i*5 +1] != 0u) 
                    {   uTemp0 = i + iFSTART[iRANK];
                        MPI_File_write_at(fp_STF, (lSTF_start[iRANK]+lTemp0),  &uTemp0,           1, MPI_UNSIGNED,  &STATUS);           lTemp0 += 1*sizeof(unsigned int);
                        MPI_File_write_at(fp_STF, (lSTF_start[iRANK]+lTemp0),  &uFEvent[i*5 +3],  1, MPI_UNSIGNED,  &STATUS);           lTemp0 += 1*sizeof(unsigned int);
                        
                        MPI_File_write_at(fp_STF, (lSTF_start[iRANK]+lTemp0),  &fFEvent[i*17 +13],1, MPI_FLOAT,     &STATUS);           lTemp0 += 1*sizeof(float);
                        MPI_File_write_at(fp_STF, (lSTF_start[iRANK]+lTemp0),  &uFEvent[i*5 +4],  1, MPI_UNSIGNED,  &STATUS);           lTemp0 += 1*sizeof(unsigned int);
                        MPI_File_write_at(fp_STF, (lSTF_start[iRANK]+lTemp0),  &fFEvent[i*17 +14],1, MPI_FLOAT,     &STATUS);           lTemp0 += 1*sizeof(float);
                        MPI_File_write_at(fp_STF, (lSTF_start[iRANK]+lTemp0),  &uFEvent[i*5 +4],  1, MPI_UNSIGNED,  &STATUS);           lTemp0 += 1*sizeof(unsigned int);
                        
                        for (j = 0; j < uFEvent[i*5 +4]; j++)       {   uTemp0 = 2u*i*MAXMOMRATEFUNCLENGTH + 2u*j;          fSTFvals0[j] = fSTFslip[uTemp0 +0];         fSTFvals1[j] = fSTFslip[uTemp0 +1];         }
                        MPI_File_write_at(fp_STF, (lSTF_start[iRANK]+lTemp0), fSTFvals0, uFEvent[i*5 +4], MPI_FLOAT, &STATUS);          lTemp0 += uFEvent[i*5 +4]*sizeof(float);
                        MPI_File_write_at(fp_STF, (lSTF_start[iRANK]+lTemp0), fSTFvals1, uFEvent[i*5 +4], MPI_FLOAT, &STATUS);          lTemp0 += uFEvent[i*5 +4]*sizeof(float);
            }   }   } 
            
            
            fPSeisTime      = 0.0f;
            
            
            {   if (uChgBtwEQs == 1u)
                {   for (i = 0u; i < uActElmL; i++)  
                    {   
                        uTemp1 = uFActivEl[i]; 
                        fTemp0 = ran0(&lSeed);      fTemp0 = fTemp0*2.0f -1.0f;
                        fFEvent[uTemp1*17 +7] = fFRef[uTemp1*14 +7] + (fFRef[uTemp1*14 +8] * fTemp0); 
                            
                        fTemp0 = ran0(&lSeed);      fTemp0 = fTemp0*2.0f -1.0f;
                        fFFric[uTemp1*6 +0] = fFRef[uTemp1*14 +3] + (fFRef[uTemp1*14 +4] * fTemp0); 
                        
                        fTemp0 = ran0(&lSeed);      fTemp0 = fTemp0*2.0f -1.0f;
                        fFFric[uTemp1*6 +1] = fFRef[uTemp1*14 +5] + (fFRef[uTemp1*14 +6] * fTemp0); 
                        
                        fTemp0 =(fFFric[uTemp1*6 +0] - fFFric[uTemp1*6 +1]) *-1.0f*fFRef[uTemp1*14 +1]; 
                        fTemp1 = fTemp0/fFRef[uTemp1*14 +2];
                        uFEvent[uTemp1*5 +0] = (fTemp0 >= 0.0f)*2u + (fTemp0 < 0.0f)*3u; 
                        uFEvent[uTemp1*5 +0] = (fTemp1 > fFEvent[uTemp1*17 +7])*1u + (fTemp1 <= fFEvent[uTemp1*17 +7])*uFEvent[uTemp1*5 +0]; 
                        
                        fTemp0 = fFFric[uTemp1*6 +1]*-1.0f*fFRef[uTemp1*14 +1] + fFRef[uTemp1*14 +2]*fFEvent[uTemp1*17 +7]; 
                        fTemp1 = fFFric[uTemp1*6 +0]*-1.0f*fFRef[uTemp1*14 +1]; 
                        fTemp2 = MAX(fTemp0, fTemp1)/(-1.0f*fFRef[uTemp1*14 +1]); 
                        fFFric[uTemp1*6 +2]  = fTemp2 - fOvershootFract*(fTemp2-fFFric[uTemp1*6 +1]); 
                        
                        fFEvent[uTemp1*17 +10] = (fFFric[uTemp1*6 +1] - fFFric[uTemp1*6 +2])*-1.0f*fFRef[uTemp1*14 +1];
                        
                }   }
                
                float fTempCurrFric;
                for (i = 0u; i < iFOFFSET[iRANK]; i++)
                {   fFRef[i*14 +13] += sqrtf(fFEvent[i*17 +13]*fFEvent[i*17 +13] +fFEvent[i*17 +14]*fFEvent[i*17 +14]);
                    fFEvent[i*17 +2] = fFRef[i*14 +1];
                    fTempCurrFric    = fFFric[i*6 +0];
                    
                    fTemp0           = sqrtf(fFEvent[i*17 +0]*fFEvent[i*17 +0] + fFEvent[i*17 +1]*fFEvent[i*17 +1]);
                    fTemp1           = fTempCurrFric *-1.0f*fFEvent[i*17 +2]; 
                    fTemp2           = MAX((fTemp0/fTemp1), 1.0f);
                    fFEvent[i*17 +0] = (uFEvent[i*5 +2] == 1u)*(fFEvent[i*17 +0]/fTemp2) + (uFEvent[i*5 +2] != 1u)*fFEvent[i*17 +0];
                    fFEvent[i*17 +1] = (uFEvent[i*5 +2] == 1u)*(fFEvent[i*17 +1]/fTemp2) + (uFEvent[i*5 +2] != 1u)*fFEvent[i*17 +1];
                    
                    fTemp0           = sqrtf(fFEvent[i*17 +0]*fFEvent[i*17 +0] + fFEvent[i*17 +1]*fFEvent[i*17 +1]);
                    uTemp0           = (fTemp0 > fTemp1)*1u           + (fTemp0 <= fTemp1)*0u;
                    fTempCurrFric    = (uTemp0 == 1u)*(fTemp0/(-1.0f*fFEvent[i*17 +2])) + (uTemp0 != 1u)*fTempCurrFric ; 
                    fTempCurrFric    = MIN(fTempCurrFric , INTERNALFRICTION);
                    
                    fTemp1           = fTempCurrFric *-1.0f*fFEvent[i*17 +2]; 
                    fTemp2           = MAX((fTemp0/fTemp1), 1.0f);
                    fFEvent[i*17 +0] = (uTemp0 == 1u)*(fFEvent[i*17 +0]/fTemp2)           + (uTemp0 != 1u)*fFEvent[i*17 +0];
                    fFEvent[i*17 +1] = (uTemp0 == 1u)*(fFEvent[i*17 +1]/fTemp2)           + (uTemp0 != 1u)*fFEvent[i*17 +1];
                    
                    fFFric[i*6 +4]   = fTempCurrFric - fFFric[i*6 +0];
                    
                    uFEvent[i*5 +1]   = 0u;         uFEvent[i*5 +2]   = 0u;         uFEvent[i*5 +3]   = 0u;         uFEvent[i*5 +4]   = 0u;
                    fFEvent[i*17 +11] = 0.0f;       fFEvent[i*17 +12] = 0.0f;       fFEvent[i*17 +13] = 0.0f;       fFEvent[i*17 +14] = 0.0f;       fFEvent[i*17 +15] = 0.0f;
                    
                } 
            }
            
            
            timerC = clock() - timerC;        time_CoSeis += (((double)timerC)/CLOCKS_PER_SEC);
            timerI = clock();
        }
    }
    
    MPI_Barrier( MPI_COMM_WORLD );
    
    
    unsigned int *uFLTemp   = calloc( iFOFFSET[iRANK],  sizeof *uFLTemp  );
    unsigned int *uFvTemp   = calloc( uFPNum,  sizeof *uFvTemp  );
    float        *fFLTemp0  = calloc( iFOFFSET[iRANK],  sizeof *fFLTemp0 );
    float        *fFvTemp0  = calloc( uFPNum,  sizeof *fFvTemp0 );
    float        *fFLTemp1  = calloc( iFOFFSET[iRANK],  sizeof *fFLTemp1 );
    float        *fFvTemp1  = calloc( uFPNum,  sizeof *fFvTemp1 );
    float        *fFLTemp2  = calloc( iFOFFSET[iRANK],  sizeof *fFLTemp2 );
    float        *fFvTemp2  = calloc( uFPNum,  sizeof *fFvTemp2 );
    float        *fFLTemp3  = calloc( iFOFFSET[iRANK],  sizeof *fFLTemp3 );
    float        *fFvTemp3  = calloc( uFPNum,  sizeof *fFvTemp3 );
    float        *fFLTemp4  = calloc( iFOFFSET[iRANK],  sizeof *fFLTemp4 );
    float        *fFvTemp4  = calloc( uFPNum,  sizeof *fFvTemp4 );
    float        *fFLTemp5  = calloc( iFOFFSET[iRANK],  sizeof *fFLTemp5 );
    float        *fFvTemp5  = calloc( uFPNum,  sizeof *fFvTemp5 );
    float        *fFLTemp6  = calloc( iFOFFSET[iRANK],  sizeof *fFLTemp6 );
    float        *fFvTemp6  = calloc( uFPNum,  sizeof *fFvTemp6 );
    float        *fFLTemp7  = calloc( iFOFFSET[iRANK],  sizeof *fFLTemp7 );
    float        *fFvTemp7  = calloc( uFPNum,  sizeof *fFvTemp7 );
    float        *fFLTemp8  = calloc( iFOFFSET[iRANK],  sizeof *fFLTemp8 );
    float        *fFvTemp8  = calloc( uFPNum,  sizeof *fFvTemp8 );
    float        *fFLTemp9  = calloc( iFOFFSET[iRANK],  sizeof *fFLTemp9 );
    float        *fFvTemp9  = calloc( uFPNum,  sizeof *fFvTemp9 );
    float        *fFLTemp10 = calloc( iFOFFSET[iRANK],  sizeof *fFLTemp10);
    float        *fFvTemp10 = calloc( uFPNum,  sizeof *fFvTemp10 );
    float        *fFLTemp11 = calloc( iFOFFSET[iRANK],  sizeof *fFLTemp11);
    float        *fFvTemp11 = calloc( uFPNum,  sizeof *fFvTemp11 );
    float        *fBLTemp0  = calloc( iBOFFSET[iRANK],  sizeof *fBLTemp0 );
    float        *fBvTemp0  = calloc( uBPNum,  sizeof *fBvTemp0 );
    float        *fBLTemp1  = calloc( iBOFFSET[iRANK],  sizeof *fBLTemp1 );
    float        *fBvTemp1  = calloc( uBPNum,  sizeof *fBvTemp1 );
    float        *fBLTemp2  = calloc( iBOFFSET[iRANK],  sizeof *fBLTemp2 );
    float        *fBvTemp2  = calloc( uBPNum,  sizeof *fBvTemp2 );
    float        *fBLTemp3  = calloc( iBOFFSET[iRANK],  sizeof *fBLTemp3 );
    float        *fBvTemp3  = calloc( uBPNum,  sizeof *fBvTemp3 );
    float        *fBLTemp4  = calloc( iBOFFSET[iRANK],  sizeof *fBLTemp4 );
    float        *fBvTemp4  = calloc( uBPNum,  sizeof *fBvTemp4 );
    
    for (i = 0u; i < iFOFFSET[iRANK]; i++)
    {   uFLTemp[i]  = uFEvent[i*5  +0];
        fFLTemp0[i] = fFEvent[i*17 +0];      fFLTemp1[i] = fFEvent[i*17 +1];        fFLTemp2[i] = fFEvent[i*17 +2];      fFLTemp3[i] = fFEvent[i*17 +3];
        fFLTemp4[i] = fFEvent[i*17 +7];      fFLTemp5[i] = fFEvent[i*17 +10];       fFLTemp6[i] = fFFric[i*6 +0];        fFLTemp7[i] = fFFric[i*6 +1];
        fFLTemp8[i] = fFFric[i*6 +2];        fFLTemp9[i] = fFFric[i*6 +4];          fFLTemp10[i]= fFRef[i*14 +12];       fFLTemp11[i]= fFRef[i*14 +13];
    }  
    for (i = 0u; i < iBOFFSET[iRANK]; i++)      
    {   fBLTemp0[i] = fBEvent[i*9 +0];       fBLTemp1[i] = fBEvent[i*9 +1];         fBLTemp2[i] = fBEvent[i*9 +2];       fBLTemp3[i] = fBEvent[i*9 +7];
        fBLTemp4[i] = fBEvent[i*9 +8];
    }
    MPI_Gatherv(uFLTemp,  iFOFFSET[iRANK], MPI_UNSIGNED, uFvTemp,  iFOFFSET, iFSTART, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    MPI_Gatherv(fFLTemp0, iFOFFSET[iRANK], MPI_FLOAT,    fFvTemp0, iFOFFSET, iFSTART, MPI_FLOAT,    0, MPI_COMM_WORLD);
    MPI_Gatherv(fFLTemp1, iFOFFSET[iRANK], MPI_FLOAT,    fFvTemp1, iFOFFSET, iFSTART, MPI_FLOAT,    0, MPI_COMM_WORLD);
    MPI_Gatherv(fFLTemp2, iFOFFSET[iRANK], MPI_FLOAT,    fFvTemp2, iFOFFSET, iFSTART, MPI_FLOAT,    0, MPI_COMM_WORLD);
    MPI_Gatherv(fFLTemp3, iFOFFSET[iRANK], MPI_FLOAT,    fFvTemp3, iFOFFSET, iFSTART, MPI_FLOAT,    0, MPI_COMM_WORLD);
    MPI_Gatherv(fFLTemp4, iFOFFSET[iRANK], MPI_FLOAT,    fFvTemp4, iFOFFSET, iFSTART, MPI_FLOAT,    0, MPI_COMM_WORLD);
    MPI_Gatherv(fFLTemp5, iFOFFSET[iRANK], MPI_FLOAT,    fFvTemp5, iFOFFSET, iFSTART, MPI_FLOAT,    0, MPI_COMM_WORLD);
    MPI_Gatherv(fFLTemp6, iFOFFSET[iRANK], MPI_FLOAT,    fFvTemp6, iFOFFSET, iFSTART, MPI_FLOAT,    0, MPI_COMM_WORLD);
    MPI_Gatherv(fFLTemp7, iFOFFSET[iRANK], MPI_FLOAT,    fFvTemp7, iFOFFSET, iFSTART, MPI_FLOAT,    0, MPI_COMM_WORLD);
    MPI_Gatherv(fFLTemp8, iFOFFSET[iRANK], MPI_FLOAT,    fFvTemp8, iFOFFSET, iFSTART, MPI_FLOAT,    0, MPI_COMM_WORLD);
    MPI_Gatherv(fFLTemp9, iFOFFSET[iRANK], MPI_FLOAT,    fFvTemp9, iFOFFSET, iFSTART, MPI_FLOAT,    0, MPI_COMM_WORLD);
    MPI_Gatherv(fFLTemp10,iFOFFSET[iRANK], MPI_FLOAT,    fFvTemp10,iFOFFSET, iFSTART, MPI_FLOAT,    0, MPI_COMM_WORLD);
    MPI_Gatherv(fFLTemp11,iFOFFSET[iRANK], MPI_FLOAT,    fFvTemp11,iFOFFSET, iFSTART, MPI_FLOAT,    0, MPI_COMM_WORLD);
    MPI_Gatherv(fBLTemp0, iBOFFSET[iRANK], MPI_FLOAT,    fBvTemp0, iBOFFSET, iBSTART, MPI_FLOAT,    0, MPI_COMM_WORLD);
    MPI_Gatherv(fBLTemp1, iBOFFSET[iRANK], MPI_FLOAT,    fBvTemp1, iBOFFSET, iBSTART, MPI_FLOAT,    0, MPI_COMM_WORLD);
    MPI_Gatherv(fBLTemp2, iBOFFSET[iRANK], MPI_FLOAT,    fBvTemp2, iBOFFSET, iBSTART, MPI_FLOAT,    0, MPI_COMM_WORLD);
    MPI_Gatherv(fBLTemp3, iBOFFSET[iRANK], MPI_FLOAT,    fBvTemp3, iBOFFSET, iBSTART, MPI_FLOAT,    0, MPI_COMM_WORLD);
    MPI_Gatherv(fBLTemp4, iBOFFSET[iRANK], MPI_FLOAT,    fBvTemp4, iBOFFSET, iBSTART, MPI_FLOAT,    0, MPI_COMM_WORLD);
    
    if (iRANK == 0)
    {   unsigned int uDummy = 3u;
        
        if ((fpPrePost = fopen(cPrevStateName,"wb")) == NULL)   {   printf("Error -cant open %s  PrevStressFile...\n",cPrevStateName);      exit(10);     }
        
        fwrite(&dCurrTime,  sizeof(double),       1,    fpPrePost);
        fwrite(&uEQcntr,    sizeof(unsigned int), 1,    fpPrePost);
        fwrite(&uDummy,     sizeof(unsigned int), 1,    fpPrePost);
        fwrite(&fPSeisTime, sizeof(float),        1,    fpPrePost);

        fwrite(uFvTemp,  sizeof(unsigned int), uFPNum, fpPrePost);              fwrite(fFvTemp0, sizeof(float),        uFPNum, fpPrePost);
        fwrite(fFvTemp1, sizeof(float),        uFPNum, fpPrePost);              fwrite(fFvTemp2, sizeof(float),        uFPNum, fpPrePost);
        fwrite(fFvTemp3, sizeof(float),        uFPNum, fpPrePost);              fwrite(fFvTemp4, sizeof(float),        uFPNum, fpPrePost);
        fwrite(fFvTemp5, sizeof(float),        uFPNum, fpPrePost);              fwrite(fFvTemp6, sizeof(float),        uFPNum, fpPrePost);
        fwrite(fFvTemp7, sizeof(float),        uFPNum, fpPrePost);              fwrite(fFvTemp8, sizeof(float),        uFPNum, fpPrePost);
        fwrite(fFvTemp9, sizeof(float),        uFPNum, fpPrePost);              fwrite(fFvTemp10,sizeof(float),        uFPNum, fpPrePost);
        fwrite(fFvTemp11,sizeof(float),        uFPNum, fpPrePost);
        fwrite(fBvTemp0, sizeof(float),        uBPNum, fpPrePost);              fwrite(fBvTemp1, sizeof(float),        uBPNum, fpPrePost);
        fwrite(fBvTemp2, sizeof(float),        uBPNum, fpPrePost);              fwrite(fBvTemp3, sizeof(float),        uBPNum, fpPrePost);
        fwrite(fBvTemp4, sizeof(float),        uBPNum, fpPrePost);
        
        fclose(fpPrePost);
    }
    MPI_Barrier( MPI_COMM_WORLD );
    
    
    timer0 = clock() - timer0;        time_taken  = ((double)timer0)/CLOCKS_PER_SEC;
    MPI_Allreduce(MPI_IN_PLACE, &time_taken, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); 
    time_taken /=(double)iSIZE;
    
    free(uKh_FFcnt); free(uKh_FFps2); free(fKh_FFvalstk); free(fKh_FFvaldip); free(fKh_FFvalnrm);
    free(uKh_FBcnt); free(uKh_FBps2); free(fKh_FBvalstk); free(fKh_FBvaldip);
    free(uKh_BBcnt); free(uKh_BBps2); free(fKh_BBvalstk); free(fKh_BBvaldip); free(fKh_BBvalnrm);
    free(uKh_BFcnt); free(uKh_BFps2); free(fKh_BFvalstk); free(fKh_BFvaldip); free(fKh_BFvalnrm);
    
    if (iRANK == 0)
    {   
        dMeanRecurTime /= (double)uEQcntr;
        fprintf(stdout,"\nTotal RunTime for EQcycle (minutes): %6.4f       and mean inter-event time (days): %3.2lf\n",time_taken/60.0f, dMeanRecurTime*365.25 );
        fprintf(stdout,"Time spent in interseismic: %4.2lf   Time spent in coseismic: %4.2lf\n",time_InterSeis/60.0f, time_CoSeis/60.0f);
        fprintf(stdout,"PostSeis iterations per event %u\n", (uPSeisSteps/uEQtickr) );
        
        float fSecsPerYear = 31622400.0f,   fNeighborDist = 0.0f,   fP2P1[3],   fP1P3[3];
        
        unsigned int *uTempF  = calloc(uFPNum,   sizeof *uTempF  );
        unsigned int *uFtemp  = calloc(5*uFPNum, sizeof *uFtemp  );
        float *fTempF         = calloc(uFPNum,   sizeof *fTempF  );
        float *fTempFV        = calloc(uFVNum,   sizeof *fTempFV );
        float *fFtemp         = calloc(3*uFPNum, sizeof *fFtemp  );
        float *fFVtemp        = calloc(4*uFVNum, sizeof *fFVtemp );
        float *fF_Area        = calloc(uFPNum,   sizeof *fF_Area );
        
        fseek(fp_CATALOG, (4*sizeof(unsigned int)+2*sizeof(float)), SEEK_SET);
        
        if (fread(uTempF,    sizeof(unsigned int),    uFPNum, fp_CATALOG) != uFPNum)       {   exit(10);   }           for (i = 0u; i < uFPNum; i++)   {   uFtemp[i*5 +0]  = uTempF[i];   }
        if (fread(uTempF,    sizeof(unsigned int),    uFPNum, fp_CATALOG) != uFPNum)       {   exit(10);   }           for (i = 0u; i < uFPNum; i++)   {   uFtemp[i*5 +1]  = uTempF[i];   }
        if (fread(uTempF,    sizeof(unsigned int),    uFPNum, fp_CATALOG) != uFPNum)       {   exit(10);   }           for (i = 0u; i < uFPNum; i++)   {   uFtemp[i*5 +2]  = uTempF[i];   }
        if (fread(uTempF,    sizeof(unsigned int),    uFPNum, fp_CATALOG) != uFPNum)       {   exit(10);   }           for (i = 0u; i < uFPNum; i++)   {   uFtemp[i*5 +3]  = uTempF[i];   }
        if (fread(uTempF,    sizeof(unsigned int),    uFPNum, fp_CATALOG) != uFPNum)       {   exit(10);   }           for (i = 0u; i < uFPNum; i++)   {   uFtemp[i*5 +4]  = uTempF[i];   }

        if (fread(fTempFV,   sizeof(float),           uFVNum, fp_CATALOG) != uFVNum)       {   exit(10);   }           for (i = 0u; i < uFVNum; i++)   {   fFVtemp[i*4 +0] = fTempFV[i];  }
        if (fread(fTempFV,   sizeof(float),           uFVNum, fp_CATALOG) != uFVNum)       {   exit(10);   }           for (i = 0u; i < uFVNum; i++)   {   fFVtemp[i*4 +1] = fTempFV[i];  }
        if (fread(fTempFV,   sizeof(float),           uFVNum, fp_CATALOG) != uFVNum)       {   exit(10);   }           for (i = 0u; i < uFVNum; i++)   {   fFVtemp[i*4 +2] = fTempFV[i];  }
        if (fread(fTempFV,   sizeof(float),           uFVNum, fp_CATALOG) != uFVNum)       {   exit(10);   }           for (i = 0u; i < uFVNum; i++)   {   fFVtemp[i*4 +3] = fTempFV[i];  }

        if (fread(fTempF,    sizeof(float),           uFPNum, fp_CATALOG) != uFPNum)       {   exit(10);   }           for (i = 0u; i < uFPNum; i++)   {   fFtemp[i*3 +0]  = fTempF[i];   }
        if (fread(fTempF,    sizeof(float),           uFPNum, fp_CATALOG) != uFPNum)       {   exit(10);   }           for (i = 0u; i < uFPNum; i++)   {   fFtemp[i*3 +1]  = fTempF[i];   }
        if (fread(fTempF,    sizeof(float),           uFPNum, fp_CATALOG) != uFPNum)       {   exit(10);   }           for (i = 0u; i < uFPNum; i++)   {   fFtemp[i*3 +2]  = fTempF[i];   }
        
        if ((fp_CATALOGCUT = fopen(cOutputNameCUT,"wb")) == NULL)         {   printf("Error -cant open %s  ParsedCatalogFile...\n",cOutputNameCUT);      exit(10);     }
        
        fwrite(&uEQcntr,     sizeof(unsigned int), 1, fp_CATALOGCUT);
        fwrite(&uCatType,    sizeof(unsigned int), 1, fp_CATALOGCUT);
        fwrite(&fModPara[2], sizeof(float),        1, fp_CATALOGCUT);
        fwrite(&fModPara[7], sizeof(float),        1, fp_CATALOGCUT);
        
        if (uCatType == 3u)
        {   
            fwrite(&uFPNum,      sizeof(unsigned int), 1, fp_CATALOGCUT);
            fwrite(&uFVNum,      sizeof(unsigned int), 1, fp_CATALOGCUT);
            
            for (i = 0u; i < uFPNum; i++)   {   uTempF[i] = uFtemp[i*5 +0];   }     fwrite(uTempF,  sizeof(unsigned int), uFPNum, fp_CATALOGCUT);
            for (i = 0u; i < uFPNum; i++)   {   uTempF[i] = uFtemp[i*5 +1];   }     fwrite(uTempF,  sizeof(unsigned int), uFPNum, fp_CATALOGCUT);
            for (i = 0u; i < uFPNum; i++)   {   uTempF[i] = uFtemp[i*5 +2];   }     fwrite(uTempF,  sizeof(unsigned int), uFPNum, fp_CATALOGCUT);
            for (i = 0u; i < uFPNum; i++)   {   uTempF[i] = uFtemp[i*5 +3];   }     fwrite(uTempF,  sizeof(unsigned int), uFPNum, fp_CATALOGCUT);
            for (i = 0u; i < uFPNum; i++)   {   uTempF[i] = uFtemp[i*5 +4];   }     fwrite(uTempF,  sizeof(unsigned int), uFPNum, fp_CATALOGCUT);

            for (i = 0u; i < uFVNum; i++)   {   fTempFV[i]= fFVtemp[i*4 +0];  }     fwrite(fTempFV, sizeof(float),        uFVNum, fp_CATALOGCUT);
            for (i = 0u; i < uFVNum; i++)   {   fTempFV[i]= fFVtemp[i*4 +1];  }     fwrite(fTempFV, sizeof(float),        uFVNum, fp_CATALOGCUT);
            for (i = 0u; i < uFVNum; i++)   {   fTempFV[i]= fFVtemp[i*4 +2];  }     fwrite(fTempFV, sizeof(float),        uFVNum, fp_CATALOGCUT);
            for (i = 0u; i < uFVNum; i++)   {   fTempFV[i]= fFVtemp[i*4 +3];  }     fwrite(fTempFV, sizeof(float),        uFVNum, fp_CATALOGCUT);

            for (i = 0u; i < uFPNum; i++)   {   fTempF[i] = fFtemp[i*3 +0];   }     fwrite(fTempF,  sizeof(float),        uFPNum, fp_CATALOGCUT);
            for (i = 0u; i < uFPNum; i++)   {   fTempF[i] = fFtemp[i*3 +1];   }     fwrite(fTempF,  sizeof(float),        uFPNum, fp_CATALOGCUT);
            for (i = 0u; i < uFPNum; i++)   {   fTempF[i] = fFtemp[i*3 +2];   }     fwrite(fTempF,  sizeof(float),        uFPNum, fp_CATALOGCUT);
        }
        
        {   unsigned int uVert1,   uVert2,   uVert3;
            float fVert1[3],   fVert2[3],   fVert3[3],   fDist[3];
            
            for (i = 0u; i < uFPNum; i++)
            {   uVert1    = uFtemp[i*5 +0];            uVert2    = uFtemp[i*5 +1];            uVert3    = uFtemp[i*5 +2];
                fVert1[0] = fFVtemp[uVert1*4 +0];      fVert1[1] = fFVtemp[uVert1*4 +1];      fVert1[2] = fFVtemp[uVert1*4 +2];
                fVert2[0] = fFVtemp[uVert2*4 +0];      fVert2[1] = fFVtemp[uVert2*4 +1];      fVert2[2] = fFVtemp[uVert2*4 +2];
                fVert3[0] = fFVtemp[uVert3*4 +0];      fVert3[1] = fFVtemp[uVert3*4 +1];      fVert3[2] = fFVtemp[uVert3*4 +2];
                fDist[0]  = sqrtf((fVert2[0] -fVert1[0])*(fVert2[0] -fVert1[0]) + (fVert2[1] -fVert1[1])*(fVert2[1] -fVert1[1]) + (fVert2[2] -fVert1[2])*(fVert2[2] -fVert1[2]));
                fDist[1]  = sqrtf((fVert3[0] -fVert2[0])*(fVert3[0] -fVert2[0]) + (fVert3[1] -fVert2[1])*(fVert3[1] -fVert2[1]) + (fVert3[2] -fVert2[2])*(fVert3[2] -fVert2[2]));
                fDist[2]  = sqrtf((fVert1[0] -fVert3[0])*(fVert1[0] -fVert3[0]) + (fVert1[1] -fVert3[1])*(fVert1[1] -fVert3[1]) + (fVert1[2] -fVert3[2])*(fVert1[2] -fVert3[2]));
                
                fNeighborDist += ((fDist[0] + fDist[1] +fDist[2])/3.0f); 
                
                fP2P1[0]   = fVert2[0] - fVert1[0];      fP2P1[1]  = fVert2[1] - fVert1[1];      fP2P1[2]  = fVert2[2] - fVert1[2];
                fP1P3[0]   = fVert1[0] - fVert3[0];      fP1P3[1]  = fVert1[1] - fVert3[1];      fP1P3[2]  = fVert1[2] - fVert3[2];
                fF_Area[i] =  0.5f*sqrtf( (fP2P1[1]*fP1P3[2]-fP2P1[2]*fP1P3[1])*(fP2P1[1]*fP1P3[2]-fP2P1[2]*fP1P3[1]) + (fP2P1[2]*fP1P3[0]-fP2P1[0]*fP1P3[2])*(fP2P1[2]*fP1P3[0]-fP2P1[0]*fP1P3[2]) + (fP2P1[0]*fP1P3[1]-fP2P1[1]*fP1P3[0])*(fP2P1[0]*fP1P3[1]-fP2P1[1]*fP1P3[0]) );
        }   }
        fNeighborDist *= (NeighborFactor/(float)uFPNum);
        fprintf(stdout,"max. neighbor distance: %f\n",fNeighborDist);
        
        unsigned int **uElNeighbors = NULL;
        
        uElNeighbors  = calloc(uFPNum, sizeof *uElNeighbors );
        
        for (i = 0; i < uFPNum; i++)
        {   uElNeighbors[i] = calloc( uFPNum, sizeof *uElNeighbors[i] );
            for (j = 0; j < uFPNum; j++)
            {   
                fTemp0 = (i == j)*(2.0f*fNeighborDist)  +  (i != j)*sqrtf((fFtemp[i*3 +0]-fFtemp[j*3 +0])*(fFtemp[i*3 +0]-fFtemp[j*3 +0]) + (fFtemp[i*3 +1]-fFtemp[j*3 +1])*(fFtemp[i*3 +1]-fFtemp[j*3 +1]) + (fFtemp[i*3 +2]-fFtemp[j*3 +2])*(fFtemp[i*3 +2]-fFtemp[j*3 +2]));
                
                uElNeighbors[i][0]                    = (fTemp0 >= fNeighborDist)*uElNeighbors[i][0]                     +  (fTemp0 < fNeighborDist)*(uElNeighbors[i][0] +1u);
                uElNeighbors[i][uElNeighbors[i][0]] = (fTemp0 >= fNeighborDist)*uElNeighbors[i][uElNeighbors[i][0]]  +  (fTemp0 < fNeighborDist)*j;
            }
            uElNeighbors[i]  = realloc(uElNeighbors[i], (uElNeighbors[i][0]+1) *sizeof *uElNeighbors[i] );
        }
        
        
        double dEQtime, dEQtimeMod;
        float fMagn,   fHypoLoc[3],   fMeanVals[6];
        unsigned int uMRFlgth,   uElNum;
        
        float *fMRF           = calloc(MAXMOMRATEFUNCLENGTH, sizeof *fMRF );
        unsigned int *uElID   = calloc(uFPNum,       sizeof *uElID );
        
        unsigned int *uStrtT  = calloc(uFPNum,       sizeof *uStrtT );
        unsigned int *uStabT  = calloc(uFPNum,       sizeof *uStabT );
        float *fDtau          = calloc(uFPNum,       sizeof *fDtau );
        float *fSlip_H        = calloc(uFPNum,       sizeof *fSlip_H );
        float *fSlip_V        = calloc(uFPNum,       sizeof *fSlip_V );
        
        unsigned int *ugStrtT = calloc(uFPNum,       sizeof *ugStrtT );
        unsigned int *ugStabT = calloc(uFPNum,       sizeof *ugStabT );
        float *fgDtau         = calloc(uFPNum,       sizeof *fgDtau );
        float *fgSlip_H       = calloc(uFPNum,       sizeof *fgSlip_H );
        float *fgSlip_V       = calloc(uFPNum,       sizeof *fgSlip_V );
        
        unsigned int *ug2Write = calloc(uFPNum,      sizeof *ug2Write );
        float *fg2Write        = calloc(uFPNum,      sizeof *fg2Write );
        
        unsigned int uSelElem,   uTstElem,   uEl2TestCnt,   uNx2TestCnt,   uSubEQElCnt,   uEQcntrNEW = 0u;
        unsigned int uMinStartTime,   uMinStartTPos,   uElLeft;
        unsigned int *uElVisited   = calloc(uFPNum, sizeof *uElVisited );
        unsigned int *uElSlipped   = calloc(uFPNum, sizeof *uElSlipped );
        unsigned int *uEl2Test     = calloc(uFPNum, sizeof *uEl2Test );
        unsigned int *uNx2Test     = calloc(uFPNum, sizeof *uNx2Test );
        unsigned int *uSubEQEl     = calloc(uFPNum, sizeof *uSubEQEl );
        
        for (i = 0; i < uEQcntr; i++)
        {   if (fread(&dEQtime,    sizeof(double),          1,     fp_CATALOG) != 1)              {   exit(10);   }
            if (fread(&fMagn,      sizeof(float),           1,     fp_CATALOG) != 1)              {   exit(10);   }
            if (fread(fHypoLoc,    sizeof(float),           3,     fp_CATALOG) != 3)              {   exit(10);   }
            if (fread(fMeanVals,   sizeof(float),           6,     fp_CATALOG) != 6)              {   exit(10);   }
            if (fread(&uMRFlgth,   sizeof(unsigned int),    1,     fp_CATALOG) != 1)              {   exit(10);   }
            if (fread(fMRF,        sizeof(float),        uMRFlgth, fp_CATALOG) != uMRFlgth)       {   exit(10);   }
            if (fread(&uElNum,     sizeof(unsigned int),    1,     fp_CATALOG) != 1)              {   exit(10);   }
            if (fread(uElID,       sizeof(unsigned int), uElNum,   fp_CATALOG) != uElNum)         {   exit(10);   }
            if (fread(uStrtT,      sizeof(unsigned int), uElNum,   fp_CATALOG) != uElNum)         {   exit(10);   }
            if (fread(fDtau,       sizeof(float),        uElNum,   fp_CATALOG) != uElNum)         {   exit(10);   }
            if (fread(fSlip_H,     sizeof(float),        uElNum,   fp_CATALOG) != uElNum)         {   exit(10);   }
            if (fread(fSlip_V,     sizeof(float),        uElNum,   fp_CATALOG) != uElNum)         {   exit(10);   }
            if (fread(uStabT,      sizeof(unsigned int), uElNum,   fp_CATALOG) != uElNum)         {   exit(10);   }
            
            
            memset(uElSlipped, 0, uFPNum*sizeof(unsigned int) );
            memset(uElVisited, 0, uFPNum*sizeof(unsigned int) );
            
            for (j = 0; j < uElNum; j++)
            {   uElSlipped[uElID[j]] = 1u;                  ugStrtT[uElID[j]]    = uStrtT[j];
                ugStabT[uElID[j]]    = uStabT[j];           fgDtau[uElID[j]]     = fDtau[j];
                fgSlip_H[uElID[j]]   = fSlip_H[j];          fgSlip_V[uElID[j]]   = fSlip_V[j];
            }
            
            uElLeft = uElNum;
            
            while (uElLeft > 0u)
            {   uMinStartTime = UINT_MAX;
                uMinStartTPos = 0;
                uSubEQElCnt   = 0;
                uEQcntrNEW   += 1u;
                
                for (j = 0; j < uElNum; j++)      {   if ((uStrtT[j] < uMinStartTime) && (uElVisited[uElID[j]]) == 0u)        {       uMinStartTime = uStrtT[j];        uMinStartTPos = uElID[j];      }        }
                
                dEQtimeMod   = dEQtime + (double)(((float)uMinStartTime *fModPara[7])/fSecsPerYear);
                uEl2Test[0]  = uMinStartTPos;                     uEl2TestCnt    = 1u;
                uSubEQEl[0]  = uMinStartTPos;                     uSubEQElCnt    = 1u;
                uElLeft     -= 1u;
                uElVisited[uMinStartTPos] = uEQcntrNEW;
                
                while (uEl2TestCnt > 0u)
                {   
                    uNx2TestCnt = 0u;
                    for (j = 0; j < uEl2TestCnt; j++)
                    {   
                        uSelElem = uEl2Test[j];
                        for (k = 1; k <= uElNeighbors[uSelElem][0]; k++)
                        {   uTstElem = uElNeighbors[uSelElem][k];
                            if ((uElSlipped[uTstElem] == 1u) && (uElVisited[uTstElem] == 0u))
                            {   uNx2Test[uNx2TestCnt] = uTstElem;           uNx2TestCnt += 1u;
                                uSubEQEl[uSubEQElCnt] = uTstElem;           uSubEQElCnt += 1u;
                                uElVisited[uTstElem]  = uEQcntrNEW;
                                uElLeft              -= 1u;
                    }   }   }
                    
                    if (uNx2TestCnt > 0u)
                    {   memcpy(uEl2Test, uNx2Test,     uNx2TestCnt* sizeof(unsigned int) );
                        uEl2TestCnt = uNx2TestCnt;
                    }
                    else    {   break;      }
                }
                
                memset(fMeanVals, 0, 6*sizeof(float) );
                fTemp2      = 0.0f;
                fHypoLoc[0] = fFtemp[uSubEQEl[0]*3 +0];            fHypoLoc[1] = fFtemp[uSubEQEl[0]*3 +1];            fHypoLoc[2] = fFtemp[uSubEQEl[0]*3 +2];
                
                for (j = 0; j < uSubEQElCnt; j++)
                {   uTemp0        = uSubEQEl[j]; 
                    fTemp0        = sqrtf( fgSlip_H[uTemp0]*fgSlip_H[uTemp0] + fgSlip_V[uTemp0]*fgSlip_V[uTemp0] );
                    fTemp1        = fTemp0*fF_Area[uTemp0];
                    fTemp2       += fTemp1;
                    
                    fMeanVals[0] += (fFtemp[uTemp0*3 +0] *fTemp1);
                    fMeanVals[1] += (fFtemp[uTemp0*3 +1] *fTemp1);
                    fMeanVals[2] += (fFtemp[uTemp0*3 +2] *fTemp1);
                    fMeanVals[3] += fF_Area[uTemp0];
                    fMeanVals[4] += fTemp0;
                    fMeanVals[5] += fgDtau[uTemp0];
                }
                
                fMeanVals[0] /= fTemp2;                 fMeanVals[1] /= fTemp2;             fMeanVals[2] /= fTemp2;
                fMeanVals[4] /= (float)uSubEQElCnt;     fMeanVals[5] /= (float)uSubEQElCnt;
                fMagn = (log10f(fTemp2*fModPara[2]) -9.1f)/1.5f;
                
                fwrite(&dEQtimeMod, sizeof(double),       1,   fp_CATALOGCUT);
                fwrite(&fMagn,      sizeof(float),        1,   fp_CATALOGCUT);
                fwrite(&fHypoLoc,   sizeof(float),        3,   fp_CATALOGCUT);
                fwrite(&fMeanVals,  sizeof(float),        6,   fp_CATALOGCUT);
                if (uCatType > 1u)
                {   fwrite(&uMRFlgth,   sizeof(unsigned int), 1,   fp_CATALOGCUT);
                    fwrite(fMRF,        sizeof(float),   uMRFlgth, fp_CATALOGCUT);
                }
                if (uCatType > 2u)
                {   fwrite(&uSubEQElCnt,sizeof(unsigned int), 1,   fp_CATALOGCUT);
                    for (j = 0; j < uSubEQElCnt; j++)   {   ug2Write[j] = uSubEQEl[j];            }      fwrite(ug2Write,  sizeof(unsigned int), uSubEQElCnt, fp_CATALOGCUT);
                    for (j = 0; j < uSubEQElCnt; j++)   {   ug2Write[j] = ugStrtT[uSubEQEl[j]];   }      fwrite(ug2Write,  sizeof(unsigned int), uSubEQElCnt, fp_CATALOGCUT);
                    for (j = 0; j < uSubEQElCnt; j++)   {   fg2Write[j] = fgDtau[uSubEQEl[j]];    }      fwrite(fg2Write,  sizeof(float),        uSubEQElCnt, fp_CATALOGCUT);
                    for (j = 0; j < uSubEQElCnt; j++)   {   fg2Write[j] = fgSlip_H[uSubEQEl[j]];  }      fwrite(fg2Write,  sizeof(float),        uSubEQElCnt, fp_CATALOGCUT);
                    for (j = 0; j < uSubEQElCnt; j++)   {   fg2Write[j] = fgSlip_V[uSubEQEl[j]];  }      fwrite(fg2Write,  sizeof(float),        uSubEQElCnt, fp_CATALOGCUT);
                    for (j = 0; j < uSubEQElCnt; j++)   {   ug2Write[j] = ugStabT[uSubEQEl[j]];   }      fwrite(ug2Write,  sizeof(unsigned int), uSubEQElCnt, fp_CATALOGCUT);
        }   }   }
        
        fseek(fp_CATALOGCUT, 0, SEEK_SET);
        fwrite(&uEQcntrNEW, sizeof(unsigned int), 1, fp_CATALOGCUT);
        fprintf(stdout,"EQ number after parsing:  %u\n",uEQcntrNEW);
        
        fclose(fp_CATALOG);
        fclose(fp_CATALOGCUT);
        
    }
    
    
    
    MPI_Finalize( );
    return 0;
 }




 float ran0(long *idum)
{   long k;
    float ans;
    
    *idum ^= MASK;
    k = (*idum)/IQ;
    *idum = IA*(*idum-k*IQ) -IR*k;
    if (*idum < 0) *idum += IM;
    ans = AM*(*idum);
    *idum ^= MASK;

    return ans;
}


float FloatPow(float Base, unsigned int Exp)
{   unsigned int i;
    float p = 1.0f;
    for (i = 1u; i <= Exp; i++)       {   p *= Base;   }
    return p;
}


void   ScaleStress6(float fStress[6], float fScaleFact)
{   fStress[0] *= fScaleFact;           fStress[1] *= fScaleFact;           fStress[2] *= fScaleFact;
    fStress[3] *= fScaleFact;           fStress[4] *= fScaleFact;           fStress[5] *= fScaleFact;
    return;
}


void RotStressG2L(float fOutTens[3], float fR[9], float fT[6])
{   fOutTens[2] = (fR[0]*fT[0] + fR[1]*fT[1] + fR[2]*fT[2])*fR[0] + (fR[0]*fT[1] + fR[1]*fT[3] + fR[2]*fT[4])*fR[1] + (fR[0]*fT[2] + fR[1]*fT[4] + fR[2]*fT[5])*fR[2];
    fOutTens[0] = (fR[3]*fT[0] + fR[4]*fT[1] + fR[5]*fT[2])*fR[0] + (fR[3]*fT[1] + fR[4]*fT[3] + fR[5]*fT[4])*fR[1] + (fR[3]*fT[2] + fR[4]*fT[4] + fR[5]*fT[5])*fR[2];
    fOutTens[1] = (fR[6]*fT[0] + fR[7]*fT[1] + fR[8]*fT[2])*fR[0] + (fR[6]*fT[1] + fR[7]*fT[3] + fR[8]*fT[4])*fR[1] + (fR[6]*fT[2] + fR[7]*fT[4] + fR[8]*fT[5])*fR[2];
    return;
}


void  SubtractVect3(float fDiff[3], float fMinu[3], float fSubt[3])
{   fDiff[0] = fMinu[0] - fSubt[0];
    fDiff[1] = fMinu[1] - fSubt[1];
    fDiff[2] = fMinu[2] - fSubt[2];
    return;
}


void  CrossProduct(float fXProd[3], float fVect1[3], float fVect2[3])
{   fXProd[0] = fVect1[1]*fVect2[2] - fVect1[2]*fVect2[1];
    fXProd[1] = fVect1[2]*fVect2[0] - fVect1[0]*fVect2[2];
    fXProd[2] = fVect1[0]*fVect2[1] - fVect1[1]*fVect2[0];
    return;
}


float VectorLength(float fTempVect[3])
{   float fVectLgth = 0.0f;
    fVectLgth = sqrtf(fTempVect[0]*fTempVect[0] + fTempVect[1]*fTempVect[1] + fTempVect[2]*fTempVect[2]);
    return fVectLgth;
}


void    NrmlzeVect(float fTempVect[3])
{   float fVectLgth = 0.0f;
    fVectLgth = sqrtf(fTempVect[0]*fTempVect[0] + fTempVect[1]*fTempVect[1] + fTempVect[2]*fTempVect[2]);
    fTempVect[0] /= fVectLgth;              fTempVect[1] /= fVectLgth;              fTempVect[2] /= fVectLgth; 
    return;
}


void FindMainStressDir(float fSMaxVect[3], float fStrIn[6])
{   
    int i, j, ip, iq, maxIndx;
    float sm, theta, thres, g, h, t, c, tau, s, maxEig;
    float b[3], z[3];
    float fStress[9],   fEigVect[9],   fEigVal[3];
    fStress[0] = fStrIn[0];         fStress[1] = fStrIn[1];         fStress[2] = fStrIn[2];
    fStress[3] = fStrIn[1];         fStress[4] = fStrIn[3];         fStress[5] = fStrIn[4];
    fStress[6] = fStrIn[2];         fStress[7] = fStrIn[4];         fStress[8] = fStrIn[5];
    
    for (i = 0; i < 3; i++)
    {   for (j = 0; j < 3; j++)             {        fEigVect[i*3 +j] = 0.0f;           }
        fEigVect[i*3 +i] = 1.0f;
        fEigVal[i]       = fStress[i*3 +i];
        b[i]             = fStress[i*3 +i];
        z[i]             = 0.0f;  
    }
    
    for (i = 0; i < 50; i++)
    {   sm = 0.0f;
        for (ip = 0; ip < 2; ip++)
        {   for (iq = (ip+1); iq < 3; iq++)
            {   sm += fabs(fStress[ip*3 +iq]);
        }   }
        if (sm <= FLT_EPSILON)
        {  
            
            maxIndx = 0;
            maxEig  = fabs(fEigVal[0]);
            for (j = 1; j < 3; j++)
            { if (fabs(fEigVal[j]) > maxEig)    {       maxEig = fabs(fEigVal[j]);      maxIndx = j;        }
            }
            for (j = 0; j < 3; j++)
            {   fSMaxVect[j] = fEigVect[j*3 +maxIndx] *fEigVal[maxIndx];
            }
            
            return;
        }
        
        if (i < 4)              {   thres = 0.2f*sm/9.0f;       }
        else                    {   thres = 0.0f;               }
        
        for (ip = 0; ip < 2; ip++)
        {   for (iq = (ip+1); iq < 3; iq++)
            {   g = 100.0f*fabs(fStress[ip*3 +iq]);
        
                if ( (i > 4) && ( fabs(fEigVal[ip]) +g == fabs(fEigVal[ip]) ) &&  ( fabs(fEigVal[iq]) +g == fabs(fEigVal[iq]) ) )
                {   fStress[ip*3 +iq] = 0.0f;   
                }
                else if (fabs(fStress[ip*3 +iq]) > thres  )
                {   h = fEigVal[iq] - fEigVal[ip];
                    if (fabs(h) +g == fabs(h))
                    {   t = fStress[ip*3 +iq]/h;
                    }
                    else
                    {   theta = 0.5f*h/fStress[ip*3 +iq];
                        t     = 1.0f/(fabs(theta) +sqrtf(1.0f +theta*theta) );
                        if (theta < 0.0f)   {   t *= -1.0f;       }
                    }
                    c   = 1.0f/sqrtf(1.0f +t*t);
                    s   = t*c;
                    tau = s/(1.0f+c);
                    h   = t*fStress[ip*3 +iq];
                    z[ip      ] -= h;
                    z[iq]       += h;
                    fEigVal[ip] -= h;
                    fEigVal[iq] += h;
                    
                    fStress[ip*3 +iq] = 0.0f;
                    
                    for (j = 0; j < ip; j++)
                    {   g = fStress[j*3 +ip];                   h = fStress[j*3 +iq];
                        fStress[j*3 +ip] = g -s*(h+g*tau);      fStress[j*3 +iq] = h +s*(g-h*tau);
                    }
                    for (j = (ip+1); j < iq; j++)
                    {   g = fStress[ip*3 +j];                   h = fStress[j*3 +iq];
                        fStress[ip*3 +j] = g -s*(h+g*tau);      fStress[j*3 +iq] = h +s*(g-h*tau);
                    }
                    for (j = (iq+1); j < 3; j++)
                    {   g = fStress[ip*3 +j];                   h = fStress[iq*3 +j];
                        fStress[ip*3 +j] = g -s*(h+g*tau);      fStress[iq*3 +j] = h +s*(g-h*tau);
                    }
                    
                    for (j = 0; j < 3; j++)
                    {   g = fEigVect[j*3 +ip];                  h = fEigVect[j*3 +iq];
                        fEigVect[j*3 +ip] = g -s*(h+g*tau);     fEigVect[j*3 +iq] = h +s*(g-h*tau);
        }   }   }   }
        
        for (ip = 0; ip < 2; ip++)
        {   b[ip]      += z[ip];
            fEigVal[ip] = b[ip];
            z[ip]       = 0.0f;
        }
        
    }
    return;
}


void GetStkAndDipVect(float fNrm[3], float fStk[3], float fDip[3])
{   
    float fVectLgth, feZ[3];
    feZ[0] = 0.0f;                          feZ[1] = 0.0f;                              feZ[2] = 1.0f;
    CrossProduct(fStk, feZ, fNrm);          fVectLgth = VectorLength(fStk);
    
    if (fVectLgth <= FLT_EPSILON)
    {   fStk[0]    = 0.0f;                  fStk[1] = 1.0f;                             fStk[2] = 0.0f;                   }
    else
    {   fStk[0]   /= fVectLgth;             fStk[1]/= fVectLgth;                        fStk[2]/= fVectLgth;              }
    
    CrossProduct(fDip, fNrm, fStk);         NrmlzeVect(fDip);
    
    return;
}


void GetGlobVertsForRectangle(float fP1[3], float fP2[3], float fP3[3], float fP4[3], float *fKh_Brch, unsigned int uStart)
{   
    float fMinSL,   fMaxSL,   fMinDL,   fMaxDL;
    float fPt1[3],   fPt2[3],   fPt3[3],   fPt4[3],   fPt5[3],   fPt6[3],   fPt7[3],   fPt8[3];
    float fPt1L[3],   fPt2L[3],   fPt3L[3],   fPt4L[3],   fPt5L[3],   fPt6L[3],   fPt7L[3],   fPt8L[3];
    float fRM_G2L[3][3];
    float fRM_L2G[3][3];
    float fEmid = 0.5f*(fKh_Brch[uStart +15] + fKh_Brch[uStart +16]);
    float fNmid = 0.5f*(fKh_Brch[uStart +17] + fKh_Brch[uStart +18]);
    float fZmid = 0.5f*(fKh_Brch[uStart +19] + fKh_Brch[uStart +20]);
    float fAspRatio, fSLdist, fDLdist;
    
    fRM_L2G[0][0] = fKh_Brch[uStart +7];        fRM_L2G[0][1] = fKh_Brch[uStart +10];       fRM_L2G[0][2] = fKh_Brch[uStart +4];
    fRM_L2G[1][0] = fKh_Brch[uStart +8];        fRM_L2G[1][1] = fKh_Brch[uStart +11];       fRM_L2G[1][2] = fKh_Brch[uStart +5];
    fRM_L2G[2][0] = fKh_Brch[uStart +9];        fRM_L2G[2][1] = fKh_Brch[uStart +12];       fRM_L2G[2][2] = fKh_Brch[uStart +6];
    
    fRM_G2L[0][0] = fKh_Brch[uStart +7];        fRM_G2L[0][1] = fKh_Brch[uStart  +8];       fRM_G2L[0][2] = fKh_Brch[uStart +9];
    fRM_G2L[1][0] = fKh_Brch[uStart+10];        fRM_G2L[1][1] = fKh_Brch[uStart +11];       fRM_G2L[1][2] = fKh_Brch[uStart+12];
    fRM_G2L[2][0] = fKh_Brch[uStart +4];        fRM_G2L[2][1] = fKh_Brch[uStart  +5];       fRM_G2L[2][2] = fKh_Brch[uStart +6];
    
    fPt1[0] = fKh_Brch[uStart +15] - fEmid;     fPt1[1] = fKh_Brch[uStart +17] - fNmid;     fPt1[2] = fKh_Brch[uStart +19] - fZmid; 
    fPt2[0] = fKh_Brch[uStart +15] - fEmid;     fPt2[1] = fKh_Brch[uStart +17] - fNmid;     fPt2[2] = fKh_Brch[uStart +20] - fZmid; 
    fPt3[0] = fKh_Brch[uStart +15] - fEmid;     fPt3[1] = fKh_Brch[uStart +18] - fNmid;     fPt3[2] = fKh_Brch[uStart +19] - fZmid; 
    fPt4[0] = fKh_Brch[uStart +16] - fEmid;     fPt4[1] = fKh_Brch[uStart +17] - fNmid;     fPt4[2] = fKh_Brch[uStart +19] - fZmid; 
    fPt5[0] = fKh_Brch[uStart +15] - fEmid;     fPt5[1] = fKh_Brch[uStart +18] - fNmid;     fPt5[2] = fKh_Brch[uStart +20] - fZmid; 
    fPt6[0] = fKh_Brch[uStart +16] - fEmid;     fPt6[1] = fKh_Brch[uStart +17] - fNmid;     fPt6[2] = fKh_Brch[uStart +20] - fZmid; 
    fPt7[0] = fKh_Brch[uStart +16] - fEmid;     fPt7[1] = fKh_Brch[uStart +18] - fNmid;     fPt7[2] = fKh_Brch[uStart +19] - fZmid; 
    fPt8[0] = fKh_Brch[uStart +16] - fEmid;     fPt8[1] = fKh_Brch[uStart +18] - fNmid;     fPt8[2] = fKh_Brch[uStart +20] - fZmid; 
    
    fPt1L[0] = fRM_G2L[0][0]*fPt1[0] + fRM_G2L[0][1]*fPt1[1] +  fRM_G2L[0][2]*fPt1[2];
    fPt1L[1] = fRM_G2L[1][0]*fPt1[0] + fRM_G2L[1][1]*fPt1[1] +  fRM_G2L[1][2]*fPt1[2];
    fPt1L[2] = fRM_G2L[2][0]*fPt1[0] + fRM_G2L[2][1]*fPt1[1] +  fRM_G2L[2][2]*fPt1[2];
    
    fPt2L[0] = fRM_G2L[0][0]*fPt2[0] + fRM_G2L[0][1]*fPt2[1] +  fRM_G2L[0][2]*fPt2[2];
    fPt2L[1] = fRM_G2L[1][0]*fPt2[0] + fRM_G2L[1][1]*fPt2[1] +  fRM_G2L[1][2]*fPt2[2];
    fPt2L[2] = fRM_G2L[2][0]*fPt2[0] + fRM_G2L[2][1]*fPt2[1] +  fRM_G2L[2][2]*fPt2[2];
    
    fPt3L[0] = fRM_G2L[0][0]*fPt3[0] + fRM_G2L[0][1]*fPt3[1] +  fRM_G2L[0][2]*fPt3[2];
    fPt3L[1] = fRM_G2L[1][0]*fPt3[0] + fRM_G2L[1][1]*fPt3[1] +  fRM_G2L[1][2]*fPt3[2];
    fPt3L[2] = fRM_G2L[2][0]*fPt3[0] + fRM_G2L[2][1]*fPt3[1] +  fRM_G2L[2][2]*fPt3[2];
    
    fPt4L[0] = fRM_G2L[0][0]*fPt4[0] + fRM_G2L[0][1]*fPt4[1] +  fRM_G2L[0][2]*fPt4[2];
    fPt4L[1] = fRM_G2L[1][0]*fPt4[0] + fRM_G2L[1][1]*fPt4[1] +  fRM_G2L[1][2]*fPt4[2];
    fPt4L[2] = fRM_G2L[2][0]*fPt4[0] + fRM_G2L[2][1]*fPt4[1] +  fRM_G2L[2][2]*fPt4[2];
    
    fPt5L[0] = fRM_G2L[0][0]*fPt5[0] + fRM_G2L[0][1]*fPt5[1] +  fRM_G2L[0][2]*fPt5[2];
    fPt5L[1] = fRM_G2L[1][0]*fPt5[0] + fRM_G2L[1][1]*fPt5[1] +  fRM_G2L[1][2]*fPt5[2];
    fPt5L[2] = fRM_G2L[2][0]*fPt5[0] + fRM_G2L[2][1]*fPt5[1] +  fRM_G2L[2][2]*fPt5[2];
    
    fPt6L[0] = fRM_G2L[0][0]*fPt6[0] + fRM_G2L[0][1]*fPt6[1] +  fRM_G2L[0][2]*fPt6[2];
    fPt6L[1] = fRM_G2L[1][0]*fPt6[0] + fRM_G2L[1][1]*fPt6[1] +  fRM_G2L[1][2]*fPt6[2];
    fPt6L[2] = fRM_G2L[2][0]*fPt6[0] + fRM_G2L[2][1]*fPt6[1] +  fRM_G2L[2][2]*fPt6[2];
    
    fPt7L[0] = fRM_G2L[0][0]*fPt7[0] + fRM_G2L[0][1]*fPt7[1] +  fRM_G2L[0][2]*fPt7[2];
    fPt7L[1] = fRM_G2L[1][0]*fPt7[0] + fRM_G2L[1][1]*fPt7[1] +  fRM_G2L[1][2]*fPt7[2];
    fPt7L[2] = fRM_G2L[2][0]*fPt7[0] + fRM_G2L[2][1]*fPt7[1] +  fRM_G2L[2][2]*fPt7[2];
    
    fPt8L[0] = fRM_G2L[0][0]*fPt8[0] + fRM_G2L[0][1]*fPt8[1] +  fRM_G2L[0][2]*fPt8[2];
    fPt8L[1] = fRM_G2L[1][0]*fPt8[0] + fRM_G2L[1][1]*fPt8[1] +  fRM_G2L[1][2]*fPt8[2];
    fPt8L[2] = fRM_G2L[2][0]*fPt8[0] + fRM_G2L[2][1]*fPt8[1] +  fRM_G2L[2][2]*fPt8[2];
    
    fMinSL = fPt1L[0];       fMaxSL = fPt1L[0];       fMinDL = fPt1L[1];       fMaxDL = fPt1L[1];
    
    fMinSL = MIN(fMinSL, fPt2L[0]);   fMinSL = MIN(fMinSL, fPt3L[0]);   fMinSL = MIN(fMinSL, fPt4L[0]);   fMinSL = MIN(fMinSL, fPt5L[0]);   fMinSL = MIN(fMinSL, fPt6L[0]);   fMinSL = MIN(fMinSL, fPt7L[0]);   fMinSL = MIN(fMinSL, fPt8L[0]);
    fMaxSL = MAX(fMaxSL, fPt2L[0]);   fMaxSL = MAX(fMaxSL, fPt3L[0]);   fMaxSL = MAX(fMaxSL, fPt4L[0]);   fMaxSL = MAX(fMaxSL, fPt5L[0]);   fMaxSL = MAX(fMaxSL, fPt6L[0]);   fMaxSL = MAX(fMaxSL, fPt7L[0]);   fMaxSL = MAX(fMaxSL, fPt8L[0]);
    
    fMinDL = MIN(fMinDL, fPt2L[1]);   fMinDL = MIN(fMinDL, fPt3L[1]);   fMinDL = MIN(fMinDL, fPt4L[1]);   fMinDL = MIN(fMinDL, fPt5L[1]);   fMinDL = MIN(fMinDL, fPt6L[1]);   fMinDL = MIN(fMinDL, fPt7L[1]);   fMinDL = MIN(fMinDL, fPt8L[1]);
    fMaxDL = MAX(fMaxDL, fPt2L[1]);   fMaxDL = MAX(fMaxDL, fPt3L[1]);   fMaxDL = MAX(fMaxDL, fPt4L[1]);   fMaxDL = MAX(fMaxDL, fPt5L[1]);   fMaxDL = MAX(fMaxDL, fPt6L[1]);   fMaxDL = MAX(fMaxDL, fPt7L[1]);   fMaxDL = MAX(fMaxDL, fPt8L[1]);
    
    fAspRatio = (fMaxSL - fMinSL)/(fMaxDL - fMinDL);
    fSLdist   = sqrtf(fKh_Brch[uStart +3] *fAspRatio);
    fDLdist   = fSLdist/fAspRatio;
    
    if (fabs(fKh_Brch[uStart+12]) > FLT_EPSILON ) 
    {   fDLdist = MIN(fDLdist, 2.0f*fabs(fZmid/fKh_Brch[uStart +12]));
        fSLdist = fKh_Brch[uStart +3]/fDLdist;
    }
    
    fPt1L[0]  = -0.5f*fSLdist;          fPt1L[1]  = -0.5f*fDLdist;              fPt1L[2]  = 0.0f;
    fPt2L[0]  =  0.5f*fSLdist;          fPt2L[1]  = -0.5f*fDLdist;              fPt2L[2]  = 0.0f;
    fPt3L[0]  =  0.5f*fSLdist;          fPt3L[1]  =  0.5f*fDLdist;              fPt3L[2]  = 0.0f;
    fPt4L[0]  = -0.5f*fSLdist;          fPt4L[1]  =  0.5f*fDLdist;              fPt4L[2]  = 0.0f;
    
    fP1[0] = fRM_L2G[0][0]*fPt1L[0] + fRM_L2G[0][1]*fPt1L[1] + fRM_L2G[0][2]*fPt1L[2] + fEmid;
    fP1[1] = fRM_L2G[1][0]*fPt1L[0] + fRM_L2G[1][1]*fPt1L[1] + fRM_L2G[1][2]*fPt1L[2] + fNmid;
    fP1[2] = fRM_L2G[2][0]*fPt1L[0] + fRM_L2G[2][1]*fPt1L[1] + fRM_L2G[2][2]*fPt1L[2] + fZmid;

    fP2[0] = fRM_L2G[0][0]*fPt2L[0] + fRM_L2G[0][1]*fPt2L[1] + fRM_L2G[0][2]*fPt2L[2] + fEmid;
    fP2[1] = fRM_L2G[1][0]*fPt2L[0] + fRM_L2G[1][1]*fPt2L[1] + fRM_L2G[1][2]*fPt2L[2] + fNmid;
    fP2[2] = fRM_L2G[2][0]*fPt2L[0] + fRM_L2G[2][1]*fPt2L[1] + fRM_L2G[2][2]*fPt2L[2] + fZmid;
    
    fP3[0] = fRM_L2G[0][0]*fPt3L[0] + fRM_L2G[0][1]*fPt3L[1] + fRM_L2G[0][2]*fPt3L[2] + fEmid;
    fP3[1] = fRM_L2G[1][0]*fPt3L[0] + fRM_L2G[1][1]*fPt3L[1] + fRM_L2G[1][2]*fPt3L[2] + fNmid;
    fP3[2] = fRM_L2G[2][0]*fPt3L[0] + fRM_L2G[2][1]*fPt3L[1] + fRM_L2G[2][2]*fPt3L[2] + fZmid;
    
    fP4[0] = fRM_L2G[0][0]*fPt4L[0] + fRM_L2G[0][1]*fPt4L[1] + fRM_L2G[0][2]*fPt4L[2] + fEmid;
    fP4[1] = fRM_L2G[1][0]*fPt4L[0] + fRM_L2G[1][1]*fPt4L[1] + fRM_L2G[1][2]*fPt4L[2] + fNmid;
    fP4[2] = fRM_L2G[2][0]*fPt4L[0] + fRM_L2G[2][1]*fPt4L[1] + fRM_L2G[2][2]*fPt4L[2] + fZmid;
    
    fP1[2] = (fP1[2] > 0.0f)*0.0f + (fP1[2] <= 0.0f)*fP1[2];        fP2[2] = (fP2[2] > 0.0f)*0.0f + (fP2[2] <= 0.0f)*fP2[2];
    fP3[2] = (fP3[2] > 0.0f)*0.0f + (fP3[2] <= 0.0f)*fP3[2];        fP4[2] = (fP4[2] > 0.0f)*0.0f + (fP4[2] <= 0.0f)*fP4[2];
    return;
}


