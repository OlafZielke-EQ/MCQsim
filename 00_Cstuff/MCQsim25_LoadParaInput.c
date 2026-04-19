#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#define ALIGN64 64 
#define VersionNumber 202511


struct _5r5uqpg 
{   long _0cdk36m; 
    double _0ur5l8a; 
    float _cmktlzk; 
    float _iwvhdk8; 
    float _kkrzd5e; 
    float _wrb7f6f; 
    float _6snglh6; 
    float _douadvw; 
    float _mzwffhh; 
    float _14ve92m; 
    float _rtj5pht; 
    float _dswm31y; 
    unsigned int _at9mx2x; 
};
struct _xjuy0bv
{   float _rdihsa7; 
    float _lw2lm34; 
    float _e1g983e; 
    float _3jczucq; 
    float _x5qjicj; 
    float _7odkef8; 
    unsigned int _3ru7myo; 
    unsigned int _5izolnn; 
    unsigned int _e6zvlhy; 
    unsigned int _dn1th1x; 
    unsigned int _w2maiet; 
    unsigned short int _wk8temp; 
    unsigned short int _qdbxvh7; 
    unsigned short int _pwh0mx1; 
};
struct _j49itn7
{   float _b68rerv; 
    float _gfu8iog; 
    float _ycz04cd; 
    float _cb1kjtg; 
    float _xf7jxjn; 
    float _z2c9ros; 
    float _dhgzz0j; 
    float _a11dpg8; 
    float _sbi4q2y; 
    float _9v86wzs; 
    unsigned int _fu7bfxq; 
    unsigned int _ryf0992; 
    unsigned short int _6wmdw0f; 
    unsigned short int _rp9dxoq; 
    unsigned short int _jcmrsez; 
    unsigned short int _qy8yb7w; 
    unsigned short int _o7vbgr8; 
    unsigned short int _s832ejx; 
};


unsigned int _ktsg7fx(unsigned int _zdo74v1) {
    volatile unsigned int _ocvitwu      = 0xF3A1C2D4;
    volatile unsigned int _163yk03  = 0xF3A00384;
    unsigned int _soyg9zq = _ocvitwu ^ _163yk03;
    if (_zdo74v1 >= _soyg9zq)
    {   return _ocvitwu ^ (_zdo74v1 * 0x1F);  
    }
    return _zdo74v1;
}


void   _148p0tv(int _fsqg8x8, const char *restrict _ut076g5, char *restrict _g79or78, char *restrict _dygflij,  struct _5r5uqpg *restrict _atj9tze , struct _xjuy0bv *restrict _vxpj4hn, struct _j49itn7 *restrict _s55j2o7)
{   FILE        *_yh1fkj4,                 *_yk1b3d5,                      *_050ui22;  
    char        _u1g463s[512],        _kzmpjrh[512],             _75aud9k[512];
    char        _ycelos6[512];
    unsigned int _8p74mjd[2];
    
    if ((_yh1fkj4 = fopen(_ut076g5,"r"))    == NULL)         {   perror("Error: Cant open input file in function 'LoadInputParaFile'\n");        exit(EXIT_FAILURE);                 }
    {
        if (fgets(_ycelos6, 512, _yh1fkj4)    != NULL)      {   sscanf(_ycelos6,"%*s %u", &_8p74mjd[0]);                 }
        if (_8p74mjd[0] != VersionNumber)                 {   perror("VersionNumber of code and parameter file are not identical. This might cause undefined behavior.\n");       }
        if (fgets(_ycelos6, 512, _yh1fkj4)    != NULL)      {                                                           }
        if (fgets(_ycelos6, 512, _yh1fkj4)    != NULL)      {                                                           }
        if (fgets(_ycelos6, 512, _yh1fkj4)    != NULL)      {   sscanf(_ycelos6,"%*s %s", _dygflij);                   }

        if (fgets(_ycelos6, 512, _yh1fkj4)    != NULL)      {   sscanf(_ycelos6,"%*s %hu",&_vxpj4hn->_qdbxvh7);           }
        if (fgets(_ycelos6, 512, _yh1fkj4)    != NULL)      {   sscanf(_ycelos6,"%*s %ld",&_atj9tze->_0cdk36m);            }
        if (fgets(_ycelos6, 512, _yh1fkj4)    != NULL)      {   sscanf(_ycelos6,"%*s %u", &_vxpj4hn->_e6zvlhy);         }



        if (fgets(_ycelos6, 512, _yh1fkj4)    != NULL)      {                                                           }
        if (fgets(_ycelos6, 512, _yh1fkj4)    != NULL)      {   sscanf(_ycelos6,"%*s %hu",&_vxpj4hn->_wk8temp);        }
        if (fgets(_ycelos6, 512, _yh1fkj4)    != NULL)      {   sscanf(_ycelos6,"%*s %f", &_vxpj4hn->_rdihsa7);        }
        if (fgets(_ycelos6, 512, _yh1fkj4)    != NULL)      {                                                           }
        if (fgets(_ycelos6, 512, _yh1fkj4)    != NULL)      {   sscanf(_ycelos6,"%*s %f", &_vxpj4hn->_lw2lm34);         }

        if (fgets(_ycelos6, 512, _yh1fkj4)    != NULL)      {                                                           }
        if (fgets(_ycelos6, 512, _yh1fkj4)    != NULL)      {   sscanf(_ycelos6,"%*s %f", &_vxpj4hn->_3jczucq);        }
        if (fgets(_ycelos6, 512, _yh1fkj4)    != NULL)      {   sscanf(_ycelos6,"%*s %f", &_vxpj4hn->_x5qjicj);       }


        if (fgets(_ycelos6, 512, _yh1fkj4)    != NULL)      {                                                           }
        if (fgets(_ycelos6, 512, _yh1fkj4)    != NULL)      {   sscanf(_ycelos6,"%*s %lf",&_atj9tze->_0ur5l8a);         }
        
        fclose(_yh1fkj4);
    }
    

    _s55j2o7->_jcmrsez       = 1;
    _s55j2o7->_qy8yb7w      = 0;
    _s55j2o7->_o7vbgr8       = 0;
    _vxpj4hn->_e1g983e      = 20.0f; 
    _s55j2o7->_dhgzz0j = 0.85f; 
    _vxpj4hn->_7odkef8    = 0.01f; 

    _atj9tze->_0cdk36m        += (long)_fsqg8x8;
    _atj9tze->_14ve92m    = 10000.0f; 
    _atj9tze->_rtj5pht        = 5.0E+6f; 
    _atj9tze->_dswm31y       = 0.75f; 
    _atj9tze->_at9mx2x     = 20000u; 
    _s55j2o7->_s832ejx  = 1u;
    _s55j2o7->_a11dpg8   = 12.0f; 
    _s55j2o7->_sbi4q2y   = 36.0f; 
    _s55j2o7->_9v86wzs    =  0.3f;
    
    
    strcpy(_u1g463s, _dygflij);     snprintf(_75aud9k,sizeof(_75aud9k),"_%u.rgh",_s55j2o7->_jcmrsez);   strcat(_u1g463s,_75aud9k);
    
    if ((_yk1b3d5 = fopen(_u1g463s,"rb"))   == NULL)     {    perror("Error:cant open *Roughn.dat file in function 'LoadInputParaFile'\n");  exit(EXIT_FAILURE);                 }
    {   
        if (fread(_8p74mjd,            sizeof(unsigned int),       1, _yk1b3d5) != 1)       {   exit(EXIT_FAILURE);   }
        if (_8p74mjd[0] != VersionNumber)                {   perror("VersionNumber of code and parameter file are not identical. This might cause undefined behavior.\n");       }
        if (fread(&_vxpj4hn->_3ru7myo,  sizeof(unsigned int),       1, _yk1b3d5) != 1)       {   exit(EXIT_FAILURE);   }
        if (fread(&_s55j2o7->_fu7bfxq,  sizeof(unsigned int),       1, _yk1b3d5) != 1)       {   exit(EXIT_FAILURE);   }
        if (fread(&_atj9tze->_cmktlzk,   sizeof(float),              1, _yk1b3d5) != 1)       {   exit(EXIT_FAILURE);   }
        if (fread(&_atj9tze->_iwvhdk8,sizeof(float),              1, _yk1b3d5) != 1)       {   exit(EXIT_FAILURE);   }
        if (fread(&_atj9tze->_kkrzd5e, sizeof(float),              1, _yk1b3d5) != 1)       {   exit(EXIT_FAILURE);   }
        if (fread(&_atj9tze->_wrb7f6f, sizeof(float),              1, _yk1b3d5) != 1)       {   exit(EXIT_FAILURE);   }
        if (fread(&_vxpj4hn->_pwh0mx1,sizeof(unsigned short int), 1, _yk1b3d5) != 1)       {   exit(EXIT_FAILURE);   }
        
        fclose(_yk1b3d5);
    }
    _vxpj4hn->_3ru7myo = _ktsg7fx(_vxpj4hn->_3ru7myo);
    
    strcpy(_kzmpjrh, _dygflij);     strcat(_kzmpjrh, ".btrg");
    
    if ((_050ui22 = fopen(_kzmpjrh,"rb"))   != NULL)
    {   if (fread(_8p74mjd,           sizeof(unsigned int),   2, _050ui22) != 2)       {   exit(EXIT_FAILURE);   }
        if (_8p74mjd[0] != VersionNumber)                {   perror("VersionNumber of code and parameter file are not identical. This might cause undefined behavior.\n");       }
        if (fread(&_vxpj4hn->_5izolnn, sizeof(unsigned int),   1, _050ui22) != 1)       {   exit(EXIT_FAILURE);   }
        if (fread(&_s55j2o7->_ryf0992, sizeof(unsigned int),   1, _050ui22) != 1)       {   exit(EXIT_FAILURE);   }
        fclose(_050ui22);
    }
    
    return;
}

