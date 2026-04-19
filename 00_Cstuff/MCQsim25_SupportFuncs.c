#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <mpi.h>
#include <limits.h>
#define VersionNumber 202511


#define IA 16807 
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0 - EPS)
#define MASK 123459876

#define ALIGN64 64 

#define MAX(A,B) ( ((A) > (B))  *(A)  +  ((A) <= (B))  *(B)    )
#define MIN(A,B) ( ((A) < (B))  *(A)  +  ((A) >= (B))  *(B)    )


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




long int _z0yncxl(int _uuijnwp, int _qtkcd0r, int _tgxnx4s)
{   long int _f1t2tal;
    long int _a753lc5;
    _a753lc5    = _qtkcd0r *_tgxnx4s;
    _f1t2tal = (_a753lc5%_uuijnwp == 0)*_a753lc5  +  (_a753lc5%_uuijnwp != 0)*((long int)(_a753lc5/_uuijnwp) +1)*_uuijnwp;

    return _f1t2tal;
}


void  _id4n1u3(float *restrict _j49yd82, const float *restrict _28xg8y7, const float *restrict _kf1622v)
{   _j49yd82[0] = _28xg8y7[0] - _kf1622v[0];
    _j49yd82[1] = _28xg8y7[1] - _kf1622v[1];
    _j49yd82[2] = _28xg8y7[2] - _kf1622v[2];
    return;
}


void  _j2qp89u(float *restrict _iquwrxz, const float *restrict _f445cri, const float *restrict _8bixoiq)
{   _iquwrxz[0] = _f445cri[1]*_8bixoiq[2] - _f445cri[2]*_8bixoiq[1];
    _iquwrxz[1] = _f445cri[2]*_8bixoiq[0] - _f445cri[0]*_8bixoiq[2];
    _iquwrxz[2] = _f445cri[0]*_8bixoiq[1] - _f445cri[1]*_8bixoiq[0];
    return;
}


float _mabsi9n(const float *restrict _n42x2ez)
{   float _2q4qerx = 0.0f;
    _2q4qerx = sqrtf(_n42x2ez[0]*_n42x2ez[0] + _n42x2ez[1]*_n42x2ez[1] + _n42x2ez[2]*_n42x2ez[2]);
    return _2q4qerx;
}


void    _bfuk9e1(float *restrict _n42x2ez)
{   float _2q4qerx = 0.0f;
    _2q4qerx = sqrtf(_n42x2ez[0]*_n42x2ez[0] + _n42x2ez[1]*_n42x2ez[1] + _n42x2ez[2]*_n42x2ez[2]);
    _n42x2ez[0] /= _2q4qerx;              _n42x2ez[1] /= _2q4qerx;              _n42x2ez[2] /= _2q4qerx; 
    return;
}


void _ktakpux(const float *restrict _088o98b, float *restrict _j3swebn, float *restrict _795e5iz)
{   
    float _2q4qerx, _74rs7gu[3];
    _74rs7gu[0] = 0.0f;                          _74rs7gu[1] = 0.0f;                              _74rs7gu[2] = 1.0f;
    _j2qp89u(_j3swebn, _74rs7gu, _088o98b);          _2q4qerx = _mabsi9n(_j3swebn);
    
    if (_2q4qerx <= FLT_EPSILON)
    {   _j3swebn[0]    = 0.0f;                  _j3swebn[1] = 1.0f;                             _j3swebn[2] = 0.0f;                   }
    else
    {   _j3swebn[0]   /= _2q4qerx;             _j3swebn[1]/= _2q4qerx;                        _j3swebn[2]/= _2q4qerx;              }
    
    _j2qp89u(_795e5iz, _088o98b, _j3swebn);         _bfuk9e1(_795e5iz);
    
    return;
}


void _ftu2gc5(float *restrict _0i9g5dz, const float *restrict _mt8vt43, const float *restrict _xyzx9jp)
{   _0i9g5dz[2] = (_mt8vt43[0]*_xyzx9jp[0] + _mt8vt43[1]*_xyzx9jp[1] + _mt8vt43[2]*_xyzx9jp[2])*_mt8vt43[0] + (_mt8vt43[0]*_xyzx9jp[1] + _mt8vt43[1]*_xyzx9jp[3] + _mt8vt43[2]*_xyzx9jp[4])*_mt8vt43[1] + (_mt8vt43[0]*_xyzx9jp[2] + _mt8vt43[1]*_xyzx9jp[4] + _mt8vt43[2]*_xyzx9jp[5])*_mt8vt43[2];
    _0i9g5dz[0] = (_mt8vt43[3]*_xyzx9jp[0] + _mt8vt43[4]*_xyzx9jp[1] + _mt8vt43[5]*_xyzx9jp[2])*_mt8vt43[0] + (_mt8vt43[3]*_xyzx9jp[1] + _mt8vt43[4]*_xyzx9jp[3] + _mt8vt43[5]*_xyzx9jp[4])*_mt8vt43[1] + (_mt8vt43[3]*_xyzx9jp[2] + _mt8vt43[4]*_xyzx9jp[4] + _mt8vt43[5]*_xyzx9jp[5])*_mt8vt43[2];
    _0i9g5dz[1] = (_mt8vt43[6]*_xyzx9jp[0] + _mt8vt43[7]*_xyzx9jp[1] + _mt8vt43[8]*_xyzx9jp[2])*_mt8vt43[0] + (_mt8vt43[6]*_xyzx9jp[1] + _mt8vt43[7]*_xyzx9jp[3] + _mt8vt43[8]*_xyzx9jp[4])*_mt8vt43[1] + (_mt8vt43[6]*_xyzx9jp[2] + _mt8vt43[7]*_xyzx9jp[4] + _mt8vt43[8]*_xyzx9jp[5])*_mt8vt43[2];
    return;
}


void   _yyv37di(float *restrict _0yknfu8, float _ktmonrv)
{   _0yknfu8[0] *= _ktmonrv;           _0yknfu8[1] *= _ktmonrv;           _0yknfu8[2] *= _ktmonrv;
    _0yknfu8[3] *= _ktmonrv;           _0yknfu8[4] *= _ktmonrv;           _0yknfu8[5] *= _ktmonrv;
    return;
}


float _53ayeqa(long *_rrhf4ab)
{   long k;
    float _764jwte;
    
    *_rrhf4ab ^= MASK;
    k = (*_rrhf4ab)/IQ;
    *_rrhf4ab = IA*(*_rrhf4ab-k*IQ) -IR*k;
    if (*_rrhf4ab < 0) *_rrhf4ab += IM;
    _764jwte = AM*(*_rrhf4ab);
    *_rrhf4ab ^= MASK;

    return _764jwte;
}


float _hthsjr0(long *_rrhf4ab)
{   int j;
    long k;
    static long _k5ogz2u = 0;
    static long _1zderi4[NTAB];
    float _6parsg7;

    if (*_rrhf4ab <= 0 || !_k5ogz2u )
    {   if (-(*_rrhf4ab) < 1)   {   *_rrhf4ab = 1;          }
        else                {   *_rrhf4ab = -(*_rrhf4ab);   }
        for (j = NTAB +7; j >= 0; j--)
        {   k = (*_rrhf4ab)/IQ;
            *_rrhf4ab = IA*(*_rrhf4ab -k*IQ) -IR*k;
            if (*_rrhf4ab < 0)  {   *_rrhf4ab += IM;        }
            if (j < NTAB)   {   _1zderi4[j]  = *_rrhf4ab;     }
        }
        _k5ogz2u = _1zderi4[0];
    }

    k     = (*_rrhf4ab)/IQ;
    *_rrhf4ab = IA*(*_rrhf4ab -k*IQ) -IR*k;
    if (*_rrhf4ab < 0)  {   *_rrhf4ab += IM;        }
    j     = _k5ogz2u/NDIV;
    _k5ogz2u    = _1zderi4[j];
    _1zderi4[j] = *_rrhf4ab;
    if ((_6parsg7 = AM*_k5ogz2u) > RNMX)  {   return RNMX;        }
    else                        {   return _6parsg7;        }
}




void _myx444n(struct _5r5uqpg *restrict _atj9tze, const struct _xjuy0bv *restrict _vxpj4hn, unsigned int i, unsigned int *restrict _p422ta8, float *restrict _6fokzix, float *restrict _0dl5e11, unsigned int _urvqt3g)
{   float _nbzug0i[4];
    
    _nbzug0i[0]      = _53ayeqa(&_atj9tze->_0cdk36m);      _nbzug0i[0] = _nbzug0i[0]*2.0f -1.0f;
    _0dl5e11[i*16 +2] = _0dl5e11[i*16 +8] + (_0dl5e11[i*16 +9] * _nbzug0i[0]); 
    _nbzug0i[0]      = _53ayeqa(&_atj9tze->_0cdk36m);      _nbzug0i[0] = _nbzug0i[0]*2.0f -1.0f;
    _0dl5e11[i*16 +3] = MAX(0.0f, (_0dl5e11[i*16 +10] + (_0dl5e11[i*16 +11] * _nbzug0i[0])) ); 
    _nbzug0i[0]      = _53ayeqa(&_atj9tze->_0cdk36m);      _nbzug0i[0] = _nbzug0i[0]*2.0f -1.0f;
    _6fokzix[i*16 +4] = _0dl5e11[i*16 +12] + (_0dl5e11[i*16 +13] * _nbzug0i[0]); 
    
    _6fokzix[i*16 +3] = (_urvqt3g == 1u)*_0dl5e11[i*16 +2]  +  (_urvqt3g != 1u)*_6fokzix[i*16 +3]; 
    
    _nbzug0i[0]      = (_0dl5e11[i*16 +2] - _0dl5e11[i*16 +3]) *-1.0f*_0dl5e11[i*16 +0]; 
    _nbzug0i[1]      = _nbzug0i[0]/_0dl5e11[i*16 +1];
    _p422ta8[i*4 +0]  = (_nbzug0i[0]         >= 0.0f)*2u  +  (_nbzug0i[0]           < 0.0f)*3u;
    _p422ta8[i*4 +0]  = (_nbzug0i[1] > _6fokzix[i*16 +4])*1u  +  (_nbzug0i[1] <= _6fokzix[i*16 +4])*_p422ta8[i*4 +0];
    
    _nbzug0i[0]      = _0dl5e11[i*16 +3]*-1.0f*_0dl5e11[i*16 +0] + _6fokzix[i*16 +4]*_0dl5e11[i*16 +1]; 
    _nbzug0i[1]      = _0dl5e11[i*16 +2]*-1.0f*_0dl5e11[i*16 +0]; 
    _nbzug0i[2]      = MAX(_nbzug0i[0], _nbzug0i[1])/(-1.0f*_0dl5e11[i*16 +0]); 
    _0dl5e11[i*16 +4] = MAX(0.0f, (_nbzug0i[2] -  _vxpj4hn->_x5qjicj*(_nbzug0i[2]-_0dl5e11[i*16 +3])) ); 
    
    return;
}


void _790pdw5(int _fsqg8x8, const int *restrict _lexdw6l, struct _5r5uqpg *restrict _atj9tze, unsigned int *restrict _p422ta8, float *restrict _6fokzix, float *restrict _0dl5e11, float *restrict _z39pecn)
{   int i;
    unsigned int _xkc5i8x;
    float _nbzug0i[4];
    for (i = 0u; i < _lexdw6l[_fsqg8x8]; i++)
    {   _nbzug0i[0]        = sqrtf(_6fokzix[i*16 +13]*_6fokzix[i*16 +13] +_6fokzix[i*16 +14]*_6fokzix[i*16 +14]);
        _z39pecn[i*8 +3] += _nbzug0i[0];
        _z39pecn[i*8 +4] += _nbzug0i[0];
        _6fokzix[i*16 +2]   = _0dl5e11[i*16 +0]; 
        
        _nbzug0i[0]        = sqrtf(_6fokzix[i*16 +0]*_6fokzix[i*16 +0] + _6fokzix[i*16 +1]*_6fokzix[i*16 +1]); 
        _nbzug0i[1]        = _0dl5e11[i*16 +2] *-1.0f*_0dl5e11[i*16 +0]; 
        _xkc5i8x           = (_nbzug0i[0] > _nbzug0i[1])*1u                          +  (_nbzug0i[0]<= _nbzug0i[1])*0u; 
        _nbzug0i[3]        = (_xkc5i8x  == 1u)*(_nbzug0i[0]/(-1.0f*_6fokzix[i*16 +2]))   +  (_xkc5i8x != 1u)*_0dl5e11[i*16 +2]; 
        _6fokzix[i*16 +3]   = MIN(_nbzug0i[3], _atj9tze->_dswm31y);
        
        _nbzug0i[1]        = _6fokzix[i*16 +3] *-1.0f*_6fokzix[i*16 +2]; 
        _nbzug0i[2]        = MAX((_nbzug0i[0]/_nbzug0i[1]), 1.0f); 
        _6fokzix[i*16 +0]  /= _nbzug0i[2]; 
        _6fokzix[i*16 +1]  /= _nbzug0i[2];
        
        _0dl5e11[i*16 +7]   = _6fokzix[i*16 +3] - _0dl5e11[i*16 +2]; 
        
        _p422ta8[i*4 +1]    = 0u;          _p422ta8[i*4 +2]   = 0u;         _p422ta8[i*4 +3]   = 0u;
        _6fokzix[i*16 +11]  = 0.0f;        _6fokzix[i*16 +12] = 0.0f;       _6fokzix[i*16 +13] = 0.0f;       _6fokzix[i*16 +14] = 0.0f;       _6fokzix[i*16 +15] = 0.0f;
        
    }
    return;
}


void   _4lmct64(int _fsqg8x8, int _aytxw1v, int _zc7d2kd, const char *restrict _zle6kwh, unsigned int _6zjq7w7, const int *restrict _lexdw6l, const unsigned int *restrict _1nlc8jp, const unsigned int *restrict _vmefozt, const float *restrict _4xkss0k, const float *restrict _zzf8s2z, const float *restrict _0dl5e11, const float *restrict _6fokzix, const float *restrict _z39pecn)
{   unsigned int i;
    float  *_5d3go8j  = malloc(_6zjq7w7 *sizeof (float));             memset(_5d3go8j, 0, _6zjq7w7 *sizeof (float));
    
    FILE   *_yh1fkj4;
    if ((_yh1fkj4 = fopen(_zle6kwh,"w")) == NULL)      {   fprintf(stdout,"Error -cant open %s file in MCQsim25_SmallFuncs.c function \n",_zle6kwh);      exit(EXIT_FAILURE);     }
    
    if (_fsqg8x8 == 0)
    {   fprintf(_yh1fkj4,"FromCodeVersion#:     %d\n", VersionNumber);
        if      (_aytxw1v == 1)                                                 {       fprintf(_yh1fkj4,"output written for:   uF_t \n");                                   }
        else if (_aytxw1v == 2)                                                 {       fprintf(_yh1fkj4,"output written for:   fF_t \n");                                   }
        else if (_aytxw1v == 3)                                                 {       fprintf(_yh1fkj4,"output written for:   fFRF \n");                                   }
        else if (_aytxw1v == 4)                                                 {       fprintf(_yh1fkj4,"output written for:   fFEv \n");                                   }
        else if (_aytxw1v == 5)                                                 {       fprintf(_yh1fkj4,"output written for:   fFLcPl \n");                                 }
        else                                                                        {       fprintf(stdout,"incorrect array selection\n");        exit(EXIT_FAILURE);       }
        if    ( (_aytxw1v == 1) && ((_zc7d2kd < 0) || (_zc7d2kd >  5)) )  {       fprintf(stdout,"incorrect array entry selection\n");  exit(EXIT_FAILURE);       }
        if    ( (_aytxw1v == 2) && ((_zc7d2kd < 0) || (_zc7d2kd > 15)) )  {       fprintf(stdout,"incorrect array entry selection\n");  exit(EXIT_FAILURE);       }
        if    ( (_aytxw1v == 3) && ((_zc7d2kd < 0) || (_zc7d2kd > 13)) )  {       fprintf(stdout,"incorrect array entry selection\n");  exit(EXIT_FAILURE);       }
        if    ( (_aytxw1v == 4) && ((_zc7d2kd < 0) || (_zc7d2kd > 15)) )  {       fprintf(stdout,"incorrect array entry selection\n");  exit(EXIT_FAILURE);       }
        if    ( (_aytxw1v == 5) && ((_zc7d2kd < 0) || (_zc7d2kd >  5)) )  {       fprintf(stdout,"incorrect array entry selection\n");  exit(EXIT_FAILURE);       }
    }
    
    if  (_aytxw1v == 3)
    {   for (i = 0u; i < _lexdw6l[_fsqg8x8]; i++)                                      {       _5d3go8j[_1nlc8jp[i]]  = _0dl5e11[i*16 +_zc7d2kd];                      }
        MPI_Reduce(MPI_IN_PLACE,  _5d3go8j,  _6zjq7w7,  MPI_FLOAT,  MPI_SUM,  0, MPI_COMM_WORLD);
    }
    else if (_aytxw1v == 4)
    {   for (i = 0u; i < _lexdw6l[_fsqg8x8]; i++)                                      {       _5d3go8j[_1nlc8jp[i]]  = _6fokzix[i*16 +_zc7d2kd];                      }
        MPI_Reduce(MPI_IN_PLACE,  _5d3go8j,  _6zjq7w7,  MPI_FLOAT,  MPI_SUM,  0, MPI_COMM_WORLD);
    }
    else if (_aytxw1v == 5)
    {   for (i = 0u; i < _lexdw6l[_fsqg8x8]; i++)                                      {       _5d3go8j[_1nlc8jp[i]]  = _z39pecn[i*8 +_zc7d2kd];                     }
        MPI_Reduce(MPI_IN_PLACE,  _5d3go8j,  _6zjq7w7,  MPI_FLOAT,  MPI_SUM,  0, MPI_COMM_WORLD);
    }
    
    if (_fsqg8x8 == 0)
    {   
        fprintf(_yh1fkj4,"ArrayEntryPosition:   %d\n", _zc7d2kd);
        fprintf(_yh1fkj4,"ElementNumber:        %d\n", _6zjq7w7);
        
        if       (_aytxw1v == 1)
        {   for (i = 0u; i < _6zjq7w7; i++)
            {   fprintf(_yh1fkj4,"ElemID:     %d\n",i);
                fprintf(_yh1fkj4,"ElemValue:  %d\n",_vmefozt[i*8 +_zc7d2kd]);
                fprintf(_yh1fkj4,"Vert1:      %e    %e    %e\n", _zzf8s2z[_vmefozt[i*8 +0] *4 +0], _zzf8s2z[_vmefozt[i*8 +0] *4 +1], _zzf8s2z[_vmefozt[i*8 +0] *4 +2]);
                fprintf(_yh1fkj4,"Vert2:      %e    %e    %e\n", _zzf8s2z[_vmefozt[i*8 +1] *4 +0], _zzf8s2z[_vmefozt[i*8 +1] *4 +1], _zzf8s2z[_vmefozt[i*8 +1] *4 +2]);
                fprintf(_yh1fkj4,"Vert3:      %e    %e    %e\n", _zzf8s2z[_vmefozt[i*8 +2] *4 +0], _zzf8s2z[_vmefozt[i*8 +2] *4 +1], _zzf8s2z[_vmefozt[i*8 +2] *4 +2]);
        }   }
        else if  (_aytxw1v == 2)
        {   for (i = 0u; i < _6zjq7w7; i++)
            {   fprintf(_yh1fkj4,"ElemID:     %d\n",i);
                fprintf(_yh1fkj4,"ElemValue:  %e\n",_4xkss0k[i*16 +_zc7d2kd]);
                fprintf(_yh1fkj4,"Vert1:      %e    %e    %e\n", _zzf8s2z[_vmefozt[i*8 +0] *4 +0], _zzf8s2z[_vmefozt[i*8 +0] *4 +1], _zzf8s2z[_vmefozt[i*8 +0] *4 +2]);
                fprintf(_yh1fkj4,"Vert2:      %e    %e    %e\n", _zzf8s2z[_vmefozt[i*8 +1] *4 +0], _zzf8s2z[_vmefozt[i*8 +1] *4 +1], _zzf8s2z[_vmefozt[i*8 +1] *4 +2]);
                fprintf(_yh1fkj4,"Vert3:      %e    %e    %e\n", _zzf8s2z[_vmefozt[i*8 +2] *4 +0], _zzf8s2z[_vmefozt[i*8 +2] *4 +1], _zzf8s2z[_vmefozt[i*8 +2] *4 +2]);
        }   }
        else if  (_aytxw1v >= 3)
        {   for (i = 0u; i < _6zjq7w7; i++)
            {   fprintf(_yh1fkj4,"ElemID:     %d\n",i);
                fprintf(_yh1fkj4,"ElemValue:  %e\n",_5d3go8j[i]);
                fprintf(_yh1fkj4,"Vert1:      %e    %e    %e\n", _zzf8s2z[_vmefozt[i*8 +0] *4 +0], _zzf8s2z[_vmefozt[i*8 +0] *4 +1], _zzf8s2z[_vmefozt[i*8 +0] *4 +2]);
                fprintf(_yh1fkj4,"Vert2:      %e    %e    %e\n", _zzf8s2z[_vmefozt[i*8 +1] *4 +0], _zzf8s2z[_vmefozt[i*8 +1] *4 +1], _zzf8s2z[_vmefozt[i*8 +1] *4 +2]);
                fprintf(_yh1fkj4,"Vert3:      %e    %e    %e\n", _zzf8s2z[_vmefozt[i*8 +2] *4 +0], _zzf8s2z[_vmefozt[i*8 +2] *4 +1], _zzf8s2z[_vmefozt[i*8 +2] *4 +2]);
    }   }   }
    
    free(_5d3go8j);
    fclose(_yh1fkj4);
    return;
}


void _z7w2seh(int _fsqg8x8, const char *restrict _zle6kwh, const struct _xjuy0bv *restrict _vxpj4hn, const struct _j49itn7 *restrict _s55j2o7, const int *restrict _lexdw6l, const unsigned int *restrict _1nlc8jp, const unsigned int *restrict _vmefozt, const float *restrict _4xkss0k, const float *restrict _zzf8s2z, const unsigned int *restrict _p422ta8, float *restrict _6fokzix, const float *restrict _0dl5e11)
{   unsigned int i;
    unsigned int *_9ufu6py  = malloc(_vxpj4hn->_3ru7myo *sizeof(unsigned int));    memset(_9ufu6py, 0, _vxpj4hn->_3ru7myo *sizeof(unsigned int));
    float *_4yfds0d        = malloc(_vxpj4hn->_3ru7myo *sizeof(float));           memset(_4yfds0d,0, _vxpj4hn->_3ru7myo *sizeof(float));
    float *_g87oi0e         = malloc(_vxpj4hn->_3ru7myo *sizeof(float));           memset(_g87oi0e, 0, _vxpj4hn->_3ru7myo *sizeof(float));
    float *_jtb9xxa         = malloc(_vxpj4hn->_3ru7myo *sizeof(float));           memset(_jtb9xxa, 0, _vxpj4hn->_3ru7myo *sizeof(float));
    float *_ypntrfz         = malloc(_vxpj4hn->_3ru7myo *sizeof(float));           memset(_ypntrfz, 0, _vxpj4hn->_3ru7myo *sizeof(float));
    float *_wkxxjve         = malloc(_vxpj4hn->_3ru7myo *sizeof(float));           memset(_wkxxjve, 0, _vxpj4hn->_3ru7myo *sizeof(float));
    float *_3dzeedh         = malloc(_vxpj4hn->_3ru7myo *sizeof(float));           memset(_3dzeedh, 0, _vxpj4hn->_3ru7myo *sizeof(float));
    float *_hjshzkb         = malloc(_vxpj4hn->_3ru7myo *sizeof(float));           memset(_hjshzkb, 0, _vxpj4hn->_3ru7myo *sizeof(float));
    
    for (i = 0u; i < _lexdw6l[_fsqg8x8]; i++)
    {   _4yfds0d[_1nlc8jp[i]]= _6fokzix[i*16 +11];                      _g87oi0e[_1nlc8jp[i]] = _6fokzix[i*16 +12];
        _jtb9xxa[_1nlc8jp[i]] = _6fokzix[i*16 +9];                       _ypntrfz[_1nlc8jp[i]] = _6fokzix[i*16 +10];
        _wkxxjve[_1nlc8jp[i]] = _0dl5e11[i*16 +2] - _0dl5e11[i*16 +3];       _3dzeedh[_1nlc8jp[i]] = (_0dl5e11[i*16 +2] - _0dl5e11[i*16 +3])* -1.0*_0dl5e11[i*16 +0]; 
        _hjshzkb[_1nlc8jp[i]] = _6fokzix[i*16 +4]; 
        _9ufu6py[_1nlc8jp[i]] = _p422ta8[i*4 +0]; 
        
        _6fokzix[i*16 +0]  = 0.0f;          _6fokzix[i*16 +1]  = 0.0f;          _6fokzix[i*16 +11] = 0.0f;          _6fokzix[i*16 +12] = 0.0f;
        _6fokzix[i*16 +13] = 0.0f;          _6fokzix[i*16 +14] = 0.0f;
    }
    MPI_Allreduce(MPI_IN_PLACE, _4yfds0d, _vxpj4hn->_3ru7myo, MPI_FLOAT,    MPI_SUM,  MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, _g87oi0e,  _vxpj4hn->_3ru7myo, MPI_FLOAT,    MPI_SUM,  MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, _jtb9xxa,  _vxpj4hn->_3ru7myo, MPI_FLOAT,    MPI_SUM,  MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, _ypntrfz,  _vxpj4hn->_3ru7myo, MPI_FLOAT,    MPI_SUM,  MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, _wkxxjve,  _vxpj4hn->_3ru7myo, MPI_FLOAT,    MPI_SUM,  MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, _3dzeedh,  _vxpj4hn->_3ru7myo, MPI_FLOAT,    MPI_SUM,  MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, _hjshzkb,  _vxpj4hn->_3ru7myo, MPI_FLOAT,    MPI_SUM,  MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, _9ufu6py,  _vxpj4hn->_3ru7myo, MPI_UNSIGNED, MPI_SUM,  MPI_COMM_WORLD);
    
    if (_fsqg8x8 == 0)
    {   unsigned int *_09lzrk7 = malloc(_vxpj4hn->_3ru7myo *sizeof(unsigned int));      memset(_09lzrk7,  0, _vxpj4hn->_3ru7myo *sizeof(unsigned int));
        float *_rh7lay6        = malloc(_s55j2o7->_fu7bfxq *sizeof(float));             memset(_rh7lay6,  0, _s55j2o7->_fu7bfxq *sizeof(float));
        float *_0vgnxb5       = malloc(_vxpj4hn->_3ru7myo *sizeof(float));             memset(_0vgnxb5, 0, _vxpj4hn->_3ru7myo *sizeof(float));
        
        FILE  *_kx7n89l;
        char _efqjg02[512],   _75aud9k[512];
        strcpy(_efqjg02, _zle6kwh);      strcat(_efqjg02,"_");      snprintf(_75aud9k,sizeof(_75aud9k), "%u",_s55j2o7->_jcmrsez);     strcat(_efqjg02,_75aud9k);     strcat(_efqjg02,".pre");
        
        if ((_kx7n89l = fopen(_efqjg02,"wb")) == NULL)        {   printf("Error -cant open %s  PrePostFile...\n", _efqjg02);      exit(EXIT_FAILURE);     }
        
        int _3zezqac = (int)VersionNumber;
        fwrite( &_3zezqac,      sizeof(int),          1, _kx7n89l);
        fwrite( &_vxpj4hn->_3ru7myo,  sizeof(unsigned int), 1, _kx7n89l);
        fwrite( &_s55j2o7->_fu7bfxq,  sizeof(unsigned int), 1, _kx7n89l);
        for (i = 0u; i < _vxpj4hn->_3ru7myo; i++)   {   _09lzrk7[i]  = _vmefozt[i*8 +0];   }     fwrite( _09lzrk7,     sizeof(unsigned int),  _vxpj4hn->_3ru7myo, _kx7n89l);
        for (i = 0u; i < _vxpj4hn->_3ru7myo; i++)   {   _09lzrk7[i]  = _vmefozt[i*8 +1];   }     fwrite( _09lzrk7,     sizeof(unsigned int),  _vxpj4hn->_3ru7myo, _kx7n89l);
        for (i = 0u; i < _vxpj4hn->_3ru7myo; i++)   {   _09lzrk7[i]  = _vmefozt[i*8 +2];   }     fwrite( _09lzrk7,     sizeof(unsigned int),  _vxpj4hn->_3ru7myo, _kx7n89l);
        for (i = 0u; i < _vxpj4hn->_3ru7myo; i++)   {   _09lzrk7[i]  = _vmefozt[i*8 +3];   }     fwrite( _09lzrk7,     sizeof(unsigned int),  _vxpj4hn->_3ru7myo, _kx7n89l);
        for (i = 0u; i < _vxpj4hn->_3ru7myo; i++)   {   _09lzrk7[i]  = _vmefozt[i*8 +4];   }     fwrite( _09lzrk7,     sizeof(unsigned int),  _vxpj4hn->_3ru7myo, _kx7n89l);

        for (i = 0u; i < _s55j2o7->_fu7bfxq; i++)   {   _rh7lay6[i]  = _zzf8s2z[i*4 +0];  }     fwrite( _rh7lay6,     sizeof(float),         _s55j2o7->_fu7bfxq, _kx7n89l);
        for (i = 0u; i < _s55j2o7->_fu7bfxq; i++)   {   _rh7lay6[i]  = _zzf8s2z[i*4 +1];  }     fwrite( _rh7lay6,     sizeof(float),         _s55j2o7->_fu7bfxq, _kx7n89l);
        for (i = 0u; i < _s55j2o7->_fu7bfxq; i++)   {   _rh7lay6[i]  = _zzf8s2z[i*4 +2];  }     fwrite( _rh7lay6,     sizeof(float),         _s55j2o7->_fu7bfxq, _kx7n89l);
        for (i = 0u; i < _s55j2o7->_fu7bfxq; i++)   {   _rh7lay6[i]  = _zzf8s2z[i*4 +3];  }     fwrite( _rh7lay6,     sizeof(float),         _s55j2o7->_fu7bfxq, _kx7n89l);

        for (i = 0u; i < _vxpj4hn->_3ru7myo; i++)   {   _0vgnxb5[i] = _4xkss0k[i*16 +0];  }     fwrite( _0vgnxb5,    sizeof(float),         _vxpj4hn->_3ru7myo, _kx7n89l);
        for (i = 0u; i < _vxpj4hn->_3ru7myo; i++)   {   _0vgnxb5[i] = _4xkss0k[i*16 +1];  }     fwrite( _0vgnxb5,    sizeof(float),         _vxpj4hn->_3ru7myo, _kx7n89l);
        for (i = 0u; i < _vxpj4hn->_3ru7myo; i++)   {   _0vgnxb5[i] = _4xkss0k[i*16 +2];  }     fwrite( _0vgnxb5,    sizeof(float),         _vxpj4hn->_3ru7myo, _kx7n89l);

        fwrite(_4yfds0d,sizeof(float), _vxpj4hn->_3ru7myo, _kx7n89l);
        fwrite(_g87oi0e, sizeof(float), _vxpj4hn->_3ru7myo, _kx7n89l);
        fwrite(_jtb9xxa, sizeof(float), _vxpj4hn->_3ru7myo, _kx7n89l);
        fwrite(_ypntrfz, sizeof(float), _vxpj4hn->_3ru7myo, _kx7n89l);
        fwrite(_wkxxjve, sizeof(float), _vxpj4hn->_3ru7myo, _kx7n89l); 
        fwrite(_3dzeedh, sizeof(float), _vxpj4hn->_3ru7myo, _kx7n89l); 
        fwrite(_hjshzkb, sizeof(float), _vxpj4hn->_3ru7myo, _kx7n89l); 
        fwrite(_9ufu6py, sizeof(unsigned int), _vxpj4hn->_3ru7myo, _kx7n89l);
        
        fclose(_kx7n89l); 
        
        free(_09lzrk7);   free(_rh7lay6);   free(_0vgnxb5);
    } 
    
    free(_9ufu6py);     free(_4yfds0d);   free(_g87oi0e);   free(_jtb9xxa);   free(_ypntrfz);   free(_wkxxjve);   free(_3dzeedh);   free(_hjshzkb);
    return;
}


void _qpki7p5(int _fsqg8x8,  const char *restrict _dygflij, struct _5r5uqpg *restrict _atj9tze,  const struct _xjuy0bv *restrict _vxpj4hn, const struct _j49itn7 *restrict _s55j2o7, unsigned int *restrict _xo100vh, double *restrict _ofn1r4k, const int *restrict _lexdw6l, const int *restrict _1poxrte, const unsigned int *restrict _1nlc8jp, const unsigned int *restrict _fftouhy, const unsigned int *restrict _vmefozt, const float *restrict _4xkss0k, const float *restrict _zzf8s2z, unsigned int *restrict _p422ta8, float *restrict _6fokzix, float *restrict _0dl5e11, float *restrict _z39pecn, float *restrict _7yt2m6t)
{   
    if (_s55j2o7->_qy8yb7w == 1u)
    {   
        MPI_Status  _0get8t6;
        MPI_Offset  _9kwv08h;
        MPI_File    _8kgyggs;
        
        int _n6sm5rr;
        unsigned int i;
        unsigned int *_p7bd4w1 = NULL;
        float *_pqvl00e = NULL,   *_q01tee1 = NULL;
        
        _p7bd4w1   = malloc(_vxpj4hn->_3ru7myo *sizeof(unsigned int));              memset(_p7bd4w1, 0, _vxpj4hn->_3ru7myo *sizeof(unsigned int));
        _pqvl00e   = malloc(_vxpj4hn->_3ru7myo *sizeof(float));                     memset(_pqvl00e, 0, _vxpj4hn->_3ru7myo *sizeof(float));
        if (_vxpj4hn->_5izolnn > 0u)
        {   _q01tee1   = malloc(_vxpj4hn->_5izolnn *sizeof(float));                 memset(_q01tee1, 0, _vxpj4hn->_5izolnn *sizeof(float));
        }
        
        char _efqjg02[512],   _75aud9k[512];
        strcpy(_efqjg02, _dygflij);      snprintf(_75aud9k,sizeof(_75aud9k), "_%u.post",_s55j2o7->_jcmrsez);        strcat(_efqjg02, _75aud9k);
        
        MPI_File_open(MPI_COMM_WORLD, _efqjg02, MPI_MODE_RDONLY, MPI_INFO_NULL, &_8kgyggs);
        

        MPI_File_read(_8kgyggs, &_n6sm5rr, 1, MPI_INT,            &_0get8t6);
        MPI_File_read(_8kgyggs, _ofn1r4k,   1, MPI_DOUBLE,         &_0get8t6);
        MPI_File_read(_8kgyggs, _xo100vh,      1, MPI_UNSIGNED,       &_0get8t6);
        MPI_File_read(_8kgyggs, _p7bd4w1,      2, MPI_UNSIGNED,       &_0get8t6);
        if (_n6sm5rr != VersionNumber)       {   perror("Version number of code and *.post file are not identical. This can cause undefined behavior\n");                            }
        if (_p7bd4w1[0] != _vxpj4hn->_3ru7myo)      {   exit(EXIT_FAILURE);     }
        if (_p7bd4w1[1] != _vxpj4hn->_5izolnn)      {   exit(EXIT_FAILURE);     }
        
        _9kwv08h = sizeof(int) + sizeof(double) +3*sizeof(unsigned int);
        
        MPI_File_read_at(_8kgyggs,   _9kwv08h, _p7bd4w1, _vxpj4hn->_3ru7myo, MPI_UNSIGNED, &_0get8t6);   _9kwv08h += (_vxpj4hn->_3ru7myo*sizeof(unsigned int));
        for (i = 0u; i < _lexdw6l[_fsqg8x8]; i++)      {   _p422ta8[i*4 +0] = _p7bd4w1[_1nlc8jp[i]];        }
        
        MPI_File_read_at(_8kgyggs,   _9kwv08h, _pqvl00e, _vxpj4hn->_3ru7myo, MPI_FLOAT, &_0get8t6);      _9kwv08h += (_vxpj4hn->_3ru7myo*sizeof(float));
        for (i = 0u; i < _lexdw6l[_fsqg8x8]; i++)      {   _6fokzix[i*16 +0] = _pqvl00e[_1nlc8jp[i]];       }
        MPI_File_read_at(_8kgyggs,   _9kwv08h, _pqvl00e, _vxpj4hn->_3ru7myo, MPI_FLOAT, &_0get8t6);      _9kwv08h += (_vxpj4hn->_3ru7myo*sizeof(float));
        for (i = 0u; i < _lexdw6l[_fsqg8x8]; i++)      {   _6fokzix[i*16 +1] = _pqvl00e[_1nlc8jp[i]];       }
        MPI_File_read_at(_8kgyggs,   _9kwv08h, _pqvl00e, _vxpj4hn->_3ru7myo, MPI_FLOAT, &_0get8t6);      _9kwv08h += (_vxpj4hn->_3ru7myo*sizeof(float));
        for (i = 0u; i < _lexdw6l[_fsqg8x8]; i++)      {   _6fokzix[i*16 +2] = _pqvl00e[_1nlc8jp[i]];       }
        MPI_File_read_at(_8kgyggs,   _9kwv08h, _pqvl00e, _vxpj4hn->_3ru7myo, MPI_FLOAT, &_0get8t6);      _9kwv08h += (_vxpj4hn->_3ru7myo*sizeof(float));
        for (i = 0u; i < _lexdw6l[_fsqg8x8]; i++)      {   _6fokzix[i*16 +3] = _pqvl00e[_1nlc8jp[i]];       }
        MPI_File_read_at(_8kgyggs,   _9kwv08h, _pqvl00e, _vxpj4hn->_3ru7myo, MPI_FLOAT, &_0get8t6);      _9kwv08h += (_vxpj4hn->_3ru7myo*sizeof(float));
        for (i = 0u; i < _lexdw6l[_fsqg8x8]; i++)      {   _6fokzix[i*16 +4] = _pqvl00e[_1nlc8jp[i]];       }
        
        MPI_File_read_at(_8kgyggs,   _9kwv08h, _pqvl00e, _vxpj4hn->_3ru7myo, MPI_FLOAT, &_0get8t6);      _9kwv08h += (_vxpj4hn->_3ru7myo*sizeof(float));
        for (i = 0u; i < _lexdw6l[_fsqg8x8]; i++)      {   _0dl5e11[i*16 +2] = _pqvl00e[_1nlc8jp[i]];       }
        MPI_File_read_at(_8kgyggs,   _9kwv08h, _pqvl00e, _vxpj4hn->_3ru7myo, MPI_FLOAT, &_0get8t6);      _9kwv08h += (_vxpj4hn->_3ru7myo*sizeof(float));
        for (i = 0u; i < _lexdw6l[_fsqg8x8]; i++)      {   _0dl5e11[i*16 +3] = _pqvl00e[_1nlc8jp[i]];       }
        MPI_File_read_at(_8kgyggs,   _9kwv08h, _pqvl00e, _vxpj4hn->_3ru7myo, MPI_FLOAT, &_0get8t6);      _9kwv08h += (_vxpj4hn->_3ru7myo*sizeof(float));
        for (i = 0u; i < _lexdw6l[_fsqg8x8]; i++)      {   _0dl5e11[i*16 +4] = _pqvl00e[_1nlc8jp[i]];       }
        MPI_File_read_at(_8kgyggs,   _9kwv08h, _pqvl00e, _vxpj4hn->_3ru7myo, MPI_FLOAT, &_0get8t6);      _9kwv08h += (_vxpj4hn->_3ru7myo*sizeof(float));
        for (i = 0u; i < _lexdw6l[_fsqg8x8]; i++)      {   _0dl5e11[i*16 +7] = _pqvl00e[_1nlc8jp[i]];       }
        
        MPI_File_read_at(_8kgyggs,   _9kwv08h, _pqvl00e, _vxpj4hn->_3ru7myo, MPI_FLOAT, &_0get8t6);      _9kwv08h += (_vxpj4hn->_3ru7myo*sizeof(float));
        for (i = 0u; i < _lexdw6l[_fsqg8x8]; i++)      { _z39pecn[i*8 +3] = _pqvl00e[_1nlc8jp[i]];        }
        MPI_File_read_at(_8kgyggs,   _9kwv08h, _pqvl00e, _vxpj4hn->_3ru7myo, MPI_FLOAT, &_0get8t6);      _9kwv08h += (_vxpj4hn->_3ru7myo*sizeof(float));
        for (i = 0u; i < _lexdw6l[_fsqg8x8]; i++)      { _z39pecn[i*8 +4] = _pqvl00e[_1nlc8jp[i]];        }
        
        if (_vxpj4hn->_5izolnn > 0u)
        {   MPI_File_read_at(_8kgyggs,   _9kwv08h, _q01tee1, _vxpj4hn->_5izolnn, MPI_FLOAT, &_0get8t6);      _9kwv08h += (_vxpj4hn->_5izolnn*sizeof(float));
            for (i = 0u; i < _1poxrte[_fsqg8x8]; i++)      {   _7yt2m6t[i*16 +0] = _q01tee1[_fftouhy[i]];       }
            MPI_File_read_at(_8kgyggs,   _9kwv08h, _q01tee1, _vxpj4hn->_5izolnn, MPI_FLOAT, &_0get8t6);      _9kwv08h += (_vxpj4hn->_5izolnn*sizeof(float));
            for (i = 0u; i < _1poxrte[_fsqg8x8]; i++)      {   _7yt2m6t[i*16 +1] = _q01tee1[_fftouhy[i]];       }
            MPI_File_read_at(_8kgyggs,   _9kwv08h, _q01tee1, _vxpj4hn->_5izolnn, MPI_FLOAT, &_0get8t6);      _9kwv08h += (_vxpj4hn->_5izolnn*sizeof(float));
            for (i = 0u; i < _1poxrte[_fsqg8x8]; i++)      {   _7yt2m6t[i*16 +2] = _q01tee1[_fftouhy[i]];       }
        }
        
        MPI_File_close(&_8kgyggs);
        MPI_Barrier( MPI_COMM_WORLD );
        
        free(_p7bd4w1);   free(_pqvl00e);   free(_q01tee1);
       
    } 
    else
    {   
        if (_fsqg8x8 == 0)
        {   FILE *_tajw253;
            
            unsigned int i;
            unsigned int *_7orfqla = malloc(_vxpj4hn->_3ru7myo *sizeof(unsigned int));                    memset(_7orfqla, 0, _vxpj4hn->_3ru7myo *sizeof(unsigned int));
            float *_vm59m6m        = malloc(_vxpj4hn->_3ru7myo *sizeof(float));                           memset(_vm59m6m, 0, _vxpj4hn->_3ru7myo *sizeof(float));
            float *_rh7lay6       = malloc(_s55j2o7->_fu7bfxq *sizeof(float));                           memset(_vm59m6m, 0, _s55j2o7->_fu7bfxq *sizeof(float));
            
            char _9tr74qg[512],    _75aud9k[512];
            strcpy(_9tr74qg,    _dygflij);      snprintf(_75aud9k,sizeof(_75aud9k), "_%u.cat",_s55j2o7->_jcmrsez);          strcat(_9tr74qg,    _75aud9k);
            
            if ((_tajw253 = fopen(_9tr74qg,"wb")) == NULL)     {   perror("Error -cannot open/write   RawCatalogFile...\n");      exit(EXIT_FAILURE);     }
            
            int _3zezqac = (int)VersionNumber;
            fwrite(_xo100vh,          sizeof(unsigned int), 1, _tajw253);
            fwrite(&_3zezqac,     sizeof(int),          1, _tajw253);
            fwrite(&_atj9tze->_kkrzd5e,sizeof(float),        1, _tajw253);
            fwrite(&_atj9tze->_mzwffhh, sizeof(float),        1, _tajw253);
            fwrite(&_vxpj4hn->_3ru7myo, sizeof(unsigned int), 1, _tajw253);
            fwrite(&_s55j2o7->_fu7bfxq, sizeof(unsigned int), 1, _tajw253);
            
            for (i = 0u; i < _vxpj4hn->_3ru7myo; i++)   {   _7orfqla[i] = _vmefozt[i*8 +0];   }     fwrite(_7orfqla,  sizeof(unsigned int), _vxpj4hn->_3ru7myo, _tajw253);
            for (i = 0u; i < _vxpj4hn->_3ru7myo; i++)   {   _7orfqla[i] = _vmefozt[i*8 +1];   }     fwrite(_7orfqla,  sizeof(unsigned int), _vxpj4hn->_3ru7myo, _tajw253);
            for (i = 0u; i < _vxpj4hn->_3ru7myo; i++)   {   _7orfqla[i] = _vmefozt[i*8 +2];   }     fwrite(_7orfqla,  sizeof(unsigned int), _vxpj4hn->_3ru7myo, _tajw253);
            for (i = 0u; i < _vxpj4hn->_3ru7myo; i++)   {   _7orfqla[i] = _vmefozt[i*8 +3];   }     fwrite(_7orfqla,  sizeof(unsigned int), _vxpj4hn->_3ru7myo, _tajw253);
            for (i = 0u; i < _vxpj4hn->_3ru7myo; i++)   {   _7orfqla[i] = _vmefozt[i*8 +4];   }     fwrite(_7orfqla,  sizeof(unsigned int), _vxpj4hn->_3ru7myo, _tajw253);
            
            for (i = 0u; i < _s55j2o7->_fu7bfxq; i++)   {   _rh7lay6[i]= _zzf8s2z[i*4 +0];  }     fwrite(_rh7lay6, sizeof(float),        _s55j2o7->_fu7bfxq, _tajw253); 
            for (i = 0u; i < _s55j2o7->_fu7bfxq; i++)   {   _rh7lay6[i]= _zzf8s2z[i*4 +1];  }     fwrite(_rh7lay6, sizeof(float),        _s55j2o7->_fu7bfxq, _tajw253);
            for (i = 0u; i < _s55j2o7->_fu7bfxq; i++)   {   _rh7lay6[i]= _zzf8s2z[i*4 +2];  }     fwrite(_rh7lay6, sizeof(float),        _s55j2o7->_fu7bfxq, _tajw253);
            for (i = 0u; i < _s55j2o7->_fu7bfxq; i++)   {   _rh7lay6[i]= _zzf8s2z[i*4 +3];  }     fwrite(_rh7lay6, sizeof(float),        _s55j2o7->_fu7bfxq, _tajw253);
            
            for (i = 0u; i < _vxpj4hn->_3ru7myo; i++)   {   _vm59m6m[i] = _4xkss0k[i*16 +0];   }    fwrite(_vm59m6m,  sizeof(float),        _vxpj4hn->_3ru7myo, _tajw253);
            for (i = 0u; i < _vxpj4hn->_3ru7myo; i++)   {   _vm59m6m[i] = _4xkss0k[i*16 +1];   }    fwrite(_vm59m6m,  sizeof(float),        _vxpj4hn->_3ru7myo, _tajw253);
            for (i = 0u; i < _vxpj4hn->_3ru7myo; i++)   {   _vm59m6m[i] = _4xkss0k[i*16 +2];   }    fwrite(_vm59m6m,  sizeof(float),        _vxpj4hn->_3ru7myo, _tajw253);
            
            free(_7orfqla);   free(_vm59m6m);   free(_rh7lay6);
            fclose(_tajw253);
        }
        
        unsigned int i;
        float _nbzug0i[3];
        for (i = 0u; i < _lexdw6l[_fsqg8x8]; i++)
        {   _nbzug0i[0]  = sqrtf(_6fokzix[i*16 +9]*_6fokzix[i*16 +9] + _6fokzix[i*16 +10]*_6fokzix[i*16 +10]); 
            if (_nbzug0i[0] > 0.0f)
            {   _nbzug0i[1]      = _6fokzix[i*16 +3]*-1.0f*_6fokzix[i*16 +2]; 
                _nbzug0i[2]      = 0.5f*(1.0f - _s55j2o7->_dhgzz0j)*(_53ayeqa(&_atj9tze->_0cdk36m)*2.0f -1.0f);
                _6fokzix[i*16 +0] = ((_nbzug0i[1]*(_s55j2o7->_dhgzz0j+_nbzug0i[2]))/_nbzug0i[0]) * _6fokzix[i*16 +9];
                _6fokzix[i*16 +1] = ((_nbzug0i[1]*(_s55j2o7->_dhgzz0j+_nbzug0i[2]))/_nbzug0i[0]) * _6fokzix[i*16 +10];
        }   } 
        
    }
    
    return;
}

