#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <mpi.h>
#define VersionNumber 202511

#define MAX(A,B) ( ((A) > (B))  *(A)  +  ((A) <= (B))  *(B)    )


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


void   _fqewaou(int _fsqg8x8, const char *restrict _dygflij, const struct _xjuy0bv *restrict _vxpj4hn, struct _j49itn7 *restrict _s55j2o7, const int *restrict _lexdw6l, const unsigned int *restrict _1nlc8jp, unsigned int *restrict _vmefozt,  unsigned int *restrict _rgb93wk, float *restrict _4xkss0k, float *restrict _rrpzlk3, float *restrict _zzf8s2z, float *restrict _3n3bzt9, float *restrict _0dl5e11)
{   
    unsigned int i,   _ytkw283;
    
    unsigned int *_7orfqla = NULL,   *_blmntfk = NULL;
    float *_vm59m6m = NULL,  *_ji9kqdg = NULL,    *_gywiyme = NULL,  *_nvg3kk0 = NULL;
    
    _7orfqla      = malloc(_vxpj4hn->_3ru7myo*sizeof(unsigned int));         memset(_7orfqla,  0, _vxpj4hn->_3ru7myo*sizeof(unsigned int));
    _vm59m6m      = malloc(_vxpj4hn->_3ru7myo*sizeof(float));                memset(_vm59m6m,  0, _vxpj4hn->_3ru7myo*sizeof(float));
    _ji9kqdg     = malloc(_s55j2o7->_fu7bfxq*sizeof(float));                memset(_ji9kqdg, 0, _s55j2o7->_fu7bfxq*sizeof(float));
    
    if (_vxpj4hn->_5izolnn > 0u)
    {   _blmntfk  = malloc(_vxpj4hn->_5izolnn*sizeof(unsigned int));         memset(_blmntfk,  0, _vxpj4hn->_5izolnn*sizeof(unsigned int));
        _gywiyme  = malloc(_vxpj4hn->_5izolnn*sizeof(float));                memset(_gywiyme,  0, _vxpj4hn->_5izolnn*sizeof(float));
        _nvg3kk0 = malloc(_s55j2o7->_ryf0992*sizeof(float));                memset(_nvg3kk0, 0, _s55j2o7->_ryf0992*sizeof(float));
    }
    
    char   _yn9382k[512],   _ytwbozd[512],   _1at7o75[512],   _rh8lx6d[512];
    FILE   *_yh1fkj4,           *_yk1b3d5,              *_050ui22;
    
    
    strcpy(_ytwbozd, _dygflij);      snprintf(_yn9382k,sizeof(_yn9382k),"_%u.rgh",_s55j2o7->_jcmrsez);     strcat(_ytwbozd,_yn9382k);
    strcpy(_1at7o75, _dygflij);      snprintf(_yn9382k,sizeof(_yn9382k),"_%u.stg",_s55j2o7->_jcmrsez);     strcat(_1at7o75,_yn9382k);
    strcpy(_rh8lx6d, _dygflij);      strcat(_rh8lx6d,  ".btrg");
    
    if ((_yh1fkj4 = fopen(_ytwbozd,"rb"))   == NULL)    {   perror("Error: -cant open *.rgh file in function 'LoadRghStgthBnd'\n");      exit(EXIT_FAILURE);      }
    {   
        if (fread(_7orfqla, sizeof(unsigned int), 1, _yh1fkj4) != 1)       {   exit(EXIT_FAILURE); } 
        if (_7orfqla[0] != VersionNumber)              {   perror("\n\nCode version number and input file version not equal (while reading *.rgh).   This might cause undefined behavior.\n");      }

        fseek(_yh1fkj4, (4*sizeof(float)+3*sizeof(unsigned int)+1*sizeof(unsigned short int)), SEEK_SET); 
        if (fread(_7orfqla, sizeof(unsigned int), _vxpj4hn->_3ru7myo, _yh1fkj4) != _vxpj4hn->_3ru7myo)       {   exit(EXIT_FAILURE); }
        for (i = 0u; i < _vxpj4hn->_3ru7myo; i++)       {   _vmefozt[i*8 +0] = _7orfqla[i] -1u;                        }
        if (fread(_7orfqla, sizeof(unsigned int), _vxpj4hn->_3ru7myo, _yh1fkj4) != _vxpj4hn->_3ru7myo)       {   exit(EXIT_FAILURE); }
        for (i = 0u; i < _vxpj4hn->_3ru7myo; i++)       {   _vmefozt[i*8 +1] = _7orfqla[i] -1u;                        }
        if (fread(_7orfqla, sizeof(unsigned int), _vxpj4hn->_3ru7myo, _yh1fkj4) != _vxpj4hn->_3ru7myo)       {   exit(EXIT_FAILURE); }
        for (i = 0u; i < _vxpj4hn->_3ru7myo; i++)       {   _vmefozt[i*8 +2] = _7orfqla[i] -1u;                        }
        if (fread(_7orfqla, sizeof(unsigned int), _vxpj4hn->_3ru7myo, _yh1fkj4) != _vxpj4hn->_3ru7myo)       {   exit(EXIT_FAILURE); }
        for (i = 0u; i < _vxpj4hn->_3ru7myo; i++)       {   _vmefozt[i*8 +3] = _7orfqla[i] -1u;                        }
        if (fread(_7orfqla, sizeof(unsigned int), _vxpj4hn->_3ru7myo, _yh1fkj4) != _vxpj4hn->_3ru7myo)       {   exit(EXIT_FAILURE); }
        for (i = 0u; i < _vxpj4hn->_3ru7myo; i++)       {   _vmefozt[i*8 +4] = _7orfqla[i] -1u;                        }
        if (fread(_7orfqla, sizeof(unsigned int), _vxpj4hn->_3ru7myo,_yh1fkj4) != _vxpj4hn->_3ru7myo)        {   exit(EXIT_FAILURE); }
        for (i = 0u; i < _vxpj4hn->_3ru7myo; i++)       {   _vmefozt[i*8 +5] = _7orfqla[i];                            }
        
        if (fread(_vm59m6m, sizeof(float), _vxpj4hn->_3ru7myo,_yh1fkj4) != _vxpj4hn->_3ru7myo)               {   exit(EXIT_FAILURE); }
        for (i = 0u; i < _vxpj4hn->_3ru7myo; i++)       {   _4xkss0k[i*16 +4] = _vm59m6m[i];                           }
        if (fread(_vm59m6m, sizeof(float), _vxpj4hn->_3ru7myo,_yh1fkj4) != _vxpj4hn->_3ru7myo)               {   exit(EXIT_FAILURE); }
        for (i = 0u; i < _vxpj4hn->_3ru7myo; i++)       {   _4xkss0k[i*16 +5] = _vm59m6m[i];                           }
        
        if (fread(_vm59m6m, sizeof(float), _vxpj4hn->_3ru7myo,_yh1fkj4) != _vxpj4hn->_3ru7myo)               {   exit(EXIT_FAILURE); }
        for (i = 0u; i < _vxpj4hn->_3ru7myo; i++)       {   _4xkss0k[i*16 +0] = _vm59m6m[i]*1000.0f;                   }
        if (fread(_vm59m6m, sizeof(float), _vxpj4hn->_3ru7myo,_yh1fkj4) != _vxpj4hn->_3ru7myo)               {   exit(EXIT_FAILURE); }
        for (i = 0u; i < _vxpj4hn->_3ru7myo; i++)       {   _4xkss0k[i*16 +1] = _vm59m6m[i]*1000.0f;                   }
        if (fread(_vm59m6m, sizeof(float), _vxpj4hn->_3ru7myo,_yh1fkj4) != _vxpj4hn->_3ru7myo)               {   exit(EXIT_FAILURE); }
        for (i = 0u; i < _vxpj4hn->_3ru7myo; i++)       {   _4xkss0k[i*16 +2] = _vm59m6m[i]*1000.0f;                   }
        
        if (fread(_ji9kqdg, sizeof(float), _s55j2o7->_fu7bfxq,_yh1fkj4) != _s55j2o7->_fu7bfxq)              {   exit(EXIT_FAILURE); }
        for (i = 0u; i < _s55j2o7->_fu7bfxq; i++)       {   _zzf8s2z[i*4 +0] = _ji9kqdg[i]*1000.0f;                  }
        if (fread(_ji9kqdg, sizeof(float), _s55j2o7->_fu7bfxq,_yh1fkj4) != _s55j2o7->_fu7bfxq)              {   exit(EXIT_FAILURE); }
        for (i = 0u; i < _s55j2o7->_fu7bfxq; i++)       {   _zzf8s2z[i*4 +1] = _ji9kqdg[i]*1000.0f;                  }
        if (fread(_ji9kqdg, sizeof(float), _s55j2o7->_fu7bfxq,_yh1fkj4) != _s55j2o7->_fu7bfxq)              {   exit(EXIT_FAILURE); }
        for (i = 0u; i < _s55j2o7->_fu7bfxq; i++)       {   _zzf8s2z[i*4 +2] = _ji9kqdg[i]*1000.0f;                  }
        if (fread(_ji9kqdg, sizeof(float), _s55j2o7->_fu7bfxq,_yh1fkj4) != _s55j2o7->_fu7bfxq)              {   exit(EXIT_FAILURE); }
        for (i = 0u; i < _s55j2o7->_fu7bfxq; i++)       {   _zzf8s2z[i*4 +3] = _ji9kqdg[i]*1000.0f;                  }
        
        fclose(_yh1fkj4);
    }
    
    if ((_yk1b3d5 = fopen(_1at7o75,"rb"))   == NULL)     {   perror("Error: -cant open *.stg file in function 'LoadRghStgthBnd'\n");      exit(EXIT_FAILURE);     }
    {   
        if (fread(_7orfqla, sizeof(unsigned int), 1, _yk1b3d5) != 1)       {   exit(EXIT_FAILURE); } 
        if (_7orfqla[0] != VersionNumber)              {   perror("\n\nCode version number and input file version not equal (while reading *.stg).   This might cause undefined behavior.\n");      }

        fseek(_yk1b3d5, 2*sizeof(unsigned int), SEEK_SET); 
        if (fread(_vm59m6m, sizeof(float), _vxpj4hn->_3ru7myo,_yk1b3d5) != _vxpj4hn->_3ru7myo) {   exit(EXIT_FAILURE);               }
        for (i = 0u; i < _lexdw6l[_fsqg8x8]; i++)       {      _0dl5e11[i*16 +8]  = _vm59m6m[_1nlc8jp[i]];            }
        if (fread(_vm59m6m, sizeof(float), _vxpj4hn->_3ru7myo,_yk1b3d5) != _vxpj4hn->_3ru7myo) {   exit(EXIT_FAILURE);               }
        for (i = 0u; i < _lexdw6l[_fsqg8x8]; i++)       {      _0dl5e11[i*16 +10] = _vm59m6m[_1nlc8jp[i]];            }
        if (fread(_vm59m6m, sizeof(float), _vxpj4hn->_3ru7myo,_yk1b3d5) != _vxpj4hn->_3ru7myo) {   exit(EXIT_FAILURE);               }
        for (i = 0u; i < _lexdw6l[_fsqg8x8]; i++)       {      _0dl5e11[i*16 +0]  = _vm59m6m[_1nlc8jp[i]];            }
        if (fread(_vm59m6m, sizeof(float), _vxpj4hn->_3ru7myo,_yk1b3d5) != _vxpj4hn->_3ru7myo) {   exit(EXIT_FAILURE);               }
        for (i = 0u; i < _lexdw6l[_fsqg8x8]; i++)       {      _0dl5e11[i*16 +12] = _vm59m6m[_1nlc8jp[i]];            }
        if (fread(_vm59m6m, sizeof(float), _vxpj4hn->_3ru7myo,_yk1b3d5) != _vxpj4hn->_3ru7myo) {   exit(EXIT_FAILURE);               }
        for (i = 0u; i < _lexdw6l[_fsqg8x8]; i++)       {      _0dl5e11[i*16 +9]  = 0.01f*_vm59m6m[_1nlc8jp[i]]*_0dl5e11[i*16 +8];  }
        if (fread(_vm59m6m, sizeof(float), _vxpj4hn->_3ru7myo,_yk1b3d5) != _vxpj4hn->_3ru7myo) {   exit(EXIT_FAILURE);                         }
        for (i = 0u; i < _lexdw6l[_fsqg8x8]; i++)       {      _0dl5e11[i*16 +11] = 0.01f*_vm59m6m[_1nlc8jp[i]]*_0dl5e11[i*16 +10]; }
        if (fread(_vm59m6m, sizeof(float), _vxpj4hn->_3ru7myo,_yk1b3d5) != _vxpj4hn->_3ru7myo) {   exit(EXIT_FAILURE);                        }
        for (i = 0u; i < _lexdw6l[_fsqg8x8]; i++)       {      _0dl5e11[i*16 +13] = 0.01f*_vm59m6m[_1nlc8jp[i]]*_0dl5e11[i*16 +12]; }
        if (fread(_vm59m6m, sizeof(float), _vxpj4hn->_3ru7myo,_yk1b3d5) != _vxpj4hn->_3ru7myo) {   exit(EXIT_FAILURE);                         }
        for (i = 0u; i < _lexdw6l[_fsqg8x8]; i++)       {      _0dl5e11[i*16 +8] += _vm59m6m[_1nlc8jp[i]];            }
        if (fread(_vm59m6m, sizeof(float), _vxpj4hn->_3ru7myo,_yk1b3d5) != _vxpj4hn->_3ru7myo) {   exit(EXIT_FAILURE);               }
        for (i = 0u; i < _lexdw6l[_fsqg8x8]; i++)       {      _0dl5e11[i*16 +10]+= _vm59m6m[_1nlc8jp[i]];            }
        if (fread(_vm59m6m, sizeof(float), _vxpj4hn->_3ru7myo,_yk1b3d5) != _vxpj4hn->_3ru7myo) {   exit(EXIT_FAILURE);               }
        for (i = 0u; i < _lexdw6l[_fsqg8x8]; i++)       {      _0dl5e11[i*16 +0] += _vm59m6m[_1nlc8jp[i]];            }
        if (fread(_vm59m6m, sizeof(float), _vxpj4hn->_3ru7myo,_yk1b3d5) != _vxpj4hn->_3ru7myo) {   exit(EXIT_FAILURE);               }
        for (i = 0u; i < _lexdw6l[_fsqg8x8]; i++)       {      _0dl5e11[i*16 +12]+= _vm59m6m[_1nlc8jp[i]];            } 
        
        fclose(_yk1b3d5);
    }
    
    if ((_050ui22 = fopen(_rh8lx6d,"rb"))   != NULL)
    {   
        if (fread(_7orfqla, sizeof(unsigned int), 1, _050ui22) != 1)       {   exit(EXIT_FAILURE); } 
        if (_7orfqla[0] != VersionNumber)              {   perror("\n\nCode version number and input file version not equal (while reading *.btrg).   This might cause undefined behavior.\n");      }

        fseek(_050ui22, (4*sizeof(unsigned int)), SEEK_SET);
        
        if (fread(_blmntfk,    sizeof(unsigned int), _vxpj4hn->_5izolnn,_050ui22) != _vxpj4hn->_5izolnn)  {   exit(EXIT_FAILURE); }
        for (i = 0u; i < _vxpj4hn->_5izolnn; i++)          {   _rgb93wk[i*8 +0] = _blmntfk[i] -1u;                  }
        if (fread(_blmntfk,    sizeof(unsigned int), _vxpj4hn->_5izolnn,_050ui22) != _vxpj4hn->_5izolnn)  {   exit(EXIT_FAILURE); }
        for (i = 0u; i < _vxpj4hn->_5izolnn; i++)          {   _rgb93wk[i*8 +1] = _blmntfk[i] -1u;                  }
        if (fread(_blmntfk,    sizeof(unsigned int), _vxpj4hn->_5izolnn,_050ui22) != _vxpj4hn->_5izolnn)  {   exit(EXIT_FAILURE); }
        for (i = 0u; i < _vxpj4hn->_5izolnn; i++)          {   _rgb93wk[i*8 +2] = _blmntfk[i] -1u;                  }
        if (fread(_blmntfk,    sizeof(unsigned int), _vxpj4hn->_5izolnn,_050ui22) != _vxpj4hn->_5izolnn)  {   exit(EXIT_FAILURE); }
        for (i = 0u; i < _vxpj4hn->_5izolnn; i++)          {   _rgb93wk[i*8 +3] = _blmntfk[i] -1u;                  }
        
        if (fread(_gywiyme,    sizeof(float), _vxpj4hn->_5izolnn,_050ui22) != _vxpj4hn->_5izolnn)  {   exit(EXIT_FAILURE);        }
        for (i = 0u; i < _vxpj4hn->_5izolnn; i++)          {   _rrpzlk3[i*16 +0] = _gywiyme[i]*1000.0f;             }
        if (fread(_gywiyme,    sizeof(float), _vxpj4hn->_5izolnn,_050ui22) != _vxpj4hn->_5izolnn)  {   exit(EXIT_FAILURE);        }
        for (i = 0u; i < _vxpj4hn->_5izolnn; i++)          {   _rrpzlk3[i*16 +1] = _gywiyme[i]*1000.0f;             }
        if (fread(_gywiyme,    sizeof(float), _vxpj4hn->_5izolnn,_050ui22) != _vxpj4hn->_5izolnn)  {   exit(EXIT_FAILURE);        }
        for (i = 0u; i < _vxpj4hn->_5izolnn; i++)          {   _rrpzlk3[i*16 +2] = _gywiyme[i]*1000.0f;             }
        
        if (fread(_nvg3kk0,    sizeof(float), _s55j2o7->_ryf0992,_050ui22) != _s55j2o7->_ryf0992)  {   exit(EXIT_FAILURE);       }
        for (i = 0u; i < _s55j2o7->_ryf0992; i++)          {   _3n3bzt9[i*4 +0] = _nvg3kk0[i]*1000.0f;            }
        if (fread(_nvg3kk0,    sizeof(float), _s55j2o7->_ryf0992,_050ui22) != _s55j2o7->_ryf0992)  {   exit(EXIT_FAILURE);       }
        for (i = 0u; i < _s55j2o7->_ryf0992; i++)          {   _3n3bzt9[i*4 +1] = _nvg3kk0[i]*1000.0f;            }
        if (fread(_nvg3kk0,    sizeof(float), _s55j2o7->_ryf0992,_050ui22) != _s55j2o7->_ryf0992)  {   exit(EXIT_FAILURE);       }
        for (i = 0u; i < _s55j2o7->_ryf0992; i++)          {   _3n3bzt9[i*4 +2] = _nvg3kk0[i]*1000.0f;            }
        
        if (fread(_gywiyme,    sizeof(float), _vxpj4hn->_5izolnn,_050ui22) != _vxpj4hn->_5izolnn)  {   exit(EXIT_FAILURE);        }
        for (i = 0u; i < _vxpj4hn->_5izolnn; i++)          {   _rrpzlk3[i*16 +3] = _gywiyme[i]*0.001f;              }
        if (fread(_gywiyme,    sizeof(float), _vxpj4hn->_5izolnn,_050ui22) != _vxpj4hn->_5izolnn)  {   exit(EXIT_FAILURE);        }
        for (i = 0u; i < _vxpj4hn->_5izolnn; i++)          {   _rrpzlk3[i*16 +4] = _gywiyme[i]*0.001f;              }
        if (fread(_gywiyme,    sizeof(float), _vxpj4hn->_5izolnn,_050ui22) != _vxpj4hn->_5izolnn)  {   exit(EXIT_FAILURE);        }
        for (i = 0u; i < _vxpj4hn->_5izolnn; i++)          {   _rrpzlk3[i*16 +5] = _gywiyme[i]*0.001f;              }
        
        fclose(_050ui22);
    }
    
    _ytkw283 = 0u;
    for (i = 0u; i < _vxpj4hn->_3ru7myo; i++)          {   _ytkw283 = MAX(_ytkw283, _vmefozt[i*8 +3]);           }
    _s55j2o7->_6wmdw0f = _ytkw283 +1u;
    
    _ytkw283 = 0u;
    for (i = 0u; i < _vxpj4hn->_5izolnn; i++)          {   _ytkw283 = MAX(_ytkw283, _rgb93wk[i*8 +3]);           }
    _s55j2o7->_rp9dxoq = (_vxpj4hn->_5izolnn > 0u)*(_ytkw283 +1u)  +  (_vxpj4hn->_5izolnn <= 0u)*0u;
    
    free(_7orfqla);   free(_vm59m6m);   free(_ji9kqdg);   free(_blmntfk);   free(_gywiyme);   free(_nvg3kk0);
    
    return;
}