#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <mpi.h>
#include <time.h>
#include <cblas.h>
#define VersionNumber 202511

#define ALIGN64 64 

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


extern long int  _z0yncxl(int _uuijnwp, int _qtkcd0r, int _tgxnx4s);

void   _rpj1epk(int _0rmhrcn, int _fsqg8x8, const struct _5r5uqpg *restrict _atj9tze, const struct _xjuy0bv *restrict _vxpj4hn, const int *restrict _lexdw6l, const int *restrict _1poxrte, const unsigned int *restrict _1nlc8jp, const unsigned int *restrict _fftouhy, float *restrict _6fokzix, float *restrict _z39pecn, float *restrict _7yt2m6t, const float *restrict _rrpzlk3, const unsigned short int *restrict _qnsn2s8, const unsigned short int *restrict _qpwbsf3, const unsigned short int *restrict _iymyrsv, const unsigned short int *restrict _iob0nlw, const unsigned short int **restrict _5gzpl1q, const unsigned short int **restrict _oqgvzhq, const unsigned short int **restrict _1ajo091, const unsigned short int **restrict _3b9so0u, const float **restrict _ckeh5sw, const float **restrict _3ac6pe6, const float **restrict _obyghgk, const float **restrict _yj9boq0, const float **restrict _ljefy3z, const float **restrict _h61iaoi, const float **restrict _plozjcn, const float **restrict _qk4ltah, const float **restrict _0ljch2v, const float **restrict _0dtirc9, const float **restrict _9dde3cu)
{   
    struct timespec _glmy62x, _1zocfi4;
    double _kicwmpe;
    
    unsigned int _18hafey = 100u;
    unsigned int i,   j,   k,   _5v2wwg8,   _nzrf5a7,   _yce4qgm,   _egok01t;
    float _d31k9qm,   _3p91i3v,   _v71448n = 0.0f,   _hez74la = 0.0f;
    float _e4xyf9i[3] = {0.0f, 0.0f, 0.0f};
    
    int __attribute__((aligned(ALIGN64))) _9lzs34j[_0rmhrcn];     memset(_9lzs34j, 0, _0rmhrcn*sizeof(int));
    int __attribute__((aligned(ALIGN64))) _sosx43w[_0rmhrcn];     memset(_sosx43w, 0, _0rmhrcn*sizeof(int));
    int __attribute__((aligned(ALIGN64))) _bh6nqqb[_0rmhrcn];     memset(_bh6nqqb, 0, _0rmhrcn*sizeof(int));
    int __attribute__((aligned(ALIGN64))) _1ni10kq[_0rmhrcn];     memset(_1ni10kq, 0, _0rmhrcn*sizeof(int));
    unsigned int  __attribute__((aligned(ALIGN64))) _xkc5i8x[16];
    float         __attribute__((aligned(ALIGN64))) _nbzug0i[16];
    long int _moxsfh4;

    float *_ozctwjl = NULL,   *_ny181ac = NULL,   *_t997e4u = NULL,   *_lqbmae3 = NULL,   *_gyspeng = NULL,   *_zkwdeek = NULL;
    
    _moxsfh4     = _z0yncxl(ALIGN64, 3*_lexdw6l[_fsqg8x8], sizeof(float));             _ozctwjl   = aligned_alloc(ALIGN64, _moxsfh4);                memset(_ozctwjl, 0, _moxsfh4);
    _moxsfh4     = _z0yncxl(ALIGN64, 3*_vxpj4hn->_3ru7myo, sizeof(float));             _ny181ac   = aligned_alloc(ALIGN64, _moxsfh4);                memset(_ny181ac, 0, _moxsfh4);
    _moxsfh4     = _z0yncxl(ALIGN64, 2*_vxpj4hn->_dn1th1x, sizeof(float));           _t997e4u    = aligned_alloc(ALIGN64, _moxsfh4);                memset(_t997e4u,  0, _moxsfh4);
    
    if (_vxpj4hn->_5izolnn > 0u)
    {   _moxsfh4 = _z0yncxl(ALIGN64, 4*_1poxrte[_fsqg8x8], sizeof(float));             _lqbmae3   = aligned_alloc(ALIGN64, _moxsfh4);                memset(_lqbmae3, 0, _moxsfh4);
        _moxsfh4 = _z0yncxl(ALIGN64, 4*_vxpj4hn->_5izolnn, sizeof(float));             _gyspeng   = aligned_alloc(ALIGN64, _moxsfh4);                memset(_gyspeng, 0, _moxsfh4);
        _moxsfh4 = _z0yncxl(ALIGN64, 3*_vxpj4hn->_w2maiet, sizeof(float));           _zkwdeek    = aligned_alloc(ALIGN64, _moxsfh4);                memset(_zkwdeek,  0, _moxsfh4);
    }
    const float  *_jbc2sws = (const float *)_ny181ac;
    const float  *_wyx9trf = (const float *)_gyspeng;
    const float  *_sh3tow4  = (const float *)_t997e4u;
    const float  *_8sryo4r  = (const float *)_zkwdeek;

    clock_gettime(CLOCK_REALTIME, &_glmy62x);
    
    
    unsigned int _8p6fvdh = 0u;
    
    if ( _8p6fvdh == 1u)
    {   _5v2wwg8 = 0u;
        for (i = 0u; i < _lexdw6l[_fsqg8x8]; i++)
        {   _ozctwjl[_5v2wwg8*3 +0] = (float)(_1nlc8jp[i]);
            _ozctwjl[_5v2wwg8*3 +1] = -1.0f*_6fokzix[i*16 +13]*_6fokzix[i*16 +8];
            _ozctwjl[_5v2wwg8*3 +2] = -1.0f*_6fokzix[i*16 +14]*_6fokzix[i*16 +8];
            _5v2wwg8              += 1u;
        }
        
        MPI_Allgather(&_5v2wwg8, 1, MPI_UNSIGNED, _sosx43w, 1, MPI_INT, MPI_COMM_WORLD);
        _9lzs34j[0] = 0u;
        _nzrf5a7    = _sosx43w[0];
        for (i = 1; i < _0rmhrcn; i++)         {        _nzrf5a7 += _sosx43w[i];        _sosx43w[i-1] *= 3;           _9lzs34j[i] = _9lzs34j[i-1] + _sosx43w[i-1];         }
        _sosx43w[_0rmhrcn-1] *= 3;
        MPI_Allgatherv(_ozctwjl, _sosx43w[_fsqg8x8], MPI_FLOAT, _ny181ac, _sosx43w, _9lzs34j, MPI_FLOAT, MPI_COMM_WORLD);
        
        for (i = 0u; i < _lexdw6l[_fsqg8x8]; i++)
        {   memset(_t997e4u, 0, 2*_qnsn2s8[i]*sizeof(float) );
            for (j = 0u; j < _nzrf5a7; j++) 
            {   _xkc5i8x[0] = (unsigned int)_jbc2sws[j*3 +0]; 
                _xkc5i8x[1] = _5gzpl1q[i][_xkc5i8x[0]] *2; 
                _xkc5i8x[0] = j*3;
                _t997e4u[_xkc5i8x[1] +0] += _jbc2sws[_xkc5i8x[0] +1];
                _t997e4u[_xkc5i8x[1] +1] += _jbc2sws[_xkc5i8x[0] +2];
            }
            _6fokzix[i*16 +9]  = cblas_sdot(2*_qnsn2s8[i], _sh3tow4, 1, _ckeh5sw[i], 1);
            _6fokzix[i*16 +10] = cblas_sdot(2*_qnsn2s8[i], _sh3tow4, 1, _3ac6pe6[i], 1);
            
            _z39pecn[i*8 +5] = sqrtf(_6fokzix[i*16 +13]*_6fokzix[i*16 +13] +_6fokzix[i*16 +14]*_6fokzix[i*16 +14]);
            _6fokzix[i*16 +11] = _6fokzix[i*16 +13];
            _6fokzix[i*16 +12] = _6fokzix[i*16 +14];
        }
    }
    else
    {   for (i = 0u; i < _lexdw6l[_fsqg8x8]; i++) 
        {   _nbzug0i[2]     = sqrtf(_6fokzix[i*16 +13]*_6fokzix[i*16 +13] +_6fokzix[i*16 +14]*_6fokzix[i*16 +14]);
            _e4xyf9i[0]  += (0.5*(_6fokzix[i*16 +5] + _6fokzix[i*16 +6]) );
            _e4xyf9i[1]  += (_nbzug0i[2] <= 0.0f)*0.0f  +  (_nbzug0i[2] > 0.0f)*_nbzug0i[2];
            _e4xyf9i[2]  += (_nbzug0i[2] <= 0.0f)*0.0f  +  (_nbzug0i[2 ]> 0.0f)*1.0f; 
            
            _6fokzix[i*16 +9] += (-1.0f*_6fokzix[i*16 +13] *_6fokzix[i*16 +5]); 
            _6fokzix[i*16 +10]+= (-1.0f*_6fokzix[i*16 +14] *_6fokzix[i*16 +6]); 
            
        }
        MPI_Allreduce(MPI_IN_PLACE, _e4xyf9i,  3, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD); 
        _v71448n = fabsf(1000.0f *_atj9tze->_14ve92m / (_e4xyf9i[0]/(float)_vxpj4hn->_3ru7myo));
        _nbzug0i[0] = (_e4xyf9i[2] <= 0.0f)*-1.0f  +  (_e4xyf9i[2] > 0.0f)*_e4xyf9i[2];
        _d31k9qm = (_e4xyf9i[2] <= 0.0f)*0.0f   +  (_e4xyf9i[2] > 0.0f)* (_e4xyf9i[1]/_nbzug0i[0]); 
        
        if (_vxpj4hn->_5izolnn > 0u)
        {   
            _yce4qgm = 0u;
            {   for (i = 0u; i < _1poxrte[_fsqg8x8]; i++)
                {   
                    _hez74la += ((_7yt2m6t[i*16 +5] + _7yt2m6t[i*16 +6] + _7yt2m6t[i*16 +7])/3.0f);         _nbzug0i[0]  = -1.0f*_7yt2m6t[i*16 +9];
                    _nbzug0i[1]  = -1.0f*_7yt2m6t[i*16 +10];                                           _nbzug0i[2]  = -1.0f*_7yt2m6t[i*16 +11];
                    
                    if ((fabsf(_nbzug0i[0]) + fabsf(_nbzug0i[1]) + fabsf(_nbzug0i[2])) > 0.0f)
                    {   _lqbmae3[_yce4qgm*4 +0] = (float)(_fftouhy[i]);
                        _lqbmae3[_yce4qgm*4 +1] = _nbzug0i[0]*_7yt2m6t[i*16 +8];
                        _lqbmae3[_yce4qgm*4 +2] = _nbzug0i[1]*_7yt2m6t[i*16 +8];
                        _lqbmae3[_yce4qgm*4 +3] = _nbzug0i[2]*_7yt2m6t[i*16 +8];
                        _yce4qgm              += 1u;
            }   }   }
            
            MPI_Allreduce(MPI_IN_PLACE, &_hez74la,  1, MPI_FLOAT,    MPI_SUM, MPI_COMM_WORLD); 
            _hez74la = fabsf(1000.0f *_atj9tze->_14ve92m / (_hez74la/(float)_vxpj4hn->_5izolnn)); 
            
            MPI_Allgather(&_yce4qgm, 1, MPI_UNSIGNED, _1ni10kq, 1, MPI_INT, MPI_COMM_WORLD);
            _bh6nqqb[0]  = 0u;
            _egok01t     = _1ni10kq[0];
            for (i = 1; i < _0rmhrcn; i++)         {       _egok01t  += _1ni10kq[i];        _1ni10kq[i-1] *= 4;           _bh6nqqb[i] = _bh6nqqb[i-1] + _1ni10kq[i-1];         }
            _1ni10kq[_0rmhrcn-1] *= 4;
            MPI_Allgatherv(_lqbmae3, _1ni10kq[_fsqg8x8], MPI_FLOAT, _gyspeng, _1ni10kq, _bh6nqqb, MPI_FLOAT, MPI_COMM_WORLD);
            
            for (i = 0u; i < _lexdw6l[_fsqg8x8]; i++)
            {   memset(_zkwdeek, 0, 3*_qpwbsf3[i]*sizeof(float) );
                for (j = 0u; j < _egok01t; j++) 
                {   _xkc5i8x[0] =  (unsigned int)_wyx9trf[j*4 +0]; 
                    _xkc5i8x[1] = _oqgvzhq[i][_xkc5i8x[0]] *3; 
                    _xkc5i8x[0] = j*4;
                    _zkwdeek[_xkc5i8x[1] +0] += _wyx9trf[_xkc5i8x[0] +1];
                    _zkwdeek[_xkc5i8x[1] +1] += _wyx9trf[_xkc5i8x[0] +2];
                    _zkwdeek[_xkc5i8x[1] +2] += _wyx9trf[_xkc5i8x[0] +3];
                }
                _6fokzix[i*16 +9] += cblas_sdot(3*_qpwbsf3[i], _8sryo4r, 1, _yj9boq0[i], 1);
                _6fokzix[i*16 +10]+= cblas_sdot(3*_qpwbsf3[i], _8sryo4r, 1, _ljefy3z[i], 1);
            }
        } 
        
        {   unsigned int _8cd31nr = 0u;
            for (i = 0u; i < _lexdw6l[_fsqg8x8]; i++)
            {   _6fokzix[i*16 +0]  = _6fokzix[i*16 +9];              _6fokzix[i*16 +1]  = _6fokzix[i*16 +10];
                _6fokzix[i*16 +11] = 0.0f;                       _6fokzix[i*16 +12] = 0.0f;
                 
                _nbzug0i[0] = sqrtf(_6fokzix[i*16 +9]*_6fokzix[i*16 +9] + _6fokzix[i*16 +10]*_6fokzix[i*16 +10]);
                _8cd31nr = (_nbzug0i[0] <= FLT_EPSILON)*_8cd31nr  +  (_nbzug0i[0] > FLT_EPSILON)*1u;
            }
            MPI_Allreduce(MPI_IN_PLACE, &_8cd31nr, 1, MPI_UNSIGNED, MPI_MAX, MPI_COMM_WORLD);
            if (_8cd31nr == 0u) {       fprintf(stdout,"no load was defined. abort \n");    exit(123);      }
        } 
        
        for (k = 0u; k < _18hafey; k++)
        {   
            _5v2wwg8     = 0u;
            _nzrf5a7 = 0u;
            for (i = 0u; i < _lexdw6l[_fsqg8x8]; i++)
            {   _nbzug0i[0]            = -_6fokzix[i*16 +0]/_6fokzix[i*16 +5];         _nbzug0i[1]        = -_6fokzix[i*16 +1]/_6fokzix[i*16 +6]; 
                _6fokzix[i*16 +11]     += _nbzug0i[0];                             _6fokzix[i*16 +12] += _nbzug0i[1];
                
                _ozctwjl[_5v2wwg8*3 +0] = (float)(_1nlc8jp[i]);
                _ozctwjl[_5v2wwg8*3 +1] = _nbzug0i[0]*_6fokzix[i*16 +8];
                _ozctwjl[_5v2wwg8*3 +2] = _nbzug0i[1]*_6fokzix[i*16 +8];
                _5v2wwg8              += 1u;
            }
            
            _yce4qgm     = 0u;
            _egok01t = 0u;
            for (i = 0u; i < _1poxrte[_fsqg8x8]; i++)
            {   _nbzug0i[0]            = -_7yt2m6t[i*16 +0]/_7yt2m6t[i*16 +5];         _nbzug0i[1]        = -_7yt2m6t[i*16 +1]/_7yt2m6t[i*16 +6];         _nbzug0i[2]        = -_7yt2m6t[i*16 +2]/_7yt2m6t[i*16 +7]; 
                
                _lqbmae3[_yce4qgm*4 +0] = (float)(_fftouhy[i]);
                _lqbmae3[_yce4qgm*4 +1] = _nbzug0i[0]*_7yt2m6t[i*16 +8];
                _lqbmae3[_yce4qgm*4 +2] = _nbzug0i[1]*_7yt2m6t[i*16 +8];
                _lqbmae3[_yce4qgm*4 +3] = _nbzug0i[2]*_7yt2m6t[i*16 +8];
                _yce4qgm              += 1u;
            }
            
            MPI_Allgather(&_5v2wwg8, 1, MPI_UNSIGNED, _sosx43w, 1, MPI_INT, MPI_COMM_WORLD);
            _9lzs34j[0] = 0u;
            _nzrf5a7    = _sosx43w[0];
            for (i = 1; i < _0rmhrcn; i++)         {        _nzrf5a7 += _sosx43w[i];        _sosx43w[i-1] *= 3;           _9lzs34j[i] = _9lzs34j[i-1] + _sosx43w[i-1];         }
            _sosx43w[_0rmhrcn-1] *= 3;
            MPI_Allgatherv(_ozctwjl, _sosx43w[_fsqg8x8], MPI_FLOAT, _ny181ac, _sosx43w, _9lzs34j, MPI_FLOAT, MPI_COMM_WORLD);
            
            if (_vxpj4hn->_5izolnn > 0u)
            {   MPI_Allgather(&_yce4qgm, 1, MPI_UNSIGNED, _1ni10kq, 1, MPI_INT, MPI_COMM_WORLD);
                _bh6nqqb[0] = 0u;
                _egok01t    = _1ni10kq[0];
                for (i = 1; i < _0rmhrcn; i++)         {        _egok01t += _1ni10kq[i];        _1ni10kq[i-1] *= 4;           _bh6nqqb[i] = _bh6nqqb[i-1] + _1ni10kq[i-1];         }
                _1ni10kq[_0rmhrcn-1] *= 4;
                MPI_Allgatherv(_lqbmae3, _1ni10kq[_fsqg8x8], MPI_FLOAT, _gyspeng, _1ni10kq, _bh6nqqb, MPI_FLOAT, MPI_COMM_WORLD);
            }
            
            for (i = 0u; i < _lexdw6l[_fsqg8x8]; i++)
            {   
                if (_vxpj4hn->_5izolnn > 0u)
                {   memset(_zkwdeek, 0, 3*_qpwbsf3[i]*sizeof(float) );
                    for (j = 0u; j < _egok01t; j++) 
                    {   _xkc5i8x[2]   = (unsigned int)_wyx9trf[j*4 +0]; 
                        _xkc5i8x[1]   = _oqgvzhq[i][_xkc5i8x[2]] *3; 
                        _xkc5i8x[2]   = j*4;
                        _zkwdeek[_xkc5i8x[1] +0] += _wyx9trf[_xkc5i8x[2] +1];
                        _zkwdeek[_xkc5i8x[1] +1] += _wyx9trf[_xkc5i8x[2] +2];
                        _zkwdeek[_xkc5i8x[1] +2] += _wyx9trf[_xkc5i8x[2] +3];
                    }
                    _6fokzix[i*16 +0] += cblas_sdot(3*_qpwbsf3[i], _8sryo4r, 1, _yj9boq0[i], 1);
                    _6fokzix[i*16 +1] += cblas_sdot(3*_qpwbsf3[i], _8sryo4r, 1, _ljefy3z[i], 1);
                }
                
                memset(_t997e4u, 0, 2*_qnsn2s8[i]*sizeof(float) );
                for (j = 0u; j < _nzrf5a7; j++) 
                {   _xkc5i8x[2] = (unsigned int)_jbc2sws[j*3 +0]; 
                    _xkc5i8x[1] = _5gzpl1q[i][_xkc5i8x[2]] *2; 
                    _xkc5i8x[2] = j*3;
                    _t997e4u[_xkc5i8x[1] +0] += _jbc2sws[_xkc5i8x[2] +1];
                    _t997e4u[_xkc5i8x[1] +1] += _jbc2sws[_xkc5i8x[2] +2];
                }
                _6fokzix[i*16 +0] += cblas_sdot(2*_qnsn2s8[i], _sh3tow4, 1, _ckeh5sw[i], 1);
                _6fokzix[i*16 +1] += cblas_sdot(2*_qnsn2s8[i], _sh3tow4, 1, _3ac6pe6[i], 1);
            } 
            
            for (i = 0u; i < _1poxrte[_fsqg8x8]; i++)
            {   
                memset(_zkwdeek, 0, 3*_iymyrsv[i]*sizeof(float) );
                for (j = 0u; j < _egok01t; j++) 
                {   _xkc5i8x[2]  = (unsigned int)_wyx9trf[j*4 +0]; 
                    _xkc5i8x[1]  = _1ajo091[i][_xkc5i8x[2]] *3; 
                    _xkc5i8x[2]  = j*4;
                    _zkwdeek[_xkc5i8x[1] +0] += _wyx9trf[_xkc5i8x[2] +1];
                    _zkwdeek[_xkc5i8x[1] +1] += _wyx9trf[_xkc5i8x[2] +2];
                    _zkwdeek[_xkc5i8x[1] +2] += _wyx9trf[_xkc5i8x[2] +3];
                }
                _7yt2m6t[i*16 +0] += cblas_sdot(3*_iymyrsv[i], _8sryo4r, 1, _h61iaoi[i], 1);
                _7yt2m6t[i*16 +1] += cblas_sdot(3*_iymyrsv[i], _8sryo4r, 1, _plozjcn[i], 1);
                _7yt2m6t[i*16 +2] += cblas_sdot(3*_iymyrsv[i], _8sryo4r, 1, _qk4ltah[i], 1);
                
                memset(_t997e4u, 0, 2*_iob0nlw[i]*sizeof(float) );
                for (j = 0u; j < _nzrf5a7; j++) 
                {   _xkc5i8x[2]  = (unsigned int)_jbc2sws[j*3 +0]; 
                    _xkc5i8x[1]  = _3b9so0u[i][_xkc5i8x[2]] *2; 
                    _xkc5i8x[2]  = j*3;
                    _t997e4u[_xkc5i8x[1] +0] += _jbc2sws[_xkc5i8x[2] +1];
                    _t997e4u[_xkc5i8x[1] +1] += _jbc2sws[_xkc5i8x[2] +2];
                }
                _7yt2m6t[i*16 +0] += cblas_sdot(2*_iob0nlw[i], _sh3tow4, 1, _0ljch2v[i], 1);
                _7yt2m6t[i*16 +1] += cblas_sdot(2*_iob0nlw[i], _sh3tow4, 1, _0dtirc9[i], 1);
                _7yt2m6t[i*16 +2] += cblas_sdot(2*_iob0nlw[i], _sh3tow4, 1, _9dde3cu[i], 1);
            }
            
        }
        
        {   memset(_e4xyf9i, 0, 3*sizeof(float));
            
            for (i = 0u; i < _lexdw6l[_fsqg8x8]; i++) 
            {   _nbzug0i[2]  = fabsf(_6fokzix[i*16 +13]) + fabsf(_6fokzix[i*16 +14]); 
                _nbzug0i[3]  = sqrtf(_6fokzix[i*16 +11]*_6fokzix[i*16 +11] +_6fokzix[i*16 +12]*_6fokzix[i*16 +12]);
                _e4xyf9i[1] += (_nbzug0i[2] <= 0.0f)*0.0f  +  (_nbzug0i[2] > 0.0f)*_nbzug0i[3];
                _e4xyf9i[2] += (_nbzug0i[2] <= 0.0f)*0.0f  +  (_nbzug0i[2] > 0.0f)*1.0f;
                
                _z39pecn[i*8 +5] = _nbzug0i[3]; 
            }
            
            MPI_Allreduce(MPI_IN_PLACE, &_e4xyf9i, 3, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
            _nbzug0i[0] = (_e4xyf9i[2] <= 0.0f)*-1.0f  +  (_e4xyf9i[2] > 0.0f)*_e4xyf9i[2];
            _3p91i3v = (_e4xyf9i[2] <= 0.0f)*0.0f   +  (_e4xyf9i[2] > 0.0f)*(_e4xyf9i[1]/_nbzug0i[0]);
            
            if ((_d31k9qm > FLT_EPSILON) && (_3p91i3v > FLT_EPSILON))   {   _nbzug0i[0] = _d31k9qm/_3p91i3v;     }
            else                                                        {   _nbzug0i[0] = 1.0f;                  }
            
            for (i = 0u; i < _lexdw6l[_fsqg8x8]; i++) 
            {   _6fokzix[i*16 +11] *= _nbzug0i[0];       _6fokzix[i*16 +12] *= _nbzug0i[0];           _z39pecn[i*8 +5] *= _nbzug0i[0]; 
            }
        }
        
        {   _5v2wwg8 = 0u;
            
            for (i = 0u; i < _lexdw6l[_fsqg8x8]; i++)
            {   _ozctwjl[_5v2wwg8*3 +0] = (float)(_1nlc8jp[i]);
                _ozctwjl[_5v2wwg8*3 +1] = -_6fokzix[i*16 +11]*_6fokzix[i*16 +8];
                _ozctwjl[_5v2wwg8*3 +2] = -_6fokzix[i*16 +12]*_6fokzix[i*16 +8];
                _5v2wwg8              += 1u;
            }
            
            MPI_Allgather(&_5v2wwg8, 1, MPI_UNSIGNED, _sosx43w, 1, MPI_INT, MPI_COMM_WORLD);
            _9lzs34j[0] = 0u;
            _nzrf5a7    = _sosx43w[0];
            for (i = 1; i < _0rmhrcn; i++)         {       _nzrf5a7  += _sosx43w[i];        _sosx43w[i-1] *= 3;           _9lzs34j[i] = _9lzs34j[i-1] + _sosx43w[i-1];         }
            _sosx43w[_0rmhrcn-1] *= 3;
            MPI_Allgatherv(_ozctwjl, _sosx43w[_fsqg8x8], MPI_FLOAT, _ny181ac, _sosx43w, _9lzs34j, MPI_FLOAT, MPI_COMM_WORLD);
            
            for (i = 0u; i < _lexdw6l[_fsqg8x8]; i++)
            {   memset(_t997e4u, 0, 2*_qnsn2s8[i]*sizeof(float) );
                for (j = 0u; j < _nzrf5a7; j++) 
                {   _xkc5i8x[0] = (unsigned int)_jbc2sws[j*3 +0]; 
                    _xkc5i8x[1] = _5gzpl1q[i][_xkc5i8x[0]] *2; 
                    _xkc5i8x[0] = j*3;
                    _t997e4u[_xkc5i8x[1] +0] += _jbc2sws[_xkc5i8x[0] +1];
                    _t997e4u[_xkc5i8x[1] +1] += _jbc2sws[_xkc5i8x[0] +2];
                }
                _6fokzix[i*16 +9]  = cblas_sdot(2*_qnsn2s8[i], _sh3tow4, 1, _ckeh5sw[i], 1);
                _6fokzix[i*16 +10] = cblas_sdot(2*_qnsn2s8[i], _sh3tow4, 1, _3ac6pe6[i], 1);
            }
        } 
        
        {   float _g34iq7z,   _uug95rl,   _nbzug0i;
            
            _nbzug0i = 0.0f;
            for (i = 0u; i < _lexdw6l[_fsqg8x8]; i++)          {           _nbzug0i += _z39pecn[i*8 +5];                 }
            MPI_Allreduce(MPI_IN_PLACE, &_nbzug0i,  1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
            _g34iq7z = _nbzug0i/(float)_vxpj4hn->_3ru7myo;
            
            for (i = 0u; i < _lexdw6l[_fsqg8x8]; i++)
            {   _6fokzix[i*16  +9] /= MAX(1.0f, _z39pecn[i*8 +5]/_g34iq7z); 
                _6fokzix[i*16 +10] /= MAX(1.0f, _z39pecn[i*8 +5]/_g34iq7z);
                _6fokzix[i*16 +11] /= MAX(1.0f, _z39pecn[i*8 +5]/_g34iq7z);
                _6fokzix[i*16 +12] /= MAX(1.0f, _z39pecn[i*8 +5]/_g34iq7z);
                _z39pecn[i*8 +5] /= MAX(1.0f, _z39pecn[i*8 +5]/_g34iq7z);
            }
            
            _nbzug0i = 0.0f;
            for (i = 0u; i < _lexdw6l[_fsqg8x8]; i++)          {           _nbzug0i += _z39pecn[i*8 +5];                 }
            MPI_Allreduce(MPI_IN_PLACE, &_nbzug0i,  1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
            _uug95rl = _nbzug0i/(float)_vxpj4hn->_3ru7myo;
            
            for (i = 0u; i < _lexdw6l[_fsqg8x8]; i++)
            {   _6fokzix[i*16  +9] *= (_g34iq7z/_uug95rl);
                _6fokzix[i*16 +10] *= (_g34iq7z/_uug95rl);
                _6fokzix[i*16 +11] *= (_g34iq7z/_uug95rl); 
                _6fokzix[i*16 +12] *= (_g34iq7z/_uug95rl);
                _z39pecn[i*8 +5] *= (_g34iq7z/_uug95rl);
            }
            
        } 
        




    }
    MPI_Barrier( MPI_COMM_WORLD );
    clock_gettime(CLOCK_REALTIME, &_1zocfi4);
    _kicwmpe = (_1zocfi4.tv_sec - _glmy62x.tv_sec) + (_1zocfi4.tv_nsec - _glmy62x.tv_nsec)/1.0E+9;
    if (_fsqg8x8 == 0)
    {   fprintf(stdout,"Time taken to determine loading function (in seconds): %lf\n",_kicwmpe);
        fprintf(stdout,"Minimum AfterSlip (mm): %3.3f\n", _v71448n);
        fprintf(stdout,"\nResulting average tectonic slip-rate on faults (mm/yr): %5.5f   and  %5.5f\n",_d31k9qm*1000.0f, _3p91i3v*1000.0f);
    }
    
    free (_ozctwjl);   free(_ny181ac);   free (_lqbmae3);   free(_gyspeng);   free(_t997e4u);   free(_zkwdeek);
    
    return;
}