#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <mpi.h>
#include <time.h>
#define VersionNumber 202511

#define ALIGN64 64 
#define fMinDistFact 1.75f 

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


extern void    _id4n1u3(float *restrict _j49yd82, const float *restrict _28xg8y7, const float *restrict _kf1622v);
extern float    _mabsi9n(const float *restrict _n42x2ez);
extern void       _bfuk9e1(float *restrict _n42x2ez);
extern void _ktakpux(const float *restrict _088o98b, float *restrict _j3swebn, float *restrict _795e5iz);
extern void     _ftu2gc5(float *restrict _0i9g5dz, const float *restrict _mt8vt43, const float *restrict _xyzx9jp);
extern void     _yyv37di(float *restrict _0yknfu8, float _ktmonrv);
extern float            _53ayeqa(long *_rrhf4ab);
extern long int  _z0yncxl(int _uuijnwp, int _qtkcd0r, int _tgxnx4s);

extern void _myx444n(struct _5r5uqpg *restrict _atj9tze, const struct _xjuy0bv *restrict _vxpj4hn, unsigned int i, unsigned int *restrict _p422ta8, float *restrict _6fokzix, float *restrict _0dl5e11, unsigned int _urvqt3g);

extern void _oyr6ccv(float _cjwtl9t[6], float _5lf6a4f[6], float X, float Y, float Z, float _vkizzqy[3], float _x63cfl0[3], float _rcbpljd[3], float _xlpu0dk, float _zfk8ut0, float _fcd6mzs, const float _mu7jae5, const float _odck828);
extern void _25sq8cr(float _cjwtl9t[6], float _5lf6a4f[6], float X, float Y, float Z, float _vkizzqy[3], float _x63cfl0[3], float _rcbpljd[3], float _xlpu0dk, float _zfk8ut0, float _fcd6mzs, const float _mu7jae5, const float _odck828);

void   _nb4gh2m(float _gnmby49, unsigned int _7c6qc6k, const float *restrict _bjfvyb6, unsigned int _6zjq7w7,  const unsigned int *restrict _r59c81d, const float *restrict _yn4u9bx, unsigned int ***restrict _h9ier7u, unsigned int ***restrict _xdrkz40, unsigned int *restrict _katqyds);
void _h75lxo6(unsigned int _6zjq7w7,  const unsigned int *restrict _r59c81d, const float *restrict _yn4u9bx, const float *restrict _mszdxim, unsigned int **restrict _nwngucz, unsigned int _katqyds,  float ***restrict _vwvm4w3);
void  _e4q6qdy(const char *restrict _dygflij, unsigned int _jz6us7r, const float *restrict _jdtvv2y, unsigned int **restrict _1xgbngz, unsigned int _gzvrctl, float **restrict _4pa8hye, unsigned int _ntbz3m5, const float *restrict _8e7flqh, unsigned int **restrict _f9yecaz, unsigned int _2vpgstq, float **restrict _32orwix);



void  _43d5iqv(int _0rmhrcn, int _fsqg8x8, const char *restrict _dygflij, struct _5r5uqpg *_atj9tze, struct _xjuy0bv *restrict _vxpj4hn, const struct _j49itn7 *restrict _s55j2o7, const int *restrict _lexdw6l, const int *restrict _1poxrte, const unsigned int *restrict _1nlc8jp, const unsigned int *restrict _fftouhy, const float *restrict _e9ftaxp, const float *restrict _9yo18pa, const unsigned int *restrict _vmefozt, const unsigned int *restrict _rgb93wk, const float *restrict _4xkss0k, const float *restrict _rrpzlk3, const float *restrict _zzf8s2z, const float *restrict _3n3bzt9, unsigned int *restrict _p422ta8, float *restrict _6fokzix, float *restrict _0dl5e11, float *restrict _7yt2m6t, unsigned short int *restrict _qnsn2s8, unsigned short int *restrict _qpwbsf3, unsigned short int *restrict _iymyrsv, unsigned short int *restrict _iob0nlw, unsigned short int ***restrict _5gzpl1q, unsigned short int ***restrict _oqgvzhq, unsigned short int ***restrict _1ajo091, unsigned short int ***restrict _3b9so0u, float ***restrict _ckeh5sw, float ***restrict _3ac6pe6, float ***restrict _obyghgk, float ***restrict _yj9boq0, float ***restrict _ljefy3z, float ***restrict _h61iaoi, float ***restrict _plozjcn, float ***restrict _qk4ltah, float ***restrict _0ljch2v, float ***restrict _0dtirc9, float ***restrict _9dde3cu)
{   
    struct timespec _glmy62x, _1zocfi4;
    double _kicwmpe;
    
    unsigned int _l9ue421 = 0u; 
    unsigned int _kb700uc = 0u; 
    unsigned int **_2h2p39y = NULL,   **_qsq0p82 = NULL,   **_2qgv3dr = NULL,   **_cw2lbm2 = NULL;
    float       **_e3nwxuj = NULL,   **_nf74t60 = NULL;
    
    clock_gettime(CLOCK_REALTIME, &_glmy62x);
    
    _nb4gh2m(    _s55j2o7->_b68rerv, _s55j2o7->_6wmdw0f, _e9ftaxp, _vxpj4hn->_3ru7myo, _vmefozt, _4xkss0k, &_2h2p39y, &_2qgv3dr, &_l9ue421);
    if (_vxpj4hn->_5izolnn > 0u)
    {   _nb4gh2m(_s55j2o7->_gfu8iog, _s55j2o7->_rp9dxoq, _9yo18pa, _vxpj4hn->_5izolnn, _rgb93wk, _rrpzlk3, &_qsq0p82, &_cw2lbm2, &_kb700uc);
    }
    
    MPI_Barrier( MPI_COMM_WORLD );
    clock_gettime(CLOCK_REALTIME, &_1zocfi4);
    _kicwmpe = (_1zocfi4.tv_sec - _glmy62x.tv_sec) + (_1zocfi4.tv_nsec - _glmy62x.tv_nsec)/1.0E+9;

    if (_fsqg8x8 == 0)         {   fprintf(stdout,"BranchCountF: %7u\nTime taken (in seconds): %lf\n", _l9ue421, _kicwmpe);             }
    
    _h75lxo6(    _vxpj4hn->_3ru7myo,  _vmefozt, _4xkss0k, _zzf8s2z, _2h2p39y, _l9ue421, &_e3nwxuj);
    if (_vxpj4hn->_5izolnn > 0u)
    {   _h75lxo6(_vxpj4hn->_5izolnn,  _rgb93wk, _rrpzlk3, _3n3bzt9, _qsq0p82, _kb700uc, &_nf74t60);
    }
    

    
    MPI_Barrier( MPI_COMM_WORLD );
    
    clock_gettime(CLOCK_REALTIME, &_glmy62x);
    
    
    {   unsigned int i,   j;
        unsigned int _d0p92k3,   _qw4c6lb,   _aq15o2n,   _4nnkm95,   _8h0imjy,   _8djzr5s,   _omxlvz4,   _qke2gjr,   _81bmqf7;
        
        float _t26ys44,   _5gdej31,   _gymdm0l,   _qw5dsk7,   _31yrujo;
        float _pt5r5rc[3],   _9p3jpxg[3],   _dan085r[3],   _ntdr5xy[3],   _fc4kl8c[3],   _gn0t6dv[3],   _n42x2ez[3];
        float _x4nxvkd[6],   _4xyjoa2[6],   _2k8xpnm[6],   _c9x5jz6[6],   _1a6xvaj[6],   _jqy3ods[6],   _vg5vmf6[6],   _p9p9w4j[9],   _kibu4fm[9],   _zy1uu9y[9];
        long int _moxsfh4;
        
        unsigned int *_be3yqy5 = NULL,   *_t4q9bc2 = NULL,   *_otn9uze = NULL,   *_mgvu4is = NULL;
        float *_2292mtr = NULL,   *_3lbqgu3 = NULL,   *_5bwikot = NULL,   *_iq2fa39 = NULL,   *_b4ldims = NULL;
        float *_m458ax4 = NULL,   *_s6g298c = NULL,   *_ukwqvr1 = NULL,   *_ffqot8h = NULL,   *_zd76szn = NULL,   *_p8ijl2z = NULL;
        
        _be3yqy5         = malloc(_vxpj4hn->_3ru7myo *sizeof(unsigned int)); 
        _t4q9bc2       = malloc(_l9ue421    *sizeof(unsigned int)); 
        _moxsfh4       = _z0yncxl(ALIGN64, 2*_vxpj4hn->_3ru7myo, sizeof(float));
        _2292mtr     = aligned_alloc(ALIGN64, _moxsfh4);
        _3lbqgu3     = aligned_alloc(ALIGN64, _moxsfh4);
        _5bwikot     = aligned_alloc(ALIGN64, _moxsfh4);
        
        if (_vxpj4hn->_5izolnn > 0u)
        {   _otn9uze     = malloc(_vxpj4hn->_5izolnn *sizeof(unsigned int)); 
            _mgvu4is   = malloc(_kb700uc    *sizeof(unsigned int));
            _moxsfh4   = _z0yncxl(ALIGN64, 3*_vxpj4hn->_5izolnn, sizeof(float));
            _iq2fa39 = aligned_alloc(ALIGN64, _moxsfh4);
            _b4ldims = aligned_alloc(ALIGN64, _moxsfh4);
        }
        
        _moxsfh4              = _z0yncxl(ALIGN64, _lexdw6l[_fsqg8x8], sizeof(unsigned short int *));
        *_5gzpl1q              = aligned_alloc(ALIGN64, _moxsfh4);
        _moxsfh4              = _z0yncxl(ALIGN64, _vxpj4hn->_3ru7myo, sizeof(unsigned short int));
        for (i = 0u; i < _lexdw6l[_fsqg8x8]; i++)
        {   (*_5gzpl1q)[i]     = aligned_alloc(ALIGN64, _moxsfh4);         memset((*_5gzpl1q)[i], 0, _moxsfh4);
        }
        
        if (_vxpj4hn->_5izolnn > 0u)
        {   _moxsfh4          = _z0yncxl(ALIGN64, _lexdw6l[_fsqg8x8], sizeof(unsigned short int *));
            *_oqgvzhq          = aligned_alloc(ALIGN64, _moxsfh4);
            _moxsfh4          = _z0yncxl(ALIGN64, _vxpj4hn->_5izolnn, sizeof(unsigned short int));
            for (i = 0u; i < _lexdw6l[_fsqg8x8]; i++)
            {   (*_oqgvzhq)[i] = aligned_alloc(ALIGN64, _moxsfh4);         memset((*_oqgvzhq)[i], 0, _moxsfh4);
            }
        }
        
        _moxsfh4       = _z0yncxl(ALIGN64, _lexdw6l[_fsqg8x8], sizeof(float *));
        *_ckeh5sw      = aligned_alloc(ALIGN64, _moxsfh4);
        *_3ac6pe6      = aligned_alloc(ALIGN64, _moxsfh4);
        *_obyghgk      = aligned_alloc(ALIGN64, _moxsfh4);
        if (_vxpj4hn->_5izolnn > 0u)
        {   *_yj9boq0  = aligned_alloc(ALIGN64, _moxsfh4);
            *_ljefy3z  = aligned_alloc(ALIGN64, _moxsfh4);
        }
        
        for (i = 0; i < _lexdw6l[_fsqg8x8]; i++) 
        {   _8h0imjy = _1nlc8jp[i];
            
            _moxsfh4   = _z0yncxl(ALIGN64, 2*_vxpj4hn->_3ru7myo, sizeof(float));
            memset(_2292mtr, 0, _moxsfh4);
            memset(_3lbqgu3, 0, _moxsfh4);
            memset(_5bwikot, 0, _moxsfh4);
            
            memcpy(_pt5r5rc,    &_4xkss0k[_8h0imjy*16 +0],  3*sizeof(float));
            memcpy(_p9p9w4j, &_4xkss0k[_8h0imjy*16 +6],  9*sizeof(float)); 
            
            memset(_be3yqy5,   0, _vxpj4hn->_3ru7myo*sizeof(unsigned int));
            memset(_t4q9bc2, 0, _l9ue421*sizeof(unsigned int));
            _d0p92k3  = 0u;
            
            while (_d0p92k3 < _vxpj4hn->_3ru7myo) 
            {   
                for (j = 0u; j < _vxpj4hn->_3ru7myo; j++)               {       if (_be3yqy5[j] == 0u)    {      _8djzr5s = j;       break;          }           } 
                
                memcpy(_9p3jpxg, &_4xkss0k[_8djzr5s*16 +0], 3*sizeof(float));
                
                _id4n1u3(_n42x2ez, _pt5r5rc, _9p3jpxg);               _31yrujo  = _mabsi9n(_n42x2ez);
                
                _qw5dsk7 = (_31yrujo     >= fMinDistFact*_s55j2o7->_b68rerv)*1.0f  +  (_31yrujo      < fMinDistFact*_s55j2o7->_b68rerv)*0.0f; 
                _qw5dsk7 = (_vmefozt[_8h0imjy*8 +4] == _vmefozt[_8djzr5s*8 +4])*1.0f  +  (_vmefozt[_8h0imjy*8 +4] != _vmefozt[_8djzr5s*8 +4])*_qw5dsk7;
                
                _omxlvz4 = _2h2p39y[_8djzr5s][0]; 
                
                _qw4c6lb = 0;
                _qke2gjr   = 1u;
                while (_qke2gjr == 1u)
                {   
                    _81bmqf7 = 0u;
                    while (_qw4c6lb < _omxlvz4-1)
                    {   _aq15o2n = _2h2p39y[_8djzr5s][16 -_omxlvz4 +_qw4c6lb];
                        
                        memcpy(_9p3jpxg, &_e3nwxuj[_aq15o2n][0], 3*sizeof(float));
                        _id4n1u3(_n42x2ez, _pt5r5rc, _9p3jpxg);               _31yrujo  = _mabsi9n(_n42x2ez);
                        
                        _81bmqf7  = (_t4q9bc2[_aq15o2n]                                                == 0u)*1u      +  (_t4q9bc2[_aq15o2n]                                                  != 0u)*0u;
                        _81bmqf7  = (_2qgv3dr[_aq15o2n][0]                                             > 1u)*_81bmqf7  +  (_2qgv3dr[_aq15o2n][0]                                              <= 1u)*0u;
                        _81bmqf7  = (_2h2p39y[_8h0imjy][16 -_omxlvz4 +_qw4c6lb]                    != _aq15o2n)*_81bmqf7  +  (_2h2p39y[_8h0imjy][16 -_omxlvz4 +_qw4c6lb]                      == _aq15o2n)*0u;
                        _81bmqf7  = (_31yrujo > (_e3nwxuj[_aq15o2n][5] +_s55j2o7->_a11dpg8*_s55j2o7->_b68rerv))*_81bmqf7  +  (_31yrujo  <= (_e3nwxuj[_aq15o2n][5] +_s55j2o7->_a11dpg8*_s55j2o7->_b68rerv))*0u;
                        
                        if (_81bmqf7 == 1u)         {       break;      }
                        
                        _qw4c6lb += 1u;
                    } 
                    
                    if (_81bmqf7 == 0u) 
                    {   memcpy(_dan085r, &_zzf8s2z[_vmefozt[_8djzr5s*8 +0]*4 +0], 3*sizeof(float));
                        memcpy(_ntdr5xy, &_zzf8s2z[_vmefozt[_8djzr5s*8 +1]*4 +0], 3*sizeof(float));
                        memcpy(_fc4kl8c, &_zzf8s2z[_vmefozt[_8djzr5s*8 +2]*4 +0], 3*sizeof(float));
                        
                        if (_s55j2o7->_s832ejx == 1u)
                        {   _oyr6ccv(_x4nxvkd, _vg5vmf6, _pt5r5rc[0], _pt5r5rc[1], _pt5r5rc[2], _dan085r, _ntdr5xy, _fc4kl8c, _s55j2o7->_ycz04cd, 0.0f, 0.0f, _atj9tze->_kkrzd5e, _atj9tze->_6snglh6);
                            _oyr6ccv(_4xyjoa2, _vg5vmf6, _pt5r5rc[0], _pt5r5rc[1], _pt5r5rc[2], _dan085r, _ntdr5xy, _fc4kl8c, 0.0f, _s55j2o7->_ycz04cd, 0.0f, _atj9tze->_kkrzd5e, _atj9tze->_6snglh6);
                            _oyr6ccv(_2k8xpnm, _vg5vmf6, _pt5r5rc[0], _pt5r5rc[1], _pt5r5rc[2], _dan085r, _ntdr5xy, _fc4kl8c, 0.0f, 0.0f, _s55j2o7->_ycz04cd, _atj9tze->_kkrzd5e, _atj9tze->_6snglh6);
                        }
                        else
                        {   _25sq8cr(_x4nxvkd, _vg5vmf6, _pt5r5rc[0], _pt5r5rc[1], _pt5r5rc[2], _dan085r, _ntdr5xy, _fc4kl8c, _s55j2o7->_ycz04cd, 0.0f, 0.0f, _atj9tze->_kkrzd5e, _atj9tze->_6snglh6);
                            _25sq8cr(_4xyjoa2, _vg5vmf6, _pt5r5rc[0], _pt5r5rc[1], _pt5r5rc[2], _dan085r, _ntdr5xy, _fc4kl8c, 0.0f, _s55j2o7->_ycz04cd, 0.0f, _atj9tze->_kkrzd5e, _atj9tze->_6snglh6);
                            _25sq8cr(_2k8xpnm, _vg5vmf6, _pt5r5rc[0], _pt5r5rc[1], _pt5r5rc[2], _dan085r, _ntdr5xy, _fc4kl8c, 0.0f, 0.0f, _s55j2o7->_ycz04cd, _atj9tze->_kkrzd5e, _atj9tze->_6snglh6);
                        }
                        
                        _yyv37di(_x4nxvkd, (_qw5dsk7/(_s55j2o7->_ycz04cd*_4xkss0k[_8djzr5s*16 +15])));      _ftu2gc5(&_kibu4fm[0], _p9p9w4j, _x4nxvkd);
                        _yyv37di(_4xyjoa2, (_qw5dsk7/(_s55j2o7->_ycz04cd*_4xkss0k[_8djzr5s*16 +15])));      _ftu2gc5(&_kibu4fm[3], _p9p9w4j, _4xyjoa2);
                        _yyv37di(_2k8xpnm, (_qw5dsk7/(_s55j2o7->_ycz04cd*_4xkss0k[_8djzr5s*16 +15])));      _ftu2gc5(&_kibu4fm[6], _p9p9w4j, _2k8xpnm);
                        
                        if (_8h0imjy == _8djzr5s)
                        {   _6fokzix[i*16 +5] = _kibu4fm[0]*_4xkss0k[_8djzr5s*16 +15];
                            _6fokzix[i*16 +6] = _kibu4fm[4]*_4xkss0k[_8djzr5s*16 +15];
                            _6fokzix[i*16 +7] = _kibu4fm[8]*_4xkss0k[_8djzr5s*16 +15];
                            _0dl5e11[i*16 +1] = MAX(fabsf(_6fokzix[i*16 +5]),fabsf(_6fokzix[i*16 +6])); 
                            _0dl5e11[i*16 +1] = MAX(_0dl5e11[i*16 +1],      fabsf(_6fokzix[i*16 +7])); 
                            
                            _myx444n(_atj9tze, _vxpj4hn, i, _p422ta8, _6fokzix, _0dl5e11, 1);
                            
                        }
                        
                        _2292mtr[_qnsn2s8[i]*2 +0] = _kibu4fm[0];         _2292mtr[_qnsn2s8[i]*2 +1] = _kibu4fm[3];
                        _3lbqgu3[_qnsn2s8[i]*2 +0] = _kibu4fm[1];         _3lbqgu3[_qnsn2s8[i]*2 +1] = _kibu4fm[4];
                        _5bwikot[_qnsn2s8[i]*2 +0] = _kibu4fm[2];         _5bwikot[_qnsn2s8[i]*2 +1] = _kibu4fm[5];
                        (*_5gzpl1q)[i][_8djzr5s]       = _qnsn2s8[i];
                        _qnsn2s8[i]                  += 1u;
                        
                        _be3yqy5[_8djzr5s]     = 1u;
                        
                        _d0p92k3 = 0u;
                        for (j = 0u; j < _vxpj4hn->_3ru7myo;  j++)   {       _d0p92k3 +=  _be3yqy5[j];                                      }
                        
                        for (j = 1u; j < _omxlvz4; j++)             {       _t4q9bc2[ _2h2p39y[_8djzr5s][16 -_omxlvz4 +j] ] = 1u;        }
                        _qke2gjr = 0u;
                    }
                    else
                    {   _aq15o2n = _2h2p39y[_8djzr5s][16 -_omxlvz4 +_qw4c6lb];
                        
                        memcpy(_dan085r, &_e3nwxuj[_aq15o2n][9],  3*sizeof(float));
                        memcpy(_ntdr5xy, &_e3nwxuj[_aq15o2n][12], 3*sizeof(float));
                        memcpy(_fc4kl8c, &_e3nwxuj[_aq15o2n][15], 3*sizeof(float));
                        memcpy(_gn0t6dv, &_e3nwxuj[_aq15o2n][18], 3*sizeof(float));
                        
                        if (_s55j2o7->_s832ejx == 1u)
                        {   _oyr6ccv(_c9x5jz6, _vg5vmf6, _pt5r5rc[0], _pt5r5rc[1], _pt5r5rc[2], _dan085r, _ntdr5xy, _gn0t6dv, _s55j2o7->_ycz04cd, 0.0f, 0.0f, _atj9tze->_kkrzd5e, _atj9tze->_6snglh6);
                            _oyr6ccv(_1a6xvaj, _vg5vmf6, _pt5r5rc[0], _pt5r5rc[1], _pt5r5rc[2], _dan085r, _ntdr5xy, _gn0t6dv, 0.0f, _s55j2o7->_ycz04cd, 0.0f, _atj9tze->_kkrzd5e, _atj9tze->_6snglh6);
                            
                            
                            _oyr6ccv(_x4nxvkd,  _vg5vmf6, _pt5r5rc[0], _pt5r5rc[1], _pt5r5rc[2], _ntdr5xy, _fc4kl8c, _gn0t6dv, _s55j2o7->_ycz04cd, 0.0f, 0.0f, _atj9tze->_kkrzd5e, _atj9tze->_6snglh6);
                            _oyr6ccv(_4xyjoa2,  _vg5vmf6, _pt5r5rc[0], _pt5r5rc[1], _pt5r5rc[2], _ntdr5xy, _fc4kl8c, _gn0t6dv, 0.0f, _s55j2o7->_ycz04cd, 0.0f, _atj9tze->_kkrzd5e, _atj9tze->_6snglh6);
                            
                            
                            _x4nxvkd[0] += _c9x5jz6[0];        _x4nxvkd[1] += _c9x5jz6[1];        _x4nxvkd[2] += _c9x5jz6[2];
                            _x4nxvkd[3] += _c9x5jz6[3];        _x4nxvkd[4] += _c9x5jz6[4];        _x4nxvkd[5] += _c9x5jz6[5];
                            _4xyjoa2[0] += _1a6xvaj[0];        _4xyjoa2[1] += _1a6xvaj[1];        _4xyjoa2[2] += _1a6xvaj[2];
                            _4xyjoa2[3] += _1a6xvaj[3];        _4xyjoa2[4] += _1a6xvaj[4];        _4xyjoa2[5] += _1a6xvaj[5];
                            
                            
                        }
                        else
                        {   _25sq8cr(_c9x5jz6, _vg5vmf6, _pt5r5rc[0], _pt5r5rc[1], _pt5r5rc[2], _dan085r, _ntdr5xy, _gn0t6dv, _s55j2o7->_ycz04cd, 0.0f, 0.0f, _atj9tze->_kkrzd5e, _atj9tze->_6snglh6);
                            _25sq8cr(_1a6xvaj, _vg5vmf6, _pt5r5rc[0], _pt5r5rc[1], _pt5r5rc[2], _dan085r, _ntdr5xy, _gn0t6dv, 0.0f, _s55j2o7->_ycz04cd, 0.0f, _atj9tze->_kkrzd5e, _atj9tze->_6snglh6);
                            
                            
                            _25sq8cr(_x4nxvkd,  _vg5vmf6, _pt5r5rc[0], _pt5r5rc[1], _pt5r5rc[2], _ntdr5xy, _fc4kl8c, _gn0t6dv, _s55j2o7->_ycz04cd, 0.0f, 0.0f, _atj9tze->_kkrzd5e, _atj9tze->_6snglh6);
                            _25sq8cr(_4xyjoa2,  _vg5vmf6, _pt5r5rc[0], _pt5r5rc[1], _pt5r5rc[2], _ntdr5xy, _fc4kl8c, _gn0t6dv, 0.0f, _s55j2o7->_ycz04cd, 0.0f, _atj9tze->_kkrzd5e, _atj9tze->_6snglh6);
                            
                            
                            _x4nxvkd[0] += _c9x5jz6[0];        _x4nxvkd[1] += _c9x5jz6[1];        _x4nxvkd[2] += _c9x5jz6[2];
                            _x4nxvkd[3] += _c9x5jz6[3];        _x4nxvkd[4] += _c9x5jz6[4];        _x4nxvkd[5] += _c9x5jz6[5];
                            _4xyjoa2[0] += _1a6xvaj[0];        _4xyjoa2[1] += _1a6xvaj[1];        _4xyjoa2[2] += _1a6xvaj[2];
                            _4xyjoa2[3] += _1a6xvaj[3];        _4xyjoa2[4] += _1a6xvaj[4];        _4xyjoa2[5] += _1a6xvaj[5];
                            
                            
                        }
                        
                        _yyv37di(_x4nxvkd, (_qw5dsk7/(_s55j2o7->_ycz04cd*_e3nwxuj[_aq15o2n][3])));          _ftu2gc5(&_kibu4fm[0], _p9p9w4j, _x4nxvkd);
                        _yyv37di(_4xyjoa2, (_qw5dsk7/(_s55j2o7->_ycz04cd*_e3nwxuj[_aq15o2n][3])));          _ftu2gc5(&_kibu4fm[3], _p9p9w4j, _4xyjoa2);
                        
                        
                        _81bmqf7 = 1u;
                        for (j = 0; j < _2qgv3dr[_aq15o2n][0]; j++)
                        {   _4nnkm95 = _2qgv3dr[_aq15o2n][(j+1)];
                            
                            memcpy(_dan085r, &_e3nwxuj[_4nnkm95][9],  3*sizeof(float));
                            memcpy(_ntdr5xy, &_e3nwxuj[_4nnkm95][12], 3*sizeof(float));
                            memcpy(_fc4kl8c, &_e3nwxuj[_4nnkm95][15], 3*sizeof(float));
                            memcpy(_gn0t6dv, &_e3nwxuj[_4nnkm95][18], 3*sizeof(float));
                            
                            if (_s55j2o7->_s832ejx == 1u)
                            {   _oyr6ccv(_c9x5jz6, _vg5vmf6, _pt5r5rc[0], _pt5r5rc[1], _pt5r5rc[2], _dan085r, _ntdr5xy, _gn0t6dv, _s55j2o7->_ycz04cd, 0.0f, 0.0f, _atj9tze->_kkrzd5e, _atj9tze->_6snglh6);
                                _oyr6ccv(_1a6xvaj, _vg5vmf6, _pt5r5rc[0], _pt5r5rc[1], _pt5r5rc[2], _dan085r, _ntdr5xy, _gn0t6dv, 0.0f, _s55j2o7->_ycz04cd, 0.0f, _atj9tze->_kkrzd5e, _atj9tze->_6snglh6);
                                
                                
                                _oyr6ccv(_x4nxvkd,  _vg5vmf6, _pt5r5rc[0], _pt5r5rc[1], _pt5r5rc[2], _ntdr5xy, _fc4kl8c, _gn0t6dv, _s55j2o7->_ycz04cd, 0.0f, 0.0f, _atj9tze->_kkrzd5e, _atj9tze->_6snglh6);
                                _oyr6ccv(_4xyjoa2,  _vg5vmf6, _pt5r5rc[0], _pt5r5rc[1], _pt5r5rc[2], _ntdr5xy, _fc4kl8c, _gn0t6dv, 0.0f, _s55j2o7->_ycz04cd, 0.0f, _atj9tze->_kkrzd5e, _atj9tze->_6snglh6);
                                
                                
                                _x4nxvkd[0] += _c9x5jz6[0];        _x4nxvkd[1] += _c9x5jz6[1];        _x4nxvkd[2] += _c9x5jz6[2];
                                _x4nxvkd[3] += _c9x5jz6[3];        _x4nxvkd[4] += _c9x5jz6[4];        _x4nxvkd[5] += _c9x5jz6[5];
                                _4xyjoa2[0] += _1a6xvaj[0];        _4xyjoa2[1] += _1a6xvaj[1];        _4xyjoa2[2] += _1a6xvaj[2];
                                _4xyjoa2[3] += _1a6xvaj[3];        _4xyjoa2[4] += _1a6xvaj[4];        _4xyjoa2[5] += _1a6xvaj[5];
                                
                                
                            }
                            else
                            {   _25sq8cr(_c9x5jz6, _vg5vmf6, _pt5r5rc[0], _pt5r5rc[1], _pt5r5rc[2], _dan085r, _ntdr5xy, _gn0t6dv, _s55j2o7->_ycz04cd, 0.0f, 0.0f, _atj9tze->_kkrzd5e, _atj9tze->_6snglh6);
                                _25sq8cr(_1a6xvaj, _vg5vmf6, _pt5r5rc[0], _pt5r5rc[1], _pt5r5rc[2], _dan085r, _ntdr5xy, _gn0t6dv, 0.0f, _s55j2o7->_ycz04cd, 0.0f, _atj9tze->_kkrzd5e, _atj9tze->_6snglh6);
                                
                                
                                _25sq8cr(_x4nxvkd,  _vg5vmf6, _pt5r5rc[0], _pt5r5rc[1], _pt5r5rc[2], _ntdr5xy, _fc4kl8c, _gn0t6dv, _s55j2o7->_ycz04cd, 0.0f, 0.0f, _atj9tze->_kkrzd5e, _atj9tze->_6snglh6);
                                _25sq8cr(_4xyjoa2,  _vg5vmf6, _pt5r5rc[0], _pt5r5rc[1], _pt5r5rc[2], _ntdr5xy, _fc4kl8c, _gn0t6dv, 0.0f, _s55j2o7->_ycz04cd, 0.0f, _atj9tze->_kkrzd5e, _atj9tze->_6snglh6);
                                
                                
                                _x4nxvkd[0] += _c9x5jz6[0];        _x4nxvkd[1] += _c9x5jz6[1];        _x4nxvkd[2] += _c9x5jz6[2];
                                _x4nxvkd[3] += _c9x5jz6[3];        _x4nxvkd[4] += _c9x5jz6[4];        _x4nxvkd[5] += _c9x5jz6[5];
                                _4xyjoa2[0] += _1a6xvaj[0];        _4xyjoa2[1] += _1a6xvaj[1];        _4xyjoa2[2] += _1a6xvaj[2];
                                _4xyjoa2[3] += _1a6xvaj[3];        _4xyjoa2[4] += _1a6xvaj[4];        _4xyjoa2[5] += _1a6xvaj[5];
                                
                                
                            }
                            
                            _yyv37di(_x4nxvkd, (_qw5dsk7/(_s55j2o7->_ycz04cd*_e3nwxuj[_4nnkm95][3])));       _ftu2gc5(&_zy1uu9y[0], _p9p9w4j, _x4nxvkd);
                            _yyv37di(_4xyjoa2, (_qw5dsk7/(_s55j2o7->_ycz04cd*_e3nwxuj[_4nnkm95][3])));       _ftu2gc5(&_zy1uu9y[3], _p9p9w4j, _4xyjoa2);
                            
                            
                            memcpy(_9p3jpxg, &_e3nwxuj[_4nnkm95][0], 3*sizeof(float));
                            _id4n1u3(_n42x2ez, _pt5r5rc, _9p3jpxg);                              _31yrujo  = _mabsi9n(_n42x2ez);
                            _gymdm0l = (_31yrujo/_s55j2o7->_b68rerv - _s55j2o7->_a11dpg8)/(_s55j2o7->_sbi4q2y - _s55j2o7->_a11dpg8);
                            _gymdm0l = (_gymdm0l >= 0.0f)*_gymdm0l  +  (_gymdm0l < 0.0f)*0.0f;
                            _gymdm0l = (_gymdm0l <= 1.0f)*_gymdm0l  +  (_gymdm0l > 1.0f)*1.0f;
                            _gymdm0l = _s55j2o7->_9v86wzs*_gymdm0l;
                            
                            _t26ys44 = sqrtf( (_zy1uu9y[0]-_kibu4fm[0])*(_zy1uu9y[0]-_kibu4fm[0]) + (_zy1uu9y[1]-_kibu4fm[1])*(_zy1uu9y[1]-_kibu4fm[1]) + (_zy1uu9y[2]-_kibu4fm[2])*(_zy1uu9y[2]-_kibu4fm[2])  + (_zy1uu9y[3]-_kibu4fm[3])*(_zy1uu9y[3]-_kibu4fm[3]) + (_zy1uu9y[4]-_kibu4fm[4])*(_zy1uu9y[4]-_kibu4fm[4]) + (_zy1uu9y[5]-_kibu4fm[5])*(_zy1uu9y[5]-_kibu4fm[5]) );
                            _5gdej31 = sqrtf(                _kibu4fm[0] *_kibu4fm[0]                 +                _kibu4fm[1] *_kibu4fm[1]                 +                _kibu4fm[2] *_kibu4fm[2]                  +                _kibu4fm[3] *_kibu4fm[3]                 +                _kibu4fm[4] *_kibu4fm[4]                 +                _kibu4fm[5] *_kibu4fm[5] );
                            
                            _81bmqf7 = ((_t26ys44/_5gdej31) <= _gymdm0l)*_81bmqf7  +  ((_t26ys44/_5gdej31) > _gymdm0l)*0u; 
                        }
                        
                        if (_81bmqf7 == 1u)
                        {   _2292mtr[_qnsn2s8[i]*2 +0] = _kibu4fm[0];         _2292mtr[_qnsn2s8[i]*2 +1] = _kibu4fm[3];
                            _3lbqgu3[_qnsn2s8[i]*2 +0] = _kibu4fm[1];         _3lbqgu3[_qnsn2s8[i]*2 +1] = _kibu4fm[4];
                            _5bwikot[_qnsn2s8[i]*2 +0] = _kibu4fm[2];         _5bwikot[_qnsn2s8[i]*2 +1] = _kibu4fm[5];
                            
                            _d0p92k3 = 0u;
                            for (j = 0u; j < _vxpj4hn->_3ru7myo;  j++)
                            {   _be3yqy5[j]        = (_2h2p39y[j][16 -_omxlvz4 +_qw4c6lb] == _aq15o2n)*1u            +  (_2h2p39y[j][16 -_omxlvz4 +_qw4c6lb] != _aq15o2n)*_be3yqy5[j];
                                (*_5gzpl1q)[i][j] = (_2h2p39y[j][16 -_omxlvz4 +_qw4c6lb] == _aq15o2n)*_qnsn2s8[i]  +  (_2h2p39y[j][16 -_omxlvz4 +_qw4c6lb] != _aq15o2n)*(*_5gzpl1q)[i][j];
                                _d0p92k3          += _be3yqy5[j];
                            }
                            _qnsn2s8[i] += 1u;
                            
                            for (j = 0u; j < _omxlvz4; j++)       {       _t4q9bc2[ _2h2p39y[_8djzr5s][16 -_omxlvz4 +j] ] = 1u;       }
                            _qke2gjr = 0u;
                        }
                    _qw4c6lb += 1u;
            }   }   }
            
            _moxsfh4       = _z0yncxl(ALIGN64, 2*_qnsn2s8[i], sizeof(float));
            (*_ckeh5sw)[i] = aligned_alloc(ALIGN64, _moxsfh4);       memcpy((*_ckeh5sw)[i], _2292mtr, 2*_qnsn2s8[i]* sizeof(float));
            (*_3ac6pe6)[i] = aligned_alloc(ALIGN64, _moxsfh4);       memcpy((*_3ac6pe6)[i], _3lbqgu3, 2*_qnsn2s8[i]* sizeof(float));
            (*_obyghgk)[i] = aligned_alloc(ALIGN64, _moxsfh4);       memcpy((*_obyghgk)[i], _5bwikot, 2*_qnsn2s8[i]* sizeof(float));
            
            if ((_vxpj4hn->_qdbxvh7 == 1)&&(_fsqg8x8 == 0))     {    fprintf(stdout,"\nElement %d of %d    with     %10u resolved sources", i, _lexdw6l[_fsqg8x8],  _qnsn2s8[i]);                  }
            
            
            if (_vxpj4hn->_5izolnn > 0u)
            {   _moxsfh4   = _z0yncxl(ALIGN64, 3*_vxpj4hn->_5izolnn, sizeof(float));
                memset(_iq2fa39, 0, _moxsfh4);
                memset(_b4ldims, 0, _moxsfh4);
                
                memset(_otn9uze,   0, _vxpj4hn->_5izolnn*sizeof(unsigned int));
                memset(_mgvu4is, 0, _kb700uc*sizeof(unsigned int));
                _d0p92k3  = 0u;
                
                while (_d0p92k3 < _vxpj4hn->_5izolnn) 
                {   
                    for (j = 0u; j < _vxpj4hn->_5izolnn; j++)               {       if (_otn9uze[j] == 0u)    {      _8djzr5s = j;       break;          }           } 
                    
                    memcpy(_9p3jpxg, &_rrpzlk3[_8djzr5s*16 +0], 3u*sizeof(float));
                    
                    _id4n1u3(_n42x2ez, _pt5r5rc, _9p3jpxg);               _31yrujo  = _mabsi9n(_n42x2ez);
                    
                    _qw5dsk7 = (_31yrujo >= fMinDistFact*_s55j2o7->_gfu8iog)*1.0f  +  (_31yrujo < fMinDistFact*_s55j2o7->_gfu8iog)*0.0f; 
                    
                    _omxlvz4 = _qsq0p82[_8djzr5s][0]; 
                    
                    _qw4c6lb = 0;
                    _qke2gjr   = 1u;
                    while (_qke2gjr == 1u)
                    {   
                        _81bmqf7 = 0u;
                        while (_qw4c6lb < _omxlvz4-1)
                        {   _aq15o2n = _qsq0p82[_8djzr5s][16 -_omxlvz4 +_qw4c6lb];
                            
                            memcpy(_9p3jpxg, &_nf74t60[_aq15o2n][0], 3*sizeof(float));
                            _id4n1u3(_n42x2ez, _pt5r5rc, _9p3jpxg);               _31yrujo  = _mabsi9n(_n42x2ez);
                            
                            _81bmqf7  = (_mgvu4is[_aq15o2n]                                                 == 0u)*1u      +  (_mgvu4is[_aq15o2n]                                                  != 0u)*0u;
                            _81bmqf7  = (_cw2lbm2[_aq15o2n][0]                                              > 1u)*_81bmqf7  +  (_cw2lbm2[_aq15o2n][0]                                              <= 1u)*0u;
                            _81bmqf7  = (_31yrujo  > (_nf74t60[_aq15o2n][5] +_s55j2o7->_a11dpg8*_s55j2o7->_gfu8iog))*_81bmqf7  +  (_31yrujo  <= (_nf74t60[_aq15o2n][5] +_s55j2o7->_a11dpg8*_s55j2o7->_gfu8iog))*0u;
                            
                            if (_81bmqf7 == 1u)         {       break;      }
                            
                            _qw4c6lb += 1u;
                        } 
                        
                        if (_81bmqf7 == 0u) 
                        {   memcpy(_dan085r, &_3n3bzt9[_rgb93wk[_8djzr5s*8 +0]*4 +0], 3*sizeof(float));
                            memcpy(_ntdr5xy, &_3n3bzt9[_rgb93wk[_8djzr5s*8 +1]*4 +0], 3*sizeof(float));
                            memcpy(_fc4kl8c, &_3n3bzt9[_rgb93wk[_8djzr5s*8 +2]*4 +0], 3*sizeof(float));
                            
                            if (_s55j2o7->_s832ejx == 1u)
                            {   _oyr6ccv(_x4nxvkd, _vg5vmf6, _pt5r5rc[0], _pt5r5rc[1], _pt5r5rc[2], _dan085r, _ntdr5xy, _fc4kl8c, _s55j2o7->_cb1kjtg, 0.0f, 0.0f, _atj9tze->_kkrzd5e, _atj9tze->_6snglh6);
                                _oyr6ccv(_4xyjoa2, _vg5vmf6, _pt5r5rc[0], _pt5r5rc[1], _pt5r5rc[2], _dan085r, _ntdr5xy, _fc4kl8c, 0.0f, _s55j2o7->_cb1kjtg, 0.0f, _atj9tze->_kkrzd5e, _atj9tze->_6snglh6);
                                _oyr6ccv(_2k8xpnm, _vg5vmf6, _pt5r5rc[0], _pt5r5rc[1], _pt5r5rc[2], _dan085r, _ntdr5xy, _fc4kl8c, 0.0f, 0.0f, _s55j2o7->_cb1kjtg, _atj9tze->_kkrzd5e, _atj9tze->_6snglh6);
                            }
                            else
                            {   _25sq8cr(_x4nxvkd, _vg5vmf6, _pt5r5rc[0], _pt5r5rc[1], _pt5r5rc[2], _dan085r, _ntdr5xy, _fc4kl8c, _s55j2o7->_cb1kjtg, 0.0f, 0.0f, _atj9tze->_kkrzd5e, _atj9tze->_6snglh6);
                                _25sq8cr(_4xyjoa2, _vg5vmf6, _pt5r5rc[0], _pt5r5rc[1], _pt5r5rc[2], _dan085r, _ntdr5xy, _fc4kl8c, 0.0f, _s55j2o7->_cb1kjtg, 0.0f, _atj9tze->_kkrzd5e, _atj9tze->_6snglh6);
                                _25sq8cr(_2k8xpnm, _vg5vmf6, _pt5r5rc[0], _pt5r5rc[1], _pt5r5rc[2], _dan085r, _ntdr5xy, _fc4kl8c, 0.0f, 0.0f, _s55j2o7->_cb1kjtg, _atj9tze->_kkrzd5e, _atj9tze->_6snglh6);
                            }
                            
                            _yyv37di(_x4nxvkd, (_qw5dsk7/(_s55j2o7->_cb1kjtg*_rrpzlk3[_8djzr5s*16 +15])));      _ftu2gc5(&_kibu4fm[0], _p9p9w4j, _x4nxvkd);
                            _yyv37di(_4xyjoa2, (_qw5dsk7/(_s55j2o7->_cb1kjtg*_rrpzlk3[_8djzr5s*16 +15])));      _ftu2gc5(&_kibu4fm[3], _p9p9w4j, _4xyjoa2);
                            _yyv37di(_2k8xpnm, (_qw5dsk7/(_s55j2o7->_cb1kjtg*_rrpzlk3[_8djzr5s*16 +15])));      _ftu2gc5(&_kibu4fm[6], _p9p9w4j, _2k8xpnm);
                            
                            _iq2fa39[_qpwbsf3[i]*3 +0] = _kibu4fm[0];         _iq2fa39[_qpwbsf3[i]*3 +1] = _kibu4fm[3];         _iq2fa39[_qpwbsf3[i]*3 +2] = _kibu4fm[6];
                            _b4ldims[_qpwbsf3[i]*3 +0] = _kibu4fm[1];         _b4ldims[_qpwbsf3[i]*3 +1] = _kibu4fm[4];         _b4ldims[_qpwbsf3[i]*3 +2] = _kibu4fm[7];
                            (*_oqgvzhq)[i][_8djzr5s]       = _qpwbsf3[i];
                            _qpwbsf3[i]                  += 1u;
                            
                            _otn9uze[_8djzr5s]     = 1u;
                            
                            _d0p92k3 = 0u;
                            for (j = 0u; j < _vxpj4hn->_5izolnn;  j++)   {       _d0p92k3 +=  _otn9uze[j];                                     }
                            
                            for (j = 0u; j < _omxlvz4; j++)             {       _mgvu4is[ _qsq0p82[_8djzr5s][16 -_omxlvz4 +j] ] = 1u;       }
                            _qke2gjr = 0u;
                        }
                        else
                        {   _aq15o2n = _qsq0p82[_8djzr5s][16 -_omxlvz4 +_qw4c6lb];
                            
                            memcpy(_dan085r, &_nf74t60[_aq15o2n][9],  3*sizeof(float));
                            memcpy(_ntdr5xy, &_nf74t60[_aq15o2n][12], 3*sizeof(float));
                            memcpy(_fc4kl8c, &_nf74t60[_aq15o2n][15], 3*sizeof(float));
                            memcpy(_gn0t6dv, &_nf74t60[_aq15o2n][18], 3*sizeof(float));
                            
                            if (_s55j2o7->_s832ejx == 1u)
                            {   _oyr6ccv(_c9x5jz6, _vg5vmf6, _pt5r5rc[0], _pt5r5rc[1], _pt5r5rc[2], _dan085r, _ntdr5xy, _gn0t6dv, _s55j2o7->_cb1kjtg, 0.0f, 0.0f, _atj9tze->_kkrzd5e, _atj9tze->_6snglh6);
                                _oyr6ccv(_1a6xvaj, _vg5vmf6, _pt5r5rc[0], _pt5r5rc[1], _pt5r5rc[2], _dan085r, _ntdr5xy, _gn0t6dv, 0.0f, _s55j2o7->_cb1kjtg, 0.0f, _atj9tze->_kkrzd5e, _atj9tze->_6snglh6);
                                _oyr6ccv(_jqy3ods, _vg5vmf6, _pt5r5rc[0], _pt5r5rc[1], _pt5r5rc[2], _dan085r, _ntdr5xy, _gn0t6dv, 0.0f, 0.0f, _s55j2o7->_cb1kjtg, _atj9tze->_kkrzd5e, _atj9tze->_6snglh6);
                                
                                _oyr6ccv(_x4nxvkd,  _vg5vmf6, _pt5r5rc[0], _pt5r5rc[1], _pt5r5rc[2], _ntdr5xy, _fc4kl8c, _gn0t6dv, _s55j2o7->_cb1kjtg, 0.0f, 0.0f, _atj9tze->_kkrzd5e, _atj9tze->_6snglh6);
                                _oyr6ccv(_4xyjoa2,  _vg5vmf6, _pt5r5rc[0], _pt5r5rc[1], _pt5r5rc[2], _ntdr5xy, _fc4kl8c, _gn0t6dv, 0.0f, _s55j2o7->_cb1kjtg, 0.0f, _atj9tze->_kkrzd5e, _atj9tze->_6snglh6);
                                _oyr6ccv(_2k8xpnm,  _vg5vmf6, _pt5r5rc[0], _pt5r5rc[1], _pt5r5rc[2], _ntdr5xy, _fc4kl8c, _gn0t6dv, 0.0f, 0.0f, _s55j2o7->_cb1kjtg, _atj9tze->_kkrzd5e, _atj9tze->_6snglh6);
                                
                                _x4nxvkd[0] += _c9x5jz6[0];        _x4nxvkd[1] += _c9x5jz6[1];        _x4nxvkd[2] += _c9x5jz6[2];
                                _x4nxvkd[3] += _c9x5jz6[3];        _x4nxvkd[4] += _c9x5jz6[4];        _x4nxvkd[5] += _c9x5jz6[5];
                                _4xyjoa2[0] += _1a6xvaj[0];        _4xyjoa2[1] += _1a6xvaj[1];        _4xyjoa2[2] += _1a6xvaj[2];
                                _4xyjoa2[3] += _1a6xvaj[3];        _4xyjoa2[4] += _1a6xvaj[4];        _4xyjoa2[5] += _1a6xvaj[5];
                                _2k8xpnm[0] += _jqy3ods[0];        _2k8xpnm[1] += _jqy3ods[1];        _2k8xpnm[2] += _jqy3ods[2];
                                _2k8xpnm[3] += _jqy3ods[3];        _2k8xpnm[4] += _jqy3ods[4];        _2k8xpnm[5] += _jqy3ods[5];
                            }
                            else
                            {   _25sq8cr(_c9x5jz6, _vg5vmf6, _pt5r5rc[0], _pt5r5rc[1], _pt5r5rc[2], _dan085r, _ntdr5xy, _gn0t6dv, _s55j2o7->_cb1kjtg, 0.0f, 0.0f, _atj9tze->_kkrzd5e, _atj9tze->_6snglh6);
                                _25sq8cr(_1a6xvaj, _vg5vmf6, _pt5r5rc[0], _pt5r5rc[1], _pt5r5rc[2], _dan085r, _ntdr5xy, _gn0t6dv, 0.0f, _s55j2o7->_cb1kjtg, 0.0f, _atj9tze->_kkrzd5e, _atj9tze->_6snglh6);
                                _25sq8cr(_jqy3ods, _vg5vmf6, _pt5r5rc[0], _pt5r5rc[1], _pt5r5rc[2], _dan085r, _ntdr5xy, _gn0t6dv, 0.0f, 0.0f, _s55j2o7->_cb1kjtg, _atj9tze->_kkrzd5e, _atj9tze->_6snglh6);
                                
                                _25sq8cr(_x4nxvkd,  _vg5vmf6, _pt5r5rc[0], _pt5r5rc[1], _pt5r5rc[2], _ntdr5xy, _fc4kl8c, _gn0t6dv, _s55j2o7->_cb1kjtg, 0.0f, 0.0f, _atj9tze->_kkrzd5e, _atj9tze->_6snglh6);
                                _25sq8cr(_4xyjoa2,  _vg5vmf6, _pt5r5rc[0], _pt5r5rc[1], _pt5r5rc[2], _ntdr5xy, _fc4kl8c, _gn0t6dv, 0.0f, _s55j2o7->_cb1kjtg, 0.0f, _atj9tze->_kkrzd5e, _atj9tze->_6snglh6);
                                _25sq8cr(_2k8xpnm,  _vg5vmf6, _pt5r5rc[0], _pt5r5rc[1], _pt5r5rc[2], _ntdr5xy, _fc4kl8c, _gn0t6dv, 0.0f, 0.0f, _s55j2o7->_cb1kjtg, _atj9tze->_kkrzd5e, _atj9tze->_6snglh6);
                                
                                _x4nxvkd[0] += _c9x5jz6[0];        _x4nxvkd[1] += _c9x5jz6[1];        _x4nxvkd[2] += _c9x5jz6[2];
                                _x4nxvkd[3] += _c9x5jz6[3];        _x4nxvkd[4] += _c9x5jz6[4];        _x4nxvkd[5] += _c9x5jz6[5];
                                _4xyjoa2[0] += _1a6xvaj[0];        _4xyjoa2[1] += _1a6xvaj[1];        _4xyjoa2[2] += _1a6xvaj[2];
                                _4xyjoa2[3] += _1a6xvaj[3];        _4xyjoa2[4] += _1a6xvaj[4];        _4xyjoa2[5] += _1a6xvaj[5];
                                _2k8xpnm[0] += _jqy3ods[0];        _2k8xpnm[1] += _jqy3ods[1];        _2k8xpnm[2] += _jqy3ods[2];
                                _2k8xpnm[3] += _jqy3ods[3];        _2k8xpnm[4] += _jqy3ods[4];        _2k8xpnm[5] += _jqy3ods[5];
                            }
                            
                            _yyv37di(_x4nxvkd, (_qw5dsk7/(_s55j2o7->_cb1kjtg*_nf74t60[_aq15o2n][3])));          _ftu2gc5(&_kibu4fm[0], _p9p9w4j, _x4nxvkd);
                            _yyv37di(_4xyjoa2, (_qw5dsk7/(_s55j2o7->_cb1kjtg*_nf74t60[_aq15o2n][3])));          _ftu2gc5(&_kibu4fm[3], _p9p9w4j, _4xyjoa2);
                            _yyv37di(_2k8xpnm, (_qw5dsk7/(_s55j2o7->_cb1kjtg*_nf74t60[_aq15o2n][3])));          _ftu2gc5(&_kibu4fm[6], _p9p9w4j, _2k8xpnm);
                            
                            _81bmqf7 = 1u;
                            for (j = 0; j < _cw2lbm2[_aq15o2n][0]; j++)
                            {   _4nnkm95 = _cw2lbm2[_aq15o2n][(j+1)];
                                
                                memcpy(_dan085r, &_nf74t60[_4nnkm95][9],  3*sizeof(float));
                                memcpy(_ntdr5xy, &_nf74t60[_4nnkm95][12], 3*sizeof(float));
                                memcpy(_fc4kl8c, &_nf74t60[_4nnkm95][15], 3*sizeof(float));
                                memcpy(_gn0t6dv, &_nf74t60[_4nnkm95][18], 3*sizeof(float));
                                
                                if (_s55j2o7->_s832ejx == 1u)
                                {   _oyr6ccv(_c9x5jz6, _vg5vmf6, _pt5r5rc[0], _pt5r5rc[1], _pt5r5rc[2], _dan085r, _ntdr5xy, _gn0t6dv, _s55j2o7->_cb1kjtg, 0.0f, 0.0f, _atj9tze->_kkrzd5e, _atj9tze->_6snglh6);
                                    _oyr6ccv(_1a6xvaj, _vg5vmf6, _pt5r5rc[0], _pt5r5rc[1], _pt5r5rc[2], _dan085r, _ntdr5xy, _gn0t6dv, 0.0f, _s55j2o7->_cb1kjtg, 0.0f, _atj9tze->_kkrzd5e, _atj9tze->_6snglh6);
                                    _oyr6ccv(_jqy3ods, _vg5vmf6, _pt5r5rc[0], _pt5r5rc[1], _pt5r5rc[2], _dan085r, _ntdr5xy, _gn0t6dv, 0.0f, 0.0f, _s55j2o7->_cb1kjtg, _atj9tze->_kkrzd5e, _atj9tze->_6snglh6);
                                    
                                    _oyr6ccv(_x4nxvkd,  _vg5vmf6, _pt5r5rc[0], _pt5r5rc[1], _pt5r5rc[2], _ntdr5xy, _fc4kl8c, _gn0t6dv, _s55j2o7->_cb1kjtg, 0.0f, 0.0f, _atj9tze->_kkrzd5e, _atj9tze->_6snglh6);
                                    _oyr6ccv(_4xyjoa2,  _vg5vmf6, _pt5r5rc[0], _pt5r5rc[1], _pt5r5rc[2], _ntdr5xy, _fc4kl8c, _gn0t6dv, 0.0f, _s55j2o7->_cb1kjtg, 0.0f, _atj9tze->_kkrzd5e, _atj9tze->_6snglh6);
                                    _oyr6ccv(_2k8xpnm,  _vg5vmf6, _pt5r5rc[0], _pt5r5rc[1], _pt5r5rc[2], _ntdr5xy, _fc4kl8c, _gn0t6dv, 0.0f, 0.0f, _s55j2o7->_cb1kjtg, _atj9tze->_kkrzd5e, _atj9tze->_6snglh6);
                                    
                                    _x4nxvkd[0] += _c9x5jz6[0];        _x4nxvkd[1] += _c9x5jz6[1];        _x4nxvkd[2] += _c9x5jz6[2];
                                    _x4nxvkd[3] += _c9x5jz6[3];        _x4nxvkd[4] += _c9x5jz6[4];        _x4nxvkd[5] += _c9x5jz6[5];
                                    _4xyjoa2[0] += _1a6xvaj[0];        _4xyjoa2[1] += _1a6xvaj[1];        _4xyjoa2[2] += _1a6xvaj[2];
                                    _4xyjoa2[3] += _1a6xvaj[3];        _4xyjoa2[4] += _1a6xvaj[4];        _4xyjoa2[5] += _1a6xvaj[5];
                                    _2k8xpnm[0] += _jqy3ods[0];        _2k8xpnm[1] += _jqy3ods[1];        _2k8xpnm[2] += _jqy3ods[2];
                                    _2k8xpnm[3] += _jqy3ods[3];        _2k8xpnm[4] += _jqy3ods[4];        _2k8xpnm[5] += _jqy3ods[5];
                                }
                                else
                                {   _25sq8cr(_c9x5jz6, _vg5vmf6, _pt5r5rc[0], _pt5r5rc[1], _pt5r5rc[2], _dan085r, _ntdr5xy, _gn0t6dv, _s55j2o7->_cb1kjtg, 0.0f, 0.0f, _atj9tze->_kkrzd5e, _atj9tze->_6snglh6);
                                    _25sq8cr(_1a6xvaj, _vg5vmf6, _pt5r5rc[0], _pt5r5rc[1], _pt5r5rc[2], _dan085r, _ntdr5xy, _gn0t6dv, 0.0f, _s55j2o7->_cb1kjtg, 0.0f, _atj9tze->_kkrzd5e, _atj9tze->_6snglh6);
                                    _25sq8cr(_jqy3ods, _vg5vmf6, _pt5r5rc[0], _pt5r5rc[1], _pt5r5rc[2], _dan085r, _ntdr5xy, _gn0t6dv, 0.0f, 0.0f, _s55j2o7->_cb1kjtg, _atj9tze->_kkrzd5e, _atj9tze->_6snglh6);
                                    
                                    _25sq8cr(_x4nxvkd,  _vg5vmf6, _pt5r5rc[0], _pt5r5rc[1], _pt5r5rc[2], _ntdr5xy, _fc4kl8c, _gn0t6dv, _s55j2o7->_cb1kjtg, 0.0f, 0.0f, _atj9tze->_kkrzd5e, _atj9tze->_6snglh6);
                                    _25sq8cr(_4xyjoa2,  _vg5vmf6, _pt5r5rc[0], _pt5r5rc[1], _pt5r5rc[2], _ntdr5xy, _fc4kl8c, _gn0t6dv, 0.0f, _s55j2o7->_cb1kjtg, 0.0f, _atj9tze->_kkrzd5e, _atj9tze->_6snglh6);
                                    _25sq8cr(_2k8xpnm,  _vg5vmf6, _pt5r5rc[0], _pt5r5rc[1], _pt5r5rc[2], _ntdr5xy, _fc4kl8c, _gn0t6dv, 0.0f, 0.0f, _s55j2o7->_cb1kjtg, _atj9tze->_kkrzd5e, _atj9tze->_6snglh6);
                                    
                                    _x4nxvkd[0] += _c9x5jz6[0];        _x4nxvkd[1] += _c9x5jz6[1];        _x4nxvkd[2] += _c9x5jz6[2];
                                    _x4nxvkd[3] += _c9x5jz6[3];        _x4nxvkd[4] += _c9x5jz6[4];        _x4nxvkd[5] += _c9x5jz6[5];
                                    _4xyjoa2[0] += _1a6xvaj[0];        _4xyjoa2[1] += _1a6xvaj[1];        _4xyjoa2[2] += _1a6xvaj[2];
                                    _4xyjoa2[3] += _1a6xvaj[3];        _4xyjoa2[4] += _1a6xvaj[4];        _4xyjoa2[5] += _1a6xvaj[5];
                                    _2k8xpnm[0] += _jqy3ods[0];        _2k8xpnm[1] += _jqy3ods[1];        _2k8xpnm[2] += _jqy3ods[2];
                                    _2k8xpnm[3] += _jqy3ods[3];        _2k8xpnm[4] += _jqy3ods[4];        _2k8xpnm[5] += _jqy3ods[5];
                                }
                                
                                _yyv37di(_x4nxvkd, (_qw5dsk7/(_s55j2o7->_cb1kjtg*_nf74t60[_4nnkm95][3])));      _ftu2gc5(&_zy1uu9y[0], _p9p9w4j, _x4nxvkd);
                                _yyv37di(_4xyjoa2, (_qw5dsk7/(_s55j2o7->_cb1kjtg*_nf74t60[_4nnkm95][3])));      _ftu2gc5(&_zy1uu9y[3], _p9p9w4j, _4xyjoa2);
                                _yyv37di(_2k8xpnm, (_qw5dsk7/(_s55j2o7->_cb1kjtg*_nf74t60[_4nnkm95][3])));      _ftu2gc5(&_zy1uu9y[6], _p9p9w4j, _2k8xpnm);
                                
                                memcpy(_9p3jpxg, &_nf74t60[_4nnkm95][0], 3*sizeof(float));
                                _id4n1u3(_n42x2ez, _pt5r5rc, _9p3jpxg);                              _31yrujo  = _mabsi9n(_n42x2ez);
                                _gymdm0l = (_31yrujo/_s55j2o7->_gfu8iog - _s55j2o7->_a11dpg8)/(_s55j2o7->_sbi4q2y - _s55j2o7->_a11dpg8);
                                _gymdm0l = (_gymdm0l >= 0.0f)*_gymdm0l + (_gymdm0l < 0.0f)*0.0f;
                                _gymdm0l = (_gymdm0l <= 1.0f)*_gymdm0l + (_gymdm0l > 1.0f)*1.0f;
                                _gymdm0l = _s55j2o7->_9v86wzs*_gymdm0l;
                                
                                _t26ys44 = sqrtf( (_zy1uu9y[0]-_kibu4fm[0])*(_zy1uu9y[0]-_kibu4fm[0]) + (_zy1uu9y[1]-_kibu4fm[1])*(_zy1uu9y[1]-_kibu4fm[1]) + (_zy1uu9y[3]-_kibu4fm[3])*(_zy1uu9y[3]-_kibu4fm[3]) + (_zy1uu9y[4]-_kibu4fm[4])*(_zy1uu9y[4]-_kibu4fm[4]) + (_zy1uu9y[6]-_kibu4fm[6])*(_zy1uu9y[6]-_kibu4fm[6]) + (_zy1uu9y[7]-_kibu4fm[7])*(_zy1uu9y[7]-_kibu4fm[7]) );
                                _5gdej31 = sqrtf(                _kibu4fm[0] *_kibu4fm[0]                 +                _kibu4fm[1] *_kibu4fm[1]                 +                _kibu4fm[3] *_kibu4fm[3]                 +                _kibu4fm[4] *_kibu4fm[4]                 +                _kibu4fm[6] *_kibu4fm[6]                 +                _kibu4fm[7] *_kibu4fm[7] );
                                
                                _81bmqf7 = ((_t26ys44/_5gdej31) <= _gymdm0l)*_81bmqf7  +  ((_t26ys44/_5gdej31) > _gymdm0l)*0u; 
                            }
                            
                            if (_81bmqf7 == 1u)
                            {   _iq2fa39[_qpwbsf3[i]*3 +0] = _kibu4fm[0];         _iq2fa39[_qpwbsf3[i]*3 +1] = _kibu4fm[3];         _iq2fa39[_qpwbsf3[i]*3 +2] = _kibu4fm[6];
                                _b4ldims[_qpwbsf3[i]*3 +0] = _kibu4fm[1];         _b4ldims[_qpwbsf3[i]*3 +1] = _kibu4fm[4];         _b4ldims[_qpwbsf3[i]*3 +2] = _kibu4fm[7];
                                
                                _d0p92k3 = 0u;
                                for (j = 0u; j < _vxpj4hn->_5izolnn;  j++)
                                {   _otn9uze[j]        = (_qsq0p82[j][16 -_omxlvz4 +_qw4c6lb] == _aq15o2n)*1u            +  (_qsq0p82[j][16 -_omxlvz4 +_qw4c6lb] != _aq15o2n)*_otn9uze[j];
                                    (*_oqgvzhq)[i][j] = (_qsq0p82[j][16 -_omxlvz4 +_qw4c6lb] == _aq15o2n)*_qpwbsf3[i]  +  (_qsq0p82[j][16 -_omxlvz4 +_qw4c6lb] != _aq15o2n)*(*_oqgvzhq)[i][j];
                                    _d0p92k3          +=  _otn9uze[j];
                                }
                                _qpwbsf3[i] += 1u;
                                
                                for (j = 0u; j < _omxlvz4; j++)       {       _mgvu4is[ _qsq0p82[_8djzr5s][16 -_omxlvz4 +j] ] = 1u;       }
                                _qke2gjr = 0u;
                            }
                            _qw4c6lb += 1u;
                }   }   }
                
                _moxsfh4       = _z0yncxl(ALIGN64, 3*_qpwbsf3[i], sizeof(float));
                (*_yj9boq0)[i] = aligned_alloc(ALIGN64, _moxsfh4);      memcpy((*_yj9boq0)[i], _iq2fa39, 3*_qpwbsf3[i]* sizeof(float));
                (*_ljefy3z)[i] = aligned_alloc(ALIGN64, _moxsfh4);      memcpy((*_ljefy3z)[i], _b4ldims, 3*_qpwbsf3[i]* sizeof(float));
                

            }
            else
            {

            }
        }
        
        free(_2292mtr);   free(_3lbqgu3);   free(_5bwikot);   free(_iq2fa39);   free(_b4ldims);
        
        
        if (_vxpj4hn->_5izolnn > 0u)
        {   _moxsfh4   = _z0yncxl(ALIGN64, 3*_vxpj4hn->_5izolnn, sizeof(float));
            _m458ax4 = aligned_alloc(ALIGN64, _moxsfh4);
            _s6g298c = aligned_alloc(ALIGN64, _moxsfh4);
            _ukwqvr1 = aligned_alloc(ALIGN64, _moxsfh4);
            _moxsfh4   = _z0yncxl(ALIGN64, 2*_vxpj4hn->_3ru7myo, sizeof(float));
            _ffqot8h = aligned_alloc(ALIGN64, _moxsfh4);
            _zd76szn = aligned_alloc(ALIGN64, _moxsfh4);
            _p8ijl2z = aligned_alloc(ALIGN64, _moxsfh4);
            
            _moxsfh4   = _z0yncxl(ALIGN64, _1poxrte[_fsqg8x8], sizeof(unsigned short int *));
            *_1ajo091   = aligned_alloc(ALIGN64, _moxsfh4);
            _moxsfh4   = _z0yncxl(ALIGN64, _vxpj4hn->_5izolnn, sizeof(unsigned short int));
            for (i = 0u; i < _1poxrte[_fsqg8x8]; i++)
            {   (*_1ajo091)[i] = aligned_alloc(ALIGN64, _moxsfh4);             memset((*_1ajo091)[i], 0, _moxsfh4);
            }
            
            _moxsfh4   = _z0yncxl(ALIGN64, _1poxrte[_fsqg8x8], sizeof(unsigned short int *));
            *_3b9so0u   = aligned_alloc(ALIGN64, _moxsfh4);
            _moxsfh4   = _z0yncxl(ALIGN64, _vxpj4hn->_3ru7myo, sizeof(unsigned short int));
            for (i = 0u; i < _1poxrte[_fsqg8x8]; i++)
            {   (*_3b9so0u)[i] = aligned_alloc(ALIGN64, _moxsfh4);             memset((*_3b9so0u)[i], 0, _moxsfh4);
            }
            
            _moxsfh4   = _z0yncxl(ALIGN64, _1poxrte[_fsqg8x8], sizeof(float *));
            *_h61iaoi  = aligned_alloc(ALIGN64, _moxsfh4);
            *_plozjcn  = aligned_alloc(ALIGN64, _moxsfh4);
            *_qk4ltah  = aligned_alloc(ALIGN64, _moxsfh4);
            *_0ljch2v  = aligned_alloc(ALIGN64, _moxsfh4);
            *_0dtirc9  = aligned_alloc(ALIGN64, _moxsfh4);
            *_9dde3cu  = aligned_alloc(ALIGN64, _moxsfh4);
            
            for (i = 0u; i < _1poxrte[_fsqg8x8]; i++) 
            {   _8h0imjy = _fftouhy[i];
                
                _moxsfh4   = _z0yncxl(ALIGN64, 3*_vxpj4hn->_5izolnn, sizeof(float));
                memset(_m458ax4, 0, _moxsfh4);
                memset(_s6g298c, 0, _moxsfh4);
                memset(_ukwqvr1, 0, _moxsfh4);
                
                memcpy(_pt5r5rc,    &_rrpzlk3[_8h0imjy*16 +0], 3*sizeof(float));
                memcpy(_p9p9w4j, &_rrpzlk3[_8h0imjy*16 +6], 9*sizeof(float));
                
                memset(_otn9uze,   0, _vxpj4hn->_5izolnn*sizeof(unsigned int));
                memset(_mgvu4is, 0, _kb700uc*sizeof(unsigned int));
                _d0p92k3  = 0u;
                
                while (_d0p92k3 <  _vxpj4hn->_5izolnn) 
                {   
                    for (j = 0u; j <  _vxpj4hn->_5izolnn; j++)               {       if (_otn9uze[j] == 0u)    {      _8djzr5s = j;       break;          }           } 
                    
                    memcpy(_9p3jpxg, &_rrpzlk3[_8djzr5s*16 +0], 3*sizeof(float));
                    
                    _id4n1u3(_n42x2ez, _pt5r5rc, _9p3jpxg);               _31yrujo  = _mabsi9n(_n42x2ez);
                    
                    _qw5dsk7 = (_31yrujo     >= fMinDistFact*_s55j2o7->_gfu8iog)*1.0f  +  (_31yrujo      < fMinDistFact*_s55j2o7->_gfu8iog)*0.0f; 
                    _qw5dsk7 = (_rgb93wk[_8h0imjy*8 +3] == _rgb93wk[_8djzr5s*8 +3])*1.0f  +  (_rgb93wk[_8h0imjy*8 +3] != _rgb93wk[_8djzr5s*8 +3])*_qw5dsk7;
                    
                    _omxlvz4 = _qsq0p82[_8djzr5s][0]; 
                    
                    _qw4c6lb = 0;
                    _qke2gjr   = 1u;
                    while (_qke2gjr == 1u)
                    {   
                        _81bmqf7 = 0u;
                        while (_qw4c6lb < _omxlvz4-1)
                        {   _aq15o2n = _qsq0p82[_8djzr5s][16 -_omxlvz4 +_qw4c6lb];
                            
                            memcpy(_9p3jpxg, &_nf74t60[_aq15o2n][0], 3*sizeof(float));
                            _id4n1u3(_n42x2ez, _pt5r5rc, _9p3jpxg);               _31yrujo  = _mabsi9n(_n42x2ez);
                            
                            _81bmqf7  = (_mgvu4is[_aq15o2n]                                                  == 0u)*1u      +  (_mgvu4is[_aq15o2n]                                                  != 0u)*0u;
                            _81bmqf7  = (_cw2lbm2[_aq15o2n][0]                                               > 1u)*_81bmqf7  +  (_cw2lbm2[_aq15o2n][0]                                              <= 1u)*0u;
                            _81bmqf7  = (_qsq0p82[_8h0imjy][16 -_omxlvz4 +_qw4c6lb]                      != _aq15o2n)*_81bmqf7  +  (_qsq0p82[_8h0imjy][16 -_omxlvz4 +_qw4c6lb]                      == _aq15o2n)*0u;
                            _81bmqf7  = (_31yrujo   > (_nf74t60[_aq15o2n][5] +_s55j2o7->_a11dpg8*_s55j2o7->_gfu8iog))*_81bmqf7  +  (_31yrujo  <= (_nf74t60[_aq15o2n][5] +_s55j2o7->_a11dpg8*_s55j2o7->_gfu8iog))*0u;
                            
                            if (_81bmqf7 == 1u)         {       break;      }
                            
                            _qw4c6lb += 1u;
                        } 
                        
                        if (_81bmqf7 == 0u)
                        {   memcpy(_dan085r, &_3n3bzt9[_rgb93wk[_8djzr5s*8 +0]*4 +0], 3*sizeof(float));
                            memcpy(_ntdr5xy, &_3n3bzt9[_rgb93wk[_8djzr5s*8 +1]*4 +0], 3*sizeof(float));
                            memcpy(_fc4kl8c, &_3n3bzt9[_rgb93wk[_8djzr5s*8 +2]*4 +0], 3*sizeof(float));
                            
                            if (_s55j2o7->_s832ejx == 1u)
                            {   _oyr6ccv(_x4nxvkd, _vg5vmf6, _pt5r5rc[0], _pt5r5rc[1], _pt5r5rc[2], _dan085r, _ntdr5xy, _fc4kl8c, _s55j2o7->_cb1kjtg, 0.0f, 0.0f, _atj9tze->_kkrzd5e, _atj9tze->_6snglh6);
                                _oyr6ccv(_4xyjoa2, _vg5vmf6, _pt5r5rc[0], _pt5r5rc[1], _pt5r5rc[2], _dan085r, _ntdr5xy, _fc4kl8c, 0.0f, _s55j2o7->_cb1kjtg, 0.0f, _atj9tze->_kkrzd5e, _atj9tze->_6snglh6);
                                _oyr6ccv(_2k8xpnm, _vg5vmf6, _pt5r5rc[0], _pt5r5rc[1], _pt5r5rc[2], _dan085r, _ntdr5xy, _fc4kl8c, 0.0f, 0.0f, _s55j2o7->_cb1kjtg, _atj9tze->_kkrzd5e, _atj9tze->_6snglh6);
                            }
                            else
                            {   _25sq8cr(_x4nxvkd, _vg5vmf6, _pt5r5rc[0], _pt5r5rc[1], _pt5r5rc[2], _dan085r, _ntdr5xy, _fc4kl8c, _s55j2o7->_cb1kjtg, 0.0f, 0.0f, _atj9tze->_kkrzd5e, _atj9tze->_6snglh6);
                                _25sq8cr(_4xyjoa2, _vg5vmf6, _pt5r5rc[0], _pt5r5rc[1], _pt5r5rc[2], _dan085r, _ntdr5xy, _fc4kl8c, 0.0f, _s55j2o7->_cb1kjtg, 0.0f, _atj9tze->_kkrzd5e, _atj9tze->_6snglh6);
                                _25sq8cr(_2k8xpnm, _vg5vmf6, _pt5r5rc[0], _pt5r5rc[1], _pt5r5rc[2], _dan085r, _ntdr5xy, _fc4kl8c, 0.0f, 0.0f, _s55j2o7->_cb1kjtg, _atj9tze->_kkrzd5e, _atj9tze->_6snglh6);
                            }
                            
                            _yyv37di(_x4nxvkd, (_qw5dsk7/(_s55j2o7->_cb1kjtg*_rrpzlk3[_8djzr5s*16 +15])));      _ftu2gc5(&_kibu4fm[0], _p9p9w4j, _x4nxvkd);
                            _yyv37di(_4xyjoa2, (_qw5dsk7/(_s55j2o7->_cb1kjtg*_rrpzlk3[_8djzr5s*16 +15])));      _ftu2gc5(&_kibu4fm[3], _p9p9w4j, _4xyjoa2);
                            _yyv37di(_2k8xpnm, (_qw5dsk7/(_s55j2o7->_cb1kjtg*_rrpzlk3[_8djzr5s*16 +15])));      _ftu2gc5(&_kibu4fm[6], _p9p9w4j, _2k8xpnm);
                            
                            if (_8h0imjy == _8djzr5s)
                            {   
                                _7yt2m6t[i*16 +5] = _kibu4fm[0]*_rrpzlk3[_8djzr5s*16 +15];
                                _7yt2m6t[i*16 +6] = _kibu4fm[4]*_rrpzlk3[_8djzr5s*16 +15];
                                _7yt2m6t[i*16 +7] = _kibu4fm[8]*_rrpzlk3[_8djzr5s*16 +15];
                            }
                            
                            _m458ax4[_iymyrsv[i]*3 +0] = _kibu4fm[0];         _m458ax4[_iymyrsv[i]*3 +1] = _kibu4fm[3];         _m458ax4[_iymyrsv[i]*3 +2] = _kibu4fm[6];
                            _s6g298c[_iymyrsv[i]*3 +0] = _kibu4fm[1];         _s6g298c[_iymyrsv[i]*3 +1] = _kibu4fm[4];         _s6g298c[_iymyrsv[i]*3 +2] = _kibu4fm[7];
                            _ukwqvr1[_iymyrsv[i]*3 +0] = _kibu4fm[2];         _ukwqvr1[_iymyrsv[i]*3 +1] = _kibu4fm[5];         _ukwqvr1[_iymyrsv[i]*3 +2] = _kibu4fm[8];
                            (*_1ajo091)[i][_8djzr5s]       = _iymyrsv[i];
                            _iymyrsv[i]                  += 1u;
                            
                            _otn9uze[_8djzr5s]     = 1u;
                            
                            _d0p92k3 = 0u;
                            for (j = 0u; j <  _vxpj4hn->_5izolnn;  j++)     {       _d0p92k3 +=  _otn9uze[j];                                      }
                            
                            for (j = 0u; j < _omxlvz4; j++)           {       _mgvu4is[ _qsq0p82[_8djzr5s][16 -_omxlvz4 +j] ] = 1u;   }
                            _qke2gjr = 0u;
                        }
                        else
                        {   _aq15o2n = _qsq0p82[_8djzr5s][16 -_omxlvz4 +_qw4c6lb];
                            
                            memcpy(_dan085r, &_nf74t60[_aq15o2n][9],  3*sizeof(float));
                            memcpy(_ntdr5xy, &_nf74t60[_aq15o2n][12], 3*sizeof(float));
                            memcpy(_fc4kl8c, &_nf74t60[_aq15o2n][15], 3*sizeof(float));
                            memcpy(_gn0t6dv, &_nf74t60[_aq15o2n][18], 3*sizeof(float));
                            
                            if (_s55j2o7->_s832ejx == 1u)
                            {   _oyr6ccv(_c9x5jz6, _vg5vmf6, _pt5r5rc[0], _pt5r5rc[1], _pt5r5rc[2], _dan085r, _ntdr5xy, _gn0t6dv, _s55j2o7->_cb1kjtg, 0.0f, 0.0f, _atj9tze->_kkrzd5e, _atj9tze->_6snglh6);
                                _oyr6ccv(_1a6xvaj, _vg5vmf6, _pt5r5rc[0], _pt5r5rc[1], _pt5r5rc[2], _dan085r, _ntdr5xy, _gn0t6dv, 0.0f, _s55j2o7->_cb1kjtg, 0.0f, _atj9tze->_kkrzd5e, _atj9tze->_6snglh6);
                                _oyr6ccv(_jqy3ods, _vg5vmf6, _pt5r5rc[0], _pt5r5rc[1], _pt5r5rc[2], _dan085r, _ntdr5xy, _gn0t6dv, 0.0f, 0.0f, _s55j2o7->_cb1kjtg, _atj9tze->_kkrzd5e, _atj9tze->_6snglh6);
                                
                                _oyr6ccv(_x4nxvkd,  _vg5vmf6, _pt5r5rc[0], _pt5r5rc[1], _pt5r5rc[2], _ntdr5xy, _fc4kl8c, _gn0t6dv, _s55j2o7->_cb1kjtg, 0.0f, 0.0f, _atj9tze->_kkrzd5e, _atj9tze->_6snglh6);
                                _oyr6ccv(_4xyjoa2,  _vg5vmf6, _pt5r5rc[0], _pt5r5rc[1], _pt5r5rc[2], _ntdr5xy, _fc4kl8c, _gn0t6dv, 0.0f, _s55j2o7->_cb1kjtg, 0.0f, _atj9tze->_kkrzd5e, _atj9tze->_6snglh6);
                                _oyr6ccv(_2k8xpnm,  _vg5vmf6, _pt5r5rc[0], _pt5r5rc[1], _pt5r5rc[2], _ntdr5xy, _fc4kl8c, _gn0t6dv, 0.0f, 0.0f, _s55j2o7->_cb1kjtg, _atj9tze->_kkrzd5e, _atj9tze->_6snglh6);
                                
                                _x4nxvkd[0] += _c9x5jz6[0];        _x4nxvkd[1] += _c9x5jz6[1];        _x4nxvkd[2] += _c9x5jz6[2];
                                _x4nxvkd[3] += _c9x5jz6[3];        _x4nxvkd[4] += _c9x5jz6[4];        _x4nxvkd[5] += _c9x5jz6[5];
                                _4xyjoa2[0] += _1a6xvaj[0];        _4xyjoa2[1] += _1a6xvaj[1];        _4xyjoa2[2] += _1a6xvaj[2];
                                _4xyjoa2[3] += _1a6xvaj[3];        _4xyjoa2[4] += _1a6xvaj[4];        _4xyjoa2[5] += _1a6xvaj[5];
                                _2k8xpnm[0] += _jqy3ods[0];        _2k8xpnm[1] += _jqy3ods[1];        _2k8xpnm[2] += _jqy3ods[2];
                                _2k8xpnm[3] += _jqy3ods[3];        _2k8xpnm[4] += _jqy3ods[4];        _2k8xpnm[5] += _jqy3ods[5];
                            }
                            else
                            {   _25sq8cr(_c9x5jz6, _vg5vmf6, _pt5r5rc[0], _pt5r5rc[1], _pt5r5rc[2], _dan085r, _ntdr5xy, _gn0t6dv, _s55j2o7->_cb1kjtg, 0.0f, 0.0f, _atj9tze->_kkrzd5e, _atj9tze->_6snglh6);
                                _25sq8cr(_1a6xvaj, _vg5vmf6, _pt5r5rc[0], _pt5r5rc[1], _pt5r5rc[2], _dan085r, _ntdr5xy, _gn0t6dv, 0.0f, _s55j2o7->_cb1kjtg, 0.0f, _atj9tze->_kkrzd5e, _atj9tze->_6snglh6);
                                _25sq8cr(_jqy3ods, _vg5vmf6, _pt5r5rc[0], _pt5r5rc[1], _pt5r5rc[2], _dan085r, _ntdr5xy, _gn0t6dv, 0.0f, 0.0f, _s55j2o7->_cb1kjtg, _atj9tze->_kkrzd5e, _atj9tze->_6snglh6);
                                
                                _25sq8cr(_x4nxvkd,  _vg5vmf6, _pt5r5rc[0], _pt5r5rc[1], _pt5r5rc[2], _ntdr5xy, _fc4kl8c, _gn0t6dv, _s55j2o7->_cb1kjtg, 0.0f, 0.0f, _atj9tze->_kkrzd5e, _atj9tze->_6snglh6);
                                _25sq8cr(_4xyjoa2,  _vg5vmf6, _pt5r5rc[0], _pt5r5rc[1], _pt5r5rc[2], _ntdr5xy, _fc4kl8c, _gn0t6dv, 0.0f, _s55j2o7->_cb1kjtg, 0.0f, _atj9tze->_kkrzd5e, _atj9tze->_6snglh6);
                                _25sq8cr(_2k8xpnm,  _vg5vmf6, _pt5r5rc[0], _pt5r5rc[1], _pt5r5rc[2], _ntdr5xy, _fc4kl8c, _gn0t6dv, 0.0f, 0.0f, _s55j2o7->_cb1kjtg, _atj9tze->_kkrzd5e, _atj9tze->_6snglh6);
                                
                                _x4nxvkd[0] += _c9x5jz6[0];        _x4nxvkd[1] += _c9x5jz6[1];        _x4nxvkd[2] += _c9x5jz6[2];
                                _x4nxvkd[3] += _c9x5jz6[3];        _x4nxvkd[4] += _c9x5jz6[4];        _x4nxvkd[5] += _c9x5jz6[5];
                                _4xyjoa2[0] += _1a6xvaj[0];        _4xyjoa2[1] += _1a6xvaj[1];        _4xyjoa2[2] += _1a6xvaj[2];
                                _4xyjoa2[3] += _1a6xvaj[3];        _4xyjoa2[4] += _1a6xvaj[4];        _4xyjoa2[5] += _1a6xvaj[5];
                                _2k8xpnm[0] += _jqy3ods[0];        _2k8xpnm[1] += _jqy3ods[1];        _2k8xpnm[2] += _jqy3ods[2];
                                _2k8xpnm[3] += _jqy3ods[3];        _2k8xpnm[4] += _jqy3ods[4];        _2k8xpnm[5] += _jqy3ods[5];
                            }
                            
                            _yyv37di(_x4nxvkd, (_qw5dsk7/(_s55j2o7->_cb1kjtg*_nf74t60[_aq15o2n][3])));          _ftu2gc5(&_kibu4fm[0], _p9p9w4j, _x4nxvkd);
                            _yyv37di(_4xyjoa2, (_qw5dsk7/(_s55j2o7->_cb1kjtg*_nf74t60[_aq15o2n][3])));          _ftu2gc5(&_kibu4fm[3], _p9p9w4j, _4xyjoa2);
                            _yyv37di(_2k8xpnm, (_qw5dsk7/(_s55j2o7->_cb1kjtg*_nf74t60[_aq15o2n][3])));          _ftu2gc5(&_kibu4fm[6], _p9p9w4j, _2k8xpnm);
                            
                            _81bmqf7 = 1u;
                            for (j = 0; j < _cw2lbm2[_aq15o2n][0]; j++)
                            {   _4nnkm95 = _cw2lbm2[_aq15o2n][(j+1)];
                                
                                memcpy(_dan085r, &_nf74t60[_4nnkm95][9],  3*sizeof(float));
                                memcpy(_ntdr5xy, &_nf74t60[_4nnkm95][12], 3*sizeof(float));
                                memcpy(_fc4kl8c, &_nf74t60[_4nnkm95][15], 3*sizeof(float));
                                memcpy(_gn0t6dv, &_nf74t60[_4nnkm95][18], 3*sizeof(float));
                                
                                if (_s55j2o7->_s832ejx == 1u)
                                {   _oyr6ccv(_c9x5jz6, _vg5vmf6, _pt5r5rc[0], _pt5r5rc[1], _pt5r5rc[2], _dan085r, _ntdr5xy, _gn0t6dv, _s55j2o7->_cb1kjtg, 0.0f, 0.0f, _atj9tze->_kkrzd5e, _atj9tze->_6snglh6);
                                    _oyr6ccv(_1a6xvaj, _vg5vmf6, _pt5r5rc[0], _pt5r5rc[1], _pt5r5rc[2], _dan085r, _ntdr5xy, _gn0t6dv, 0.0f, _s55j2o7->_cb1kjtg, 0.0f, _atj9tze->_kkrzd5e, _atj9tze->_6snglh6);
                                    _oyr6ccv(_jqy3ods, _vg5vmf6, _pt5r5rc[0], _pt5r5rc[1], _pt5r5rc[2], _dan085r, _ntdr5xy, _gn0t6dv, 0.0f, 0.0f, _s55j2o7->_cb1kjtg, _atj9tze->_kkrzd5e, _atj9tze->_6snglh6);
                                    
                                    _oyr6ccv(_x4nxvkd,  _vg5vmf6, _pt5r5rc[0], _pt5r5rc[1], _pt5r5rc[2], _ntdr5xy, _fc4kl8c, _gn0t6dv, _s55j2o7->_cb1kjtg, 0.0f, 0.0f, _atj9tze->_kkrzd5e, _atj9tze->_6snglh6);
                                    _oyr6ccv(_4xyjoa2,  _vg5vmf6, _pt5r5rc[0], _pt5r5rc[1], _pt5r5rc[2], _ntdr5xy, _fc4kl8c, _gn0t6dv, 0.0f, _s55j2o7->_cb1kjtg, 0.0f, _atj9tze->_kkrzd5e, _atj9tze->_6snglh6);
                                    _oyr6ccv(_2k8xpnm,  _vg5vmf6, _pt5r5rc[0], _pt5r5rc[1], _pt5r5rc[2], _ntdr5xy, _fc4kl8c, _gn0t6dv, 0.0f, 0.0f, _s55j2o7->_cb1kjtg, _atj9tze->_kkrzd5e, _atj9tze->_6snglh6);
                                    
                                    _x4nxvkd[0] += _c9x5jz6[0];        _x4nxvkd[1] += _c9x5jz6[1];        _x4nxvkd[2] += _c9x5jz6[2];
                                    _x4nxvkd[3] += _c9x5jz6[3];        _x4nxvkd[4] += _c9x5jz6[4];        _x4nxvkd[5] += _c9x5jz6[5];
                                    _4xyjoa2[0] += _1a6xvaj[0];        _4xyjoa2[1] += _1a6xvaj[1];        _4xyjoa2[2] += _1a6xvaj[2];
                                    _4xyjoa2[3] += _1a6xvaj[3];        _4xyjoa2[4] += _1a6xvaj[4];        _4xyjoa2[5] += _1a6xvaj[5];
                                    _2k8xpnm[0] += _jqy3ods[0];        _2k8xpnm[1] += _jqy3ods[1];        _2k8xpnm[2] += _jqy3ods[2];
                                    _2k8xpnm[3] += _jqy3ods[3];        _2k8xpnm[4] += _jqy3ods[4];        _2k8xpnm[5] += _jqy3ods[5];
                                }
                                else
                                {   _25sq8cr(_c9x5jz6, _vg5vmf6, _pt5r5rc[0], _pt5r5rc[1], _pt5r5rc[2], _dan085r, _ntdr5xy, _gn0t6dv, _s55j2o7->_cb1kjtg, 0.0f, 0.0f, _atj9tze->_kkrzd5e, _atj9tze->_6snglh6);
                                    _25sq8cr(_1a6xvaj, _vg5vmf6, _pt5r5rc[0], _pt5r5rc[1], _pt5r5rc[2], _dan085r, _ntdr5xy, _gn0t6dv, 0.0f, _s55j2o7->_cb1kjtg, 0.0f, _atj9tze->_kkrzd5e, _atj9tze->_6snglh6);
                                    _25sq8cr(_jqy3ods, _vg5vmf6, _pt5r5rc[0], _pt5r5rc[1], _pt5r5rc[2], _dan085r, _ntdr5xy, _gn0t6dv, 0.0f, 0.0f, _s55j2o7->_cb1kjtg, _atj9tze->_kkrzd5e, _atj9tze->_6snglh6);
                                    
                                    _25sq8cr(_x4nxvkd,  _vg5vmf6, _pt5r5rc[0], _pt5r5rc[1], _pt5r5rc[2], _ntdr5xy, _fc4kl8c, _gn0t6dv, _s55j2o7->_cb1kjtg, 0.0f, 0.0f, _atj9tze->_kkrzd5e, _atj9tze->_6snglh6);
                                    _25sq8cr(_4xyjoa2,  _vg5vmf6, _pt5r5rc[0], _pt5r5rc[1], _pt5r5rc[2], _ntdr5xy, _fc4kl8c, _gn0t6dv, 0.0f, _s55j2o7->_cb1kjtg, 0.0f, _atj9tze->_kkrzd5e, _atj9tze->_6snglh6);
                                    _25sq8cr(_2k8xpnm,  _vg5vmf6, _pt5r5rc[0], _pt5r5rc[1], _pt5r5rc[2], _ntdr5xy, _fc4kl8c, _gn0t6dv, 0.0f, 0.0f, _s55j2o7->_cb1kjtg, _atj9tze->_kkrzd5e, _atj9tze->_6snglh6);
                                    
                                    _x4nxvkd[0] += _c9x5jz6[0];        _x4nxvkd[1] += _c9x5jz6[1];        _x4nxvkd[2] += _c9x5jz6[2];
                                    _x4nxvkd[3] += _c9x5jz6[3];        _x4nxvkd[4] += _c9x5jz6[4];        _x4nxvkd[5] += _c9x5jz6[5];
                                    _4xyjoa2[0] += _1a6xvaj[0];        _4xyjoa2[1] += _1a6xvaj[1];        _4xyjoa2[2] += _1a6xvaj[2];
                                    _4xyjoa2[3] += _1a6xvaj[3];        _4xyjoa2[4] += _1a6xvaj[4];        _4xyjoa2[5] += _1a6xvaj[5];
                                    _2k8xpnm[0] += _jqy3ods[0];        _2k8xpnm[1] += _jqy3ods[1];        _2k8xpnm[2] += _jqy3ods[2];
                                    _2k8xpnm[3] += _jqy3ods[3];        _2k8xpnm[4] += _jqy3ods[4];        _2k8xpnm[5] += _jqy3ods[5];
                                }
                                
                                _yyv37di(_x4nxvkd, (_qw5dsk7/(_s55j2o7->_cb1kjtg*_nf74t60[_4nnkm95][3])));      _ftu2gc5(&_zy1uu9y[0], _p9p9w4j, _x4nxvkd);
                                _yyv37di(_4xyjoa2, (_qw5dsk7/(_s55j2o7->_cb1kjtg*_nf74t60[_4nnkm95][3])));      _ftu2gc5(&_zy1uu9y[3], _p9p9w4j, _4xyjoa2);
                                _yyv37di(_2k8xpnm, (_qw5dsk7/(_s55j2o7->_cb1kjtg*_nf74t60[_4nnkm95][3])));      _ftu2gc5(&_zy1uu9y[6], _p9p9w4j, _2k8xpnm);
                                
                                memcpy(_9p3jpxg, &_nf74t60[_4nnkm95][0], 3*sizeof(float));
                                _id4n1u3(_n42x2ez, _pt5r5rc, _9p3jpxg);                              _31yrujo  = _mabsi9n(_n42x2ez);
                                _gymdm0l = (_31yrujo/_s55j2o7->_gfu8iog - _s55j2o7->_a11dpg8)/(_s55j2o7->_sbi4q2y - _s55j2o7->_a11dpg8);
                                _gymdm0l = (_gymdm0l >= 0.0f)*_gymdm0l  +  (_gymdm0l < 0.0f)*0.0f;
                                _gymdm0l = (_gymdm0l <= 1.0f)*_gymdm0l  +  (_gymdm0l > 1.0f)*1.0f;
                                _gymdm0l = _s55j2o7->_9v86wzs*_gymdm0l;
                                
                                _t26ys44 = sqrtf( (_zy1uu9y[0]-_kibu4fm[0])*(_zy1uu9y[0]-_kibu4fm[0]) + (_zy1uu9y[1]-_kibu4fm[1])*(_zy1uu9y[1]-_kibu4fm[1]) + (_zy1uu9y[2]-_kibu4fm[2])*(_zy1uu9y[2]-_kibu4fm[2])  + (_zy1uu9y[3]-_kibu4fm[3])*(_zy1uu9y[3]-_kibu4fm[3]) + (_zy1uu9y[4]-_kibu4fm[4])*(_zy1uu9y[4]-_kibu4fm[4]) + (_zy1uu9y[5]-_kibu4fm[5])*(_zy1uu9y[5]-_kibu4fm[5])   + (_zy1uu9y[6]-_kibu4fm[6])*(_zy1uu9y[6]-_kibu4fm[6]) + (_zy1uu9y[7]-_kibu4fm[7])*(_zy1uu9y[7]-_kibu4fm[7]) + (_zy1uu9y[8]-_kibu4fm[8])*(_zy1uu9y[8]-_kibu4fm[8]));
                                _5gdej31 = sqrtf(                _kibu4fm[0] *_kibu4fm[0]                 +                _kibu4fm[1] *_kibu4fm[1]                 +                _kibu4fm[2] *_kibu4fm[2]                  +                _kibu4fm[3] *_kibu4fm[3]                 +                _kibu4fm[4] *_kibu4fm[4]                 +                _kibu4fm[5] *_kibu4fm[5]                   +                _kibu4fm[6] *_kibu4fm[6]                 +                _kibu4fm[7] *_kibu4fm[7]                 +                _kibu4fm[8] *_kibu4fm[8]);
                                
                                _81bmqf7 = ((_t26ys44/_5gdej31) <= _gymdm0l)*_81bmqf7  +  ((_t26ys44/_5gdej31) > _gymdm0l)*0u; 
                            }
                            
                            if (_81bmqf7 == 1u)
                            {   _m458ax4[_iymyrsv[i]*3 +0] = _kibu4fm[0];         _m458ax4[_iymyrsv[i]*3 +1] = _kibu4fm[3];            _m458ax4[_iymyrsv[i]*3 +2] = _kibu4fm[6];
                                _s6g298c[_iymyrsv[i]*3 +0] = _kibu4fm[1];         _s6g298c[_iymyrsv[i]*3 +1] = _kibu4fm[4];            _s6g298c[_iymyrsv[i]*3 +2] = _kibu4fm[7];
                                _ukwqvr1[_iymyrsv[i]*3 +0] = _kibu4fm[2];         _ukwqvr1[_iymyrsv[i]*3 +1] = _kibu4fm[5];            _ukwqvr1[_iymyrsv[i]*3 +2] = _kibu4fm[8];
                                
                                _d0p92k3 = 0u;
                                for (j = 0u; j <  _vxpj4hn->_5izolnn;  j++)
                                {   _otn9uze[j]        = (_qsq0p82[j][16 -_omxlvz4 +_qw4c6lb] == _aq15o2n)*1u            +  (_qsq0p82[j][16 -_omxlvz4 +_qw4c6lb] != _aq15o2n)*_otn9uze[j];
                                    (*_1ajo091)[i][j] = (_qsq0p82[j][16 -_omxlvz4 +_qw4c6lb] == _aq15o2n)*_iymyrsv[i]  +  (_qsq0p82[j][16 -_omxlvz4 +_qw4c6lb] != _aq15o2n)*(*_1ajo091)[i][j];
                                    _d0p92k3          +=  _otn9uze[j];
                                }
                                _iymyrsv[i] += 1u;
                                
                                for (j = 0u; j < _omxlvz4; j++)       {       _mgvu4is[ _qsq0p82[_8djzr5s][16 -_omxlvz4 +j] ] = 1u;       }
                                _qke2gjr = 0u;
                            }
                            _qw4c6lb += 1u;
                }   }   }
                
                _moxsfh4       = _z0yncxl(ALIGN64, 3*_iymyrsv[i], sizeof(float));
                (*_h61iaoi)[i] = aligned_alloc(ALIGN64, _moxsfh4);        memcpy((*_h61iaoi)[i], _m458ax4, 3*_iymyrsv[i]* sizeof(float));
                (*_plozjcn)[i] = aligned_alloc(ALIGN64, _moxsfh4);        memcpy((*_plozjcn)[i], _s6g298c, 3*_iymyrsv[i]* sizeof(float));
                (*_qk4ltah)[i] = aligned_alloc(ALIGN64, _moxsfh4);        memcpy((*_qk4ltah)[i], _ukwqvr1, 3*_iymyrsv[i]* sizeof(float));
                

                
                
                _moxsfh4   = _z0yncxl(ALIGN64, 2*_vxpj4hn->_3ru7myo, sizeof(float));
                memset(_ffqot8h, 0, _moxsfh4);
                memset(_zd76szn, 0, _moxsfh4);
                memset(_p8ijl2z, 0, _moxsfh4);
                
                memset(_be3yqy5,   0, _vxpj4hn->_3ru7myo*sizeof(unsigned int));
                memset(_t4q9bc2, 0, _l9ue421*sizeof(unsigned int));
                _d0p92k3  = 0u;
                
                while (_d0p92k3 < _vxpj4hn->_3ru7myo) 
                {   
                    for (j = 0u; j < _vxpj4hn->_3ru7myo; j++)               {       if (_be3yqy5[j] == 0u)    {      _8djzr5s = j;       break;          }           } 
                    
                    memcpy(_9p3jpxg, &_4xkss0k[_8djzr5s*16 +0], 3*sizeof(float));
                    
                    _id4n1u3(_n42x2ez, _pt5r5rc, _9p3jpxg);               _31yrujo  = _mabsi9n(_n42x2ez);
                    
                    _qw5dsk7 = (_31yrujo >= fMinDistFact*_s55j2o7->_b68rerv)*1.0f  +  (_31yrujo < fMinDistFact*_s55j2o7->_b68rerv)*0.0f; 
                    
                    _omxlvz4 = _2h2p39y[_8djzr5s][0]; 
                    
                    _qw4c6lb = 0;
                    _qke2gjr   = 1u;
                    while (_qke2gjr == 1u)
                    {   
                        _81bmqf7 = 0u;
                        while (_qw4c6lb < _omxlvz4-1)
                        {   _aq15o2n = _2h2p39y[_8djzr5s][16 -_omxlvz4 +_qw4c6lb];
                            
                            memcpy(_9p3jpxg, &_e3nwxuj[_aq15o2n][0], 3*sizeof(float));
                            _id4n1u3(_n42x2ez, _pt5r5rc, _9p3jpxg);               _31yrujo  = _mabsi9n(_n42x2ez);
                            
                            _81bmqf7  = (_t4q9bc2[_aq15o2n]                                                 == 0u)*1u      +  (_t4q9bc2[_aq15o2n]                                                  != 0u)*0u;
                            _81bmqf7  = (_2qgv3dr[_aq15o2n][0]                                              > 1u)*_81bmqf7  +  (_2qgv3dr[_aq15o2n][0]                                              <= 1u)*0u;
                            _81bmqf7  = (_31yrujo  > (_e3nwxuj[_aq15o2n][5] +_s55j2o7->_a11dpg8*_s55j2o7->_b68rerv))*_81bmqf7  +  (_31yrujo  <= (_e3nwxuj[_aq15o2n][5] +_s55j2o7->_a11dpg8*_s55j2o7->_b68rerv))*0u;
                            
                            if (_81bmqf7 == 1u)         {       break;      }
                            
                            _qw4c6lb += 1u;
                        } 
                        
                        if (_81bmqf7 == 0u) 
                        {   memcpy(_dan085r, &_zzf8s2z[_vmefozt[_8djzr5s*8 +0]*4 +0], 3*sizeof(float));
                            memcpy(_ntdr5xy, &_zzf8s2z[_vmefozt[_8djzr5s*8 +1]*4 +0], 3*sizeof(float));
                            memcpy(_fc4kl8c, &_zzf8s2z[_vmefozt[_8djzr5s*8 +2]*4 +0], 3*sizeof(float));
                            
                            if (_s55j2o7->_s832ejx == 1u)
                            {   _oyr6ccv(_x4nxvkd, _vg5vmf6, _pt5r5rc[0], _pt5r5rc[1], _pt5r5rc[2], _dan085r, _ntdr5xy, _fc4kl8c, _s55j2o7->_ycz04cd, 0.0f, 0.0f, _atj9tze->_kkrzd5e, _atj9tze->_6snglh6);
                                _oyr6ccv(_4xyjoa2, _vg5vmf6, _pt5r5rc[0], _pt5r5rc[1], _pt5r5rc[2], _dan085r, _ntdr5xy, _fc4kl8c, 0.0f, _s55j2o7->_ycz04cd, 0.0f, _atj9tze->_kkrzd5e, _atj9tze->_6snglh6);
                           }
                            else
                            {   _25sq8cr(_x4nxvkd, _vg5vmf6, _pt5r5rc[0], _pt5r5rc[1], _pt5r5rc[2], _dan085r, _ntdr5xy, _fc4kl8c, _s55j2o7->_ycz04cd, 0.0f, 0.0f, _atj9tze->_kkrzd5e, _atj9tze->_6snglh6);
                                _25sq8cr(_4xyjoa2, _vg5vmf6, _pt5r5rc[0], _pt5r5rc[1], _pt5r5rc[2], _dan085r, _ntdr5xy, _fc4kl8c, 0.0f, _s55j2o7->_ycz04cd, 0.0f, _atj9tze->_kkrzd5e, _atj9tze->_6snglh6);
                            }
                            
                            _yyv37di(_x4nxvkd, (_qw5dsk7/(_s55j2o7->_ycz04cd*_4xkss0k[_8djzr5s*16 +15])));      _ftu2gc5(&_kibu4fm[0], _p9p9w4j, _x4nxvkd);
                            _yyv37di(_4xyjoa2, (_qw5dsk7/(_s55j2o7->_ycz04cd*_4xkss0k[_8djzr5s*16 +15])));      _ftu2gc5(&_kibu4fm[3], _p9p9w4j, _4xyjoa2);
                            
                            _ffqot8h[_iob0nlw[i]*2 +0] = _kibu4fm[0];         _ffqot8h[_iob0nlw[i]*2 +1] = _kibu4fm[3];
                            _zd76szn[_iob0nlw[i]*2 +0] = _kibu4fm[1];         _zd76szn[_iob0nlw[i]*2 +1] = _kibu4fm[4];
                            _p8ijl2z[_iob0nlw[i]*2 +0] = _kibu4fm[2];         _p8ijl2z[_iob0nlw[i]*2 +1] = _kibu4fm[5];
                            (*_3b9so0u)[i][_8djzr5s]       = _iob0nlw[i];
                            _iob0nlw[i]                  += 1u;
                            
                            _be3yqy5[_8djzr5s]     = 1u;
                            
                            _d0p92k3 = 0u;
                            for (j = 0u; j < _vxpj4hn->_3ru7myo;  j++)      {       _d0p92k3 +=  _be3yqy5[j];                                      }
                            
                            for (j = 0u; j < _omxlvz4; j++)                {       _t4q9bc2[ _2h2p39y[_8djzr5s][16 -_omxlvz4 +j] ] = 1u;        }
                            _qke2gjr = 0u;
                        }
                        else
                        {   _aq15o2n = _2h2p39y[_8djzr5s][16 -_omxlvz4 +_qw4c6lb];
                            
                            memcpy(_dan085r, &_e3nwxuj[_aq15o2n][9],  3*sizeof(float));
                            memcpy(_ntdr5xy, &_e3nwxuj[_aq15o2n][12], 3*sizeof(float));
                            memcpy(_fc4kl8c, &_e3nwxuj[_aq15o2n][15], 3*sizeof(float));
                            memcpy(_gn0t6dv, &_e3nwxuj[_aq15o2n][18], 3*sizeof(float));
                            
                            if (_s55j2o7->_s832ejx == 1u)
                            {   _oyr6ccv(_c9x5jz6, _vg5vmf6, _pt5r5rc[0], _pt5r5rc[1], _pt5r5rc[2], _dan085r, _ntdr5xy, _gn0t6dv, _s55j2o7->_ycz04cd, 0.0f, 0.0f, _atj9tze->_kkrzd5e, _atj9tze->_6snglh6);
                                _oyr6ccv(_1a6xvaj, _vg5vmf6, _pt5r5rc[0], _pt5r5rc[1], _pt5r5rc[2], _dan085r, _ntdr5xy, _gn0t6dv, 0.0f, _s55j2o7->_ycz04cd, 0.0f, _atj9tze->_kkrzd5e, _atj9tze->_6snglh6);
                                
                                _oyr6ccv(_x4nxvkd,  _vg5vmf6, _pt5r5rc[0], _pt5r5rc[1], _pt5r5rc[2], _ntdr5xy, _fc4kl8c, _gn0t6dv, _s55j2o7->_ycz04cd, 0.0f, 0.0f, _atj9tze->_kkrzd5e, _atj9tze->_6snglh6);
                                _oyr6ccv(_4xyjoa2,  _vg5vmf6, _pt5r5rc[0], _pt5r5rc[1], _pt5r5rc[2], _ntdr5xy, _fc4kl8c, _gn0t6dv, 0.0f, _s55j2o7->_ycz04cd, 0.0f, _atj9tze->_kkrzd5e, _atj9tze->_6snglh6);
                                
                                _x4nxvkd[0] += _c9x5jz6[0];        _x4nxvkd[1] += _c9x5jz6[1];        _x4nxvkd[2] += _c9x5jz6[2];
                                _x4nxvkd[3] += _c9x5jz6[3];        _x4nxvkd[4] += _c9x5jz6[4];        _x4nxvkd[5] += _c9x5jz6[5];
                                _4xyjoa2[0] += _1a6xvaj[0];        _4xyjoa2[1] += _1a6xvaj[1];        _4xyjoa2[2] += _1a6xvaj[2];
                                _4xyjoa2[3] += _1a6xvaj[3];        _4xyjoa2[4] += _1a6xvaj[4];        _4xyjoa2[5] += _1a6xvaj[5];
                            }
                            else
                            {   _25sq8cr(_c9x5jz6, _vg5vmf6, _pt5r5rc[0], _pt5r5rc[1], _pt5r5rc[2], _dan085r, _ntdr5xy, _gn0t6dv, _s55j2o7->_ycz04cd, 0.0f, 0.0f, _atj9tze->_kkrzd5e, _atj9tze->_6snglh6);
                                _25sq8cr(_1a6xvaj, _vg5vmf6, _pt5r5rc[0], _pt5r5rc[1], _pt5r5rc[2], _dan085r, _ntdr5xy, _gn0t6dv, 0.0f, _s55j2o7->_ycz04cd, 0.0f, _atj9tze->_kkrzd5e, _atj9tze->_6snglh6);
                                
                                _25sq8cr(_x4nxvkd,  _vg5vmf6, _pt5r5rc[0], _pt5r5rc[1], _pt5r5rc[2], _ntdr5xy, _fc4kl8c, _gn0t6dv, _s55j2o7->_ycz04cd, 0.0f, 0.0f, _atj9tze->_kkrzd5e, _atj9tze->_6snglh6);
                                _25sq8cr(_4xyjoa2,  _vg5vmf6, _pt5r5rc[0], _pt5r5rc[1], _pt5r5rc[2], _ntdr5xy, _fc4kl8c, _gn0t6dv, 0.0f, _s55j2o7->_ycz04cd, 0.0f, _atj9tze->_kkrzd5e, _atj9tze->_6snglh6);
                                
                                _x4nxvkd[0] += _c9x5jz6[0];        _x4nxvkd[1] += _c9x5jz6[1];        _x4nxvkd[2] += _c9x5jz6[2];
                                _x4nxvkd[3] += _c9x5jz6[3];        _x4nxvkd[4] += _c9x5jz6[4];        _x4nxvkd[5] += _c9x5jz6[5];
                                _4xyjoa2[0] += _1a6xvaj[0];        _4xyjoa2[1] += _1a6xvaj[1];        _4xyjoa2[2] += _1a6xvaj[2];
                                _4xyjoa2[3] += _1a6xvaj[3];        _4xyjoa2[4] += _1a6xvaj[4];        _4xyjoa2[5] += _1a6xvaj[5];
                            }
                            
                            _yyv37di(_x4nxvkd, (_qw5dsk7/(_s55j2o7->_ycz04cd*_e3nwxuj[_aq15o2n][3])));          _ftu2gc5(&_kibu4fm[0], _p9p9w4j, _x4nxvkd);
                            _yyv37di(_4xyjoa2, (_qw5dsk7/(_s55j2o7->_ycz04cd*_e3nwxuj[_aq15o2n][3])));          _ftu2gc5(&_kibu4fm[3], _p9p9w4j, _4xyjoa2);
                            
                            _81bmqf7 = 1u;
                            for (j = 0; j < _2qgv3dr[_aq15o2n][0]; j++)
                            {   _4nnkm95 = _2qgv3dr[_aq15o2n][(j+1)];
                                
                                memcpy(_dan085r, &_e3nwxuj[_4nnkm95][9],  3*sizeof(float));
                                memcpy(_ntdr5xy, &_e3nwxuj[_4nnkm95][12], 3*sizeof(float));
                                memcpy(_fc4kl8c, &_e3nwxuj[_4nnkm95][15], 3*sizeof(float));
                                memcpy(_gn0t6dv, &_e3nwxuj[_4nnkm95][18], 3*sizeof(float));
                                
                                if (_s55j2o7->_s832ejx == 1u)
                                {   _oyr6ccv(_c9x5jz6, _vg5vmf6, _pt5r5rc[0], _pt5r5rc[1], _pt5r5rc[2], _dan085r, _ntdr5xy, _gn0t6dv, _s55j2o7->_ycz04cd, 0.0f, 0.0f, _atj9tze->_kkrzd5e, _atj9tze->_6snglh6);
                                    _oyr6ccv(_1a6xvaj, _vg5vmf6, _pt5r5rc[0], _pt5r5rc[1], _pt5r5rc[2], _dan085r, _ntdr5xy, _gn0t6dv, 0.0f, _s55j2o7->_ycz04cd, 0.0f, _atj9tze->_kkrzd5e, _atj9tze->_6snglh6);
                                    
                                    _oyr6ccv(_x4nxvkd,  _vg5vmf6, _pt5r5rc[0], _pt5r5rc[1], _pt5r5rc[2], _ntdr5xy, _fc4kl8c, _gn0t6dv, _s55j2o7->_ycz04cd, 0.0f, 0.0f, _atj9tze->_kkrzd5e, _atj9tze->_6snglh6);
                                    _oyr6ccv(_4xyjoa2,  _vg5vmf6, _pt5r5rc[0], _pt5r5rc[1], _pt5r5rc[2], _ntdr5xy, _fc4kl8c, _gn0t6dv, 0.0f, _s55j2o7->_ycz04cd, 0.0f, _atj9tze->_kkrzd5e, _atj9tze->_6snglh6);
                                    
                                    _x4nxvkd[0] += _c9x5jz6[0];        _x4nxvkd[1] += _c9x5jz6[1];        _x4nxvkd[2] += _c9x5jz6[2];
                                    _x4nxvkd[3] += _c9x5jz6[3];        _x4nxvkd[4] += _c9x5jz6[4];        _x4nxvkd[5] += _c9x5jz6[5];
                                    _4xyjoa2[0] += _1a6xvaj[0];        _4xyjoa2[1] += _1a6xvaj[1];        _4xyjoa2[2] += _1a6xvaj[2];
                                    _4xyjoa2[3] += _1a6xvaj[3];        _4xyjoa2[4] += _1a6xvaj[4];        _4xyjoa2[5] += _1a6xvaj[5];
                                }
                                else
                                {   _25sq8cr(_c9x5jz6, _vg5vmf6, _pt5r5rc[0], _pt5r5rc[1], _pt5r5rc[2], _dan085r, _ntdr5xy, _gn0t6dv, _s55j2o7->_ycz04cd, 0.0f, 0.0f, _atj9tze->_kkrzd5e, _atj9tze->_6snglh6);
                                    _25sq8cr(_1a6xvaj, _vg5vmf6, _pt5r5rc[0], _pt5r5rc[1], _pt5r5rc[2], _dan085r, _ntdr5xy, _gn0t6dv, 0.0f, _s55j2o7->_ycz04cd, 0.0f, _atj9tze->_kkrzd5e, _atj9tze->_6snglh6);
                                    
                                    _25sq8cr(_x4nxvkd,  _vg5vmf6, _pt5r5rc[0], _pt5r5rc[1], _pt5r5rc[2], _ntdr5xy, _fc4kl8c, _gn0t6dv, _s55j2o7->_ycz04cd, 0.0f, 0.0f, _atj9tze->_kkrzd5e, _atj9tze->_6snglh6);
                                    _25sq8cr(_4xyjoa2,  _vg5vmf6, _pt5r5rc[0], _pt5r5rc[1], _pt5r5rc[2], _ntdr5xy, _fc4kl8c, _gn0t6dv, 0.0f, _s55j2o7->_ycz04cd, 0.0f, _atj9tze->_kkrzd5e, _atj9tze->_6snglh6);
                                    
                                    _x4nxvkd[0] += _c9x5jz6[0];        _x4nxvkd[1] += _c9x5jz6[1];        _x4nxvkd[2] += _c9x5jz6[2];
                                    _x4nxvkd[3] += _c9x5jz6[3];        _x4nxvkd[4] += _c9x5jz6[4];        _x4nxvkd[5] += _c9x5jz6[5];
                                    _4xyjoa2[0] += _1a6xvaj[0];        _4xyjoa2[1] += _1a6xvaj[1];        _4xyjoa2[2] += _1a6xvaj[2];
                                    _4xyjoa2[3] += _1a6xvaj[3];        _4xyjoa2[4] += _1a6xvaj[4];        _4xyjoa2[5] += _1a6xvaj[5];
                                }
                                
                                _yyv37di(_x4nxvkd, (_qw5dsk7/(_s55j2o7->_ycz04cd*_e3nwxuj[_4nnkm95][3])));      _ftu2gc5(&_zy1uu9y[0], _p9p9w4j, _x4nxvkd);
                                _yyv37di(_4xyjoa2, (_qw5dsk7/(_s55j2o7->_ycz04cd*_e3nwxuj[_4nnkm95][3])));      _ftu2gc5(&_zy1uu9y[3], _p9p9w4j, _4xyjoa2);
                                
                                memcpy(_9p3jpxg, &_e3nwxuj[_4nnkm95][0], 3*sizeof(float));
                                _id4n1u3(_n42x2ez, _pt5r5rc, _9p3jpxg);                              _31yrujo  = _mabsi9n(_n42x2ez);
                                _gymdm0l = (_31yrujo/_s55j2o7->_b68rerv - _s55j2o7->_a11dpg8)/(_s55j2o7->_sbi4q2y - _s55j2o7->_a11dpg8);
                                _gymdm0l = (_gymdm0l >= 0.0f)*_gymdm0l  +  (_gymdm0l < 0.0f)*0.0f;
                                _gymdm0l = (_gymdm0l <= 1.0f)*_gymdm0l  +  (_gymdm0l > 1.0f)*1.0f;
                                _gymdm0l = _s55j2o7->_9v86wzs*_gymdm0l;
                                
                                _t26ys44 = sqrtf( (_zy1uu9y[0]-_kibu4fm[0])*(_zy1uu9y[0]-_kibu4fm[0]) + (_zy1uu9y[1]-_kibu4fm[1])*(_zy1uu9y[1]-_kibu4fm[1]) + (_zy1uu9y[2]-_kibu4fm[2])*(_zy1uu9y[2]-_kibu4fm[2])  + (_zy1uu9y[3]-_kibu4fm[3])*(_zy1uu9y[3]-_kibu4fm[3]) + (_zy1uu9y[4]-_kibu4fm[4])*(_zy1uu9y[4]-_kibu4fm[4]) + (_zy1uu9y[5]-_kibu4fm[5])*(_zy1uu9y[5]-_kibu4fm[5]) );
                                _5gdej31 = sqrtf(                _kibu4fm[0] *_kibu4fm[0]                 +                _kibu4fm[1] *_kibu4fm[1]                 +                _kibu4fm[2] *_kibu4fm[2]                  +                _kibu4fm[3] *_kibu4fm[3]                 +                _kibu4fm[4] *_kibu4fm[4]                 +                _kibu4fm[5] *_kibu4fm[5] );
                                
                                _81bmqf7 = ((_t26ys44/_5gdej31) <= _gymdm0l)*_81bmqf7  +  ((_t26ys44/_5gdej31) > _gymdm0l)*0u; 
                            }
                            
                            if (_81bmqf7 == 1u)
                            {   _ffqot8h[_iob0nlw[i]*2 +0] = _kibu4fm[0];         _ffqot8h[_iob0nlw[i]*2 +1] = _kibu4fm[3];
                                _zd76szn[_iob0nlw[i]*2 +0] = _kibu4fm[1];         _zd76szn[_iob0nlw[i]*2 +1] = _kibu4fm[4];
                                _p8ijl2z[_iob0nlw[i]*2 +0] = _kibu4fm[2];         _p8ijl2z[_iob0nlw[i]*2 +1] = _kibu4fm[5];
                                
                                _d0p92k3 = 0u;
                                for (j = 0u; j < _vxpj4hn->_3ru7myo;  j++)
                                {   _be3yqy5[j]        = (_2h2p39y[j][16 -_omxlvz4 +_qw4c6lb] == _aq15o2n)*1u            +  (_2h2p39y[j][16 -_omxlvz4 +_qw4c6lb] != _aq15o2n)*_be3yqy5[j];
                                    (*_3b9so0u)[i][j] = (_2h2p39y[j][16 -_omxlvz4 +_qw4c6lb] == _aq15o2n)*_iob0nlw[i]  +  (_2h2p39y[j][16 -_omxlvz4 +_qw4c6lb] != _aq15o2n)*(*_3b9so0u)[i][j];
                                    _d0p92k3          +=  _be3yqy5[j];
                                }
                                _iob0nlw[i] += 1u;
                                
                                for (j = 0u; j < _omxlvz4; j++)       {       _t4q9bc2[ _2h2p39y[_8djzr5s][16 -_omxlvz4 +j] ] = 1u;       }
                                _qke2gjr = 0u;
                            }
                            _qw4c6lb += 1u;
                }   }   }
                
                _moxsfh4       = _z0yncxl(ALIGN64, 2*_iob0nlw[i], sizeof(float));
                (*_0ljch2v)[i] = aligned_alloc(ALIGN64, _moxsfh4);        memcpy((*_0ljch2v)[i], _ffqot8h, 2*_iob0nlw[i]* sizeof(float));
                (*_0dtirc9)[i] = aligned_alloc(ALIGN64, _moxsfh4);        memcpy((*_0dtirc9)[i], _zd76szn, 2*_iob0nlw[i]* sizeof(float));
                (*_9dde3cu)[i] = aligned_alloc(ALIGN64, _moxsfh4);        memcpy((*_9dde3cu)[i], _p8ijl2z, 2*_iob0nlw[i]* sizeof(float));
                

            }
        }
        
        free(_be3yqy5);   free(_t4q9bc2);   free(_otn9uze);   free(_mgvu4is);
        free(_m458ax4);   free(_s6g298c);   free(_ukwqvr1);  free(_ffqot8h);  free(_zd76szn);  free(_p8ijl2z);
        
        
    } 
    
    
    MPI_Barrier( MPI_COMM_WORLD );
    clock_gettime(CLOCK_REALTIME, &_1zocfi4);
    _kicwmpe = (_1zocfi4.tv_sec - _glmy62x.tv_sec) + (_1zocfi4.tv_nsec - _glmy62x.tv_nsec)/1.0E+9;
    if (_fsqg8x8 == 0)         {   fprintf(stdout,"\nTime taken for Kh_matrix(in seconds): %lf\n", _kicwmpe);             }
    
    for (int i = 0; i < _l9ue421; i++)          {       free(_e3nwxuj[i]);            }
    for (int i = 0; i < _vxpj4hn->_3ru7myo; i++)       {       free(_2h2p39y[i]);             }
    for (int i = 0; i < _vxpj4hn->_3ru7myo*4; i++)     {       free(_2qgv3dr[i]);            } 
    for (int i = 0; i < _kb700uc; i++)          {       free(_nf74t60[i]);            }
    for (int i = 0; i < _vxpj4hn->_5izolnn; i++)       {       free(_qsq0p82[i]);             }
    for (int i = 0; i < _vxpj4hn->_5izolnn*4; i++)     {       free(_cw2lbm2[i]);            } 
    free(_e3nwxuj);   free(_nf74t60);   free(_2h2p39y);   free(_qsq0p82);   free(_2qgv3dr);   free(_cw2lbm2);
    
    
    unsigned int i, _1b07g2r[4],   _l2ld05n[4],   _ybxkpm2[4],   _jg2y8rg[4];
    memset(_1b07g2r,     0, 4*sizeof(unsigned int));
    memset(_l2ld05n,     0, 4*sizeof(unsigned int));
    memset(_ybxkpm2,     0, 4*sizeof(unsigned int));
    memset(_jg2y8rg, 0, 4*sizeof(unsigned int));
    
    for (i = 0u; i < _lexdw6l[_fsqg8x8]; i++)
    {   _1b07g2r[0] += _qnsn2s8[i];             _jg2y8rg[0] = (_jg2y8rg[0] >= _qnsn2s8[i])*_jg2y8rg[0]  +  (_jg2y8rg[0] < _qnsn2s8[i])*_qnsn2s8[i];
    }
    if (_vxpj4hn->_5izolnn > 0u)
    {   for (i = 0u; i < _lexdw6l[_fsqg8x8]; i++)
        {   _1b07g2r[1] += _qpwbsf3[i];         _jg2y8rg[1] = (_jg2y8rg[1] >= _qpwbsf3[i])*_jg2y8rg[1]  +  (_jg2y8rg[1] < _qpwbsf3[i])*_qpwbsf3[i];
        }
    }
    for (i = 0u; i < _1poxrte[_fsqg8x8]; i++)      
    {   _1b07g2r[2] += _iymyrsv[i];             _jg2y8rg[2] = (_jg2y8rg[2] >= _iymyrsv[i])*_jg2y8rg[2]  +  (_jg2y8rg[2] < _iymyrsv[i])*_iymyrsv[i];
        _1b07g2r[3] += _iob0nlw[i];             _jg2y8rg[3] = (_jg2y8rg[3] >= _iob0nlw[i])*_jg2y8rg[3]  +  (_jg2y8rg[3] < _iob0nlw[i])*_iob0nlw[i];
    }
    
    MPI_Allreduce(MPI_IN_PLACE, _jg2y8rg, 4, MPI_UNSIGNED, MPI_MAX, MPI_COMM_WORLD);
    _vxpj4hn->_dn1th1x = MAX(_jg2y8rg[0], _jg2y8rg[3]);
    _vxpj4hn->_w2maiet = MAX(_jg2y8rg[1], _jg2y8rg[2]);
    
    _1b07g2r[0] =                             _1b07g2r[0]/_lexdw6l[_fsqg8x8];
    _1b07g2r[1] =                             _1b07g2r[1]/_lexdw6l[_fsqg8x8];
    _1b07g2r[2] = (_1poxrte[_fsqg8x8] > 0u)  ?  (_1b07g2r[2]/_1poxrte[_fsqg8x8])  :  0u;
    _1b07g2r[3] = (_1poxrte[_fsqg8x8] > 0u)  ?  (_1b07g2r[3]/_1poxrte[_fsqg8x8])  :  0u;
    
    memcpy(_l2ld05n, _1b07g2r,  4*sizeof(unsigned int));
    memcpy(_ybxkpm2, _1b07g2r,  4*sizeof(unsigned int));
    MPI_Allreduce(MPI_IN_PLACE, _l2ld05n, 4, MPI_UNSIGNED, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, _ybxkpm2, 4, MPI_UNSIGNED, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, _1b07g2r, 4, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
    
    _1b07g2r[0] /= _0rmhrcn;        _1b07g2r[1] /= _0rmhrcn;        _1b07g2r[2] /= _0rmhrcn;        _1b07g2r[3] /= _0rmhrcn;
    
    if (_fsqg8x8 == 0)
    {   fprintf(stdout,"Kh_lengths (average and range) across RANKs\n");



        fprintf(stdout,"Minimum:         %u\n", _l2ld05n[0]);
        fprintf(stdout,"Average:         %u\n", _1b07g2r[0]);
        fprintf(stdout,"Maximum:         %u\n", _ybxkpm2[0]);
    }
    
    
    MPI_Barrier( MPI_COMM_WORLD );
    
    return;
}



void   _82n2atq(int _fsqg8x8, const char *restrict _g79or78, const struct _xjuy0bv *restrict _vxpj4hn, const struct _j49itn7 *restrict _s55j2o7, const int *restrict _paih0p4, const int *restrict _20n8tte, const int *restrict _lexdw6l, const int *restrict _1poxrte, const unsigned int *restrict _1nlc8jp, const unsigned int *restrict _fftouhy, const float *restrict _6fokzix, const float *restrict _7yt2m6t, unsigned short int *restrict _qnsn2s8, unsigned short int *restrict _qpwbsf3, unsigned short int *restrict _iymyrsv, unsigned short int *restrict _iob0nlw, unsigned short int **restrict _5gzpl1q, unsigned short int **restrict _oqgvzhq, unsigned short int **restrict _1ajo091, unsigned short int **restrict _3b9so0u, float **restrict _ckeh5sw, float **restrict _3ac6pe6, float **restrict _obyghgk, float **restrict _yj9boq0, float **restrict _ljefy3z, float **restrict _h61iaoi, float **restrict _plozjcn, float **restrict _qk4ltah, float **restrict _0ljch2v, float **restrict _0dtirc9, float **restrict _9dde3cu)
{   
    MPI_Status  _0get8t6;
    MPI_Offset  _9kwv08h;
    MPI_File  _6zx1tv1;
    
    unsigned int i;
    unsigned long long _9wszv83;
    unsigned short int *_dg1wjui = NULL,   *_op767j2 = NULL,   *_9f77lu8 = NULL,   *_el5015o = NULL;
    unsigned long long *_7pchcbn = NULL,   *_zzf8adf = NULL,   *_lmz04x0 = NULL,   *_qzfhomt = NULL;
    float *_0wfoy4d = NULL,   *_jssuuqo = NULL;
    
    _dg1wjui      = malloc(_vxpj4hn->_3ru7myo *sizeof(unsigned short int));        memset(_dg1wjui, 0, _vxpj4hn->_3ru7myo *sizeof(unsigned short int));
    _7pchcbn      = malloc(_vxpj4hn->_3ru7myo *sizeof(unsigned long long));        memset(_7pchcbn, 0, _vxpj4hn->_3ru7myo *sizeof(unsigned long long));
    _0wfoy4d     = malloc((3*_vxpj4hn->_3ru7myo) *sizeof(float));                 memset(_0wfoy4d,0,(3*_vxpj4hn->_3ru7myo) *sizeof(float));
    
    if (_vxpj4hn->_5izolnn > 0u)
    {   _op767j2  = malloc(_vxpj4hn->_3ru7myo *sizeof(unsigned short int));        memset(_op767j2, 0, _vxpj4hn->_3ru7myo *sizeof(unsigned short int));
        _zzf8adf  = malloc(_vxpj4hn->_3ru7myo *sizeof(unsigned long long));        memset(_zzf8adf, 0, _vxpj4hn->_3ru7myo *sizeof(unsigned long long));
        _9f77lu8  = malloc(_vxpj4hn->_5izolnn *sizeof(unsigned short int));        memset(_9f77lu8, 0, _vxpj4hn->_5izolnn *sizeof(unsigned short int));
        _lmz04x0  = malloc(_vxpj4hn->_5izolnn *sizeof(unsigned long long));        memset(_lmz04x0, 0, _vxpj4hn->_5izolnn *sizeof(unsigned long long));
        _el5015o  = malloc(_vxpj4hn->_5izolnn *sizeof(unsigned short int));        memset(_el5015o, 0, _vxpj4hn->_5izolnn *sizeof(unsigned short int));
        _qzfhomt  = malloc(_vxpj4hn->_5izolnn *sizeof(unsigned long long));        memset(_qzfhomt, 0, _vxpj4hn->_5izolnn *sizeof(unsigned long long));
        _jssuuqo = malloc((3*_vxpj4hn->_5izolnn) *sizeof(float));                 memset(_jssuuqo,0,(3*_vxpj4hn->_5izolnn) *sizeof(float));
    }
    
    for (i = 0u; i < _lexdw6l[_fsqg8x8]; i++)          {   _dg1wjui[_1nlc8jp[i]] = _qnsn2s8[i];       }
    MPI_Allreduce(MPI_IN_PLACE, _dg1wjui, _vxpj4hn->_3ru7myo,  MPI_UNSIGNED_SHORT, MPI_SUM, MPI_COMM_WORLD);
    
    for (i = 1u; i < _vxpj4hn->_3ru7myo; i++)
    {   _7pchcbn[i] = _7pchcbn[i-1] + (_vxpj4hn->_3ru7myo*sizeof(unsigned short int) + _dg1wjui[i-1]*6*sizeof(float));
    }
    
    for (i = 0u; i < _lexdw6l[_fsqg8x8]; i++)          {   _0wfoy4d[_1nlc8jp[i] +0*_vxpj4hn->_3ru7myo] = _6fokzix[i*16 +5];      _0wfoy4d[_1nlc8jp[i] +1*_vxpj4hn->_3ru7myo] = _6fokzix[i*16 +6];          _0wfoy4d[_1nlc8jp[i] +2*_vxpj4hn->_3ru7myo] = _6fokzix[i*16 +7];            }
    MPI_Allreduce(MPI_IN_PLACE, _0wfoy4d, (3*_vxpj4hn->_3ru7myo), MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    
    if (_vxpj4hn->_5izolnn > 0u)
    {   for (i = 0u; i < _lexdw6l[_fsqg8x8]; i++)      {   _op767j2[_1nlc8jp[i]] = _qpwbsf3[i];       }
        for (i = 0u; i < _1poxrte[_fsqg8x8]; i++)      {   _9f77lu8[_fftouhy[i]] = _iymyrsv[i];       _el5015o[_fftouhy[i]] = _iob0nlw[i];          }
        MPI_Allreduce(MPI_IN_PLACE, _op767j2, _vxpj4hn->_3ru7myo,  MPI_UNSIGNED_SHORT, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, _9f77lu8, _vxpj4hn->_5izolnn,  MPI_UNSIGNED_SHORT, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, _el5015o, _vxpj4hn->_5izolnn,  MPI_UNSIGNED_SHORT, MPI_SUM, MPI_COMM_WORLD);
        
        for (i = 1u; i < _vxpj4hn->_3ru7myo; i++)
        {   _zzf8adf[i] = _zzf8adf[i-1] + (_vxpj4hn->_5izolnn*sizeof(unsigned short int) + _op767j2[i-1]*6*sizeof(float));
        }
        for (i = 1u; i < _vxpj4hn->_5izolnn; i++)
        {   _lmz04x0[i] = _lmz04x0[i-1] + (_vxpj4hn->_5izolnn*sizeof(unsigned short int) + _9f77lu8[i-1]*9*sizeof(float));
            _qzfhomt[i] = _qzfhomt[i-1] + (_vxpj4hn->_3ru7myo*sizeof(unsigned short int) + _el5015o[i-1]*6*sizeof(float));
        }
        
        for (i = 0u; i < _1poxrte[_fsqg8x8]; i++)      {   _jssuuqo[_fftouhy[i] +0*_vxpj4hn->_5izolnn] = _7yt2m6t[i*16 +5];      _jssuuqo[_fftouhy[i] +1*_vxpj4hn->_5izolnn] = _7yt2m6t[i*16 +6];          _jssuuqo[_fftouhy[i] +2*_vxpj4hn->_5izolnn] = _7yt2m6t[i*16 +7];            }
        MPI_Allreduce(MPI_IN_PLACE, _jssuuqo, (3*_vxpj4hn->_5izolnn), MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    }
    
    MPI_File_delete(_g79or78, MPI_INFO_NULL);
    MPI_File_open(MPI_COMM_WORLD, _g79or78, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &_6zx1tv1);  
    
    if (_fsqg8x8 == 0)
    {   MPI_File_write(_6zx1tv1,  &_vxpj4hn->_3ru7myo,     1, MPI_UNSIGNED,  &_0get8t6);      MPI_File_write(_6zx1tv1,  &_s55j2o7->_fu7bfxq,      1, MPI_UNSIGNED,  &_0get8t6);
        MPI_File_write(_6zx1tv1,  &_vxpj4hn->_5izolnn,     1, MPI_UNSIGNED,  &_0get8t6);      MPI_File_write(_6zx1tv1,  &_s55j2o7->_ryf0992,      1, MPI_UNSIGNED,  &_0get8t6);
    } 
    _9kwv08h = (4 * sizeof(unsigned int));
    
    MPI_File_write_at(_6zx1tv1,      (_9kwv08h + _paih0p4[_fsqg8x8]*sizeof(unsigned short int)),       &_dg1wjui[_paih0p4[_fsqg8x8]], _lexdw6l[_fsqg8x8], MPI_UNSIGNED_SHORT,  &_0get8t6);        _9kwv08h += (_vxpj4hn->_3ru7myo * sizeof(unsigned short int));
    if (_vxpj4hn->_5izolnn > 0u)
    {   MPI_File_write_at(_6zx1tv1,  (_9kwv08h + _paih0p4[_fsqg8x8]*sizeof(unsigned short int)),       &_op767j2[_paih0p4[_fsqg8x8]], _lexdw6l[_fsqg8x8], MPI_UNSIGNED_SHORT,  &_0get8t6);        _9kwv08h += (_vxpj4hn->_3ru7myo * sizeof(unsigned short int));
        MPI_File_write_at(_6zx1tv1,  (_9kwv08h + _20n8tte[_fsqg8x8]*sizeof(unsigned short int)),       &_9f77lu8[_20n8tte[_fsqg8x8]], _1poxrte[_fsqg8x8], MPI_UNSIGNED_SHORT,  &_0get8t6);        _9kwv08h += (_vxpj4hn->_5izolnn * sizeof(unsigned short int));
        MPI_File_write_at(_6zx1tv1,  (_9kwv08h + _20n8tte[_fsqg8x8]*sizeof(unsigned short int)),       &_el5015o[_20n8tte[_fsqg8x8]], _1poxrte[_fsqg8x8], MPI_UNSIGNED_SHORT,  &_0get8t6);        _9kwv08h += (_vxpj4hn->_5izolnn * sizeof(unsigned short int));
    }
    
    MPI_File_write_at(_6zx1tv1,      (_9kwv08h + _paih0p4[_fsqg8x8]*sizeof(float)), &_0wfoy4d[0*_vxpj4hn->_3ru7myo +_paih0p4[_fsqg8x8]], _lexdw6l[_fsqg8x8], MPI_FLOAT,          &_0get8t6);        _9kwv08h += (_vxpj4hn->_3ru7myo * sizeof(float));
    MPI_File_write_at(_6zx1tv1,      (_9kwv08h + _paih0p4[_fsqg8x8]*sizeof(float)), &_0wfoy4d[1*_vxpj4hn->_3ru7myo +_paih0p4[_fsqg8x8]], _lexdw6l[_fsqg8x8], MPI_FLOAT,          &_0get8t6);        _9kwv08h += (_vxpj4hn->_3ru7myo * sizeof(float));
    MPI_File_write_at(_6zx1tv1,      (_9kwv08h + _paih0p4[_fsqg8x8]*sizeof(float)), &_0wfoy4d[2*_vxpj4hn->_3ru7myo +_paih0p4[_fsqg8x8]], _lexdw6l[_fsqg8x8], MPI_FLOAT,          &_0get8t6);        _9kwv08h += (_vxpj4hn->_3ru7myo * sizeof(float));
    if (_vxpj4hn->_5izolnn > 0u)
    {   MPI_File_write_at(_6zx1tv1,  (_9kwv08h + _20n8tte[_fsqg8x8]*sizeof(float)), &_jssuuqo[0*_vxpj4hn->_5izolnn +_20n8tte[_fsqg8x8]], _1poxrte[_fsqg8x8], MPI_FLOAT,          &_0get8t6);        _9kwv08h += (_vxpj4hn->_5izolnn * sizeof(float));
        MPI_File_write_at(_6zx1tv1,  (_9kwv08h + _20n8tte[_fsqg8x8]*sizeof(float)), &_jssuuqo[1*_vxpj4hn->_5izolnn +_20n8tte[_fsqg8x8]], _1poxrte[_fsqg8x8], MPI_FLOAT,          &_0get8t6);        _9kwv08h += (_vxpj4hn->_5izolnn * sizeof(float));
        MPI_File_write_at(_6zx1tv1,  (_9kwv08h + _20n8tte[_fsqg8x8]*sizeof(float)), &_jssuuqo[2*_vxpj4hn->_5izolnn +_20n8tte[_fsqg8x8]], _1poxrte[_fsqg8x8], MPI_FLOAT,          &_0get8t6);        _9kwv08h += (_vxpj4hn->_5izolnn * sizeof(float));
    }
    
    for (i = 0u; i < _lexdw6l[_fsqg8x8]; i++)  
    {   _9wszv83 = 0;
        MPI_File_write_at(_6zx1tv1,  (_9kwv08h +_7pchcbn[_1nlc8jp[i]] +_9wszv83),  _5gzpl1q[i],  _vxpj4hn->_3ru7myo, MPI_UNSIGNED_SHORT,                         &_0get8t6);         _9wszv83  += _vxpj4hn->_3ru7myo *sizeof(unsigned short int);
        MPI_File_write_at(_6zx1tv1,  (_9kwv08h +_7pchcbn[_1nlc8jp[i]] +_9wszv83),  _ckeh5sw[i], 2*_qnsn2s8[i], MPI_FLOAT,                                   &_0get8t6);         _9wszv83  += (_qnsn2s8[i]*2*sizeof(float));
        MPI_File_write_at(_6zx1tv1,  (_9kwv08h +_7pchcbn[_1nlc8jp[i]] +_9wszv83),  _3ac6pe6[i], 2*_qnsn2s8[i], MPI_FLOAT,                                   &_0get8t6);         _9wszv83  += (_qnsn2s8[i]*2*sizeof(float));
        MPI_File_write_at(_6zx1tv1,  (_9kwv08h +_7pchcbn[_1nlc8jp[i]] +_9wszv83),  _obyghgk[i], 2*_qnsn2s8[i], MPI_FLOAT,                                   &_0get8t6);         _9wszv83  += (_qnsn2s8[i]*2*sizeof(float));
    }
    _9kwv08h += (_7pchcbn[_vxpj4hn->_3ru7myo-1] + _vxpj4hn->_3ru7myo*sizeof(unsigned short int) + 6*_dg1wjui[_vxpj4hn->_3ru7myo-1]*sizeof(float)  );
    
    if (_vxpj4hn->_5izolnn > 0u)
    {   for (i = 0u; i < _lexdw6l[_fsqg8x8]; i++)  
        {   _9wszv83 = 0;
            MPI_File_write_at(_6zx1tv1,  (_9kwv08h +_zzf8adf[_1nlc8jp[i]] +_9wszv83),  _oqgvzhq[i],  _vxpj4hn->_5izolnn, MPI_UNSIGNED_SHORT,                         &_0get8t6);        _9wszv83  += _vxpj4hn->_5izolnn *sizeof(unsigned short int);
            MPI_File_write_at(_6zx1tv1,  (_9kwv08h +_zzf8adf[_1nlc8jp[i]] +_9wszv83),  _yj9boq0[i], 3*_qpwbsf3[i], MPI_FLOAT,                                   &_0get8t6);        _9wszv83  += (_qpwbsf3[i]*3*sizeof(float));
            MPI_File_write_at(_6zx1tv1,  (_9kwv08h +_zzf8adf[_1nlc8jp[i]] +_9wszv83),  _ljefy3z[i], 3*_qpwbsf3[i], MPI_FLOAT,                                   &_0get8t6);        _9wszv83  += (_qpwbsf3[i]*3*sizeof(float));
        }
        _9kwv08h += (_zzf8adf[_vxpj4hn->_3ru7myo-1] + _vxpj4hn->_5izolnn*sizeof(unsigned short int) + 6*_op767j2[_vxpj4hn->_3ru7myo-1]*sizeof(float)  );
        
        if (_vxpj4hn->_5izolnn > 0)
        {   for (i = 0u; i < _1poxrte[_fsqg8x8]; i++)  
            {   _9wszv83 = 0;
                MPI_File_write_at(_6zx1tv1,  (_9kwv08h +_lmz04x0[_fftouhy[i]] +_9wszv83),  _1ajo091[i],  _vxpj4hn->_5izolnn, MPI_UNSIGNED_SHORT,                     &_0get8t6);        _9wszv83  += _vxpj4hn->_5izolnn *sizeof(unsigned short  int);
                MPI_File_write_at(_6zx1tv1,  (_9kwv08h +_lmz04x0[_fftouhy[i]] +_9wszv83),  _h61iaoi[i], 3*_iymyrsv[i], MPI_FLOAT,                               &_0get8t6);        _9wszv83  += (_iymyrsv[i]*3*sizeof(float));
                MPI_File_write_at(_6zx1tv1,  (_9kwv08h +_lmz04x0[_fftouhy[i]] +_9wszv83),  _plozjcn[i], 3*_iymyrsv[i], MPI_FLOAT,                               &_0get8t6);        _9wszv83  += (_iymyrsv[i]*3*sizeof(float));
                MPI_File_write_at(_6zx1tv1,  (_9kwv08h +_lmz04x0[_fftouhy[i]] +_9wszv83),  _qk4ltah[i], 3*_iymyrsv[i], MPI_FLOAT,                               &_0get8t6);        _9wszv83  += (_iymyrsv[i]*3*sizeof(float));
            }
            _9kwv08h += (_lmz04x0[_vxpj4hn->_5izolnn-1] + _vxpj4hn->_5izolnn*sizeof(unsigned short int) + 9*_9f77lu8[_vxpj4hn->_5izolnn-1]*sizeof(float)  );
            
            for (i = 0u; i < _1poxrte[_fsqg8x8]; i++)  
            {   _9wszv83 = 0;
                MPI_File_write_at(_6zx1tv1,  (_9kwv08h +_qzfhomt[_fftouhy[i]] +_9wszv83),  _3b9so0u[i],  _vxpj4hn->_3ru7myo, MPI_UNSIGNED_SHORT,                     &_0get8t6);        _9wszv83  += _vxpj4hn->_3ru7myo *sizeof(unsigned short  int);
                MPI_File_write_at(_6zx1tv1,  (_9kwv08h +_qzfhomt[_fftouhy[i]] +_9wszv83),  _0ljch2v[i], 2*_iob0nlw[i], MPI_FLOAT,                               &_0get8t6);        _9wszv83  += (_iob0nlw[i]*2*sizeof(float));
                MPI_File_write_at(_6zx1tv1,  (_9kwv08h +_qzfhomt[_fftouhy[i]] +_9wszv83),  _0dtirc9[i], 2*_iob0nlw[i], MPI_FLOAT,                               &_0get8t6);        _9wszv83  += (_iob0nlw[i]*2*sizeof(float));
                MPI_File_write_at(_6zx1tv1,  (_9kwv08h +_qzfhomt[_fftouhy[i]] +_9wszv83),  _9dde3cu[i], 2*_iob0nlw[i], MPI_FLOAT,                               &_0get8t6);        _9wszv83  += (_iob0nlw[i]*2*sizeof(float));
            }
        }
    }
    
    MPI_Barrier( MPI_COMM_WORLD );
    MPI_File_close(&_6zx1tv1);
    
    free(_dg1wjui);   free(_op767j2);   free(_9f77lu8);   free(_el5015o);   free(_7pchcbn);   free(_zzf8adf);   free(_lmz04x0);   free(_qzfhomt);   free(_0wfoy4d);   free(_jssuuqo);
    
    
    return;
}



void   _l1n13jc(int _0rmhrcn, int _fsqg8x8, const char *restrict _g79or78, struct _5r5uqpg *restrict _atj9tze, struct _xjuy0bv *restrict _vxpj4hn, const struct _j49itn7 *restrict _s55j2o7, const int *restrict _paih0p4, const int *restrict _20n8tte, const int *restrict _lexdw6l, const int *restrict _1poxrte, const unsigned int *restrict _1nlc8jp, unsigned int *restrict _p422ta8, float *restrict _6fokzix, float *restrict _0dl5e11, const unsigned int *restrict _fftouhy, float *restrict _7yt2m6t, unsigned short int *restrict _qnsn2s8, unsigned short int *restrict _qpwbsf3, unsigned short int *restrict _iymyrsv, unsigned short int *restrict _iob0nlw, unsigned short int ***restrict _5gzpl1q, unsigned short int ***restrict _oqgvzhq, unsigned short int ***restrict _1ajo091, unsigned short int ***restrict _3b9so0u, float ***restrict _ckeh5sw, float ***restrict _3ac6pe6, float ***restrict _obyghgk, float ***restrict _yj9boq0, float ***restrict _ljefy3z, float ***restrict _h61iaoi, float ***restrict _plozjcn, float ***restrict _qk4ltah, float ***restrict _0ljch2v, float ***restrict _0dtirc9, float ***restrict _9dde3cu)
{   
    unsigned int i;
    unsigned int _w8dt8cq,   _txrgn74,   _7bzegbl,   _qc6dr08;
    long int _moxsfh4;
    
    MPI_Status  _0get8t6;
    MPI_Offset  _9kwv08h;
    MPI_File  _6zx1tv1;
    
    MPI_File_open(MPI_COMM_WORLD, _g79or78, MPI_MODE_RDONLY, MPI_INFO_NULL, &_6zx1tv1);  
    
    MPI_File_read(_6zx1tv1, &_w8dt8cq,   1, MPI_UNSIGNED, &_0get8t6);       MPI_File_read(_6zx1tv1, &_7bzegbl,    1, MPI_UNSIGNED, &_0get8t6);
    MPI_File_read(_6zx1tv1, &_txrgn74,   1, MPI_UNSIGNED, &_0get8t6);       MPI_File_read(_6zx1tv1, &_qc6dr08,    1, MPI_UNSIGNED, &_0get8t6);
    if ((_w8dt8cq != _vxpj4hn->_3ru7myo) || (_txrgn74 != _vxpj4hn->_5izolnn) || (_7bzegbl != _s55j2o7->_fu7bfxq) || (_qc6dr08 != _s55j2o7->_ryf0992))
    {   perror("Kh_matrix file was openend, but fault/boundary element numbers not matching => exit\n");        exit(10);    }
    
    _moxsfh4   = _z0yncxl(ALIGN64, _lexdw6l[_fsqg8x8], sizeof(unsigned short int *));
    *_5gzpl1q   = aligned_alloc(ALIGN64, _moxsfh4);
    _moxsfh4   = _z0yncxl(ALIGN64, _lexdw6l[_fsqg8x8], sizeof(float *));
    *_ckeh5sw  = aligned_alloc(ALIGN64, _moxsfh4);
    *_3ac6pe6  = aligned_alloc(ALIGN64, _moxsfh4);
    *_obyghgk  = aligned_alloc(ALIGN64, _moxsfh4);
    
    if (_vxpj4hn->_5izolnn > 0u)
    {   _moxsfh4   = _z0yncxl(ALIGN64, _lexdw6l[_fsqg8x8], sizeof(unsigned short int *));
        *_oqgvzhq   = aligned_alloc(ALIGN64, _moxsfh4);
        _moxsfh4   = _z0yncxl(ALIGN64, _1poxrte[_fsqg8x8], sizeof(unsigned short int *));
        *_1ajo091   = aligned_alloc(ALIGN64, _moxsfh4);
        *_3b9so0u   = aligned_alloc(ALIGN64, _moxsfh4);
        _moxsfh4   = _z0yncxl(ALIGN64, _lexdw6l[_fsqg8x8], sizeof(float *));
        *_yj9boq0  = aligned_alloc(ALIGN64, _moxsfh4);
        *_ljefy3z  = aligned_alloc(ALIGN64, _moxsfh4);
        _moxsfh4   = _z0yncxl(ALIGN64, _1poxrte[_fsqg8x8], sizeof(float *));
        *_h61iaoi  = aligned_alloc(ALIGN64, _moxsfh4);
        *_plozjcn  = aligned_alloc(ALIGN64, _moxsfh4);
        *_qk4ltah  = aligned_alloc(ALIGN64, _moxsfh4);
        *_0ljch2v  = aligned_alloc(ALIGN64, _moxsfh4);
        *_0dtirc9  = aligned_alloc(ALIGN64, _moxsfh4);
        *_9dde3cu  = aligned_alloc(ALIGN64, _moxsfh4);
    }
    
    
    _9kwv08h = (4 * sizeof(unsigned int)) ;
    
    unsigned long long _9wszv83;
    unsigned short int *_dg1wjui = NULL,   *_op767j2 = NULL,   *_9f77lu8 = NULL,   *_el5015o = NULL;
    unsigned long long *_7pchcbn = NULL,   *_zzf8adf = NULL,   *_lmz04x0 = NULL,   *_qzfhomt = NULL;
    float *_0wfoy4d = NULL,   *_jssuuqo = NULL;
    
    _dg1wjui  = malloc(_vxpj4hn->_3ru7myo *sizeof(unsigned short int));        memset(_dg1wjui, 0, _vxpj4hn->_3ru7myo *sizeof(unsigned short int));
    _7pchcbn  = malloc(_vxpj4hn->_3ru7myo *sizeof(unsigned long long));        memset(_7pchcbn, 0, _vxpj4hn->_3ru7myo *sizeof(unsigned long long));
    _0wfoy4d = malloc((3*_vxpj4hn->_3ru7myo) *sizeof(float));                 memset(_0wfoy4d,0,(3*_vxpj4hn->_3ru7myo) *sizeof(float));
    
    if (_vxpj4hn->_5izolnn > 0u)
    {   _op767j2  = malloc(_vxpj4hn->_3ru7myo *sizeof(unsigned short int));        memset(_op767j2, 0, _vxpj4hn->_3ru7myo *sizeof(unsigned short int));
        _zzf8adf  = malloc(_vxpj4hn->_3ru7myo *sizeof(unsigned long long));        memset(_zzf8adf, 0, _vxpj4hn->_3ru7myo *sizeof(unsigned long long));
        _9f77lu8  = malloc(_vxpj4hn->_5izolnn *sizeof(unsigned short int));        memset(_9f77lu8, 0, _vxpj4hn->_5izolnn *sizeof(unsigned short int));
        _lmz04x0  = malloc(_vxpj4hn->_5izolnn *sizeof(unsigned long long));        memset(_lmz04x0, 0, _vxpj4hn->_5izolnn *sizeof(unsigned long long));
        _el5015o  = malloc(_vxpj4hn->_5izolnn *sizeof(unsigned short int));        memset(_el5015o, 0, _vxpj4hn->_5izolnn *sizeof(unsigned short int));
        _qzfhomt  = malloc(_vxpj4hn->_5izolnn *sizeof(unsigned long long));        memset(_qzfhomt, 0, _vxpj4hn->_5izolnn *sizeof(unsigned long long));
        _jssuuqo = malloc((3*_vxpj4hn->_5izolnn) *sizeof(float));                 memset(_jssuuqo,0,(3*_vxpj4hn->_5izolnn) *sizeof(float));
    }
    
    MPI_File_read_at(_6zx1tv1,      (_9kwv08h + _paih0p4[_fsqg8x8]*sizeof(unsigned short int)),   &_dg1wjui[_paih0p4[_fsqg8x8]], _lexdw6l[_fsqg8x8], MPI_UNSIGNED_SHORT,  &_0get8t6);        _9kwv08h += (_vxpj4hn->_3ru7myo * sizeof(unsigned short int));
    if (_vxpj4hn->_5izolnn > 0u)
    {   MPI_File_read_at(_6zx1tv1,  (_9kwv08h + _paih0p4[_fsqg8x8]*sizeof(unsigned short int)),   &_op767j2[_paih0p4[_fsqg8x8]], _lexdw6l[_fsqg8x8], MPI_UNSIGNED_SHORT,  &_0get8t6);        _9kwv08h += (_vxpj4hn->_3ru7myo * sizeof(unsigned short int));
        MPI_File_read_at(_6zx1tv1,  (_9kwv08h + _20n8tte[_fsqg8x8]*sizeof(unsigned short int)),   &_9f77lu8[_20n8tte[_fsqg8x8]], _1poxrte[_fsqg8x8], MPI_UNSIGNED_SHORT,  &_0get8t6);        _9kwv08h += (_vxpj4hn->_5izolnn * sizeof(unsigned short int));
        MPI_File_read_at(_6zx1tv1,  (_9kwv08h + _20n8tte[_fsqg8x8]*sizeof(unsigned short int)),   &_el5015o[_20n8tte[_fsqg8x8]], _1poxrte[_fsqg8x8], MPI_UNSIGNED_SHORT,  &_0get8t6);        _9kwv08h += (_vxpj4hn->_5izolnn * sizeof(unsigned short int));
    }
    
    MPI_Allreduce(MPI_IN_PLACE,     _dg1wjui, _vxpj4hn->_3ru7myo,  MPI_UNSIGNED_SHORT, MPI_SUM, MPI_COMM_WORLD);
    if (_vxpj4hn->_5izolnn > 0u)
    {   MPI_Allreduce(MPI_IN_PLACE, _op767j2, _vxpj4hn->_3ru7myo,  MPI_UNSIGNED_SHORT, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, _9f77lu8, _vxpj4hn->_5izolnn,  MPI_UNSIGNED_SHORT, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, _el5015o, _vxpj4hn->_5izolnn,  MPI_UNSIGNED_SHORT, MPI_SUM, MPI_COMM_WORLD);
    }
    
    for (i = 1u; i < _vxpj4hn->_3ru7myo; i++)
    {   _7pchcbn[i] = _7pchcbn[i-1] + (_vxpj4hn->_3ru7myo*sizeof(unsigned short int) + _dg1wjui[i-1]*6*sizeof(float));
    }
    if (_vxpj4hn->_5izolnn > 0u)
    {   for (i = 1u; i < _vxpj4hn->_3ru7myo; i++)
        {   _zzf8adf[i] = _zzf8adf[i-1] + (_vxpj4hn->_5izolnn*sizeof(unsigned short int) + _op767j2[i-1]*6*sizeof(float));
        }
        for (i = 1u; i < _vxpj4hn->_5izolnn; i++)
        {   _lmz04x0[i] = _lmz04x0[i-1] + (_vxpj4hn->_5izolnn*sizeof(unsigned short int) + _9f77lu8[i-1]*9*sizeof(float));
            _qzfhomt[i] = _qzfhomt[i-1] + (_vxpj4hn->_3ru7myo*sizeof(unsigned short int) + _el5015o[i-1]*6*sizeof(float));
        }
    }
    
    MPI_File_read_at(_6zx1tv1,      (_9kwv08h + _paih0p4[_fsqg8x8]*sizeof(float)),     &_0wfoy4d[0*_vxpj4hn->_3ru7myo +_paih0p4[_fsqg8x8]], _lexdw6l[_fsqg8x8], MPI_FLOAT,     &_0get8t6);        _9kwv08h += (_vxpj4hn->_3ru7myo * sizeof(float));
    MPI_File_read_at(_6zx1tv1,      (_9kwv08h + _paih0p4[_fsqg8x8]*sizeof(float)),     &_0wfoy4d[1*_vxpj4hn->_3ru7myo +_paih0p4[_fsqg8x8]], _lexdw6l[_fsqg8x8], MPI_FLOAT,     &_0get8t6);        _9kwv08h += (_vxpj4hn->_3ru7myo * sizeof(float));
    MPI_File_read_at(_6zx1tv1,      (_9kwv08h + _paih0p4[_fsqg8x8]*sizeof(float)),     &_0wfoy4d[2*_vxpj4hn->_3ru7myo +_paih0p4[_fsqg8x8]], _lexdw6l[_fsqg8x8], MPI_FLOAT,     &_0get8t6);        _9kwv08h += (_vxpj4hn->_3ru7myo * sizeof(float));
    MPI_Allreduce(MPI_IN_PLACE, _0wfoy4d, (3*_vxpj4hn->_3ru7myo), MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);

    if (_vxpj4hn->_5izolnn > 0u)
    {   MPI_File_read_at(_6zx1tv1,  (_9kwv08h + _20n8tte[_fsqg8x8]*sizeof(float)),     &_jssuuqo[0*_vxpj4hn->_5izolnn +_20n8tte[_fsqg8x8]], _1poxrte[_fsqg8x8], MPI_FLOAT,     &_0get8t6);        _9kwv08h += (_vxpj4hn->_5izolnn * sizeof(float));
        MPI_File_read_at(_6zx1tv1,  (_9kwv08h + _20n8tte[_fsqg8x8]*sizeof(float)),     &_jssuuqo[1*_vxpj4hn->_5izolnn +_20n8tte[_fsqg8x8]], _1poxrte[_fsqg8x8], MPI_FLOAT,     &_0get8t6);        _9kwv08h += (_vxpj4hn->_5izolnn * sizeof(float));
        MPI_File_read_at(_6zx1tv1,  (_9kwv08h + _20n8tte[_fsqg8x8]*sizeof(float)),     &_jssuuqo[2*_vxpj4hn->_5izolnn +_20n8tte[_fsqg8x8]], _1poxrte[_fsqg8x8], MPI_FLOAT,     &_0get8t6);        _9kwv08h += (_vxpj4hn->_5izolnn * sizeof(float));
        MPI_Allreduce(MPI_IN_PLACE, _jssuuqo, (3*_vxpj4hn->_5izolnn), MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    }
    
    for (i = 0u; i < _lexdw6l[_fsqg8x8]; i++)
    {   _qnsn2s8[i]  = _dg1wjui[_1nlc8jp[i]];
        _6fokzix[i*16 +5] = _0wfoy4d[_1nlc8jp[i] +0*_vxpj4hn->_3ru7myo];     _6fokzix[i*16 +6] = _0wfoy4d[_1nlc8jp[i] +1*_vxpj4hn->_3ru7myo];         _6fokzix[i*16 +7] = _0wfoy4d[_1nlc8jp[i] +2*_vxpj4hn->_3ru7myo];
    }
    if (_vxpj4hn->_5izolnn > 0u)
    {   for (i = 0u; i < _lexdw6l[_fsqg8x8]; i++)
        {   _qpwbsf3[i]  = _op767j2[_1nlc8jp[i]];      }
        for (i = 0u; i < _1poxrte[_fsqg8x8]; i++)
        {   _iymyrsv[i]  = _9f77lu8[_fftouhy[i]];                         _iob0nlw[i]  = _el5015o[_fftouhy[i]];
            _7yt2m6t[i*16 +5] = _jssuuqo[_fftouhy[i] +0*_vxpj4hn->_5izolnn];     _7yt2m6t[i*16 +6] = _jssuuqo[_fftouhy[i] +1*_vxpj4hn->_5izolnn];         _7yt2m6t[i*16 +7] = _jssuuqo[_fftouhy[i] +2*_vxpj4hn->_5izolnn];
        }
    }
    
    for (i = 0u; i < _lexdw6l[_fsqg8x8]; i++)  
    {   _9wszv83 = 0;
        _moxsfh4       = _z0yncxl(ALIGN64, _vxpj4hn->_3ru7myo, sizeof(unsigned short int));
        (*_5gzpl1q)[i]  = aligned_alloc(ALIGN64, _moxsfh4);
        _moxsfh4       = _z0yncxl(ALIGN64, 2*_qnsn2s8[i], sizeof(float));
        (*_ckeh5sw)[i] = aligned_alloc(ALIGN64, _moxsfh4);
        (*_3ac6pe6)[i] = aligned_alloc(ALIGN64, _moxsfh4);
        (*_obyghgk)[i] = aligned_alloc(ALIGN64, _moxsfh4);
        MPI_File_read_at(_6zx1tv1,  (_9kwv08h +_7pchcbn[_1nlc8jp[i]] +_9wszv83),  (*_5gzpl1q)[i],    _vxpj4hn->_3ru7myo, MPI_UNSIGNED_SHORT,  &_0get8t6);          _9wszv83  += _vxpj4hn->_3ru7myo *sizeof(unsigned short int);
        MPI_File_read_at(_6zx1tv1,  (_9kwv08h +_7pchcbn[_1nlc8jp[i]] +_9wszv83),  (*_ckeh5sw)[i], 2*_qnsn2s8[i], MPI_FLOAT,              &_0get8t6);          _9wszv83  += (_qnsn2s8[i]*2*sizeof(float));
        MPI_File_read_at(_6zx1tv1,  (_9kwv08h +_7pchcbn[_1nlc8jp[i]] +_9wszv83),  (*_3ac6pe6)[i], 2*_qnsn2s8[i], MPI_FLOAT,              &_0get8t6);          _9wszv83  += (_qnsn2s8[i]*2*sizeof(float));
        MPI_File_read_at(_6zx1tv1,  (_9kwv08h +_7pchcbn[_1nlc8jp[i]] +_9wszv83),  (*_obyghgk)[i], 2*_qnsn2s8[i], MPI_FLOAT,              &_0get8t6);          _9wszv83  += (_qnsn2s8[i]*2*sizeof(float));
    }
    _9kwv08h += (_7pchcbn[_vxpj4hn->_3ru7myo-1] + _vxpj4hn->_3ru7myo*sizeof(unsigned short int) + 6*_dg1wjui[_vxpj4hn->_3ru7myo-1]*sizeof(float)  );
    
    if (_vxpj4hn->_5izolnn > 0u)
    {   for (i = 0u; i < _lexdw6l[_fsqg8x8]; i++)  
        {   _9wszv83 = 0;
            _moxsfh4       = _z0yncxl(ALIGN64, _vxpj4hn->_5izolnn, sizeof(unsigned short int));
            (*_oqgvzhq)[i]  = aligned_alloc(ALIGN64, _moxsfh4);
            _moxsfh4       = _z0yncxl(ALIGN64, 3*_qpwbsf3[i], sizeof(float));
            (*_yj9boq0)[i] = aligned_alloc(ALIGN64, _moxsfh4);
            (*_ljefy3z)[i] = aligned_alloc(ALIGN64, _moxsfh4);
            MPI_File_read_at(_6zx1tv1,  (_9kwv08h +_zzf8adf[_1nlc8jp[i]] +_9wszv83),  (*_oqgvzhq)[i],    _vxpj4hn->_5izolnn, MPI_UNSIGNED_SHORT,  &_0get8t6);          _9wszv83  += _vxpj4hn->_5izolnn *sizeof(unsigned short int);
            MPI_File_read_at(_6zx1tv1,  (_9kwv08h +_zzf8adf[_1nlc8jp[i]] +_9wszv83),  (*_yj9boq0)[i], 3*_qpwbsf3[i], MPI_FLOAT,              &_0get8t6);          _9wszv83  += (_qpwbsf3[i]*3*sizeof(float));
            MPI_File_read_at(_6zx1tv1,  (_9kwv08h +_zzf8adf[_1nlc8jp[i]] +_9wszv83),  (*_ljefy3z)[i], 3*_qpwbsf3[i], MPI_FLOAT,              &_0get8t6);          _9wszv83  += (_qpwbsf3[i]*3*sizeof(float));
        }
        _9kwv08h += (_zzf8adf[_vxpj4hn->_3ru7myo-1] + _vxpj4hn->_5izolnn*sizeof(unsigned short int) + 6*_op767j2[_vxpj4hn->_3ru7myo-1]*sizeof(float)  );
        
        for (i = 0u; i < _1poxrte[_fsqg8x8]; i++)  
        {   _9wszv83 = 0;
            _moxsfh4       = _z0yncxl(ALIGN64, _vxpj4hn->_5izolnn, sizeof(unsigned short int));
            (*_1ajo091)[i]  = aligned_alloc(ALIGN64, _moxsfh4);
            _moxsfh4       = _z0yncxl(ALIGN64, 3*_iymyrsv[i], sizeof(float));
            (*_h61iaoi)[i] = aligned_alloc(ALIGN64, _moxsfh4);
            (*_plozjcn)[i] = aligned_alloc(ALIGN64, _moxsfh4);
            (*_qk4ltah)[i] = aligned_alloc(ALIGN64, _moxsfh4);
            MPI_File_read_at(_6zx1tv1,  (_9kwv08h +_lmz04x0[_fftouhy[i]] +_9wszv83),  (*_1ajo091)[i],    _vxpj4hn->_5izolnn, MPI_UNSIGNED_SHORT,  &_0get8t6);          _9wszv83  += _vxpj4hn->_5izolnn *sizeof(unsigned short int);
            MPI_File_read_at(_6zx1tv1,  (_9kwv08h +_lmz04x0[_fftouhy[i]] +_9wszv83),  (*_h61iaoi)[i], 3*_iymyrsv[i], MPI_FLOAT,              &_0get8t6);          _9wszv83  += (_iymyrsv[i]*3*sizeof(float));
            MPI_File_read_at(_6zx1tv1,  (_9kwv08h +_lmz04x0[_fftouhy[i]] +_9wszv83),  (*_plozjcn)[i], 3*_iymyrsv[i], MPI_FLOAT,              &_0get8t6);          _9wszv83  += (_iymyrsv[i]*3*sizeof(float));
            MPI_File_read_at(_6zx1tv1,  (_9kwv08h +_lmz04x0[_fftouhy[i]] +_9wszv83),  (*_qk4ltah)[i], 3*_iymyrsv[i], MPI_FLOAT,              &_0get8t6);          _9wszv83  += (_iymyrsv[i]*3*sizeof(float));
        }
        _9kwv08h += (_lmz04x0[_vxpj4hn->_5izolnn-1] + _vxpj4hn->_5izolnn*sizeof(unsigned short int) + 9*_9f77lu8[_vxpj4hn->_5izolnn-1]*sizeof(float)  );
        
        for (i = 0u; i < _1poxrte[_fsqg8x8]; i++)  
        {   _9wszv83 = 0;
            _moxsfh4       = _z0yncxl(ALIGN64, _vxpj4hn->_3ru7myo, sizeof(unsigned short int));
            (*_3b9so0u)[i]  = aligned_alloc(ALIGN64, _moxsfh4);
            _moxsfh4       = _z0yncxl(ALIGN64, 2*_iob0nlw[i], sizeof(float));
            (*_0ljch2v)[i] = aligned_alloc(ALIGN64, _moxsfh4);
            (*_0dtirc9)[i] = aligned_alloc(ALIGN64, _moxsfh4);
            (*_9dde3cu)[i] = aligned_alloc(ALIGN64, _moxsfh4);
            MPI_File_read_at(_6zx1tv1,  (_9kwv08h +_qzfhomt[_fftouhy[i]] +_9wszv83),  (*_3b9so0u)[i],    _vxpj4hn->_3ru7myo, MPI_UNSIGNED_SHORT,  &_0get8t6);          _9wszv83  += _vxpj4hn->_3ru7myo *sizeof(unsigned short int);
            MPI_File_read_at(_6zx1tv1,  (_9kwv08h +_qzfhomt[_fftouhy[i]] +_9wszv83),  (*_0ljch2v)[i], 2*_iob0nlw[i], MPI_FLOAT,              &_0get8t6);          _9wszv83  += (_iob0nlw[i]*2*sizeof(float));
            MPI_File_read_at(_6zx1tv1,  (_9kwv08h +_qzfhomt[_fftouhy[i]] +_9wszv83),  (*_0dtirc9)[i], 2*_iob0nlw[i], MPI_FLOAT,              &_0get8t6);          _9wszv83  += (_iob0nlw[i]*2*sizeof(float));
            MPI_File_read_at(_6zx1tv1,  (_9kwv08h +_qzfhomt[_fftouhy[i]] +_9wszv83),  (*_9dde3cu)[i], 2*_iob0nlw[i], MPI_FLOAT,              &_0get8t6);          _9wszv83  += (_iob0nlw[i]*2*sizeof(float));
        }
    }
    
    MPI_Barrier( MPI_COMM_WORLD );
    MPI_File_close(&_6zx1tv1);
    
    free(_dg1wjui);   free(_op767j2);   free(_9f77lu8);   free(_el5015o);   free(_7pchcbn);   free(_zzf8adf);   free(_lmz04x0);   free(_qzfhomt);   free(_0wfoy4d);   free(_jssuuqo);
    
    
    for (i = 0u; i < _lexdw6l[_fsqg8x8]; i++)  
    {   
        _0dl5e11[i*16 +1] = MAX(fabsf(_6fokzix[i*16 +5]),fabsf(_6fokzix[i*16 +6])); 
        _0dl5e11[i*16 +1] = MAX(_0dl5e11[i*16 +1],      fabsf(_6fokzix[i*16 +7])); 
        
        _myx444n(_atj9tze, _vxpj4hn, i, _p422ta8, _6fokzix, _0dl5e11, 1);
        
    } 
    
    
    unsigned int _1b07g2r[4],   _l2ld05n[4],   _ybxkpm2[4],   _jg2y8rg[4];
    memset(_1b07g2r,     0, 4*sizeof(unsigned int));
    memset(_l2ld05n,     0, 4*sizeof(unsigned int));
    memset(_ybxkpm2,     0, 4*sizeof(unsigned int));
    memset(_jg2y8rg, 0, 4*sizeof(unsigned int));
    
    for (i = 0u; i < _lexdw6l[_fsqg8x8]; i++)
    {   _1b07g2r[0] += _qnsn2s8[i];             _jg2y8rg[0] = (_jg2y8rg[0] >= _qnsn2s8[i])*_jg2y8rg[0]  +  (_jg2y8rg[0] < _qnsn2s8[i])*_qnsn2s8[i];
    }
    if (_vxpj4hn->_5izolnn > 0u)
    {   for (i = 0u; i < _lexdw6l[_fsqg8x8]; i++)
        {   _1b07g2r[1] += _qpwbsf3[i];             _jg2y8rg[1] = (_jg2y8rg[1] >= _qpwbsf3[i])*_jg2y8rg[1]  +  (_jg2y8rg[1] < _qpwbsf3[i])*_qpwbsf3[i];
        }
        for (i = 0u; i < _1poxrte[_fsqg8x8]; i++)      
        {   _1b07g2r[2] += _iymyrsv[i];             _jg2y8rg[2] = (_jg2y8rg[2] >= _iymyrsv[i])*_jg2y8rg[2]  +  (_jg2y8rg[2] < _iymyrsv[i])*_iymyrsv[i];
            _1b07g2r[3] += _iob0nlw[i];             _jg2y8rg[3] = (_jg2y8rg[3] >= _iob0nlw[i])*_jg2y8rg[3]  +  (_jg2y8rg[3] < _iob0nlw[i])*_iob0nlw[i];
        }
    }
    
    MPI_Allreduce(MPI_IN_PLACE, _jg2y8rg, 4, MPI_UNSIGNED, MPI_MAX, MPI_COMM_WORLD);
    _vxpj4hn->_dn1th1x = MAX(_jg2y8rg[0], _jg2y8rg[3]);
    _vxpj4hn->_w2maiet = MAX(_jg2y8rg[1], _jg2y8rg[2]);
    
    _1b07g2r[0] =                             _1b07g2r[0]/_lexdw6l[_fsqg8x8];
    _1b07g2r[1] =                             _1b07g2r[1]/_lexdw6l[_fsqg8x8];
    _1b07g2r[2] = (_1poxrte[_fsqg8x8] > 0u)  ?  (_1b07g2r[2]/_1poxrte[_fsqg8x8])  :  0u;
    _1b07g2r[3] = (_1poxrte[_fsqg8x8] > 0u)  ?  (_1b07g2r[3]/_1poxrte[_fsqg8x8])  :  0u;
    
    memcpy(_l2ld05n, _1b07g2r,  4*sizeof(unsigned int));
    memcpy(_ybxkpm2, _1b07g2r,  4*sizeof(unsigned int));
    MPI_Allreduce(MPI_IN_PLACE, _l2ld05n, 4, MPI_UNSIGNED, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, _ybxkpm2, 4, MPI_UNSIGNED, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, _1b07g2r, 4, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
    
    _1b07g2r[0] /= _0rmhrcn;        _1b07g2r[1] /= _0rmhrcn;        _1b07g2r[2] /= _0rmhrcn;        _1b07g2r[3] /= _0rmhrcn;
    
    if (_fsqg8x8 == 0)
    {   fprintf(stdout,"Kh_lengths (average and range) across RANKs\n");



        fprintf(stdout,"Minimum:         %u\n", _l2ld05n[0]);
        fprintf(stdout,"Average:         %u\n", _1b07g2r[0]);
        fprintf(stdout,"Maximum:         %u\n", _ybxkpm2[0]);
    }
    
    
    return;
}






void   _nb4gh2m(float _gnmby49, unsigned int _7c6qc6k, const float *restrict _bjfvyb6, unsigned int _6zjq7w7,  const unsigned int *restrict _r59c81d, const float *restrict _yn4u9bx, unsigned int ***restrict _h9ier7u, unsigned int ***restrict _xdrkz40, unsigned int *restrict _katqyds)
{   
    float _glqj5sx = 3.0f; 
    float _8kzruin     = 0.666f;
    float _3oeuvty     = 1.0f/_8kzruin;
    long int _moxsfh4;
    
    
    _moxsfh4   = _z0yncxl(ALIGN64, _6zjq7w7, sizeof(unsigned int *));
    *_h9ier7u    = aligned_alloc(ALIGN64, _moxsfh4);
    _moxsfh4   = _z0yncxl(ALIGN64, 16, sizeof(unsigned int));
    for (int i = 0u; i <  _6zjq7w7; i++)
    {   (*_h9ier7u)[i]  = aligned_alloc(ALIGN64, _moxsfh4);            memset((*_h9ier7u)[i], 0, _moxsfh4);
    }
    
    _moxsfh4   = _z0yncxl(ALIGN64, _6zjq7w7*4, sizeof(unsigned int *));
    *_xdrkz40   = aligned_alloc(ALIGN64, _moxsfh4); 
    _moxsfh4   = _z0yncxl(ALIGN64, 8, sizeof(unsigned int));
    for (int i = 0u; i < (_6zjq7w7*4); i++)
    {   (*_xdrkz40)[i] = aligned_alloc(ALIGN64, _moxsfh4);            memset((*_xdrkz40)[i],0, _moxsfh4);
    }
    
    unsigned int i,   j,   k,   _f7fr91i,   _jzr25md,   _gntb34k,   _d7x48wz,   _cs4kek4,   _omxlvz4;
    float _xpbwhvn,   _t8qzpky;
    unsigned int _hqr6f7q[2][4];
    float _088o98b[3],   _j3swebn[3],   _795e5iz[3],   _v3z6bz8[6];
    
    unsigned int *_e4ex6u8  = malloc(_7c6qc6k *sizeof(unsigned int)); 
    unsigned int **_ddqhkga = malloc(_7c6qc6k *sizeof(unsigned int *)); 
    memset(_e4ex6u8, 0, _7c6qc6k*sizeof(unsigned int));
    for (i = 0u; i < _6zjq7w7; i++)      {      _e4ex6u8[_r59c81d[i*8 +3]] += 1u;                                 }
    for (i = 0u; i < _7c6qc6k; i++)       {      _ddqhkga[i]  = malloc(_e4ex6u8[i]*sizeof(unsigned int));        }
    
    memset(_e4ex6u8, 0, _7c6qc6k*sizeof(unsigned int));
    for (i = 0u; i < _6zjq7w7; i++)
    {   _ddqhkga[_r59c81d[i*8 +3]][_e4ex6u8[_r59c81d[i*8 +3]]] = i;
        _e4ex6u8[_r59c81d[i*8 +3]]                            += 1u;
    }
    
    for (i = 0u; i < _7c6qc6k; i++)
    {   float _li076a3[3],   _h42m44g[4] = {FLT_MAX, -FLT_MAX, FLT_MAX, -FLT_MAX};
        
        unsigned int *_aegvfq3   = malloc(_e4ex6u8[i]* sizeof(unsigned int));        memset(_aegvfq3,  0, _e4ex6u8[i]* sizeof(unsigned int));
        unsigned int *_o6ouz1e  = malloc(_e4ex6u8[i]* sizeof(unsigned int));        memset(_o6ouz1e, 0, _e4ex6u8[i]* sizeof(unsigned int));
        unsigned int *_5wryesu   = malloc(_e4ex6u8[i]* sizeof(unsigned int));        memset(_5wryesu,  0, _e4ex6u8[i]* sizeof(unsigned int));
        
        unsigned int **_nwngucz = malloc(_e4ex6u8[i] *sizeof(unsigned int *));
        float **_zh13cd5  = malloc(_e4ex6u8[i] *sizeof(float *));
        
        for (j = 0u; j < _e4ex6u8[i]; j++)
        {   _moxsfh4    = _z0yncxl(ALIGN64, 16, sizeof(unsigned int));
            _nwngucz[j]       = aligned_alloc(ALIGN64, _moxsfh4);                       memset(_nwngucz[j],       0, _moxsfh4);
            _moxsfh4    = _z0yncxl(ALIGN64, 16, sizeof(float));
            _zh13cd5[j] = aligned_alloc(ALIGN64, _moxsfh4);                       memset(_zh13cd5[j], 0, _moxsfh4); 
        }
        
        memcpy(_088o98b, &_bjfvyb6[i*8 +3], 3*sizeof(float));
        _ktakpux(_088o98b, _j3swebn, _795e5iz);
        
        _v3z6bz8[0] = FLT_MAX;            _v3z6bz8[1] = 0.0f;            _v3z6bz8[2] = -FLT_MAX;
        _v3z6bz8[3] = FLT_MAX;            _v3z6bz8[4] = 0.0f;            _v3z6bz8[5] = -FLT_MAX;
        
        for (j = 0; j < _e4ex6u8[i]; j++) 
        {   _li076a3[0]    = _yn4u9bx[ _ddqhkga[i][j]*16 +0] - _bjfvyb6[i*8 +0];
            _li076a3[1]    = _yn4u9bx[ _ddqhkga[i][j]*16 +1] - _bjfvyb6[i*8 +1];
            _li076a3[2]    = _yn4u9bx[ _ddqhkga[i][j]*16 +2] - _bjfvyb6[i*8 +2];
            
            _zh13cd5[j][8] = _088o98b[0]*_li076a3[0] + _088o98b[1]*_li076a3[1] + _088o98b[2]*_li076a3[2];
            _zh13cd5[j][9] = _j3swebn[0]*_li076a3[0] + _j3swebn[1]*_li076a3[1] + _j3swebn[2]*_li076a3[2];
            _zh13cd5[j][10]= _795e5iz[0]*_li076a3[0] + _795e5iz[1]*_li076a3[1] + _795e5iz[2]*_li076a3[2];
            
            _h42m44g[0]      = MIN(_h42m44g[0], _zh13cd5[j][9]);     _h42m44g[1]      = MAX(_h42m44g[1], _zh13cd5[j][9]);
            _h42m44g[2]      = MIN(_h42m44g[2], _zh13cd5[j][10]);    _h42m44g[3]      = MAX(_h42m44g[3], _zh13cd5[j][10]);
        }
        
        _v3z6bz8[0]    = MIN(_v3z6bz8[0], _h42m44g[0]);        _v3z6bz8[2] = MAX(_v3z6bz8[2], _h42m44g[1]);
        _v3z6bz8[3]    = MIN(_v3z6bz8[3], _h42m44g[2]);        _v3z6bz8[5] = MAX(_v3z6bz8[5], _h42m44g[3]);
        
        _xpbwhvn         = MAX((_v3z6bz8[2] -_v3z6bz8[0]), (_v3z6bz8[5] -_v3z6bz8[3]));
        _zh13cd5[0][0] = _v3z6bz8[0];              _zh13cd5[0][1] = _v3z6bz8[2];
        _zh13cd5[0][2] = _v3z6bz8[3];              _zh13cd5[0][3] = _v3z6bz8[5];
        
        _jzr25md = 0u;
        _cs4kek4 = 0u;
        _f7fr91i   = 1u;
        
        while (0.5f*_xpbwhvn >= _glqj5sx*_gnmby49) 
        {   _jzr25md += 1u; 
            _xpbwhvn  *= 0.5f; 
            _gntb34k   = 0u; 
            
            for (j = 0u; j < _f7fr91i; j++)
            {   _v3z6bz8[0] = _zh13cd5[j][0];      _v3z6bz8[2] = _zh13cd5[j][1];       _v3z6bz8[1] = 0.5f*(_v3z6bz8[0] +_v3z6bz8[2]);
                _v3z6bz8[3] = _zh13cd5[j][2];      _v3z6bz8[5] = _zh13cd5[j][3];       _v3z6bz8[4] = 0.5f*(_v3z6bz8[3] +_v3z6bz8[5]);
                
                _d7x48wz = 0u; 
                memset(_5wryesu,  0, _e4ex6u8[i]* sizeof(unsigned int)); 
                _hqr6f7q[0][0] = 0u;      _hqr6f7q[0][1] = 0u;      _hqr6f7q[0][2] = 0u;      _hqr6f7q[0][3] = 0u;
                _hqr6f7q[1][0] = 0u;      _hqr6f7q[1][1] = 0u;      _hqr6f7q[1][2] = 0u;      _hqr6f7q[1][3] = 0u;
                
                for (k = 0u; k < _e4ex6u8[i]; k++)
                {   _omxlvz4        = (_nwngucz[k][_jzr25md-1] == _aegvfq3[j])*1u + (_nwngucz[k][_jzr25md-1] != _aegvfq3[j])*0u;
                    _5wryesu[_d7x48wz] = (_omxlvz4 == 1u)*k                     + (_omxlvz4 != 1u)*_5wryesu[_d7x48wz];
                    _d7x48wz         = (_omxlvz4 == 1u)*(_d7x48wz+1u)            + (_omxlvz4 != 1u)*_d7x48wz;
                }
                
                _t8qzpky = (_v3z6bz8[2]-_v3z6bz8[0]+_gnmby49)/(_v3z6bz8[5]-_v3z6bz8[3]+_gnmby49);
                
                if      (_t8qzpky >= _3oeuvty)
                {   for (k = 0u; k < _d7x48wz; k++)
                    {   if      (_zh13cd5[_5wryesu[k]][9] < _v3z6bz8[1])
                        {   if (_hqr6f7q[0][0] == 0u)
                            {   _cs4kek4     += 1u;
                                _hqr6f7q[0][0]  = 1u;                            _hqr6f7q[1][0]          = _cs4kek4;
                                _zh13cd5[_gntb34k][4] = _v3z6bz8[0];        _zh13cd5[_gntb34k][5] = _v3z6bz8[1];
                                _zh13cd5[_gntb34k][6] = _v3z6bz8[3];        _zh13cd5[_gntb34k][7] = _v3z6bz8[5];
                                _o6ouz1e[_gntb34k]      = _cs4kek4;              _gntb34k              += 1u;
                                (*_xdrkz40)[(_aegvfq3[j]+ *_katqyds)][0] += 1u;
                                (*_xdrkz40)[(_aegvfq3[j]+ *_katqyds)][(*_xdrkz40)[(_aegvfq3[j]+ *_katqyds)][0]] = (_cs4kek4 + *_katqyds);
                            }
                            _nwngucz[_5wryesu[k]][_jzr25md] = _hqr6f7q[1][0];
                        }
                        else if (_zh13cd5[_5wryesu[k]][9] >= _v3z6bz8[1])
                        {   if (_hqr6f7q[0][1] == 0u)
                            {   _cs4kek4     += 1u;
                                _hqr6f7q[0][1]  = 1u;                            _hqr6f7q[1][1]          = _cs4kek4;
                                _zh13cd5[_gntb34k][4] = _v3z6bz8[1];        _zh13cd5[_gntb34k][5] = _v3z6bz8[2];
                                _zh13cd5[_gntb34k][6] = _v3z6bz8[3];        _zh13cd5[_gntb34k][7] = _v3z6bz8[5];
                                _o6ouz1e[_gntb34k]      = _cs4kek4;              _gntb34k              += 1u;
                                (*_xdrkz40)[(_aegvfq3[j]+ *_katqyds)][0] += 1u;
                                (*_xdrkz40)[(_aegvfq3[j]+ *_katqyds)][(*_xdrkz40)[(_aegvfq3[j]+ *_katqyds)][0]] = (_cs4kek4 + *_katqyds);
                            }
                            _nwngucz[_5wryesu[k]][_jzr25md] = _hqr6f7q[1][1];
                }   }   }
                else if (_t8qzpky <= _8kzruin)
                {   for (k = 0u; k < _d7x48wz; k++)
                    {   if      (_zh13cd5[_5wryesu[k]][10] < _v3z6bz8[4])
                        {   if (_hqr6f7q[0][0] == 0u)
                            {   _cs4kek4     += 1u;
                                _hqr6f7q[0][0]  = 1u;                            _hqr6f7q[1][0]          = _cs4kek4;
                                _zh13cd5[_gntb34k][4] = _v3z6bz8[0];        _zh13cd5[_gntb34k][5] = _v3z6bz8[2];
                                _zh13cd5[_gntb34k][6] = _v3z6bz8[3];        _zh13cd5[_gntb34k][7] = _v3z6bz8[4];
                                _o6ouz1e[_gntb34k]      = _cs4kek4;              _gntb34k              += 1u;
                                (*_xdrkz40)[(_aegvfq3[j]+ *_katqyds)][0] += 1u;
                                (*_xdrkz40)[(_aegvfq3[j]+ *_katqyds)][(*_xdrkz40)[(_aegvfq3[j]+ *_katqyds)][0]] = (_cs4kek4 + *_katqyds);
                                                                   }
                            _nwngucz[_5wryesu[k]][_jzr25md] = _hqr6f7q[1][0];
                        }
                        else if (_zh13cd5[_5wryesu[k]][10] >= _v3z6bz8[4])
                        {   if (_hqr6f7q[0][1] == 0u)
                            {   _cs4kek4     += 1u;
                                _hqr6f7q[0][1]  = 1u;                            _hqr6f7q[1][1]          = _cs4kek4;
                                _zh13cd5[_gntb34k][4] = _v3z6bz8[0];        _zh13cd5[_gntb34k][5] = _v3z6bz8[2];
                                _zh13cd5[_gntb34k][6] = _v3z6bz8[4];        _zh13cd5[_gntb34k][7] = _v3z6bz8[5];
                                _o6ouz1e[_gntb34k]      = _cs4kek4;              _gntb34k              += 1u;
                                (*_xdrkz40)[(_aegvfq3[j]+ *_katqyds)][0] += 1u;
                                (*_xdrkz40)[(_aegvfq3[j]+ *_katqyds)][(*_xdrkz40)[(_aegvfq3[j]+ *_katqyds)][0]] = (_cs4kek4 + *_katqyds);
                            }
                            _nwngucz[_5wryesu[k]][_jzr25md] = _hqr6f7q[1][1];
                }   }   }
                else
                {   for (k = 0u; k < _d7x48wz; k++)
                    {   if      ( (_zh13cd5[_5wryesu[k]][9] >= _v3z6bz8[1])  &&  (_zh13cd5[_5wryesu[k]][10] >= _v3z6bz8[4]) )
                        {   if (_hqr6f7q[0][0] == 0u)
                            {   _cs4kek4     += 1u;
                                _hqr6f7q[0][0]  = 1u;                            _hqr6f7q[1][0]          = _cs4kek4;
                                _zh13cd5[_gntb34k][4] = _v3z6bz8[1];        _zh13cd5[_gntb34k][5] = _v3z6bz8[2];
                                _zh13cd5[_gntb34k][6] = _v3z6bz8[4];        _zh13cd5[_gntb34k][7] = _v3z6bz8[5];
                                _o6ouz1e[_gntb34k]      = _cs4kek4;              _gntb34k              += 1u;
                                (*_xdrkz40)[(_aegvfq3[j]+ *_katqyds)][0] +=1u;
                                (*_xdrkz40)[(_aegvfq3[j]+ *_katqyds)][(*_xdrkz40)[(_aegvfq3[j]+ *_katqyds)][0]] = (_cs4kek4 + *_katqyds);
                            }
                            _nwngucz[_5wryesu[k]][_jzr25md] = _hqr6f7q[1][0];
                        }
                        else if ( (_zh13cd5[_5wryesu[k]][9] <  _v3z6bz8[1])  &&  (_zh13cd5[_5wryesu[k]][10] >= _v3z6bz8[4]) )
                        {   if (_hqr6f7q[0][1] == 0u)
                            {   _cs4kek4     += 1u;
                                _hqr6f7q[0][1]  = 1u;                            _hqr6f7q[1][1]          = _cs4kek4;
                                _zh13cd5[_gntb34k][4] = _v3z6bz8[0];        _zh13cd5[_gntb34k][5] = _v3z6bz8[1];
                                _zh13cd5[_gntb34k][6] = _v3z6bz8[4];        _zh13cd5[_gntb34k][7] = _v3z6bz8[5];
                                _o6ouz1e[_gntb34k]      = _cs4kek4;              _gntb34k              += 1u;
                                (*_xdrkz40)[(_aegvfq3[j]+ *_katqyds)][0] +=1u;
                                (*_xdrkz40)[(_aegvfq3[j]+ *_katqyds)][(*_xdrkz40)[(_aegvfq3[j]+ *_katqyds)][0]] = (_cs4kek4 + *_katqyds);
                            }
                            _nwngucz[_5wryesu[k]][_jzr25md] = _hqr6f7q[1][1];
                        }
                        else if ( (_zh13cd5[_5wryesu[k]][9] <  _v3z6bz8[1])  &&  (_zh13cd5[_5wryesu[k]][10] <  _v3z6bz8[4]) )
                        {   if (_hqr6f7q[0][2] == 0u)
                            {   _cs4kek4     += 1u;
                                _hqr6f7q[0][2]  = 1u;                            _hqr6f7q[1][2]          = _cs4kek4;
                                _zh13cd5[_gntb34k][4] = _v3z6bz8[0];        _zh13cd5[_gntb34k][5] = _v3z6bz8[1];
                                _zh13cd5[_gntb34k][6] = _v3z6bz8[3];        _zh13cd5[_gntb34k][7] = _v3z6bz8[4];
                                _o6ouz1e[_gntb34k]      = _cs4kek4;              _gntb34k              += 1u;
                                (*_xdrkz40)[(_aegvfq3[j]+ *_katqyds)][0] +=1u;
                                (*_xdrkz40)[(_aegvfq3[j]+ *_katqyds)][(*_xdrkz40)[(_aegvfq3[j]+ *_katqyds)][0]] = (_cs4kek4 + *_katqyds);
                            }
                            _nwngucz[_5wryesu[k]][_jzr25md] = _hqr6f7q[1][2];
                        }
                        else if ( (_zh13cd5[_5wryesu[k]][9] >= _v3z6bz8[1])  &&  (_zh13cd5[_5wryesu[k]][10] <  _v3z6bz8[4]) )
                        {   if (_hqr6f7q[0][3] == 0u)
                            {   _cs4kek4     += 1u;
                                _hqr6f7q[0][3]  = 1u;                            _hqr6f7q[1][3]          = _cs4kek4;
                                _zh13cd5[_gntb34k][4] = _v3z6bz8[1];        _zh13cd5[_gntb34k][5] = _v3z6bz8[2];
                                _zh13cd5[_gntb34k][6] = _v3z6bz8[3];        _zh13cd5[_gntb34k][7] = _v3z6bz8[4];
                                _o6ouz1e[_gntb34k]      = _cs4kek4;              _gntb34k              += 1u;
                                (*_xdrkz40)[(_aegvfq3[j]+ *_katqyds)][0] +=1u;
                                (*_xdrkz40)[(_aegvfq3[j]+ *_katqyds)][(*_xdrkz40)[(_aegvfq3[j]+ *_katqyds)][0]] = (_cs4kek4 + *_katqyds);
                            }
                            _nwngucz[_5wryesu[k]][_jzr25md] = _hqr6f7q[1][3];
                }   }   }
                
            }
            memcpy(_aegvfq3, _o6ouz1e, _gntb34k* sizeof(unsigned int) );
            for (j = 0u; j < _gntb34k; j++)  {   memcpy(&_zh13cd5[j][0], &_zh13cd5[j][4], 4* sizeof(float) );         }
            _f7fr91i = _gntb34k;
        }
        
        _cs4kek4   += 1u;
        _jzr25md   += 2u;
        for (j = 0; j < _e4ex6u8[i]; j++ )
        {   _nwngucz[j][_jzr25md -1]        = (_cs4kek4 +j);
            (*_h9ier7u)[_ddqhkga[i][j]][0] = _jzr25md; 
            
            for (k = 0; k < _jzr25md; k++)
            {   (*_h9ier7u)[_ddqhkga[i][j]][(16 -_jzr25md +k)] = (_nwngucz[j][k] + *_katqyds);
        }   }
        
        *_katqyds += (_cs4kek4 + _e4ex6u8[i]);
        
        for (j = 0; j < _e4ex6u8[i]; j++)     {      free(_nwngucz[j]);   free(_zh13cd5[j]);     }
        free(_aegvfq3);   free(_o6ouz1e);   free(_5wryesu);   free(_nwngucz);   free(_zh13cd5);
    }
    
    for (i = 0u; i < _7c6qc6k; i++)       {       free(_ddqhkga[i]);        }
    free(_e4ex6u8);   free(_ddqhkga);
    
} 



void _h75lxo6(unsigned int _6zjq7w7,  const unsigned int *restrict _r59c81d, const float *restrict _yn4u9bx, const float *restrict _mszdxim, unsigned int **restrict _nwngucz, unsigned int _katqyds,  float ***restrict _vwvm4w3)
{   
    
    
    long int _moxsfh4;
    
    _moxsfh4   = _z0yncxl(ALIGN64, _katqyds, sizeof(float *));
    *_vwvm4w3 = aligned_alloc(ALIGN64, _moxsfh4);
    _moxsfh4   = _z0yncxl(ALIGN64, 21, sizeof(float));
    for (int i = 0; i < _katqyds; i++)
    {   (*_vwvm4w3)[i] = aligned_alloc(ALIGN64, _moxsfh4);         memset((*_vwvm4w3)[i], 0, _moxsfh4);
    }
    
    int i,   j;
    unsigned int _omxlvz4,   _qke2gjr;
    float _5gdej31,   _gymdm0l,   _gwaj1g1,   _r1dwwo8,   _nbhnal3,   _31yrujo;
    float _dan085r[3],   _ntdr5xy[3],   _fc4kl8c[3],   _gn0t6dv[3],   _9p3jpxg[3],   _088o98b[3],   _j3swebn[3],   _795e5iz[3],   _bk2l0d1[3],   _b4xgmpw[3],   _n42x2ez[3];
    
    for (i = 0u; i < _6zjq7w7; i++)
    {   _omxlvz4= _nwngucz[i][0];
        for (j = 0u; j < _omxlvz4; j++)
        {   _qke2gjr = _nwngucz[i][16 -_omxlvz4 +j];
            (*_vwvm4w3)[_qke2gjr][0] += _yn4u9bx[i*16 +0];      (*_vwvm4w3)[_qke2gjr][1] += _yn4u9bx[i*16 +1];      (*_vwvm4w3)[_qke2gjr][2] += _yn4u9bx[i*16 +2];
            (*_vwvm4w3)[_qke2gjr][3] += _yn4u9bx[i*16 +15];     (*_vwvm4w3)[_qke2gjr][4] += 1.0f;
            (*_vwvm4w3)[_qke2gjr][6] += _yn4u9bx[i*16 +6];      (*_vwvm4w3)[_qke2gjr][7] += _yn4u9bx[i*16 +7];      (*_vwvm4w3)[_qke2gjr][8] += _yn4u9bx[i*16 +8];
    }   }
    
    for (i = 0u; i < _katqyds; i++)
    {   (*_vwvm4w3)[i][0] /= (*_vwvm4w3)[i][4];               (*_vwvm4w3)[i][1] /= (*_vwvm4w3)[i][4];           (*_vwvm4w3)[i][2] /= (*_vwvm4w3)[i][4];
        _bfuk9e1( &(*_vwvm4w3)[i][6] );
        (*_vwvm4w3)[i][5]  = -FLT_MAX;
        (*_vwvm4w3)[i][9]  =  FLT_MAX;                       (*_vwvm4w3)[i][10] = -FLT_MAX;
        (*_vwvm4w3)[i][11] =  FLT_MAX;                       (*_vwvm4w3)[i][12] = -FLT_MAX;
    }
    
    for (i = 0u; i < _6zjq7w7; i++)
    {   memcpy(_dan085r, &_mszdxim[_r59c81d[i*8 +0]*4 +0],  3*sizeof(float));
        memcpy(_ntdr5xy, &_mszdxim[_r59c81d[i*8 +1]*4 +0],  3*sizeof(float));
        memcpy(_fc4kl8c, &_mszdxim[_r59c81d[i*8 +2]*4 +0],  3*sizeof(float));
        _omxlvz4 = _nwngucz[i][0];
        for (j = 0u; j < _omxlvz4; j++)
        {   _qke2gjr = _nwngucz[i][16 -_omxlvz4 +j];
            
            memcpy(_9p3jpxg, &(*_vwvm4w3)[_qke2gjr][0], 3*sizeof(float));
            memcpy(_088o98b, &(*_vwvm4w3)[_qke2gjr][6], 3*sizeof(float));
            _ktakpux(_088o98b, _j3swebn, _795e5iz);
            
            _id4n1u3(_bk2l0d1, _dan085r, _9p3jpxg);
            _b4xgmpw[0] = _j3swebn[0]*_bk2l0d1[0] + _j3swebn[1]*_bk2l0d1[1] +  _j3swebn[2]*_bk2l0d1[2]; 
            _b4xgmpw[1] = _795e5iz[0]*_bk2l0d1[0] + _795e5iz[1]*_bk2l0d1[1] +  _795e5iz[2]*_bk2l0d1[2]; 
            
            (*_vwvm4w3)[_qke2gjr][9]  = MIN((*_vwvm4w3)[_qke2gjr][9],  _b4xgmpw[0]);
            (*_vwvm4w3)[_qke2gjr][10] = MAX((*_vwvm4w3)[_qke2gjr][10], _b4xgmpw[0]);
            (*_vwvm4w3)[_qke2gjr][11] = MIN((*_vwvm4w3)[_qke2gjr][11], _b4xgmpw[1]);
            (*_vwvm4w3)[_qke2gjr][12] = MAX((*_vwvm4w3)[_qke2gjr][12], _b4xgmpw[1]);
            
            _id4n1u3(_bk2l0d1, _ntdr5xy, _9p3jpxg);
            _b4xgmpw[0] = _j3swebn[0]*_bk2l0d1[0] + _j3swebn[1]*_bk2l0d1[1] +  _j3swebn[2]*_bk2l0d1[2]; 
            _b4xgmpw[1] = _795e5iz[0]*_bk2l0d1[0] + _795e5iz[1]*_bk2l0d1[1] +  _795e5iz[2]*_bk2l0d1[2]; 
            
            (*_vwvm4w3)[_qke2gjr][9]  = MIN((*_vwvm4w3)[_qke2gjr][9],  _b4xgmpw[0]);
            (*_vwvm4w3)[_qke2gjr][10] = MAX((*_vwvm4w3)[_qke2gjr][10], _b4xgmpw[0]);
            (*_vwvm4w3)[_qke2gjr][11] = MIN((*_vwvm4w3)[_qke2gjr][11], _b4xgmpw[1]);
            (*_vwvm4w3)[_qke2gjr][12] = MAX((*_vwvm4w3)[_qke2gjr][12], _b4xgmpw[1]);
            
            _id4n1u3(_bk2l0d1, _fc4kl8c, _9p3jpxg);
            _b4xgmpw[0] = _j3swebn[0]*_bk2l0d1[0] + _j3swebn[1]*_bk2l0d1[1] +  _j3swebn[2]*_bk2l0d1[2]; 
            _b4xgmpw[1] = _795e5iz[0]*_bk2l0d1[0] + _795e5iz[1]*_bk2l0d1[1] +  _795e5iz[2]*_bk2l0d1[2]; 
            
            (*_vwvm4w3)[_qke2gjr][9]  = MIN((*_vwvm4w3)[_qke2gjr][9],  _b4xgmpw[0]);
            (*_vwvm4w3)[_qke2gjr][10] = MAX((*_vwvm4w3)[_qke2gjr][10], _b4xgmpw[0]);
            (*_vwvm4w3)[_qke2gjr][11] = MIN((*_vwvm4w3)[_qke2gjr][11], _b4xgmpw[1]);
            (*_vwvm4w3)[_qke2gjr][12] = MAX((*_vwvm4w3)[_qke2gjr][12], _b4xgmpw[1]);
    }   }
    
    for (i = 0u; i < _katqyds; i++)
    {   _dan085r[0] = (*_vwvm4w3)[i][9];              _dan085r[1] = (*_vwvm4w3)[i][11];          _dan085r[2] = 0.0f;
        _ntdr5xy[0] = (*_vwvm4w3)[i][10];             _ntdr5xy[1] = (*_vwvm4w3)[i][11];          _ntdr5xy[2] = 0.0f;
        _fc4kl8c[0] = (*_vwvm4w3)[i][10];             _fc4kl8c[1] = (*_vwvm4w3)[i][12];          _fc4kl8c[2] = 0.0f;
        _gn0t6dv[0] = (*_vwvm4w3)[i][9];              _gn0t6dv[1] = (*_vwvm4w3)[i][12];          _gn0t6dv[2] = 0.0f;

        _5gdej31 = fabsf(_ntdr5xy[0]-_dan085r[0])/fabsf(_fc4kl8c[1]-_ntdr5xy[1]); 
        _gymdm0l = sqrtf((*_vwvm4w3)[i][3]/_5gdej31); 
        _gwaj1g1 = (*_vwvm4w3)[i][3]/_gymdm0l; 

        _r1dwwo8 = fabsf(_ntdr5xy[0]-_dan085r[0])/_gwaj1g1;
        _nbhnal3 = fabsf(_fc4kl8c[1]-_ntdr5xy[1])/_gymdm0l;
        
        _dan085r[0] /= _r1dwwo8;              _dan085r[1] /= _nbhnal3;
        _ntdr5xy[0] /= _r1dwwo8;              _ntdr5xy[1] /= _nbhnal3;
        _fc4kl8c[0] /= _r1dwwo8;              _fc4kl8c[1] /= _nbhnal3;
        _gn0t6dv[0] /= _r1dwwo8;              _gn0t6dv[1] /= _nbhnal3;
        
        memcpy(_088o98b, &(*_vwvm4w3)[i][6], 3*sizeof(float));
        
        _ktakpux(_088o98b, _j3swebn, _795e5iz);

        _bk2l0d1[0] = _j3swebn[0]*_dan085r[0] + _795e5iz[0]*_dan085r[1] + _088o98b[0]*_dan085r[2] + (*_vwvm4w3)[i][0];
        _bk2l0d1[1] = _j3swebn[1]*_dan085r[0] + _795e5iz[1]*_dan085r[1] + _088o98b[1]*_dan085r[2] + (*_vwvm4w3)[i][1];
        _bk2l0d1[2] = _j3swebn[2]*_dan085r[0] + _795e5iz[2]*_dan085r[1] + _088o98b[2]*_dan085r[2] + (*_vwvm4w3)[i][2];
        _bk2l0d1[2] = MIN(_bk2l0d1[2], -0.0f);
        memcpy(&(*_vwvm4w3)[i][9], _bk2l0d1, 3*sizeof(float));

        _bk2l0d1[0] = _j3swebn[0]*_ntdr5xy[0] + _795e5iz[0]*_ntdr5xy[1] + _088o98b[0]*_ntdr5xy[2] + (*_vwvm4w3)[i][0];
        _bk2l0d1[1] = _j3swebn[1]*_ntdr5xy[0] + _795e5iz[1]*_ntdr5xy[1] + _088o98b[1]*_ntdr5xy[2] + (*_vwvm4w3)[i][1];
        _bk2l0d1[2] = _j3swebn[2]*_ntdr5xy[0] + _795e5iz[2]*_ntdr5xy[1] + _088o98b[2]*_ntdr5xy[2] + (*_vwvm4w3)[i][2];
        _bk2l0d1[2] = MIN(_bk2l0d1[2], -0.0f);
        memcpy(&(*_vwvm4w3)[i][12], _bk2l0d1, 3*sizeof(float));

        _bk2l0d1[0] = _j3swebn[0]*_fc4kl8c[0] + _795e5iz[0]*_fc4kl8c[1] + _088o98b[0]*_fc4kl8c[2] + (*_vwvm4w3)[i][0];
        _bk2l0d1[1] = _j3swebn[1]*_fc4kl8c[0] + _795e5iz[1]*_fc4kl8c[1] + _088o98b[1]*_fc4kl8c[2] + (*_vwvm4w3)[i][1];
        _bk2l0d1[2] = _j3swebn[2]*_fc4kl8c[0] + _795e5iz[2]*_fc4kl8c[1] + _088o98b[2]*_fc4kl8c[2] + (*_vwvm4w3)[i][2];
        _bk2l0d1[2] = MIN(_bk2l0d1[2], -0.0f);
        memcpy(&(*_vwvm4w3)[i][15], _bk2l0d1, 3*sizeof(float));

        _bk2l0d1[0] = _j3swebn[0]*_gn0t6dv[0] + _795e5iz[0]*_gn0t6dv[1] + _088o98b[0]*_gn0t6dv[2] + (*_vwvm4w3)[i][0];
        _bk2l0d1[1] = _j3swebn[1]*_gn0t6dv[0] + _795e5iz[1]*_gn0t6dv[1] + _088o98b[1]*_gn0t6dv[2] + (*_vwvm4w3)[i][1];
        _bk2l0d1[2] = _j3swebn[2]*_gn0t6dv[0] + _795e5iz[2]*_gn0t6dv[1] + _088o98b[2]*_gn0t6dv[2] + (*_vwvm4w3)[i][2];
        _bk2l0d1[2] = MIN(_bk2l0d1[2], -0.0f);
        memcpy(&(*_vwvm4w3)[i][18], _bk2l0d1, 3*sizeof(float));

        _id4n1u3(_n42x2ez,  &(*_vwvm4w3)[i][0], &(*_vwvm4w3)[i][9]);                _31yrujo  = _mabsi9n(_n42x2ez);
        (*_vwvm4w3)[i][5] =  MAX(_31yrujo, (*_vwvm4w3)[i][5]);
        _id4n1u3(_n42x2ez,  &(*_vwvm4w3)[i][0], &(*_vwvm4w3)[i][12]);               _31yrujo  = _mabsi9n(_n42x2ez);
        (*_vwvm4w3)[i][5] =  MAX(_31yrujo, (*_vwvm4w3)[i][5]);
        _id4n1u3(_n42x2ez,  &(*_vwvm4w3)[i][0], &(*_vwvm4w3)[i][15]);               _31yrujo  = _mabsi9n(_n42x2ez);
        (*_vwvm4w3)[i][5] =  MAX(_31yrujo, (*_vwvm4w3)[i][5]);
        _id4n1u3(_n42x2ez,  &(*_vwvm4w3)[i][0], &(*_vwvm4w3)[i][18]);               _31yrujo  = _mabsi9n(_n42x2ez);
        (*_vwvm4w3)[i][5] =  MAX(_31yrujo, (*_vwvm4w3)[i][5]);
    }
    
    return;
}



void _e4q6qdy(const char *restrict _dygflij, unsigned int _jz6us7r, const float *restrict _jdtvv2y, unsigned int **restrict _1xgbngz, unsigned int _gzvrctl, float **restrict _4pa8hye, unsigned int _ntbz3m5, const float *restrict _8e7flqh, unsigned int **restrict _f9yecaz, unsigned int _2vpgstq, float **restrict _32orwix)
{   
    int i, _5ccasbu  = 16; 
    char _t7zewvi[512];
    float *_zareiya = NULL,   *_wd1xmxb = NULL;
    
    strcpy(_t7zewvi,_dygflij);   strcat(_t7zewvi,".brch");
    FILE *_yk1b3d5;
    
    _zareiya     = malloc(3*_jz6us7r*sizeof(float));
    for (i = 0; i < _jz6us7r; i++)      {       memcpy(&_zareiya[i*3 +0], &_jdtvv2y[i*16 +0],3*sizeof(float));     }
    
    if (_ntbz3m5 > 0u)
    {   _wd1xmxb = malloc(3*_ntbz3m5*sizeof(float));
        for (i = 0; i < _ntbz3m5; i++)      {       memcpy(&_wd1xmxb[i*3 +0], &_8e7flqh[i*16 +0],3*sizeof(float));     }
    }
    
    if ((_yk1b3d5 = fopen(_t7zewvi,"wb")) == NULL)     {    printf("Error -cant open BranchInfo file for writing   %s\n", _t7zewvi);      exit(10);     }
    
    int _3zezqac = (int)VersionNumber;
    fwrite( &_3zezqac,          sizeof(int),   1, _yk1b3d5);
    fwrite( &_jz6us7r,            sizeof(unsigned int),   1, _yk1b3d5); 
    fwrite( &_gzvrctl,          sizeof(unsigned int),   1, _yk1b3d5); 
    for (i = 0; i < _jz6us7r; i++)      {   fwrite(_1xgbngz[i],        sizeof(unsigned int),  _5ccasbu, _yk1b3d5);   }
    for (i = 0; i < _gzvrctl; i++)    {   fwrite(_4pa8hye[i], sizeof(float),         21,        _yk1b3d5);   }
    fwrite(_zareiya,         sizeof(float), 3*_jz6us7r, _yk1b3d5);
    
    if (_ntbz3m5 > 0u)
    {   fwrite( &_ntbz3m5,        sizeof(unsigned int),   1, _yk1b3d5); 
        fwrite( &_2vpgstq,      sizeof(unsigned int),   1, _yk1b3d5); 
        for (i = 0; i < _ntbz3m5; i++)      {   fwrite(_f9yecaz[i],        sizeof(unsigned int),  _5ccasbu, _yk1b3d5);   }
        for (i = 0; i < _2vpgstq; i++)    {   fwrite(_32orwix[i], sizeof(float),         21,        _yk1b3d5);   }
        fwrite(_wd1xmxb,         sizeof(float), 3*_ntbz3m5, _yk1b3d5);
    }
    
    fclose(_yk1b3d5);
    free(_zareiya);   free(_wd1xmxb);
    
    return;
}