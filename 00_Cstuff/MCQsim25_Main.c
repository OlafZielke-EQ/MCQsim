#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <mpi.h>
#include <limits.h>
#include <time.h>
#include <cblas.h>
#define VersionNumber 202511

#define ALIGN64 64 
#define ADDEDSTRGTHVALUE -50000.0f 

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


extern void  _148p0tv(int _fsqg8x8, const char *restrict _ut076g5, char *restrict _g79or78, char *restrict _dygflij,  struct _5r5uqpg *restrict _atj9tze , struct _xjuy0bv *restrict _vxpj4hn, struct _j49itn7 *restrict _s55j2o7);
extern void     _fqewaou(int _fsqg8x8, const char *restrict _dygflij, const struct _xjuy0bv *restrict _vxpj4hn, struct _j49itn7 *restrict _s55j2o7, const int *restrict _lexdw6l, const unsigned int *restrict _1nlc8jp, unsigned int *restrict _vmefozt,  unsigned int *restrict _rgb93wk, float *restrict _4xkss0k, float *restrict _rrpzlk3, float *restrict _zzf8s2z, float *restrict _3n3bzt9, float *restrict _0dl5e11);
extern void    _qjiujv6(int _fsqg8x8, struct _5r5uqpg *restrict _atj9tze, struct _xjuy0bv *restrict _vxpj4hn, struct _j49itn7 *restrict _s55j2o7, const int *restrict _lexdw6l, const int *restrict _1poxrte, const unsigned int *restrict _1nlc8jp, const unsigned int *restrict _fftouhy, float *restrict _e9ftaxp, float *restrict _9yo18pa, unsigned int *restrict _vmefozt, unsigned int *restrict _rgb93wk, float *restrict _4xkss0k, float *restrict _rrpzlk3, const float *restrict _zzf8s2z, const float *restrict _3n3bzt9, float *restrict _6fokzix, float *restrict _0dl5e11, float *restrict _z39pecn, float *restrict _7yt2m6t);

extern void  _43d5iqv(int _0rmhrcn, int _fsqg8x8, const char *restrict _dygflij, struct _5r5uqpg *_atj9tze, struct _xjuy0bv *restrict _vxpj4hn, const struct _j49itn7 *restrict _s55j2o7, const int *restrict _lexdw6l, const int *restrict _1poxrte, const unsigned int *restrict _1nlc8jp, const unsigned int *restrict _fftouhy, const float *restrict _e9ftaxp, const float *restrict _9yo18pa, const unsigned int *restrict _vmefozt, const unsigned int *restrict _rgb93wk, const float *restrict _4xkss0k, const float *restrict _rrpzlk3, const float *restrict _zzf8s2z, const float *restrict _3n3bzt9, unsigned int *restrict _p422ta8, float *restrict _6fokzix, float *restrict _0dl5e11, float *restrict _7yt2m6t, unsigned short int *restrict _qnsn2s8, unsigned short int *restrict _qpwbsf3, unsigned short int *restrict _iymyrsv, unsigned short int *restrict _iob0nlw, unsigned short int ***restrict _5gzpl1q, unsigned short int ***restrict _oqgvzhq, unsigned short int ***restrict _1ajo091, unsigned short int ***restrict _3b9so0u, float ***restrict _ckeh5sw, float ***restrict _3ac6pe6, float ***restrict _obyghgk, float ***restrict _yj9boq0, float ***restrict _ljefy3z, float ***restrict _h61iaoi, float ***restrict _plozjcn, float ***restrict _qk4ltah, float ***restrict _0ljch2v, float ***restrict _0dtirc9, float ***restrict _9dde3cu);
extern void     _82n2atq(int _fsqg8x8, const char *restrict _g79or78, const struct _xjuy0bv *restrict _vxpj4hn, const struct _j49itn7 *restrict _s55j2o7, const int *restrict _paih0p4, const int *restrict _20n8tte, const int *restrict _lexdw6l, const int *restrict _1poxrte, const unsigned int *restrict _1nlc8jp, const unsigned int *restrict _fftouhy, const float *restrict _6fokzix, const float *restrict _7yt2m6t, unsigned short int *restrict _qnsn2s8, unsigned short int *restrict _qpwbsf3, unsigned short int *restrict _iymyrsv, unsigned short int *restrict _iob0nlw, unsigned short int **restrict _5gzpl1q, unsigned short int **restrict _oqgvzhq, unsigned short int **restrict _1ajo091, unsigned short int **restrict _3b9so0u, float **restrict _ckeh5sw, float **restrict _3ac6pe6, float **restrict _obyghgk, float **restrict _yj9boq0, float **restrict _ljefy3z, float **restrict _h61iaoi, float **restrict _plozjcn, float **restrict _qk4ltah, float **restrict _0ljch2v, float **restrict _0dtirc9, float **restrict _9dde3cu);
extern void     _l1n13jc(int _0rmhrcn, int _fsqg8x8, const char *restrict _g79or78, struct _5r5uqpg *restrict _atj9tze, struct _xjuy0bv *restrict _vxpj4hn, const struct _j49itn7 *restrict _s55j2o7, const int *restrict _paih0p4, const int *restrict _20n8tte, const int *restrict _lexdw6l, const int *restrict _1poxrte, const unsigned int *restrict _1nlc8jp, unsigned int *restrict _p422ta8, float *restrict _6fokzix, float *restrict _0dl5e11, const unsigned int *restrict _fftouhy, float *restrict _7yt2m6t, unsigned short int *restrict _qnsn2s8, unsigned short int *restrict _qpwbsf3, unsigned short int *restrict _iymyrsv, unsigned short int *restrict _iob0nlw, unsigned short int ***restrict _5gzpl1q, unsigned short int ***restrict _oqgvzhq, unsigned short int ***restrict _1ajo091, unsigned short int ***restrict _3b9so0u, float ***restrict _ckeh5sw, float ***restrict _3ac6pe6, float ***restrict _obyghgk, float ***restrict _yj9boq0, float ***restrict _ljefy3z, float ***restrict _h61iaoi, float ***restrict _plozjcn, float ***restrict _qk4ltah, float ***restrict _0ljch2v, float ***restrict _0dtirc9, float ***restrict _9dde3cu);

extern void   _rpj1epk(int _0rmhrcn, int _fsqg8x8, const struct _5r5uqpg *restrict _atj9tze, const struct _xjuy0bv *restrict _vxpj4hn, const int *restrict _lexdw6l, const int *restrict _1poxrte, const unsigned int *restrict _1nlc8jp, const unsigned int *restrict _fftouhy, float *restrict _6fokzix, float *restrict _z39pecn, const float *restrict _7yt2m6t, float *restrict _rrpzlk3, const unsigned short int *restrict _qnsn2s8, const unsigned short int *restrict _qpwbsf3, const unsigned short int *restrict _iymyrsv, const unsigned short int *restrict _iob0nlw, const unsigned short int **restrict _5gzpl1q, const unsigned short int **restrict _oqgvzhq, const unsigned short int **restrict _1ajo091, const unsigned short int **restrict _3b9so0u, const float **restrict _ckeh5sw, const float **restrict _3ac6pe6, const float **restrict _obyghgk, const float **restrict _yj9boq0, const float **restrict _ljefy3z, const float **restrict _h61iaoi, const float **restrict _plozjcn, const float **restrict _qk4ltah, const float **restrict _0ljch2v, const float **restrict _0dtirc9, const float **restrict _9dde3cu);

extern long int  _z0yncxl(int _uuijnwp, int _qtkcd0r, int _tgxnx4s);
extern void    _id4n1u3(float *restrict _j49yd82, const float *restrict _28xg8y7, const float *restrict _kf1622v);
extern void     _j2qp89u(float *restrict _iquwrxz, const float *restrict _f445cri, const float *restrict _8bixoiq);
extern float    _mabsi9n(const float *restrict _n42x2ez);
extern void       _bfuk9e1(float *restrict _n42x2ez);
extern void _ktakpux(const float *restrict _088o98b, float *restrict _j3swebn, float *restrict _795e5iz);
extern void     _ftu2gc5(float *restrict _0i9g5dz, const float *restrict _mt8vt43, const float *restrict _xyzx9jp);
extern void     _yyv37di(float *restrict _0yknfu8, float _ktmonrv);
extern float            _53ayeqa(long *_rrhf4ab);
extern float            _hthsjr0(long *_rrhf4ab);

extern void _myx444n(struct _5r5uqpg *restrict _atj9tze, const struct _xjuy0bv *restrict _vxpj4hn, unsigned int i, unsigned int *restrict _p422ta8, float *restrict _6fokzix, float *restrict _0dl5e11, unsigned int _urvqt3g);
extern void  _790pdw5(int _fsqg8x8, const int *restrict _lexdw6l, struct _5r5uqpg *restrict _atj9tze, unsigned int *restrict _p422ta8, float *restrict _6fokzix, float *restrict _0dl5e11, float *restrict _z39pecn);
extern void     _4lmct64(int _fsqg8x8, int _aytxw1v, int _zc7d2kd, const char *restrict _zle6kwh, unsigned int _6zjq7w7, const int *restrict _lexdw6l, const unsigned int *restrict _1nlc8jp, const unsigned int *restrict _vmefozt, const float *restrict _4xkss0k, const float *restrict _zzf8s2z, const float *restrict _0dl5e11, const float *restrict _6fokzix, const float *restrict _z39pecn);
extern void  _z7w2seh(int _fsqg8x8, const char *restrict _zle6kwh, const struct _xjuy0bv *restrict _vxpj4hn, const struct _j49itn7 *restrict _s55j2o7, const int *restrict _lexdw6l, const unsigned int *restrict _1nlc8jp, const unsigned int *restrict _vmefozt, const float *restrict _4xkss0k, const float *restrict _zzf8s2z, const unsigned int *restrict _p422ta8, float *restrict _6fokzix, const float *restrict _0dl5e11);
extern void  _qpki7p5(int _fsqg8x8,  const char *restrict _dygflij, struct _5r5uqpg *restrict _atj9tze,  const struct _xjuy0bv *restrict _vxpj4hn, const struct _j49itn7 *restrict _s55j2o7, unsigned int *restrict _xo100vh, double *restrict _ofn1r4k, const int *restrict _lexdw6l, const int *restrict _1poxrte, const unsigned int *restrict _1nlc8jp, const unsigned int *restrict _fftouhy, const unsigned int *restrict _vmefozt, const float *restrict _4xkss0k, const float *restrict _zzf8s2z, unsigned int *restrict _p422ta8, float *restrict _6fokzix, float *restrict _0dl5e11, float *restrict _z39pecn, float *restrict _7yt2m6t);




int main(int _wk8r2yk, char **_ccprhvl)
{   if ( (_wk8r2yk  > 3 ) || (_wk8r2yk < 2) ) {   fprintf(stdout,"Input Error\n Please start the code in the following way:\n\n mpirun -np 4 ./MCQsim25 RunParaFile.txt");   exit(EXIT_FAILURE);    }
    
    int _fsqg8x8,   _0rmhrcn;

    MPI_Init(&_wk8r2yk, &_ccprhvl);
    MPI_Comm_rank(MPI_COMM_WORLD, &_fsqg8x8);
    MPI_Comm_size(MPI_COMM_WORLD, &_0rmhrcn);
    
    struct _5r5uqpg  __attribute__((aligned(ALIGN64))) _atj9tze;       memset(&_atj9tze, 0, sizeof(_atj9tze));
    struct _xjuy0bv    __attribute__((aligned(ALIGN64))) _vxpj4hn;       memset(&_vxpj4hn, 0, sizeof(_vxpj4hn));
    struct _j49itn7   __attribute__((aligned(ALIGN64))) _s55j2o7;       memset(&_s55j2o7, 0, sizeof(_s55j2o7));
    
    int            *_e05qmkt      = NULL,    *_avt8k06 = NULL,  *_qwlvunl = NULL,  *_o2t8fy6 = NULL;
    unsigned int   *_r7vgrf1    = NULL,  *_lu6k7p9 = NULL;
    const int      *_paih0p4       = NULL,    *_20n8tte  = NULL,  *_lexdw6l  = NULL,  *_1poxrte  = NULL;
    const unsigned int *_1nlc8jp = NULL,  *_fftouhy  = NULL;
    unsigned int   *_vmefozt          = NULL,        *_rgb93wk = NULL,       *_p422ta8 = NULL;
    float          *_4xkss0k          = NULL,       *_zzf8s2z = NULL,       *_rrpzlk3 = NULL,      *_3n3bzt9 = NULL,       *_e9ftaxp     = NULL,      *_9yo18pa = NULL;
    float          *_0dl5e11          = NULL,        *_6fokzix = NULL,       *_7yt2m6t = NULL,     *_z39pecn = NULL;

    unsigned short int      *_nn16uxf = NULL,     *_q0e27ze = NULL,    *_zn4regb = NULL,    *_0vs52sq = NULL;
    unsigned short int     **_weehwjq  = NULL,    **_jy62zee  = NULL,   **_c795rlq  = NULL,   **_swbkok3  = NULL;
    float                  **_2292mtr = NULL,    **_iq2fa39 = NULL,   **_m458ax4 = NULL,   **_ffqot8h = NULL;
    float                  **_3lbqgu3 = NULL,    **_b4ldims = NULL,   **_s6g298c = NULL,   **_zd76szn = NULL;
    float                  **_5bwikot = NULL,                            **_ukwqvr1 = NULL,   **_p8ijl2z = NULL;
    const unsigned short int  *_qnsn2s8 = NULL,     *_qpwbsf3   = NULL,    *_iymyrsv   = NULL,    *_iob0nlw   = NULL;
    const unsigned short int **_5gzpl1q  = NULL,    **_oqgvzhq    = NULL,   **_1ajo091    = NULL,   **_3b9so0u    = NULL;
    const float              **_ckeh5sw = NULL,    **_yj9boq0   = NULL,   **_h61iaoi   = NULL,   **_0ljch2v   = NULL;
    const float              **_3ac6pe6 = NULL,    **_ljefy3z   = NULL,   **_plozjcn   = NULL,   **_0dtirc9   = NULL;
    const float              **_obyghgk = NULL,                            **_qk4ltah   = NULL,   **_9dde3cu   = NULL;

    float       **_66cnolj = NULL;
    unsigned int *_pj8t22c  = NULL,   *_50b8y5d = NULL,   *_zk94yxf = NULL,   *_iblejgu = NULL,   *_hx6zdcu = NULL;
    float        *_dtg5nqj= NULL,  *_kc4lorx = NULL,    *_pkxvhag = NULL,     *_qijc1yd = NULL,   *_2ljk4q5= NULL,   *_v7c6s5n = NULL,   *_kb4ehou = NULL,   *_6qlz2rw = NULL,   *_d1vn8lf = NULL,   *_1vsmado = NULL;
    float         *_ozctwjl = NULL,    *_ny181ac = NULL,       *_t997e4u = NULL,      *_lqbmae3 = NULL,     *_gyspeng = NULL,      *_zkwdeek = NULL;
    const float  *_jbc2sws = NULL,   *_wyx9trf = NULL,      *_sh3tow4 = NULL,      *_8sryo4r = NULL;

    long int _moxsfh4;
    
    char _ut076g5[512],   _g79or78[512],   _dygflij[512];
    strcpy(_ut076g5,  _ccprhvl[1]);
    
    _148p0tv(_fsqg8x8, _ut076g5, _g79or78, _dygflij, &_atj9tze, &_vxpj4hn, &_s55j2o7);
    
    {   unsigned int i;
        int  _jq6fx7x = (int)(_vxpj4hn._3ru7myo/_0rmhrcn); 
        int  _jp61lfc  = (int)(_vxpj4hn._3ru7myo%_0rmhrcn);
        int  _mda9oek = (int)(_vxpj4hn._5izolnn/_0rmhrcn); 
        int  _xd356wa  = (int)(_vxpj4hn._5izolnn%_0rmhrcn);

        _moxsfh4   = _z0yncxl(ALIGN64, _0rmhrcn, sizeof(int));
        _e05qmkt    = aligned_alloc(ALIGN64, _moxsfh4);
        _avt8k06    = aligned_alloc(ALIGN64, _moxsfh4);
        _qwlvunl   = aligned_alloc(ALIGN64, _moxsfh4);
        _o2t8fy6   = aligned_alloc(ALIGN64, _moxsfh4);
        _e05qmkt[0] = 0;
        _avt8k06[0] = 0;
        
        for (i = 0u; i < _0rmhrcn;     i++)      {   _qwlvunl[i]     = _jq6fx7x;                            }
        for (i = 0u; i < _jp61lfc; i++)      {   _qwlvunl[i]    += 1;                                     }
        for (i = 1u; i < _0rmhrcn;     i++)      {   _e05qmkt[i]      = _e05qmkt[i-1] + _qwlvunl[i-1];        }

        for (i = 0u; i < _0rmhrcn;     i++)      {   _o2t8fy6[i]     = _mda9oek;                            }
        for (i = 0u; i < _xd356wa; i++)      {   _o2t8fy6[i]    += 1;                                     }
        for (i = 1u; i < _0rmhrcn;     i++)      {   _avt8k06[i]      = _avt8k06[i-1] + _o2t8fy6[i-1];        }
        
        _moxsfh4   = _z0yncxl(ALIGN64, _qwlvunl[_fsqg8x8], sizeof(unsigned int));
        _r7vgrf1  = aligned_alloc(ALIGN64, _moxsfh4);
        for (i = 0u; i < _qwlvunl[_fsqg8x8]; i++) {   _r7vgrf1[i] = i*_0rmhrcn +_fsqg8x8;                         }
        _paih0p4 = (const int *)_e05qmkt;            _lexdw6l = (const int *)_qwlvunl;          _1nlc8jp = (const unsigned int *)_r7vgrf1;
        
        if (_vxpj4hn._5izolnn > 0u)
        {   _moxsfh4   = _z0yncxl(ALIGN64, _o2t8fy6[_fsqg8x8], sizeof(unsigned int));
            _lu6k7p9  = aligned_alloc(ALIGN64, _moxsfh4);
            for (i = 0u; i < _o2t8fy6[_fsqg8x8]; i++) {   _lu6k7p9[i] = i*_0rmhrcn +_fsqg8x8;                          }
        }
        _20n8tte = (const int *)_avt8k06;            _1poxrte = (const int *)_o2t8fy6;          _fftouhy = (const unsigned int *)_lu6k7p9;
        
        _moxsfh4   = _z0yncxl(ALIGN64, 8*_vxpj4hn._3ru7myo,  sizeof(unsigned int));            _vmefozt        = aligned_alloc(ALIGN64, _moxsfh4);                memset(_vmefozt,  0, _moxsfh4);  
        _moxsfh4   = _z0yncxl(ALIGN64, 16*_vxpj4hn._3ru7myo, sizeof(float));                   _4xkss0k        = aligned_alloc(ALIGN64, _moxsfh4);                memset(_4xkss0k,  0, _moxsfh4);  
        _moxsfh4   = _z0yncxl(ALIGN64, 4*_s55j2o7._fu7bfxq,  sizeof(float));                   _zzf8s2z       = aligned_alloc(ALIGN64, _moxsfh4);                memset(_zzf8s2z, 0, _moxsfh4); 
        
        _moxsfh4   = _z0yncxl(ALIGN64, 16*_lexdw6l[_fsqg8x8], sizeof(float));                  _0dl5e11        = aligned_alloc(ALIGN64, _moxsfh4);                memset(_0dl5e11,   0, _moxsfh4);
        _moxsfh4   = _z0yncxl(ALIGN64, 8*_lexdw6l[_fsqg8x8], sizeof(float));                   _z39pecn      = aligned_alloc(ALIGN64, _moxsfh4);                memset(_z39pecn, 0, _moxsfh4);
        _moxsfh4   = _z0yncxl(ALIGN64, 4*_lexdw6l[_fsqg8x8], sizeof(unsigned int));            _p422ta8        = aligned_alloc(ALIGN64, _moxsfh4);                memset(_p422ta8,   0, _moxsfh4);
        _moxsfh4   = _z0yncxl(ALIGN64, 16*_lexdw6l[_fsqg8x8], sizeof(float));                  _6fokzix        = aligned_alloc(ALIGN64, _moxsfh4);                memset(_6fokzix,   0, _moxsfh4);
        
        _moxsfh4   = _z0yncxl(ALIGN64, _lexdw6l[_fsqg8x8], sizeof(unsigned short int));        _nn16uxf = aligned_alloc(ALIGN64, _moxsfh4);                memset(_nn16uxf, 0, _moxsfh4);
        
        if (_vxpj4hn._5izolnn > 0u)
        {   
            _moxsfh4   = _z0yncxl(ALIGN64, 8*_vxpj4hn._5izolnn, sizeof(unsigned int));         _rgb93wk        = aligned_alloc(ALIGN64, _moxsfh4);                memset(_rgb93wk, 0, _moxsfh4);   
            _moxsfh4   = _z0yncxl(ALIGN64, 16*_vxpj4hn._5izolnn, sizeof(float));               _rrpzlk3        = aligned_alloc(ALIGN64, _moxsfh4);                memset(_rrpzlk3, 0, _moxsfh4);  
            _moxsfh4   = _z0yncxl(ALIGN64, 4*_s55j2o7._ryf0992,  sizeof(float));               _3n3bzt9       = aligned_alloc(ALIGN64, _moxsfh4);                memset(_3n3bzt9,0, _moxsfh4);  
            
            _moxsfh4   = _z0yncxl(ALIGN64, 16*_1poxrte[_fsqg8x8], sizeof(float));              _7yt2m6t        = aligned_alloc(ALIGN64, _moxsfh4);                memset(_7yt2m6t,   0, _moxsfh4);
            
            _moxsfh4   = _z0yncxl(ALIGN64, _lexdw6l[_fsqg8x8], sizeof(unsigned short int));    _q0e27ze = aligned_alloc(ALIGN64, _moxsfh4);                memset(_q0e27ze, 0, _moxsfh4);
            _moxsfh4   = _z0yncxl(ALIGN64, _1poxrte[_fsqg8x8], sizeof(unsigned short int));    _zn4regb = aligned_alloc(ALIGN64, _moxsfh4);                memset(_zn4regb, 0, _moxsfh4);
                                                                                                _0vs52sq = aligned_alloc(ALIGN64, _moxsfh4);                memset(_0vs52sq, 0, _moxsfh4);
            
        }
    } 
    
    _fqewaou(_fsqg8x8, _dygflij, &_vxpj4hn, &_s55j2o7, _lexdw6l, _1nlc8jp, _vmefozt, _rgb93wk, _4xkss0k, _rrpzlk3, _zzf8s2z, _3n3bzt9, _0dl5e11);
    
    {   _moxsfh4 = _z0yncxl(ALIGN64, 8*_s55j2o7._6wmdw0f, sizeof(float));                      _e9ftaxp     = (float *) aligned_alloc(ALIGN64, _moxsfh4);        memset(_e9ftaxp, 0, _moxsfh4);
        
        if (_vxpj4hn._5izolnn > 0u)
        {   _moxsfh4 = _z0yncxl(ALIGN64, 8*_s55j2o7._rp9dxoq, sizeof(float));                  _9yo18pa     = (float *) aligned_alloc(ALIGN64, _moxsfh4);        memset(_9yo18pa, 0, _moxsfh4);
        }
    } 
    
    _qjiujv6(_fsqg8x8, &_atj9tze, &_vxpj4hn, &_s55j2o7, _lexdw6l, _1poxrte, _1nlc8jp, _fftouhy, _e9ftaxp, _9yo18pa, _vmefozt, _rgb93wk, _4xkss0k, _rrpzlk3, _zzf8s2z, _3n3bzt9, _6fokzix, _0dl5e11, _z39pecn, _7yt2m6t);
    
    if (_fsqg8x8 == 0) 
    {   float _nbzug0i; 
        fprintf(stdout,"\n\n\n------------------------------------------------------------------\n");
        fprintf(stdout,"MCQsim version number:  %u\n",VersionNumber);
        fprintf(stdout,"----------------------\n");
        fprintf(stdout,"Number of RANKS: %u\n",_0rmhrcn);
        fprintf(stdout,"----------------------\n");
        fprintf(stdout,"System info: Byte Size for FLOAT: %lu     unsigned INT: %lu    long INT: %lu\n", sizeof(float), sizeof(unsigned int), sizeof(long int));
        fprintf(stdout,"----------------------\n");
        fprintf(stdout,"FileName:               %s\n",_dygflij);                        fprintf(stdout,"RunNumber:              %u\n",_s55j2o7._jcmrsez);
        fprintf(stdout,"PlotCat2Screen:         %u\n",_vxpj4hn._qdbxvh7);
        fprintf(stdout,"MinElem4Catalog:        %u\n",_vxpj4hn._e6zvlhy);                fprintf(stdout,"ContPreviousRun:        %u\n",_s55j2o7._qy8yb7w);
        fprintf(stdout,"LoadPrev_Kh_mat:        %u\n",_s55j2o7._o7vbgr8);                 fprintf(stdout,"Kh_mat file name:       %s\n",_g79or78);
        fprintf(stdout,"----------------------\n");
        fprintf(stdout,"SaveSTF4LargeEQ:        %u\n",_vxpj4hn._wk8temp);               fprintf(stdout,"MinMag2SaveSTF:         %f\n",_vxpj4hn._rdihsa7); 
        fprintf(stdout,"UseHalfSpace:           %u\n",_s55j2o7._s832ejx);
        fprintf(stdout,"----------------------\n"); 
        fprintf(stdout,"ViscAftSlipTime:        %fyears\n",_vxpj4hn._lw2lm34);

        fprintf(stdout,"----------------------\n"); 
        fprintf(stdout,"CoSeisHealFract:        %f\n",_vxpj4hn._3jczucq);               fprintf(stdout,"OvershootFract:         %f\n",_vxpj4hn._x5qjicj);
        fprintf(stdout,"PreStressFract:         %f\n",_s55j2o7._dhgzz0j);           fprintf(stdout,"MinCoSeisSlipRate(m/s): %f\n",(_vxpj4hn._7odkef8/_atj9tze._mzwffhh));
        fprintf(stdout,"----------------------\n"); 
        fprintf(stdout,"RecLength:              %lf\n",_atj9tze._0ur5l8a);
        fprintf(stdout,"----------------------\n");
        fprintf(stdout,"Elastic properties\n");
        fprintf(stdout,"Medium density (kg/m^3):   %e\n",_atj9tze._cmktlzk);                fprintf(stdout,"AddedNormalStress (MPa):   %e\n",_atj9tze._iwvhdk8);
        fprintf(stdout,"ShearModulus (Pa):         %e\n",_atj9tze._kkrzd5e);              fprintf(stdout,"PoissonRatio:              %e\n",_atj9tze._wrb7f6f);
        fprintf(stdout,"ChangeFricBtwEQs:          %u\n",_vxpj4hn._pwh0mx1);             fprintf(stdout,"Lambda (Pa):               %e\n",_atj9tze._6snglh6);
        fprintf(stdout,"S-waveVelocity (m/s):      %e\n",_atj9tze._douadvw);
        fprintf(stdout,"unit slip (m, fault):      %e\n",_s55j2o7._ycz04cd);

        fprintf(stdout,"min. slip 2 start EQ (m):  %e\n",_vxpj4hn._7odkef8);           fprintf(stdout,"coseis. timestep (s):      %e\n",_atj9tze._mzwffhh);
        fprintf(stdout,"------------------------------------------------------------------\n\n");
        _nbzug0i = 100.0f - expf(-1.0f/_vxpj4hn._lw2lm34)*100.0f; 
        fprintf(stdout,"\nFractional post-seismic change during first year of after-slip:      %2.3f%% released \n",_nbzug0i);


        



       fprintf(stdout,"FaultCent2EdgeLgth and ElementArea: %5.1fm  and %2.4fkm^2\n\n",   _s55j2o7._b68rerv, _s55j2o7._xf7jxjn*1.0E-6f);

    } 
    
    int _vsffcxh = 0;
    if (_vsffcxh == 1)
    {   char _zle6kwh[] = "ArrayEntry2Check.txt";
        int  _aytxw1v = 4; 
        int  _zc7d2kd  = 8; 
        _4lmct64(_fsqg8x8, _aytxw1v, _zc7d2kd, _zle6kwh,  _vxpj4hn._3ru7myo,  _lexdw6l, _1nlc8jp,  _vmefozt, _4xkss0k, _zzf8s2z, _0dl5e11, _6fokzix, _z39pecn);
    } 
    
    
    if      (_s55j2o7._o7vbgr8 == 0) 
    {   
        _43d5iqv(_0rmhrcn, _fsqg8x8, _dygflij, &_atj9tze, &_vxpj4hn, &_s55j2o7, _lexdw6l, _1poxrte, _1nlc8jp, _fftouhy, _e9ftaxp, _9yo18pa, _vmefozt, _rgb93wk, _4xkss0k, _rrpzlk3, _zzf8s2z, _3n3bzt9, _p422ta8, _6fokzix, _0dl5e11, _7yt2m6t, _nn16uxf, _q0e27ze, _zn4regb, _0vs52sq, &_weehwjq, &_jy62zee, &_c795rlq, &_swbkok3, &_2292mtr, &_3lbqgu3, &_5bwikot, &_iq2fa39, &_b4ldims, &_m458ax4, &_s6g298c, &_ukwqvr1, &_ffqot8h, &_zd76szn, &_p8ijl2z);
    }
    else if (_s55j2o7._o7vbgr8 == 1) 
    {   
        _43d5iqv(_0rmhrcn, _fsqg8x8, _dygflij, &_atj9tze, &_vxpj4hn, &_s55j2o7, _lexdw6l, _1poxrte, _1nlc8jp, _fftouhy, _e9ftaxp, _9yo18pa, _vmefozt, _rgb93wk, _4xkss0k, _rrpzlk3, _zzf8s2z, _3n3bzt9, _p422ta8, _6fokzix, _0dl5e11, _7yt2m6t, _nn16uxf, _q0e27ze, _zn4regb, _0vs52sq, &_weehwjq, &_jy62zee, &_c795rlq, &_swbkok3, &_2292mtr, &_3lbqgu3, &_5bwikot, &_iq2fa39, &_b4ldims, &_m458ax4, &_s6g298c, &_ukwqvr1, &_ffqot8h, &_zd76szn, &_p8ijl2z);
        _82n2atq(_fsqg8x8, _g79or78, &_vxpj4hn, &_s55j2o7, _paih0p4, _20n8tte, _lexdw6l, _1poxrte, _1nlc8jp, _fftouhy, _6fokzix, _7yt2m6t, _nn16uxf, _q0e27ze, _zn4regb, _0vs52sq, _weehwjq, _jy62zee, _c795rlq, _swbkok3, _2292mtr, _3lbqgu3, _5bwikot, _iq2fa39, _b4ldims, _m458ax4, _s6g298c, _ukwqvr1, _ffqot8h, _zd76szn, _p8ijl2z);
    }
    else if (_s55j2o7._o7vbgr8 == 2) 
    {    _l1n13jc(_0rmhrcn, _fsqg8x8, _g79or78, &_atj9tze, &_vxpj4hn, &_s55j2o7, _paih0p4, _20n8tte, _lexdw6l, _1poxrte, _1nlc8jp, _p422ta8, _6fokzix, _0dl5e11, _fftouhy, _7yt2m6t, _nn16uxf, _q0e27ze, _zn4regb, _0vs52sq, &_weehwjq, &_jy62zee, &_c795rlq, &_swbkok3, &_2292mtr, &_3lbqgu3, &_5bwikot, &_iq2fa39, &_b4ldims, &_m458ax4, &_s6g298c, &_ukwqvr1, &_ffqot8h, &_zd76szn, &_p8ijl2z);
    }
    else if (_s55j2o7._o7vbgr8 == 3) 
    {   
        
    }
    _qnsn2s8 = (const unsigned short int  *)_nn16uxf;           _qpwbsf3 = (const unsigned short int  *)_q0e27ze;           _iymyrsv = (const unsigned short int  *)_zn4regb;           _iob0nlw = (const unsigned short int  *)_0vs52sq;
    _5gzpl1q  = (const unsigned short int **)_weehwjq;            _oqgvzhq  = (const unsigned short int **)_jy62zee;            _1ajo091  = (const unsigned short int **)_c795rlq;            _3b9so0u  = (const unsigned short int **)_swbkok3;
    _ckeh5sw = (const float              **)_2292mtr;           _3ac6pe6 = (const float              **)_3lbqgu3;           _obyghgk = (const float              **)_5bwikot;
    _yj9boq0 = (const float              **)_iq2fa39;           _ljefy3z = (const float              **)_b4ldims;
    _h61iaoi = (const float              **)_m458ax4;           _plozjcn = (const float              **)_s6g298c;           _qk4ltah = (const float              **)_ukwqvr1;
    _0ljch2v = (const float              **)_ffqot8h;           _0dtirc9 = (const float              **)_zd76szn;           _9dde3cu = (const float              **)_p8ijl2z;
    
    
    _rpj1epk(_0rmhrcn, _fsqg8x8, &_atj9tze, &_vxpj4hn, _lexdw6l, _1poxrte, _1nlc8jp, _fftouhy, _6fokzix, _z39pecn, _7yt2m6t, _rrpzlk3, _qnsn2s8, _qpwbsf3, _iymyrsv, _iob0nlw, _5gzpl1q, _oqgvzhq, _1ajo091, _3b9so0u, _ckeh5sw, _3ac6pe6, _obyghgk, _yj9boq0, _ljefy3z, _h61iaoi, _plozjcn, _qk4ltah, _0ljch2v, _0dtirc9, _9dde3cu);
    _z7w2seh(_fsqg8x8,_dygflij, &_vxpj4hn, &_s55j2o7, _lexdw6l, _1nlc8jp, _vmefozt, _4xkss0k, _zzf8s2z, _p422ta8, _6fokzix, _0dl5e11);
    
    float _flse45b = FLT_MAX;
    float _gl58kxo        = 1.0E-6f;
    {   unsigned int i;
        float _6zs62yw,   _2ivglye,   _nbzug0i[3] = {0.0f, 0.0f, 0.0f};
        
        for (i = 0u; i < _lexdw6l[_fsqg8x8]; i++)
        {   _6zs62yw = (-1.0f*_6fokzix[i*16 +2]*_6fokzix[i*16 +3]) *_gl58kxo;
            _2ivglye       = _6zs62yw /(sqrtf(_6fokzix[i*16 +9]*_6fokzix[i*16 +9] + _6fokzix[i*16 +10]*_6fokzix[i*16 +10])) *365.25f;
            _flse45b  = MIN(_2ivglye,_flse45b); 
            _nbzug0i[0]       = (_p422ta8[i*4 +0] == 1u)*(_nbzug0i[0] +1.0f)  +  (_p422ta8[i*4 +0] != 1u)*_nbzug0i[0];
            _nbzug0i[1]       = (_p422ta8[i*4 +0] == 2u)*(_nbzug0i[1] +1.0f)  +  (_p422ta8[i*4 +0] != 2u)*_nbzug0i[1];
            _nbzug0i[2]       = (_p422ta8[i*4 +0] == 3u)*(_nbzug0i[2] +1.0f)  +  (_p422ta8[i*4 +0] != 3u)*_nbzug0i[2];
        } 
        
        MPI_Allreduce(MPI_IN_PLACE,     _nbzug0i,       3, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, &_flse45b, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);
        _nbzug0i[0] /= (float)_vxpj4hn._3ru7myo;          _nbzug0i[1] /= (float)_vxpj4hn._3ru7myo;          _nbzug0i[2] /= (float)_vxpj4hn._3ru7myo;
        if (_fsqg8x8 == 0)         {       fprintf(stdout,"\nStability at catalog start:  unstable %2.2f%%     cond. stable %2.2f%%     stable %2.2f%%\n\n",_nbzug0i[0]*100.0f, _nbzug0i[1]*100.0f, _nbzug0i[2]*100.0f);        }
        _flse45b *= 5.0f;
    } 
    
    FILE *_tajw253;
    char _9tr74qg[512], _75aud9k[512];
    strcpy(_9tr74qg,    _dygflij);      snprintf(_75aud9k,sizeof(_75aud9k), "_%u.cat",_s55j2o7._jcmrsez);          strcat(_9tr74qg,    _75aud9k);
    
    unsigned int _xo100vh = 0u;
    double _ofn1r4k    = 0.0; 
    _qpki7p5(_fsqg8x8, _dygflij, &_atj9tze, &_vxpj4hn, &_s55j2o7, &_xo100vh, &_ofn1r4k, _lexdw6l, _1poxrte, _1nlc8jp, _fftouhy, _vmefozt, _4xkss0k, _zzf8s2z, _p422ta8, _6fokzix, _0dl5e11, _z39pecn, _7yt2m6t);
    _atj9tze._0ur5l8a    += _ofn1r4k;
    
    MPI_Barrier( MPI_COMM_WORLD );
    if (_fsqg8x8 == 0)
    {   if ((_tajw253 = fopen(_9tr74qg,"rb+"))     == NULL)     {   perror("Error -cant open/write RawCatalogFile...\n");      exit(EXIT_FAILURE);     }
    }
    
    free(_vmefozt);          free(_zzf8s2z);         free(_4xkss0k);      free(_e9ftaxp);
    free(_rgb93wk);          free(_3n3bzt9);         free(_rrpzlk3);      free(_9yo18pa);
    
    
    
    MPI_Status  _0get8t6;
    MPI_File    _npuevjl;
    
    unsigned int i,   j,   k,   _xj4pnjy,   _1o2rntj,   _aragdx5,   _o31yhxs,   _5v2wwg8,   _yce4qgm,   _nzrf5a7,   _egok01t,   _p2tixc7,   _t9nm3v2,   _mygty3a,   _hlv76mf;
    float _v78hu7f,   _qxmcthp,   _0isb5iq,   _mybd2qe,   _dxj37c9;

    double _z36u8ue   = _ofn1r4k;
    double _ihx74d0 = _ofn1r4k;
    
    _dxj37c9       = -FLT_MAX;
    
    double _qydbpgv = 0.0,   _wshmlfz = 0.0,   _jlpefja = 0.0,   _gjyv5g1 = 0.0,   _h6t1ntm;
    struct timespec   _glmy62x,   _1zocfi4,   _y25m8kb,   _8o1bj9r,   _rw6tm2o,   _o6l1f90,   _04k0lth,   _s9bp4m1;
    
    int           __attribute__((aligned(ALIGN64))) _9lzs34j[_0rmhrcn];     memset(_9lzs34j,  0, _0rmhrcn*sizeof(int));
    int           __attribute__((aligned(ALIGN64))) _sosx43w[_0rmhrcn];     memset(_sosx43w,  0, _0rmhrcn*sizeof(int));
    int           __attribute__((aligned(ALIGN64))) _bh6nqqb[_0rmhrcn];     memset(_bh6nqqb,  0, _0rmhrcn*sizeof(int));
    int           __attribute__((aligned(ALIGN64))) _1ni10kq[_0rmhrcn];     memset(_1ni10kq,  0, _0rmhrcn*sizeof(int));
    long long int __attribute__((aligned(ALIGN64))) _vsigbia[_0rmhrcn];    memset(_vsigbia, 0, _0rmhrcn*sizeof(long long int));
    long long int __attribute__((aligned(ALIGN64))) _f92movs[_0rmhrcn];     memset(_f92movs,  0, _0rmhrcn*sizeof(long long int));
    
    unsigned int  __attribute__((aligned(ALIGN64))) _xkc5i8x[16];
    float         __attribute__((aligned(ALIGN64))) _nbzug0i[16];
    float         __attribute__((aligned(ALIGN64))) _mrunjdu[3];
    float         __attribute__((aligned(ALIGN64))) _hj3muur[6];
    
    {   _moxsfh4 = _z0yncxl(ALIGN64, (_lexdw6l[_fsqg8x8]+1), sizeof(unsigned int));        _hx6zdcu = aligned_alloc(ALIGN64, _moxsfh4);      memset(_hx6zdcu, 0, _moxsfh4);

        _moxsfh4 = _z0yncxl(ALIGN64, 4*_lexdw6l[_fsqg8x8], sizeof(float));                 _d1vn8lf    = aligned_alloc(ALIGN64,_moxsfh4);       memset(_d1vn8lf,    0, _moxsfh4);
        _moxsfh4 = _z0yncxl(ALIGN64, 3*_lexdw6l[_fsqg8x8], sizeof(float));                 _ozctwjl   = aligned_alloc(ALIGN64, _moxsfh4);      memset(_ozctwjl,   0, _moxsfh4);
        _moxsfh4 = _z0yncxl(ALIGN64, 3*_vxpj4hn._3ru7myo, sizeof(float));                  _ny181ac   = aligned_alloc(ALIGN64, _moxsfh4);      memset(_ny181ac,   0, _moxsfh4);
        _moxsfh4 = _z0yncxl(ALIGN64, 2*_vxpj4hn._dn1th1x, sizeof(float));                _t997e4u    = aligned_alloc(ALIGN64, _moxsfh4);      memset(_t997e4u,    0, _moxsfh4);

        _moxsfh4 = _z0yncxl(ALIGN64, _lexdw6l[_fsqg8x8], sizeof(unsigned int));            _pj8t22c   = aligned_alloc(ALIGN64, _moxsfh4);      memset(_pj8t22c,   0, _moxsfh4);
                                                                                            _50b8y5d  = aligned_alloc(ALIGN64, _moxsfh4);      memset(_50b8y5d,  0, _moxsfh4);
        _moxsfh4 = _z0yncxl(ALIGN64, _lexdw6l[_fsqg8x8], sizeof(float));                   _dtg5nqj = aligned_alloc(ALIGN64, _moxsfh4);      memset(_dtg5nqj, 0, _moxsfh4);
                                                                                            _kc4lorx = aligned_alloc(ALIGN64, _moxsfh4);      memset(_kc4lorx, 0, _moxsfh4);
                                                                                            _pkxvhag = aligned_alloc(ALIGN64, _moxsfh4);      memset(_pkxvhag, 0, _moxsfh4);
                                                                                            _qijc1yd  = aligned_alloc(ALIGN64, _moxsfh4);      memset(_qijc1yd,  0, _moxsfh4);

        _moxsfh4 = _z0yncxl(ALIGN64, _vxpj4hn._3ru7myo, sizeof(unsigned int));             _zk94yxf= aligned_alloc(ALIGN64, _moxsfh4);      memset(_zk94yxf,0, _moxsfh4);
                                                                                            _iblejgu= aligned_alloc(ALIGN64, _moxsfh4);      memset(_iblejgu,0, _moxsfh4);
        _moxsfh4 = _z0yncxl(ALIGN64, _vxpj4hn._3ru7myo, sizeof(float));                    _2ljk4q5= aligned_alloc(ALIGN64, _moxsfh4);      memset(_2ljk4q5,0, _moxsfh4);
                                                                                            _v7c6s5n= aligned_alloc(ALIGN64, _moxsfh4);      memset(_v7c6s5n,0, _moxsfh4);
                                                                                            _kb4ehou= aligned_alloc(ALIGN64, _moxsfh4);      memset(_kb4ehou,0, _moxsfh4);
                                                                                            _6qlz2rw= aligned_alloc(ALIGN64, _moxsfh4);      memset(_6qlz2rw,0, _moxsfh4);

        _moxsfh4 = _z0yncxl(ALIGN64, _lexdw6l[_fsqg8x8], sizeof(float *));                 _66cnolj  = aligned_alloc(ALIGN64, _moxsfh4);
        _moxsfh4 = _z0yncxl(ALIGN64, 2*_atj9tze._at9mx2x, sizeof(float));
        for (i = 0u; i < _lexdw6l[_fsqg8x8]; i++)                                      {       _66cnolj[i] = aligned_alloc(ALIGN64, _moxsfh4);    memset(_66cnolj[i],0, _moxsfh4);       }
        
        if (_vxpj4hn._5izolnn > 0u)
        {   _moxsfh4 = _z0yncxl(ALIGN64, 3*_1poxrte[_fsqg8x8], sizeof(float));             _1vsmado    = aligned_alloc(ALIGN64, _moxsfh4);      memset(_1vsmado,  0, _moxsfh4);
            _moxsfh4 = _z0yncxl(ALIGN64, 4*_1poxrte[_fsqg8x8], sizeof(float));             _lqbmae3   = aligned_alloc(ALIGN64, _moxsfh4);      memset(_lqbmae3, 0, _moxsfh4);
            _moxsfh4 = _z0yncxl(ALIGN64, 4*_vxpj4hn._5izolnn, sizeof(float));              _gyspeng   = aligned_alloc(ALIGN64, _moxsfh4);      memset(_gyspeng, 0, _moxsfh4);
            _moxsfh4 = _z0yncxl(ALIGN64, 3*_vxpj4hn._w2maiet, sizeof(float));            _zkwdeek    = aligned_alloc(ALIGN64, _moxsfh4);      memset(_zkwdeek,  0, _moxsfh4);
        }
        
        _jbc2sws  = (const float *)_ny181ac;                 _wyx9trf  = (const float *)_gyspeng;
        _sh3tow4   = (const float *)_t997e4u;                  _8sryo4r   = (const float *)_zkwdeek;
    } 
    
    if (_fsqg8x8 == 0)     {   fprintf(stdout,"Starting the catalog....\n");           }
    
    for (i = 0u; i < _1poxrte[_fsqg8x8]; i++)      {   _7yt2m6t[i*16 +0] = 0.0f;        _7yt2m6t[i*16 +1] = 0.0f;         _7yt2m6t[i*16 +2] = 0.0f;          }
    
    
    
    clock_gettime(CLOCK_REALTIME, &_glmy62x);
    unsigned int _w2pcn9t  = 0u;
    unsigned int _gvyie6i = 0u;
    unsigned int _owgzrix       = 1u;

    while (_z36u8ue <= _atj9tze._0ur5l8a)
    {   
        clock_gettime(CLOCK_REALTIME, &_y25m8kb);
        
        for (i = 0u; i < _lexdw6l[_fsqg8x8]; i++)              {       memcpy(&_d1vn8lf[i*4 +0], &_6fokzix[i*16 +0], 4*sizeof(float));           }
        for (i = 0u; i < _1poxrte[_fsqg8x8]; i++)              {       memcpy(&_1vsmado[i*3 +0], &_7yt2m6t[i*16 +0], 3*sizeof(float));           }
        
        _v78hu7f  = _flse45b; 
        
        k   = 1;
        _1o2rntj = 0;
        while (k == 1)
        {   _1o2rntj++;
            _nbzug0i[12] = MAX(0.0f, expf(-1.0f*_v78hu7f/_vxpj4hn._lw2lm34)); 
            _nbzug0i[13] = MAX(0.0f, expf(-1.0f*_v78hu7f/_vxpj4hn._e1g983e));
            _nbzug0i[14] = (_atj9tze._kkrzd5e/(1.0f*_atj9tze._douadvw)) *(_vxpj4hn._7odkef8/_atj9tze._mzwffhh); 
            _nbzug0i[15] = FLT_MAX; 
            
            for (i = 0u; i < _lexdw6l[_fsqg8x8]; i++)
            {   _h6t1ntm         = (_z36u8ue + (double)_v78hu7f)*(double)_z39pecn[i*8 +5]; 
                _6fokzix[i*16 +0]  = (_z39pecn[i*8 +3] <= _h6t1ntm)*(_d1vn8lf[i*4 +0] +_v78hu7f*_6fokzix[i*16 +9])   +  (_z39pecn[i*8 +3] > _h6t1ntm)*_d1vn8lf[i*4 +0];
                _6fokzix[i*16 +1]  = (_z39pecn[i*8 +3] <= _h6t1ntm)*(_d1vn8lf[i*4 +1] +_v78hu7f*_6fokzix[i*16 +10])  +  (_z39pecn[i*8 +3] > _h6t1ntm)*_d1vn8lf[i*4 +1];
                _6fokzix[i*16 +2]  = _d1vn8lf[i*4 +2];
                _6fokzix[i*16 +3]  = _d1vn8lf[i*4 +3] - (1.0f - _nbzug0i[12])*_0dl5e11[i*16 +7];
                _6fokzix[i*16 +15] = 0.0f;
                
                
                
            }
            for (i = 0u; i < _1poxrte[_fsqg8x8]; i++)
            {   memcpy(&_7yt2m6t[i*16 +0], &_1vsmado[i*3 +0], 3*sizeof(float));
            }
            
            for (_xj4pnjy = 0; _xj4pnjy < _owgzrix; _xj4pnjy++)
            {   
                _5v2wwg8     = 0u;
                _nzrf5a7 = 0u;
                for (i = 0u; i < _lexdw6l[_fsqg8x8]; i++)
                {   
                    _nbzug0i[0]        = sqrtf( _6fokzix[i*16 +0]*_6fokzix[i*16 +0] + _6fokzix[i*16 +1]*_6fokzix[i*16 +1] );
                    _nbzug0i[0]        = (_p422ta8[i*4 +0] == 1u)*0.0f  +  (_p422ta8[i*4 +0] != 1u)*_nbzug0i[0]; 
                    _nbzug0i[1]        = _6fokzix[i*16 +3] *-1.0f*_6fokzix[i*16 +2] ; 
                    
                    _nbzug0i[3]        = MAX(0.0f, (_nbzug0i[0] -_nbzug0i[1])); 
                    
                    _nbzug0i[0]        = (_nbzug0i[0] <= 0.0f)*-1.0f + (_nbzug0i[0] > 0.0f)*_nbzug0i[0]; 
                    _nbzug0i[1]        = -1.0f*(((_nbzug0i[3]/_nbzug0i[0])*_6fokzix[i*16 +0])/_6fokzix[i*16 +5]);
                    _nbzug0i[2]        = -1.0f*(((_nbzug0i[3]/_nbzug0i[0])*_6fokzix[i*16 +1])/_6fokzix[i*16 +6]);
                    
                    _6fokzix[i*16 +15]       = (_nbzug0i[3] <= _atj9tze._14ve92m)*_6fokzix[i*16 +15]        +  (_nbzug0i[3] > _atj9tze._14ve92m)*(_6fokzix[i*16 +15] + sqrtf(_nbzug0i[1]*_nbzug0i[1] + _nbzug0i[2]*_nbzug0i[2])); 
                    _ozctwjl[_5v2wwg8*3 +0] = (_nbzug0i[3] <= _atj9tze._14ve92m)*_ozctwjl[_5v2wwg8*3 +0]  +  (_nbzug0i[3] > _atj9tze._14ve92m)*(float)(_1nlc8jp[i]);
                    _ozctwjl[_5v2wwg8*3 +1] = (_nbzug0i[3] <= _atj9tze._14ve92m)*_ozctwjl[_5v2wwg8*3 +1]  +  (_nbzug0i[3] > _atj9tze._14ve92m)*_nbzug0i[1]*_6fokzix[i*16 +8];
                    _ozctwjl[_5v2wwg8*3 +2] = (_nbzug0i[3] <= _atj9tze._14ve92m)*_ozctwjl[_5v2wwg8*3 +2]  +  (_nbzug0i[3] > _atj9tze._14ve92m)*_nbzug0i[2]*_6fokzix[i*16 +8];
                    _5v2wwg8               = (_nbzug0i[3] <= _atj9tze._14ve92m)*_5v2wwg8                +  (_nbzug0i[3] > _atj9tze._14ve92m)*(_5v2wwg8 + 1u);
                }
                
                _yce4qgm     = 0u;
                _egok01t = 0u;
                for (i = 0u; i < _1poxrte[_fsqg8x8]; i++)
                {   _nbzug0i[0] = sqrtf( _7yt2m6t[i*16 +0]*_7yt2m6t[i*16 +0] + _7yt2m6t[i*16 +1]*_7yt2m6t[i*16 +1] ); 
                    _nbzug0i[1] = fabsf(_7yt2m6t[i*16 +2]); 
                    _nbzug0i[2] = MAX(0.0f, (1.0f - _nbzug0i[13])*_nbzug0i[0]); 
                    _nbzug0i[3] = MAX(0.0f, (1.0f - _nbzug0i[13])*_nbzug0i[1]); 
                    _nbzug0i[4] = sqrtf(_nbzug0i[2]*_nbzug0i[2] + _nbzug0i[3]*_nbzug0i[3]);
                    
                    _nbzug0i[0] = (_nbzug0i[0] <= 0.0f)*-1.0f  +  (_nbzug0i[0] > 0.0f)*_nbzug0i[0]; 
                    _nbzug0i[1] = (_nbzug0i[1] <= 0.0f)*-1.0f  +  (_nbzug0i[1] > 0.0f)*_nbzug0i[1];

                    _nbzug0i[5] = (_nbzug0i[0] <= 0.0f)*0.0f   +  (_nbzug0i[0] > 0.0f)*(-1.0f*((_nbzug0i[2]/_nbzug0i[0])*_7yt2m6t[i*16 +0])/_7yt2m6t[i*16 +5]);
                    _nbzug0i[6] = (_nbzug0i[0] <= 0.0f)*0.0f   +  (_nbzug0i[0] > 0.0f)*(-1.0f*((_nbzug0i[2]/_nbzug0i[0])*_7yt2m6t[i*16 +1])/_7yt2m6t[i*16 +6]);
                    _nbzug0i[7] = (_nbzug0i[1] <= 0.0f)*0.0f   +  (_nbzug0i[1] > 0.0f)*(-1.0f*((_nbzug0i[3]/_nbzug0i[1])*_7yt2m6t[i*16 +2])/_7yt2m6t[i*16 +7]);
                    
                    _lqbmae3[_yce4qgm*4 +0] = (_nbzug0i[4] <= _atj9tze._14ve92m)*_lqbmae3[_yce4qgm*4 +0]  +  (_nbzug0i[4] > _atj9tze._14ve92m)*(float)(_fftouhy[i]);
                    _lqbmae3[_yce4qgm*4 +1] = (_nbzug0i[4] <= _atj9tze._14ve92m)*_lqbmae3[_yce4qgm*4 +1]  +  (_nbzug0i[4] > _atj9tze._14ve92m)*_nbzug0i[5]*_7yt2m6t[i*16 +8];
                    _lqbmae3[_yce4qgm*4 +2] = (_nbzug0i[4] <= _atj9tze._14ve92m)*_lqbmae3[_yce4qgm*4 +2]  +  (_nbzug0i[4] > _atj9tze._14ve92m)*_nbzug0i[6]*_7yt2m6t[i*16 +8];
                    _lqbmae3[_yce4qgm*4 +3] = (_nbzug0i[4] <= _atj9tze._14ve92m)*_lqbmae3[_yce4qgm*4 +3]  +  (_nbzug0i[4] > _atj9tze._14ve92m)*_nbzug0i[7]*_7yt2m6t[i*16 +8];
                    _yce4qgm               = (_nbzug0i[4] <= _atj9tze._14ve92m)*_yce4qgm                +  (_nbzug0i[4] > _atj9tze._14ve92m)*(_yce4qgm + 1u);
                }
                
                MPI_Allgather(&_5v2wwg8, 1, MPI_UNSIGNED, _sosx43w, 1, MPI_INT, MPI_COMM_WORLD);
                _9lzs34j[0] = 0;
                _nzrf5a7    = _sosx43w[0];
                for (i = 1; i < _0rmhrcn; i++)         {       _nzrf5a7 += _sosx43w[i];        _sosx43w[i-1] *= 3;           _9lzs34j[i] = _9lzs34j[i-1] + _sosx43w[i-1];         }
                _sosx43w[_0rmhrcn-1] *= 3;
                MPI_Allgatherv(_ozctwjl, _sosx43w[_fsqg8x8], MPI_FLOAT, _ny181ac, _sosx43w, _9lzs34j, MPI_FLOAT, MPI_COMM_WORLD);
                
                if (_vxpj4hn._5izolnn > 0u)
                {   MPI_Allgather(&_yce4qgm, 1, MPI_UNSIGNED, _1ni10kq, 1, MPI_INT, MPI_COMM_WORLD);
                    _bh6nqqb[0] = 0u;
                    _egok01t    = _1ni10kq[0];
                    for (i = 1; i < _0rmhrcn; i++)     {       _egok01t += _1ni10kq[i];        _1ni10kq[i-1] *= 4;           _bh6nqqb[i] = _bh6nqqb[i-1] + _1ni10kq[i-1];         }
                    _1ni10kq[_0rmhrcn-1] *= 4;
                    MPI_Allgatherv(_lqbmae3, _1ni10kq[_fsqg8x8], MPI_FLOAT, _gyspeng, _1ni10kq, _bh6nqqb, MPI_FLOAT, MPI_COMM_WORLD);
                }
                
                for (i = 0u; i < _lexdw6l[_fsqg8x8]; i++)
                {   
                    if (_vxpj4hn._5izolnn > 0u)
                    {   _xkc5i8x[0] = (_egok01t <= _qpwbsf3[i])*1u  +  (_egok01t > _qpwbsf3[i])*2u;
                        _xkc5i8x[0] = (_egok01t <= 0u)*0u            +  (_egok01t > 0u)*_xkc5i8x[0];
                        
                        if      (_xkc5i8x[0] == 1u)
                        {   _nbzug0i[0] = 0.0f;   _nbzug0i[1] = 0.0f;
                            for (j = 0u; j < _egok01t; j++) 
                            {   _xkc5i8x[2]   = (unsigned int)_wyx9trf[j*4 +0]; 
                                _xkc5i8x[1]   = _oqgvzhq[i][_xkc5i8x[2]] *3; 
                                _xkc5i8x[2]   = j*4;
                                _nbzug0i[0]  += (_wyx9trf[_xkc5i8x[2] +1] *_yj9boq0[i][_xkc5i8x[1] +0] + _wyx9trf[_xkc5i8x[2] +2] *_yj9boq0[i][_xkc5i8x[1] +1] + _wyx9trf[_xkc5i8x[2]+3] *_yj9boq0[i][_xkc5i8x[1] +2]);
                                _nbzug0i[1]  += (_wyx9trf[_xkc5i8x[2] +1] *_ljefy3z[i][_xkc5i8x[1] +0] + _wyx9trf[_xkc5i8x[2] +2] *_ljefy3z[i][_xkc5i8x[1] +1] + _wyx9trf[_xkc5i8x[2]+3] *_ljefy3z[i][_xkc5i8x[1] +2]);
                            }
                            _6fokzix[i*16 +0] += _nbzug0i[0];
                            _6fokzix[i*16 +1] += _nbzug0i[1];
                        }
                        else if (_xkc5i8x[0] == 2u)
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
                    }
                    
                    _xkc5i8x[0] = (_nzrf5a7 <= _qnsn2s8[i])*1u  +  (_nzrf5a7 > _qnsn2s8[i])*2u;
                    _xkc5i8x[0] = (_nzrf5a7 <= 0u)*0u            +  (_nzrf5a7 > 0u)*_xkc5i8x[0];
                    
                    if      (_xkc5i8x[0] == 1u)
                    {   _nbzug0i[0] = 0.0f;   _nbzug0i[1] = 0.0f;
                        for (j = 0u; j < _nzrf5a7; j++) 
                        {   _xkc5i8x[2]   = (unsigned int)_jbc2sws[j*3 +0]; 
                            _xkc5i8x[1]   = _5gzpl1q[i][_xkc5i8x[2]] *2; 
                            _xkc5i8x[2]   = j*3;
                            _nbzug0i[0]  += (_jbc2sws[_xkc5i8x[2] +1]*_ckeh5sw[i][_xkc5i8x[1] +0] + _jbc2sws[_xkc5i8x[2] +2]*_ckeh5sw[i][_xkc5i8x[1] +1] );
                            _nbzug0i[1]  += (_jbc2sws[_xkc5i8x[2] +1]*_3ac6pe6[i][_xkc5i8x[1] +0] + _jbc2sws[_xkc5i8x[2] +2]*_3ac6pe6[i][_xkc5i8x[1] +1] );
                        }
                        _6fokzix[i*16 +0] += _nbzug0i[0];
                        _6fokzix[i*16 +1] += _nbzug0i[1];
                    }
                    else if (_xkc5i8x[0] == 2u)
                    {   memset(_t997e4u, 0, 2*_qnsn2s8[i]*sizeof(float) );
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
                } 
                
                for (i = 0u; i < _1poxrte[_fsqg8x8]; i++)
                {   _xkc5i8x[0] = (_egok01t <= _iymyrsv[i])*1u  +  (_egok01t > _iymyrsv[i])*2u;
                    _xkc5i8x[0] =           (_egok01t <= 0u)*0u  +  (_egok01t > 0u)*_xkc5i8x[0];
                    
                    if      ( _xkc5i8x[0] == 1)
                    {   _nbzug0i[0] = 0.0f;   _nbzug0i[1] = 0.0f;   _nbzug0i[2] = 0.0f;
                        for (j = 0u; j < _egok01t; j++) 
                        {   _xkc5i8x[2]  = (unsigned int)_wyx9trf[j*4 +0]; 
                            _xkc5i8x[1]  = _1ajo091[i][_xkc5i8x[2]] *3; 
                            _xkc5i8x[2]  = j*4;
                            _nbzug0i[0] += (_wyx9trf[_xkc5i8x[2] +1]*_h61iaoi[i][_xkc5i8x[1] +0] + _wyx9trf[_xkc5i8x[2] +2]*_h61iaoi[i][_xkc5i8x[1] +1] + _wyx9trf[_xkc5i8x[2] +3]*_h61iaoi[i][_xkc5i8x[1] +2]);
                            _nbzug0i[1] += (_wyx9trf[_xkc5i8x[2] +1]*_plozjcn[i][_xkc5i8x[1] +0] + _wyx9trf[_xkc5i8x[2] +2]*_plozjcn[i][_xkc5i8x[1] +1] + _wyx9trf[_xkc5i8x[2] +3]*_plozjcn[i][_xkc5i8x[1] +2]);
                            _nbzug0i[2] += (_wyx9trf[_xkc5i8x[2] +1]*_qk4ltah[i][_xkc5i8x[1] +0] + _wyx9trf[_xkc5i8x[2] +2]*_qk4ltah[i][_xkc5i8x[1] +1] + _wyx9trf[_xkc5i8x[2] +3]*_qk4ltah[i][_xkc5i8x[1] +2]);
                        }
                        _7yt2m6t[i*16 +0] += _nbzug0i[0];
                        _7yt2m6t[i*16 +1] += _nbzug0i[1];
                        _7yt2m6t[i*16 +2] += _nbzug0i[2];
                    }
                    else if ( _xkc5i8x[0] == 2)
                    {   memset(_zkwdeek, 0, 3*_iymyrsv[i]*sizeof(float) );
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
                    }
                    
                    _xkc5i8x[0] = (_nzrf5a7 <= _iob0nlw[i])*1u  +  (_nzrf5a7 > _iob0nlw[i])*2u;
                    _xkc5i8x[0] =           (_nzrf5a7 <= 0u)*0u  +  (_nzrf5a7 > 0u)*_xkc5i8x[0];
                    
                    if      ( _xkc5i8x[0] == 1)
                    {   _nbzug0i[0] = 0.0f;   _nbzug0i[1] = 0.0f;   _nbzug0i[2] = 0.0f;
                        for (j = 0u; j < _nzrf5a7; j++) 
                        {   _xkc5i8x[2]  = (unsigned int)_jbc2sws[j*3 +0]; 
                            _xkc5i8x[1]  = _3b9so0u[i][_xkc5i8x[2]] *2; 
                            _xkc5i8x[2]  = j*3;
                            _nbzug0i[0] += (_jbc2sws[_xkc5i8x[2] +1]*_0ljch2v[i][_xkc5i8x[1] +0] + _jbc2sws[_xkc5i8x[2] +2]*_0ljch2v[i][_xkc5i8x[1] +1]);
                            _nbzug0i[1] += (_jbc2sws[_xkc5i8x[2] +1]*_0dtirc9[i][_xkc5i8x[1] +0] + _jbc2sws[_xkc5i8x[2] +2]*_0dtirc9[i][_xkc5i8x[1] +1]);
                            _nbzug0i[2] += (_jbc2sws[_xkc5i8x[2] +1]*_9dde3cu[i][_xkc5i8x[1] +0] + _jbc2sws[_xkc5i8x[2] +2]*_9dde3cu[i][_xkc5i8x[1] +1]);
                        }
                        _7yt2m6t[i*16 +0] += _nbzug0i[0];
                        _7yt2m6t[i*16 +1] += _nbzug0i[1];
                        _7yt2m6t[i*16 +2] += _nbzug0i[2];
                    }
                    else if ( _xkc5i8x[0] == 2)
                    {   memset(_t997e4u, 0, 2*_iob0nlw[i]*sizeof(float) );
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
                
            }
            
            for (i = 0u; i < _lexdw6l[_fsqg8x8]; i++)
            {   
                
                
                
                if (_p422ta8[i*4 +0]  == 1u)
                {   _nbzug0i[0] = sqrtf(_d1vn8lf[i*4 +0]*_d1vn8lf[i*4 +0]  + _d1vn8lf[i*4 +1]*_d1vn8lf[i*4 +1]); 
                    _nbzug0i[1] = _d1vn8lf[i*4 +3]* -1.0f*_d1vn8lf[i*4 +2] + 1.05f*_nbzug0i[14]; 
                    _nbzug0i[2] = sqrtf(_6fokzix[i*16 + 0]*_6fokzix[i*16 + 0]  + _6fokzix[i*16 + 1]*_6fokzix[i*16 + 1]); 
                    _nbzug0i[3] = _6fokzix[i*16 +3]*  -1.0f*_6fokzix[i*16 +2]  + 1.05f*_nbzug0i[14]; 
                    
                    _nbzug0i[4] = _nbzug0i[0] - _nbzug0i[1]; 
                    _nbzug0i[5] = _nbzug0i[2] - _nbzug0i[3]; 
                    
                    _nbzug0i[8] = sqrtf(_6fokzix[i*16 +9]*_6fokzix[i*16 +9] + _6fokzix[i*16 +10]*_6fokzix[i*16 +10]); 
                    _nbzug0i[8] = (_nbzug0i[8] > 0.0f)*_nbzug0i[8] + (_nbzug0i[8] > 0.0f)*-1.0f;
                    _nbzug0i[6] = (MAX(_nbzug0i[0], _nbzug0i[1])*_gl58kxo)/_nbzug0i[8]; 
                    _nbzug0i[6] = (_nbzug0i[8] > 0.0f)*_nbzug0i[6]  +  (_nbzug0i[8] <= 0.0f)*-FLT_MAX; 
                    
                    _nbzug0i[8] = (_nbzug0i[5] > _nbzug0i[4])*(_nbzug0i[5] - _nbzug0i[4])                      +  (_nbzug0i[5] <= _nbzug0i[4])*-1.0f; 
                    _nbzug0i[7] = (_nbzug0i[5] > _nbzug0i[4])*((-1.0f*_nbzug0i[4])/_nbzug0i[8] *_v78hu7f)  +  (_nbzug0i[5] <= _nbzug0i[4])*FLT_MAX; 

                    _nbzug0i[7] = MAX(_nbzug0i[7], 0.0f); 

                    _xkc5i8x[0] = (_nbzug0i[7] > _v78hu7f)*1u                     +  (_nbzug0i[7] <= _v78hu7f)*0u;
                    _xkc5i8x[1] = (fabsf(_nbzug0i[7] - _v78hu7f) <= _nbzug0i[6])*1u  +  (fabsf(_nbzug0i[7] - _v78hu7f) > _nbzug0i[6])*0u;
                    _nbzug0i[7] = (_xkc5i8x[0] == 1u)*MAX(_v78hu7f+_nbzug0i[6], (_nbzug0i[7]))  +  (_xkc5i8x[0] != 1u)*((_xkc5i8x[1] == 1u)*_v78hu7f  +  (_xkc5i8x[1] != 1u)*_nbzug0i[7]);

                    _nbzug0i[7] = (_nbzug0i[0] >= _nbzug0i[1])*0.0f  +  (_nbzug0i[0] < _nbzug0i[1])*_nbzug0i[7];
                    _nbzug0i[15]= MIN(_nbzug0i[7], _nbzug0i[15]);
                 }
            }
            
            MPI_Allreduce(MPI_IN_PLACE, &_nbzug0i[15], 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);
            
            if (_nbzug0i[15] == 0.0f)
            {   k = 0;
                _v78hu7f = 0.0f;
                for (i = 0u; i < _lexdw6l[_fsqg8x8]; i++)      {       memcpy(&_6fokzix[i*16 +0], &_d1vn8lf[i*4 +0], 4*sizeof(float));           _6fokzix[i*16 +15]= 0.0f;       }
                for (i = 0u; i < _1poxrte[_fsqg8x8]; i++)      {       memcpy(&_7yt2m6t[i*16 +0], &_1vsmado[i*3 +0], 3*sizeof(float));                                       }
            }
            else
            {
                if      (_nbzug0i[15]  > _v78hu7f) 
                {   _z36u8ue   += (double)_v78hu7f;
                    _v78hu7f = (_nbzug0i[15] -_v78hu7f); 
                    for (i = 0u; i < _lexdw6l[_fsqg8x8]; i++)  {       memcpy(&_d1vn8lf[i*4 +0], &_6fokzix[i*16 +0], 4*sizeof(float));           _0dl5e11[i*16 +7] = _6fokzix[i*16 +3] - _0dl5e11[i*16+2];           _z39pecn[i*8 +3] += _6fokzix[i*16 +15];                }
                    for (i = 0u; i < _1poxrte[_fsqg8x8]; i++)  {       memcpy(&_1vsmado[i*3 +0], &_7yt2m6t[i*16 +0], 3*sizeof(float));                                                                                                                    }
                    _w2pcn9t += 1u;
                }
                else if (_nbzug0i[15] == _v78hu7f)         {       k = 0;                                          }
                else if (_nbzug0i[15]  < _v78hu7f)         {       _v78hu7f = _nbzug0i[15];                       }
            }
            
        }
        
        _z36u8ue  += (double)_v78hu7f;
        for (i = 0u; i < _lexdw6l[_fsqg8x8]; i++)              {       _z39pecn[i*8 +3] += _6fokzix[i*16 +15];               }
        
        _o31yhxs  = 0u;
        _aragdx5  = 0u;
        _0isb5iq  = 0.0f;
        memset(_mrunjdu, 0, 3*sizeof(float));
        for (i = 0u; i < _lexdw6l[_fsqg8x8]; i++)
        {   _nbzug0i[0] = sqrtf( _6fokzix[i*16 +0]*_6fokzix[i*16 +0] + _6fokzix[i*16 +1]*_6fokzix[i*16 +1] );
            _nbzug0i[0] = (_p422ta8[i*4 +0] == 1u)*_nbzug0i[0]  +  (_p422ta8[i*4 +0] != 1u)*0.0f; 
            _nbzug0i[1] = _6fokzix[i*16 +3] *-1.0f*_6fokzix[i*16 +2]; 
            _nbzug0i[2] =  ((1.0f*_atj9tze._douadvw)/_atj9tze._kkrzd5e) *MAX(0.0f,(_nbzug0i[0] -_nbzug0i[1])) *_atj9tze._mzwffhh; 
            
            _xkc5i8x[0] = (_nbzug0i[2] >= _vxpj4hn._7odkef8)*1u       +  (_nbzug0i[2]  < _vxpj4hn._7odkef8)*0u;
            _nbzug0i[2] =                (_xkc5i8x[0] == 1u)*_nbzug0i[2]  +  (_xkc5i8x[0] != 1u)*0.0f; 
            
            _mrunjdu[0] = (_nbzug0i[2] > _0isb5iq)*_z39pecn[i*8 +0]  +  (_nbzug0i[2] <= _0isb5iq)*_mrunjdu[0];
            _mrunjdu[1] = (_nbzug0i[2] > _0isb5iq)*_z39pecn[i*8 +1]  +  (_nbzug0i[2] <= _0isb5iq)*_mrunjdu[1];
            _mrunjdu[2] = (_nbzug0i[2] > _0isb5iq)*_z39pecn[i*8 +2]  +  (_nbzug0i[2] <= _0isb5iq)*_mrunjdu[2];
            _aragdx5   = (_nbzug0i[2] > _0isb5iq)*i               +  (_nbzug0i[2] <= _0isb5iq)*_aragdx5;
            _0isb5iq   = (_nbzug0i[2] > _0isb5iq)*_nbzug0i[2]        +  (_nbzug0i[2] <= _0isb5iq)*_0isb5iq;
        }
        
        {   struct { float _4fieecy;  int _k2wbxan;  } _g8vuazx[1], _yllwb7s[1];
            _g8vuazx[0]._4fieecy  = _0isb5iq;         _g8vuazx[0]._k2wbxan = _fsqg8x8;
            
            MPI_Allreduce(_g8vuazx, _yllwb7s, 1, MPI_FLOAT_INT, MPI_MAXLOC, MPI_COMM_WORLD);
            
            MPI_Bcast(&_yllwb7s[0]._k2wbxan, 1, MPI_INT,   _yllwb7s[0]._k2wbxan, MPI_COMM_WORLD);
            MPI_Bcast(_mrunjdu,     3, MPI_FLOAT, _yllwb7s[0]._k2wbxan, MPI_COMM_WORLD);
            _o31yhxs = _yllwb7s[0]._k2wbxan;
        } 
        
        _mygty3a = (_0isb5iq > 0.0f)*1u + (_0isb5iq <= 0.0f)*0u;
        MPI_Allreduce(MPI_IN_PLACE, &_mygty3a, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
        if (_mygty3a == 0u)   {   perror("something went wrong, left Load2NextEQ without finding next event");       exit(EXIT_FAILURE);      }
        
        clock_gettime(CLOCK_REALTIME, &_8o1bj9r);
        _wshmlfz += (_8o1bj9r.tv_sec - _y25m8kb.tv_sec) + (_8o1bj9r.tv_nsec - _y25m8kb.tv_nsec)/1.0E+9;
        
        
        for (i = 0u; i < _lexdw6l[_fsqg8x8]; i++)
        {   memcpy(&_d1vn8lf[i*4 +0], &_6fokzix[i*16 +0], 2*sizeof(float));
            
            _nbzug0i[0]      = (_0dl5e11[i*16 +3] *-1.0f*_6fokzix[i*16 +2]) + _0dl5e11[i*16 +1]*_6fokzix[i*16 +4]; 
            _nbzug0i[0]     /= (-1.0f *_6fokzix[i*16 +2]);
            _6fokzix[i*16 +3] = MAX(_6fokzix[i*16 +3], _nbzug0i[0]);
            _0dl5e11[i*16 +5] = _6fokzix[i*16 +3]; 
            _0dl5e11[i*16 +6] = _6fokzix[i*16 +3]; 
            
            _6fokzix[i*16 +2] = (_fsqg8x8 != _o31yhxs)*((i != _aragdx5)*(_6fokzix[i*16 +2] + ADDEDSTRGTHVALUE)  +  (i == _aragdx5)*_6fokzix[i*16 +2])  +  (_fsqg8x8 == _o31yhxs)*_6fokzix[i*16 +2];
            
        }
        
        clock_gettime(CLOCK_REALTIME, &_rw6tm2o);
        clock_gettime(CLOCK_REALTIME, &_04k0lth);
        
        _nzrf5a7 = 1u;
        _p2tixc7   = 0u;
        _hlv76mf = 0u;
        _mybd2qe  = 0.0f; 

        while (_nzrf5a7 > 0u)
        {   
            _5v2wwg8     = 0u;
            _nzrf5a7 = 0u;
            
            for (i = 0u; i < _lexdw6l[_fsqg8x8]; i++)
            {   _6fokzix[i*16 +12] = (_6fokzix[i*16 +11] <= 0.0f)*0.0f                               +  (_6fokzix[i*16 +11] > 0.0f)*_6fokzix[i*16 +12]; 
                _0dl5e11[i*16 +5]  = (_6fokzix[i*16 +11] <= 0.0f)*MAX(_0dl5e11[i*16 +3], _0dl5e11[i*16 +5])  +  (_6fokzix[i*16 +11] > 0.0f)*_0dl5e11[i*16 +5];  
                _6fokzix[i*16 +3]  = (_6fokzix[i*16 +11] <= 0.0f)*MAX(_0dl5e11[i*16 +3], _6fokzix[i*16 +3])  +  (_6fokzix[i*16 +11] > 0.0f)*_6fokzix[i*16 +3];  
                
                _nbzug0i[1]       = (_0dl5e11[i*16 +3] - _0dl5e11[i*16 +5])/(1.0f*_6fokzix[i*16 +4]); 
                _nbzug0i[2]       = (_0dl5e11[i*16 +4] - _0dl5e11[i*16 +3])/(0.5f*_6fokzix[i*16 +4]); 
                _nbzug0i[3]       = _0dl5e11[i*16 +3] - _nbzug0i[2]*_6fokzix[i*16 +4]; 
                _nbzug0i[4]       = (_6fokzix[i*16 +12] <= _6fokzix[i*16 +4]) *(_0dl5e11[i*16 +5] + _nbzug0i[1]*_6fokzix[i*16 +12])  +  (_6fokzix[i*16 +12] > _6fokzix[i*16 +4])*(_nbzug0i[3] + _nbzug0i[2]*_6fokzix[i*16 +12]);
                _nbzug0i[4]       = MAX(_nbzug0i[4], _0dl5e11[i*16 +4]); 
                
                _nbzug0i[5]       = _6fokzix[i*16 +3] + _vxpj4hn._3jczucq*(_0dl5e11[i*16 +6] - _6fokzix[i*16 +3]); 
                
                _6fokzix[i*16 +3]  = (_p422ta8[i*4 +1]   ==   0u)*_6fokzix[i*16 +3] +  (_p422ta8[i*4 +1] != 0u)*(   (_6fokzix[i*16 +11] <= 0.0f)*_nbzug0i[5]  + (_6fokzix[i*16 +11] > 0.0f)*_nbzug0i[4]   ); 
                _0dl5e11[i*16 +5]  = (_6fokzix[i*16 +11] <= 0.0f)*_nbzug0i[5]      +  (_6fokzix[i*16 +11] > 0.0f)*_0dl5e11[i*16 +5];
                _6fokzix[i*16 +11] = 0.0f;
                
                _nbzug0i[0]       = sqrtf( _6fokzix[i*16 +0]*_6fokzix[i*16 +0] + _6fokzix[i*16 +1]*_6fokzix[i*16 +1] ); 
                
                _6fokzix[i*16 +2]  = MIN(0.0f, _6fokzix[i*16 +2]);
                _nbzug0i[15]      = _atj9tze._rtj5pht + _atj9tze._dswm31y*-1.0f*_6fokzix[i*16 +2]; 
                _nbzug0i[1]       = _6fokzix[i*16 +3] *-1.0f*_6fokzix[i*16 +2]; 
                _nbzug0i[3]       = _nbzug0i[0] - _nbzug0i[1]; 
                _nbzug0i[5]       = ((1.0f*_atj9tze._douadvw)/_atj9tze._kkrzd5e) *_nbzug0i[3] *_atj9tze._mzwffhh; 
                
                _xkc5i8x[0]       = (_nbzug0i[5] < _vxpj4hn._7odkef8)*0u  +  (_nbzug0i[5] >= _vxpj4hn._7odkef8)*1u;
                _xkc5i8x[1]       = (_p422ta8[i*4 +1] == 2u)*2u             +  (_p422ta8[i*4 +1] != 2u)*1u; 
                _xkc5i8x[0]      *= _xkc5i8x[1]; 
                
                if      (_xkc5i8x[0] == 1u) 
                {   
                    _hx6zdcu[_p2tixc7] = (_p422ta8[i*4 +1] == 0u)*i          +  (_p422ta8[i*4 +1] != 0u)*_hx6zdcu[_p2tixc7];  
                    _p2tixc7        = (_p422ta8[i*4 +1] == 0u)*(_p2tixc7 +1u) +  (_p422ta8[i*4 +1] != 0u)*_p2tixc7;
                    
                    _p422ta8[i*4 +2]    = (_p422ta8[i*4 +1] == 0u)*_hlv76mf     +  (_p422ta8[i*4 +1] != 0u)*_p422ta8[i*4 +2];
                    _p422ta8[i*4 +1]    = (_p422ta8[i*4 +1] == 0u)*1u             +  (_p422ta8[i*4 +1] != 0u)*_p422ta8[i*4 +1];
                    
                    _p422ta8[i*4 +1]    = (_nbzug0i[0] > _nbzug0i[15])*2u           +  (_nbzug0i[0] <= _nbzug0i[15])*_p422ta8[i*4 +1]; 
                    
                    _p422ta8[i*4 +3]    = _hlv76mf +1u - _p422ta8[i*4 +2];
                    _6fokzix[i*16 +15]  = fabsf(_nbzug0i[5]);
                    
                    _nbzug0i[5]        = -1.0f*((_nbzug0i[3]/_nbzug0i[0]) *_6fokzix[i*16 +0]) / _6fokzix[i*16 +5];
                    _nbzug0i[6]        = -1.0f*((_nbzug0i[3]/_nbzug0i[0]) *_6fokzix[i*16 +1]) / _6fokzix[i*16 +6];
                    _nbzug0i[7]        = _6fokzix[i*16 +15]/sqrtf(_nbzug0i[5]*_nbzug0i[5] +_nbzug0i[6]*_nbzug0i[6]);
                    _nbzug0i[7]        = MIN(_nbzug0i[7], 1.0f);
                    _nbzug0i[5]       *= _nbzug0i[7];
                    _nbzug0i[6]       *= _nbzug0i[7];
                    
                    _6fokzix[i*16 +11]  = sqrtf(_nbzug0i[5]*_nbzug0i[5] +_nbzug0i[6]*_nbzug0i[6]);
                    _6fokzix[i*16 +12] += _6fokzix[i*16 +11];
                    _6fokzix[i*16 +13] += _nbzug0i[5];
                    _6fokzix[i*16 +14] += _nbzug0i[6];
                    _mybd2qe      += (_6fokzix[i*16 +11]*_6fokzix[i*16 +8]);
                    
                    _66cnolj[i][ (2u*(_hlv76mf - _p422ta8[i*4 +2])) +0] = _nbzug0i[5];
                    _66cnolj[i][ (2u*(_hlv76mf - _p422ta8[i*4 +2])) +1] = _nbzug0i[6];
                    
                    _ozctwjl[_5v2wwg8*3 +0] = (float)(_1nlc8jp[i]);
                    _ozctwjl[_5v2wwg8*3 +1] = _nbzug0i[5]*_6fokzix[i*16 +8];
                    _ozctwjl[_5v2wwg8*3 +2] = _nbzug0i[6]*_6fokzix[i*16 +8];
                    _5v2wwg8              += 1u;
                    
                }
                else if (_xkc5i8x[0] == 2u) 
                {   _p422ta8[i*4 +3]    = _hlv76mf +1u - _p422ta8[i*4 +2];
                    _6fokzix[i*16 +15]  = fabsf(_nbzug0i[5]);
                    
                    _nbzug0i[5]        = -1.0f*((_nbzug0i[3]/_nbzug0i[0]) *_6fokzix[i*16 +0]) / _6fokzix[i*16 +5];
                    _nbzug0i[6]        = -1.0f*((_nbzug0i[3]/_nbzug0i[0]) *_6fokzix[i*16 +1]) / _6fokzix[i*16 +6];
                    _nbzug0i[7]        = _6fokzix[i*16 +15]/sqrtf(_nbzug0i[5]*_nbzug0i[5] +_nbzug0i[6]*_nbzug0i[6]);
                    _nbzug0i[7]        = MIN(_nbzug0i[7], 1.0f);
                    _nbzug0i[5]       *= _nbzug0i[7];
                    _nbzug0i[6]       *= _nbzug0i[7];
                    
                    _6fokzix[i*16 +11]  = sqrtf(_nbzug0i[5]*_nbzug0i[5] +_nbzug0i[6]*_nbzug0i[6]);
                    _6fokzix[i*16 +12] += _6fokzix[i*16 +11];
                    _6fokzix[i*16 +13] += _nbzug0i[5];
                    _6fokzix[i*16 +14] += _nbzug0i[6];
                    _mybd2qe      += (_6fokzix[i*16 +11]*_6fokzix[i*16 +8]);
                    
                    _66cnolj[i][ (2u*(_hlv76mf - _p422ta8[i*4 +2])) +0] = _nbzug0i[5];
                    _66cnolj[i][ (2u*(_hlv76mf - _p422ta8[i*4 +2])) +1] = _nbzug0i[6];
                    
                    _6fokzix[i*16 +0]  += (_nbzug0i[5] *_6fokzix[i*16 +5]); 
                    _6fokzix[i*16 +1]  += (_nbzug0i[6] *_6fokzix[i*16 +6]); 
                    
            }   }
            
            MPI_Allgather(&_5v2wwg8, 1, MPI_UNSIGNED, _sosx43w, 1, MPI_INT, MPI_COMM_WORLD);
            _9lzs34j[0] = 0u;
            _nzrf5a7    = _sosx43w[0];
            for (i = 1; i < _0rmhrcn; i++)         {       _nzrf5a7 += _sosx43w[i];        _sosx43w[i-1] *= 3;           _9lzs34j[i] = _9lzs34j[i-1] + _sosx43w[i-1];         }
            _sosx43w[_0rmhrcn-1] *= 3;
            MPI_Allgatherv(_ozctwjl, _sosx43w[_fsqg8x8], MPI_FLOAT, _ny181ac, _sosx43w, _9lzs34j, MPI_FLOAT, MPI_COMM_WORLD);
            
            for (i = 0u; i < _lexdw6l[_fsqg8x8]; i++)
            {   
                _xkc5i8x[0] = (_nzrf5a7 <= _qnsn2s8[i])*1u  +  (_nzrf5a7 > _qnsn2s8[i])*2u;
                _xkc5i8x[0] = (_nzrf5a7 <= 0u)*0u            +  (_nzrf5a7 > 0u)*_xkc5i8x[0];
                
                if      (_xkc5i8x[0] == 1u)
                {   _nbzug0i[0] = 0.0f;   _nbzug0i[1] = 0.0f;   _nbzug0i[2] = 0.0f;
                    for (j = 0u; j < _nzrf5a7; j++) 
                    {   _xkc5i8x[2]  = (unsigned int)_jbc2sws[j*3 +0]; 
                        _xkc5i8x[1]  = _5gzpl1q[i][_xkc5i8x[2]] *2; 
                        _xkc5i8x[2]  = j*3;
                        _nbzug0i[0] += (_jbc2sws[_xkc5i8x[2] +1]*_ckeh5sw[i][_xkc5i8x[1] +0] + _jbc2sws[_xkc5i8x[2] +2]*_ckeh5sw[i][_xkc5i8x[1] +1] );
                        _nbzug0i[1] += (_jbc2sws[_xkc5i8x[2] +1]*_3ac6pe6[i][_xkc5i8x[1] +0] + _jbc2sws[_xkc5i8x[2] +2]*_3ac6pe6[i][_xkc5i8x[1] +1] );
                        _nbzug0i[2] += (_jbc2sws[_xkc5i8x[2] +1]*_obyghgk[i][_xkc5i8x[1] +0] + _jbc2sws[_xkc5i8x[2] +2]*_obyghgk[i][_xkc5i8x[1] +1] );
                    }
                    _6fokzix[i*16 +0] += _nbzug0i[0];
                    _6fokzix[i*16 +1] += _nbzug0i[1];
                    _6fokzix[i*16 +2] += _nbzug0i[2];
                }
                else if (_xkc5i8x[0] == 2u)
                {   memset(_t997e4u, 0, 2*_qnsn2s8[i]*sizeof(float) );
                    for (j = 0u; j < _nzrf5a7; j++) 
                    {   _xkc5i8x[2] = (unsigned int)_jbc2sws[j*3 +0]; 
                        _xkc5i8x[1] = _5gzpl1q[i][_xkc5i8x[2]] *2; 
                        _xkc5i8x[2]  = j*3;
                        _t997e4u[_xkc5i8x[1] +0] += _jbc2sws[_xkc5i8x[2] +1];
                        _t997e4u[_xkc5i8x[1] +1] += _jbc2sws[_xkc5i8x[2] +2];
                    }
                    _6fokzix[i*16 +0] += cblas_sdot(2*_qnsn2s8[i], _sh3tow4, 1, _ckeh5sw[i], 1);
                    _6fokzix[i*16 +1] += cblas_sdot(2*_qnsn2s8[i], _sh3tow4, 1, _3ac6pe6[i], 1);
                    _6fokzix[i*16 +2] += cblas_sdot(2*_qnsn2s8[i], _sh3tow4, 1, _obyghgk[i], 1);
            }   }
            _nzrf5a7 = (_hlv76mf < _atj9tze._at9mx2x-1)*_nzrf5a7  + (_hlv76mf >= _atj9tze._at9mx2x-1)*0; 
            
            _hlv76mf += 1u;
        }
        
        
        if (_vxpj4hn._5izolnn > 0u)
        {   
            _5v2wwg8 = 0u;
            for (i = 0u; i < _p2tixc7; i++)
            {   _xkc5i8x[0] =((fabsf(_6fokzix[_hx6zdcu[i]*16 +13]) + fabsf(_6fokzix[_hx6zdcu[i]*16 +14])) > 0.0f)*1u  +  ((fabsf(_6fokzix[_hx6zdcu[i]*16 +13]) + fabsf(_6fokzix[_hx6zdcu[i]*16 +14])) <= 0.0f)*0u;
                
                _ozctwjl[_5v2wwg8*3 +0] = (_xkc5i8x[0] == 0u)*_ozctwjl[_5v2wwg8*3 +0]  +  (_xkc5i8x[0] != 0u)*(float)(_1nlc8jp[_hx6zdcu[i]]);
                _ozctwjl[_5v2wwg8*3 +1] = (_xkc5i8x[0] == 0u)*_ozctwjl[_5v2wwg8*3 +1]  +  (_xkc5i8x[0] != 0u)*_6fokzix[_hx6zdcu[i]*16 +13]*_6fokzix[_hx6zdcu[i]*16 +8];
                _ozctwjl[_5v2wwg8*3 +2] = (_xkc5i8x[0] == 0u)*_ozctwjl[_5v2wwg8*3 +2]  +  (_xkc5i8x[0] != 0u)*_6fokzix[_hx6zdcu[i]*16 +14]*_6fokzix[_hx6zdcu[i]*16 +8];
                _5v2wwg8               = (_xkc5i8x[0] == 0u)*_5v2wwg8                +  (_xkc5i8x[0] != 0u)*(_5v2wwg8 + 1u);
            }
            
            MPI_Allgather(&_5v2wwg8, 1, MPI_UNSIGNED, _sosx43w, 1, MPI_INT, MPI_COMM_WORLD);
            _9lzs34j[0] = 0u;
            _nzrf5a7    = _sosx43w[0];
            for (i = 1; i < _0rmhrcn; i++)         {       _nzrf5a7 += _sosx43w[i];        _sosx43w[i-1] *= 3;           _9lzs34j[i] = _9lzs34j[i-1] + _sosx43w[i-1];         }
            _sosx43w[_0rmhrcn-1] *= 3;
            MPI_Allgatherv(_ozctwjl, _sosx43w[_fsqg8x8], MPI_FLOAT, _ny181ac, _sosx43w, _9lzs34j, MPI_FLOAT, MPI_COMM_WORLD);
            
            for (i = 0u; i < _1poxrte[_fsqg8x8]; i++)
            {   _xkc5i8x[0] = (_nzrf5a7 <= _iob0nlw[i])*1u  +  (_nzrf5a7 > _iob0nlw[i])*2u;
                _xkc5i8x[0] = (_nzrf5a7 <= 0u)*0u            +  (_nzrf5a7 > 0u)*_xkc5i8x[0];
                if (_xkc5i8x[0] == 1)
                {   _nbzug0i[0] = 0.0f;   _nbzug0i[1] = 0.0f;   _nbzug0i[2] = 0.0f;
                    for (j = 0u; j < _nzrf5a7; j++)
                    {   _xkc5i8x[2] = (unsigned int)_jbc2sws[j*3 +0]; 
                        _xkc5i8x[1] = _3b9so0u[i][_xkc5i8x[2]] *2; 
                        _xkc5i8x[2] = j*3;
                        
                        _nbzug0i[0] += (_jbc2sws[_xkc5i8x[2] +1]*_0ljch2v[i][_xkc5i8x[1] +0] + _jbc2sws[_xkc5i8x[2] +2]*_0ljch2v[i][_xkc5i8x[1] +1]);
                        _nbzug0i[1] += (_jbc2sws[_xkc5i8x[2]+ 1]*_0dtirc9[i][_xkc5i8x[1] +0] + _jbc2sws[_xkc5i8x[2] +2]*_0dtirc9[i][_xkc5i8x[1] +1]);
                        _nbzug0i[2] += (_jbc2sws[_xkc5i8x[2] +1]*_9dde3cu[i][_xkc5i8x[1] +0] + _jbc2sws[_xkc5i8x[2] +2]*_9dde3cu[i][_xkc5i8x[1] +1]);
                    }
                    _7yt2m6t[i*16 +0] += _nbzug0i[0];
                    _7yt2m6t[i*16 +1] += _nbzug0i[1];
                    _7yt2m6t[i*16 +2] += _nbzug0i[2];
                }
                else if (_xkc5i8x[0] == 2)
                {   memset(_t997e4u, 0, 2*_iob0nlw[i]*sizeof(float) );
                    for (j = 0u; j < _nzrf5a7; j++) 
                    {   _xkc5i8x[2] = (unsigned int)_jbc2sws[j*3 +0]; 
                        _xkc5i8x[1] = _3b9so0u[i][_xkc5i8x[2]]; 
                        _xkc5i8x[2] = j*3;
                        _t997e4u[_xkc5i8x[1] +0] += _jbc2sws[_xkc5i8x[2] +1];
                        _t997e4u[_xkc5i8x[1] +1] += _jbc2sws[_xkc5i8x[2] +2];
                    }
                    _7yt2m6t[i*16 +0] +=  cblas_sdot(2*_iob0nlw[i], _sh3tow4, 1, _0ljch2v[i], 1);
                    _7yt2m6t[i*16 +1] +=  cblas_sdot(2*_iob0nlw[i], _sh3tow4, 1, _0dtirc9[i], 1);
                    _7yt2m6t[i*16 +2] +=  cblas_sdot(2*_iob0nlw[i], _sh3tow4, 1, _9dde3cu[i], 1);
                }
            }
        } 
        
        
        _t9nm3v2 = _p2tixc7;
        
        MPI_Allreduce(MPI_IN_PLACE, &_t9nm3v2,  1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, &_mybd2qe, 1, MPI_FLOAT,    MPI_SUM, MPI_COMM_WORLD);
        
        
        _mybd2qe *= _atj9tze._kkrzd5e;
        _qxmcthp    = (log10f(_mybd2qe) -9.1f)/1.5f;
        
        
        if (_t9nm3v2 >= _vxpj4hn._e6zvlhy)
        {   int _gg24479 = 0;

            _xo100vh += 1u;
            _gvyie6i += 1u;
            memset(_hj3muur, 0, 6*sizeof(float));

            for (i = 0u; i < _p2tixc7; i++)
            {   _xkc5i8x[1] = _hx6zdcu[i];
                _nbzug0i[0] = sqrtf(_6fokzix[_hx6zdcu[i]*16 +13]*_6fokzix[_hx6zdcu[i]*16 +13] +_6fokzix[_hx6zdcu[i]*16 +14]*_6fokzix[_hx6zdcu[i]*16 +14]);
                
                _hj3muur[0] += (_z39pecn[_hx6zdcu[i]*8 +0]*_nbzug0i[0]*_6fokzix[_hx6zdcu[i]*16 +8]); 
                _hj3muur[1] += (_z39pecn[_hx6zdcu[i]*8 +1]*_nbzug0i[0]*_6fokzix[_hx6zdcu[i]*16 +8]);
                _hj3muur[2] += (_z39pecn[_hx6zdcu[i]*8 +2]*_nbzug0i[0]*_6fokzix[_hx6zdcu[i]*16 +8]);
                _hj3muur[3] += _6fokzix[_hx6zdcu[i]*16 +8]; 
                _hj3muur[4] += _nbzug0i[0];
                _hj3muur[5] += sqrtf(_d1vn8lf[_hx6zdcu[i]*4 +0]*_d1vn8lf[_hx6zdcu[i]*4 +0] + _d1vn8lf[_hx6zdcu[i]*4 +1]*_d1vn8lf[_hx6zdcu[i]*4 +1]) - sqrtf(_6fokzix[_hx6zdcu[i]*16 +0]*_6fokzix[_hx6zdcu[i]*16 +0] + _6fokzix[_hx6zdcu[i]*16 +1]*_6fokzix[_hx6zdcu[i]*16 +1]); 
                
                _pj8t22c[i]      = _1nlc8jp[_hx6zdcu[i]]; 
                _50b8y5d[i]     = _p422ta8[_hx6zdcu[i]*4 +2];
                _dtg5nqj[i]    = _6fokzix[_hx6zdcu[i]*16 +13];
                _kc4lorx[i]    = _6fokzix[_hx6zdcu[i]*16 +14];
                _pkxvhag[i]    = _z39pecn[_hx6zdcu[i]*8 +3] + sqrtf(_dtg5nqj[i]*_dtg5nqj[i] + _kc4lorx[i]*_kc4lorx[i]); 
                _qijc1yd[i]     = sqrtf(_d1vn8lf[_hx6zdcu[i]*4 +0]*_d1vn8lf[_hx6zdcu[i]*4 +0] + _d1vn8lf[_hx6zdcu[i]*4 +1]*_d1vn8lf[_hx6zdcu[i]*4 +1]) - sqrtf(_6fokzix[_hx6zdcu[i]*16 +0]*_6fokzix[_hx6zdcu[i]*16 +0] + _6fokzix[_hx6zdcu[i]*16 +1]*_6fokzix[_hx6zdcu[i]*16 +1]); 
                _gg24479   += (_p422ta8[_hx6zdcu[i]*4 +1] == 2)*1  +  (_p422ta8[_hx6zdcu[i]*4 +1] != 2)*0;
            }
            
            MPI_Allreduce(MPI_IN_PLACE, _hj3muur,      6,            MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, &_gg24479,    1,            MPI_INT,   MPI_SUM, MPI_COMM_WORLD);
            
            _hj3muur[0] /= (_mybd2qe/_atj9tze._kkrzd5e);              _hj3muur[1] /= (_mybd2qe/_atj9tze._kkrzd5e);            _hj3muur[2] /= (_mybd2qe/_atj9tze._kkrzd5e);
            _hj3muur[4] /= (float)_t9nm3v2;      _hj3muur[5] /= (float)_t9nm3v2;
            
            MPI_Allgather(&_p2tixc7, 1, MPI_UNSIGNED, _sosx43w, 1, MPI_INT, MPI_COMM_WORLD);
            
            _9lzs34j[0] = 0;
            for (i = 1; i < _0rmhrcn; i++)         {       _9lzs34j[i]  = _9lzs34j[i-1] + _sosx43w[i-1];         }
            
            MPI_Gatherv(_pj8t22c,   _sosx43w[_fsqg8x8], MPI_UNSIGNED, _zk94yxf, _sosx43w, _9lzs34j, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
            MPI_Gatherv(_50b8y5d,  _sosx43w[_fsqg8x8], MPI_UNSIGNED, _iblejgu, _sosx43w, _9lzs34j, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
            MPI_Gatherv(_dtg5nqj, _sosx43w[_fsqg8x8], MPI_FLOAT,    _2ljk4q5, _sosx43w, _9lzs34j, MPI_FLOAT,    0, MPI_COMM_WORLD);
            MPI_Gatherv(_kc4lorx, _sosx43w[_fsqg8x8], MPI_FLOAT,    _v7c6s5n, _sosx43w, _9lzs34j, MPI_FLOAT,    0, MPI_COMM_WORLD);
            MPI_Gatherv(_pkxvhag, _sosx43w[_fsqg8x8], MPI_FLOAT,    _kb4ehou, _sosx43w, _9lzs34j, MPI_FLOAT,    0, MPI_COMM_WORLD);
            MPI_Gatherv(_qijc1yd,  _sosx43w[_fsqg8x8], MPI_FLOAT,    _6qlz2rw, _sosx43w, _9lzs34j, MPI_FLOAT,    0, MPI_COMM_WORLD);
            
            if (_fsqg8x8 == 0)
            {   _dxj37c9 = (_qxmcthp > _dxj37c9)*_qxmcthp  +  (_qxmcthp <= _dxj37c9)*_dxj37c9;
                if (_vxpj4hn._qdbxvh7 == 1u) { fprintf(stdout,"%7u  Time: %9.4lf  (%7.2lf days since last)   intSeisLoops: %2u   Element#: %6u (%3u broken)    RA: %8.2f   mSlip: %7.3f      MRF %5u     M: %4.2f      M_max: %4.2f\n", _xo100vh, _z36u8ue, (_z36u8ue-_ihx74d0)*365.25, _1o2rntj, _t9nm3v2, _gg24479, (_hj3muur[3]*1.0E-6f), _hj3muur[4], (_hlv76mf-1u), _qxmcthp, _dxj37c9);       }
                
                fseek(_tajw253, 0, SEEK_END);
                fwrite(&_z36u8ue,    sizeof(double),              1, _tajw253);
                fwrite(&_mybd2qe,    sizeof(float),               1, _tajw253);
                fwrite(&_qxmcthp,      sizeof(float),               1, _tajw253);
                fwrite(&_hlv76mf,   sizeof(unsigned int),        1, _tajw253);
                fwrite(_mrunjdu,      sizeof(float),               3, _tajw253);
                fwrite(_hj3muur,   sizeof(float),               6, _tajw253);
                
                fwrite(&_t9nm3v2,     sizeof(unsigned),            1, _tajw253);
                fwrite(_zk94yxf,    sizeof(unsigned int), _t9nm3v2, _tajw253);
                fwrite(_iblejgu,    sizeof(unsigned int), _t9nm3v2, _tajw253);
                fwrite(_2ljk4q5,    sizeof(float),        _t9nm3v2, _tajw253);
                fwrite(_v7c6s5n,    sizeof(float),        _t9nm3v2, _tajw253);
                fwrite(_kb4ehou,    sizeof(float),        _t9nm3v2, _tajw253);
                fwrite(_6qlz2rw,    sizeof(float),        _t9nm3v2, _tajw253);
                fseek(_tajw253, 0, SEEK_SET);
                fwrite(&_xo100vh,      sizeof(unsigned int),        1, _tajw253);
            }
            
            _gjyv5g1 += (_xo100vh > 1u)*(_z36u8ue - _ihx74d0) + (_xo100vh <= 1u)*0.0;
            _ihx74d0     = _z36u8ue; 
        } 
        MPI_Barrier( MPI_COMM_WORLD );
        
        
        if ((_vxpj4hn._wk8temp == 1u) && (_qxmcthp >= _vxpj4hn._rdihsa7))
        {   float *_q52wspo = (float *) malloc( _atj9tze._at9mx2x *sizeof(float));
            float *_fv177ud = (float *) malloc( _atj9tze._at9mx2x *sizeof(float));
            long long int _8jvh7uj = 0;
            long long int _xphx2kb = 0;
            char _i18k0ax[512];
            
            for (i = 0; i < _lexdw6l[_fsqg8x8]; i++)
            {   _xphx2kb = (_p422ta8[i*4 +1] == 0u)*_xphx2kb  +  (_p422ta8[i*4 +1] != 0u)*(_xphx2kb + 2*sizeof(float) +4*sizeof(unsigned int) +2*_p422ta8[i*4 +3]*sizeof(float) );
            }
            
            MPI_Allgather(&_xphx2kb, 1, MPI_LONG_LONG_INT, _vsigbia, 1, MPI_LONG_LONG_INT, MPI_COMM_WORLD);
            _f92movs[0] = 0;
            for (i = 1; i < _0rmhrcn; i++)         {       _f92movs[i]  = _f92movs[i-1]  +_vsigbia[i-1];         }
            
            strcpy(_i18k0ax,_dygflij);       strcat(_i18k0ax,"_t");         snprintf(_75aud9k,sizeof(_75aud9k), "%lf",_z36u8ue);      strcat(_i18k0ax,_75aud9k);  strcat(_i18k0ax,"_M");         snprintf(_75aud9k,sizeof(_75aud9k), "%f",_qxmcthp);      strcat(_i18k0ax,_75aud9k);           strcat(_i18k0ax,".srfb");
            
            MPI_File_delete(_i18k0ax, MPI_INFO_NULL);
            MPI_File_open(MPI_COMM_WORLD, _i18k0ax, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &_npuevjl);
            
            if (_fsqg8x8 == 0)
            {   int _3zezqac = (int)VersionNumber;
                MPI_File_write_at(_npuevjl, _8jvh7uj,  &_3zezqac,   1, MPI_INT,       &_0get8t6);          _8jvh7uj += 1*sizeof(int);
                MPI_File_write_at(_npuevjl, _8jvh7uj,  &_qxmcthp,         1, MPI_FLOAT,     &_0get8t6);          _8jvh7uj += 1*sizeof(float);
                _nbzug0i[0] = (float)_z36u8ue;
                MPI_File_write_at(_npuevjl, _8jvh7uj,  &_nbzug0i[0],        1, MPI_FLOAT,     &_0get8t6);          _8jvh7uj += 1*sizeof(float);
                MPI_File_write_at(_npuevjl, _8jvh7uj,  &_t9nm3v2,        1, MPI_UNSIGNED,  &_0get8t6);          _8jvh7uj += 1*sizeof(unsigned int);
                MPI_File_write_at(_npuevjl, _8jvh7uj,  &_atj9tze._mzwffhh,  1, MPI_FLOAT,     &_0get8t6);          _8jvh7uj += 1*sizeof(float);
                MPI_File_write_at(_npuevjl, _8jvh7uj,  &_atj9tze._douadvw, 1, MPI_FLOAT,     &_0get8t6);          _8jvh7uj += 1*sizeof(float);
                MPI_File_write_at(_npuevjl, _8jvh7uj,  &_atj9tze._cmktlzk,   1, MPI_FLOAT,     &_0get8t6);          _8jvh7uj += 1*sizeof(float);
            }
            
            _8jvh7uj = 1*sizeof(int) +1*sizeof(unsigned int) +5*sizeof(float); 
            
            for (i = 0; i < _lexdw6l[_fsqg8x8]; i++)
            {   if (_p422ta8[i*4 +1] != 0u) 
                {   _xkc5i8x[0] = _1nlc8jp[i];
                    MPI_File_write_at(_npuevjl, (_f92movs[_fsqg8x8]+_8jvh7uj),  &_xkc5i8x[0],       1, MPI_UNSIGNED,  &_0get8t6);           _8jvh7uj += 1*sizeof(unsigned int);
                    MPI_File_write_at(_npuevjl, (_f92movs[_fsqg8x8]+_8jvh7uj),  &_p422ta8[i*4 +2],   1, MPI_UNSIGNED,  &_0get8t6);           _8jvh7uj += 1*sizeof(unsigned int);
                    
                    MPI_File_write_at(_npuevjl, (_f92movs[_fsqg8x8]+_8jvh7uj),  &_6fokzix[i*16 +13], 1, MPI_FLOAT,     &_0get8t6);           _8jvh7uj += 1*sizeof(float);
                    MPI_File_write_at(_npuevjl, (_f92movs[_fsqg8x8]+_8jvh7uj),  &_p422ta8[i*4 +3],   1, MPI_UNSIGNED,  &_0get8t6);           _8jvh7uj += 1*sizeof(unsigned int);
                    MPI_File_write_at(_npuevjl, (_f92movs[_fsqg8x8]+_8jvh7uj),  &_6fokzix[i*16 +14], 1, MPI_FLOAT,     &_0get8t6);           _8jvh7uj += 1*sizeof(float);
                    MPI_File_write_at(_npuevjl, (_f92movs[_fsqg8x8]+_8jvh7uj),  &_p422ta8[i*4 +3],   1, MPI_UNSIGNED,  &_0get8t6);           _8jvh7uj += 1*sizeof(unsigned int);
                    
                    for (j = 0; j < _p422ta8[i*4 +3]; j++)       {        _q52wspo[j] = _66cnolj[i][2u*j +0];         _fv177ud[j] = _66cnolj[i][2u*j +1];         }
                    MPI_File_write_at(_npuevjl, (_f92movs[_fsqg8x8]+_8jvh7uj), _q52wspo, _p422ta8[i*4 +3], MPI_FLOAT, &_0get8t6);           _8jvh7uj += _p422ta8[i*4 +3]*sizeof(float);
                    MPI_File_write_at(_npuevjl, (_f92movs[_fsqg8x8]+_8jvh7uj), _fv177ud, _p422ta8[i*4 +3], MPI_FLOAT, &_0get8t6);           _8jvh7uj += _p422ta8[i*4 +3]*sizeof(float);
            }   }
            
            free(_q52wspo);   free(_fv177ud);
            MPI_Barrier( MPI_COMM_WORLD );
            MPI_File_close(&_npuevjl);
            
            double _eimvihc;
            clock_gettime(CLOCK_REALTIME, &_s9bp4m1);
            _eimvihc = (double)_0rmhrcn *(_s9bp4m1.tv_sec - _04k0lth.tv_sec) + (_s9bp4m1.tv_nsec - _04k0lth.tv_nsec)/1.0E+9;
            if (_fsqg8x8 == 0)         {       fprintf(stdout,"CPU time (hours) for STF event %s:  %lf\n",_i18k0ax, _eimvihc/3600.0);              }
        } 
        
        
        for (i = 0u; i < _p2tixc7; i++)          {       memset(_66cnolj[_hx6zdcu[i]],  0, 2*_atj9tze._at9mx2x*sizeof(float));            }
        
        if (_vxpj4hn._pwh0mx1 == 1u)
        {   
            for (i = 0u; i < _p2tixc7; i++)
            {   
                _myx444n(&_atj9tze, &_vxpj4hn, _hx6zdcu[i], _p422ta8, _6fokzix, _0dl5e11, 0);
        }   }
        
        _790pdw5(_fsqg8x8, _lexdw6l, &_atj9tze, _p422ta8, _6fokzix, _0dl5e11, _z39pecn);
        
        
        clock_gettime(CLOCK_REALTIME, &_o6l1f90);
        _jlpefja += (_o6l1f90.tv_sec - _rw6tm2o.tv_sec) + (_o6l1f90.tv_nsec - _rw6tm2o.tv_nsec)/1.0E+9;
    }

    clock_gettime(CLOCK_REALTIME, &_1zocfi4);
    _qydbpgv = (_1zocfi4.tv_sec - _glmy62x.tv_sec) + (_1zocfi4.tv_nsec - _glmy62x.tv_nsec)/1.0E+9;
    
    if (_fsqg8x8 == 0)
    {   
        _gjyv5g1 /= (double)_xo100vh;
        fprintf(stdout,"\nTotal RunTime for EQcycle (minutes): %6.4f       and mean inter-event time (days): %3.2lf\n",_qydbpgv/60.0f, _gjyv5g1*365.25 );
        fprintf(stdout,"Time spent in interseismic: %4.2lf   Time spent in coseismic: %4.2lf\n",_wshmlfz/60.0f, _jlpefja/60.0f);
        

        fprintf(stdout,"\n average iteration number: %f\n",(float)_w2pcn9t*_owgzrix/_gvyie6i);
        fclose(_tajw253);
    }
    
    
    
    MPI_Barrier( MPI_COMM_WORLD );





























































































    
    
    free(_r7vgrf1);   free(_lu6k7p9);   free(_p422ta8);   free(_6fokzix);   free(_0dl5e11);   free(_z39pecn);   free(_7yt2m6t);   free(_t997e4u);   free(_zkwdeek);   free(_zk94yxf);   free(_iblejgu);
    free(_ozctwjl);   free(_ny181ac);   free(_lqbmae3);   free(_gyspeng);   free(_d1vn8lf);   free(_1vsmado);   free(_pj8t22c);   free(_dtg5nqj);   free(_kc4lorx);   free(_pkxvhag);
    free(_2ljk4q5);   free(_v7c6s5n);   free(_kb4ehou);   free(_6qlz2rw);   free(_hx6zdcu);   free(_50b8y5d);   free(_qijc1yd);
    
    for (i = 0; i < _lexdw6l[_fsqg8x8]; i++)
    {   free(_66cnolj[i]);
        free(_weehwjq[i]);   free(_2292mtr[i]);   free(_3lbqgu3[i]);   free(_5bwikot[i]);
    }
    if (_vxpj4hn._5izolnn > 0u)
    {   for (i = 0; i < _lexdw6l[_fsqg8x8]; i++)
        {   free(_jy62zee[i]);   free(_iq2fa39[i]);   free(_b4ldims[i]);
        }
        for (i = 0; i < _1poxrte[_fsqg8x8]; i++)
        {   free(_c795rlq[i]);   free(_swbkok3[i]);   free(_m458ax4[i]);   free(_ffqot8h[i]);   free(_s6g298c[i]);   free(_zd76szn[i]);   free(_ukwqvr1[i]);   free(_p8ijl2z[i]);
        }
    }
    
    free(_e05qmkt);    free(_avt8k06);   free(_qwlvunl);   free(_o2t8fy6);   free(_66cnolj);
    free(_nn16uxf); free(_weehwjq); free(_2292mtr); free(_3lbqgu3); free(_5bwikot);
    free(_q0e27ze); free(_jy62zee); free(_iq2fa39); free(_b4ldims);
    free(_zn4regb); free(_c795rlq); free(_m458ax4); free(_s6g298c); free(_ukwqvr1);
    free(_0vs52sq); free(_swbkok3); free(_ffqot8h); free(_zd76szn); free(_p8ijl2z);
    
    
    MPI_Barrier( MPI_COMM_WORLD );
    MPI_Finalize( );
    return 0;
}









