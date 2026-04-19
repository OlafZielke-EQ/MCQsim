#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
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


extern void    _id4n1u3(float *restrict _j49yd82, const float *restrict _28xg8y7, const float *restrict _kf1622v);
extern void     _j2qp89u(float *restrict _iquwrxz, const float *restrict _f445cri, const float *restrict _8bixoiq);
extern float    _mabsi9n(const float *restrict _n42x2ez);
extern void       _bfuk9e1(float *restrict _n42x2ez);
extern void _ktakpux(const float *restrict _088o98b, float *restrict _j3swebn, float *restrict _795e5iz);
extern void     _ftu2gc5(float *restrict _0i9g5dz, const float *restrict _mt8vt43, const float *restrict _xyzx9jp);
extern float            _53ayeqa(long *_rrhf4ab);
extern float            _hthsjr0(long *_rrhf4ab);

void _l4zkloo(unsigned int _6zjq7w7, const unsigned int *restrict _r59c81d, float *restrict _yn4u9bx, const float *restrict _mszdxim, float *restrict _9el8hl8, float *restrict _y9x8syl);
void _kep7y6w(int _fsqg8x8, unsigned int _7c6qc6k, float *restrict _bjfvyb6, unsigned int _6zjq7w7, unsigned int *restrict _r59c81d, float *restrict _yn4u9bx, const float *restrict _mszdxim, const unsigned int *restrict _lx2toow, const int *restrict _ov6tt9q);


void   _qjiujv6(int _fsqg8x8, struct _5r5uqpg *restrict _atj9tze, struct _xjuy0bv *restrict _vxpj4hn, struct _j49itn7 *restrict _s55j2o7, const int *restrict _lexdw6l, const int *restrict _1poxrte, const unsigned int *restrict _1nlc8jp, const unsigned int *restrict _fftouhy, float *restrict _e9ftaxp, float *restrict _9yo18pa, unsigned int *restrict _vmefozt, unsigned int *restrict _rgb93wk, float *restrict _4xkss0k, float *restrict _rrpzlk3, const float *restrict _zzf8s2z, const float *restrict _3n3bzt9, float *restrict _6fokzix, float *restrict _0dl5e11, float *restrict _z39pecn, float *restrict _7yt2m6t)
{  
    unsigned int i;
    float _9el8hl8,   _y9x8syl,   _nbzug0i;
    float _088o98b[3];
    float _ge2j51b[6],   _0i9g5dz[6];
    
    _9el8hl8 = 0.0f,      _y9x8syl = 0.0f;
    _l4zkloo(_vxpj4hn->_3ru7myo, _vmefozt, _4xkss0k, _zzf8s2z, &_9el8hl8, &_y9x8syl);
    _s55j2o7->_b68rerv    = _y9x8syl; 
    _s55j2o7->_xf7jxjn  = _9el8hl8; 
    
    _9el8hl8 = 0.0f,      _y9x8syl = 0.0f;
    _l4zkloo(_vxpj4hn->_5izolnn, _rgb93wk, _rrpzlk3, _3n3bzt9, &_9el8hl8, &_y9x8syl);
    _s55j2o7->_gfu8iog    = _y9x8syl;
    _s55j2o7->_z2c9ros  = _9el8hl8;
    
    _s55j2o7->_ycz04cd    = _s55j2o7->_b68rerv*1.0E-4f; 
    _s55j2o7->_cb1kjtg    = (_s55j2o7->_z2c9ros/_s55j2o7->_xf7jxjn)*_s55j2o7->_ycz04cd;
    
    
    _atj9tze->_kkrzd5e    = _atj9tze->_kkrzd5e*1.0E+9f; 
    _atj9tze->_6snglh6    = (2.0f*_atj9tze->_kkrzd5e*_atj9tze->_wrb7f6f)/(1.0f-2.0f*_atj9tze->_wrb7f6f); 
    _nbzug0i               = (_atj9tze->_cmktlzk > 0.0f)*_atj9tze->_cmktlzk + (_atj9tze->_cmktlzk <= 0.0f)*2700.0f; 
    _atj9tze->_douadvw    = sqrtf(_atj9tze->_kkrzd5e/_nbzug0i); 
    _atj9tze->_mzwffhh     = (1.0f*_s55j2o7->_b68rerv) /(sqrtf((_atj9tze->_6snglh6 +2.0f*_atj9tze->_kkrzd5e)/_nbzug0i));
    _vxpj4hn->_7odkef8 = _vxpj4hn->_7odkef8*_atj9tze->_mzwffhh;
    
    
    _kep7y6w(_fsqg8x8, _s55j2o7->_6wmdw0f, _e9ftaxp, _vxpj4hn->_3ru7myo,  _vmefozt, _4xkss0k, _zzf8s2z, _1nlc8jp, _lexdw6l);
    _kep7y6w(_fsqg8x8, _s55j2o7->_rp9dxoq, _9yo18pa, _vxpj4hn->_5izolnn,  _rgb93wk, _rrpzlk3, _3n3bzt9, _fftouhy, _1poxrte);
    
    for (i = 0u; i < _lexdw6l[_fsqg8x8]; i++) 
    {   
        _6fokzix[i*16 +8]  = _4xkss0k[_1nlc8jp[i]*16 +15];
        _6fokzix[i*16 +2]  = _atj9tze->_cmktlzk * 9.81f * _4xkss0k[_1nlc8jp[i]*16 +2] +(_atj9tze->_iwvhdk8*1.0E+6f); 
        _nbzug0i          = _53ayeqa(&_atj9tze->_0cdk36m);      _nbzug0i = _nbzug0i -0.5f; 
        _nbzug0i          = _atj9tze->_cmktlzk * 9.81f *_s55j2o7->_b68rerv*_nbzug0i; 
        _6fokzix[i*16 +2]  = _6fokzix[i*16 +2] +_nbzug0i; 

        _0dl5e11[i*16 +0]  = _6fokzix[i*16 +2]; 
        _z39pecn[i*8 +0] = _4xkss0k[_1nlc8jp[i]*16 +0];
        _z39pecn[i*8 +1] = _4xkss0k[_1nlc8jp[i]*16 +1];
        _z39pecn[i*8 +2] = _4xkss0k[_1nlc8jp[i]*16 +2];
        
        if (_vmefozt[_1nlc8jp[i]*8 +5] == 1u) 
        {   _6fokzix[i*16 +9]  = 0.0f;
            _6fokzix[i*16 +10] = 0.0f;
            _6fokzix[i*16 +13] = 0.001f*_4xkss0k[_1nlc8jp[i]*16 +4] *cosf(_4xkss0k[_1nlc8jp[i]*16 +5]*(M_PI/180.0f));
            _6fokzix[i*16 +14] = 0.001f*_4xkss0k[_1nlc8jp[i]*16 +4] *sinf(_4xkss0k[_1nlc8jp[i]*16 +5]*(M_PI/180.0f));
        } 
        else
        {   _088o98b[0]      = sinf(_4xkss0k[_1nlc8jp[i]*16 +5]*(M_PI/180.0f));
            _088o98b[1]      = cosf(_4xkss0k[_1nlc8jp[i]*16 +5]*(M_PI/180.0f));
            _088o98b[2]      = 0.0f;
            _ge2j51b[0] = 1.0E+6f*_4xkss0k[_1nlc8jp[i]*16 +4] *_088o98b[0]*_088o98b[0]; 
            _ge2j51b[1] = 1.0E+6f*_4xkss0k[_1nlc8jp[i]*16 +4] *_088o98b[0]*_088o98b[1]; 
            _ge2j51b[2] = 1.0E+6f*_4xkss0k[_1nlc8jp[i]*16 +4] *_088o98b[0]*_088o98b[2]; 
            _ge2j51b[3] = 1.0E+6f*_4xkss0k[_1nlc8jp[i]*16 +4] *_088o98b[1]*_088o98b[1]; 
            _ge2j51b[4] = 1.0E+6f*_4xkss0k[_1nlc8jp[i]*16 +4] *_088o98b[1]*_088o98b[2]; 
            _ge2j51b[5] = 1.0E+6f*_4xkss0k[_1nlc8jp[i]*16 +4] *_088o98b[2]*_088o98b[2]; 
            
            _ftu2gc5(_0i9g5dz, &_4xkss0k[_1nlc8jp[i]*16 +6], _ge2j51b);
            
            _6fokzix[i*16 +9]   = _0i9g5dz[0]; 
            _6fokzix[i*16 +10]  = _0i9g5dz[1]; 
            _6fokzix[i*16 +13]  = 0.0f;
            _6fokzix[i*16 +14]  = 0.0f;
        }
    }
    
    for (i = 0u; i < _1poxrte[_fsqg8x8]; i++)
    {   _7yt2m6t[i*16 +8]  = _rrpzlk3[_fftouhy[i]*16 +15];
        _7yt2m6t[i*16 +9]  = _rrpzlk3[_fftouhy[i]*16 +3];
        _7yt2m6t[i*16 +10] = _rrpzlk3[_fftouhy[i]*16 +4];
        _7yt2m6t[i*16 +11] = _rrpzlk3[_fftouhy[i]*16 +5];
    }
    
    return;
}



void _l4zkloo(unsigned int _6zjq7w7, const unsigned int *restrict _r59c81d, float *restrict _yn4u9bx, const float *restrict _mszdxim, float *restrict _9el8hl8, float *restrict _y9x8syl)
{   unsigned int i;
    float _nbzug0i;
    float _dan085r[3],   _ntdr5xy[3],   _fc4kl8c[3],   _xrum2lh[3],   _gcsotdh[3],   _hlmh1t9[3],   _39ayq6w[3],   _w8yl68n[3],   _088o98b[3],   _n42x2ez[3];
    
    for (i = 0u; i < _6zjq7w7; i++)
    {   
        memcpy(_dan085r, &_mszdxim[_r59c81d[i*8 +0]*4 +0], 3*sizeof(float));
        memcpy(_ntdr5xy, &_mszdxim[_r59c81d[i*8 +1]*4 +0], 3*sizeof(float));
        memcpy(_fc4kl8c, &_mszdxim[_r59c81d[i*8 +2]*4 +0], 3*sizeof(float));
        
        _id4n1u3(_xrum2lh, _ntdr5xy, _dan085r);                  _id4n1u3(_gcsotdh, _fc4kl8c, _dan085r);              _id4n1u3(_hlmh1t9, _fc4kl8c, _ntdr5xy);
        _id4n1u3(_39ayq6w, &_yn4u9bx[i*16 +0], _dan085r);    _id4n1u3(_w8yl68n, &_yn4u9bx[i*16 +0], _ntdr5xy);
        
        _j2qp89u(_n42x2ez,_xrum2lh, _39ayq6w);         _nbzug0i = 0.5f*_mabsi9n(_n42x2ez);          *_y9x8syl += ((2.0f*_nbzug0i)/_mabsi9n(_xrum2lh));
        _j2qp89u(_n42x2ez,_gcsotdh, _39ayq6w);         _nbzug0i = 0.5f*_mabsi9n(_n42x2ez);          *_y9x8syl += ((2.0f*_nbzug0i)/_mabsi9n(_gcsotdh));
        _j2qp89u(_n42x2ez,_hlmh1t9, _w8yl68n);         _nbzug0i = 0.5f*_mabsi9n(_n42x2ez);          *_y9x8syl += ((2.0f*_nbzug0i)/_mabsi9n(_hlmh1t9));
        
        _j2qp89u(_088o98b,_xrum2lh, _gcsotdh);              _yn4u9bx[i*16 +15] = 0.5f*_mabsi9n(_088o98b); 
        
        *_9el8hl8 += _yn4u9bx[i*16 +15];
        
    }
    *_y9x8syl = (_6zjq7w7 > 0u) ? *_y9x8syl/((float)_6zjq7w7*3.0f) : *_y9x8syl;
    *_9el8hl8  = (_6zjq7w7 > 0u) ? *_9el8hl8/  (float)_6zjq7w7       : *_9el8hl8;
    
    return;
}



void _kep7y6w(int _fsqg8x8, unsigned int _7c6qc6k, float *restrict _bjfvyb6, unsigned int _6zjq7w7, unsigned int *restrict _r59c81d, float *restrict _yn4u9bx, const float *restrict _mszdxim, const unsigned int *restrict _lx2toow, const int *restrict _ov6tt9q)
{   
    unsigned int i,   j,   _06ftlo4;
    unsigned int _nezw9he,   _omxlvz4,   _qke2gjr,   _4d4r9v4,   _k00ilu6,   _vfua9bq,   _1btr23q,   _iua08ek;
    float _nbzug0i = 0.0f,   _t26ys44,   _c7fov2y;
    float _dan085r[3],   _ntdr5xy[3],   _fc4kl8c[3],   _gn0t6dv[3],   _xrum2lh[3],   _gcsotdh[3],   _hlmh1t9[3],   _39ayq6w[3],   _w8yl68n[3],   _088o98b[3],   _j3swebn[3],   _795e5iz[3],   _n42x2ez[3];
    struct { float _2c9t199;  int _gv44ftj;  } _g8vuazx[1], _yllwb7s[1];
    
    for (_nezw9he = 0u; _nezw9he < _7c6qc6k; _nezw9he++)
    {   
        _4d4r9v4   = 0u;
        for (i = 0u; i < _6zjq7w7; i++)         {   if (_r59c81d[i*8 +3] == _nezw9he)       {   _4d4r9v4 = i;        break;          }              } 
        
        _g8vuazx[0]._gv44ftj = 0;                 _g8vuazx[0]._2c9t199 = 0.0f;
        _06ftlo4        = 0u;                 _c7fov2y        = 0.0f;
        for (i = 0u; i < _ov6tt9q[_fsqg8x8]; i++)
        {   if (_r59c81d[_lx2toow[i]*8 +3] == _nezw9he)
            {   for (j = 0u; j < 3u; j++)
                {   _id4n1u3(_n42x2ez, &_yn4u9bx[_4d4r9v4*16 +0], &_mszdxim[ _r59c81d[_lx2toow[i]*8 +j]*4 +0]);
                    _nbzug0i   = _mabsi9n(_n42x2ez);
                    _06ftlo4  = (_nbzug0i < _c7fov2y)*_06ftlo4   +  (_nbzug0i >= _c7fov2y)*_r59c81d[_lx2toow[i]*8 +j];
                    _c7fov2y = (_nbzug0i < _c7fov2y)*_c7fov2y  +  (_nbzug0i >= _c7fov2y)*_nbzug0i;
        }   }   }
        _g8vuazx[0]._gv44ftj = _06ftlo4;            _g8vuazx[0]._2c9t199 = _c7fov2y;
        
        MPI_Allreduce(_g8vuazx, _yllwb7s, 1, MPI_FLOAT_INT, MPI_MAXLOC, MPI_COMM_WORLD);
        
        _k00ilu6 = (unsigned int)_yllwb7s[0]._gv44ftj;
        
        _g8vuazx[0]._gv44ftj = 0;                 _g8vuazx[0]._2c9t199 = 0.0f;
        _06ftlo4         = 0u;                _c7fov2y        = 0.0f;
        for (i = 0u; i < _ov6tt9q[_fsqg8x8]; i++)
        {   if (_r59c81d[_lx2toow[i]*8 +3] == _nezw9he)
            {   for (j = 0u; j < 3u; j++)
                {   _id4n1u3(_n42x2ez, &_mszdxim[_k00ilu6*4 +0], &_mszdxim[ _r59c81d[_lx2toow[i]*8 +j]*4 +0]);
                    _nbzug0i   = _mabsi9n(_n42x2ez);
                    _06ftlo4  = (_nbzug0i < _c7fov2y)*_06ftlo4   +  (_nbzug0i >= _c7fov2y)*_r59c81d[_lx2toow[i]*8 +j];
                    _c7fov2y = (_nbzug0i < _c7fov2y)*_c7fov2y  +  (_nbzug0i >= _c7fov2y)*_nbzug0i;
        }   }   }
        _g8vuazx[0]._gv44ftj = _06ftlo4;      _g8vuazx[0]._2c9t199 = _c7fov2y;
        
        MPI_Allreduce(_g8vuazx, _yllwb7s, 1, MPI_FLOAT_INT, MPI_MAXLOC, MPI_COMM_WORLD);
        
        _vfua9bq = (unsigned int)_yllwb7s[0]._gv44ftj;
        
        _g8vuazx[0]._gv44ftj = 0;                 _g8vuazx[0]._2c9t199 = 0.0f;
        _06ftlo4         = 0u;                _c7fov2y        = 0.0f;
        for (i = 0u; i < _ov6tt9q[_fsqg8x8]; i++)
        {   if (_r59c81d[_lx2toow[i]*8 +3] == _nezw9he)
            {   for (j = 0u; j < 3u; j++)
                {   _id4n1u3(_n42x2ez, &_mszdxim[_k00ilu6*4 +0], &_mszdxim[ _r59c81d[_lx2toow[i]*8 +j]*4 +0]);
                    _nbzug0i   = _mabsi9n(_n42x2ez);
                    _id4n1u3(_n42x2ez, &_mszdxim[_vfua9bq*4 +0], &_mszdxim[ _r59c81d[_lx2toow[i]*8 +j]*4 +0]);
                    _nbzug0i  += _mabsi9n(_n42x2ez);
                    _06ftlo4  = (_nbzug0i < _c7fov2y)*_06ftlo4   +  (_nbzug0i >= _c7fov2y)*_r59c81d[_lx2toow[i]*8 +j];
                    _c7fov2y = (_nbzug0i < _c7fov2y)*_c7fov2y  +  (_nbzug0i >= _c7fov2y)*_nbzug0i;
        }   }   }
        
        _g8vuazx[0]._gv44ftj = _06ftlo4;      _g8vuazx[0]._2c9t199 = _c7fov2y;
        
        MPI_Allreduce(_g8vuazx, _yllwb7s, 1, MPI_FLOAT_INT, MPI_MAXLOC, MPI_COMM_WORLD);
        
        _1btr23q = (unsigned int)_yllwb7s[0]._gv44ftj;
        
        _g8vuazx[0]._gv44ftj = 0;                 _g8vuazx[0]._2c9t199 = 0.0f;
        _06ftlo4         = 0u;                _c7fov2y        = 0.0f;
        for (i = 0u; i < _ov6tt9q[_fsqg8x8]; i++)
        {   if (_r59c81d[_lx2toow[i]*8 +3] == _nezw9he)
            {   for (j = 0u; j < 3u; j++)
                {   _id4n1u3(_n42x2ez, &_mszdxim[_k00ilu6*4 +0], &_mszdxim[_r59c81d[ _lx2toow[i]*8 +j]*4 +0]);
                    _nbzug0i   = _mabsi9n(_n42x2ez);
                    _id4n1u3(_n42x2ez, &_mszdxim[_vfua9bq*4 +0], &_mszdxim[_r59c81d[ _lx2toow[i]*8 +j]*4 +0]);
                    _nbzug0i  += _mabsi9n(_n42x2ez);
                    _id4n1u3(_n42x2ez, &_mszdxim[_1btr23q*4 +0], &_mszdxim[_r59c81d[ _lx2toow[i]*8 +j]*4 +0]);
                    _nbzug0i  += _mabsi9n(_n42x2ez);
                    _06ftlo4  = (_nbzug0i < _c7fov2y)*_06ftlo4   +  (_nbzug0i >= _c7fov2y)*_r59c81d[_lx2toow[i]*8 +j];
                    _c7fov2y = (_nbzug0i < _c7fov2y)*_c7fov2y  +  (_nbzug0i >= _c7fov2y)*_nbzug0i;
        }   }   }
        
        _g8vuazx[0]._gv44ftj = _06ftlo4;      _g8vuazx[0]._2c9t199 = _c7fov2y;
        
        MPI_Allreduce(_g8vuazx, _yllwb7s, 1, MPI_FLOAT_INT, MPI_MAXLOC, MPI_COMM_WORLD);
       
        _iua08ek = (unsigned int)_yllwb7s[0]._gv44ftj;
        







        
        memcpy(_dan085r, &_mszdxim[_k00ilu6*4 +0], 3*sizeof(float));
        memcpy(_ntdr5xy, &_mszdxim[_vfua9bq*4 +0], 3*sizeof(float));
        memcpy(_fc4kl8c, &_mszdxim[_1btr23q*4 +0], 3*sizeof(float));
        memcpy(_gn0t6dv, &_mszdxim[_iua08ek*4 +0], 3*sizeof(float));
        _id4n1u3(_xrum2lh, _ntdr5xy, _dan085r);              _bfuk9e1(_xrum2lh);
        _id4n1u3(_gcsotdh, _fc4kl8c, _dan085r);              _bfuk9e1(_gcsotdh);
        _j2qp89u(_088o98b, _xrum2lh, _gcsotdh);             _bfuk9e1(_088o98b);
        if (_088o98b[2] < 0.0f)  {           _088o98b[0] *= -1.0f;       _088o98b[1] *= -1.0f;       _088o98b[2] *= -1.0f;}
        
        if ((_k00ilu6 == _iua08ek) || (_vfua9bq == _iua08ek) || (_1btr23q == _iua08ek))
        {   

            _bjfvyb6[_nezw9he*8 +3] = _088o98b[0];
            _bjfvyb6[_nezw9he*8 +4] = _088o98b[1];
            _bjfvyb6[_nezw9he*8 +5] = _088o98b[2];
        }
        else
        {   

            _id4n1u3(_xrum2lh, _dan085r, _ntdr5xy);              _bfuk9e1(_xrum2lh);
            _id4n1u3(_gcsotdh, _gn0t6dv, _ntdr5xy);              _bfuk9e1(_gcsotdh);
            _j2qp89u(_hlmh1t9, _xrum2lh, _gcsotdh);             _bfuk9e1(_hlmh1t9);
            _nbzug0i = _088o98b[0]*_hlmh1t9[0] + _088o98b[1]*_hlmh1t9[1] + _088o98b[2]*_hlmh1t9[2];
            if (_nbzug0i < 0.0f)  {           _hlmh1t9[0] *= -1.0f;       _hlmh1t9[1] *= -1.0f;       _hlmh1t9[2] *= -1.0f;}
            
            _id4n1u3(_xrum2lh, _dan085r, _fc4kl8c);              _bfuk9e1(_xrum2lh);
            _id4n1u3(_gcsotdh, _gn0t6dv, _fc4kl8c);              _bfuk9e1(_gcsotdh);
            _j2qp89u(_39ayq6w, _xrum2lh, _gcsotdh);             _bfuk9e1(_39ayq6w);
            _nbzug0i = _088o98b[0]*_39ayq6w[0] + _088o98b[1]*_39ayq6w[1] + _088o98b[2]*_39ayq6w[2];
            if (_nbzug0i < 0.0f)  {           _39ayq6w[0] *= -1.0f;       _39ayq6w[1] *= -1.0f;       _39ayq6w[2] *= -1.0f;}
            
            _id4n1u3(_xrum2lh, _ntdr5xy, _gn0t6dv);              _bfuk9e1(_xrum2lh);
            _id4n1u3(_gcsotdh, _fc4kl8c, _gn0t6dv);              _bfuk9e1(_gcsotdh);
            _j2qp89u(_w8yl68n, _xrum2lh, _gcsotdh);             _bfuk9e1(_w8yl68n);
            _nbzug0i = _088o98b[0]*_w8yl68n[0] + _088o98b[1]*_w8yl68n[1] + _088o98b[2]*_w8yl68n[2];
            if (_nbzug0i < 0.0f)  {           _w8yl68n[0] *= -1.0f;       _w8yl68n[1] *= -1.0f;       _w8yl68n[2] *= -1.0f;}
            
            _bjfvyb6[_nezw9he*8 +3] = _088o98b[0] + _hlmh1t9[0] + _39ayq6w[0] + _w8yl68n[0];
            _bjfvyb6[_nezw9he*8 +4] = _088o98b[1] + _hlmh1t9[1] + _39ayq6w[1] + _w8yl68n[1];
            _bjfvyb6[_nezw9he*8 +5] = _088o98b[2] + _hlmh1t9[2] + _39ayq6w[2] + _w8yl68n[2];
        }
        _bfuk9e1(&_bjfvyb6[_nezw9he*8 +3]);
        
        _nbzug0i = 0.0f;
        for (i = 0u; i < _6zjq7w7; i++)
        {   if (_r59c81d[i*8 +3] == _nezw9he)
            {   _nbzug0i               += 1.0f;
                _bjfvyb6[_nezw9he*8 +0] += _yn4u9bx[i*16 +0];
                _bjfvyb6[_nezw9he*8 +1] += _yn4u9bx[i*16 +1];
                _bjfvyb6[_nezw9he*8 +2] += _yn4u9bx[i*16 +2];
                
                memcpy(_dan085r, &_mszdxim[_r59c81d[i*8 +0]*4 +0], 3*sizeof(float));
                memcpy(_ntdr5xy, &_mszdxim[_r59c81d[i*8 +1]*4 +0], 3*sizeof(float));
                memcpy(_fc4kl8c, &_mszdxim[_r59c81d[i*8 +2]*4 +0], 3*sizeof(float));
                _id4n1u3(_xrum2lh, _ntdr5xy, _dan085r);              _bfuk9e1(_xrum2lh);
                _id4n1u3(_gcsotdh, _fc4kl8c, _dan085r);              _bfuk9e1(_gcsotdh);
                _j2qp89u(_088o98b, _xrum2lh, _gcsotdh);             _bfuk9e1(_088o98b);
                _t26ys44       = _bjfvyb6[_nezw9he*8 +3]*_088o98b[0] + _bjfvyb6[_nezw9he*8 +4]*_088o98b[1] + _bjfvyb6[_nezw9he*8 +5]*_088o98b[2];
                
                _omxlvz4       = (_t26ys44 < 0.0f)*1u + (_t26ys44 >= 0.0f)*0u; 
                _qke2gjr       = _r59c81d[i*8 +1];
                _r59c81d[i*8 +1] = (_omxlvz4 == 0u)*_r59c81d[i*8 +1] + (_omxlvz4 != 0u)*_r59c81d[i*8 +2];
                _r59c81d[i*8 +2] = (_omxlvz4 == 0u)*_r59c81d[i*8 +2] + (_omxlvz4 != 0u)*_qke2gjr;
                
                _088o98b[0]      = (_omxlvz4 == 0u)*_088o98b[0]      + (_omxlvz4 != 0u)*-_088o98b[0];
                _088o98b[1]      = (_omxlvz4 == 0u)*_088o98b[1]      + (_omxlvz4 != 0u)*-_088o98b[1];
                _088o98b[2]      = (_omxlvz4 == 0u)*_088o98b[2]      + (_omxlvz4 != 0u)*-_088o98b[2];
                _ktakpux(_088o98b, _j3swebn, _795e5iz);
                
                memcpy(&_yn4u9bx[i*16 +6],  _088o98b, 3*sizeof(float));
                memcpy(&_yn4u9bx[i*16 +9],  _j3swebn, 3*sizeof(float));
                memcpy(&_yn4u9bx[i*16 +12], _795e5iz, 3*sizeof(float));
        }   }
        _bjfvyb6[_nezw9he*8 +0] /= _nbzug0i;
        _bjfvyb6[_nezw9he*8 +1] /= _nbzug0i;
        _bjfvyb6[_nezw9he*8 +2] /= _nbzug0i;
    } 
    
    return;
}


