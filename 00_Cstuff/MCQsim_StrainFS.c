# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <float.h>





void       _tb94874(double _ikwe4bx[3], double _rwedbcy, double _mdmmmxi, double _dmc2rh2, double _py6b2kn[3][3]); 

void      _w5wfcpk(int _26q3xlr[1], double x,double y,double z,double _b2rvndf,double _dw79wcx, double _2pgyme6,double _rpkzggm,double _kjnawng,double _a41izxt);

void        _37d4ivs(double _9iy4x41[6], double _s3jonou[6], double B[3][3]);

void         _3jgzlma(double x,double y,double z,double _52xg8p6,double _1750ceh,double _lxslg6m,double _iui9kkh,double _qa78rap, double _c54wsqt[3],double _6kof8kr[3],double e[6]);

void     _t3w5uzb(double x, double y, double z, double _52xg8p6, double _1750ceh, double _lxslg6m, double _iui9kkh, double _qa78rap, double e[6]);






void _25sq8cr(float _cjwtl9t[6], float _5lf6a4f[6], float _kshy8nb, float _e1hm38k, float _szbffzz, float _dan085r[3], float _ntdr5xy[3], float _fc4kl8c[3], float _h30153s, float _n8la100, float _e7eobjl, const float _emruyp0, const float _z7whwkm)



























































































{
   
    double X  = (double)_kshy8nb;         double Y = (double)_e1hm38k;          double Z = (double)_szbffzz; 
    double _1750ceh = (double)_e7eobjl;        double _lxslg6m= (double)_h30153s;         double _iui9kkh= (double)_n8la100;  
    double _mu7jae5 = (double)_emruyp0;        double _odck828= (double)_z7whwkm; 
    double _vkizzqy[3],                   _x63cfl0[3],                          _rcbpljd[3];
    _vkizzqy[0] = (double)_dan085r[0];         _vkizzqy[1] = (double)_dan085r[1];         _vkizzqy[2] = (double)_dan085r[2];
    _x63cfl0[0] = (double)_ntdr5xy[0];         _x63cfl0[1] = (double)_ntdr5xy[1];         _x63cfl0[2] = (double)_ntdr5xy[2];
    _rcbpljd[0] = (double)_fc4kl8c[0];         _rcbpljd[1] = (double)_fc4kl8c[1];         _rcbpljd[2] = (double)_fc4kl8c[2];
    
    int      i, _ds3qsqc,          _uflsehu;
    double   A,      B,      C,  x,  y,      z,      _qa78rap,  _zpeal2w;
    
    int      _r11q7uw[1];
    double   _ja2m0w1[3],        _j6fp8hq[3],            _s9toqw9[3];
    double   _0s1f537[3],      _ye0qzsa[3];
    double   _wlgje2a[3],             _z2ei921[3],              _0j36tcf[3];
    double   _txcb48y[3],            _1o4izas[3],             _nwywrra[3];
    double   _tatj9a3[3],             _y36ddfg[6];
    double   _dntcagr[6],        _p2u4dtb[6],         _h55rbc7[6];
    double   _zmf8mdd[3][3];

    _qa78rap = _odck828/(2.0*(_odck828+_mu7jae5)); 
    _wlgje2a[0]        = 0.0;                     _wlgje2a[1]        = 0.0;                     _wlgje2a[2]        = 0.0;
    _z2ei921[0]        = 0.0;                     _z2ei921[1]        = 0.0;                     _z2ei921[2]        = 0.0;
    _0j36tcf[0]        = 0.0;                     _0j36tcf[1]        = 0.0;                     _0j36tcf[2]        = 0.0;
    _tatj9a3[0]        = 0.0;                     _tatj9a3[1]        = 0.0;                     _tatj9a3[2]        = 1.0;
    
    

    _0s1f537[0] = _x63cfl0[0] -_vkizzqy[0];            _0s1f537[1] = _x63cfl0[1] -_vkizzqy[1];            _0s1f537[2] = _x63cfl0[2] -_vkizzqy[2];
    _ye0qzsa[0] = _rcbpljd[0] -_vkizzqy[0];            _ye0qzsa[1] = _rcbpljd[1] -_vkizzqy[1];            _ye0qzsa[2] = _rcbpljd[2] -_vkizzqy[2];
      
    _s9toqw9[0]     = _0s1f537[1]*_ye0qzsa[2] - _0s1f537[2]*_ye0qzsa[1];
    _s9toqw9[1]     = _0s1f537[2]*_ye0qzsa[0] - _0s1f537[0]*_ye0qzsa[2];
    _s9toqw9[2]     = _0s1f537[0]*_ye0qzsa[1] - _0s1f537[1]*_ye0qzsa[0];
    _zpeal2w   = sqrt(  _s9toqw9[0]*_s9toqw9[0] +_s9toqw9[1]*_s9toqw9[1] + _s9toqw9[2]*_s9toqw9[2]);
    
    _s9toqw9[0]     = _s9toqw9[0]/_zpeal2w;     _s9toqw9[1]     = _s9toqw9[1]/_zpeal2w;     _s9toqw9[2]     = _s9toqw9[2]/_zpeal2w;
   
    _ja2m0w1[0]   = _tatj9a3[1]*_s9toqw9[2] - _tatj9a3[2]*_s9toqw9[1];
    _ja2m0w1[1]   = _tatj9a3[2]*_s9toqw9[0] - _tatj9a3[0]*_s9toqw9[2];
    _ja2m0w1[2]   = _tatj9a3[0]*_s9toqw9[1] - _tatj9a3[1]*_s9toqw9[0];
    
    _zpeal2w    = sqrt(  _ja2m0w1[0]*_ja2m0w1[0] +_ja2m0w1[1]*_ja2m0w1[1] + _ja2m0w1[2]*_ja2m0w1[2]);
    if (_zpeal2w < FLT_EPSILON)
    {   _ja2m0w1[0] = 0.0 ;                  _ja2m0w1[1] = 1.0;            _ja2m0w1[2] = 0.0;    
    }
    else
    {   _ja2m0w1[0]/= _zpeal2w;            _ja2m0w1[1]/= _zpeal2w;     _ja2m0w1[2]/= _zpeal2w;
    }
    _j6fp8hq[0]    = _s9toqw9[1]*_ja2m0w1[2] - _s9toqw9[2]*_ja2m0w1[1];
    _j6fp8hq[1]    = _s9toqw9[2]*_ja2m0w1[0] - _s9toqw9[0]*_ja2m0w1[2];
    _j6fp8hq[2]    = _s9toqw9[0]*_ja2m0w1[1] - _s9toqw9[1]*_ja2m0w1[0];
    _zpeal2w = sqrt(  _j6fp8hq[0]*_j6fp8hq[0] +_j6fp8hq[1]*_j6fp8hq[1] + _j6fp8hq[2]*_j6fp8hq[2]);

    _j6fp8hq[0]    = _j6fp8hq[0]/_zpeal2w;     _j6fp8hq[1] = _j6fp8hq[1]/_zpeal2w;           _j6fp8hq[2]      = _j6fp8hq[2]/_zpeal2w;
    
    _zmf8mdd[0][0] = _s9toqw9[0];               _zmf8mdd[0][1] = _s9toqw9[1];                  _zmf8mdd[0][2]   = _s9toqw9[2];
    _zmf8mdd[1][0] = _ja2m0w1[0];             _zmf8mdd[1][1] = _ja2m0w1[1];                _zmf8mdd[1][2]   = _ja2m0w1[2];
    _zmf8mdd[2][0] = _j6fp8hq[0];                _zmf8mdd[2][1] = _j6fp8hq[1];                   _zmf8mdd[2][2]   = _j6fp8hq[2];

    _tb94874(_0s1f537, (X-_x63cfl0[0]),  (Y-_x63cfl0[1]), (Z-_x63cfl0[2]), _zmf8mdd); 
    x          = _0s1f537[0];                  y          = _0s1f537[1];                  z          = _0s1f537[2];
    _tb94874(_0s1f537, (_vkizzqy[0]-_x63cfl0[0]),  (_vkizzqy[1]-_x63cfl0[1]), (_vkizzqy[2]-_x63cfl0[2]), _zmf8mdd); 
    _wlgje2a[0]      = _0s1f537[0];                   _wlgje2a[1]      = _0s1f537[1];                  _wlgje2a[2]      = _0s1f537[2];
    _tb94874(_ye0qzsa, (_rcbpljd[0]-_x63cfl0[0]),  (_rcbpljd[1]-_x63cfl0[1]), (_rcbpljd[2]-_x63cfl0[2]), _zmf8mdd); 
    _0j36tcf[0]      = _ye0qzsa[0];                   _0j36tcf[1]      = _ye0qzsa[1];                  _0j36tcf[2]      = _ye0qzsa[2];
    
    _zpeal2w = sqrt(  (_z2ei921[0]-_wlgje2a[0])*(_z2ei921[0]-_wlgje2a[0]) +(_z2ei921[1]-_wlgje2a[1])*(_z2ei921[1]-_wlgje2a[1]) + (_z2ei921[2]-_wlgje2a[2])*(_z2ei921[2]-_wlgje2a[2]));
    
    _txcb48y[0]     = (_z2ei921[0]-_wlgje2a[0])/_zpeal2w;       _txcb48y[1]      = (_z2ei921[1]-_wlgje2a[1])/_zpeal2w;       _txcb48y[2]     = (_z2ei921[2]-_wlgje2a[2])/_zpeal2w;     
    _zpeal2w = sqrt(  (_0j36tcf[0]-_wlgje2a[0])*(_0j36tcf[0]-_wlgje2a[0]) +(_0j36tcf[1]-_wlgje2a[1])*(_0j36tcf[1]-_wlgje2a[1]) + (_0j36tcf[2]-_wlgje2a[2])*(_0j36tcf[2]-_wlgje2a[2]));
    _1o4izas[0]     = (_0j36tcf[0]-_wlgje2a[0])/_zpeal2w;       _1o4izas[1]      = (_0j36tcf[1]-_wlgje2a[1])/_zpeal2w;       _1o4izas[2]     = (_0j36tcf[2]-_wlgje2a[2])/_zpeal2w;
    _zpeal2w = sqrt(  (_0j36tcf[0]-_z2ei921[0])*(_0j36tcf[0]-_z2ei921[0]) +(_0j36tcf[1]-_z2ei921[1])*(_0j36tcf[1]-_z2ei921[1]) + (_0j36tcf[2]-_z2ei921[2])*(_0j36tcf[2]-_z2ei921[2]));
    _nwywrra[0]     = (_0j36tcf[0]-_z2ei921[0])/_zpeal2w;       _nwywrra[1]      = (_0j36tcf[1]-_z2ei921[1])/_zpeal2w;       _nwywrra[2]     = (_0j36tcf[2]-_z2ei921[2])/_zpeal2w; 
    
    _zpeal2w  = _txcb48y[0]*_1o4izas[0] +_txcb48y[1]*_1o4izas[1] +_txcb48y[2]*_1o4izas[2];
    A = acos(_zpeal2w) ;
    _zpeal2w  = -1.0*_txcb48y[0]*_nwywrra[0] + -1.0*_txcb48y[1]*_nwywrra[1] + -1.0*_txcb48y[2]*_nwywrra[2];
    B = acos(_zpeal2w) ;
    _zpeal2w  = _nwywrra[0]*_1o4izas[0] +_nwywrra[1]*_1o4izas[1] +_nwywrra[2]*_1o4izas[2];
    C = acos(_zpeal2w) ;
     
    

    _w5wfcpk(_r11q7uw,y,z,x,_wlgje2a[1],_wlgje2a[2], _z2ei921[1], _z2ei921[2], _0j36tcf[1], _0j36tcf[2]);

    if (_r11q7uw[0] == 1)       {       _ds3qsqc = 1;   _uflsehu = 0;          } 
    if (_r11q7uw[0] ==-1)       {       _ds3qsqc = 0;   _uflsehu = 1;          } 
    if (_r11q7uw[0] == 0)       {       _ds3qsqc = 0;   _uflsehu = 0;          } 

    if (_ds3qsqc == 1) 
    {   
        _0s1f537[0] = -1.0*_1o4izas[0];         _0s1f537[1] = -1.0*_1o4izas[1];         _0s1f537[2] = -1.0*_1o4izas[2];              
        _3jgzlma(x, y, z, A, _1750ceh, _lxslg6m, _iui9kkh, _qa78rap, _wlgje2a, _0s1f537, _dntcagr);
        
        _3jgzlma(x, y, z, B, _1750ceh, _lxslg6m, _iui9kkh, _qa78rap, _z2ei921, _txcb48y,       _p2u4dtb);
        
        _3jgzlma(x, y, z, C, _1750ceh, _lxslg6m, _iui9kkh, _qa78rap, _0j36tcf, _nwywrra,       _h55rbc7);  
    }
    if (_uflsehu == 1) 
    {  
        _3jgzlma(x, y, z, A, _1750ceh, _lxslg6m, _iui9kkh, _qa78rap, _wlgje2a, _1o4izas,       _dntcagr);
        
        _0s1f537[0] = -1.0*_txcb48y[0];         _0s1f537[1] = -1.0*_txcb48y[1];         _0s1f537[2] = -1.0*_txcb48y[2];      
        _3jgzlma(x, y, z, B, _1750ceh, _lxslg6m, _iui9kkh, _qa78rap, _z2ei921, _0s1f537, _p2u4dtb);
        
        _0s1f537[0] = -1.0*_nwywrra[0];         _0s1f537[1] = -1.0*_nwywrra[1];         _0s1f537[2] = -1.0*_nwywrra[2]; 
        _3jgzlma(x, y, z, C, _1750ceh, _lxslg6m, _iui9kkh, _qa78rap, _0j36tcf, _0s1f537, _h55rbc7);     
    }
    if ((_uflsehu == 1) || (_ds3qsqc == 1))
    {
        _y36ddfg[0]       = _dntcagr[0]+_p2u4dtb[0]+_h55rbc7[0]; 
        _y36ddfg[1]       = _dntcagr[1]+_p2u4dtb[1]+_h55rbc7[1]; 
        _y36ddfg[2]       = _dntcagr[2]+_p2u4dtb[2]+_h55rbc7[2]; 
        _y36ddfg[3]       = _dntcagr[3]+_p2u4dtb[3]+_h55rbc7[3]; 
        _y36ddfg[4]       = _dntcagr[4]+_p2u4dtb[4]+_h55rbc7[4]; 
        _y36ddfg[5]       = _dntcagr[5]+_p2u4dtb[5]+_h55rbc7[5]; 
    }
    else
    {   
        _y36ddfg[0]       = NAN; 
        _y36ddfg[1]       = NAN; 
        _y36ddfg[2]       = NAN; 
        _y36ddfg[3]       = NAN; 
        _y36ddfg[4]       = NAN; 
        _y36ddfg[5]       = NAN; 
    }

    _zmf8mdd[0][0] = _s9toqw9[0];         _zmf8mdd[0][1] = _ja2m0w1[0];       _zmf8mdd[0][2] = _j6fp8hq[0];
    _zmf8mdd[1][0] = _s9toqw9[1];         _zmf8mdd[1][1] = _ja2m0w1[1];       _zmf8mdd[1][2] = _j6fp8hq[1];
    _zmf8mdd[2][0] = _s9toqw9[2];         _zmf8mdd[2][1] = _ja2m0w1[2];       _zmf8mdd[2][2] = _j6fp8hq[2];
    

    double _78e9riv[6];
    _37d4ivs(_78e9riv, _y36ddfg, _zmf8mdd);
    for (i = 0; i < 6; i++) {   _5lf6a4f[i] = (float)_78e9riv[i];     }

    
    _cjwtl9t[0] = 2.0*_mu7jae5* _5lf6a4f[0]+_odck828*(_5lf6a4f[0]+_5lf6a4f[3]+_5lf6a4f[5]);       
    _cjwtl9t[3] = 2.0*_mu7jae5* _5lf6a4f[3]+_odck828*(_5lf6a4f[0]+_5lf6a4f[3]+_5lf6a4f[5]);       
    _cjwtl9t[5] = 2.0*_mu7jae5* _5lf6a4f[5]+_odck828*(_5lf6a4f[0]+_5lf6a4f[3]+_5lf6a4f[5]);       
    _cjwtl9t[1] = 2.0*_mu7jae5* _5lf6a4f[1];                                           
    _cjwtl9t[2] = 2.0*_mu7jae5* _5lf6a4f[2];                                           
    _cjwtl9t[4] = 2.0*_mu7jae5* _5lf6a4f[4];                                           
    
    return;
}






void  _tb94874(double _ikwe4bx[3], double _mv9e40k, double _zuj2xor, double _vork1dj, double A[3][3]) 
{
    
    
    
    
    _ikwe4bx[0] = A[0][0]*_mv9e40k + A[0][1]*_zuj2xor + A[0][2]*_vork1dj;     _ikwe4bx[1] = A[1][0]*_mv9e40k + A[1][1]*_zuj2xor + A[1][2]*_vork1dj;       _ikwe4bx[2] = A[2][0]*_mv9e40k + A[2][1]*_zuj2xor + A[2][2]*_vork1dj;

    return;
}



void _37d4ivs(double _9iy4x41[6], double _s3jonou[6], double B[3][3])
{
    
    
    
    double _gippsfn,         _dyk0x65,           _v7ln23y,           _ngh8rie;
    double _m1enf31,         _k0qtl3q,           _komv2xc,           _mufm8zo;
    double _30ysxpn,         _c4d145b,           _3fpoyip,           _4n89qkn;
    double A[9];
    
    _gippsfn = _s3jonou[0];         _dyk0x65 = _s3jonou[1];         _v7ln23y = _s3jonou[2];
    _ngh8rie = _s3jonou[3];         _m1enf31 = _s3jonou[4];         _k0qtl3q = _s3jonou[5];

    A[0] = B[0][0];         A[1] = B[1][0];         A[2] = B[2][0];
    A[3] = B[0][1];         A[4] = B[1][1];         A[5] = B[2][1];
    A[6] = B[0][2];         A[7] = B[1][2];         A[8] = B[2][2];

    _komv2xc = A[0]*A[0]*_gippsfn +         2.0*A[0]*A[3] *_dyk0x65 +          2.0*A[0]*A[6] *_v7ln23y +          2.0*A[3]*A[6] *_m1enf31 + A[3]*A[3]*_ngh8rie + A[6]*A[6]*_k0qtl3q;
    _c4d145b = A[1]*A[1]*_gippsfn +         2.0*A[1]*A[4] *_dyk0x65 +          2.0*A[1]*A[7] *_v7ln23y +          2.0*A[4]*A[7] *_m1enf31 + A[4]*A[4]*_ngh8rie + A[7]*A[7]*_k0qtl3q;
    _4n89qkn = A[2]*A[2]*_gippsfn +         2.0*A[2]*A[5] *_dyk0x65 +          2.0*A[2]*A[8] *_v7ln23y +          2.0*A[5]*A[8] *_m1enf31 + A[5]*A[5]*_ngh8rie + A[8]*A[8]*_k0qtl3q;
    _mufm8zo = A[0]*A[1]*_gippsfn + (A[0]*A[4]+ A[1]*A[3])*_dyk0x65 + (A[0]*A[7] + A[1]*A[6])*_v7ln23y + (A[7]*A[3] + A[6]*A[4])*_m1enf31 + A[4]*A[3]*_ngh8rie + A[6]*A[7]*_k0qtl3q;
    _30ysxpn = A[0]*A[2]*_gippsfn + (A[0]*A[5]+ A[2]*A[3])*_dyk0x65 + (A[0]*A[8] + A[2]*A[6])*_v7ln23y + (A[8]*A[3] + A[6]*A[5])*_m1enf31 + A[5]*A[3]*_ngh8rie + A[6]*A[8]*_k0qtl3q;
    _3fpoyip = A[1]*A[2]*_gippsfn + (A[2]*A[4]+ A[1]*A[5])*_dyk0x65 + (A[2]*A[7] + A[1]*A[8])*_v7ln23y + (A[7]*A[5] + A[8]*A[4])*_m1enf31 + A[4]*A[5]*_ngh8rie + A[7]*A[8]*_k0qtl3q;
    
    _9iy4x41[0] = _komv2xc;        _9iy4x41[1] = _mufm8zo;      _9iy4x41[2] = _30ysxpn;
    _9iy4x41[3] = _c4d145b;        _9iy4x41[4] = _3fpoyip;      _9iy4x41[5] = _4n89qkn;
    
    return;
}



void      _w5wfcpk(int _r11q7uw[1], double x,double y,double z,double _b2rvndf,double _dw79wcx, double _2pgyme6,double _rpkzggm,double _kjnawng,double _a41izxt)
{
    
    
    
    
    
    double a,            b,          c;

    a = ((_rpkzggm-_a41izxt)*(x-_kjnawng) +(_kjnawng-_2pgyme6)*(y-_a41izxt)) /  ((_rpkzggm-_a41izxt)*(_b2rvndf-_kjnawng) +(_kjnawng-_2pgyme6)*(_dw79wcx-_a41izxt));
    b = ((_a41izxt-_dw79wcx)*(x-_kjnawng) +(_b2rvndf-_kjnawng)*(y-_a41izxt)) /  ((_rpkzggm-_a41izxt)*(_b2rvndf-_kjnawng) +(_kjnawng-_2pgyme6)*(_dw79wcx-_a41izxt));
    c = 1.0 -a -b;

    _r11q7uw[0] = 1;
    if       ((a <= 0.0) && (b > c) && (c > a))             {   _r11q7uw[0] = -1;                  }
    else if  ((b <= 0.0) && (c > a) && (a > b))             {   _r11q7uw[0] = -1;                  }
    else if  ((c <= 0.0) && (a > b) && (b > c))             {   _r11q7uw[0] = -1;                  }

    else if  ((a == 0.0) && (b >= 0.0) && (c >= 0.0))       {   _r11q7uw[0] = 0;                   }
    else if  ((a >= 0.0) && (b == 0.0) && (c >= 0.0))       {   _r11q7uw[0] = 0;                   }
    else if  ((a >= 0.0) && (b >= 0.0) && (c == 0.0))       {   _r11q7uw[0] = 0;                   }
    if       ((_r11q7uw[0] == 0) && (z != 0.0))              {   _r11q7uw[0] = 1;                   } 

    return;
}



void _3jgzlma(double x,double y,double z,double _52xg8p6,double _1750ceh,double _lxslg6m,double _iui9kkh,double _qa78rap, double _c54wsqt[3],double _6kof8kr[3],double _9iy4x41[6])
{   
    
    double A[2][2];
    double B[3][3];
    double _vk0ucc1[2];
    double _m3o55ix[2];
    double _8nyfwdc;
    double _v9nou3d;
    double _m1zxyvg;
    double _zlm72je;
    double _0s1f537[2];
    double e[6];
    
    
    A[0][0] = _6kof8kr[2];                   A[0][1]      = -1.0*_6kof8kr[1];
    A[1][0] = _6kof8kr[1];                   A[1][1]      =      _6kof8kr[2]; 
    
    _0s1f537[0] = y -_c54wsqt[1];         _0s1f537[1] = z -_c54wsqt[2];
    _vk0ucc1[0]        = A[0][0]*_0s1f537[0] +A[0][1]*_0s1f537[1];
    _vk0ucc1[1]        = A[1][0]*_0s1f537[0] +A[1][1]*_0s1f537[1];
    _8nyfwdc           = _vk0ucc1[0];                   _v9nou3d           = _vk0ucc1[1];
    
    _0s1f537[0] = _lxslg6m;                      _0s1f537[1] = _iui9kkh;
    _m3o55ix[0]        = A[0][0]*_0s1f537[0] +A[0][1]*_0s1f537[1];
    _m3o55ix[1]        = A[1][0]*_0s1f537[0] +A[1][1]*_0s1f537[1];
    _m1zxyvg          = _m3o55ix[0];                   _zlm72je          = _m3o55ix[1];
    
    
    _t3w5uzb(x, _8nyfwdc, _v9nou3d, (-1.0*M_PI+_52xg8p6), _1750ceh, _m1zxyvg, _zlm72je, _qa78rap, e); 
    
    B[0][0] = 1.0;          B[0][1] = 0.0;      B[0][2] = 0.0;
    B[1][0] = 0.0;          B[1][1] = A[0][0];  B[1][2] = A[1][0];
    B[2][0] = 0.0;          B[2][1] = A[0][1];  B[2][2] = A[1][1];
    
    _37d4ivs(_9iy4x41, e, B); 




    return;
}



void _t3w5uzb(double x, double y, double z, double _52xg8p6, double _1750ceh, double _lxslg6m, double _iui9kkh, double _qa78rap, double e[6])

{   double       _vyv2gqu,           _3p15ahd,           _ldalmje,            _xk2ptho;
    double       _zuj2xor,             _m4t8bp6,             _pr424tk,             _m3o55ix;
    double       r,              _fshb0oq,             _25lzse9,             _va1nkbt;
    double       _vshkoec,            W,              _qzw3c1q,             _xie0rnu;
    double       _kbhb7s9,            _jjq9k6r,            _tj1gm25,           C;
    double       S,              _ubods1e,         _eurxhw0,         _ic1l7dy;
    double       _0tpac56,            _t4ufel6,            _246hid2,            _5owf759;
    double       _35a0l72,            _53j2d54;  
    

    _vyv2gqu = sin(_52xg8p6);                          _3p15ahd = cos(_52xg8p6);
    _ldalmje = y*_3p15ahd - z*_vyv2gqu;                      _xk2ptho = y*_vyv2gqu + z*_3p15ahd;

    _zuj2xor = x*x;               _m4t8bp6   = y*y;                 _pr424tk  = z*z;
    _m3o55ix = _zuj2xor + _m4t8bp6 + _pr424tk;      r    = pow(_m3o55ix,0.5);         _fshb0oq  = r*_m3o55ix;
    _25lzse9 = r*(r-z);           _va1nkbt = _m3o55ix*pow((r-z),2.0);   _vshkoec = _fshb0oq*(r-z);

    W  = _xk2ptho-r;            _qzw3c1q   = W*W;                 _xie0rnu  = W*r;
    _kbhb7s9= _qzw3c1q*r;              _jjq9k6r  = W*_fshb0oq;                _tj1gm25= _qzw3c1q*_m3o55ix;

    C = (r*_3p15ahd-z)/_xie0rnu;      S    = (r*_vyv2gqu-y)/_xie0rnu;
    
    _ubods1e = (_ldalmje/r/(r-_xk2ptho)  -y/r/(r-z))   /4.0/M_PI;
    _eurxhw0 = (  x/r/(r-z)-_3p15ahd*x/r/(r-_xk2ptho))/4.0/M_PI;
    _ic1l7dy = (_vyv2gqu*x/r/(r-_xk2ptho))/4.0/M_PI;
    
    _0tpac56 = _1750ceh*(_ubods1e)  +_1750ceh/8.0/M_PI/(1.0-_qa78rap)*(_ldalmje/_xie0rnu+_ldalmje*_zuj2xor/_tj1gm25-_ldalmje*_zuj2xor/_jjq9k6r+y/_25lzse9-_zuj2xor*y/_va1nkbt-_zuj2xor*y/_vshkoec)
                       -_lxslg6m*x/8.0/M_PI/(1.0-_qa78rap)*(((2.0*_qa78rap+1.0)/_xie0rnu+_zuj2xor/_tj1gm25-_zuj2xor/_jjq9k6r)*_3p15ahd+(2.0*_qa78rap+1.0)/_25lzse9-_zuj2xor/_va1nkbt-_zuj2xor/_vshkoec)
                       +_iui9kkh*x*_vyv2gqu/8.0/M_PI/(1.0-_qa78rap)*((2.0*_qa78rap+1.0)/_xie0rnu+_zuj2xor/_tj1gm25-_zuj2xor/_jjq9k6r);

    _5owf759 = _lxslg6m*(_eurxhw0)  +_1750ceh/8.0/M_PI/(1.0-_qa78rap)*((1.0/_xie0rnu+S*S-_m4t8bp6/_jjq9k6r)*_ldalmje+(2.0*_qa78rap+1.0)*y/_25lzse9-pow(y,3.0)/_va1nkbt-pow(y,3.0)/_vshkoec-2.0*_qa78rap*_3p15ahd*S)
                       -_lxslg6m*x/8.0/M_PI/(1.0-_qa78rap)*(1.0/_25lzse9-_m4t8bp6/_va1nkbt-_m4t8bp6/_vshkoec+(1.0/_xie0rnu+S*S-_m4t8bp6/_jjq9k6r)*_3p15ahd)
                       +_iui9kkh*x*_vyv2gqu/8.0/M_PI/(1.0-_qa78rap)*(1.0/_xie0rnu+S*S-_m4t8bp6/_jjq9k6r);

    _53j2d54 = _iui9kkh*(_ic1l7dy)+_1750ceh/8.0/M_PI/(1.0-_qa78rap)*(_ldalmje/W/r+_ldalmje*C*C-_ldalmje*_pr424tk/_jjq9k6r+y*z/_fshb0oq+2.0*_qa78rap*_vyv2gqu*C)
                      -_lxslg6m*x/8.0/M_PI/(1.0-_qa78rap)*((1.0/_xie0rnu+C*C-_pr424tk/_jjq9k6r)*_3p15ahd+z/_fshb0oq)
                      +_iui9kkh*x*_vyv2gqu/8.0/M_PI/(1.0-_qa78rap)*(1.0/_xie0rnu+C*C-_pr424tk/_jjq9k6r);

    _t4ufel6 = _1750ceh*(_eurxhw0)/2.0+_lxslg6m*(_ubods1e)/2.0  -_1750ceh/8.0/M_PI/(1.0-_qa78rap)*(x*_m4t8bp6/_va1nkbt-_qa78rap*x/_25lzse9+x*_m4t8bp6/_vshkoec-_qa78rap*x*_3p15ahd/_xie0rnu+_ldalmje*x*S/_xie0rnu+_ldalmje*x*y/_jjq9k6r)
                                           +_lxslg6m/8.0/M_PI/(1.0-_qa78rap)*(_zuj2xor*y/_va1nkbt-_qa78rap*y/_25lzse9+_zuj2xor*y/_vshkoec+_qa78rap*_3p15ahd*S+_zuj2xor*y*_3p15ahd/_jjq9k6r+_zuj2xor*_3p15ahd*S/_xie0rnu)
                                           -_iui9kkh*_vyv2gqu/8.0/M_PI/(1.0-_qa78rap)*(_qa78rap*S+_zuj2xor*S/_xie0rnu+_zuj2xor*y/_jjq9k6r);

    _246hid2 = _1750ceh*(_ic1l7dy)/2.0+_iui9kkh*(_ubods1e)/2.0  -_1750ceh/8.0/M_PI/(1.0-_qa78rap)*(-x*y/_fshb0oq+_qa78rap*x*_vyv2gqu/_xie0rnu+_ldalmje*x*C/_xie0rnu+_ldalmje*x*z/_jjq9k6r)
                                           +_lxslg6m/8.0/M_PI/(1.0-_qa78rap)*(-_zuj2xor/_fshb0oq+_qa78rap/r+_qa78rap*_3p15ahd*C+_zuj2xor*z*_3p15ahd/_jjq9k6r+_zuj2xor*_3p15ahd*C/_xie0rnu)
                                           -_iui9kkh*_vyv2gqu/8.0/M_PI/(1.0-_qa78rap)*(_qa78rap*C+_zuj2xor*C/_xie0rnu+_zuj2xor*z/_jjq9k6r);

    _35a0l72 = _lxslg6m*(_ic1l7dy)/2.0+_iui9kkh*(_eurxhw0)/2.0  +_1750ceh/8.0/M_PI/(1.0-_qa78rap)*(_m4t8bp6/_fshb0oq-_qa78rap/r-_qa78rap*_3p15ahd*C+_qa78rap*_vyv2gqu*S+_ldalmje*_vyv2gqu*_3p15ahd/_qzw3c1q-_ldalmje*(y*_3p15ahd+z*_vyv2gqu)/_kbhb7s9+_ldalmje*y*z/_tj1gm25-_ldalmje*y*z/_jjq9k6r)
                                           -_lxslg6m*x/8.0/M_PI/(1.0-_qa78rap)*(y/_fshb0oq+_vyv2gqu*_3p15ahd*_3p15ahd/_qzw3c1q-_3p15ahd*(y*_3p15ahd+z*_vyv2gqu)/_kbhb7s9+y*z*_3p15ahd/_tj1gm25-y*z*_3p15ahd/_jjq9k6r)
                                           -_iui9kkh*x*_vyv2gqu/8.0/M_PI/(1.0-_qa78rap)*(y*z/_jjq9k6r-_vyv2gqu*_3p15ahd/_qzw3c1q+(y*_3p15ahd+z*_vyv2gqu)/_kbhb7s9-y*z/_tj1gm25);

    e[0] = _0tpac56;             e[1] = _t4ufel6;             e[2] = _246hid2;
    e[3] = _5owf759;             e[4] = _35a0l72;             e[5] = _53j2d54;

    return;
}




