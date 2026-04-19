# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <float.h>

void       _hpmnz1u(double _z8vmuhx[6], double _ieuc1hh[6],  double X, double Y, double Z, double _vkizzqy[3], double _x63cfl0[3], double _rcbpljd[3], double _p0fodro, double _zfk8ut0, double _fcd6mzs, double _mu7jae5, double _odck828);

void _dbfwyjx(double _3dd03uv[6],double _yt9moge[6], double X, double Y, double Z, double _vkizzqy[3] ,double _x63cfl0[3], double _rcbpljd[3], double _p0fodro, double _zfk8ut0, double _fcd6mzs, double _mu7jae5, double _odck828);

void       _frhocn9(double _ikwe4bx[3], double _rwedbcy, double _mdmmmxi, double _dmc2rh2, double _py6b2kn[3][3]); 

void      _jcslqc3(int _26q3xlr[1], double x,double y,double z,double _b2rvndf,double _dw79wcx, double _2pgyme6,double _rpkzggm,double _kjnawng,double _a41izxt);

void        _olmlgse(double _9iy4x41[6], double _s3jonou[6], double B[3][3]);

void         _5uofwb0(double x,double y,double z,double _52xg8p6,double _1750ceh,double _lxslg6m,double _iui9kkh,double _qa78rap, double _c54wsqt[3],double _6kof8kr[3],double e[6]);

void     _ikcidkw(double x, double y, double z, double _52xg8p6, double _1750ceh, double _lxslg6m, double _iui9kkh, double _qa78rap, double e[6]);

void    _nnhj7xv(double _qvdz4x1[6],double _ejgxaa2[6], double X,double Y,double Z,double _1nouhtc,double _s1k9ye9,double _o2m71hs,double _s3kkv9r[3], double _ftmttql[3], double _mu7jae5,double _odck828);

void  _4v30t9n(double _8nyfwdc, double _m4t8bp6, double _nlnqbu4, double _slwf90k, double _1rsrla8, double _s8hmbr5, double _xqt4p50, double _qa78rap, double a, double _5lf6a4f[6]);


void _oyr6ccv(float _cjwtl9t[6], float _5lf6a4f[6], float _kshy8nb, float _e1hm38k, float _szbffzz, float _dan085r[3], float _ntdr5xy[3], float _fc4kl8c[3], float _h30153s, float _n8la100, float _e7eobjl, const float _emruyp0, const float _z7whwkm)



























































































{

    double X  = (double)_kshy8nb;         double Y = (double)_e1hm38k;          double Z = (double)_szbffzz; 
    double _p0fodro = (double)_h30153s;        double _zfk8ut0= (double)_n8la100;         double _fcd6mzs= (double)_e7eobjl;
    double _mu7jae5 = (double)_emruyp0;        double _odck828= (double)_z7whwkm; 
    double _vkizzqy[3],                   _x63cfl0[3],                          _rcbpljd[3];
    _vkizzqy[0] = (double)_dan085r[0];         _vkizzqy[1] = (double)_dan085r[1];         _vkizzqy[2] = (double)_dan085r[2];
    _x63cfl0[0] = (double)_ntdr5xy[0];         _x63cfl0[1] = (double)_ntdr5xy[1];         _x63cfl0[2] = (double)_ntdr5xy[2];
    _rcbpljd[0] = (double)_fc4kl8c[0];         _rcbpljd[1] = (double)_fc4kl8c[1];         _rcbpljd[2] = (double)_fc4kl8c[2];

    int   i;
    double _z8vmuhx[6],                 _3dd03uv[6],                  _ulemx15[6];
    double _ieuc1hh[6],                 _yt9moge[6],                  _4qxxva2[6];
    double _72wx7pt[3],                _ukeyz7s[3],                  _84u89yj[3];
    
    _72wx7pt[0] = _vkizzqy[0];              _ukeyz7s[0] = _x63cfl0[0];          _84u89yj[0] = _rcbpljd[0];
    _72wx7pt[1] = _vkizzqy[1];              _ukeyz7s[1] = _x63cfl0[1];          _84u89yj[1] = _rcbpljd[1];
    _72wx7pt[2] = -1.0*_vkizzqy[2];         _ukeyz7s[2] = -1.0*_x63cfl0[2];     _84u89yj[2] = -1.0*_rcbpljd[2];
    
    for (i = 0; i < 6; i++)
    {   _z8vmuhx[i]   = 0.0;           _ieuc1hh[i]  = 0.0;     
        _3dd03uv[i]  = 0.0;           _yt9moge[i] = 0.0;
        _ulemx15[i]   = 0.0;           _4qxxva2[i]  = 0.0;
    }
    
    
    if ((Z > 0.0) || (_vkizzqy[2] > 0.0) || (_x63cfl0[2] > 0.0) || (_rcbpljd[2] > 0.0))
    {   fprintf(stdout,"A triangle vertex or the center of testing location have positive z-position (above half-space) => abort\n");           
        fprintf(stdout,"Z: %f    P1[2]: %f    P2[2]: %f      P3[2]: %f   \n",Z, _vkizzqy[2], _x63cfl0[2],_rcbpljd[2]);       
        exit(10);
    }
    if ((_vkizzqy[2] == 0.0) && (_x63cfl0[2] == 0.0) && (_rcbpljd[2] == 0.0))
    {   for (i = 0; i < 6; i++)
        {   _cjwtl9t[i] = 0.0;
            _5lf6a4f[i] = 0.0;
    }   }
    else
    {
        
        _hpmnz1u(       _z8vmuhx, _ieuc1hh,  X, Y, Z, _vkizzqy,     _x63cfl0,     _rcbpljd,     _p0fodro, _zfk8ut0, _fcd6mzs, _mu7jae5, _odck828);
        
        _dbfwyjx(_3dd03uv, _yt9moge, X, Y, Z, _vkizzqy,     _x63cfl0,     _rcbpljd,     _p0fodro, _zfk8ut0, _fcd6mzs, _mu7jae5, _odck828);
        
        _hpmnz1u(       _ulemx15, _4qxxva2,  X, Y, Z, _72wx7pt, _ukeyz7s, _84u89yj, _p0fodro, _zfk8ut0, _fcd6mzs, _mu7jae5, _odck828);

        if ((_72wx7pt[2] == 0.0) && (_ukeyz7s[2] == 0.0) && (_84u89yj[2] == 0.0))
        {   _ulemx15[4] = -1.0*_ulemx15[4];           _ulemx15[5] = -1.0*_ulemx15[5];
            _4qxxva2[4] = -1.0*_4qxxva2[4];           _4qxxva2[5] = -1.0*_4qxxva2[5];
        }
        
        for (i = 0; i < 6; i++)
        {   _cjwtl9t[i] = (_z8vmuhx[i] +_ulemx15[i] +_3dd03uv[i]);
            _5lf6a4f[i] = (_ieuc1hh[i] +_4qxxva2[i] +_yt9moge[i]);
    }   }
    return;
}

void _hpmnz1u(double _z8vmuhx[6], double _ieuc1hh[6],  double X, double Y, double Z, double _vkizzqy[3], double _x63cfl0[3], double _rcbpljd[3], double _lxslg6m, double _iui9kkh, double _1750ceh, double _mu7jae5, double _odck828)
{
    
    int      _ds3qsqc,          _uflsehu;
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

    _frhocn9(_0s1f537, (X-_x63cfl0[0]),  (Y-_x63cfl0[1]), (Z-_x63cfl0[2]), _zmf8mdd); 
    x          = _0s1f537[0];                  y          = _0s1f537[1];                  z          = _0s1f537[2];
    _frhocn9(_0s1f537, (_vkizzqy[0]-_x63cfl0[0]),  (_vkizzqy[1]-_x63cfl0[1]), (_vkizzqy[2]-_x63cfl0[2]), _zmf8mdd); 
    _wlgje2a[0]      = _0s1f537[0];                   _wlgje2a[1]      = _0s1f537[1];                  _wlgje2a[2]      = _0s1f537[2];
    _frhocn9(_ye0qzsa, (_rcbpljd[0]-_x63cfl0[0]),  (_rcbpljd[1]-_x63cfl0[1]), (_rcbpljd[2]-_x63cfl0[2]), _zmf8mdd); 
    _0j36tcf[0]      = _ye0qzsa[0];                   _0j36tcf[1]      = _ye0qzsa[1];                  _0j36tcf[2]      = _ye0qzsa[2];
    
    _zpeal2w = sqrt(  (_z2ei921[0]-_wlgje2a[0])*(_z2ei921[0]-_wlgje2a[0]) +(_z2ei921[1]-_wlgje2a[1])*(_z2ei921[1]-_wlgje2a[1]) + (_z2ei921[2]-_wlgje2a[2])*(_z2ei921[2]-_wlgje2a[2]));
    
    _txcb48y[0]     = (_z2ei921[0]-_wlgje2a[0])/_zpeal2w;       _txcb48y[1]      = (_z2ei921[1]-_wlgje2a[1])/_zpeal2w;       _txcb48y[2]     = (_z2ei921[2]-_wlgje2a[2])/_zpeal2w;     
    _zpeal2w = sqrt(  (_0j36tcf[0]-_wlgje2a[0])*(_0j36tcf[0]-_wlgje2a[0]) +(_0j36tcf[1]-_wlgje2a[1])*(_0j36tcf[1]-_wlgje2a[1]) + (_0j36tcf[2]-_wlgje2a[2])*(_0j36tcf[2]-_wlgje2a[2]));
    _1o4izas[0]     = (_0j36tcf[0]-_wlgje2a[0])/_zpeal2w;       _1o4izas[1]      = (_0j36tcf[1]-_wlgje2a[1])/_zpeal2w;       _1o4izas[2]     = (_0j36tcf[2]-_wlgje2a[2])/_zpeal2w;
    _zpeal2w = sqrt(  (_0j36tcf[0]-_z2ei921[0])*(_0j36tcf[0]-_z2ei921[0]) +(_0j36tcf[1]-_z2ei921[1])*(_0j36tcf[1]-_z2ei921[1]) + (_0j36tcf[2]-_z2ei921[2])*(_0j36tcf[2]-_z2ei921[2]));
    _nwywrra[0]     = (_0j36tcf[0]-_z2ei921[0])/_zpeal2w;       _nwywrra[1]      = (_0j36tcf[1]-_z2ei921[1])/_zpeal2w;       _nwywrra[2]     = (_0j36tcf[2]-_z2ei921[2])/_zpeal2w; 
    
    _zpeal2w  =      _txcb48y[0]*_1o4izas[0] +      _txcb48y[1]*_1o4izas[1] +      _txcb48y[2]*_1o4izas[2];
    A = acos(_zpeal2w) ;
    _zpeal2w  = -1.0*_txcb48y[0]*_nwywrra[0] + -1.0*_txcb48y[1]*_nwywrra[1] + -1.0*_txcb48y[2]*_nwywrra[2];
    B = acos(_zpeal2w) ;
    _zpeal2w  =      _nwywrra[0]*_1o4izas[0] +      _nwywrra[1]*_1o4izas[1] +     _nwywrra[2]*_1o4izas[2];
    C = acos(_zpeal2w) ;
     
    

    _jcslqc3(_r11q7uw,y,z,x,_wlgje2a[1],_wlgje2a[2], _z2ei921[1], _z2ei921[2], _0j36tcf[1], _0j36tcf[2]);

    if (_r11q7uw[0] == 1)       {       _ds3qsqc = 1;   _uflsehu = 0;      } 
    if (_r11q7uw[0] ==-1)       {       _ds3qsqc = 0;   _uflsehu = 1;      } 
    if (_r11q7uw[0] == 0)       {       _ds3qsqc = 0;   _uflsehu = 0;      } 

    if (_ds3qsqc == 1) 
    {   
        _0s1f537[0] = -1.0*_1o4izas[0];         _0s1f537[1] = -1.0*_1o4izas[1];         _0s1f537[2] = -1.0*_1o4izas[2];              
        _5uofwb0(x, y, z, A, _1750ceh, _lxslg6m, _iui9kkh, _qa78rap, _wlgje2a, _0s1f537, _dntcagr);
        
        _5uofwb0(x, y, z, B, _1750ceh, _lxslg6m, _iui9kkh, _qa78rap, _z2ei921, _txcb48y,       _p2u4dtb);
        
        _5uofwb0(x, y, z, C, _1750ceh, _lxslg6m, _iui9kkh, _qa78rap, _0j36tcf, _nwywrra,       _h55rbc7);  
    }
    if (_uflsehu == 1) 
    {                
        _5uofwb0(x, y, z, A, _1750ceh, _lxslg6m, _iui9kkh, _qa78rap, _wlgje2a, _1o4izas,       _dntcagr);
        
        _0s1f537[0] = -1.0*_txcb48y[0];         _0s1f537[1] = -1.0*_txcb48y[1];         _0s1f537[2] = -1.0*_txcb48y[2];      
        _5uofwb0(x, y, z, B, _1750ceh, _lxslg6m, _iui9kkh, _qa78rap, _z2ei921, _0s1f537, _p2u4dtb);
        
        _0s1f537[0] = -1.0*_nwywrra[0];         _0s1f537[1] = -1.0*_nwywrra[1];         _0s1f537[2] = -1.0*_nwywrra[2]; 
        _5uofwb0(x, y, z, C, _1750ceh, _lxslg6m, _iui9kkh, _qa78rap, _0j36tcf, _0s1f537, _h55rbc7);     
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
    
    _olmlgse(_ieuc1hh, _y36ddfg, _zmf8mdd);

    
    _z8vmuhx[0] = 2.0*_mu7jae5*_ieuc1hh[0]+_odck828*(_ieuc1hh[0]+_ieuc1hh[3]+_ieuc1hh[5]);  
    _z8vmuhx[3] = 2.0*_mu7jae5*_ieuc1hh[3]+_odck828*(_ieuc1hh[0]+_ieuc1hh[3]+_ieuc1hh[5]);  
    _z8vmuhx[5] = 2.0*_mu7jae5*_ieuc1hh[5]+_odck828*(_ieuc1hh[0]+_ieuc1hh[3]+_ieuc1hh[5]);  
    _z8vmuhx[1] = 2.0*_mu7jae5*_ieuc1hh[1];                            
    _z8vmuhx[2] = 2.0*_mu7jae5*_ieuc1hh[2];                            
    _z8vmuhx[4] = 2.0*_mu7jae5*_ieuc1hh[4];                            
    
    return;
}

void _dbfwyjx(double _3dd03uv[6], double _yt9moge[6], double X, double Y, double Z, double _vkizzqy[3], double _x63cfl0[3], double _rcbpljd[3], double _lxslg6m, double _iui9kkh, double _1750ceh, double _mu7jae5, double _odck828)
{
    
    

    
    
    int   i;          
    double _0s1f537[3],     _ye0qzsa[3];
    double _qvdz4x1[6],       _zhqimqo[6],         _dpxnf7c[6];
    double _ejgxaa2[6],       _ilgb0s0[6],         _r2sds8b[6];
    
    double _zmf8mdd[3][3];

    
    
    double _zpeal2w,                      _tatj9a3[3];
    double _s9toqw9[3],                        _ja2m0w1[3],                    _j6fp8hq[3];
    
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

    
    _zmf8mdd[0][0] = _s9toqw9[0];          _zmf8mdd[0][1] = _ja2m0w1[0];            _zmf8mdd[0][2] = _j6fp8hq[0];
    _zmf8mdd[1][0] = _s9toqw9[1];          _zmf8mdd[1][1] = _ja2m0w1[1];            _zmf8mdd[1][2] = _j6fp8hq[1];
    _zmf8mdd[2][0] = _s9toqw9[2];          _zmf8mdd[2][1] = _ja2m0w1[2];            _zmf8mdd[2][2] = _j6fp8hq[2];

    _frhocn9(_0s1f537, _1750ceh,  _lxslg6m, _iui9kkh, _zmf8mdd); 

    
    _nnhj7xv(_qvdz4x1,_ejgxaa2, X,Y,Z,_0s1f537[0],_0s1f537[1],_0s1f537[2],_vkizzqy,_x63cfl0,_mu7jae5,_odck828); 
    _nnhj7xv(_zhqimqo,_ilgb0s0, X,Y,Z,_0s1f537[0],_0s1f537[1],_0s1f537[2],_x63cfl0,_rcbpljd,_mu7jae5,_odck828); 
    _nnhj7xv(_dpxnf7c,_r2sds8b, X,Y,Z,_0s1f537[0],_0s1f537[1],_0s1f537[2],_rcbpljd,_vkizzqy,_mu7jae5,_odck828); 
    
    for (i = 0; i < 6; i++)
    {   _3dd03uv[i] = _qvdz4x1[i] + _zhqimqo[i] + _dpxnf7c[i];
        _yt9moge[i] = _ejgxaa2[i] + _ilgb0s0[i] + _r2sds8b[i];
    }
    
    return;
}

void _olmlgse(double _9iy4x41[6], double _s3jonou[6], double B[3][3])
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

void  _frhocn9(double _ikwe4bx[3], double _mv9e40k, double _zuj2xor, double _vork1dj, double A[3][3]) 
{
    
     
    
    
    _ikwe4bx[0] = A[0][0]*_mv9e40k + A[0][1]*_zuj2xor + A[0][2]*_vork1dj;     _ikwe4bx[1] = A[1][0]*_mv9e40k + A[1][1]*_zuj2xor + A[1][2]*_vork1dj;       _ikwe4bx[2] = A[2][0]*_mv9e40k + A[2][1]*_zuj2xor + A[2][2]*_vork1dj;

    return;
}

void      _jcslqc3(int _r11q7uw[1], double x,double y,double z,double _b2rvndf,double _dw79wcx, double _2pgyme6,double _rpkzggm,double _kjnawng,double _a41izxt)
{
    
    
    
    
    
    double a,            b,          c;

    a = ((_rpkzggm-_a41izxt)*(x-_kjnawng) +(_kjnawng-_2pgyme6)*(y-_a41izxt)) /  ((_rpkzggm-_a41izxt)*(_b2rvndf-_kjnawng) +(_kjnawng-_2pgyme6)*(_dw79wcx-_a41izxt));
    b = ((_a41izxt-_dw79wcx)*(x-_kjnawng) +(_b2rvndf-_kjnawng)*(y-_a41izxt)) /  ((_rpkzggm-_a41izxt)*(_b2rvndf-_kjnawng) +(_kjnawng-_2pgyme6)*(_dw79wcx-_a41izxt));
    c = 1.0 -a -b;

    _r11q7uw[0] = 1;
    if  ((a <= 0.0) && (b > c)    && (c > a))             {   _r11q7uw[0] = -1;                  }
    if  ((a > b)    && (b <= 0.0) && (c > a))             {   _r11q7uw[0] = -1;                  }
    if  ((a > b)    && (b > c)    && (c <= 0.0))          {   _r11q7uw[0] = -1;                  }
    if  ((a == 0.0) && (b >= 0.0) && (c >= 0.0))       {   _r11q7uw[0] = 0;                   }
    if  ((a >= 0.0) && (b == 0.0) && (c >= 0.0))       {   _r11q7uw[0] = 0;                   }
    if  ((a >= 0.0) && (b >= 0.0) && (c == 0.0))       {   _r11q7uw[0] = 0;                   }
    if  ((_r11q7uw[0] == 0) && (z != 0.0))              {   _r11q7uw[0] = 1;                   } 

    return;
}

void _5uofwb0(double x,double y,double z,double _52xg8p6,double _1750ceh,double _lxslg6m,double _iui9kkh,double _qa78rap, double _c54wsqt[3],double _6kof8kr[3],double _9iy4x41[6]) 
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
    
    
    _ikcidkw(x, _8nyfwdc, _v9nou3d, (-1.0*M_PI+_52xg8p6), _1750ceh, _m1zxyvg, _zlm72je, _qa78rap, e); 
    
    B[0][0] = 1.0;          B[0][1] = 0.0;      B[0][2] = 0.0;
    B[1][0] = 0.0;          B[1][1] = A[0][0];  B[1][2] = A[1][0];
    B[2][0] = 0.0;          B[2][1] = A[0][1];  B[2][2] = A[1][1];
    
    _olmlgse(_9iy4x41, e, B); 

    return;
}

void    _nnhj7xv(double _cjwtl9t[6],double _5lf6a4f[6], double X,double Y,double Z,double _1nouhtc,double _s1k9ye9,double _o2m71hs,double _b1qx6gv[3], double _wrdmly4[3],double _mu7jae5,double _odck828)
{   
    int   i,        I;
    double _qa78rap,      _slwf90k,   _yrybova,       _khejla0;
    
    double _6kof8kr[3],       _tatj9a3[3],          _vlnknd6[3],         _69ugnek[3],     _rcpprcj[3];
    double _gnzs43c[3],     _5ba75uu[3],   _or8wm9l[3],          _bb4m0a6[3];
    double _3ykyydu[6],           _3asnrm6[6],         _v1mzb3i[6];

    double A[3][3],          _xaji0y6[3][3];
   
    _qa78rap  = _odck828/(_mu7jae5+_odck828)/2.0; 
    
    _6kof8kr[0] = _wrdmly4[0]-_b1qx6gv[0];           _6kof8kr[1] = _wrdmly4[1]-_b1qx6gv[1];       _6kof8kr[2] = _wrdmly4[2]-_b1qx6gv[2];
    _tatj9a3[0]      = 0.0;                   _tatj9a3[1]      = 0.0;               _tatj9a3[2]      = 1.0;

    _yrybova   = sqrt( _6kof8kr[0]*_6kof8kr[0] +_6kof8kr[1]*_6kof8kr[1] +_6kof8kr[2]*_6kof8kr[2]);
    _khejla0   = -1.0*(_6kof8kr[0]*_tatj9a3[0]      +_6kof8kr[1]*_tatj9a3[1]      +_6kof8kr[2]*_tatj9a3[2]);    
    _slwf90k       = acos(_khejla0/_yrybova);
    if (fabs(cos(_slwf90k)/sin(_slwf90k)) > 5.0e+5*(M_PI/360.0))
   
    {   for (i = 0; i < 6; i++)         {       _cjwtl9t[i] = 0.0;        _5lf6a4f[i] = 0.0;            } 
    }
    else
    {
        _vlnknd6[0]   = _6kof8kr[0];          _vlnknd6[1]   = _6kof8kr[1];           _vlnknd6[2]   = 0.0;
        _yrybova = sqrt( _vlnknd6[0]*_vlnknd6[0] +_vlnknd6[1]*_vlnknd6[1] +_vlnknd6[2]*_vlnknd6[2]);
        
        _vlnknd6[0]  /= _yrybova;            _vlnknd6[1]  /= _yrybova;             _vlnknd6[2]  /= _yrybova;
        _rcpprcj[0]   = -1.0*_tatj9a3[0];          _rcpprcj[1]   = -1.0*_tatj9a3[1];           _rcpprcj[2]   = -1.0*_tatj9a3[2];
        
        _69ugnek[0]   = _rcpprcj[1]*_vlnknd6[2] -_rcpprcj[2]*_vlnknd6[1]; 
        _69ugnek[1]   = _rcpprcj[2]*_vlnknd6[0] -_rcpprcj[0]*_vlnknd6[2]; 
        _69ugnek[2]   = _rcpprcj[0]*_vlnknd6[1] -_rcpprcj[1]*_vlnknd6[0]; 
        _yrybova = sqrt( _69ugnek[0]*_69ugnek[0] +_69ugnek[1]*_69ugnek[1] +_69ugnek[2]*_69ugnek[2]);
        _69ugnek[0]  /= _yrybova;            _69ugnek[1]  /= _yrybova;             _69ugnek[2]  /= _yrybova;
        
        A[0][0]   = _vlnknd6[0];              A[0][1]   = _69ugnek[0];              A[0][2]   = _rcpprcj[0];
        A[1][0]   = _vlnknd6[1];              A[1][1]   = _69ugnek[1];              A[1][2]   = _rcpprcj[1];
        A[2][0]   = _vlnknd6[2];              A[2][1]   = _69ugnek[2];              A[2][2]   = _rcpprcj[2];
        
        _xaji0y6[0][0] = _vlnknd6[0];              _xaji0y6[0][1] = _vlnknd6[1];              _xaji0y6[0][2] = _vlnknd6[2];
        _xaji0y6[1][0] = _69ugnek[0];              _xaji0y6[1][1] = _69ugnek[1];              _xaji0y6[1][2] = _69ugnek[2];
        _xaji0y6[2][0] = _rcpprcj[0];              _xaji0y6[2][1] = _rcpprcj[1];              _xaji0y6[2][2] = _rcpprcj[2];
       
        
        _frhocn9(_or8wm9l, (X-_b1qx6gv[0]), (Y-_b1qx6gv[1]), (Z-_b1qx6gv[2]), A); 
        
        _frhocn9(_gnzs43c, _6kof8kr[0], _6kof8kr[1], _6kof8kr[2], A); 
        _bb4m0a6[0] = _or8wm9l[0] - _gnzs43c[0];
        _bb4m0a6[1] = _or8wm9l[1] - _gnzs43c[1];
        _bb4m0a6[2] = _or8wm9l[2] - _gnzs43c[2];
        
        _frhocn9(_5ba75uu, _1nouhtc, _s1k9ye9, _o2m71hs, A); 
        
        I = (_slwf90k*_or8wm9l[0]) >=0.0 ? 1 : 0;
        
        for (i = 0; i < 6; i++)         {       _3ykyydu[i] = 0.0;       _3asnrm6[i] = 0.0;           }
        
        if (I == 1)
        {   _4v30t9n(_or8wm9l[0],_or8wm9l[1], _or8wm9l[2], (-M_PI+_slwf90k),_5ba75uu[0],_5ba75uu[1], _5ba75uu[2], _qa78rap,(-1.0*_b1qx6gv[2]), _3ykyydu);
            _4v30t9n(_bb4m0a6[0],_bb4m0a6[1], _bb4m0a6[2], (-M_PI+_slwf90k),_5ba75uu[0],_5ba75uu[1], _5ba75uu[2], _qa78rap,(-1.0*_wrdmly4[2]), _3asnrm6);
        }
        
        else if (I == 0)
        {   _4v30t9n(_or8wm9l[0],_or8wm9l[1], _or8wm9l[2], _slwf90k,_5ba75uu[0],_5ba75uu[1], _5ba75uu[2], _qa78rap,(-1.0*_b1qx6gv[2]), _3ykyydu);
            _4v30t9n(_bb4m0a6[0],_bb4m0a6[1], _bb4m0a6[2], _slwf90k,_5ba75uu[0],_5ba75uu[1], _5ba75uu[2], _qa78rap,(-1.0*_wrdmly4[2]), _3asnrm6);
        }
        
         for (i = 0; i < 6; i++)         {           _v1mzb3i[i] = _3asnrm6[i] - _3ykyydu[i];            }
          
         
         _olmlgse(_5lf6a4f, _v1mzb3i, _xaji0y6);
               
         _cjwtl9t[0] = 2.0*_mu7jae5*_5lf6a4f[0] +_odck828*(_5lf6a4f[0] +_5lf6a4f[3] +_5lf6a4f[5]); 
              _cjwtl9t[3] = 2.0*_mu7jae5*_5lf6a4f[3] +_odck828*(_5lf6a4f[0] +_5lf6a4f[3] +_5lf6a4f[5]); 
         _cjwtl9t[5] = 2.0*_mu7jae5*_5lf6a4f[5] +_odck828*(_5lf6a4f[0] +_5lf6a4f[3] +_5lf6a4f[5]); 
         _cjwtl9t[1] = 2.0*_mu7jae5*_5lf6a4f[1]; 
         _cjwtl9t[2] = 2.0*_mu7jae5*_5lf6a4f[2]; 
          _cjwtl9t[4] = 2.0*_mu7jae5*_5lf6a4f[4]; 
    }
    return;
}

void _ikcidkw(double x, double y, double z, double _52xg8p6, double _1750ceh, double _lxslg6m, double _iui9kkh, double _qa78rap, double e[6]) 
{   double       _vyv2gqu,           _3p15ahd,           _ldalmje,            _xk2ptho;
    double       _zuj2xor,             _m4t8bp6,             _pr424tk,             _m3o55ix;
    double       r,              _fshb0oq,             _25lzse9,             _va1nkbt;
    double       _vshkoec,            W,              _qzw3c1q,             _xie0rnu;
    double       _kbhb7s9,            _jjq9k6r,            _tj1gm25,           C;
    double       S,              _ubods1e,         _eurxhw0,         _ic1l7dy;
    double       _0tpac56,            _t4ufel6,            _246hid2,            _5owf759;
    double       _35a0l72,            _53j2d54;
    double       _wyzx4w6,             _nlnqbu4,             _ufwbucj,           _7qqyzqi;
    double       _8drc11e,             _bvzb0oi;

    _vyv2gqu = sin(_52xg8p6);         _3p15ahd = cos(_52xg8p6);
    _ldalmje  = y*_3p15ahd - z*_vyv2gqu;     _xk2ptho = y*_vyv2gqu + z*_3p15ahd;
    _zuj2xor   = x*x;                 _m4t8bp6   = y*y;
    _pr424tk   = z*z;                 _m3o55ix   = _zuj2xor+_m4t8bp6+_pr424tk;
    r    = sqrt(_m3o55ix);           _fshb0oq   = r*r*r;
    _25lzse9   = r*(r-z);             _va1nkbt = _m3o55ix*(r-z)*(r-z);
    _vshkoec  = _fshb0oq*(r-z);
    
    W    = _xk2ptho-r;              _qzw3c1q   = W*W;
    _xie0rnu   = W*r;                 _kbhb7s9  = _qzw3c1q*r;
    _jjq9k6r  = W*_fshb0oq;                _tj1gm25 = _qzw3c1q*_m3o55ix;
    C    = (r*_3p15ahd -z)/_xie0rnu;      S    = (r*_vyv2gqu -y)/_xie0rnu;
    
    _wyzx4w6   = S*S;                 _nlnqbu4   = y*y*y;
    _8drc11e   = C*C;
    _ufwbucj = (1.0-_qa78rap);            _7qqyzqi = (2.0*_qa78rap+1.0);
    _bvzb0oi= _3p15ahd*_3p15ahd;

    
    _ubods1e = (_ldalmje/r/(r-_xk2ptho)-y/r/(r-z))/4.0/M_PI;
    _eurxhw0 = (x/r/(r-z)-_3p15ahd*x/r/(r-_xk2ptho))/4.0/M_PI;
    _ic1l7dy = (_vyv2gqu*x/r/(r-_xk2ptho))/4.0/M_PI;

    _0tpac56 = _1750ceh*(_ubods1e)                     + _1750ceh/8.0/M_PI/_ufwbucj*(_ldalmje/_xie0rnu+_ldalmje*_zuj2xor/_tj1gm25-_ldalmje*_zuj2xor/_jjq9k6r+y/_25lzse9-_zuj2xor*y/_va1nkbt-_zuj2xor*y/_vshkoec)                                            - _lxslg6m*x/8.0/M_PI/_ufwbucj*((_7qqyzqi/_xie0rnu+_zuj2xor/_tj1gm25-_zuj2xor/_jjq9k6r)*_3p15ahd+_7qqyzqi/_25lzse9-_zuj2xor/_va1nkbt-_zuj2xor/_vshkoec)                   + _iui9kkh*x*_vyv2gqu/8.0/M_PI/_ufwbucj*(_7qqyzqi/_xie0rnu+_zuj2xor/_tj1gm25-_zuj2xor/_jjq9k6r);

    _5owf759 = _lxslg6m*(_eurxhw0)                     + _1750ceh/8.0/M_PI/_ufwbucj*((1.0/_xie0rnu+_wyzx4w6-_m4t8bp6/_jjq9k6r)*_ldalmje+_7qqyzqi*y/_25lzse9-_nlnqbu4/_va1nkbt-_nlnqbu4/_vshkoec-2.0*_qa78rap*_3p15ahd*S)                                    - _lxslg6m*x/8.0/M_PI/_ufwbucj*(1.0/_25lzse9-_m4t8bp6/_va1nkbt-_m4t8bp6/_vshkoec+(1.0/_xie0rnu+_wyzx4w6-_m4t8bp6/_jjq9k6r)*_3p15ahd)                          + _iui9kkh*x*_vyv2gqu/8.0/M_PI/_ufwbucj*(1.0/_xie0rnu+_wyzx4w6-_m4t8bp6/_jjq9k6r);

    _53j2d54 = _iui9kkh*(_ic1l7dy)                     + _1750ceh/8.0/M_PI/_ufwbucj*(_ldalmje/W/r+_ldalmje*_8drc11e-_ldalmje*_pr424tk/_jjq9k6r+y*z/_fshb0oq+2.0*_qa78rap*_vyv2gqu*C)                                                   - _lxslg6m*x/8.0/M_PI/_ufwbucj*((1.0/_xie0rnu+_8drc11e-_pr424tk/_jjq9k6r)*_3p15ahd+z/_fshb0oq)                                           + _iui9kkh*x*_vyv2gqu/8.0/M_PI/_ufwbucj*(1.0/_xie0rnu+_8drc11e-_pr424tk/_jjq9k6r);

    _t4ufel6 = _1750ceh*(_eurxhw0)/2.0+_lxslg6m*(_ubods1e)/2.0 - _1750ceh/8.0/M_PI/_ufwbucj*(x*_m4t8bp6/_va1nkbt-_qa78rap*x/_25lzse9+x*_m4t8bp6/_vshkoec-_qa78rap*x*_3p15ahd/_xie0rnu+_ldalmje*x*S/_xie0rnu+_ldalmje*x*y/_jjq9k6r)                                   + _lxslg6m/8.0/M_PI/_ufwbucj*(_zuj2xor*y/_va1nkbt-_qa78rap*y/_25lzse9+_zuj2xor*y/_vshkoec+_qa78rap*_3p15ahd*S+_zuj2xor*y*_3p15ahd/_jjq9k6r+_zuj2xor*_3p15ahd*S/_xie0rnu)          - _iui9kkh*_vyv2gqu/8.0/M_PI/_ufwbucj*(_qa78rap*S+_zuj2xor*S/_xie0rnu+_zuj2xor*y/_jjq9k6r);

    _246hid2 = _1750ceh*(_ic1l7dy)/2.0+_iui9kkh*(_ubods1e)/2.0 - _1750ceh/8.0/M_PI/_ufwbucj*(-x*y/_fshb0oq+_qa78rap*x*_vyv2gqu/_xie0rnu+_ldalmje*x*C/_xie0rnu+_ldalmje*x*z/_jjq9k6r)                                                      + _lxslg6m/8.0/M_PI/_ufwbucj*(-_zuj2xor/_fshb0oq+_qa78rap/r+_qa78rap*_3p15ahd*C+_zuj2xor*z*_3p15ahd/_jjq9k6r+_zuj2xor*_3p15ahd*C/_xie0rnu)                         - _iui9kkh*_vyv2gqu/8.0/M_PI/_ufwbucj*(_qa78rap*C+_zuj2xor*C/_xie0rnu+_zuj2xor*z/_jjq9k6r);

    _35a0l72 = _lxslg6m*(_ic1l7dy)/2.0+_iui9kkh*(_eurxhw0)/2.0 + _1750ceh/8.0/M_PI/_ufwbucj*(_m4t8bp6/_fshb0oq-_qa78rap/r-_qa78rap*_3p15ahd*C+_qa78rap*_vyv2gqu*S+_ldalmje*_vyv2gqu*_3p15ahd/_qzw3c1q-_ldalmje*(y*_3p15ahd+z*_vyv2gqu)/_kbhb7s9+_ldalmje*y*z/_tj1gm25-_ldalmje*y*z/_jjq9k6r) - _lxslg6m*x/8.0/M_PI/_ufwbucj*(y/_fshb0oq+_vyv2gqu*_bvzb0oi/_qzw3c1q-_3p15ahd*(y*_3p15ahd+z*_vyv2gqu)/_kbhb7s9+y*z*_3p15ahd/_tj1gm25-y*z*_3p15ahd/_jjq9k6r) - _iui9kkh*x*_vyv2gqu/8.0/M_PI/_ufwbucj*(y*z/_jjq9k6r-_vyv2gqu*_3p15ahd/_qzw3c1q+(y*_3p15ahd+z*_vyv2gqu)/_kbhb7s9-y*z/_tj1gm25);

    e[0] = _0tpac56;             e[1] = _t4ufel6;             e[2] = _246hid2;
    e[3] = _5owf759;             e[4] = _35a0l72;             e[5] = _53j2d54;
    return;
}

void  _4v30t9n(double _8nyfwdc, double _m4t8bp6, double _nlnqbu4, double _slwf90k, double _1rsrla8, double _s8hmbr5, double _xqt4p50, double _qa78rap, double a, double _5lf6a4f[6])
{
    
    double _7dknncp,         _h5csk30,           _pr97y2l,           _kzrezrd,        _yms1mwn;
    double _b9oi0d4,          _gsq7ait,            _in9dl5t,             _0lbt3g2,         _qzw3c1q;
    double _len606n,           _9rr7s20,             _cv3tlcx,             _m9j98cw,         _n2y8fv2;
    double _bc2nnuv,           _y8n1uto,             _b7r3ocw,             _w9pfb3h,   _pmlvak9;
    double _n14eb87,     _4fcl5ef,            _jp6mkt7,            _s1gd1t5,        _dfhg6w3;
    double _stt1eqc,          _37rln36;
    
    double _mf9pml9,       _iwg0odw,          _6c2cu7o,               _h17j8qo,          _nz3jbck;
    double _h4a3yjt,           _mqx5ggo,          _09vs9kk,               _ln1p0vu,          _u7fbmen;
    double _88d0tgs,     _yt91qe3,         _zv4u3bb,          _112hxrj,     _mofmmwi;
    double _yc3rwqm,       _wgabtgd,          _lvpymb5,           _yxyur4s,     _jz0ifrm;
     
    _7dknncp = sin(_slwf90k);                  _h5csk30 = cos(_slwf90k);               _pr97y2l = _h5csk30/_7dknncp;
    _kzrezrd  = _nlnqbu4+a+a;                      _yms1mwn  = _8nyfwdc*_h5csk30+_kzrezrd*_7dknncp;       _b9oi0d4  = -_8nyfwdc*_7dknncp+_kzrezrd*_h5csk30;
    _gsq7ait  = _8nyfwdc*_8nyfwdc + _m4t8bp6*_m4t8bp6 + _kzrezrd*_kzrezrd;     _in9dl5t   = sqrt(_gsq7ait);              _b7r3ocw   = 1.0-_qa78rap-_qa78rap;

    _0lbt3g2   = _in9dl5t*_h5csk30+_kzrezrd;                 _qzw3c1q   = _h5csk30  +a/_in9dl5t;            _len606n   = _h5csk30 +_kzrezrd/_in9dl5t;
    _9rr7s20   = _qa78rap +a/_in9dl5t;                    _cv3tlcx   = _qa78rap+_qa78rap +a/_in9dl5t;            _m9j98cw   = _in9dl5t   +_kzrezrd;
    _n2y8fv2   = _in9dl5t +_b9oi0d4;                     _bc2nnuv   = _nlnqbu4 +a;                  _y8n1uto   = 1.0  +a/_in9dl5t/_h5csk30;

    
    _w9pfb3h =      _yms1mwn/_in9dl5t/(_in9dl5t+_b9oi0d4)      -_8nyfwdc/_in9dl5t/(_in9dl5t+_kzrezrd); 
    _pmlvak9 =       _m4t8bp6/_in9dl5t/(_in9dl5t+_kzrezrd) -_h5csk30*_m4t8bp6/_in9dl5t/(_in9dl5t+_b9oi0d4); 
    _n14eb87 = -_7dknncp*_m4t8bp6/_in9dl5t/(_in9dl5t+_b9oi0d4);                      
 
    _mf9pml9   = _in9dl5t*_in9dl5t*_in9dl5t;                        _iwg0odw = _pr97y2l*_pr97y2l;                            _6c2cu7o   = _m9j98cw*_m9j98cw;
    _h17j8qo   = _8nyfwdc*_8nyfwdc;                            _nz3jbck   = _n2y8fv2*_n2y8fv2;                                 _h4a3yjt   = _m4t8bp6*_m4t8bp6;
    _mqx5ggo  = _gsq7ait*_gsq7ait;                           _09vs9kk   = _8nyfwdc*_8nyfwdc*_8nyfwdc;                             _ln1p0vu   = _in9dl5t*_in9dl5t*_in9dl5t*_in9dl5t*_in9dl5t;
    _u7fbmen   = _m4t8bp6*_m4t8bp6*_m4t8bp6;                        _88d0tgs = _h5csk30*_h5csk30;                       _112hxrj   = (_8nyfwdc/_in9dl5t-_7dknncp);
    
   
    _yt91qe3   = (2.0-_qa78rap-_qa78rap);                  _zv4u3bb    = (a/_gsq7ait+1.0/_m9j98cw);                  _mofmmwi   = (_in9dl5t*_7dknncp-_8nyfwdc);
    _yc3rwqm  = (-2.0+_qa78rap+_qa78rap);                 _wgabtgd     = (1.0-_cv3tlcx);                        _lvpymb5     = (_cv3tlcx-1.0);
    _yxyur4s  = (_kzrezrd/_in9dl5t+1.0);                 _jz0ifrm     = (1.0-_qa78rap);

    _4fcl5ef = _1rsrla8*(0.25*(_yc3rwqm*_b7r3ocw*_pmlvak9*_iwg0odw-_b7r3ocw*_m4t8bp6/_6c2cu7o*(_wgabtgd*_pr97y2l-_8nyfwdc/_m9j98cw*_9rr7s20)/_in9dl5t*_8nyfwdc+_b7r3ocw*_m4t8bp6/_m9j98cw*(a/_mf9pml9*_8nyfwdc*_pr97y2l-1.0/_m9j98cw*_9rr7s20+_h17j8qo/_6c2cu7o*_9rr7s20/_in9dl5t+_h17j8qo/_m9j98cw*a/_mf9pml9)-_b7r3ocw*_m4t8bp6*_h5csk30*_pr97y2l/_nz3jbck*_qzw3c1q*_112hxrj-_b7r3ocw*_m4t8bp6*_h5csk30*_pr97y2l/_n2y8fv2*a/_mf9pml9*_8nyfwdc-3.0*a*_m4t8bp6*_bc2nnuv*_pr97y2l/_ln1p0vu*_8nyfwdc-_m4t8bp6*_bc2nnuv/_mf9pml9/_m9j98cw*(-_b7r3ocw*_pr97y2l+_8nyfwdc/_m9j98cw*_cv3tlcx+a*_8nyfwdc/_gsq7ait)*_8nyfwdc-_m4t8bp6*_bc2nnuv/_gsq7ait/_6c2cu7o*(-_b7r3ocw*_pr97y2l+_8nyfwdc/_m9j98cw*_cv3tlcx+a*_8nyfwdc/_gsq7ait)*_8nyfwdc+_m4t8bp6*_bc2nnuv/_in9dl5t/_m9j98cw*(1.0/_m9j98cw*_cv3tlcx-_h17j8qo/_6c2cu7o*_cv3tlcx/_in9dl5t-_h17j8qo/_m9j98cw*a/_mf9pml9+a/_gsq7ait-2.0*a*_h17j8qo/_mqx5ggo)-_m4t8bp6*_bc2nnuv/_mf9pml9/_n2y8fv2*(_h5csk30/_n2y8fv2*(_0lbt3g2*(_b7r3ocw*_h5csk30-a/_in9dl5t)*_pr97y2l+_yt91qe3*_mofmmwi*_h5csk30)-a*_kzrezrd*_h5csk30*_pr97y2l/_gsq7ait)*_8nyfwdc-_m4t8bp6*_bc2nnuv/_in9dl5t/_nz3jbck*(_h5csk30/_n2y8fv2*(_0lbt3g2*(_b7r3ocw*_h5csk30-a/_in9dl5t)*_pr97y2l+_yt91qe3*_mofmmwi*_h5csk30)-a*_kzrezrd*_h5csk30*_pr97y2l/_gsq7ait)*_112hxrj+_m4t8bp6*_bc2nnuv/_in9dl5t/_n2y8fv2*(-_h5csk30/_nz3jbck*(_0lbt3g2*(_b7r3ocw*_h5csk30-a/_in9dl5t)*_pr97y2l+_yt91qe3*_mofmmwi*_h5csk30)*_112hxrj+_h5csk30/_n2y8fv2*(1.0/_in9dl5t*_h5csk30*_8nyfwdc*(_b7r3ocw*_h5csk30-a/_in9dl5t)*_pr97y2l+_0lbt3g2*a/_mf9pml9*_8nyfwdc*_pr97y2l+_yt91qe3*(1.0/_in9dl5t*_7dknncp*_8nyfwdc-1.0)*_h5csk30)+2.0*a*_kzrezrd*_h5csk30*_pr97y2l/_mqx5ggo*_8nyfwdc))/M_PI/_jz0ifrm)
        + _s8hmbr5*(0.25*(_b7r3ocw*((_yt91qe3*_iwg0odw+_qa78rap)/_in9dl5t*_8nyfwdc/_m9j98cw-(_yt91qe3*_iwg0odw+1.0)*_h5csk30*_112hxrj/_n2y8fv2)-_b7r3ocw/_6c2cu7o*(-_b7r3ocw*_8nyfwdc*_pr97y2l+_qa78rap*_kzrezrd-a+a*_8nyfwdc*_pr97y2l/_in9dl5t+_h17j8qo/_m9j98cw*_9rr7s20)/_in9dl5t*_8nyfwdc+_b7r3ocw/_m9j98cw*(-_b7r3ocw*_pr97y2l+a*_pr97y2l/_in9dl5t-a*_h17j8qo*_pr97y2l/_mf9pml9+2.0*_8nyfwdc/_m9j98cw*_9rr7s20-_09vs9kk/_6c2cu7o*_9rr7s20/_in9dl5t-_09vs9kk/_m9j98cw*a/_mf9pml9)+_b7r3ocw*_pr97y2l/_nz3jbck*(_yms1mwn*_h5csk30-a*_mofmmwi/_in9dl5t/_h5csk30)*_112hxrj-_b7r3ocw*_pr97y2l/_n2y8fv2*(_88d0tgs-a*(1.0/_in9dl5t*_7dknncp*_8nyfwdc-1.0)/_in9dl5t/_h5csk30+a*_mofmmwi/_mf9pml9/_h5csk30*_8nyfwdc)-a*_bc2nnuv*_pr97y2l/_mf9pml9+3.0*a*_h17j8qo*_bc2nnuv*_pr97y2l/_ln1p0vu-_bc2nnuv/_6c2cu7o*(2.0*_qa78rap+1.0/_in9dl5t*(_b7r3ocw*_8nyfwdc*_pr97y2l+a)-_h17j8qo/_in9dl5t/_m9j98cw*_cv3tlcx-a*_h17j8qo/_mf9pml9)/_in9dl5t*_8nyfwdc+_bc2nnuv/_m9j98cw*(-1.0/_mf9pml9*(_b7r3ocw*_8nyfwdc*_pr97y2l+a)*_8nyfwdc+1.0/_in9dl5t*_b7r3ocw*_pr97y2l-2.0*_8nyfwdc/_in9dl5t/_m9j98cw*_cv3tlcx+_09vs9kk/_mf9pml9/_m9j98cw*_cv3tlcx+_09vs9kk/_gsq7ait/_6c2cu7o*_cv3tlcx+_09vs9kk/_mqx5ggo/_m9j98cw*a-2.0*a/_mf9pml9*_8nyfwdc+3.0*a*_09vs9kk/_ln1p0vu)-_bc2nnuv*_pr97y2l/_nz3jbck*(-_h5csk30*_7dknncp+a*_8nyfwdc*_kzrezrd/_mf9pml9/_h5csk30+_mofmmwi/_in9dl5t*(_yt91qe3*_h5csk30-_0lbt3g2/_n2y8fv2*_y8n1uto))*_112hxrj+_bc2nnuv*_pr97y2l/_n2y8fv2*(a*_kzrezrd/_mf9pml9/_h5csk30-3.0*a*_h17j8qo*_kzrezrd/_ln1p0vu/_h5csk30+(1.0/_in9dl5t*_7dknncp*_8nyfwdc-1.0)/_in9dl5t*(_yt91qe3*_h5csk30-_0lbt3g2/_n2y8fv2*_y8n1uto)-_mofmmwi/_mf9pml9*(_yt91qe3*_h5csk30-_0lbt3g2/_n2y8fv2*_y8n1uto)*_8nyfwdc+_mofmmwi/_in9dl5t*(-1.0/_in9dl5t*_h5csk30*_8nyfwdc/_n2y8fv2*_y8n1uto+_0lbt3g2/_nz3jbck*_y8n1uto*_112hxrj+_0lbt3g2/_n2y8fv2*a/_mf9pml9/_h5csk30*_8nyfwdc)))/M_PI/_jz0ifrm)
        + _xqt4p50*(0.25*(_b7r3ocw*(-_m4t8bp6/_6c2cu7o*(1.0+a/_in9dl5t)/_in9dl5t*_8nyfwdc-_m4t8bp6/_m9j98cw*a/_mf9pml9*_8nyfwdc+_m4t8bp6*_h5csk30/_nz3jbck*_qzw3c1q*_112hxrj+_m4t8bp6*_h5csk30/_n2y8fv2*a/_mf9pml9*_8nyfwdc)+_m4t8bp6*_bc2nnuv/_mf9pml9*_zv4u3bb*_8nyfwdc-_m4t8bp6*_bc2nnuv/_in9dl5t*(-2.0*a/_mqx5ggo*_8nyfwdc-1.0/_6c2cu7o/_in9dl5t*_8nyfwdc)-_m4t8bp6*_bc2nnuv*_h5csk30/_mf9pml9/_n2y8fv2*(_0lbt3g2/_n2y8fv2*_qzw3c1q+a*_kzrezrd/_gsq7ait)*_8nyfwdc-_m4t8bp6*_bc2nnuv*_h5csk30/_in9dl5t/_nz3jbck*(_0lbt3g2/_n2y8fv2*_qzw3c1q+a*_kzrezrd/_gsq7ait)*_112hxrj+_m4t8bp6*_bc2nnuv*_h5csk30/_in9dl5t/_n2y8fv2*(1.0/_in9dl5t*_h5csk30*_8nyfwdc/_n2y8fv2*_qzw3c1q-_0lbt3g2/_nz3jbck*_qzw3c1q*_112hxrj-_0lbt3g2/_n2y8fv2*a/_mf9pml9*_8nyfwdc-2.0*a*_kzrezrd/_mqx5ggo*_8nyfwdc))/M_PI/_jz0ifrm);
        
    _dfhg6w3 = _1rsrla8*(0.25*(_b7r3ocw*((_yt91qe3*_iwg0odw-_qa78rap)/_in9dl5t*_m4t8bp6/_m9j98cw-(_yt91qe3*_iwg0odw+1.0-2.0*_qa78rap)*_h5csk30/_in9dl5t*_m4t8bp6/_n2y8fv2)+_b7r3ocw/_6c2cu7o*(_8nyfwdc*_pr97y2l*_wgabtgd+_qa78rap*_kzrezrd-a+_h4a3yjt/_m9j98cw*_9rr7s20)/_in9dl5t*_m4t8bp6-_b7r3ocw/_m9j98cw*(a*_8nyfwdc*_pr97y2l/_mf9pml9*_m4t8bp6+2.0*_m4t8bp6/_m9j98cw*_9rr7s20-_u7fbmen/_6c2cu7o*_9rr7s20/_in9dl5t-_u7fbmen/_m9j98cw*a/_mf9pml9)+_b7r3ocw*_yms1mwn*_pr97y2l/_nz3jbck*_qzw3c1q/_in9dl5t*_m4t8bp6+_b7r3ocw*_yms1mwn*_pr97y2l/_n2y8fv2*a/_mf9pml9*_m4t8bp6+3.0*a*_m4t8bp6*_bc2nnuv*_pr97y2l/_ln1p0vu*_8nyfwdc-_bc2nnuv/_6c2cu7o*(-2.0*_qa78rap+1.0/_in9dl5t*(_b7r3ocw*_8nyfwdc*_pr97y2l-a)+_h4a3yjt/_in9dl5t/_m9j98cw*_cv3tlcx+a*_h4a3yjt/_mf9pml9)/_in9dl5t*_m4t8bp6+_bc2nnuv/_m9j98cw*(-1.0/_mf9pml9*(_b7r3ocw*_8nyfwdc*_pr97y2l-a)*_m4t8bp6+2.0*_m4t8bp6/_in9dl5t/_m9j98cw*_cv3tlcx-_u7fbmen/_mf9pml9/_m9j98cw*_cv3tlcx-_u7fbmen/_gsq7ait/_6c2cu7o*_cv3tlcx-_u7fbmen/_mqx5ggo/_m9j98cw*a+2.0*a/_mf9pml9*_m4t8bp6-3.0*a*_u7fbmen/_ln1p0vu)-_bc2nnuv/_nz3jbck*(_88d0tgs-1.0/_in9dl5t*(_b7r3ocw*_yms1mwn*_pr97y2l+a*_h5csk30)+a*_kzrezrd*_yms1mwn*_pr97y2l/_mf9pml9-1.0/_in9dl5t/_n2y8fv2*(_h4a3yjt*_88d0tgs-a*_yms1mwn*_pr97y2l/_in9dl5t*_0lbt3g2))/_in9dl5t*_m4t8bp6+_bc2nnuv/_n2y8fv2*(1.0/_mf9pml9*(_b7r3ocw*_yms1mwn*_pr97y2l+a*_h5csk30)*_m4t8bp6-3.0*a*_kzrezrd*_yms1mwn*_pr97y2l/_ln1p0vu*_m4t8bp6+1.0/_mf9pml9/_n2y8fv2*(_h4a3yjt*_88d0tgs-a*_yms1mwn*_pr97y2l/_in9dl5t*_0lbt3g2)*_m4t8bp6+1.0/_gsq7ait/_nz3jbck*(_h4a3yjt*_88d0tgs-a*_yms1mwn*_pr97y2l/_in9dl5t*_0lbt3g2)*_m4t8bp6-1.0/_in9dl5t/_n2y8fv2*(2.0*_m4t8bp6*_88d0tgs+a*_yms1mwn*_pr97y2l/_mf9pml9*_0lbt3g2*_m4t8bp6-a*_yms1mwn*_pr97y2l/_gsq7ait*_h5csk30*_m4t8bp6)))/M_PI/_jz0ifrm)
        + _s8hmbr5*(0.25*(_yt91qe3*_b7r3ocw*_w9pfb3h*_iwg0odw+_b7r3ocw/_m9j98cw*(_lvpymb5*_pr97y2l+_8nyfwdc/_m9j98cw*_9rr7s20)-_b7r3ocw*_h4a3yjt/_6c2cu7o*(_lvpymb5*_pr97y2l+_8nyfwdc/_m9j98cw*_9rr7s20)/_in9dl5t+_b7r3ocw*_m4t8bp6/_m9j98cw*(-a/_mf9pml9*_m4t8bp6*_pr97y2l-_8nyfwdc/_6c2cu7o*_9rr7s20/_in9dl5t*_m4t8bp6-_m4t8bp6/_m9j98cw*a/_mf9pml9*_8nyfwdc)-_b7r3ocw*_pr97y2l/_n2y8fv2*_y8n1uto+_b7r3ocw*_h4a3yjt*_pr97y2l/_nz3jbck*_y8n1uto/_in9dl5t+_b7r3ocw*_h4a3yjt*_pr97y2l/_n2y8fv2*a/_mf9pml9/_h5csk30-a*_bc2nnuv*_pr97y2l/_mf9pml9+3.0*a*_h4a3yjt*_bc2nnuv*_pr97y2l/_ln1p0vu+_bc2nnuv/_in9dl5t/_m9j98cw*(_b7r3ocw*_pr97y2l-2.0*_qa78rap*_8nyfwdc/_m9j98cw-a*_8nyfwdc/_in9dl5t*(1.0/_in9dl5t+1.0/_m9j98cw))-_h4a3yjt*_bc2nnuv/_mf9pml9/_m9j98cw*(_b7r3ocw*_pr97y2l-2.0*_qa78rap*_8nyfwdc/_m9j98cw-a*_8nyfwdc/_in9dl5t*(1.0/_in9dl5t+1.0/_m9j98cw))-_h4a3yjt*_bc2nnuv/_gsq7ait/_6c2cu7o*(_b7r3ocw*_pr97y2l-2.0*_qa78rap*_8nyfwdc/_m9j98cw-a*_8nyfwdc/_in9dl5t*(1.0/_in9dl5t+1.0/_m9j98cw))+_m4t8bp6*_bc2nnuv/_in9dl5t/_m9j98cw*(2.0*_qa78rap*_8nyfwdc/_6c2cu7o/_in9dl5t*_m4t8bp6+a*_8nyfwdc/_mf9pml9*(1.0/_in9dl5t+1.0/_m9j98cw)*_m4t8bp6-a*_8nyfwdc/_in9dl5t*(-1.0/_mf9pml9*_m4t8bp6-1.0/_6c2cu7o/_in9dl5t*_m4t8bp6))+_bc2nnuv*_pr97y2l/_in9dl5t/_n2y8fv2*(_yc3rwqm*_h5csk30+_0lbt3g2/_n2y8fv2*_y8n1uto+a*_kzrezrd/_gsq7ait/_h5csk30)-_h4a3yjt*_bc2nnuv*_pr97y2l/_mf9pml9/_n2y8fv2*(_yc3rwqm*_h5csk30+_0lbt3g2/_n2y8fv2*_y8n1uto+a*_kzrezrd/_gsq7ait/_h5csk30)-_h4a3yjt*_bc2nnuv*_pr97y2l/_gsq7ait/_nz3jbck*(_yc3rwqm*_h5csk30+_0lbt3g2/_n2y8fv2*_y8n1uto+a*_kzrezrd/_gsq7ait/_h5csk30)+_m4t8bp6*_bc2nnuv*_pr97y2l/_in9dl5t/_n2y8fv2*(1.0/_in9dl5t*_h5csk30*_m4t8bp6/_n2y8fv2*_y8n1uto-_0lbt3g2/_nz3jbck*_y8n1uto/_in9dl5t*_m4t8bp6-_0lbt3g2/_n2y8fv2*a/_mf9pml9/_h5csk30*_m4t8bp6-2.0*a*_kzrezrd/_mqx5ggo/_h5csk30*_m4t8bp6))/M_PI/_jz0ifrm)
        + _xqt4p50*(0.25*(_b7r3ocw*(-_7dknncp/_in9dl5t*_m4t8bp6/_n2y8fv2+_m4t8bp6/_6c2cu7o*(1.0+a/_in9dl5t)/_in9dl5t*_8nyfwdc+_m4t8bp6/_m9j98cw*a/_mf9pml9*_8nyfwdc-_yms1mwn/_nz3jbck*_qzw3c1q/_in9dl5t*_m4t8bp6-_yms1mwn/_n2y8fv2*a/_mf9pml9*_m4t8bp6)-_m4t8bp6*_bc2nnuv/_mf9pml9*_zv4u3bb*_8nyfwdc+_8nyfwdc*_bc2nnuv/_in9dl5t*(-2.0*a/_mqx5ggo*_m4t8bp6-1.0/_6c2cu7o/_in9dl5t*_m4t8bp6)+_bc2nnuv/_nz3jbck*(_7dknncp*(_h5csk30-a/_in9dl5t)+_yms1mwn/_in9dl5t*(1.0+a*_kzrezrd/_gsq7ait)-1.0/_in9dl5t/_n2y8fv2*(_h4a3yjt*_h5csk30*_7dknncp-a*_yms1mwn/_in9dl5t*_0lbt3g2))/_in9dl5t*_m4t8bp6-_bc2nnuv/_n2y8fv2*(_7dknncp*a/_mf9pml9*_m4t8bp6-_yms1mwn/_mf9pml9*(1.0+a*_kzrezrd/_gsq7ait)*_m4t8bp6-2.0*_yms1mwn/_ln1p0vu*a*_kzrezrd*_m4t8bp6+1.0/_mf9pml9/_n2y8fv2*(_h4a3yjt*_h5csk30*_7dknncp-a*_yms1mwn/_in9dl5t*_0lbt3g2)*_m4t8bp6+1.0/_gsq7ait/_nz3jbck*(_h4a3yjt*_h5csk30*_7dknncp-a*_yms1mwn/_in9dl5t*_0lbt3g2)*_m4t8bp6-1/_in9dl5t/_n2y8fv2*(2*_m4t8bp6*_h5csk30*_7dknncp+a*_yms1mwn/_mf9pml9*_0lbt3g2*_m4t8bp6-a*_yms1mwn/_gsq7ait*_h5csk30*_m4t8bp6)))/M_PI/_jz0ifrm);

    _37rln36 = _1rsrla8*(0.25*(_yt91qe3*(_b7r3ocw*_n14eb87*_pr97y2l-_m4t8bp6/_6c2cu7o*_cv3tlcx*_yxyur4s-0.5*_m4t8bp6/_m9j98cw*a/_mf9pml9*2.0*_kzrezrd+_m4t8bp6*_h5csk30/_nz3jbck*_qzw3c1q*_len606n+0.5*_m4t8bp6*_h5csk30/_n2y8fv2*a/_mf9pml9*2.0*_kzrezrd)+_m4t8bp6/_in9dl5t*(2.0*_qa78rap/_m9j98cw+a/_gsq7ait)-0.5*_m4t8bp6*_bc2nnuv/_mf9pml9*(2.0*_qa78rap/_m9j98cw+a/_gsq7ait)*2.0*_kzrezrd+_m4t8bp6*_bc2nnuv/_in9dl5t*(-2.0*_qa78rap/_6c2cu7o*_yxyur4s-a/_mqx5ggo*2.0*_kzrezrd)+_m4t8bp6*_h5csk30/_in9dl5t/_n2y8fv2*(1.0-2.0*_qa78rap-_0lbt3g2/_n2y8fv2*_qzw3c1q-a*_kzrezrd/_gsq7ait)-0.5*_m4t8bp6*_bc2nnuv*_h5csk30/_mf9pml9/_n2y8fv2*(1.0-2.0*_qa78rap-_0lbt3g2/_n2y8fv2*_qzw3c1q-a*_kzrezrd/_gsq7ait)*2.0*_kzrezrd-_m4t8bp6*_bc2nnuv*_h5csk30/_in9dl5t/_nz3jbck*(1.0-2.0*_qa78rap-_0lbt3g2/_n2y8fv2*_qzw3c1q-a*_kzrezrd/_gsq7ait)*_len606n+_m4t8bp6*_bc2nnuv*_h5csk30/_in9dl5t/_n2y8fv2*(-(_h5csk30*_kzrezrd/_in9dl5t+1.0)/_n2y8fv2*_qzw3c1q+_0lbt3g2/_nz3jbck*_qzw3c1q*_len606n+0.5*_0lbt3g2/_n2y8fv2*a/_mf9pml9*2.0*_kzrezrd-a/_gsq7ait+a*_kzrezrd/_mqx5ggo*2.0*_kzrezrd))/M_PI/_jz0ifrm)
        + _s8hmbr5*(0.25*(_yc3rwqm*_b7r3ocw*_pr97y2l*(_yxyur4s/_m9j98cw-_h5csk30*_len606n/_n2y8fv2)+_yt91qe3*_8nyfwdc/_6c2cu7o*_cv3tlcx*_yxyur4s+0.5*_yt91qe3*_8nyfwdc/_m9j98cw*a/_mf9pml9*2.0*_kzrezrd+_yt91qe3*_7dknncp/_n2y8fv2*_qzw3c1q-_yt91qe3*_yms1mwn/_nz3jbck*_qzw3c1q*_len606n-0.5*_yt91qe3*_yms1mwn/_n2y8fv2*a/_mf9pml9*2.0*_kzrezrd+1.0/_in9dl5t*(_b7r3ocw*_pr97y2l-2.0*_qa78rap*_8nyfwdc/_m9j98cw-a*_8nyfwdc/_gsq7ait)-0.5*_bc2nnuv/_mf9pml9*(_b7r3ocw*_pr97y2l-2.0*_qa78rap*_8nyfwdc/_m9j98cw-a*_8nyfwdc/_gsq7ait)*2.0*_kzrezrd+_bc2nnuv/_in9dl5t*(2.0*_qa78rap*_8nyfwdc/_6c2cu7o*_yxyur4s+a*_8nyfwdc/_mqx5ggo*2.0*_kzrezrd)-1.0/_n2y8fv2*(_h5csk30*_7dknncp+_0lbt3g2*_pr97y2l/_in9dl5t*(_yt91qe3*_h5csk30-_0lbt3g2/_n2y8fv2)+a/_in9dl5t*(_7dknncp-_kzrezrd*_yms1mwn/_gsq7ait-_yms1mwn*_0lbt3g2/_in9dl5t/_n2y8fv2))+_bc2nnuv/_nz3jbck*(_h5csk30*_7dknncp+_0lbt3g2*_pr97y2l/_in9dl5t*(_yt91qe3*_h5csk30-_0lbt3g2/_n2y8fv2)+a/_in9dl5t*(_7dknncp-_kzrezrd*_yms1mwn/_gsq7ait-_yms1mwn*_0lbt3g2/_in9dl5t/_n2y8fv2))*_len606n-_bc2nnuv/_n2y8fv2*((_h5csk30*_kzrezrd/_in9dl5t+1.0)*_pr97y2l/_in9dl5t*(_yt91qe3*_h5csk30-_0lbt3g2/_n2y8fv2)-0.5*_0lbt3g2*_pr97y2l/_mf9pml9*(_yt91qe3*_h5csk30-_0lbt3g2/_n2y8fv2)*2.0*_kzrezrd+_0lbt3g2*_pr97y2l/_in9dl5t*(-(_h5csk30*_kzrezrd/_in9dl5t+1.0)/_n2y8fv2+_0lbt3g2/_nz3jbck*_len606n)-0.5*a/_mf9pml9*(_7dknncp-_kzrezrd*_yms1mwn/_gsq7ait-_yms1mwn*_0lbt3g2/_in9dl5t/_n2y8fv2)*2.0*_kzrezrd+a/_in9dl5t*(-_yms1mwn/_gsq7ait-_kzrezrd*_7dknncp/_gsq7ait+_kzrezrd*_yms1mwn/_mqx5ggo*2.0*_kzrezrd-_7dknncp*_0lbt3g2/_in9dl5t/_n2y8fv2-_yms1mwn*(_h5csk30*_kzrezrd/_in9dl5t+1.0)/_in9dl5t/_n2y8fv2+0.5*_yms1mwn*_0lbt3g2/_mf9pml9/_n2y8fv2*2.0*_kzrezrd+_yms1mwn*_0lbt3g2/_in9dl5t/_nz3jbck*_len606n)))/M_PI/_jz0ifrm)
        + _xqt4p50*(0.25*(_yt91qe3*_n14eb87-_yt91qe3*_m4t8bp6*_7dknncp/_nz3jbck*_qzw3c1q*_len606n-0.5*_yt91qe3*_m4t8bp6*_7dknncp/_n2y8fv2*a/_mf9pml9*2.0*_kzrezrd+_m4t8bp6*_7dknncp/_in9dl5t/_n2y8fv2*(1.0+_0lbt3g2/_n2y8fv2*_qzw3c1q+a*_kzrezrd/_gsq7ait)-0.5*_m4t8bp6*_bc2nnuv*_7dknncp/_mf9pml9/_n2y8fv2*(1.0+_0lbt3g2/_n2y8fv2*_qzw3c1q+a*_kzrezrd/_gsq7ait)*2.0*_kzrezrd-_m4t8bp6*_bc2nnuv*_7dknncp/_in9dl5t/_nz3jbck*(1.0+_0lbt3g2/_n2y8fv2*_qzw3c1q+a*_kzrezrd/_gsq7ait)*_len606n+_m4t8bp6*_bc2nnuv*_7dknncp/_in9dl5t/_n2y8fv2*((_h5csk30*_kzrezrd/_in9dl5t+1.0)/_n2y8fv2*_qzw3c1q-_0lbt3g2/_nz3jbck*_qzw3c1q*_len606n-0.5*_0lbt3g2/_n2y8fv2*a/_mf9pml9*2.0*_kzrezrd+a/_gsq7ait-a*_kzrezrd/_mqx5ggo*2.0*_kzrezrd))/M_PI/_jz0ifrm);

    _jp6mkt7 = _1rsrla8/2.0*(0.25*(_yc3rwqm*_b7r3ocw*_w9pfb3h*_iwg0odw+_b7r3ocw/_m9j98cw*(_wgabtgd*_pr97y2l-_8nyfwdc/_m9j98cw*_9rr7s20)-_b7r3ocw*_h4a3yjt/_6c2cu7o*(_wgabtgd*_pr97y2l-_8nyfwdc/_m9j98cw*_9rr7s20)/_in9dl5t+_b7r3ocw*_m4t8bp6/_m9j98cw*(a/_mf9pml9*_m4t8bp6*_pr97y2l+_8nyfwdc/_6c2cu7o*_9rr7s20/_in9dl5t*_m4t8bp6+_m4t8bp6/_m9j98cw*a/_mf9pml9*_8nyfwdc)+_b7r3ocw*_h5csk30*_pr97y2l/_n2y8fv2*_qzw3c1q-_b7r3ocw*_h4a3yjt*_h5csk30*_pr97y2l/_nz3jbck*_qzw3c1q/_in9dl5t-_b7r3ocw*_h4a3yjt*_h5csk30*_pr97y2l/_n2y8fv2*a/_mf9pml9+a*_bc2nnuv*_pr97y2l/_mf9pml9-3.0*a*_h4a3yjt*_bc2nnuv*_pr97y2l/_ln1p0vu+_bc2nnuv/_in9dl5t/_m9j98cw*(-_b7r3ocw*_pr97y2l+_8nyfwdc/_m9j98cw*_cv3tlcx+a*_8nyfwdc/_gsq7ait)-_h4a3yjt*_bc2nnuv/_mf9pml9/_m9j98cw*(-_b7r3ocw*_pr97y2l+_8nyfwdc/_m9j98cw*_cv3tlcx+a*_8nyfwdc/_gsq7ait)-_h4a3yjt*_bc2nnuv/_gsq7ait/_6c2cu7o*(-_b7r3ocw*_pr97y2l+_8nyfwdc/_m9j98cw*_cv3tlcx+a*_8nyfwdc/_gsq7ait)+_m4t8bp6*_bc2nnuv/_in9dl5t/_m9j98cw*(-_8nyfwdc/_6c2cu7o*_cv3tlcx/_in9dl5t*_m4t8bp6-_m4t8bp6/_m9j98cw*a/_mf9pml9*_8nyfwdc-2.0*a*_8nyfwdc/_mqx5ggo*_m4t8bp6)+_bc2nnuv/_in9dl5t/_n2y8fv2*(_h5csk30/_n2y8fv2*(_0lbt3g2*(_b7r3ocw*_h5csk30-a/_in9dl5t)*_pr97y2l+_yt91qe3*_mofmmwi*_h5csk30)-a*_kzrezrd*_h5csk30*_pr97y2l/_gsq7ait)-_h4a3yjt*_bc2nnuv/_mf9pml9/_n2y8fv2*(_h5csk30/_n2y8fv2*(_0lbt3g2*(_b7r3ocw*_h5csk30-a/_in9dl5t)*_pr97y2l+_yt91qe3*_mofmmwi*_h5csk30)-a*_kzrezrd*_h5csk30*_pr97y2l/_gsq7ait)-_h4a3yjt*_bc2nnuv/_gsq7ait/_nz3jbck*(_h5csk30/_n2y8fv2*(_0lbt3g2*(_b7r3ocw*_h5csk30-a/_in9dl5t)*_pr97y2l+_yt91qe3*_mofmmwi*_h5csk30)-a*_kzrezrd*_h5csk30*_pr97y2l/_gsq7ait)+_m4t8bp6*_bc2nnuv/_in9dl5t/_n2y8fv2*(-_h5csk30/_nz3jbck*(_0lbt3g2*(_b7r3ocw*_h5csk30-a/_in9dl5t)*_pr97y2l+_yt91qe3*_mofmmwi*_h5csk30)/_in9dl5t*_m4t8bp6+_h5csk30/_n2y8fv2*(1.0/_in9dl5t*_h5csk30*_m4t8bp6*(_b7r3ocw*_h5csk30-a/_in9dl5t)*_pr97y2l+_0lbt3g2*a/_mf9pml9*_m4t8bp6*_pr97y2l+_yt91qe3/_in9dl5t*_7dknncp*_m4t8bp6*_h5csk30)+2.0*a*_kzrezrd*_h5csk30*_pr97y2l/_mqx5ggo*_m4t8bp6))/M_PI/_jz0ifrm)
        + _s8hmbr5/2.0*(0.25*(_b7r3ocw*((_yt91qe3*_iwg0odw+_qa78rap)/_in9dl5t*_m4t8bp6/_m9j98cw-(_yt91qe3*_iwg0odw+1.0)*_h5csk30/_in9dl5t*_m4t8bp6/_n2y8fv2)-_b7r3ocw/_6c2cu7o*(-_b7r3ocw*_8nyfwdc*_pr97y2l+_qa78rap*_kzrezrd-a+a*_8nyfwdc*_pr97y2l/_in9dl5t+_h17j8qo/_m9j98cw*_9rr7s20)/_in9dl5t*_m4t8bp6+_b7r3ocw/_m9j98cw*(-a*_8nyfwdc*_pr97y2l/_mf9pml9*_m4t8bp6-_h17j8qo/_6c2cu7o*_9rr7s20/_in9dl5t*_m4t8bp6-_h17j8qo/_m9j98cw*a/_mf9pml9*_m4t8bp6)+_b7r3ocw*_pr97y2l/_nz3jbck*(_yms1mwn*_h5csk30-a*_mofmmwi/_in9dl5t/_h5csk30)/_in9dl5t*_m4t8bp6-_b7r3ocw*_pr97y2l/_n2y8fv2*(-a/_gsq7ait*_7dknncp*_m4t8bp6/_h5csk30+a*_mofmmwi/_mf9pml9/_h5csk30*_m4t8bp6)+3.0*a*_m4t8bp6*_bc2nnuv*_pr97y2l/_ln1p0vu*_8nyfwdc-_bc2nnuv/_6c2cu7o*(2.0*_qa78rap+1.0/_in9dl5t*(_b7r3ocw*_8nyfwdc*_pr97y2l+a)-_h17j8qo/_in9dl5t/_m9j98cw*_cv3tlcx-a*_h17j8qo/_mf9pml9)/_in9dl5t*_m4t8bp6+_bc2nnuv/_m9j98cw*(-1.0/_mf9pml9*(_b7r3ocw*_8nyfwdc*_pr97y2l+a)*_m4t8bp6+_h17j8qo/_mf9pml9/_m9j98cw*_cv3tlcx*_m4t8bp6+_h17j8qo/_gsq7ait/_6c2cu7o*_cv3tlcx*_m4t8bp6+_h17j8qo/_mqx5ggo/_m9j98cw*a*_m4t8bp6+3.0*a*_h17j8qo/_ln1p0vu*_m4t8bp6)-_bc2nnuv*_pr97y2l/_nz3jbck*(-_h5csk30*_7dknncp+a*_8nyfwdc*_kzrezrd/_mf9pml9/_h5csk30+_mofmmwi/_in9dl5t*(_yt91qe3*_h5csk30-_0lbt3g2/_n2y8fv2*_y8n1uto))/_in9dl5t*_m4t8bp6+_bc2nnuv*_pr97y2l/_n2y8fv2*(-3.0*a*_8nyfwdc*_kzrezrd/_ln1p0vu/_h5csk30*_m4t8bp6+1.0/_gsq7ait*_7dknncp*_m4t8bp6*(_yt91qe3*_h5csk30-_0lbt3g2/_n2y8fv2*_y8n1uto)-_mofmmwi/_mf9pml9*(_yt91qe3*_h5csk30-_0lbt3g2/_n2y8fv2*_y8n1uto)*_m4t8bp6+_mofmmwi/_in9dl5t*(-1.0/_in9dl5t*_h5csk30*_m4t8bp6/_n2y8fv2*_y8n1uto+_0lbt3g2/_nz3jbck*_y8n1uto/_in9dl5t*_m4t8bp6+_0lbt3g2/_n2y8fv2*a/_mf9pml9/_h5csk30*_m4t8bp6)))/M_PI/_jz0ifrm)
        + _xqt4p50/2.0*(0.25*(_b7r3ocw*(1.0/_m9j98cw*(1.0+a/_in9dl5t)-_h4a3yjt/_6c2cu7o*(1.0+a/_in9dl5t)/_in9dl5t-_h4a3yjt/_m9j98cw*a/_mf9pml9-_h5csk30/_n2y8fv2*_qzw3c1q+_h4a3yjt*_h5csk30/_nz3jbck*_qzw3c1q/_in9dl5t+_h4a3yjt*_h5csk30/_n2y8fv2*a/_mf9pml9)-_bc2nnuv/_in9dl5t*_zv4u3bb+_h4a3yjt*_bc2nnuv/_mf9pml9*_zv4u3bb-_m4t8bp6*_bc2nnuv/_in9dl5t*(-2.0*a/_mqx5ggo*_m4t8bp6-1.0/_6c2cu7o/_in9dl5t*_m4t8bp6)+_bc2nnuv*_h5csk30/_in9dl5t/_n2y8fv2*(_0lbt3g2/_n2y8fv2*_qzw3c1q+a*_kzrezrd/_gsq7ait)-_h4a3yjt*_bc2nnuv*_h5csk30/_mf9pml9/_n2y8fv2*(_0lbt3g2/_n2y8fv2*_qzw3c1q+a*_kzrezrd/_gsq7ait)-_h4a3yjt*_bc2nnuv*_h5csk30/_gsq7ait/_nz3jbck*(_0lbt3g2/_n2y8fv2*_qzw3c1q+a*_kzrezrd/_gsq7ait)+_m4t8bp6*_bc2nnuv*_h5csk30/_in9dl5t/_n2y8fv2*(1.0/_in9dl5t*_h5csk30*_m4t8bp6/_n2y8fv2*_qzw3c1q-_0lbt3g2/_nz3jbck*_qzw3c1q/_in9dl5t*_m4t8bp6-_0lbt3g2/_n2y8fv2*a/_mf9pml9*_m4t8bp6-2.0*a*_kzrezrd/_mqx5ggo*_m4t8bp6))/M_PI/_jz0ifrm)
        + _1rsrla8/2.0*(0.25*(_b7r3ocw*((_yt91qe3*_iwg0odw-_qa78rap)/_in9dl5t*_8nyfwdc/_m9j98cw-(_yt91qe3*_iwg0odw+1.0-2.0*_qa78rap)*_h5csk30*_112hxrj/_n2y8fv2)+_b7r3ocw/_6c2cu7o*(_8nyfwdc*_pr97y2l*_wgabtgd+_qa78rap*_kzrezrd-a+_h4a3yjt/_m9j98cw*_9rr7s20)/_in9dl5t*_8nyfwdc-_b7r3ocw/_m9j98cw*(_wgabtgd*_pr97y2l+a*_h17j8qo*_pr97y2l/_mf9pml9-_h4a3yjt/_6c2cu7o*_9rr7s20/_in9dl5t*_8nyfwdc-_h4a3yjt/_m9j98cw*a/_mf9pml9*_8nyfwdc)-_b7r3ocw*_h5csk30*_pr97y2l/_n2y8fv2*_qzw3c1q+_b7r3ocw*_yms1mwn*_pr97y2l/_nz3jbck*_qzw3c1q*_112hxrj+_b7r3ocw*_yms1mwn*_pr97y2l/_n2y8fv2*a/_mf9pml9*_8nyfwdc-a*_bc2nnuv*_pr97y2l/_mf9pml9+3.0*a*_h17j8qo*_bc2nnuv*_pr97y2l/_ln1p0vu-_bc2nnuv/_6c2cu7o*(-2.0*_qa78rap+1.0/_in9dl5t*(_b7r3ocw*_8nyfwdc*_pr97y2l-a)+_h4a3yjt/_in9dl5t/_m9j98cw*_cv3tlcx+a*_h4a3yjt/_mf9pml9)/_in9dl5t*_8nyfwdc+_bc2nnuv/_m9j98cw*(-1.0/_mf9pml9*(_b7r3ocw*_8nyfwdc*_pr97y2l-a)*_8nyfwdc+1.0/_in9dl5t*_b7r3ocw*_pr97y2l-_h4a3yjt/_mf9pml9/_m9j98cw*_cv3tlcx*_8nyfwdc-_h4a3yjt/_gsq7ait/_6c2cu7o*_cv3tlcx*_8nyfwdc-_h4a3yjt/_mqx5ggo/_m9j98cw*a*_8nyfwdc-3.0*a*_h4a3yjt/_ln1p0vu*_8nyfwdc)-_bc2nnuv/_nz3jbck*(_88d0tgs-1.0/_in9dl5t*(_b7r3ocw*_yms1mwn*_pr97y2l+a*_h5csk30)+a*_kzrezrd*_yms1mwn*_pr97y2l/_mf9pml9-1.0/_in9dl5t/_n2y8fv2*(_h4a3yjt*_88d0tgs-a*_yms1mwn*_pr97y2l/_in9dl5t*_0lbt3g2))*_112hxrj+_bc2nnuv/_n2y8fv2*(1.0/_mf9pml9*(_b7r3ocw*_yms1mwn*_pr97y2l+a*_h5csk30)*_8nyfwdc-1.0/_in9dl5t*_b7r3ocw*_h5csk30*_pr97y2l+a*_kzrezrd*_h5csk30*_pr97y2l/_mf9pml9-3.0*a*_kzrezrd*_yms1mwn*_pr97y2l/_ln1p0vu*_8nyfwdc+1.0/_mf9pml9/_n2y8fv2*(_h4a3yjt*_88d0tgs-a*_yms1mwn*_pr97y2l/_in9dl5t*_0lbt3g2)*_8nyfwdc+1.0/_in9dl5t/_nz3jbck*(_h4a3yjt*_88d0tgs-a*_yms1mwn*_pr97y2l/_in9dl5t*_0lbt3g2)*_112hxrj-1.0/_in9dl5t/_n2y8fv2*(-a*_h5csk30*_pr97y2l/_in9dl5t*_0lbt3g2+a*_yms1mwn*_pr97y2l/_mf9pml9*_0lbt3g2*_8nyfwdc-a*_yms1mwn*_pr97y2l/_gsq7ait*_h5csk30*_8nyfwdc)))/M_PI/_jz0ifrm)
        + _s8hmbr5/2.0*(0.25*(_yt91qe3*_b7r3ocw*_pmlvak9*_iwg0odw-_b7r3ocw*_m4t8bp6/_6c2cu7o*(_lvpymb5*_pr97y2l+_8nyfwdc/_m9j98cw*_9rr7s20)/_in9dl5t*_8nyfwdc+_b7r3ocw*_m4t8bp6/_m9j98cw*(-a/_mf9pml9*_8nyfwdc*_pr97y2l+1.0/_m9j98cw*_9rr7s20-_h17j8qo/_6c2cu7o*_9rr7s20/_in9dl5t-_h17j8qo/_m9j98cw*a/_mf9pml9)+_b7r3ocw*_m4t8bp6*_pr97y2l/_nz3jbck*_y8n1uto*_112hxrj+_b7r3ocw*_m4t8bp6*_pr97y2l/_n2y8fv2*a/_mf9pml9/_h5csk30*_8nyfwdc+3.0*a*_m4t8bp6*_bc2nnuv*_pr97y2l/_ln1p0vu*_8nyfwdc-_m4t8bp6*_bc2nnuv/_mf9pml9/_m9j98cw*(_b7r3ocw*_pr97y2l-2.0*_qa78rap*_8nyfwdc/_m9j98cw-a*_8nyfwdc/_in9dl5t*(1.0/_in9dl5t+1.0/_m9j98cw))*_8nyfwdc-_m4t8bp6*_bc2nnuv/_gsq7ait/_6c2cu7o*(_b7r3ocw*_pr97y2l-2.0*_qa78rap*_8nyfwdc/_m9j98cw-a*_8nyfwdc/_in9dl5t*(1.0/_in9dl5t+1.0/_m9j98cw))*_8nyfwdc+_m4t8bp6*_bc2nnuv/_in9dl5t/_m9j98cw*(-2.0*_qa78rap/_m9j98cw+2.0*_qa78rap*_h17j8qo/_6c2cu7o/_in9dl5t-a/_in9dl5t*(1.0/_in9dl5t+1.0/_m9j98cw)+a*_h17j8qo/_mf9pml9*(1.0/_in9dl5t+1.0/_m9j98cw)-a*_8nyfwdc/_in9dl5t*(-1.0/_mf9pml9*_8nyfwdc-1.0/_6c2cu7o/_in9dl5t*_8nyfwdc))-_m4t8bp6*_bc2nnuv*_pr97y2l/_mf9pml9/_n2y8fv2*(_yc3rwqm*_h5csk30+_0lbt3g2/_n2y8fv2*_y8n1uto+a*_kzrezrd/_gsq7ait/_h5csk30)*_8nyfwdc-_m4t8bp6*_bc2nnuv*_pr97y2l/_in9dl5t/_nz3jbck*(_yc3rwqm*_h5csk30+_0lbt3g2/_n2y8fv2*_y8n1uto+a*_kzrezrd/_gsq7ait/_h5csk30)*_112hxrj+_m4t8bp6*_bc2nnuv*_pr97y2l/_in9dl5t/_n2y8fv2*(1.0/_in9dl5t*_h5csk30*_8nyfwdc/_n2y8fv2*_y8n1uto-_0lbt3g2/_nz3jbck*_y8n1uto*_112hxrj-_0lbt3g2/_n2y8fv2*a/_mf9pml9/_h5csk30*_8nyfwdc-2.0*a*_kzrezrd/_mqx5ggo/_h5csk30*_8nyfwdc))/M_PI/_jz0ifrm)
        + _xqt4p50/2.0*(0.25*(_b7r3ocw*(-_7dknncp*_112hxrj/_n2y8fv2-1.0/_m9j98cw*(1.0+a/_in9dl5t)+_h17j8qo/_6c2cu7o*(1.0+a/_in9dl5t)/_in9dl5t+_h17j8qo/_m9j98cw*a/_mf9pml9+_h5csk30/_n2y8fv2*_qzw3c1q-_yms1mwn/_nz3jbck*_qzw3c1q*_112hxrj-_yms1mwn/_n2y8fv2*a/_mf9pml9*_8nyfwdc)+_bc2nnuv/_in9dl5t*_zv4u3bb-_h17j8qo*_bc2nnuv/_mf9pml9*_zv4u3bb+_8nyfwdc*_bc2nnuv/_in9dl5t*(-2.0*a/_mqx5ggo*_8nyfwdc-1.0/_6c2cu7o/_in9dl5t*_8nyfwdc)+_bc2nnuv/_nz3jbck*(_7dknncp*(_h5csk30-a/_in9dl5t)+_yms1mwn/_in9dl5t*(1.0+a*_kzrezrd/_gsq7ait)-1.0/_in9dl5t/_n2y8fv2*(_h4a3yjt*_h5csk30*_7dknncp-a*_yms1mwn/_in9dl5t*_0lbt3g2))*_112hxrj-_bc2nnuv/_n2y8fv2*(_7dknncp*a/_mf9pml9*_8nyfwdc+_h5csk30/_in9dl5t*(1.0+a*_kzrezrd/_gsq7ait)-_yms1mwn/_mf9pml9*(1.0+a*_kzrezrd/_gsq7ait)*_8nyfwdc-2.0*_yms1mwn/_ln1p0vu*a*_kzrezrd*_8nyfwdc+1.0/_mf9pml9/_n2y8fv2*(_h4a3yjt*_h5csk30*_7dknncp-a*_yms1mwn/_in9dl5t*_0lbt3g2)*_8nyfwdc+1.0/_in9dl5t/_nz3jbck*(_h4a3yjt*_h5csk30*_7dknncp-a*_yms1mwn/_in9dl5t*_0lbt3g2)*_112hxrj-1.0/_in9dl5t/_n2y8fv2*(-a*_h5csk30/_in9dl5t*_0lbt3g2+a*_yms1mwn/_mf9pml9*_0lbt3g2*_8nyfwdc-a*_yms1mwn/_gsq7ait*_h5csk30*_8nyfwdc)))/M_PI/_jz0ifrm);

    _s1gd1t5 = _1rsrla8/2.0*(0.25*(_yc3rwqm*_b7r3ocw*_n14eb87*_iwg0odw-_b7r3ocw*_m4t8bp6/_6c2cu7o*(_wgabtgd*_pr97y2l-_8nyfwdc/_m9j98cw*_9rr7s20)*_yxyur4s+_b7r3ocw*_m4t8bp6/_m9j98cw*(0.5*a/_mf9pml9*2.0*_kzrezrd*_pr97y2l+_8nyfwdc/_6c2cu7o*_9rr7s20*_yxyur4s+0.5*_8nyfwdc/_m9j98cw*a/_mf9pml9*2.0*_kzrezrd)-_b7r3ocw*_m4t8bp6*_h5csk30*_pr97y2l/_nz3jbck*_qzw3c1q*_len606n-0.5*_b7r3ocw*_m4t8bp6*_h5csk30*_pr97y2l/_n2y8fv2*a/_mf9pml9*2.0*_kzrezrd+a/_mf9pml9*_m4t8bp6*_pr97y2l-1.5*a*_m4t8bp6*_bc2nnuv*_pr97y2l/_ln1p0vu*2.0*_kzrezrd+_m4t8bp6/_in9dl5t/_m9j98cw*(-_b7r3ocw*_pr97y2l+_8nyfwdc/_m9j98cw*_cv3tlcx+a*_8nyfwdc/_gsq7ait)-0.5*_m4t8bp6*_bc2nnuv/_mf9pml9/_m9j98cw*(-_b7r3ocw*_pr97y2l+_8nyfwdc/_m9j98cw*_cv3tlcx+a*_8nyfwdc/_gsq7ait)*2.0*_kzrezrd-_m4t8bp6*_bc2nnuv/_in9dl5t/_6c2cu7o*(-_b7r3ocw*_pr97y2l+_8nyfwdc/_m9j98cw*_cv3tlcx+a*_8nyfwdc/_gsq7ait)*_yxyur4s+_m4t8bp6*_bc2nnuv/_in9dl5t/_m9j98cw*(-_8nyfwdc/_6c2cu7o*_cv3tlcx*_yxyur4s-0.5*_8nyfwdc/_m9j98cw*a/_mf9pml9*2.0*_kzrezrd-a*_8nyfwdc/_mqx5ggo*2.0*_kzrezrd)+_m4t8bp6/_in9dl5t/_n2y8fv2*(_h5csk30/_n2y8fv2*(_0lbt3g2*(_b7r3ocw*_h5csk30-a/_in9dl5t)*_pr97y2l+_yt91qe3*_mofmmwi*_h5csk30)-a*_kzrezrd*_h5csk30*_pr97y2l/_gsq7ait)-0.5*_m4t8bp6*_bc2nnuv/_mf9pml9/_n2y8fv2*(_h5csk30/_n2y8fv2*(_0lbt3g2*(_b7r3ocw*_h5csk30-a/_in9dl5t)*_pr97y2l+_yt91qe3*_mofmmwi*_h5csk30)-a*_kzrezrd*_h5csk30*_pr97y2l/_gsq7ait)*2.0*_kzrezrd-_m4t8bp6*_bc2nnuv/_in9dl5t/_nz3jbck*(_h5csk30/_n2y8fv2*(_0lbt3g2*(_b7r3ocw*_h5csk30-a/_in9dl5t)*_pr97y2l+_yt91qe3*_mofmmwi*_h5csk30)-a*_kzrezrd*_h5csk30*_pr97y2l/_gsq7ait)*_len606n+_m4t8bp6*_bc2nnuv/_in9dl5t/_n2y8fv2*(-_h5csk30/_nz3jbck*(_0lbt3g2*(_b7r3ocw*_h5csk30-a/_in9dl5t)*_pr97y2l+_yt91qe3*_mofmmwi*_h5csk30)*_len606n+_h5csk30/_n2y8fv2*((_h5csk30*_kzrezrd/_in9dl5t+1.0)*(_b7r3ocw*_h5csk30-a/_in9dl5t)*_pr97y2l+0.5*_0lbt3g2*a/_mf9pml9*2.0*_kzrezrd*_pr97y2l+0.5*_yt91qe3/_in9dl5t*_7dknncp*2.0*_kzrezrd*_h5csk30)-a*_h5csk30*_pr97y2l/_gsq7ait+a*_kzrezrd*_h5csk30*_pr97y2l/_mqx5ggo*2.0*_kzrezrd))/M_PI/_jz0ifrm)
        + _s8hmbr5/2.0*(0.25*(_b7r3ocw*((_yt91qe3*_iwg0odw+_qa78rap)*_yxyur4s/_m9j98cw-(_yt91qe3*_iwg0odw+1.0)*_h5csk30*_len606n/_n2y8fv2)-_b7r3ocw/_6c2cu7o*(-_b7r3ocw*_8nyfwdc*_pr97y2l+_qa78rap*_kzrezrd-a+a*_8nyfwdc*_pr97y2l/_in9dl5t+_h17j8qo/_m9j98cw*_9rr7s20)*_yxyur4s+_b7r3ocw/_m9j98cw*(_qa78rap-0.5*a*_8nyfwdc*_pr97y2l/_mf9pml9*2.0*_kzrezrd-_h17j8qo/_6c2cu7o*_9rr7s20*_yxyur4s-0.5*_h17j8qo/_m9j98cw*a/_mf9pml9*2.0*_kzrezrd)+_b7r3ocw*_pr97y2l/_nz3jbck*(_yms1mwn*_h5csk30-a*_mofmmwi/_in9dl5t/_h5csk30)*_len606n-_b7r3ocw*_pr97y2l/_n2y8fv2*(_h5csk30*_7dknncp-0.5*a/_gsq7ait*_7dknncp*2.0*_kzrezrd/_h5csk30+0.5*a*_mofmmwi/_mf9pml9/_h5csk30*2.0*_kzrezrd)-a/_mf9pml9*_8nyfwdc*_pr97y2l+1.5*a*_8nyfwdc*_bc2nnuv*_pr97y2l/_ln1p0vu*2.0*_kzrezrd+1.0/_m9j98cw*(2.0*_qa78rap+1.0/_in9dl5t*(_b7r3ocw*_8nyfwdc*_pr97y2l+a)-_h17j8qo/_in9dl5t/_m9j98cw*_cv3tlcx-a*_h17j8qo/_mf9pml9)-_bc2nnuv/_6c2cu7o*(2.0*_qa78rap+1.0/_in9dl5t*(_b7r3ocw*_8nyfwdc*_pr97y2l+a)-_h17j8qo/_in9dl5t/_m9j98cw*_cv3tlcx-a*_h17j8qo/_mf9pml9)*_yxyur4s+_bc2nnuv/_m9j98cw*(-0.5/_mf9pml9*(_b7r3ocw*_8nyfwdc*_pr97y2l+a)*2.0*_kzrezrd+0.5*_h17j8qo/_mf9pml9/_m9j98cw*_cv3tlcx*2.0*_kzrezrd+_h17j8qo/_in9dl5t/_6c2cu7o*_cv3tlcx*_yxyur4s+0.5*_h17j8qo/_mqx5ggo/_m9j98cw*a*2.0*_kzrezrd+1.5*a*_h17j8qo/_ln1p0vu*2.0*_kzrezrd)+_pr97y2l/_n2y8fv2*(-_h5csk30*_7dknncp+a*_8nyfwdc*_kzrezrd/_mf9pml9/_h5csk30+_mofmmwi/_in9dl5t*(_yt91qe3*_h5csk30-_0lbt3g2/_n2y8fv2*_y8n1uto))-_bc2nnuv*_pr97y2l/_nz3jbck*(-_h5csk30*_7dknncp+a*_8nyfwdc*_kzrezrd/_mf9pml9/_h5csk30+_mofmmwi/_in9dl5t*(_yt91qe3*_h5csk30-_0lbt3g2/_n2y8fv2*_y8n1uto))*_len606n+_bc2nnuv*_pr97y2l/_n2y8fv2*(a/_mf9pml9/_h5csk30*_8nyfwdc-1.5*a*_8nyfwdc*_kzrezrd/_ln1p0vu/_h5csk30*2.0*_kzrezrd+0.5/_gsq7ait*_7dknncp*2.0*_kzrezrd*(_yt91qe3*_h5csk30-_0lbt3g2/_n2y8fv2*_y8n1uto)-0.5*_mofmmwi/_mf9pml9*(_yt91qe3*_h5csk30-_0lbt3g2/_n2y8fv2*_y8n1uto)*2.0*_kzrezrd+_mofmmwi/_in9dl5t*(-(_h5csk30*_kzrezrd/_in9dl5t+1.0)/_n2y8fv2*_y8n1uto+_0lbt3g2/_nz3jbck*_y8n1uto*_len606n+0.5*_0lbt3g2/_n2y8fv2*a/_mf9pml9/_h5csk30*2.0*_kzrezrd)))/M_PI/_jz0ifrm)
        + _xqt4p50/2.0*(0.25*(_b7r3ocw*(-_m4t8bp6/_6c2cu7o*(1.0+a/_in9dl5t)*_yxyur4s-0.5*_m4t8bp6/_m9j98cw*a/_mf9pml9*2.0*_kzrezrd+_m4t8bp6*_h5csk30/_nz3jbck*_qzw3c1q*_len606n+0.5*_m4t8bp6*_h5csk30/_n2y8fv2*a/_mf9pml9*2.0*_kzrezrd)-_m4t8bp6/_in9dl5t*_zv4u3bb+0.5*_m4t8bp6*_bc2nnuv/_mf9pml9*_zv4u3bb*2.0*_kzrezrd-_m4t8bp6*_bc2nnuv/_in9dl5t*(-a/_mqx5ggo*2.0*_kzrezrd-1.0/_6c2cu7o*_yxyur4s)+_m4t8bp6*_h5csk30/_in9dl5t/_n2y8fv2*(_0lbt3g2/_n2y8fv2*_qzw3c1q+a*_kzrezrd/_gsq7ait)-0.5*_m4t8bp6*_bc2nnuv*_h5csk30/_mf9pml9/_n2y8fv2*(_0lbt3g2/_n2y8fv2*_qzw3c1q+a*_kzrezrd/_gsq7ait)*2.0*_kzrezrd-_m4t8bp6*_bc2nnuv*_h5csk30/_in9dl5t/_nz3jbck*(_0lbt3g2/_n2y8fv2*_qzw3c1q+a*_kzrezrd/_gsq7ait)*_len606n+_m4t8bp6*_bc2nnuv*_h5csk30/_in9dl5t/_n2y8fv2*((_h5csk30*_kzrezrd/_in9dl5t+1.0)/_n2y8fv2*_qzw3c1q-_0lbt3g2/_nz3jbck*_qzw3c1q*_len606n-0.5*_0lbt3g2/_n2y8fv2*a/_mf9pml9*2.0*_kzrezrd+a/_gsq7ait-a*_kzrezrd/_mqx5ggo*2.0*_kzrezrd))/M_PI/_jz0ifrm)
        + _1rsrla8/2.0*(0.25*(_yt91qe3*(_b7r3ocw*_pmlvak9*_pr97y2l-_8nyfwdc/_6c2cu7o*_cv3tlcx/_in9dl5t*_m4t8bp6-_m4t8bp6/_m9j98cw*a/_mf9pml9*_8nyfwdc+_m4t8bp6*_h5csk30/_nz3jbck*_qzw3c1q*_112hxrj+_m4t8bp6*_h5csk30/_n2y8fv2*a/_mf9pml9*_8nyfwdc)-_m4t8bp6*_bc2nnuv/_mf9pml9*(2.0*_qa78rap/_m9j98cw+a/_gsq7ait)*_8nyfwdc+_m4t8bp6*_bc2nnuv/_in9dl5t*(-2.0*_qa78rap/_6c2cu7o/_in9dl5t*_8nyfwdc-2.0*a/_mqx5ggo*_8nyfwdc)-_m4t8bp6*_bc2nnuv*_h5csk30/_mf9pml9/_n2y8fv2*(1.0-2.0*_qa78rap-_0lbt3g2/_n2y8fv2*_qzw3c1q-a*_kzrezrd/_gsq7ait)*_8nyfwdc-_m4t8bp6*_bc2nnuv*_h5csk30/_in9dl5t/_nz3jbck*(1.0-2.0*_qa78rap-_0lbt3g2/_n2y8fv2*_qzw3c1q-a*_kzrezrd/_gsq7ait)*_112hxrj+_m4t8bp6*_bc2nnuv*_h5csk30/_in9dl5t/_n2y8fv2*(-1.0/_in9dl5t*_h5csk30*_8nyfwdc/_n2y8fv2*_qzw3c1q+_0lbt3g2/_nz3jbck*_qzw3c1q*_112hxrj+_0lbt3g2/_n2y8fv2*a/_mf9pml9*_8nyfwdc+2.0*a*_kzrezrd/_mqx5ggo*_8nyfwdc))/M_PI/_jz0ifrm)
        + _s8hmbr5/2.0*(0.25*(_yc3rwqm*_b7r3ocw*_pr97y2l*(1.0/_in9dl5t*_8nyfwdc/_m9j98cw-_h5csk30*_112hxrj/_n2y8fv2)-_yt91qe3/_m9j98cw*_cv3tlcx+_yt91qe3*_h17j8qo/_6c2cu7o*_cv3tlcx/_in9dl5t+_yt91qe3*_h17j8qo/_m9j98cw*a/_mf9pml9+_yt91qe3*_h5csk30/_n2y8fv2*_qzw3c1q-_yt91qe3*_yms1mwn/_nz3jbck*_qzw3c1q*_112hxrj-_yt91qe3*_yms1mwn/_n2y8fv2*a/_mf9pml9*_8nyfwdc-_bc2nnuv/_mf9pml9*(_b7r3ocw*_pr97y2l-2.0*_qa78rap*_8nyfwdc/_m9j98cw-a*_8nyfwdc/_gsq7ait)*_8nyfwdc+_bc2nnuv/_in9dl5t*(-2.0*_qa78rap/_m9j98cw+2.0*_qa78rap*_h17j8qo/_6c2cu7o/_in9dl5t-a/_gsq7ait+2.0*a*_h17j8qo/_mqx5ggo)+_bc2nnuv/_nz3jbck*(_h5csk30*_7dknncp+_0lbt3g2*_pr97y2l/_in9dl5t*(_yt91qe3*_h5csk30-_0lbt3g2/_n2y8fv2)+a/_in9dl5t*(_7dknncp-_kzrezrd*_yms1mwn/_gsq7ait-_yms1mwn*_0lbt3g2/_in9dl5t/_n2y8fv2))*_112hxrj-_bc2nnuv/_n2y8fv2*(1.0/_gsq7ait*_h5csk30*_8nyfwdc*_pr97y2l*(_yt91qe3*_h5csk30-_0lbt3g2/_n2y8fv2)-_0lbt3g2*_pr97y2l/_mf9pml9*(_yt91qe3*_h5csk30-_0lbt3g2/_n2y8fv2)*_8nyfwdc+_0lbt3g2*_pr97y2l/_in9dl5t*(-1.0/_in9dl5t*_h5csk30*_8nyfwdc/_n2y8fv2+_0lbt3g2/_nz3jbck*_112hxrj)-a/_mf9pml9*(_7dknncp-_kzrezrd*_yms1mwn/_gsq7ait-_yms1mwn*_0lbt3g2/_in9dl5t/_n2y8fv2)*_8nyfwdc+a/_in9dl5t*(-_kzrezrd*_h5csk30/_gsq7ait+2.0*_kzrezrd*_yms1mwn/_mqx5ggo*_8nyfwdc-_h5csk30*_0lbt3g2/_in9dl5t/_n2y8fv2-_yms1mwn/_gsq7ait*_h5csk30*_8nyfwdc/_n2y8fv2+_yms1mwn*_0lbt3g2/_mf9pml9/_n2y8fv2*_8nyfwdc+_yms1mwn*_0lbt3g2/_in9dl5t/_nz3jbck*_112hxrj)))/M_PI/_jz0ifrm)
        + _xqt4p50/2.0*(0.25*(_yt91qe3*_pmlvak9-_yt91qe3*_m4t8bp6*_7dknncp/_nz3jbck*_qzw3c1q*_112hxrj-_yt91qe3*_m4t8bp6*_7dknncp/_n2y8fv2*a/_mf9pml9*_8nyfwdc-_m4t8bp6*_bc2nnuv*_7dknncp/_mf9pml9/_n2y8fv2*(1.0+_0lbt3g2/_n2y8fv2*_qzw3c1q+a*_kzrezrd/_gsq7ait)*_8nyfwdc-_m4t8bp6*_bc2nnuv*_7dknncp/_in9dl5t/_nz3jbck*(1.0+_0lbt3g2/_n2y8fv2*_qzw3c1q+a*_kzrezrd/_gsq7ait)*_112hxrj+_m4t8bp6*_bc2nnuv*_7dknncp/_in9dl5t/_n2y8fv2*(1.0/_in9dl5t*_h5csk30*_8nyfwdc/_n2y8fv2*_qzw3c1q-_0lbt3g2/_nz3jbck*_qzw3c1q*_112hxrj-_0lbt3g2/_n2y8fv2*a/_mf9pml9*_8nyfwdc-2.0*a*_kzrezrd/_mqx5ggo*_8nyfwdc))/M_PI/_jz0ifrm);

    _stt1eqc = _1rsrla8/2.0*(0.25*(_b7r3ocw*((_yt91qe3*_iwg0odw-_qa78rap)*_yxyur4s/_m9j98cw-(_yt91qe3*_iwg0odw+1.0-2.0*_qa78rap)*_h5csk30*_len606n/_n2y8fv2)+_b7r3ocw/_6c2cu7o*(_8nyfwdc*_pr97y2l*_wgabtgd+_qa78rap*_kzrezrd-a+_h4a3yjt/_m9j98cw*_9rr7s20)*_yxyur4s-_b7r3ocw/_m9j98cw*(0.5*a*_8nyfwdc*_pr97y2l/_mf9pml9*2.0*_kzrezrd+_qa78rap-_h4a3yjt/_6c2cu7o*_9rr7s20*_yxyur4s-0.5*_h4a3yjt/_m9j98cw*a/_mf9pml9*2.0*_kzrezrd)-_b7r3ocw*_7dknncp*_pr97y2l/_n2y8fv2*_qzw3c1q+_b7r3ocw*_yms1mwn*_pr97y2l/_nz3jbck*_qzw3c1q*_len606n+0.5*_b7r3ocw*_yms1mwn*_pr97y2l/_n2y8fv2*a/_mf9pml9*2.0*_kzrezrd-a/_mf9pml9*_8nyfwdc*_pr97y2l+1.5*a*_8nyfwdc*_bc2nnuv*_pr97y2l/_ln1p0vu*2.0*_kzrezrd+1.0/_m9j98cw*(-2.0*_qa78rap+1.0/_in9dl5t*(_b7r3ocw*_8nyfwdc*_pr97y2l-a)+_h4a3yjt/_in9dl5t/_m9j98cw*_cv3tlcx+a*_h4a3yjt/_mf9pml9)-_bc2nnuv/_6c2cu7o*(-2.0*_qa78rap+1.0/_in9dl5t*(_b7r3ocw*_8nyfwdc*_pr97y2l-a)+_h4a3yjt/_in9dl5t/_m9j98cw*_cv3tlcx+a*_h4a3yjt/_mf9pml9)*_yxyur4s+_bc2nnuv/_m9j98cw*(-0.5/_mf9pml9*(_b7r3ocw*_8nyfwdc*_pr97y2l-a)*2.0*_kzrezrd-0.5*_h4a3yjt/_mf9pml9/_m9j98cw*_cv3tlcx*2.0*_kzrezrd-_h4a3yjt/_in9dl5t/_6c2cu7o*_cv3tlcx*_yxyur4s-0.5*_h4a3yjt/_mqx5ggo/_m9j98cw*a*2.0*_kzrezrd-1.5*a*_h4a3yjt/_ln1p0vu*2.0*_kzrezrd)+1.0/_n2y8fv2*(_88d0tgs-1.0/_in9dl5t*(_b7r3ocw*_yms1mwn*_pr97y2l+a*_h5csk30)+a*_kzrezrd*_yms1mwn*_pr97y2l/_mf9pml9-1.0/_in9dl5t/_n2y8fv2*(_h4a3yjt*_88d0tgs-a*_yms1mwn*_pr97y2l/_in9dl5t*_0lbt3g2))-_bc2nnuv/_nz3jbck*(_88d0tgs-1.0/_in9dl5t*(_b7r3ocw*_yms1mwn*_pr97y2l+a*_h5csk30)+a*_kzrezrd*_yms1mwn*_pr97y2l/_mf9pml9-1.0/_in9dl5t/_n2y8fv2*(_h4a3yjt*_88d0tgs-a*_yms1mwn*_pr97y2l/_in9dl5t*_0lbt3g2))*_len606n+_bc2nnuv/_n2y8fv2*(0.5/_mf9pml9*(_b7r3ocw*_yms1mwn*_pr97y2l+a*_h5csk30)*2.0*_kzrezrd-1.0/_in9dl5t*_b7r3ocw*_7dknncp*_pr97y2l+a*_yms1mwn*_pr97y2l/_mf9pml9+a*_kzrezrd*_7dknncp*_pr97y2l/_mf9pml9-1.5*a*_kzrezrd*_yms1mwn*_pr97y2l/_ln1p0vu*2.0*_kzrezrd+0.5/_mf9pml9/_n2y8fv2*(_h4a3yjt*_88d0tgs-a*_yms1mwn*_pr97y2l/_in9dl5t*_0lbt3g2)*2.0*_kzrezrd+1.0/_in9dl5t/_nz3jbck*(_h4a3yjt*_88d0tgs-a*_yms1mwn*_pr97y2l/_in9dl5t*_0lbt3g2)*_len606n-1.0/_in9dl5t/_n2y8fv2*(-a*_7dknncp*_pr97y2l/_in9dl5t*_0lbt3g2+0.5*a*_yms1mwn*_pr97y2l/_mf9pml9*_0lbt3g2*2.0*_kzrezrd-a*_yms1mwn*_pr97y2l/_in9dl5t*(_h5csk30*_kzrezrd/_in9dl5t+1.0))))/M_PI/_jz0ifrm)
        + _s8hmbr5/2.0*(0.25*(_yt91qe3*_b7r3ocw*_n14eb87*_iwg0odw-_b7r3ocw*_m4t8bp6/_6c2cu7o*(_lvpymb5*_pr97y2l+_8nyfwdc/_m9j98cw*_9rr7s20)*_yxyur4s+_b7r3ocw*_m4t8bp6/_m9j98cw*(-0.5*a/_mf9pml9*2.0*_kzrezrd*_pr97y2l-_8nyfwdc/_6c2cu7o*_9rr7s20*_yxyur4s-0.5*_8nyfwdc/_m9j98cw*a/_mf9pml9*2.0*_kzrezrd)+_b7r3ocw*_m4t8bp6*_pr97y2l/_nz3jbck*_y8n1uto*_len606n+0.5*_b7r3ocw*_m4t8bp6*_pr97y2l/_n2y8fv2*a/_mf9pml9/_h5csk30*2.0*_kzrezrd-a/_mf9pml9*_m4t8bp6*_pr97y2l+1.5*a*_m4t8bp6*_bc2nnuv*_pr97y2l/_ln1p0vu*2.0*_kzrezrd+_m4t8bp6/_in9dl5t/_m9j98cw*(_b7r3ocw*_pr97y2l-2.0*_qa78rap*_8nyfwdc/_m9j98cw-a*_8nyfwdc/_in9dl5t*(1.0/_in9dl5t+1.0/_m9j98cw))-0.5*_m4t8bp6*_bc2nnuv/_mf9pml9/_m9j98cw*(_b7r3ocw*_pr97y2l-2.0*_qa78rap*_8nyfwdc/_m9j98cw-a*_8nyfwdc/_in9dl5t*(1.0/_in9dl5t+1.0/_m9j98cw))*2.0*_kzrezrd-_m4t8bp6*_bc2nnuv/_in9dl5t/_6c2cu7o*(_b7r3ocw*_pr97y2l-2.0*_qa78rap*_8nyfwdc/_m9j98cw-a*_8nyfwdc/_in9dl5t*(1.0/_in9dl5t+1.0/_m9j98cw))*_yxyur4s+_m4t8bp6*_bc2nnuv/_in9dl5t/_m9j98cw*(2.0*_qa78rap*_8nyfwdc/_6c2cu7o*_yxyur4s+0.5*a*_8nyfwdc/_mf9pml9*(1.0/_in9dl5t+1.0/_m9j98cw)*2.0*_kzrezrd-a*_8nyfwdc/_in9dl5t*(-0.5/_mf9pml9*2.0*_kzrezrd-1.0/_6c2cu7o*_yxyur4s))+_m4t8bp6*_pr97y2l/_in9dl5t/_n2y8fv2*(_yc3rwqm*_h5csk30+_0lbt3g2/_n2y8fv2*_y8n1uto+a*_kzrezrd/_gsq7ait/_h5csk30)-0.5*_m4t8bp6*_bc2nnuv*_pr97y2l/_mf9pml9/_n2y8fv2*(_yc3rwqm*_h5csk30+_0lbt3g2/_n2y8fv2*_y8n1uto+a*_kzrezrd/_gsq7ait/_h5csk30)*2.0*_kzrezrd-_m4t8bp6*_bc2nnuv*_pr97y2l/_in9dl5t/_nz3jbck*(_yc3rwqm*_h5csk30+_0lbt3g2/_n2y8fv2*_y8n1uto+a*_kzrezrd/_gsq7ait/_h5csk30)*_len606n+_m4t8bp6*_bc2nnuv*_pr97y2l/_in9dl5t/_n2y8fv2*((_h5csk30*_kzrezrd/_in9dl5t+1.0)/_n2y8fv2*_y8n1uto-_0lbt3g2/_nz3jbck*_y8n1uto*_len606n-0.5*_0lbt3g2/_n2y8fv2*a/_mf9pml9/_h5csk30*2.0*_kzrezrd+a/_gsq7ait/_h5csk30-a*_kzrezrd/_mqx5ggo/_h5csk30*2.0*_kzrezrd))/M_PI/_jz0ifrm)
        + _xqt4p50/2.0*(0.25*(_b7r3ocw*(-_7dknncp*_len606n/_n2y8fv2+_8nyfwdc/_6c2cu7o*(1.0+a/_in9dl5t)*_yxyur4s+0.5*_8nyfwdc/_m9j98cw*a/_mf9pml9*2.0*_kzrezrd+_7dknncp/_n2y8fv2*_qzw3c1q-_yms1mwn/_nz3jbck*_qzw3c1q*_len606n-0.5*_yms1mwn/_n2y8fv2*a/_mf9pml9*2.0*_kzrezrd)+_8nyfwdc/_in9dl5t*_zv4u3bb-0.5*_8nyfwdc*_bc2nnuv/_mf9pml9*_zv4u3bb*2.0*_kzrezrd+_8nyfwdc*_bc2nnuv/_in9dl5t*(-a/_mqx5ggo*2.0*_kzrezrd-1.0/_6c2cu7o*_yxyur4s)-1.0/_n2y8fv2*(_7dknncp*(_h5csk30-a/_in9dl5t)+_yms1mwn/_in9dl5t*(1.0+a*_kzrezrd/_gsq7ait)-1.0/_in9dl5t/_n2y8fv2*(_h4a3yjt*_h5csk30*_7dknncp-a*_yms1mwn/_in9dl5t*_0lbt3g2))+_bc2nnuv/_nz3jbck*(_7dknncp*(_h5csk30-a/_in9dl5t)+_yms1mwn/_in9dl5t*(1.0+a*_kzrezrd/_gsq7ait)-1.0/_in9dl5t/_n2y8fv2*(_h4a3yjt*_h5csk30*_7dknncp-a*_yms1mwn/_in9dl5t*_0lbt3g2))*_len606n-_bc2nnuv/_n2y8fv2*(0.5*_7dknncp*a/_mf9pml9*2.0*_kzrezrd+_7dknncp/_in9dl5t*(1.0+a*_kzrezrd/_gsq7ait)-0.5*_yms1mwn/_mf9pml9*(1.0+a*_kzrezrd/_gsq7ait)*2.0*_kzrezrd+_yms1mwn/_in9dl5t*(a/_gsq7ait-a*_kzrezrd/_mqx5ggo*2.0*_kzrezrd)+0.5/_mf9pml9/_n2y8fv2*(_h4a3yjt*_h5csk30*_7dknncp-a*_yms1mwn/_in9dl5t*_0lbt3g2)*2.0*_kzrezrd+1.0/_in9dl5t/_nz3jbck*(_h4a3yjt*_h5csk30*_7dknncp-a*_yms1mwn/_in9dl5t*_0lbt3g2)*_len606n-1.0/_in9dl5t/_n2y8fv2*(-a*_7dknncp/_in9dl5t*_0lbt3g2+0.5*a*_yms1mwn/_mf9pml9*_0lbt3g2*2.0*_kzrezrd-a*_yms1mwn/_in9dl5t*(_h5csk30*_kzrezrd/_in9dl5t+1.0))))/M_PI/_jz0ifrm)
        + _1rsrla8/2.0*(0.25*(_yt91qe3*(_b7r3ocw*_w9pfb3h*_pr97y2l+1.0/_m9j98cw*_cv3tlcx-_h4a3yjt/_6c2cu7o*_cv3tlcx/_in9dl5t-_h4a3yjt/_m9j98cw*a/_mf9pml9-_h5csk30/_n2y8fv2*_qzw3c1q+_h4a3yjt*_h5csk30/_nz3jbck*_qzw3c1q/_in9dl5t+_h4a3yjt*_h5csk30/_n2y8fv2*a/_mf9pml9)+_bc2nnuv/_in9dl5t*(2.0*_qa78rap/_m9j98cw+a/_gsq7ait)-_h4a3yjt*_bc2nnuv/_mf9pml9*(2.0*_qa78rap/_m9j98cw+a/_gsq7ait)+_m4t8bp6*_bc2nnuv/_in9dl5t*(-2.0*_qa78rap/_6c2cu7o/_in9dl5t*_m4t8bp6-2.0*a/_mqx5ggo*_m4t8bp6)+_bc2nnuv*_h5csk30/_in9dl5t/_n2y8fv2*(1.0-2.0*_qa78rap-_0lbt3g2/_n2y8fv2*_qzw3c1q-a*_kzrezrd/_gsq7ait)-_h4a3yjt*_bc2nnuv*_h5csk30/_mf9pml9/_n2y8fv2*(1.0-2.0*_qa78rap-_0lbt3g2/_n2y8fv2*_qzw3c1q-a*_kzrezrd/_gsq7ait)-_h4a3yjt*_bc2nnuv*_h5csk30/_gsq7ait/_nz3jbck*(1.0-2.0*_qa78rap-_0lbt3g2/_n2y8fv2*_qzw3c1q-a*_kzrezrd/_gsq7ait)+_m4t8bp6*_bc2nnuv*_h5csk30/_in9dl5t/_n2y8fv2*(-1.0/_in9dl5t*_h5csk30*_m4t8bp6/_n2y8fv2*_qzw3c1q+_0lbt3g2/_nz3jbck*_qzw3c1q/_in9dl5t*_m4t8bp6+_0lbt3g2/_n2y8fv2*a/_mf9pml9*_m4t8bp6+2.0*a*_kzrezrd/_mqx5ggo*_m4t8bp6))/M_PI/_jz0ifrm)
        + _s8hmbr5/2.0*(0.25*(_yc3rwqm*_b7r3ocw*_pr97y2l*(1.0/_in9dl5t*_m4t8bp6/_m9j98cw-_h5csk30/_in9dl5t*_m4t8bp6/_n2y8fv2)+_yt91qe3*_8nyfwdc/_6c2cu7o*_cv3tlcx/_in9dl5t*_m4t8bp6+_yt91qe3*_8nyfwdc/_m9j98cw*a/_mf9pml9*_m4t8bp6-_yt91qe3*_yms1mwn/_nz3jbck*_qzw3c1q/_in9dl5t*_m4t8bp6-_yt91qe3*_yms1mwn/_n2y8fv2*a/_mf9pml9*_m4t8bp6-_bc2nnuv/_mf9pml9*(_b7r3ocw*_pr97y2l-2.0*_qa78rap*_8nyfwdc/_m9j98cw-a*_8nyfwdc/_gsq7ait)*_m4t8bp6+_bc2nnuv/_in9dl5t*(2.0*_qa78rap*_8nyfwdc/_6c2cu7o/_in9dl5t*_m4t8bp6+2.0*a*_8nyfwdc/_mqx5ggo*_m4t8bp6)+_bc2nnuv/_nz3jbck*(_h5csk30*_7dknncp+_0lbt3g2*_pr97y2l/_in9dl5t*(_yt91qe3*_h5csk30-_0lbt3g2/_n2y8fv2)+a/_in9dl5t*(_7dknncp-_kzrezrd*_yms1mwn/_gsq7ait-_yms1mwn*_0lbt3g2/_in9dl5t/_n2y8fv2))/_in9dl5t*_m4t8bp6-_bc2nnuv/_n2y8fv2*(1.0/_gsq7ait*_h5csk30*_m4t8bp6*_pr97y2l*(_yt91qe3*_h5csk30-_0lbt3g2/_n2y8fv2)-_0lbt3g2*_pr97y2l/_mf9pml9*(_yt91qe3*_h5csk30-_0lbt3g2/_n2y8fv2)*_m4t8bp6+_0lbt3g2*_pr97y2l/_in9dl5t*(-_h5csk30/_in9dl5t*_m4t8bp6/_n2y8fv2+_0lbt3g2/_nz3jbck/_in9dl5t*_m4t8bp6)-a/_mf9pml9*(_7dknncp-_kzrezrd*_yms1mwn/_gsq7ait-_yms1mwn*_0lbt3g2/_in9dl5t/_n2y8fv2)*_m4t8bp6+a/_in9dl5t*(2.0*_kzrezrd*_yms1mwn/_mqx5ggo*_m4t8bp6-_yms1mwn/_gsq7ait*_h5csk30*_m4t8bp6/_n2y8fv2+_yms1mwn*_0lbt3g2/_mf9pml9/_n2y8fv2*_m4t8bp6+_yms1mwn*_0lbt3g2/_gsq7ait/_nz3jbck*_m4t8bp6)))/M_PI/_jz0ifrm)
        + _xqt4p50/2.0*(0.25*(_yt91qe3*_w9pfb3h+_yt91qe3*_7dknncp/_n2y8fv2*_qzw3c1q-_yt91qe3*_h4a3yjt*_7dknncp/_nz3jbck*_qzw3c1q/_in9dl5t-_yt91qe3*_h4a3yjt*_7dknncp/_n2y8fv2*a/_mf9pml9+_bc2nnuv*_7dknncp/_in9dl5t/ _n2y8fv2*(1.0+_0lbt3g2/_n2y8fv2*_qzw3c1q+a*_kzrezrd/_gsq7ait)-_h4a3yjt*_bc2nnuv*_7dknncp/_mf9pml9/_n2y8fv2*(1.0+_0lbt3g2/_n2y8fv2*_qzw3c1q+a*_kzrezrd/_gsq7ait)-_h4a3yjt*_bc2nnuv*_7dknncp/_gsq7ait/_nz3jbck*(1.0+_0lbt3g2/_n2y8fv2*_qzw3c1q+a* _kzrezrd/_gsq7ait)+_m4t8bp6*_bc2nnuv*_7dknncp/_in9dl5t/_n2y8fv2*(1.0/_in9dl5t*_h5csk30*_m4t8bp6/_n2y8fv2*_qzw3c1q-_0lbt3g2/_nz3jbck* _qzw3c1q/_in9dl5t*_m4t8bp6-_0lbt3g2/_n2y8fv2*a/_mf9pml9*_m4t8bp6-2.0*a*_kzrezrd/_mqx5ggo*_m4t8bp6))/M_PI/_jz0ifrm);

    _5lf6a4f[0] = _4fcl5ef;           _5lf6a4f[1] = _jp6mkt7;         _5lf6a4f[2] = _s1gd1t5;
    _5lf6a4f[3] = _dfhg6w3;           _5lf6a4f[4] = _stt1eqc;         _5lf6a4f[5] = _37rln36;
    return;
}
