# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <float.h>





void       CoordTrans_inStrnFS(double newVal[3], double x_shift, double y_shift, double z_shift, double RotMat[3][3]); 

void      TriModeFind_inStrnFS(int TrimMode[1], double x,double y,double z,double p1_a,double p1_b, double p2_a,double p2_b,double p3_a,double p3_b);

void        TensTrans_inStrnFS(double e_out[6], double e_in[6], double B[3][3]);

void         TDSetupS_inStrnFS(double x,double y,double z,double alpha,double bx,double by,double bz,double nu, double TriVertex[3],double SideVec[3],double e[6]);

void     AngDisStrain_inStrnFS(double x, double y, double z, double alpha, double bx, double by, double bz, double nu, double e[6]);






void StrainFS_Nikkhoo(float Stress[6], float Strain[6], float fX, float fY, float fZ, float fP1[3], float fP2[3], float fP3[3], float fSs, float fDs, float fTs, const float fmu, const float flambda)
//
// this function is translated by Olaf Zielke from Matlab to C
// TDstressFS 
// Calculates stresses and strains associated with a triangular dislocation 
// in an elastic full-space.
//
// TD: Triangular Dislocation
// EFCS: Earth-Fixed Coordinate System
// TDCS: Triangular Dislocation Coordinate System
// ADCS: Angular Dislocation Coordinate System
// 
// INPUTS
// X, Y and Z: 
// Coordinates of calculation points in EFCS (East, North, Up). X, Y and Z 
// must have the same size.
//
// P1,P2 and P3:
// Coordinates of TD vertices in EFCS.
// 
// Ss, Ds and Ts:
// TD slip vector components (Strike-slip, Dip-slip, Tensile-slip).
//
// mu and lambda:
// Lame constants.
//
// OUTPUTS
// Stress:
// Calculated stress tensor components in EFCS. The six columns of Stress 
// are Sxx, Syy, Szz, Sxy, Sxz and Syz, respectively. The stress components 
// have the same unit as Lame constants.
//
// Strain:
// Calculated strain tensor components in EFCS. The six columns of Strain 
// are Exx, Eyy, Ezz, Exy, Exz and Eyz, respectively. The strain components 
// are dimensionless.
// 
// 
// Example: Calculate and plot the first component of stress tensor on a  
// regular grid.
// 
// [X,Y,Z] = meshgrid(-3:.02:3,-3:.02:3,2);
// [Stress,Strain] = TDstressFS(X,Y,Z,[-1 0 0],[1 -1 -1],[0 1.5 .5],...
// -1,2,3,.33e11,.33e11);
// h = surf(X,Y,reshape(Stress(:,1),size(X)),'edgecolor','none');
// view(2)
// axis equal
// axis tight
// set(gcf,'renderer','painters')

// Reference journal article: 
// Nikkhoo M. and Walter T.R., 2015. Triangular dislocation: An analytical, 
// artefact-free solution. 
// Submitted to Geophysical Journal International 

// Copyright (c) 2014 Mehdi Nikkhoo
// 
// Permission is hereby granted, free of charge, to any person obtaining a 
// copy of this software and associated documentation files 
// (the "Software"), to deal in the Software without restriction, including 
// without limitation the rights to use, copy, modify, merge, publish, 
// distribute, sublicense, and/or sell copies of the Software, and to permit
// persons to whom the Software is furnished to do so, subject to the 
// following conditions:
// 
// The above copyright notice and this permission notice shall be included 
// in all copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS 
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
// NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, 
// DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR 
// OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
// USE OR OTHER DEALINGS IN THE SOFTWARE.

// I appreciate any comments or bug reports.

// Mehdi Nikkhoo
// created: 2012.5.14
// Last modified: 2014.7.30
//
// VolcanoTectonics Research Group
// Section 2.1, Physics of Earthquakes and Volcanoes
// Department 2, Physics of the Earth
// Helmholtz Centre Potsdam
// German Research Centre for Geosciences (GFZ)
// 
// Email: 
// mehdi.nikkhoo@gfz-potsdam.de 
// mehdi.nikkhoo@gmail.com

{
   
    double X  = (double)fX;         double Y = (double)fY;          double Z = (double)fZ; 
    double bx = (double)fTs;        double by= (double)fSs;         double bz= (double)fDs;  /* bx = Ts; % Tensile-slip;      by = Ss; % Strike-slip;         bz = Ds; % Dip-slip */
    double mu = (double)fmu;        double lambda= (double)flambda; 
    double P1[3],                   P2[3],                          P3[3];
    P1[0] = (double)fP1[0];         P1[1] = (double)fP1[1];         P1[2] = (double)fP1[2];
    P2[0] = (double)fP2[0];         P2[1] = (double)fP2[1];         P2[2] = (double)fP2[2];
    P3[0] = (double)fP3[0];         P3[1] = (double)fP3[1];         P3[2] = (double)fP3[2];
    
    int      i, casepLog,          casenLog;
    double   A,      B,      C,  x,  y,      z,      nu,  Tempdouble;
    
    int      TriMode[1];
    double   Vstrike[3],        Vdip[3],            Vnorm[3];
    double   tempVect1[3],      tempVect2[3];
    double   p1[3],             p2[3],              p3[3];
    double   e12[3],            e13[3],             e23[3];
    double   eZ[3],             e_comb[6];
    double   e_vals1[6],        e_vals2[6],         e_vals3[6];
    double   Amat[3][3];

    nu = lambda/(2.0*(lambda+mu)); /* Poisson's ratio */
    p1[0]        = 0.0;                     p1[1]        = 0.0;                     p1[2]        = 0.0;
    p2[0]        = 0.0;                     p2[1]        = 0.0;                     p2[2]        = 0.0;
    p3[0]        = 0.0;                     p3[1]        = 0.0;                     p3[2]        = 0.0;
    eZ[0]        = 0.0;                     eZ[1]        = 0.0;                     eZ[2]        = 1.0;
    // Calculate unit strike, dip and normal to TD vectors: For a horizontal TD as an exception, if the normal vector points upward, the strike and dip 
    // vectors point Northward and Westward, whereas if the normal vector points downward, the strike and dip vectors point Southward and Westward, respectively.

    tempVect1[0] = P2[0] -P1[0];            tempVect1[1] = P2[1] -P1[1];            tempVect1[2] = P2[2] -P1[2];
    tempVect2[0] = P3[0] -P1[0];            tempVect2[1] = P3[1] -P1[1];            tempVect2[2] = P3[2] -P1[2];
      
    Vnorm[0]     = tempVect1[1]*tempVect2[2] - tempVect1[2]*tempVect2[1];
    Vnorm[1]     = tempVect1[2]*tempVect2[0] - tempVect1[0]*tempVect2[2];
    Vnorm[2]     = tempVect1[0]*tempVect2[1] - tempVect1[1]*tempVect2[0];
    Tempdouble   = sqrt(  Vnorm[0]*Vnorm[0] +Vnorm[1]*Vnorm[1] + Vnorm[2]*Vnorm[2]);
    
    Vnorm[0]     = Vnorm[0]/Tempdouble;     Vnorm[1]     = Vnorm[1]/Tempdouble;     Vnorm[2]     = Vnorm[2]/Tempdouble;
   
    Vstrike[0]   = eZ[1]*Vnorm[2] - eZ[2]*Vnorm[1];
    Vstrike[1]   = eZ[2]*Vnorm[0] - eZ[0]*Vnorm[2];
    Vstrike[2]   = eZ[0]*Vnorm[1] - eZ[1]*Vnorm[0];
    // For horizontal elements ("Vnorm(3)" adjusts for Northward or Southward direction)
    Tempdouble    = sqrt(  Vstrike[0]*Vstrike[0] +Vstrike[1]*Vstrike[1] + Vstrike[2]*Vstrike[2]);
    //if (Tempdouble < DBL_EPSILON)
    if (Tempdouble < FLT_EPSILON)
    {   Vstrike[0] = 0.0 ;                  Vstrike[1] = 1.0;            Vstrike[2] = 0.0;    /* For horizontal elements in case of half-space calculation!!! => Correct the strike vector of image dislocation only */
    }
    else
    {   Vstrike[0]/= Tempdouble;            Vstrike[1]/= Tempdouble;     Vstrike[2]/= Tempdouble;
    }
    Vdip[0]    = Vnorm[1]*Vstrike[2] - Vnorm[2]*Vstrike[1];
    Vdip[1]    = Vnorm[2]*Vstrike[0] - Vnorm[0]*Vstrike[2];
    Vdip[2]    = Vnorm[0]*Vstrike[1] - Vnorm[1]*Vstrike[0];
    Tempdouble = sqrt(  Vdip[0]*Vdip[0] +Vdip[1]*Vdip[1] + Vdip[2]*Vdip[2]);

    Vdip[0]    = Vdip[0]/Tempdouble;     Vdip[1] = Vdip[1]/Tempdouble;           Vdip[2]      = Vdip[2]/Tempdouble;
    /* Transform coordinates and slip vector components from EFCS into TDCS */
    Amat[0][0] = Vnorm[0];               Amat[0][1] = Vnorm[1];                  Amat[0][2]   = Vnorm[2];
    Amat[1][0] = Vstrike[0];             Amat[1][1] = Vstrike[1];                Amat[1][2]   = Vstrike[2];
    Amat[2][0] = Vdip[0];                Amat[2][1] = Vdip[1];                   Amat[2][2]   = Vdip[2];

    CoordTrans_inStrnFS(tempVect1, (X-P2[0]),  (Y-P2[1]), (Z-P2[2]), Amat); 
    x          = tempVect1[0];                  y          = tempVect1[1];                  z          = tempVect1[2];
    CoordTrans_inStrnFS(tempVect1, (P1[0]-P2[0]),  (P1[1]-P2[1]), (P1[2]-P2[2]), Amat); 
    p1[0]      = tempVect1[0];                   p1[1]      = tempVect1[1];                  p1[2]      = tempVect1[2];
    CoordTrans_inStrnFS(tempVect2, (P3[0]-P2[0]),  (P3[1]-P2[1]), (P3[2]-P2[2]), Amat); 
    p3[0]      = tempVect2[0];                   p3[1]      = tempVect2[1];                  p3[2]      = tempVect2[2];
    // Calculate the unit vectors along TD sides in TDCS
    Tempdouble = sqrt(  (p2[0]-p1[0])*(p2[0]-p1[0]) +(p2[1]-p1[1])*(p2[1]-p1[1]) + (p2[2]-p1[2])*(p2[2]-p1[2]));
    
    e12[0]     = (p2[0]-p1[0])/Tempdouble;       e12[1]      = (p2[1]-p1[1])/Tempdouble;       e12[2]     = (p2[2]-p1[2])/Tempdouble;     
    Tempdouble = sqrt(  (p3[0]-p1[0])*(p3[0]-p1[0]) +(p3[1]-p1[1])*(p3[1]-p1[1]) + (p3[2]-p1[2])*(p3[2]-p1[2]));
    e13[0]     = (p3[0]-p1[0])/Tempdouble;       e13[1]      = (p3[1]-p1[1])/Tempdouble;       e13[2]     = (p3[2]-p1[2])/Tempdouble;
    Tempdouble = sqrt(  (p3[0]-p2[0])*(p3[0]-p2[0]) +(p3[1]-p2[1])*(p3[1]-p2[1]) + (p3[2]-p2[2])*(p3[2]-p2[2]));
    e23[0]     = (p3[0]-p2[0])/Tempdouble;       e23[1]      = (p3[1]-p2[1])/Tempdouble;       e23[2]     = (p3[2]-p2[2])/Tempdouble; 
    // Calculate the TD angles
    Tempdouble  = e12[0]*e13[0] +e12[1]*e13[1] +e12[2]*e13[2];
    A = acos(Tempdouble) ;
    Tempdouble  = -1.0*e12[0]*e23[0] + -1.0*e12[1]*e23[1] + -1.0*e12[2]*e23[2];
    B = acos(Tempdouble) ;
    Tempdouble  = e23[0]*e13[0] +e23[1]*e13[1] +e23[2]*e13[2];
    C = acos(Tempdouble) ;
     
    // Determine the best arteact-free configuration for each calculation point

    TriModeFind_inStrnFS(TriMode,y,z,x,p1[1],p1[2], p2[1], p2[2], p3[1], p3[2]);

    if (TriMode[0] == 1)       {       casepLog = 1;   casenLog = 0;          } 
    if (TriMode[0] ==-1)       {       casepLog = 0;   casenLog = 1;          } 
    if (TriMode[0] == 0)       {       casepLog = 0;   casenLog = 0;          } 

    if (casepLog == 1) // Configuration I
    {   // Calculate first angular dislocation contribution
        tempVect1[0] = -1.0*e13[0];         tempVect1[1] = -1.0*e13[1];         tempVect1[2] = -1.0*e13[2];              
        TDSetupS_inStrnFS(x, y, z, A, bx, by, bz, nu, p1, tempVect1, e_vals1);
        // Calculate second angular dislocation contribution
        TDSetupS_inStrnFS(x, y, z, B, bx, by, bz, nu, p2, e12,       e_vals2);
        // Calculate third angular dislocation contribution
        TDSetupS_inStrnFS(x, y, z, C, bx, by, bz, nu, p3, e23,       e_vals3);  
    }
    if (casenLog == 1) // Configuration II
    {  // Calculate first angular dislocation contribution             
        TDSetupS_inStrnFS(x, y, z, A, bx, by, bz, nu, p1, e13,       e_vals1);
        // Calculate second angular dislocation contribution
        tempVect1[0] = -1.0*e12[0];         tempVect1[1] = -1.0*e12[1];         tempVect1[2] = -1.0*e12[2];      
        TDSetupS_inStrnFS(x, y, z, B, bx, by, bz, nu, p2, tempVect1, e_vals2);
        // Calculate third angular dislocation contribution
        tempVect1[0] = -1.0*e23[0];         tempVect1[1] = -1.0*e23[1];         tempVect1[2] = -1.0*e23[2]; 
        TDSetupS_inStrnFS(x, y, z, C, bx, by, bz, nu, p3, tempVect1, e_vals3);     
    }
    if ((casenLog == 1) || (casepLog == 1))
    {
        e_comb[0]       = e_vals1[0]+e_vals2[0]+e_vals3[0]; // exx
        e_comb[1]       = e_vals1[1]+e_vals2[1]+e_vals3[1]; // exy
        e_comb[2]       = e_vals1[2]+e_vals2[2]+e_vals3[2]; // exz
        e_comb[3]       = e_vals1[3]+e_vals2[3]+e_vals3[3]; // eyy
        e_comb[4]       = e_vals1[4]+e_vals2[4]+e_vals3[4]; // eyz
        e_comb[5]       = e_vals1[5]+e_vals2[5]+e_vals3[5]; // ezz
    }
    else
    {   
        e_comb[0]       = NAN; // exx => supposed to be "NaN" => have to check
        e_comb[1]       = NAN; // exy
        e_comb[2]       = NAN; // exz
        e_comb[3]       = NAN; // eyy
        e_comb[4]       = NAN; // eyz
        e_comb[5]       = NAN; // ezz
    }

    Amat[0][0] = Vnorm[0];         Amat[0][1] = Vstrike[0];       Amat[0][2] = Vdip[0];
    Amat[1][0] = Vnorm[1];         Amat[1][1] = Vstrike[1];       Amat[1][2] = Vdip[1];
    Amat[2][0] = Vnorm[2];         Amat[2][1] = Vstrike[2];       Amat[2][2] = Vdip[2];
    // Transform the strain tensor components from TDCS into EFCS

    double dStrain[6];
    TensTrans_inStrnFS(dStrain, e_comb, Amat);
    for (i = 0; i < 6; i++) {   Strain[i] = (float)dStrain[i];     }

    // Calculate the stress tensor components in EFCS
    Stress[0] = 2.0*mu* Strain[0]+lambda*(Strain[0]+Strain[3]+Strain[5]);       // sxx
    Stress[3] = 2.0*mu* Strain[3]+lambda*(Strain[0]+Strain[3]+Strain[5]);       // syy
    Stress[5] = 2.0*mu* Strain[5]+lambda*(Strain[0]+Strain[3]+Strain[5]);       // szz
    Stress[1] = 2.0*mu* Strain[1];                                           // sxy
    Stress[2] = 2.0*mu* Strain[2];                                           // sxz
    Stress[4] = 2.0*mu* Strain[4];                                           // syz
    
    return;
}
//---------------------------------------------------------------------------------
//---------------------------------------------------------------------------------
//---------------------------------------------------------------------------------
//---------------------------------------------------------------------------------
//---------------------------------------------------------------------------------
//---------------------------------------------------------------------------------
void  CoordTrans_inStrnFS(double newVal[3], double x1, double x2, double x3, double A[3][3]) 
{
    // CoordTrans transforms the coordinates of the vectors, from x1x2x3 coordinate system to X1X2X3 coordinate system. "A" is the
    // transformation matrix, whose columns e1,e2 and e3 are the unit base vectors of the x1x2x3. The coordinates of e1,e2 and e3 in A must be given 
    // in X1X2X3. The transpose of A (i.e., A') will transform the coordinates  from X1X2X3 into x1x2x3.
    
    newVal[0] = A[0][0]*x1 + A[0][1]*x2 + A[0][2]*x3;     newVal[1] = A[1][0]*x1 + A[1][1]*x2 + A[1][2]*x3;       newVal[2] = A[2][0]*x1 + A[2][1]*x2 + A[2][2]*x3;

    return;
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
void TensTrans_inStrnFS(double e_out[6], double e_in[6], double B[3][3])
{
    // TensTrans Transforms the coordinates of tensors,from x1y1z1 coordinate system to x2y2z2 coordinate system. "A" is the transformation matrix, 
    // whose columns e1,e2 and e3 are the unit base vectors of the x1y1z1. The coordinates of e1,e2 and e3 in A must be given in x2y2z2. The transpose 
    // of A (i.e., A') does the transformation from x2y2z2 into x1y1z1.
    double Txx1,         Txy1,           Txz1,           Tyy1;
    double Tyz1,         Tzz1,           Txx2,           Txy2;
    double Txz2,         Tyy2,           Tyz2,           Tzz2;
    double A[9];
    //---------------------------------
    Txx1 = e_in[0];         Txy1 = e_in[1];         Txz1 = e_in[2];
    Tyy1 = e_in[3];         Tyz1 = e_in[4];         Tzz1 = e_in[5];

    A[0] = B[0][0];         A[1] = B[1][0];         A[2] = B[2][0];
    A[3] = B[0][1];         A[4] = B[1][1];         A[5] = B[2][1];
    A[6] = B[0][2];         A[7] = B[1][2];         A[8] = B[2][2];

    Txx2 = A[0]*A[0]*Txx1 +         2.0*A[0]*A[3] *Txy1 +          2.0*A[0]*A[6] *Txz1 +          2.0*A[3]*A[6] *Tyz1 + A[3]*A[3]*Tyy1 + A[6]*A[6]*Tzz1;
    Tyy2 = A[1]*A[1]*Txx1 +         2.0*A[1]*A[4] *Txy1 +          2.0*A[1]*A[7] *Txz1 +          2.0*A[4]*A[7] *Tyz1 + A[4]*A[4]*Tyy1 + A[7]*A[7]*Tzz1;
    Tzz2 = A[2]*A[2]*Txx1 +         2.0*A[2]*A[5] *Txy1 +          2.0*A[2]*A[8] *Txz1 +          2.0*A[5]*A[8] *Tyz1 + A[5]*A[5]*Tyy1 + A[8]*A[8]*Tzz1;
    Txy2 = A[0]*A[1]*Txx1 + (A[0]*A[4]+ A[1]*A[3])*Txy1 + (A[0]*A[7] + A[1]*A[6])*Txz1 + (A[7]*A[3] + A[6]*A[4])*Tyz1 + A[4]*A[3]*Tyy1 + A[6]*A[7]*Tzz1;
    Txz2 = A[0]*A[2]*Txx1 + (A[0]*A[5]+ A[2]*A[3])*Txy1 + (A[0]*A[8] + A[2]*A[6])*Txz1 + (A[8]*A[3] + A[6]*A[5])*Tyz1 + A[5]*A[3]*Tyy1 + A[6]*A[8]*Tzz1;
    Tyz2 = A[1]*A[2]*Txx1 + (A[2]*A[4]+ A[1]*A[5])*Txy1 + (A[2]*A[7] + A[1]*A[8])*Txz1 + (A[7]*A[5] + A[8]*A[4])*Tyz1 + A[4]*A[5]*Tyy1 + A[7]*A[8]*Tzz1;
    //---------------------------------
    e_out[0] = Txx2;        e_out[1] = Txy2;      e_out[2] = Txz2;
    e_out[3] = Tyy2;        e_out[4] = Tyz2;      e_out[5] = Tzz2;
    
    return;
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
void      TriModeFind_inStrnFS(int TriMode[1], double x,double y,double z,double p1_a,double p1_b, double p2_a,double p2_b,double p3_a,double p3_b)
{
    // trimodefinder calculates the normalized barycentric coordinates of  the points with respect to the TD vertices and specifies the appropriate
    // artefact-free configuration of the angular dislocations for the  calculations. The input matrices x, y and z share the same size and
    // correspond to the y, z and x coordinates in the TDCS, respectively. p1, p2 and p3 are two-component matrices representing the y and z coordinates
    // of the TD vertices in the TDCS, respectively. The components of the output (trimode) corresponding to each calculation 
    // points, are 1 for the first configuration, -1 for the second configuration and 0 for the calculation point that lie on the TD sides.
    double a,            b,          c;

    a = ((p2_b-p3_b)*(x-p3_a) +(p3_a-p2_a)*(y-p3_b)) /  ((p2_b-p3_b)*(p1_a-p3_a) +(p3_a-p2_a)*(p1_b-p3_b));
    b = ((p3_b-p1_b)*(x-p3_a) +(p1_a-p3_a)*(y-p3_b)) /  ((p2_b-p3_b)*(p1_a-p3_a) +(p3_a-p2_a)*(p1_b-p3_b));
    c = 1.0 -a -b;

    TriMode[0] = 1;
    if       ((a <= 0.0) && (b > c) && (c > a))             {   TriMode[0] = -1;                  }
    else if  ((b <= 0.0) && (c > a) && (a > b))             {   TriMode[0] = -1;                  }
    else if  ((c <= 0.0) && (a > b) && (b > c))             {   TriMode[0] = -1;                  }

    else if  ((a == 0.0) && (b >= 0.0) && (c >= 0.0))       {   TriMode[0] = 0;                   }
    else if  ((a >= 0.0) && (b == 0.0) && (c >= 0.0))       {   TriMode[0] = 0;                   }
    else if  ((a >= 0.0) && (b >= 0.0) && (c == 0.0))       {   TriMode[0] = 0;                   }
    if       ((TriMode[0] == 0) && (z != 0.0))              {   TriMode[0] = 1;                   } 

    return;
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
void TDSetupS_inStrnFS(double x,double y,double z,double alpha,double bx,double by,double bz,double nu, double TriVertex[3],double SideVec[3],double e_out[6])
{   // TDSetupS transforms coordinates of the calculation points as well as  slip vector components from ADCS into TDCS. 
    // It then calculates the  strains in ADCS and transforms them into TDCS.
    double A[2][2];
    double B[3][3];
    double r1[2];
    double r2[2];
    double y1;
    double z1;
    double by1;
    double bz1;
    double tempVect1[2];
    double e[6];
    //-----------------------------------
    // Transformation matrix
    A[0][0] = SideVec[2];                   A[0][1]      = -1.0*SideVec[1];
    A[1][0] = SideVec[1];                   A[1][1]      =      SideVec[2]; 
    // Transform coordinates of the calculation points from TDCS into ADCS
    tempVect1[0] = y -TriVertex[1];         tempVect1[1] = z -TriVertex[2];
    r1[0]        = A[0][0]*tempVect1[0] +A[0][1]*tempVect1[1];
    r1[1]        = A[1][0]*tempVect1[0] +A[1][1]*tempVect1[1];
    y1           = r1[0];                   z1           = r1[1];
    // Transform the in-plane slip vector components from TDCS into ADCS
    tempVect1[0] = by;                      tempVect1[1] = bz;
    r2[0]        = A[0][0]*tempVect1[0] +A[0][1]*tempVect1[1];
    r2[1]        = A[1][0]*tempVect1[0] +A[1][1]*tempVect1[1];
    by1          = r2[0];                   bz1          = r2[1];
    
    // Calculate strains associated with an angular dislocation in ADCS
    AngDisStrain_inStrnFS(x, y1, z1, (-1.0*M_PI+alpha), bx, by1, bz1, nu, e); 
    // Transform strains from ADCS into TDCS
    B[0][0] = 1.0;          B[0][1] = 0.0;      B[0][2] = 0.0;
    B[1][0] = 0.0;          B[1][1] = A[0][0];  B[1][2] = A[1][0];
    B[2][0] = 0.0;          B[2][1] = A[0][1];  B[2][2] = A[1][1];// 3x3 Transformation matrix
    
    TensTrans_inStrnFS(e_out, e, B); //the e_out is then send back from the function




    return;
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
void AngDisStrain_inStrnFS(double x, double y, double z, double alpha, double bx, double by, double bz, double nu, double e[6])
// AngDisStrain calculates the strains associated with an angular dislocation in an elastic full-space.
{   double       sinA,           cosA,           eta,            zeta;
    double       x2,             y2,             z2,             r2;
    double       r,              r3,             rz,             r2z2;
    double       r3z,            W,              W2,             Wr;
    double       W2r,            Wr3,            W2r2,           C;
    double       S,              rFi_rx,         rFi_ry,         rFi_rz;
    double       Exx,            Exy,            Exz,            Eyy;
    double       Eyz,            Ezz;  
    //----------------------------------------------------------   

    sinA = sin(alpha);                          cosA = cos(alpha);
    eta = y*cosA - z*sinA;                      zeta = y*sinA + z*cosA;

    x2 = x*x;               y2   = y*y;                 z2  = z*z;
    r2 = x2 + y2 + z2;      r    = pow(r2,0.5);         r3  = r*r2;
    rz = r*(r-z);           r2z2 = r2*pow((r-z),2.0);   r3z = r3*(r-z);

    W  = zeta-r;            W2   = W*W;                 Wr  = W*r;
    W2r= W2*r;              Wr3  = W*r3;                W2r2= W2*r2;

    C = (r*cosA-z)/Wr;      S    = (r*sinA-y)/Wr;
    // Partial derivatives of the Burgers' function
    rFi_rx = (eta/r/(r-zeta)  -y/r/(r-z))   /4.0/M_PI;
    rFi_ry = (  x/r/(r-z)-cosA*x/r/(r-zeta))/4.0/M_PI;
    rFi_rz = (sinA*x/r/(r-zeta))/4.0/M_PI;
    //----------------------------------------------------------
    Exx = bx*(rFi_rx)  +bx/8.0/M_PI/(1.0-nu)*(eta/Wr+eta*x2/W2r2-eta*x2/Wr3+y/rz-x2*y/r2z2-x2*y/r3z)
                       -by*x/8.0/M_PI/(1.0-nu)*(((2.0*nu+1.0)/Wr+x2/W2r2-x2/Wr3)*cosA+(2.0*nu+1.0)/rz-x2/r2z2-x2/r3z)
                       +bz*x*sinA/8.0/M_PI/(1.0-nu)*((2.0*nu+1.0)/Wr+x2/W2r2-x2/Wr3);

    Eyy = by*(rFi_ry)  +bx/8.0/M_PI/(1.0-nu)*((1.0/Wr+S*S-y2/Wr3)*eta+(2.0*nu+1.0)*y/rz-pow(y,3.0)/r2z2-pow(y,3.0)/r3z-2.0*nu*cosA*S)
                       -by*x/8.0/M_PI/(1.0-nu)*(1.0/rz-y2/r2z2-y2/r3z+(1.0/Wr+S*S-y2/Wr3)*cosA)
                       +bz*x*sinA/8.0/M_PI/(1.0-nu)*(1.0/Wr+S*S-y2/Wr3);

    Ezz = bz*(rFi_rz)+bx/8.0/M_PI/(1.0-nu)*(eta/W/r+eta*C*C-eta*z2/Wr3+y*z/r3+2.0*nu*sinA*C)
                      -by*x/8.0/M_PI/(1.0-nu)*((1.0/Wr+C*C-z2/Wr3)*cosA+z/r3)
                      +bz*x*sinA/8.0/M_PI/(1.0-nu)*(1.0/Wr+C*C-z2/Wr3);

    Exy = bx*(rFi_ry)/2.0+by*(rFi_rx)/2.0  -bx/8.0/M_PI/(1.0-nu)*(x*y2/r2z2-nu*x/rz+x*y2/r3z-nu*x*cosA/Wr+eta*x*S/Wr+eta*x*y/Wr3)
                                           +by/8.0/M_PI/(1.0-nu)*(x2*y/r2z2-nu*y/rz+x2*y/r3z+nu*cosA*S+x2*y*cosA/Wr3+x2*cosA*S/Wr)
                                           -bz*sinA/8.0/M_PI/(1.0-nu)*(nu*S+x2*S/Wr+x2*y/Wr3);

    Exz = bx*(rFi_rz)/2.0+bz*(rFi_rx)/2.0  -bx/8.0/M_PI/(1.0-nu)*(-x*y/r3+nu*x*sinA/Wr+eta*x*C/Wr+eta*x*z/Wr3)
                                           +by/8.0/M_PI/(1.0-nu)*(-x2/r3+nu/r+nu*cosA*C+x2*z*cosA/Wr3+x2*cosA*C/Wr)
                                           -bz*sinA/8.0/M_PI/(1.0-nu)*(nu*C+x2*C/Wr+x2*z/Wr3);

    Eyz = by*(rFi_rz)/2.0+bz*(rFi_ry)/2.0  +bx/8.0/M_PI/(1.0-nu)*(y2/r3-nu/r-nu*cosA*C+nu*sinA*S+eta*sinA*cosA/W2-eta*(y*cosA+z*sinA)/W2r+eta*y*z/W2r2-eta*y*z/Wr3)
                                           -by*x/8.0/M_PI/(1.0-nu)*(y/r3+sinA*cosA*cosA/W2-cosA*(y*cosA+z*sinA)/W2r+y*z*cosA/W2r2-y*z*cosA/Wr3)
                                           -bz*x*sinA/8.0/M_PI/(1.0-nu)*(y*z/Wr3-sinA*cosA/W2+(y*cosA+z*sinA)/W2r-y*z/W2r2);

    e[0] = Exx;             e[1] = Exy;             e[2] = Exz;
    e[3] = Eyy;             e[4] = Eyz;             e[5] = Ezz;

    return;
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

