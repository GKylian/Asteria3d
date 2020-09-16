#pragma once
#include <iostream>
#define _USE_MATH_DEFINES
#include <math.h>
using namespace std;

#define M_rt4PI 3.54490770181103205460
const double gamma = 5.0/3.0;


//TODO://////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//TODO://////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//TODO://////////////////---  EIGENSYSTEM  ---///////////////////////////////////////////////////////////////////////////////
//TODO://////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



bool getEigenHD(long double *w, long double *eigen, long double L[][5], long double R[][5]) {
    long double rh = w[0];  long double vn = w[1];  long double P = w[4];
    if (rh <= 0) {
        cout << "RANGE-ERROR::getEigenHD:: rho <= 0" << endl; return false;
    }
    if (P <= 1e-6) {
        cout << "RANGE-ERROR::getEigenHD:: P <= 1e-6" << endl; return false;
    }
    long double a2 = gamma*P/rh; long double a = sqrtl(a2);

    eigen[0] = vn - a;
    eigen[1] = vn; eigen[2] = vn; eigen[3] = vn;
    eigen[4] = vn + a;

    R[0][0] = 1.0; R[1][0] = -a/rh; R[4][0] = a2;
    R[0][1] = 1.0;
    R[2][2] = 1.0;
    R[3][3] = 1.0;
    R[0][4] = 1.0; R[1][4] = a/rh; R[4][4] = a2;

    L[0][1] = -0.5*rh/a; L[0][4] = 0.5/a2;
    L[1][0] = 1.0; L[1][4] = -1.0/a2;
    L[2][2] = 1.0;
    L[3][3] = 1.0;
    L[4][1] = 0.5*rh/a; L[4][4] = 0.5/a2;


    return true;
}





//TODO:----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//TODO:----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//TODO:----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//TODO:----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//TODO:----------------------------------------------------------------------------------------------------------------------------------------------------------------------------




/// <summary>
/// Computes the eigensystem (eigenvalues as array, left/right eigenvectors as matrices) for the vector of primitive variables w in the n-direction.
/// </summary>
/// <param name="w">Primitive variables in order (rho, vn, v1, v2, P, B1, B2, Bn).</param>
/// <param name="eigen">The vector of eigenvalues.</param>
/// <param name="L">The matrix whose rows are the left eigenvectors.</param>
/// <param name="R">The matrix whose columns are the right eigenvectors.</param>
/// <returns>Has the computation worked ?</returns>
bool getEigen(long double *w, long double *eigen, long double L[][7], long double R[][7]) {

    //0) Shortcuts for the variables.
    long double rh = w[0]; long double vn = w[1]; long double P = w[4]; long double rtrh = sqrtl(rh);
    long double B1 = w[5]; long double B2 = w[6]; long double Bn = w[7];
    long double b1 = B1/rtrh; long double b2 = B2/rtrh; long double bn = b1/rtrh;


    //1) Compute some reused values.
    long double Sn = Sgn(Bn);

    long double a2 = gamma*P/rh; long double a = sqrtl(a2); if (isnan(a)) cout << "NaN::EIGENVALUES.H::getEigen:: wavespeed a is NaN: a2 = " << a2 << endl;
    long double CAn2 = Bn*Bn/rh; long double CAn = sqrtl(CAn2); if (isnan(CAn)) cout << "NaN::EIGENVALUES.H::getEigen:: wavespeed CAn is NaN: CAn2 = " << CAn2 << endl;
    long double CA2 = (Bn*Bn+B1*B1+B2*B2)/rh; long double CA = sqrtl(CA2); if (isnan(CA)) cout << "NaN::EIGENVALUES.H::getEigen:: wavespeed CA is NaN: CA2 = " << CA2 << endl;
    long double Cf2 = 0.5*(  (a2+CA2) + sqrtl(powl(a2+CA2, 2) - 4*a2*CAn2)  ); long double Cf = sqrtl(Cf2); if (isnan(Cf)) cout << "NaN::EIGENVALUES.H::getEigen:: wavespeed Cf is NaN: Cf2 = " << Cf2 << endl;
    long double Cs2 = 0.5*(  (a2+CA2) - sqrtl(powl(a2+CA2, 2) - 4*a2*CAn2)  ); long double Cs = sqrtl(Cs2); if (isnan(Cs)) cout << "NaN::EIGENVALUES.H::getEigen:: wavespeed Cs is NaN: Cs2 = " << Cs2 << endl;

    long double alf = sqrtl((a2 - Cs2)/(Cf2 - Cs2)); long double als = sqrtl((Cf2 - a2)/(Cf2 - Cs2));
    //cout << alf*alf+als*als-1.0 << endl;

    long double Bperp = sqrtl(B1*B1+B2*B2); long double bperp = Bperp/rtrh;
    long double beta1 = B1/Bperp; long double beta2 = B2/Bperp;

    if (null(Bperp) && !null(Bn)) {
        beta1 = 1.0/sqrtl(2.0); beta2 = 1.0/sqrtl(2.0);
    }
    if (null(Bn)) Sn = 1.0;



    //Handle special cases:
    if (null(Bn)) { //Bn = 0
        Cf2 = a2+Bperp*Bperp; Cf = sqrtl(Cf2); Cs2 = a2*Bn*Bn/(a2+Bperp*Bperp); Cs = sqrtl(Cs2);
        alf = sqrtl(a2/(a2+Bperp*Bperp)); als = sqrtl(Bperp*Bperp/(a2+Bperp*Bperp));
    }

    if (null(Bn) && null(Bperp)) { //HD limit
        //HD limit:
        Cf2 = a2; Cs2 = 0.0; Cf = a; Cs = 0.0;
        alf = 1.0; als = 0.0;
        beta1 = 1.0/sqrtl(2.0); beta2 = 1.0/sqrtl(2.0);
    }

    if (fabsl(bn)-a < -1e-16) { //fast waves = acoustic waves
        Cf2 = a2; Cf = a; Cs2 = Bn*Bn; Cs = Bn;
        alf = 1.0; als = 0.0;
    }
    
    if (fabsl(bn)-a > 1e-16) { //slow waves = acoustic waves
        Cf2 = Bn*Bn; Cf = Bn; Cs2 = a2; Cs = a;
        alf = 0.0; als = 1.0;
    }

    if (null(fabsl(bn)-a) && null(Bperp)) { //magnetosonic case
        Cf2 = a2; Cf = a; Cs2 = a2; Cs = a;
        long double phi = atanl( Bperp/(fabsl(bn)-a) );
        //alf = sinl(phi/2.0) + 
        cout << "Magnetosonic case" << endl;
    }


    if (isnan(alf)) cout << "NaN::EIGENVALUES.H::getEigen:: alpha_f is NaN; (alpha_f)^2 = " << (a2 - Cs2)/(Cf2 - Cs2) << endl;
    if (isnan(als)) cout << "NaN::EIGENVALUES.H::getEigen:: alpha_s is NaN; (alpha_s)^2 = " << (Cf2 - a2)/(Cf2 - Cs2) << endl;
    if (isnan(als)) {
        cout << "B(Bn, B1, B2) = " << Bn << ", " << B1 << ", " << B2 << endl;
    }

    long double Cff = Cf*alf; long double Css = Cs*als;
    long double Qf = Cf*alf*Sn; long double Qs = Cs*als*Sn;
    long double Af = a*alf*rtrh; long double As = a*als*rtrh;
    long double N = 0.5/a2;


    //2) Compute the eigenvalues
    eigen[0] = vn - Cf; eigen[1] = vn - CA; eigen[2] = vn - Cs;
    eigen[3] = vn;
    eigen[4] = vn + Cs; eigen[5] = vn + CA; eigen[6] = vn + Cf;


    //3) Compute the left eigenvectors (as rows in a 7*7 matrix)

    L[0][1] = -N*Cff; L[0][2] = N*Qs*beta1; L[0][3] = N*Qs*beta2; L[0][4] = N*alf/rh; L[0][5] = N*As*beta1/rh; L[0][6] = N*As*beta2/rh;
    L[1][2] = -0.5*beta2; L[1][3] = 0.5*beta1; L[1][5] = -0.5*beta2*Sn/rtrh; L[1][6] = 0.5*beta1*Sn/rtrh;
    L[2][1] = -N*Css; L[2][2] = -N*Qf*beta1; L[2][3] = -N*Qf*beta2; L[2][4] = N*als/rh; L[2][5] = -N*Af*beta1/rh; L[2][6] = -N*Af*beta2/rh;
    L[3][0] = 1.0; L[3][4] = -1.0/a2;
    L[4][1] = N*Css; L[4][2] = N*Qf*beta1; L[4][3] = N*Qf*beta2; L[4][4] = N*als/rh; L[4][5] = -N*Af*beta1/rh; L[4][6] = -N*Af*beta2/rh;
    L[5][2] = 0.5*beta2; L[5][3] = -0.5*beta1; L[5][5] = -0.5*beta2*Sn/rtrh; L[5][6] = 0.5*beta1*Sn/rtrh;
    L[6][1] = N*Cff; L[6][2] = -N*Qs*beta1; L[6][3] = -N*Qs*beta2; L[6][4] = N*alf/rh; L[6][5] = N*As*beta1/rh; L[6][6] = N*As*beta2/rh;


    R[0][0] = rh*alf; R[1][0] = -Cff; R[2][0] = Qs*beta1; R[3][0] = Qs*beta2; R[4][0] = rh*a2*alf; R[5][0] = As*beta1; R[6][0] = As*beta2;
    R[2][1] = -beta2; R[3][1] = beta1; R[5][1] = -beta2*Sn*rtrh; R[6][1] = beta1*Sn*rtrh;
    R[0][2] = rh*als; R[1][2] = -Css; R[2][2] = -Qf*beta1; R[3][2] = -Qf*beta2; R[4][2] = rh*a2*als; R[5][2] = -Af*beta1; R[6][2] = -Af*beta2;
    R[0][3] = 1.0;
    R[0][4] = rh*als; R[1][4] = Css; R[2][4] = Qf*beta1; R[3][4] = Qf*beta2; R[4][4] = rh*a2*als; R[5][4] = -Af*beta1; R[6][4] = -Af*beta2;
    R[2][5] = beta2; R[3][5] = -beta1; R[5][5] = -beta2*Sn*rtrh; R[6][5] = beta1*Sn*rtrh;
    R[0][6] = rh*alf; R[1][6] = Cff; R[2][6] = -Qs*beta1; R[3][6] = -Qs*beta2; R[4][6] = rh*a2*alf; R[5][6] = As*beta1; R[6][6] = As*beta2;





    //4) Check that nothing is NaN
    for (int m = 0; m < 7; m++)
    {
        if (isnan(eigen[m])) {
            cout << "NaN::EIGENVALUES.h::getEigen:: the eigenvalue eigen[" << m << "] is NaN !" << endl; return false;
        }
        for (int n = 0; n < 7; n++)
        {
            if (isnan(R[m][n])) {
                cout << "NaN::EIGENVALUES.h::getEigen:: the right eigenvector component R[" << m << "][" << n << "] is NaN !" << endl; return false;
            }
            if (isnan(L[m][n])) {
                cout << "NaN::EIGENVALUES.h::getEigen:: the left eigenvector component L[" << m << "][" << n << "] is NaN !" << endl; return false;
            }

            /*if (L[m][n] != 0.0 && null(L[m][n])) cout << "L[" << m << "][" << n << "] = " << L[m][n] << endl;
            if (R[m][n] != 0.0 && null(R[m][n])) cout << "R[" << m << "][" << n << "] = " << R[m][n] << endl;*/
            if (null(L[m][n])) L[m][n] = 0.0;
            if (null(R[m][n])) R[m][n] = 0.0;
        }
    }


    return true;
}



//TODO://////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//TODO://////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//TODO://////////////////---  EIGENVALUES  ---////////////////////////////////////////////////////////////////////////////////
//TODO://////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



bool getXEigenvalues(long double u_L[8], long double u_R[8], long double *eigen) // --OK--
{
    


    /// 0. Variables used later
    //long double rt4PI = sqrtl(4 * M_PI);

    long double rho = sqrtl(u_L[0])*sqrtl(u_R[0]);
    if (rho <= 0 || u_R[0] <= 0 || u_L[0] <= 0) {
        cout << "SOLVER::EIGENVALUES::ERROR::The density is null or negative" << endl;
        return false;
    }

    long double vx = (sqrtl(u_L[0])*u_L[1]/u_L[0] + sqrtl(u_R[0])*u_R[1]/u_R[0]) / (sqrtl(u_L[0]) + sqrtl(u_R[0]));
    long double vy = (sqrtl(u_L[0])*u_L[2]/u_L[0] + sqrtl(u_R[0])*u_R[2]/u_R[0]) / (sqrtl(u_L[0]) + sqrtl(u_R[0]));
    long double vz = (sqrtl(u_L[0])*u_L[3]/u_L[0] + sqrtl(u_R[0])*u_R[3]/u_R[0]) / (sqrtl(u_L[0]) + sqrtl(u_R[0]));
    long double By = (sqrtl(u_L[0])*u_L[6] + sqrtl(u_R[0])*u_R[6]) / (sqrtl(u_L[0]) + sqrtl(u_R[0]));
    long double Bz = (sqrtl(u_L[0])*u_L[7] + sqrtl(u_R[0])*u_R[7]) / (sqrtl(u_L[0]) + sqrtl(u_R[0]));
    long double Bx = u_L[5];
    
    long double vL = sqrtl(powl(u_L[1]/u_L[0],2) + powl(u_L[2]/u_L[0],2) + powl(u_L[3]/u_L[0],2));
    long double BL = sqrtl(u_L[5]*u_L[5] + u_L[6]*u_L[6] + u_L[7]*u_L[7]);
    long double PL = (u_L[4] - u_L[0]*vL*vL/2 - BL*BL/2) * (gamma - 1);
    if (PL < 0) {
        cout << "X-EIGENVALUES::WARNING::The pressure at left interface is negative: " << PL << endl; return false; }
    long double HL = (u_L[4] + PL + BL*BL/2) / u_L[0];

    long double vR = sqrtl(powl(u_R[1]/u_R[0], 2) + powl(u_R[2]/u_R[0], 2) + powl(u_R[3]/u_R[0], 2));
    long double BR = sqrtl(u_R[5]*u_R[5] + u_R[6]*u_R[6] + u_R[7]*u_R[7]);
    long double PR = (u_R[4] - u_R[0]*vR*vR/2 - BR*BR/2) * (gamma - 1);
    if (PR < 0) {
        cout << "X-EIGENVARUES::WARNING::The pressure at left interface is negative: " << PR << endl; return false;
    }
    long double HR = (u_R[4] + PR + BR*BR/2) / u_R[0];

    long double H = (sqrtl(u_L[0])*HL + sqrtl(u_R[0])*HR) / (sqrtl(u_L[0]) + sqrtl(u_R[0]));
    long double v = sqrtl(vx*vx + vy*vy + vz*vz); long double B = sqrtl(Bx*Bx + By*By + Bz*Bz);

    
    

    ///1. Compute in preparation for interface wavespeeds
    long double X = (powl(u_L[6]-u_R[6], 2) + powl(u_L[7]-u_R[7], 2)) / (2*(sqrtl(u_L[0]) + sqrtl(u_R[0])));
    //if (isnan(X)) cout << "Solver: X is NaN" << endl;
    long double Y = (u_L[0] + u_R[0]) / (2*rho);  if (isnan(Y)) cout << "Solver: Y is NaN" << endl;
    long double _gamma = gamma - 1;
    long double _X = (gamma - 2) * X; long double _Y = (gamma - 2) * Y; // < 0

    ///2. Compute the interface wavespeeds
    long double b_perp = sqrtl((_gamma - _Y)*(By*By+Bz*Bz)); if (isnan(b_perp)) { cout << "Solver: b_perp is NaN" << endl; return false; }
    long double CAx = sqrtl(Bx*Bx/rho);
    long double tCA = sqrtl(CAx*CAx + b_perp*b_perp/rho); if (isnan(tCA)) { cout << "Solver: tCA is NaN" << endl; return false; }
    long double ta = sqrtl(_gamma*(H-v*v/2-B*B/rho) - _X); if (isnan(ta)) { cout << "Solver: ta is NaN" << endl; return false; }

    long double Cs = sqrtl(0.5*((ta*ta+tCA*tCA) - sqrtl(powl(ta*ta+tCA*tCA, 2) - 4*ta*ta*CAx*CAx)));
    long double Cf = sqrtl(0.5*((ta*ta+tCA*tCA) + sqrtl(powl(ta*ta+tCA*tCA, 2) - 4*ta*ta*CAx*CAx)));

    //if ((ta*ta+tCA*tCA) - sqrtl(pow(ta*ta+tCA*tCA, 2) - 4*ta*ta*CAx*CAx) < 0) cout << "Solver: Cs is going to be NaN" << endl;
    //if (pow(ta*ta+tCA*tCA, 2) - 4*ta*ta*CAx*CAx < 0) cout << "Solver: Cs and Cf are going to be NaN" << endl;
    //if (isnan(Cs)) cout << "Solver: Cs is NaN" << endl; if (isnan(CAx)) cout << "Solver: CAx is NaN" << endl; if (isnan(Cf)) cout << "Solver: Cf is NaN" << endl;

    ///3. Update eigenvalues
    eigen[0] = vx - Cf; eigen[1] = vx - CAx; eigen[2] = vx - Cs;
    eigen[3] = vx;
    eigen[6] = vx + Cf; eigen[5] = vx + CAx; eigen[4] = vx + Cs;



    return true;
}


bool getYEigenvalues(long double u_L[8], long double u_R[8], long double *eigen)
{
    long double rho = sqrtl(u_L[0])*sqrtl(u_R[0]);
    if (rho <= 0 || u_R[0] <= 0 || u_L[0] <= 0) {
        cout << "SOLVER::EIGENVALUES::ERROR::The density is null or negative" << endl;
        return false;
    }

    long double vx = (sqrtl(u_L[0])*u_L[1]/u_L[0] + sqrtl(u_R[0])*u_R[1]/u_R[0]) / (sqrtl(u_L[0]) + sqrtl(u_R[0]));
    long double vy = (sqrtl(u_L[0])*u_L[2]/u_L[0] + sqrtl(u_R[0])*u_R[2]/u_R[0]) / (sqrtl(u_L[0]) + sqrtl(u_R[0]));
    long double vz = (sqrtl(u_L[0])*u_L[3]/u_L[0] + sqrtl(u_R[0])*u_R[3]/u_R[0]) / (sqrtl(u_L[0]) + sqrtl(u_R[0]));
    long double Bx = (sqrtl(u_L[0])*u_L[5] + sqrtl(u_R[0])*u_R[5]) / (sqrtl(u_L[0]) + sqrtl(u_R[0]));
    long double Bz = (sqrtl(u_L[0])*u_L[7] + sqrtl(u_R[0])*u_R[7]) / (sqrtl(u_L[0]) + sqrtl(u_R[0]));
    long double By = u_L[6];

    long double vL = sqrtl(powl(u_L[1]/u_L[0], 2) + powl(u_L[2]/u_L[0], 2) + powl(u_L[3]/u_L[0], 2));
    long double BL = sqrtl(u_L[5]*u_L[5] + u_L[6]*u_L[6] + u_L[7]*u_L[7]);
    long double PL = (u_L[4] - u_L[0]*vL*vL/2 - BL*BL/2) * (gamma - 1);
    if (PL < 0) {
        cout << "X-EIGENVALUES::WARNING::The pressure at left interface is negative: " << PL << endl; return false;
    }
    long double HL = (u_L[4] + PL + BL*BL/2) / u_L[0];

    long double vR = sqrtl(powl(u_R[1]/u_R[0], 2) + powl(u_R[2]/u_R[0], 2) + powl(u_R[3]/u_R[0], 2));
    long double BR = sqrtl(u_R[5]*u_R[5] + u_R[6]*u_R[6] + u_R[7]*u_R[7]);
    long double PR = (u_R[4] - u_R[0]*vR*vR/2 - BR*BR/2) * (gamma - 1);
    if (PR < 0) {
        cout << "X-EIGENVARUES::WARNING::The pressure at left interface is negative: " << PR << endl; return false;
    }
    long double HR = (u_R[4] + PR + BR*BR/2) / u_R[0];

    long double H = (sqrtl(u_L[0])*HL + sqrtl(u_R[0])*HR) / (sqrtl(u_L[0]) + sqrtl(u_R[0]));
    long double v = sqrtl(vx*vx + vy*vy + vz*vz); long double B = sqrtl(Bx*Bx + By*By + Bz*Bz);


    ///1. Compute in preparation for interface wavespeeds
    long double X = (powl(u_L[5]-u_R[5], 2) + powl(u_L[7]-u_R[7], 2)) / (2*(sqrtl(u_L[0]) + sqrtl(u_R[0])));
    if (isnan(X)) cout << "Solver: X is NaN" << endl;
    long double Y = (u_L[0] + u_R[0]) / (2*rho);  if (isnan(Y)) cout << "Solver: Y is NaN" << endl;
    long double _gamma = gamma - 1;
    long double _X = (gamma - 2) * X; long double _Y = (gamma - 2) * Y; // < 0

    ///2. Compute the interface wavespeeds
    long double b_perp = sqrtl((_gamma - _Y)*(Bx*Bx+Bz*Bz)); if (isnan(b_perp)) { cout << "Solver: b_perp is NaN" << endl; return false; }
    long double CAy = sqrtl(By*By/rho);
    long double tCA = sqrtl(CAy*CAy + b_perp*b_perp/rho); if (isnan(tCA)) { cout << "Solver: tCA is NaN" << endl; return false; }
    long double ta = sqrtl(_gamma*(H-v*v/2-B*B/rho) - _X); if (isnan(ta)) { cout << "Solver: ta is NaN" << endl; return false; }

    long double Cs = sqrtl(0.5*((ta*ta+tCA*tCA) - sqrtl(powl(ta*ta+tCA*tCA, 2) - 4*ta*ta*CAy*CAy)));
    long double Cf = sqrtl(0.5*((ta*ta+tCA*tCA) + sqrtl(powl(ta*ta+tCA*tCA, 2) - 4*ta*ta*CAy*CAy)));

    //if ((ta*ta+tCA*tCA) - sqrtl(pow(ta*ta+tCA*tCA, 2) - 4*ta*ta*CAy*CAy) < 0) cout << "Solver: Cs is going to be NaN" << endl;
    //if (pow(ta*ta+tCA*tCA, 2) - 4*ta*ta*CAy*CAy < 0) cout << "Solver: Cs and Cf are going to be NaN" << endl;
    //if (isnan(Cs)) cout << "Solver: Cs is NaN" << endl; if (isnan(CAy)) cout << "Solver: CAy is NaN" << endl; if (isnan(Cf)) cout << "Solver: Cf is NaN" << endl;

    ///3. Update eigenvalues
    eigen[0] = vy - Cf; eigen[1] = vy - CAy; eigen[2] = vy - Cs;
    eigen[3] = vy;
    eigen[6] = vy + Cf; eigen[5] = vy + CAy; eigen[4] = vy + Cs;



    return true;
}