#pragma once
#include "../Defs/Arrays.h"
#include "../prob.h"
#include "../Defs/phymaths.h"

#ifdef MHD
#ifdef ISO


/* Get eigenvalues for isometric MHD system */
void getEigen(long double *w, long double *eigen, long double L[][NWAVE], long double R[][NWAVE]) {

}



#else // Not ISO //

/* Get eigenvalues for adiabatic MHD system */
bool getEigen(long double *w, long double *eigen, long double L[][NWAVE], long double R[][NWAVE]) {

    //0) Shortcuts for the variables.
    long double rh = w[0]; long double vn = w[1]; long double P = w[4]; long double rtrh = sqrtl(rh);
    long double B1 = w[5]; long double B2 = w[6]; long double Bn = w[7];
    long double b1 = B1/rtrh; long double b2 = B2/rtrh; long double bn = b1/rtrh;

    if (rh <= 0) {
        cout << "Density is null or negative: " << rh << endl; return false;
    }
    if (P <= 0) {
        cout << "Pressure is null or negative: " << P << endl; return false;
    }
    

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



#endif // ISO  //
#else // Not MHD //
#ifdef ISO


/* Get eigenvalues for isometric HD system */
void getEigen() {

}



#else // Not ISO //

/* Get eigenvalues for adiabatic HD system */
void getEigen() {
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



#endif // ISO  //

#endif // MHD //   
