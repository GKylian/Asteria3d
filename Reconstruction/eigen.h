#pragma once
#include "../Defs/Arrays.h"
#include "../prob.h"
#include "../Defs/phymaths.h"

#ifdef MHD
#ifdef ISO


/* Get eigenvalues for isometric MHD system */
void getEigen(ld *w, ld *eigen, ld L[][NWAVE], ld R[][NWAVE]) {

}
 


#else // Not ISO //

/* Get eigenvalues for adiabatic MHD system */
bool getEigen(ld *w, ld *eigen, ld L[][NWAVE], ld R[][NWAVE], ld gamma) {

    //0) Shortcuts for the variables.
    ld rh = w[0]; ld vn = w[1]; ld P = w[4]; ld rtrh = sqrtl(rh);
    ld B1 = w[5]; ld B2 = w[6]; ld Bn = w[7];
    ld b1 = B1/rtrh; ld b2 = B2/rtrh; ld bn = b1/rtrh;

    if (rh <= 0) {
        cout << "Density is null or negative: " << rh << endl; return false;
    }
    if (P <= 0) {
        cout << "Pressure is null or negative: " << P << endl; return false;
    }
    

    //1) Compute some reused values.
    ld Sn = Sgn(Bn);

    ld a2 = gamma*P/rh; ld a = sqrtl(a2); if (isnan(a)) cout << "NaN::EIGENVALUES.H::getEigen:: wavespeed a is NaN: a2 = " << a2 << endl;
    ld CAn2 = Bn*Bn/rh; ld CAn = sqrtl(CAn2); if (isnan(CAn)) cout << "NaN::EIGENVALUES.H::getEigen:: wavespeed CAn is NaN: CAn2 = " << CAn2 << endl;
    ld CA2 = (Bn*Bn+B1*B1+B2*B2)/rh; ld CA = sqrtl(CA2); if (isnan(CA)) cout << "NaN::EIGENVALUES.H::getEigen:: wavespeed CA is NaN: CA2 = " << CA2 << endl;
    ld Cf2 = 0.5*(  (a2+CA2) + sqrtl(powl(a2+CA2, 2) - 4*a2*CAn2)  ); ld Cf = sqrtl(Cf2); if (isnan(Cf)) cout << "NaN::EIGENVALUES.H::getEigen:: wavespeed Cf is NaN: Cf2 = " << Cf2 << endl;
    ld Cs2 = 0.5*(  (a2+CA2) - sqrtl(powl(a2+CA2, 2) - 4*a2*CAn2)  ); ld Cs = sqrtl(Cs2); if (isnan(Cs)) cout << "NaN::EIGENVALUES.H::getEigen:: wavespeed Cs is NaN: Cs2 = " << Cs2 << endl;

    ld alf = sqrtl((a2 - Cs2)/(Cf2 - Cs2)); ld als = sqrtl((Cf2 - a2)/(Cf2 - Cs2));
    //cout << alf*alf+als*als-1.0 << endl;

    ld Bperp = sqrtl(B1*B1+B2*B2); ld bperp = Bperp/rtrh;
    ld beta1 = B1/Bperp; ld beta2 = B2/Bperp;

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
        ld phi = atanl( Bperp/(fabsl(bn)-a) );
        //alf = sinl(phi/2.0) + 
        cout << "Magnetosonic case" << endl;
    }

    if (isnan(beta1)) {
        cout << "NaN:::eigen.h::getEigen:: beta 1 is NaN !" << endl;
        cout << "Magnetic field: " << w[5] << ", " << w[6] << ", " << w[7] << endl;
        cout << "Bperp: " << Bperp << ", Bn: " << Bn << endl;
        cout << "null(Bperp): " << null(Bperp) << ", null(Bn): " << null(Bn) << " --> null(Bperp) && !null(Bn): " << (null(Bperp) && !null(Bn)) << endl;
    }


    if (isnan(alf)) cout << "NaN::EIGENVALUES.H::getEigen:: alpha_f is NaN; (alpha_f)^2 = " << (a2 - Cs2)/(Cf2 - Cs2) << endl;
    if (isnan(als)) cout << "NaN::EIGENVALUES.H::getEigen:: alpha_s is NaN; (alpha_s)^2 = " << (Cf2 - a2)/(Cf2 - Cs2) << endl;
    if (isnan(als)) {
        cout << "B(Bn, B1, B2) = " << Bn << ", " << B1 << ", " << B2 << endl;
    }

    ld Cff = Cf*alf; ld Css = Cs*als; if (isnan(Cff)) cout << "NaN:::eigen.h::getEigen:: Cff is NaN !" << endl;   if (isnan(Css)) cout << "NaN:::eigen.h::getEigen:: Css is NaN !" << endl;
    ld Qf = Cf*alf*Sn; ld Qs = Cs*als*Sn; if (isnan(Qf)) cout << "NaN:::eigen.h::getEigen:: Qf is NaN !" << endl;   if (isnan(Qs)) cout << "NaN:::eigen.h::getEigen:: Qs is NaN !" << endl;
    ld Af = a*alf*rtrh; ld As = a*als*rtrh; if (isnan(Af)) cout << "NaN:::eigen.h::getEigen:: Af is NaN !" << endl;   if (isnan(As)) cout << "NaN:::eigen.h::getEigen:: Af is NaN !" << endl;
    ld N = 0.5/a2; if (isnan(N)) cout << "NaN:::eigen.h::getEigen:: a is NaN,   a2 = " << a2 << endl;


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
bool getEigen(ld *w, ld *eigen, ld L[][NWAVE], ld R[][NWAVE], ld gamma) {
    ld rh = w[0];  ld vn = w[1];  ld P = w[4];
    if (rh <= 0) {
        cout << "RANGE-ERROR::getEigenHD:: rho <= 0" << endl; return false;
    }
    if (P <= 1e-6) {
        cout << "RANGE-ERROR::getEigenHD:: P <= 1e-6" << endl; return false;
    }
    ld a2 = gamma*P/rh; ld a = sqrtl(a2);

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
