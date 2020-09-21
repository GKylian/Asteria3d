#pragma once
#include "../Defs/Arrays.h"
#include "../Defs/phymaths.h"


/*
<-------------------- Functions -------------------->
<->
---> bool BCenter(Arrays *uP, Arrays *uC, int i, int j, int k)
     <-> Computes B at cell center (in uP) from B at interfaces (in uC) for cell (i, j, k)
--->  - Arrays *u: the primitive variables, with B at cell center and the conserved variables, with B at interfaces
--->  - int i, j, k: the cell indices
<->
<->
---> bool toConserved(Arrays *u, int i, int j, int k)
     <-> Computes the conserved variables (in uC) from the primitive variables (in uP) at (i, j, k)
     <-> !!! We assume B is already at cell center in uP
--->  - Arrays *u: the primitive variables, with B at cell center and the conserved variables, with B at interfaces
<->
<->
---> bool toPrimitive(Arrays *uP, Arrays *uC, int i, int j, int k)
     <-> Computes the primitive variables (in uC) from the conserved variables (in uP) at (i, j, k)
     <-> !!! We assume B is already at cell center in uP
--->  - Arrays *u: the primitive variables, with B at cell center and the conserved variables, with B at interfaces
<->
*/


/*  uP to uC  */
bool toConserved(Arrays *u, int i, int j, int k, long double gamma) {
    long double v2 = SQ(u->uP(1, i, j, k))+SQ(u->uP(2, i, j, k))+SQ(u->uP(3, i, j, k));

    long double rh = u->uP(0, i, j, k); u->uC(0, i, j, k) = rh;
    u->uC(1, i, j, k) = rh*u->uP(1, i, j, k); u->uC(2, i, j, k) = rh*u->uP(2, i, j, k); u->uC(3, i, j, k) = rh*u->uP(3, i, j, k);

    u->uC(4, i, j, k) = u->uP(4, i, j, k)/(gamma-1) + rh*v2/2.0;

    if (isnan(u->uC(4, i, j, k))) return false;
    return true;
}

bool toConserved(long double *w, long double *u, long double gamma) {
    long double v2 = SQ(w[1])+SQ(w[2])+SQ(w[3]);

    long double rh = w[0]; u[0] = rh;
    u[1] = rh*w[1]; u[2] = rh*w[2]; u[3] = rh*u[3];

    u[4] = w[4]/(gamma-1) + rh*v2/2.0;

    if (isnan(u[4]) || null(u[4])) {
        std::cout << "ERROR:::AdiabHD.h::toConserved:: E is NaN or null/negative: " << u[4] << std::endl;
        return false;
    }

}

/*  uC to uP  */
bool toPrimitive(Arrays *u, int i, int j, int k, long double gamma) {
    long double rh = u->uC(0, i, j, k); u->uP(0, i, j, k) = rh;
    u->uP(1, i, j, k) = u->uC(1, i, j, k)/rh; u->uP(2, i, j, k) = u->uC(2, i, j, k)/rh; u->uP(3, i, j, k) = u->uC(3, i, j, k)/rh;

    long double v2 = SQ(u->uP(1, i, j, k))+SQ(u->uP(2, i, j, k))+SQ(u->uP(3, i, j, k));

    u->uP(4, i, j, k) = (u->uC(4, i, j, k) - rh*v2/2.0) * (gamma-1);

    if (isnan(u->uP(4, i, j, k))) return false;
    return true;
}

bool toPrimitive(long double *u, long double *w, long double gamma) {
    if (null(u[0])) {
        std::cout << "ERROR:::AdiabHD.h::toPrimitive:: the density is null !" << std::endl;
        return false;
    }
    long double rh = u[0]; w[0] = rh;
    w[1] = u[1]/rh; w[2] = u[2]/rh; w[3] = u[3]/rh;

    long double v2 = SQ(w[1])+SQ(w[2])+SQ(w[3]);

    w[4] = (u[4] - rh*v2/2.0)*(gamma-1);

    if (isnan(w[4]) || null(w[4])) {
        std::cout << "ERROR:::AdiabHD.h::toPrimitive:: P is NaN or null/negative: " << w[4] << std::endl;
        std::cout << "--> rho = " << rh << ", E: " << u[4] << ", v2: " << v2 << std::endl;

        return false;
    }

    return true;

}

/* Transform the (i, j, k) interface values from primitive to conserved variables (used at the end of the reconstruction) */
bool toConserved(Arrays *u, int dim, int i, int j, int k, long double gamma) {
    if (dim == 0) {
        long double v2 = SQ(u->ix(1, i, j, k))+SQ(u->ix(2, i, j, k))+SQ(u->ix(3, i, j, k));

        long double rh = u->ix(0, i, j, k); u->ix(0, i, j, k) = rh;
        u->ix(1, i, j, k) = rh*u->ix(1, i, j, k); u->ix(2, i, j, k) = rh*u->ix(2, i, j, k); u->ix(3, i, j, k) = rh*u->ix(3, i, j, k);

        u->ix(4, i, j, k) = u->ix(4, i, j, k)/(gamma-1) + rh*v2/2.0;

        if (isnan(u->ix(4, i, j, k))) return false;
        return true;
    }
    if (dim == 1) {
        long double v2 = SQ(u->iy(1, i, j, k))+SQ(u->iy(2, i, j, k))+SQ(u->iy(3, i, j, k));

        long double rh = u->iy(0, i, j, k); u->iy(0, i, j, k) = rh;
        u->iy(1, i, j, k) = rh*u->iy(1, i, j, k); u->iy(2, i, j, k) = rh*u->iy(2, i, j, k); u->iy(3, i, j, k) = rh*u->iy(3, i, j, k);

        u->iy(4, i, j, k) = u->iy(4, i, j, k)/(gamma-1) + rh*v2/2.0;

        if (isnan(u->iy(4, i, j, k))) return false;
        return true;
    }
    if (dim == 2) {
        long double v2 = SQ(u->iz(1, i, j, k))+SQ(u->iz(2, i, j, k))+SQ(u->iz(3, i, j, k));

        long double rh = u->iz(0, i, j, k); u->iz(0, i, j, k) = rh;
        u->iz(1, i, j, k) = rh*u->iz(1, i, j, k); u->iz(2, i, j, k) = rh*u->iz(2, i, j, k); u->iz(3, i, j, k) = rh*u->iz(3, i, j, k);

        u->iz(4, i, j, k) = u->iz(4, i, j, k)/(gamma-1) + rh*v2/2.0;

        if (isnan(u->iz(4, i, j, k))) return false;
        return true;
    }

}



bool F(long double uP[8], long double uC[8], long double *flux) {
    long double v = sqrtl(SQ(uP[1])+SQ(uP[2])+SQ(uP[3]));
    long double Pstar = uP[4];

    flux[0] = uC[1];
    flux[1] = uC[1]*uP[1] + Pstar;
    flux[2] = uC[1]*uP[2];
    flux[3] = uC[1]*uP[3];
    flux[4] = (uC[4]+Pstar)*uP[1];

    return true;
}

bool G(long double uP[8], long double uC[8], long double *flux) {
    long double v = sqrtl(SQ(uP[1])+SQ(uP[2])+SQ(uP[3]));
    long double Pstar = uP[4];

    flux[0] = uC[2];
    flux[1] = uC[2]*uP[1];
    flux[2] = uC[2]*uP[2] + Pstar;
    flux[3] = uC[2]*uP[3];
    flux[4] = (uC[4]+Pstar)*uP[2];

    return true;
}

bool H(long double uP[8], long double uC[8], long double *flux) {
    long double v = sqrtl(SQ(uP[1])+SQ(uP[2])+SQ(uP[3]));
    long double Pstar = uP[4];

    flux[0] = uC[3];
    flux[1] = uC[3]*uP[1];
    flux[2] = uC[3]*uP[2];
    flux[3] = uC[3]*uP[3] + Pstar;
    flux[4] = (uC[4]+Pstar)*uP[3];

    return true;
}


#ifndef ISO
void getWavespeeds(Arrays *u, int i, int j, int k, long double gamma, long double *wvx, long double *wvy, long double *wvz) {

    long double v = sqrtl(SQ(u->uP(1, i, j, k))+SQ(u->uP(2, i, j, k))+SQ(u->uP(3, i, j, k)));


    long double rh = u->uC(0, i, j, k); if (rh <= 0) { cout << "ERROR:::AdiabHD.h::getWavespeeds:: The density is null or negative !" << endl; return; }
    long double a = sqrtl(gamma*u->uP(4, i, j, k)/rh);

    if (isnan(a)) {
        cout << "ERROR:::AdiabHD.h::getWavespeeds:: a is NaN" << endl;
        cout << "d = " << rh << ", P = " << u->uP(4, i, j, k) << endl;
        return;
    }


    wvx[0] = a;
    wvx[1] = 0.0;

    wvy[0] = a;
    wvy[1] = 0.0;

    wvz[0] = a;
    wvz[1] = 0.0;

}

void getWavespeeds(long double *w, long double gamma, long double *wvn) {

    /* Order: (rho, vn, v1, v2, P) */
    long double v = sqrtl(SQ(w[1])+SQ(w[2])+SQ(w[3]));


    long double rh = w[0]; if (rh <= 0) { cout << "ERROR:::AdiabHD.h::getWavespeeds:: The density is null or negative !" << endl; return; }
    long double a = sqrtl(gamma*w[4]/rh);

    if (isnan(a)) {
        cout << "ERROR:::AdiabHD.h::getWavespeeds:: a is NaN" << endl;
        cout << "d = " << rh << ", P = " << w[4] << endl;
        return;
    }
    wvn[0] = a;
    wvn[1] = 0.0;

}
#endif // !ISO


