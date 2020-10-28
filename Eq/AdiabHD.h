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
bool toConserved(Arrays *u, int i, int j, int k, ld gamma) {
    ld v2 = SQ(u->uP(1, i, j, k))+SQ(u->uP(2, i, j, k))+SQ(u->uP(3, i, j, k));

    ld rh = u->uP(0, i, j, k); u->uC(0, i, j, k) = rh;
    u->uC(1, i, j, k) = rh*u->uP(1, i, j, k); u->uC(2, i, j, k) = rh*u->uP(2, i, j, k); u->uC(3, i, j, k) = rh*u->uP(3, i, j, k);

    u->uC(4, i, j, k) = u->uP(4, i, j, k)/(gamma-1) + rh*v2/2.0;

    if (isnan(u->uC(4, i, j, k))) return false;
    return true;
}

bool toConserved(ld *w, ld *u, ld gamma) {
    ld v2 = SQ(w[1])+SQ(w[2])+SQ(w[3]);

    ld rh = w[0]; u[0] = rh;
    u[1] = rh*w[1]; u[2] = rh*w[2]; u[3] = rh*u[3];

    u[4] = w[4]/(gamma-1) + rh*v2/2.0;

    if (isnan(u[4]) || null(u[4])) {
        std::cout << "ERROR:::AdiabHD.h::toConserved:: E is NaN or null/negative: " << u[4] << std::endl;
        return false;
    }

}

/*  uC to uP  */
bool toPrimitive(Arrays *u, int i, int j, int k, ld gamma) {
    ld rh = u->uC(0, i, j, k); u->uP(0, i, j, k) = rh;
    u->uP(1, i, j, k) = u->uC(1, i, j, k)/rh; u->uP(2, i, j, k) = u->uC(2, i, j, k)/rh; u->uP(3, i, j, k) = u->uC(3, i, j, k)/rh;

    ld v2 = SQ(u->uP(1, i, j, k))+SQ(u->uP(2, i, j, k))+SQ(u->uP(3, i, j, k));

    u->uP(4, i, j, k) = (u->uC(4, i, j, k) - rh*v2/2.0) * (gamma-1);

    if (isnan(u->uP(4, i, j, k))) return false;
    return true;
}

bool toPrimitive(ld *u, ld *w, ld gamma) {
    if (null(u[0])) {
        std::cout << "ERROR:::AdiabHD.h::toPrimitive:: the density is null !" << std::endl;
        return false;
    }
    ld rh = u[0]; w[0] = rh;
    w[1] = u[1]/rh; w[2] = u[2]/rh; w[3] = u[3]/rh;

    ld v2 = SQ(w[1])+SQ(w[2])+SQ(w[3]);

    w[4] = (u[4] - rh*v2/2.0)*(gamma-1);

    if (isnan(w[4]) || null(w[4])) {
        std::cout << "ERROR:::AdiabHD.h::toPrimitive:: P is NaN or null/negative: " << w[4] << std::endl;
        std::cout << "--> rho = " << rh << ", E: " << u[4] << ", v2: " << v2 << std::endl;

        return false;
    }

    return true;

}

/* Transform the (i, j, k) interface values from primitive to conserved variables (used at the end of the reconstruction) */
bool toConserved(Arrays *u, int dim, int i, int j, int k, ld gamma) {
    if (dim == 0) {
        ld v2 = SQ(u->ix(1, i, j, k))+SQ(u->ix(2, i, j, k))+SQ(u->ix(3, i, j, k));

        ld rh = u->ix(0, i, j, k); u->ix(0, i, j, k) = rh;
        u->ix(1, i, j, k) = rh*u->ix(1, i, j, k); u->ix(2, i, j, k) = rh*u->ix(2, i, j, k); u->ix(3, i, j, k) = rh*u->ix(3, i, j, k);

        u->ix(4, i, j, k) = u->ix(4, i, j, k)/(gamma-1) + rh*v2/2.0;

        if (isnan(u->ix(4, i, j, k))) return false;
        return true;
    }
    if (dim == 1) {
        ld v2 = SQ(u->iy(1, i, j, k))+SQ(u->iy(2, i, j, k))+SQ(u->iy(3, i, j, k));

        ld rh = u->iy(0, i, j, k); u->iy(0, i, j, k) = rh;
        u->iy(1, i, j, k) = rh*u->iy(1, i, j, k); u->iy(2, i, j, k) = rh*u->iy(2, i, j, k); u->iy(3, i, j, k) = rh*u->iy(3, i, j, k);

        u->iy(4, i, j, k) = u->iy(4, i, j, k)/(gamma-1) + rh*v2/2.0;

        if (isnan(u->iy(4, i, j, k))) return false;
        return true;
    }
    if (dim == 2) {
        ld v2 = SQ(u->iz(1, i, j, k))+SQ(u->iz(2, i, j, k))+SQ(u->iz(3, i, j, k));

        ld rh = u->iz(0, i, j, k); u->iz(0, i, j, k) = rh;
        u->iz(1, i, j, k) = rh*u->iz(1, i, j, k); u->iz(2, i, j, k) = rh*u->iz(2, i, j, k); u->iz(3, i, j, k) = rh*u->iz(3, i, j, k);

        u->iz(4, i, j, k) = u->iz(4, i, j, k)/(gamma-1) + rh*v2/2.0;

        if (isnan(u->iz(4, i, j, k))) return false;
        return true;
    }

}



bool F(ld uP[8], ld uC[8], ld *flux) {
    ld v = sqrtl(SQ(uP[1])+SQ(uP[2])+SQ(uP[3]));
    ld Pstar = uP[4];

    flux[0] = uC[1];
    flux[1] = uC[1]*uP[1] + Pstar;
    flux[2] = uC[1]*uP[2];
    flux[3] = uC[1]*uP[3];
    flux[4] = (uC[4]+Pstar)*uP[1];

    return true;
}

bool G(ld uP[8], ld uC[8], ld *flux) {
    ld v = sqrtl(SQ(uP[1])+SQ(uP[2])+SQ(uP[3]));
    ld Pstar = uP[4];

    flux[0] = uC[2];
    flux[1] = uC[2]*uP[1];
    flux[2] = uC[2]*uP[2] + Pstar;
    flux[3] = uC[2]*uP[3];
    flux[4] = (uC[4]+Pstar)*uP[2];

    return true;
}

bool H(ld uP[8], ld uC[8], ld *flux) {
    ld v = sqrtl(SQ(uP[1])+SQ(uP[2])+SQ(uP[3]));
    ld Pstar = uP[4];

    flux[0] = uC[3];
    flux[1] = uC[3]*uP[1];
    flux[2] = uC[3]*uP[2];
    flux[3] = uC[3]*uP[3] + Pstar;
    flux[4] = (uC[4]+Pstar)*uP[3];

    return true;
}


#ifndef ISO
void getWavespeeds(Arrays *u, int i, int j, int k, ld gamma, ld *wvx, ld *wvy, ld *wvz) {

    ld v = sqrtl(SQ(u->uP(1, i, j, k))+SQ(u->uP(2, i, j, k))+SQ(u->uP(3, i, j, k)));


    ld rh = u->uC(0, i, j, k); if (rh <= 0) { cout << "ERROR:::AdiabHD.h::getWavespeeds:: The density is null or negative !" << endl; return; }
    ld a = sqrtl(gamma*u->uP(4, i, j, k)/rh);

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

void getWavespeeds(ld *w, ld gamma, ld *wvn) {

    /* Order: (rho, vn, v1, v2, P) */
    ld v = sqrtl(SQ(w[1])+SQ(w[2])+SQ(w[3]));


    ld rh = w[0]; if (rh <= 0) { cout << "ERROR:::AdiabHD.h::getWavespeeds:: The density is null or negative !" << endl; return; }
    ld a = sqrtl(gamma*w[4]/rh);

    if (isnan(a)) {
        cout << "ERROR:::AdiabHD.h::getWavespeeds:: a is NaN" << endl;
        cout << "d = " << rh << ", P = " << w[4] << endl;
        return;
    }
    wvn[0] = a;
    wvn[1] = 0.0;

}

ld maxWavespeed(ld *w, ld gamma) {

    /* Order: (rho, vn, v1, v2, P) */
    return sqrtl(gamma*w[4]/w[0]);
}
#endif // !ISO





void minmaxRoeEigenvalue(ld *uL, ld *wL, ld *uR, ld *wR, ld gamma, ld *l) {
    //ld rho = sqrtl(wL[0])*sqrtl(wR[0]);
    ld vx = ( sqrtl(wL[0])*wL[1] + sqrtl(wR[0])*wR[1] ) / (sqrtl(wL[0])+sqrtl(wR[0]));
    ld vy = ( sqrtl(wL[0])*wL[2] + sqrtl(wR[0])*wR[2] ) / (sqrtl(wL[0])+sqrtl(wR[0]));
    ld vz = ( sqrtl(wL[0])*wL[3] + sqrtl(wR[0])*wR[3] ) / (sqrtl(wL[0])+sqrtl(wR[0]));
    ld v2 = vx*vx + vy*vy + vz*vz;

    ld HL = (uL[4]+wL[4])/wL[0];   ld HR = (uR[4]+wR[4])/wR[0];
    ld H = ( sqrtl(wL[0])*HL + sqrtl(wR[0])*HR ) / (sqrtl(wL[0])+sqrtl(wR[0]));
    ld a = sqrtl((gamma-1)*(H-v2/2));
    l[0] = vx - a;
    l[1] = vx + a;
}