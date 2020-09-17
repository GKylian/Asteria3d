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



bool BCenter(Arrays *u, int i, int j, int k) {
    if (u->Nx > 1)
        u->uP(5, i, j, k) = (u->uC(5, i, j, k)+u->uC(5, i+1, j, k))/2.0;
    if(u->Ny > 1)
        u->uP(6, i, j, k) = (u->uC(6, i, j, k)+u->uC(6, i, j+1, k))/2.0;
    if(u->Nz > 1)
        u->uP(7, i, j, k) = (u->uC(7, i, j, k)+u->uC(7, i, j, k+1))/2.0;

    return true;
}


/*  uP to uC  */
bool toConserved(Arrays *u, int i, int j, int k, long double gamma) {
    long double v2 = SQ(u->uP(1, i, j, k))+SQ(u->uP(2, i, j, k))+SQ(u->uP(3, i, j, k));
    long double B2 = SQ(u->uP(5, i, j, k))+SQ(u->uP(6, i, j, k))+SQ(u->uP(7, i, j, k));
    
    long double rh = u->uP(0, i, j, k); u->uC(0, i, j, k) = rh;
    u->uC(1, i, j, k) = rh*u->uP(1, i, j, k); u->uC(2, i, j, k) = rh*u->uP(2, i, j, k); u->uC(3, i, j, k) = rh*u->uP(3, i, j, k);
    
    u->uC(4, i, j, k) = u->uP(4, i, j, k)/(gamma-1) + rh*v2/2.0 + B2/2.0;

    if (isnan(u->uC(4, i, j, k))) return false;
    return true;
}

/*  uC to uP  */
bool toPrimitive(Arrays *u, int i, int j, int k, long double gamma) {
    long double rh = u->uC(0, i, j, k); u->uP(0, i, j, k) = rh;
    u->uP(1, i, j, k) = u->uC(1, i, j, k)/rh; u->uP(2, i, j, k) = u->uC(2, i, j, k)/rh; u->uP(3, i, j, k) = u->uC(3, i, j, k)/rh;

    long double v2 = SQ(u->uP(1, i, j, k))+SQ(u->uP(2, i, j, k))+SQ(u->uP(3, i, j, k));
    long double B2 = SQ(u->uP(5, i, j, k))+SQ(u->uP(6, i, j, k))+SQ(u->uP(7, i, j, k));

    u->uP(4, i, j, k) = ( u->uC(4, i, j, k) - rh*v2/2.0 - B2/2.0 ) * (gamma-1);

    if (isnan(u->uP(4, i, j, k))) return false;
    return true;
}

/* Transform the (i, j, k) interface values from primitive to conserved variables (used at the end of the reconstruction) */ 
bool toConserved(Arrays *u, int dim, int i, int j, int k, long double gamma) {
    if (dim == 0) {
        long double v2 = SQ(u->ix(1, i, j, k))+SQ(u->ix(2, i, j, k))+SQ(u->ix(3, i, j, k));
        long double B2 = SQ(u->ix(5, i, j, k))+SQ(u->ix(6, i, j, k))+SQ(u->ix(7, i, j, k));

        long double rh = u->ix(0, i, j, k); u->ix(0, i, j, k) = rh;
        u->ix(1, i, j, k) = rh*u->ix(1, i, j, k); u->ix(2, i, j, k) = rh*u->ix(2, i, j, k); u->ix(3, i, j, k) = rh*u->ix(3, i, j, k);

        u->ix(4, i, j, k) = u->ix(4, i, j, k)/(gamma-1) + rh*v2/2.0 + B2/2.0;

        if (isnan(u->ix(4, i, j, k))) return false;
        return true;
    }
    if (dim == 1) {
        long double v2 = SQ(u->iy(1, i, j, k))+SQ(u->iy(2, i, j, k))+SQ(u->iy(3, i, j, k));
        long double B2 = SQ(u->iy(5, i, j, k))+SQ(u->iy(6, i, j, k))+SQ(u->iy(7, i, j, k));

        long double rh = u->iy(0, i, j, k); u->iy(0, i, j, k) = rh;
        u->iy(1, i, j, k) = rh*u->iy(1, i, j, k); u->iy(2, i, j, k) = rh*u->iy(2, i, j, k); u->iy(3, i, j, k) = rh*u->iy(3, i, j, k);

        u->iy(4, i, j, k) = u->iy(4, i, j, k)/(gamma-1) + rh*v2/2.0 + B2/2.0;

        if (isnan(u->iy(4, i, j, k))) return false;
        return true;
    }
    if (dim == 2) {
        long double v2 = SQ(u->iz(1, i, j, k))+SQ(u->iz(2, i, j, k))+SQ(u->iz(3, i, j, k));
        long double B2 = SQ(u->iz(5, i, j, k))+SQ(u->iz(6, i, j, k))+SQ(u->iz(7, i, j, k));

        long double rh = u->iz(0, i, j, k); u->iz(0, i, j, k) = rh;
        u->iz(1, i, j, k) = rh*u->iz(1, i, j, k); u->iz(2, i, j, k) = rh*u->iz(2, i, j, k); u->iz(3, i, j, k) = rh*u->iz(3, i, j, k);

        u->iz(4, i, j, k) = u->iz(4, i, j, k)/(gamma-1) + rh*v2/2.0 + B2/2.0;

        if (isnan(u->iz(4, i, j, k))) return false;
        return true;
    }
    
}



bool F(long double uP[8], long double uC[8], long double *flux) {
    long double v = sqrtl(SQ(uP[1])+SQ(uP[2])+SQ(uP[3]));
    long double B = sqrtl(SQ(uP[5])+SQ(uP[6])+SQ(uP[7]));
    long double Pstar = uP[4] + B*B/2.0;

    flux[0] = uC[1];
    flux[1] = uC[1]*uP[1] + Pstar - uP[5]*uP[5];
    flux[2] = uC[1]*uP[2] - uP[5]*uP[6];
    flux[3] = uC[1]*uP[3] - uP[5]*uP[7];
    flux[4] = (uP[4]+Pstar)*uP[1] - (uP[5]*uP[1]+uP[6]*uP[2]+uP[7]*uP[3])*uP[5];
    flux[5] = 0.0;
    flux[6] = uP[6]*uP[1] - uP[5]*uP[2];
    flux[7] = uP[7]*uP[1] - uP[5]*uP[3];

    return true;
}

bool G(long double uP[8], long double uC[8], long double *flux) {
    long double v = sqrtl(SQ(uP[1])+SQ(uP[2])+SQ(uP[3]));
    long double B = sqrtl(SQ(uP[5])+SQ(uP[6])+SQ(uP[7]));
    long double Pstar = uP[4] + B*B/2.0;

    flux[0] = uC[2];
    flux[1] = uC[2]*uP[1] -         uP[6]*uP[5];
    flux[2] = uC[2]*uP[2] + Pstar - uP[6]*uP[6];
    flux[3] = uC[2]*uP[3] -         uP[6]*uP[7];
    flux[4] = (uC[4]+Pstar)*uP[2] - (uP[5]*uP[1]+uP[6]*uP[2]+uP[7]*uP[3])*uP[6];
    flux[5] = uP[5]*uP[2] - uP[6]*uP[1];
    flux[6] = 0.0;
    flux[7] = uP[7]*uP[2] - uP[6]*uP[3];

    return true;
}

bool H(long double uP[8], long double uC[8], long double *flux) {
    long double v = sqrtl(SQ(uP[1])+SQ(uP[2])+SQ(uP[3]));
    long double B = sqrtl(SQ(uP[5])+SQ(uP[6])+SQ(uP[7]));
    long double Pstar = uP[4] + B*B/2.0;

    flux[0] = uC[3];
    flux[1] = uC[3]*uP[1] -         uP[7]*uP[5];
    flux[2] = uC[3]*uP[2] -         uP[7]*uP[6];
    flux[3] = uC[3]*uP[3] + Pstar - uP[7]*uP[7];
    flux[4] = (uC[4]+Pstar)*uP[3] - (uP[5]*uP[1]+uP[6]*uP[2]+uP[7]*uP[3])*uP[7];
    flux[5] = uP[5]*uP[3] - uP[7]*uP[1];
    flux[6] = uP[6]*uP[3] - uP[7]*uP[2];
    flux[7] = 0.0;

    return true;
}


long double getDiv(Arrays *u, int i, int j, int k) {
    
    return (u->Nx>1)*(u->uC(5, i+1, j, k)-u->uC(5, i, j, k))/u->dx + (u->Ny>1)*(u->uC(6, i, j+1, k)-u->uC(6, i, j, k))/u->dy + (u->Nz>1)*(u->uC(7, i, j, k+1)-u->uC(7, i, j, k))/u->dz;
}