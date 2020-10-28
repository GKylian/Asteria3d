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
        u->uP(NVAL-3, i, j, k) = (u->uC(NVAL-3, i, j, k)+u->uC(NVAL-3, i+1, j, k))/2.0;
    else
        u->uP(NVAL-3, i, j, k) = u->uC(NVAL-3, i, j, k);

    if (u->Ny > 1)
        u->uP(NVAL-2, i, j, k) = (u->uC(NVAL-2, i, j, k)+u->uC(NVAL-2, i, j+1, k))/2.0;
    else
        u->uP(NVAL-2, i, j, k) = u->uC(NVAL-2, i, j, k);

    if (u->Nz > 1)
        u->uP(NVAL-1, i, j, k) = (u->uC(NVAL-1, i, j, k)+u->uC(NVAL-1, i, j, k+1))/2.0;
    else
        u->uP(NVAL-1, i, j, k) = u->uC(NVAL-1, i, j, k);

    return true;
}


/*  uP to uC  */
bool toConserved(Arrays *u, int i, int j, int k, ld gamma) {
    ld v2 = SQ(u->uP(1, i, j, k))+SQ(u->uP(2, i, j, k))+SQ(u->uP(3, i, j, k));
    ld B2 = SQ(u->uP(5, i, j, k))+SQ(u->uP(6, i, j, k))+SQ(u->uP(7, i, j, k));

    ld rh = u->uP(0, i, j, k); u->uC(0, i, j, k) = rh;
    u->uC(1, i, j, k) = rh*u->uP(1, i, j, k); u->uC(2, i, j, k) = rh*u->uP(2, i, j, k); u->uC(3, i, j, k) = rh*u->uP(3, i, j, k);
    
    u->uC(4, i, j, k) = u->uP(4, i, j, k)/(gamma-1.0) + rh*v2/2.0 + B2/2.0;

    if (isnan(u->uC(4, i, j, k))) return false;
    return true;
}

bool toConserved(ld *w, ld *u, ld gamma) {
    ld v2 = SQ(w[1])+SQ(w[2])+SQ(w[3]);
    ld B2 = SQ(w[5])+SQ(w[6])+SQ(w[7]);

    ld rh = w[0]; u[0] = rh;
    u[1] = rh*w[1]; u[2] = rh*w[2]; u[3] = rh*u[3];
    u[5] = w[5]; u[6] = w[6]; u[7] = w[7];

    u[4] = w[4]/(gamma-1) + rh*v2/2.0 + B2/2.0;

    if (isnan(u[4]) || null(u[4])) {
        std::cout << "ERROR:::AdiabMHD.h::toConserved:: E is NaN or null/negative: " << u[4] << std::endl;
        return false;
    }

}

/*  uC to uP  */
bool toPrimitive(Arrays *u, int i, int j, int k, ld gamma) {

    ld rh = u->uC(0, i, j, k); u->uP(0, i, j, k) = rh;
    u->uP(1, i, j, k) = u->uC(1, i, j, k)/rh; u->uP(2, i, j, k) = u->uC(2, i, j, k)/rh; u->uP(3, i, j, k) = u->uC(3, i, j, k)/rh;
    //u->uP(5, i, j, k) = u->uC(5, i, j, k); u->uP(6, i, j, k) = u->uC(6, i, j, k); u->uP(7, i, j, k) = u->uC(7, i, j, k);

    ld v2 = SQ(u->uP(1, i, j, k))+SQ(u->uP(2, i, j, k))+SQ(u->uP(3, i, j, k));
    ld B2 = SQ(u->uP(5, i, j, k))+SQ(u->uP(6, i, j, k))+SQ(u->uP(7, i, j, k));

    u->uP(4, i, j, k) = ( u->uC(4, i, j, k) - rh*v2/2.0 - B2/2.0 ) * (gamma-1.0);

    if (isnan(u->uP(4, i, j, k)) || u->uP(4, i, j, k) <= 0.0) {
        std::cout << "\ttoPrimitive()::P at (" << i << ", " << j << ", " << k << ") is NaN or null/negative: " << u->uP(4, i, j, k) << std::endl;
        std::cout << "\t\trho = " << u->uP(0, i, j, k) << ",  E = " << u->uC(4, i, j, k) << ",  v2 = " << v2 << ",  B2 = " << B2 << std::endl;
        std::cout << "\t\tBx at interfaces: " << u->uC(5, i, j, k) << " and " << u->uC(5, i+1, j, k) << std::endl;
        std::cout << "\t\tBy at interfaces: " << u->uC(6, i, j, k) << " and " << u->uC(6, i, j+1, k) << std::endl;
        std::cout << "\t\tBz at interfaces: " << u->uC(7, i, j, k) << " and " << u->uC(7, i, j, k+1) << std::endl;
        return false;
    }
    return true;
}

bool toPrimitive(ld *u, ld *w, ld gamma) {
    if (null(u[0])) {
        std::cout << "ERROR:::AdiabMHD.h::toPrimitive:: the density is null !" << std::endl;
        return false;
    }
    ld rh = u[0]; w[0] = rh;
    w[1] = u[1]/rh; w[2] = u[2]/rh; w[3] = u[3]/rh;
    w[5] = u[5]; w[6] = u[6]; w[7] = u[7];

    ld v2 = SQ(w[1])+SQ(w[2])+SQ(w[3]);
    ld B2 = SQ(w[5])+SQ(w[6])+SQ(w[7]);

    w[4] = (u[4] - rh*v2/2.0 - B2/2.0)*(gamma-1.0);

    if (isnan(w[4]) || null(w[4])) {
        std::cout << "ERROR:::AdiabMHD.h::toPrimitive:: P is NaN or null/negative: " << w[4] << std::endl;
        std::cout << "--> rho = " << rh << ", E: " << u[4] << ", v2: " << v2 << ", B2: " << B2 << std::endl;
        
        return false;
    }

    return true;

}

/* Transform the (i, j, k) interface values from primitive to conserved variables (used at the end of the reconstruction) */ 
bool toConserved(Arrays *u, int dim, int i, int j, int k, ld gamma) {
    if (dim == 0) {
        ld v2 = SQ(u->ix(1, i, j, k))+SQ(u->ix(2, i, j, k))+SQ(u->ix(3, i, j, k));
        ld B2 = SQ(u->ix(5, i, j, k))+SQ(u->ix(6, i, j, k))+SQ(u->ix(7, i, j, k));

        ld rh = u->ix(0, i, j, k); u->ix(0, i, j, k) = rh;
        u->ix(1, i, j, k) = rh*u->ix(1, i, j, k); u->ix(2, i, j, k) = rh*u->ix(2, i, j, k); u->ix(3, i, j, k) = rh*u->ix(3, i, j, k);

        u->ix(4, i, j, k) = u->ix(4, i, j, k)/(gamma-1) + rh*v2/2.0 + B2/2.0;

        if (isnan(u->ix(4, i, j, k))) return false;
        return true;
    }
    if (dim == 1) {
        ld v2 = SQ(u->iy(1, i, j, k))+SQ(u->iy(2, i, j, k))+SQ(u->iy(3, i, j, k));
        ld B2 = SQ(u->iy(5, i, j, k))+SQ(u->iy(6, i, j, k))+SQ(u->iy(7, i, j, k));

        ld rh = u->iy(0, i, j, k); u->iy(0, i, j, k) = rh;
        u->iy(1, i, j, k) = rh*u->iy(1, i, j, k); u->iy(2, i, j, k) = rh*u->iy(2, i, j, k); u->iy(3, i, j, k) = rh*u->iy(3, i, j, k);

        u->iy(4, i, j, k) = u->iy(4, i, j, k)/(gamma-1) + rh*v2/2.0 + B2/2.0;

        if (isnan(u->iy(4, i, j, k))) return false;
        return true;
    }
    if (dim == 2) {
        ld v2 = SQ(u->iz(1, i, j, k))+SQ(u->iz(2, i, j, k))+SQ(u->iz(3, i, j, k));
        ld B2 = SQ(u->iz(5, i, j, k))+SQ(u->iz(6, i, j, k))+SQ(u->iz(7, i, j, k));

        ld rh = u->iz(0, i, j, k); u->iz(0, i, j, k) = rh;
        u->iz(1, i, j, k) = rh*u->iz(1, i, j, k); u->iz(2, i, j, k) = rh*u->iz(2, i, j, k); u->iz(3, i, j, k) = rh*u->iz(3, i, j, k);

        u->iz(4, i, j, k) = u->iz(4, i, j, k)/(gamma-1) + rh*v2/2.0 + B2/2.0;

        if (isnan(u->iz(4, i, j, k))) return false;
        return true;
    }
    
}



bool F(ld uP[8], ld uC[8], ld *flux) {
    ld v = sqrtl(SQ(uP[1])+SQ(uP[2])+SQ(uP[3]));
    ld B = sqrtl(SQ(uP[5])+SQ(uP[6])+SQ(uP[7]));
    ld Pstar = uP[4] + B*B/2.0;

    flux[0] = uC[1];
    flux[1] = uC[1]*uP[1] + Pstar - uP[5]*uP[5];
    flux[2] = uC[1]*uP[2] - uP[5]*uP[6];
    flux[3] = uC[1]*uP[3] - uP[5]*uP[7];
    flux[4] = (uC[4]+Pstar)*uP[1] - (uP[5]*uP[1]+uP[6]*uP[2]+uP[7]*uP[3])*uP[5];
    flux[5] = 0.0;
    flux[6] = uP[6]*uP[1] - uP[5]*uP[2];
    flux[7] = uP[7]*uP[1] - uP[5]*uP[3];

    return true;
}

bool G(ld uP[8], ld uC[8], ld *flux) {
    ld v = sqrtl(SQ(uP[1])+SQ(uP[2])+SQ(uP[3]));
    ld B = sqrtl(SQ(uP[5])+SQ(uP[6])+SQ(uP[7]));
    ld Pstar = uP[4] + B*B/2.0;

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

bool H(ld uP[8], ld uC[8], ld *flux) {
    ld v = sqrtl(SQ(uP[1])+SQ(uP[2])+SQ(uP[3]));
    ld B = sqrtl(SQ(uP[5])+SQ(uP[6])+SQ(uP[7]));
    ld Pstar = uP[4] + B*B/2.0;

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


ld getDiv(Arrays *u, int i, int j, int k) {
    ld Div = 0;
    if (u->Nx>1) Div += (u->uC(5, i+1, j, k)-u->uC(5, i, j, k))/u->dx;
    if (u->Ny>1) Div += (u->uC(6, i, j+1, k)-u->uC(6, i, j, k))/u->dy;
    if (u->Nz>1) Div += (u->uC(7, i, j, k+1)-u->uC(7, i, j, k))/u->dz;
    return Div;
}


void getWavespeeds(Arrays *u, int i, int j, int k, ld gamma, ld *wvx, ld *wvy, ld *wvz) {
    int iBx = NVAL-3; int iBy = NVAL-2; int iBz = NVAL-1;
    ld v = sqrtl(SQ(u->uP(1, i, j, k))+SQ(u->uP(2, i, j, k))+SQ(u->uP(3, i, j, k)));
    ld B = sqrtl(SQ(u->uP(iBx, i, j, k))+SQ(u->uP(iBy, i, j, k))+SQ(u->uP(iBz, i, j, k)));
    

    ld rh = u->uC(0, i, j, k); if (rh <= 0) { cout << "ERROR:::AdiabMHD.h::getWavespeeds:: The density is null or negative !" << endl; return; }
    ld a = sqrtl(gamma*u->uP(4, i, j, k)/rh);
    ld CA = sqrtl(B*B/rh);

    if (isnan(a)) {
        cout << "ERROR:::AdiabMHD.h::getWavespeeds:: a is NaN" << endl;
        cout << "d = " << rh << ", P = " << u->uP(4, i, j, k) << endl;
        return;
    }
    

    /*wvx[0] = a;*/ wvx[1] = CA;
    ld CAx = sqrtl(SQ(u->uP(iBx, i, j, k))/rh);
    wvx[0] = sqrtl(0.5)*sqrtl(  (a*a+CA*CA) - sqrtl(SQ(a*a+CA*CA) - 4*a*a*CAx*CAx)  );// if(isnan(wvx[0])) cout << "wvx[0] is NaN" << endl;
    wvx[2] = sqrtl(0.5)*sqrtl(  (a*a+CA*CA) + sqrtl(SQ(a*a+CA*CA) - 4*a*a*CAx*CAx)  );// if(isnan(wvx[2])) cout << "wvx[2] is NaN" << endl;

    if (isnan(wvx[0])) {
        cout << "ERROR:::AdiabMHD.h::getWavespeeds(u):: wvx[0] is NaN" << endl;
        cout << "--> square: " << 0.5*((a*a+CA*CA) - sqrtl(SQ(a*a+CA*CA) - 4*a*a*CAx*CAx)) << endl;
        cout << "--> a2 = " << a*a << ", CA2 = " << CA*CA << ", CAx2 = " << CAx*CAx << endl;
        cout << "--> second square root: " << sqrtl(powl(a*a+CA*CA, 2) - 4*a*a*CAx*CAx) << endl;
        cout << "--> sqrtl(sq(a2)) = " << sqrtl(SQ(a*a)) << ", should be equal to a2 = " << a*a << endl;
        cout << "--> SQ(a2+CA2) = " << SQ(a*a+CA*CA) << " should be equal to a4 = " << SQ(a*a) << endl;
        
    }
    if (isnan(wvx[2])) {
        cout << "ERROR:::AdiabMHD.h::getWavespeeds(u):: wvx[2] is NaN" << endl;
        cout << "--> square: " << 0.5*((a*a+CA*CA) + sqrtl(SQ(a*a+CA*CA) - 4*a*a*CAx*CAx)) << endl;
        cout << "--> a2 = " << a*a << ", CA2 = " << CA*CA << ", CAx2 = " << CAx*CAx << endl;
        cout << "--> second square root: " << sqrtl(powl(a*a+CA*CA, 2) - 4*a*a*CAx*CAx) << endl;
        cout << "--> sqrtl(sq(a2)) = " << sqrtl(SQ(a*a)) << ", should be equal to a2 = " << a*a << endl;
        cout << "--> SQ(a2+CA2) = " << SQ(a*a+CA*CA) << " should be equal to a4 = " << SQ(a*a) << endl;
    }
    
    


    /*wvy[0] = a;*/ wvy[1] = CA;
    ld CAy = sqrtl(SQ(u->uP(iBy, i, j, k))/rh);
    wvy[0] = sqrtl(0.5)*sqrtl(  (a*a+CA*CA) - sqrtl(SQ(a*a+CA*CA) - 4*a*a*CAy*CAy)  );
    wvy[2] = sqrtl(0.5)*sqrtl(  (a*a+CA*CA) + sqrtl(SQ(a*a+CA*CA) - 4*a*a*CAy*CAy)  );

    /*wvz[0] = a;*/ wvz[1] = CA;
    ld CAz = sqrtl(SQ(u->uP(iBz, i, j, k))/rh);
    wvz[0] = sqrtl(0.5)*sqrtl(  (a*a+CA*CA) - sqrtl(SQ(a*a+CA*CA) - 4*a*a*CAz*CAz)  );
    wvz[2] = sqrtl(0.5)*sqrtl(  (a*a+CA*CA) + sqrtl(SQ(a*a+CA*CA) - 4*a*a*CAz*CAz)  );

}

void getWavespeeds(ld *w, ld gamma, ld *wvn) {
    int iBx = NVAL-3; int iBy = NVAL-2; int iBz = NVAL-1;
    /* Order: (rho, vn, v1, v2, P, Bn, B1, B2) */
    ld v = sqrtl(SQ(w[1])+SQ(w[2])+SQ(w[3]));
    ld B = sqrtl(SQ(w[iBx])+SQ(w[iBy])+SQ(w[iBz]));


    ld rh = w[0]; if (rh <= 0) { cout << "ERROR:::AdiabMHD.h::getWavespeeds:: The density is null or negative !" << endl; return; }
    ld a = sqrtl(gamma*w[4]/rh);
    ld CA = sqrtl(B*B/rh);


    /*wvn[0] = a;*/ wvn[1] = CA; if (isnan(wvn[1])) cout << "wvn[1] is NaN" << endl;
    ld CAn = sqrtl(SQ(w[5])/rh);
    wvn[0] = sqrtl(0.5)*sqrtl((a*a+CA*CA) - sqrtl(SQ(a*a+CA*CA) - 4*a*a*CAn*CAn));
    if (isnan(wvn[0])) {
        cout << "ERROR:::AdiabMHD.h::getWavespeeds(w):: wvn[0] is NaN" << endl;
        cout << "--> square: " << 0.5*((a*a+CA*CA) - sqrtl(SQ(a*a+CA*CA) - 4*a*a*CAn*CAn)) << endl;
        cout << "--> a2 = " << a*a << ", CA2 = " << CA*CA << ", CAn2 = " << CAn*CAn << endl;
        cout << "--> second square root: " << sqrtl(powl(a*a+CA*CA,2) - 4*a*a*CAn*CAn) << endl;
        cout << "--> sqrtl(sq(a2)) = " << sqrtl(SQ(a*a)) << ", should be equal to a2 = " << a*a << endl;
        cout << "--> SQ(a2+CA2) = " << SQ(a*a+CA*CA) << " should be equal to a4 = " << SQ(a*a) << endl;
    }
    wvn[2] = sqrtl(0.5)*sqrtl((a*a+CA*CA) + sqrtl(SQ(a*a+CA*CA) - 4*a*a*CAn*CAn));
    if (isnan(wvn[2])) {
        cout << "ERROR:::AdiabMHD.h::getWavespeeds(w):: wvn[2] is NaN" << endl;
        cout << "--> square: " << 0.5*((a*a+CA*CA) + sqrtl(SQ(a*a+CA*CA) - 4*a*a*CAn*CAn)) << endl;
        cout << "--> a2 = " << a*a << ", CA2 = " << CA*CA << ", CAn2 = " << CAn*CAn << endl;
        cout << "--> second square root: " << sqrtl(powl(a*a+CA*CA, 2) - 4*a*a*CAn*CAn) << endl;
        cout << "--> sqrtl(sq(a2)) = " << sqrtl(SQ(a*a)) << ", should be equal to a2 = " << a*a << endl;
        cout << "--> SQ(a2+CA2) = " << SQ(a*a+CA*CA) << " should be equal to a4 = " << SQ(a*a) << endl;
    }

}

ld maxWavespeed(ld *w, ld gamma) {
    int iBx = NVAL-3; int iBy = NVAL-2; int iBz = NVAL-1;

    /* Order: (rho, vn, v1, v2, P) */
    ld B = sqrtl(SQ(w[NVAL-3])+SQ(w[NVAL-2])+SQ(w[NVAL-1]));

    ld a = sqrtl(gamma*w[4]/w[0]);
    ld CA = sqrtl(B*B/w[0]);
    ld CAn = sqrtl(SQ(w[5])/w[0]);

    return sqrtl(0.5)*sqrtl((a*a+CA*CA) + sqrtl(SQ(a*a+CA*CA) - 4*a*a*CAn*CAn));
}


void minmaxRoeEigenvalue(ld *uL, ld *wL, ld *uR, ld *wR, ld gamma, ld *l) {
    ld rho = sqrtl(wL[0])*sqrtl(wR[0]);
    ld vx = (sqrtl(wL[0])*wL[1] + sqrtl(wR[0])*wR[1]) / (sqrtl(wL[0])+sqrtl(wR[0]));
    ld vy = (sqrtl(wL[0])*wL[2] + sqrtl(wR[0])*wR[2]) / (sqrtl(wL[0])+sqrtl(wR[0]));
    ld vz = (sqrtl(wL[0])*wL[3] + sqrtl(wR[0])*wR[3]) / (sqrtl(wL[0])+sqrtl(wR[0]));
    ld v2 = vx*vx + vy*vy + vz*vz;
    ld By = (sqrtl(wL[0])*wL[6] + sqrtl(wR[0])*wR[6]) / (sqrtl(wL[0])+sqrtl(wR[0]));
    ld Bz = (sqrtl(wL[0])*wL[7] + sqrtl(wR[0])*wR[7]) / (sqrtl(wL[0])+sqrtl(wR[0]));
    ld B2L = wL[5]*wL[5]+wL[6]*wL[6]+wL[7]*wL[7];   ld B2R = wR[5]*wR[5]+wR[6]*wR[6]+wR[7]*wR[7];
    ld B2 = wL[5]*wL[5]+By*By+Bz*Bz;

    ld HL = (uL[4]+wL[4]+0.5*B2L)/wL[0];   ld HR = (uR[4]+wR[4]+0.5*B2R)/wR[0];
    ld H = (sqrtl(wL[0])*HL + sqrtl(wR[0])*HR) / (sqrtl(wL[0])+sqrtl(wR[0]));
    ld X = ( SQ(wL[6]-wR[6]) + SQ(wL[7]-wR[7]) ) / ( 2*(sqrtl(wL[0])+sqrtl(wR[0])) );
    ld Y = (wL[0]+wR[0])/(2*rho);


    ld _gamma = gamma-1;   ld _X = (gamma-2)*X;   ld _Y = (gamma-2)*Y;

    ld bperp2 = (_gamma-_Y)*(By*By+Bz*Bz);
    ld CAn2 = wL[5]*wL[5]/rho;   ld CA2 = CAn2 + bperp2/rho;
    ld a2 = _gamma*(H-v2/2-B2/rho)-_X;
    
    ld Cf2 = 0.5*( (a2+CA2) + sqrtl(SQ(a2+CA2) - 4*a2*CAn2) );

    l[0] = vx - sqrtl(Cf2);
    l[1] = vx + sqrtl(Cf2);
}
