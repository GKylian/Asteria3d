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

    u->uC(5, i, j, k) = u->uP(5, i, j, k); u->uC(6, i, j, k) = u->uP(6, i, j, k); u->uC(7, i, j, k) = u->uP(7, i, j, k);
    long double rh = u->uP(0, i, j, k); u->uC(0, i, j, k) = rh;
    u->uC(1, i, j, k) = rh*u->uP(1, i, j, k); u->uC(2, i, j, k) = rh*u->uP(2, i, j, k); u->uC(3, i, j, k) = rh*u->uP(3, i, j, k);
    
    u->uC(4, i, j, k) = u->uP(4, i, j, k)/(gamma-1) + rh*v2/2.0 + B2/2.0;

    if (isnan(u->uC(4, i, j, k))) return false;
    return true;
}

bool toConserved(long double *w, long double *u, long double gamma) {
    long double v2 = SQ(w[1])+SQ(w[2])+SQ(w[3]);
    long double B2 = SQ(w[5])+SQ(w[6])+SQ(w[7]);

    long double rh = w[0]; u[0] = rh;
    u[1] = rh*w[1]; u[2] = rh*w[2]; u[3] = rh*u[3];
    u[5] = w[5]; u[6] = w[6]; u[7] = w[7];

    u[4] = w[4]/(gamma-1) + rh*v2/2.0 + B2/2.0;

    if (isnan(u[4]) || null(u[4])) {
        std::cout << "ERROR:::AdiabMHD.h::toConserved:: E is NaN or null/negative: " << u[4] << std::endl;
        return false;
    }

}

/*  uC to uP  */
bool toPrimitive(Arrays *u, int i, int j, int k, long double gamma) {
    long double rh = u->uC(0, i, j, k); u->uP(0, i, j, k) = rh;
    u->uP(1, i, j, k) = u->uC(1, i, j, k)/rh; u->uP(2, i, j, k) = u->uC(2, i, j, k)/rh; u->uP(3, i, j, k) = u->uC(3, i, j, k)/rh;
    u->uP(5, i, j, k) = u->uC(5, i, j, k); u->uP(6, i, j, k) = u->uC(6, i, j, k); u->uP(7, i, j, k) = u->uC(7, i, j, k);

    long double v2 = SQ(u->uP(1, i, j, k))+SQ(u->uP(2, i, j, k))+SQ(u->uP(3, i, j, k));
    long double B2 = SQ(u->uP(5, i, j, k))+SQ(u->uP(6, i, j, k))+SQ(u->uP(7, i, j, k));

    u->uP(4, i, j, k) = ( u->uC(4, i, j, k) - rh*v2/2.0 - B2/2.0 ) * (gamma-1);

    if (isnan(u->uP(4, i, j, k))) return false;
    return true;
}

bool toPrimitive(long double *u, long double *w, long double gamma) {
    if (null(u[0])) {
        std::cout << "ERROR:::AdiabMHD.h::toPrimitive:: the density is null !" << std::endl;
        return false;
    }
    long double rh = u[0]; w[0] = rh;
    w[1] = u[1]/rh; w[2] = u[2]/rh; w[3] = u[3]/rh;
    w[5] = u[5]; w[6] = u[6]; w[7] = u[7];

    long double v2 = SQ(w[1])+SQ(w[2])+SQ(w[3]);
    long double B2 = SQ(w[5])+SQ(w[6])+SQ(w[7]);

    w[4] = (u[4] - rh*v2/2.0 - B2/2.0)*(gamma-1);

    if (isnan(w[4]) || null(w[4])) {
        std::cout << "ERROR:::AdiabMHD.h::toPrimitive:: P is NaN or null/negative: " << w[4] << std::endl;
        std::cout << "--> rho = " << rh << ", E: " << u[4] << ", v2: " << v2 << ", B2: " << B2 << std::endl;
        
        return false;
    }

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
    flux[4] = (uC[4]+Pstar)*uP[1] - (uP[5]*uP[1]+uP[6]*uP[2]+uP[7]*uP[3])*uP[5];
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

#ifndef ISO
void getWavespeeds(Arrays *u, int i, int j, int k, long double gamma, long double *wvx, long double *wvy, long double *wvz) {
    int iBx = NVAL-3; int iBy = NVAL-2; int iBz = NVAL-1;
    long double v = sqrtl(SQ(u->uP(1, i, j, k))+SQ(u->uP(2, i, j, k))+SQ(u->uP(3, i, j, k)));
    long double B = sqrtl(SQ(u->uP(iBx, i, j, k))+SQ(u->uP(iBy, i, j, k))+SQ(u->uP(iBz, i, j, k)));
    

    long double rh = u->uC(0, i, j, k); if (rh <= 0) { cout << "ERROR:::AdiabMHD.h::getWavespeeds:: The density is null or negative !" << endl; return; }
    long double a = sqrtl(gamma*u->uP(4, i, j, k)/rh);
    long double CA = sqrtl(B*B/rh);

    if (isnan(a)) {
        cout << "ERROR:::AdiabMHD.h::getWavespeeds:: a is NaN" << endl;
        cout << "d = " << rh << ", P = " << u->uP(4, i, j, k) << endl;
        return;
    }
    

    /*wvx[0] = a;*/ wvx[1] = CA;
    long double CAx = sqrtl(SQ(u->uP(iBx, i, j, k))/rh);
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
    long double CAy = sqrtl(SQ(u->uP(iBy, i, j, k))/rh);
    wvy[0] = sqrtl(0.5)*sqrtl(  (a*a+CA*CA) - sqrtl(SQ(a*a+CA*CA) - 4*a*a*CAy*CAy)  );
    wvy[2] = sqrtl(0.5)*sqrtl(  (a*a+CA*CA) + sqrtl(SQ(a*a+CA*CA) - 4*a*a*CAy*CAy)  );

    /*wvz[0] = a;*/ wvz[1] = CA;
    long double CAz = sqrtl(SQ(u->uP(iBz, i, j, k))/rh);
    wvz[0] = sqrtl(0.5)*sqrtl(  (a*a+CA*CA) - sqrtl(SQ(a*a+CA*CA) - 4*a*a*CAz*CAz)  );
    wvz[2] = sqrtl(0.5)*sqrtl(  (a*a+CA*CA) + sqrtl(SQ(a*a+CA*CA) - 4*a*a*CAz*CAz)  );

}

void getWavespeeds(long double *w, long double gamma, long double *wvn) {
    int iBx = NVAL-3; int iBy = NVAL-2; int iBz = NVAL-1;
    /* Order: (rho, vn, v1, v2, P, Bn, B1, B2) */
    long double v = sqrtl(SQ(w[1])+SQ(w[2])+SQ(w[3]));
    long double B = sqrtl(SQ(w[iBx])+SQ(w[iBy])+SQ(w[iBz]));


    long double rh = w[0]; if (rh <= 0) { cout << "ERROR:::AdiabMHD.h::getWavespeeds:: The density is null or negative !" << endl; return; }
    long double a = sqrtl(gamma*w[4]/rh);
    long double CA = sqrtl(B*B/rh);


    /*wvn[0] = a;*/ wvn[1] = CA; if (isnan(wvn[1])) cout << "wvn[1] is NaN" << endl;
    long double CAn = sqrtl(SQ(w[5])/rh);
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
#endif // !ISO


