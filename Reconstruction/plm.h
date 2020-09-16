#pragma once
#include "../Defs/Arrays.h"
#include "eigen.h"

/* Computes the left and right interfaces of the (i, j, k) cell using PLM reconstruction
   (with characteristic variables characteritic tracing)                                 */


#ifdef MHD
bool PLM(Arrays *u, int i, int j, int k) {

    /* ----- 0. Initialize all the values, arrays and matrices ----- */

    /* Arrays holding the primitive values of the cell and its neighbors */
    long double wi[NVAL] = { 0 }; long double wip[NVAL] = { 0 }; long double wim[NVAL] = { 0 };
    long double wj[NVAL] = { 0 }; long double wjp[NVAL] = { 0 }; long double wjm[NVAL] = { 0 };
    long double wk[NVAL] = { 0 }; long double wkp[NVAL] = { 0 }; long double wkm[NVAL] = { 0 };

    /* Arrays and matrices holding the eigenvalues and the left/right eigenvectors */
    long double eigen_x[NWAVE] = { 0 }; long double eigen_y[NWAVE] = { 0 }; long double eigen_z[NWAVE] = { 0 };
    long double Lx[NWAVE][NWAVE] = { 0 }; long double Rx[NWAVE][NWAVE] = { 0 };
    long double Ly[NWAVE][NWAVE] = { 0 }; long double Ry[NWAVE][NWAVE] = { 0 };
    long double Lz[NWAVE][NWAVE] = { 0 }; long double Rz[NWAVE][NWAVE] = { 0 };

    /* The left, center and right differences in primitive variables */
    long double dwL_x[NWAVE] = { 0 }; long double dwR_x[NWAVE] = { 0 }; long double dwC_x[NWAVE] = { 0 };
    long double dwL_y[NWAVE] = { 0 }; long double dwR_y[NWAVE] = { 0 }; long double dwC_y[NWAVE] = { 0 };
    long double dwL_z[NWAVE] = { 0 }; long double dwR_z[NWAVE] = { 0 }; long double dwC_z[NWAVE] = { 0 };

    /* The differences projected onto characteristic variables */
    long double daL_x[NWAVE] = { 0 }; long double daR_x[NWAVE] = { 0 }; long double daC_x[NWAVE] = { 0 };
    long double daL_y[NWAVE] = { 0 }; long double daR_y[NWAVE] = { 0 }; long double daC_y[NWAVE] = { 0 };
    long double daL_z[NWAVE] = { 0 }; long double daR_z[NWAVE] = { 0 }; long double daC_z[NWAVE] = { 0 };

    /* The differences with monotonicity constraints applied */
    long double da_x[NWAVE] = { 0 }; long double da_y[NWAVE] = { 0 }; long double da_z[NWAVE] = { 0 };

    /* The differences projected back onto primitive variables */
    long double dw_x[NWAVE] = { 0 }; long double dw_y[NWAVE] = { 0 }; long double dw_z[NWAVE] = { 0 };



    /* ----- 1. Load the primitive values to the local arrays ----- */
    /* The order should be (rho, vn, v1, v2, P, Bn, B1, B2) */
    for (int n = 0; n < NVAL; n++) {
        wi[n] = u->uP(n, i, j, k);    wip[n] = u->uP(n, i+1, j, k);    wim[n] = u->uP(n, i-1, j, k);
        wj[n] = u->uP(n, i, j, k);    wjp[n] = u->uP(n, i, j+1, k);    wjm[n] = u->uP(n, i, j-1, k);
        wk[n] = u->uP(n, i, j, k);    wkp[n] = u->uP(n, i, j, k+1);    wkm[n] = u->uP(n, i, j, k-1);
    }
    /* Swap v and B to get the correct order */
    /* y: (vx, vy, vz) -> (vy, vx, vz)       z: (vx, vy, vz) -> (vz, vx, vy) */
    swap(wj[1], wj[2]);    swap(wjp[1], wjp[2]);    swap(wjm[1], wjm[2]);
    swap(wj[NVAL-3], wj[NVAL-2]);    swap(wjp[NVAL-3], wjp[NVAL-2]);    swap(wjm[NVAL-3], wjm[NVAL-2]);

    swap(wk[1], wk[3]); swap(wk[2], wk[3]);    swap(wkp[1], wkp[3]); swap(wkp[2], wkp[3]);    swap(wkm[1], wkm[3]); swap(wkm[2], wkm[3]);
    swap(wk[NVAL-3], wk[NVAL-1]); swap(wk[NVAL-2], wk[NVAL-1]);
    swap(wkp[NVAL-3], wkp[NVAL-1]); swap(wkp[NVAL-2], wkp[NVAL-1]);
    swap(wkm[NVAL-3], wkm[NVAL-1]); swap(wkm[NVAL-2], wkm[NVAL-1]);


    /* ----- 2. Compute the eigenvectors and eigenvalues ----- */
    if (!getEigen(wi, eigen_x, Lx, Rx)) std::cout << "plm.h::PLM:: Could not compute eigenvalues and eigenvectors in the 'x' direction." << std::endl;
    if (!getEigen(wj, eigen_y, Ly, Ry)) std::cout << "plm.h::PLM:: Could not compute eigenvalues and eigenvectors in the 'y' direction." << std::endl;
    if (!getEigen(wk, eigen_z, Lz, Rz)) std::cout << "plm.h::PLM:: Could not compute eigenvalues and eigenvectors in the 'z' direction." << std::endl;

}

#else
bool PLM(Arrays *u, int i, int j, int k) {

}

#endif // MHD