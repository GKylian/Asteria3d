#pragma once
#include "../Defs/Arrays.h"
#include "eigen.h"
#include "../prob.h"

/* Computes the left and right interfaces of the (i, j, k) cell using PLM reconstruction
   (with characteristic variables characteritic tracing)                                 */


#ifdef MHD
bool PLM(Arrays *u, int i, int j, int k) {
    int nx = u->Nx; int ny = u->Ny; int nz = u->Nz;
    int iR = 2*i; int iL = 2*i+1; int jR = 2*j; int jL = 2*j+1; int kR = 2*k; int kL = 2*k+1; //interface indices of the right (i-1/2) and left (i+1/2) interfaces (in the cell i, j, k)
    int iBx = NVAL-3; int iBy = NVAL-2; int iBz = NVAL-1;

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
    /* The order should be (rho, vn, v1, v2, P, B1, B2, Bn) */
    for (int n = 0; n < NVAL; n++) {
        wi[n] = u->uP(n, i, j, k);    wip[n] = u->uP(n, i+1, j, k);    wim[n] = u->uP(n, i-1, j, k);
        wj[n] = u->uP(n, i, j, k);    wjp[n] = u->uP(n, i, j+1, k);    wjm[n] = u->uP(n, i, j-1, k);
        wk[n] = u->uP(n, i, j, k);    wkp[n] = u->uP(n, i, j, k+1);    wkm[n] = u->uP(n, i, j, k-1);
    }
    /* Swap v and B to get the correct order */
    /* y: (vx, vy, vz) -> (vy, vx, vz)       z: (vx, vy, vz) -> (vz, vx, vy) */
    /* x: (Bx, By, Bz) -> (By, Bz, Bx)       y: (Bx, By, Bz) -> (Bx, Bz, By) */
    swap( wi[iBx],  wi[iBz]); swap( wi[iBx],  wi[iBy]);
    swap(wip[iBx], wip[iBz]); swap(wip[iBx], wip[iBy]);
    swap(wim[iBx], wim[iBz]); swap(wim[iBx], wim[iBy]);

    swap(wj[1], wj[2]);    swap(wjp[1], wjp[2]);    swap(wjm[1], wjm[2]);
    swap(wj[iBy], wj[iBz]);    swap(wjp[iBy], wjp[iBz]);    swap(wjm[iBy], wjm[iBz]);

    swap(wk[1], wk[3]); swap(wk[2], wk[3]);    swap(wkp[1], wkp[3]); swap(wkp[2], wkp[3]);    swap(wkm[1], wkm[3]); swap(wkm[2], wkm[3]);
    


    /* ----- 2. Compute and store the eigenvectors and eigenvalues, as well as the min/max eigenvalues ----- */
    if (nx > 1) { if (!getEigen(wi, eigen_x, Lx, Rx)) std::cout << "plm.h::PLM:: Could not compute eigenvalues and eigenvectors in the 'x' direction." << std::endl; return false; }
    if (ny > 1) { if (!getEigen(wj, eigen_y, Ly, Ry)) std::cout << "plm.h::PLM:: Could not compute eigenvalues and eigenvectors in the 'y' direction." << std::endl; return false; }
    if (nz > 1) { if (!getEigen(wk, eigen_z, Lz, Rz)) std::cout << "plm.h::PLM:: Could not compute eigenvalues and eigenvectors in the 'z' direction." << std::endl; return false; }

    long double eigM_x = (long double)*max_element(eigen_x, eigen_x+NWAVE); long double eig0_x = (long double)*min_element(eigen_x, eigen_x+NWAVE);
    long double eigM_y = (long double)*max_element(eigen_y, eigen_y+NWAVE); long double eig0_y = (long double)*min_element(eigen_y, eigen_y+NWAVE);
    long double eigM_z = (long double)*max_element(eigen_z, eigen_z+NWAVE); long double eig0_z = (long double)*min_element(eigen_z, eigen_z+NWAVE);




    /* ----- 3. Compute the left, right and centered differences in primitive variables. ----- */
    if (nx > 1) {    for (int n = 0; n < NWAVE; n++) { dwL_x[n] = wi[n]-wim[n];   dwR_x[n] = wip[n]-wi[n];   dwC_x[n] = (wip[n]-wim[n])/2.0;}    }
    if (ny > 1) {    for (int n = 0; n < NWAVE; n++) { dwL_y[n] = wj[n]-wjm[n];   dwR_y[n] = wjp[n]-wj[n];   dwC_y[n] = (wjp[n]-wjm[n])/2.0;}    }
    if (nz > 1) {    for (int n = 0; n < NWAVE; n++) { dwL_z[n] = wk[n]-wkm[n];   dwR_z[n] = wkp[n]-wk[n];   dwC_z[n] = (wkp[n]-wkm[n])/2.0;}    }



    /* ----- 4. Project the differences along the characteristics. ----- */
    if (nx > 1) {    for (int n = 0; n < NWAVE; n++) for (int m = 0; m < NWAVE; m++){
                         daL_x[n] += Lx[n][m]*dwL_x[m]; daR_x[n] += Lx[n][m]*dwR_x[m]; daC_x[n] += Lx[n][m]*dwC_x[m]; }
    }
    if (ny > 1) {    for (int n = 0; n < NWAVE; n++) for (int m = 0; m < NWAVE; m++){
                         daL_y[n] += Ly[n][m]*dwL_y[m]; daR_y[n] += Ly[n][m]*dwR_y[m]; daC_y[n] += Ly[n][m]*dwC_y[m]; }
    }
    if (nz > 1) {    for (int n = 0; n < NWAVE; n++) for (int m = 0; m < NWAVE; m++){
                         daL_z[n] += Lz[n][m]*dwL_z[m]; daR_z[n] += Lz[n][m]*dwR_z[m]; daC_z[n] += Lz[n][m]*dwC_z[m]; }
    }
    


    /* ----- 5. Apply the monotonicity constraints. ----- */
    if (nx > 1) {    for (int n = 0; n < NWAVE; n++) { da_x[n] = sgn(daC_x[n])*fminl(fminl(2.0*fabsl(daL_x[n]), 2.0*fabsl(daR_x[n])), fabsl(daC_x[n]));  if (null(da_x[n])) da_x[n] = 0.0; }    }
    if (ny > 1) {    for (int n = 0; n < NWAVE; n++) { da_y[n] = sgn(daC_y[n])*fminl(fminl(2.0*fabsl(daL_y[n]), 2.0*fabsl(daR_y[n])), fabsl(daC_y[n]));  if (null(da_y[n])) da_y[n] = 0.0; }    }
    if (nz > 1) {    for (int n = 0; n < NWAVE; n++) { da_z[n] = sgn(daC_z[n])*fminl(fminl(2.0*fabsl(daL_z[n]), 2.0*fabsl(daR_z[n])), fabsl(daC_z[n]));  if (null(da_z[n])) da_z[n] = 0.0; }    }


    /* ----- 6. Project back to the primitive variabkes. ----- */
    if (nx > 1){   for (int n = 0; n < NWAVE; n++) for (int m = 0; m < NWAVE; m++){
                        dw_x[n] += da_x[m]*Rx[n][m]; if (null(dw_x[n])) dw_x[n] = 0.0;}
    }
    if (ny > 1){   for (int n = 0; n < NWAVE; n++) for (int m = 0; m < NWAVE; m++){
                        dw_y[n] += da_y[m]*Ry[n][m]; if (null(dw_y[n])) dw_y[n] = 0.0;}
    }
    if (nz > 1){   for (int n = 0; n < NWAVE; n++) for (int m = 0; m < NWAVE; m++){
                        dw_z[n] += da_z[m]*Rz[n][m]; if (null(dw_z[n])) dw_z[n] = 0.0;}
    }


    
    /* ----- 7. Compute the values at R(i-1/2) and L(i+1/2). ----- */
    if (nx > 1) { for(int n = 0; n < NWAVE; n++){
        u->ix(n, iR, j, k) = wi[n] - 0.5*dw_x[n];    u->ix(n, iL, j, k) = wi[n] + 0.5*dw_x[n];
        u->ix(n, iR, j, k) = fmaxl(fminl(wi[n], wim[n]), u->ix(n, iR, j, k)); u->ix(n, iR, j, k) = fminl(fmaxl(wi[n], wim[n]), u->ix(n, iR, j, k));
        u->ix(n, iL, j, k) = fmaxl(fminl(wi[n], wip[n]), u->ix(n, iL, j, k)); u->ix(n, iL, j, k) = fminl(fmaxl(wi[n], wip[n]), u->ix(n, iL, j, k));
        u->ix(n, iR, j, k) += -u->dt/(2.0*u->dx)*fminl(eig0_x, 0.0)*dw_x[n];    u->ix(n, iL, j, k) = -u->dt/(2.0*u->dx)*fmaxl(eigM_x, 0.0)*dw_x[n];
        //if (null(u->ix(n, iR, j, k))) u->ix(n, iR, j, k) = 0.0;   if (u->ix(n, iL, j, k)) u->ix(n, iL, j, k) = 0.0;
                  }
    }
    if (ny > 1) { for(int n = 0; n < NWAVE; n++){
        u->iy(n, i, jR, k) = wj[n] - 0.5*dw_y[n];    u->iy(n, i, jL, k) = wj[n] + 0.5*dw_y[n];
        u->iy(n, i, jR, k) = fmaxl(fminl(wj[n], wjm[n]), u->iy(n, i, jR, k)); u->iy(n, i, jR, k) = fminl(fmaxl(wj[n], wjm[n]), u->iy(n, i, jR, k));
        u->iy(n, i, jL, k) = fmaxl(fminl(wj[n], wjp[n]), u->iy(n, i, jL, k)); u->iy(n, i, jL, k) = fminl(fmaxl(wj[n], wjp[n]), u->iy(n, i, jL, k));
        u->iy(n, i, jR, k) += -u->dt/(2.0*u->dy)*fminl(eig0_y, 0.0)*dw_y[n];    u->iy(n, i, jL, k) = -u->dt/(2.0*u->dy)*fmaxl(eigM_y, 0.0)*dw_y[n];
        //if (null(u->iy(n, i, jR, k))) u->iy(n, i, jR, k) = 0.0;   if (u->iy(n, i, jL, k)) u->iy(n, i, jL, k) = 0.0;
                  }
    }
    if (nz > 1) { for(int n = 0; n < NWAVE; n++){
        u->iz(n, i, j, kR) = wk[n] - 0.5*dw_z[n];    u->iz(n, i, j, kL) = wk[n] + 0.5*dw_z[n];
        u->iz(n, i, j, kR) = fmaxl(fminl(wk[n], wkm[n]), u->iz(n, i, j, kR)); u->iz(n, i, j, kR) = fminl(fmaxl(wk[n], wkm[n]), u->iz(n, i, j, kR));
        u->iz(n, i, j, kL) = fmaxl(fminl(wk[n], wkp[n]), u->iz(n, i, j, kL)); u->iz(n, i, j, kL) = fminl(fmaxl(wk[n], wkp[n]), u->iz(n, i, j, kL));
        u->iz(n, i, j, kR) += -u->dt/(2.0*u->dz)*fminl(eig0_z, 0.0)*dw_z[n];    u->iz(n, i, j, kL) = -u->dt/(2.0*u->dz)*fmaxl(eigM_z, 0.0)*dw_z[n];
        //if (null(u->iz(n, i, j, kR))) u->iz(n, i, j, kR) = 0.0;   if (u->iz(n, i, j, kL)) u->iz(n, i, j, kL) = 0.0;
                  }
    }



    /* ----- 8. Perform the characteristic tracing on all dimensions where it is needed ----- */

#ifdef TRACING

    if (nx > 1) {
        long double qa = 0.0; long double qn;
        for (int n = 0; n < 7; n++)
        {
            if (eigen_x[n] >= 0.0) {
                //Substract amount of each wave that doesn't reach the interface
                qa = 0.0; qn = 0.5*u->dt/u->dx*(eigM_x-eigen_x[n]);
                for (int m = 0; m < 7; m++)
                    qa += Lx[n][m]*qn*dw_x[m];
                if (null(qa)) qa = 0.0;
                for (int m = 0; m < 7; m++)
                    u->ix(m, iL, j, k) += qa*Rx[m][n];


                //For HLL, substract those moving AWAY from the interface
                qa = 0.0; qn = 0.5*u->dt/u->dx*(eig0_x-eigen_x[n]);
                for (int m = 0; m < 7; m++)
                    qa += Lx[n][m]*qn*dw_x[m];
                if (null(qa)) qa = 0.0;
                for (int m = 0; m < 7; m++)
                    u->ix(m, iR, j, k) += qa*Rx[m][n];
            }

            if (eigen_x[n] <= 0.0) {
                //Substract amount of each wave that doesn't reach the interface
                qa = 0.0; qn = 0.5*u->dt/u->dx*(eig0_x-eigen_x[n]);
                for (int m = 0; m < 7; m++)
                    qa += Lx[n][m]*qn*dw_x[m];
                if (null(qa)) qa = 0.0;
                for (int m = 0; m < 7; m++)
                    u->ix(m, iR, j, k) += qa*Rx[m][n];

                //For HLL, substract those moving AWAY from the interface
                qa = 0.0; qn = 0.5*u->dt/u->dx*(eigM_x-eigen_x[n]);
                for (int m = 0; m < 7; m++)
                    qa += Lx[n][m]*qn*dw_x[m];
                if (null(qa)) qa = 0.0;
                for (int m = 0; m < 7; m++)
                    u->ix(m, iL, j, k) += qa*Rx[m][n];
            }
        }
    }
    if (ny > 1) {
        long double qa = 0.0; long double qn;
        for (int n = 0; n < 7; n++)
        {
            if (eigen_y[n] >= 0.0) {
                //Substract amount of each wave that doesn't reach the interface
                qa = 0.0; qn = 0.5*u->dt/u->dy*(eigM_y-eigen_y[n]);
                for (int m = 0; m < 7; m++)
                    qa += Ly[n][m]*qn*dw_y[m];
                if (null(qa)) qa = 0.0;
                for (int m = 0; m < 7; m++)
                    u->iy(m, i, jL, k) += qa*Ry[m][n];


                //For HLL, substract those moving AWAY from the interface
                qa = 0.0; qn = 0.5*u->dt/u->dy*(eig0_y-eigen_y[n]);
                for (int m = 0; m < 7; m++)
                    qa += Ly[n][m]*qn*dw_y[m];
                if (null(qa)) qa = 0.0;
                for (int m = 0; m < 7; m++)
                    u->iy(m, i, jR, k) += qa*Ry[m][n];
            }

            if (eigen_y[n] <= 0.0) {
                //Substract amount of each wave that doesn't reach the interface
                qa = 0.0; qn = 0.5*u->dt/u->dy*(eig0_y-eigen_y[n]);
                for (int m = 0; m < 7; m++)
                    qa += Ly[n][m]*qn*dw_y[m];
                if (null(qa)) qa = 0.0;
                for (int m = 0; m < 7; m++)
                    u->iy(m, i, jR, k) += qa*Ry[m][n];

                //For HLL, substract those moving AWAY from the interface
                qa = 0.0; qn = 0.5*u->dt/u->dy*(eigM_y-eigen_y[n]);
                for (int m = 0; m < 7; m++)
                    qa += Ly[n][m]*qn*dw_y[m];
                if (null(qa)) qa = 0.0;
                for (int m = 0; m < 7; m++)
                    u->iy(m, i, jL, k) += qa*Ry[m][n];
            }
        }
    }
    if (nz > 1) {
        long double qa = 0.0; long double qn;
        for (int n = 0; n < 7; n++)
        {
            if (eigen_z[n] >= 0.0) {
                //Substract amount of each wave that doesn't reach the interface
                qa = 0.0; qn = 0.5*u->dt/u->dz*(eigM_z-eigen_z[n]);
                for (int m = 0; m < 7; m++)
                    qa += Lz[n][m]*qn*dw_z[m];
                if (null(qa)) qa = 0.0;
                for (int m = 0; m < 7; m++)
                    u->iz(m, i, j, kL) += qa*Rz[m][n];


                //For HLL, substract those moving AWAY from the interface
                qa = 0.0; qn = 0.5*u->dt/u->dz*(eig0_z-eigen_z[n]);
                for (int m = 0; m < 7; m++)
                    qa += Lz[n][m]*qn*dw_z[m];
                if (null(qa)) qa = 0.0;
                for (int m = 0; m < 7; m++)
                    u->iz(m, i, j, kR) += qa*Rz[m][n];
            }

            if (eigen_z[n] <= 0.0) {
                //Substract amount of each wave that doesn't reach the interface
                qa = 0.0; qn = 0.5*u->dt/u->dz*(eig0_z-eigen_z[n]);
                for (int m = 0; m < 7; m++)
                    qa += Lz[n][m]*qn*dw_z[m];
                if (null(qa)) qa = 0.0;
                for (int m = 0; m < 7; m++)
                    u->iz(m, i, j, kR) += qa*Rz[m][n];

                //For HLL, substract those moving AWAY from the interface
                qa = 0.0; qn = 0.5*u->dt/u->dz*(eigM_z-eigen_z[n]);
                for (int m = 0; m < 7; m++)
                    qa += Lz[n][m]*qn*dw_z[m];
                if (null(qa)) qa = 0.0;
                for (int m = 0; m < 7; m++)
                    u->iz(m, i, j, kL) += qa*Rz[m][n];
            }
        }
    }
#endif // TRACING

    /* ----- 8.2 Switch back elements to normal order ----- */
    /* y: (vy, vx, vz) -> (vx, vy, vz)       z: (vz, vx, vy) -> (vx, vy, vz) */
    /* x: (By, Bz, Bx) -> (Bx, By, Bz)       y: (Bx, Bz, By) -> (Bx, By, Bz) */
    swap(u->iy(1, i, jL, k), u->iy(2, i, jL, k));    swap(u->iy(1, i, jR, k), u->iy(2, i, jR, k));
    swap(u->iz(1, i, j, kL), u->iz(2, i, j, kL)); swap(u->iz(1, i, j, kR), u->iz(2, i, j, kR));    swap(u->iz(2, i, j, kL), u->iz(3, i, j, kL)); swap(u->iz(2, i, j, kR), u->iz(3, i, j, kR));
    
    swap(u->ix(iBx, iL, j, k), u->ix(iBy, iL, j, k)); swap(u->ix(iBx, iL, j, k), u->ix(iBz, iL, j, k));    swap(u->ix(iBx, iR, j, k), u->ix(iBy, iR, j, k)); swap(u->ix(iBx, iR, j, k), u->ix(iBz, iR, j, k));
    swap(u->iy(iBy, i, jL, k), u->iy(iBz, i, jL, k));    swap(u->iy(iBy, i, jR, k), u->iy(iBz, i, jR, k));



    /* ----- 9. Add magnetic field source terms ----- */
    long double dBx = (u->uC(iBx, i+1, j, k)-u->uC(iBx, i, j, k))/u->dx;
    long double dBy = (u->uC(iBy, i, j+1, k)-u->uC(iBy, i, j, k))/u->dy;
    long double dBz = (u->uC(iBz, i, j, k+1)-u->uC(iBz, i, j, k))/u->dz;

    long double dBy_x = -u->dt/2.0 * u->uP(2, i, j, k) * minmod(  (nz>1)*dBz, -(nx>1)*dBx  ); //Changed to By and Bz at the x interfaces
    long double dBz_x = -u->dt/2.0 * u->uP(3, i, j, k) * minmod(  (ny>1)*dBy, -(nx>1)*dBx  );

    long double dBz_y = -u->dt/2.0 * u->uP(3, i, j, k) * minmod(  (nx>1)*dBx, -(ny>1)*dBy  ); //Changed to Bx and Bz at the y interfaces
    long double dBx_y = -u->dt/2.0 * u->uP(1, i, j, k) * minmod(  (nz>1)*dBz, -(ny>1)*dBy  );

    long double dBx_z = -u->dt/2.0 * u->uP(1, i, j, k) * minmod(  (ny>1)*dBy, -(nz>1)*dBz  ); //Changed to Bx and Bz at the y interfaces
    long double dBy_z = -u->dt/2.0 * u->uP(2, i, j, k) * minmod(  (nx>1)*dBx, -(nz>1)*dBz  );


    u->ix(iBy, iR, j, k) += dBy_x;  u->ix(iBy, iL, j, k) += dBy_x;
    u->ix(iBz, iR, j, k) += dBz_x;  u->ix(iBz, iL, j, k) += dBz_x;

    u->iy(iBz, i, jR, k) += dBz_y;  u->iy(iBz, i, jL, k) += dBz_y;
    u->iy(iBx, i, jR, k) += dBx_y;  u->iy(iBx, i, jL, k) += dBx_y;

    u->iz(iBx, i, j, kR) += dBx_z;  u->iz(iBx, i, j, kL) += dBx_z;
    u->iz(iBy, i, j, kR) += dBy_z;  u->iz(iBy, i, j, kL) += dBy_z;



    /* ----- 10. Add gravity source terms (static or potential) ----- */
    if (gx != 0.0 || gy != 0.0 || gz != 0.0) {
        u->ix(1, iR, j, k) += u->dt/2.0 *gx; u->ix(1, iL, j, k) += u->dt/2.0 *gx;
        u->iy(2, i, jR, k) += u->dt/2.0 *gy; u->iy(2, i, jL, k) += u->dt/2.0 *gy;
        u->iz(3, i, j, kR) += u->dt/2.0 *gz; u->iz(3, i, j, kL) += u->dt/2.0 *gz;
    }
    else {
        long double dpdx_p = (Phii(u, i+1, j, k)-Phii(u, i, j, k))/u->dx;   long double dpdx_m = (Phii(u, i, j, k)-Phii(u, i-1, j, k))/u->dx;
        long double dpdy_p = (Phii(u, i, j+1, k)-Phii(u, i, j, k))/u->dy;   long double dpdy_m = (Phii(u, i, j, k)-Phii(u, i, j-1, k))/u->dy;
        long double dpdz_p = (Phii(u, i, j, k+1)-Phii(u, i, j, k))/u->dz;   long double dpdz_m = (Phii(u, i, j, k)-Phii(u, i, j, k-1))/u->dz;

        u->ix(1, iR, j, k) -= u->dt/2.0 *dpdx_m; u->ix(1, iL, j, k) -= u->dt/2.0 *dpdx_p;
        u->iy(2, i, jR, k) -= u->dt/2.0 *dpdy_m; u->iy(2, i, jL, k) -= u->dt/2.0 *dpdy_p;
        u->iz(3, i, j, kR) -= u->dt/2.0 *dpdz_m; u->iz(3, i, j, kL) -= u->dt/2.0 *dpdz_p;
    }


    /* ----- 11. Set Bn at the n interface from the values in the uC arrays ----- */
    u->ix(iBx, iR, j, k) = u->uC(iBx, i, j, k); u->ix(iBx, iL, j, k) = u->uC(iBx, i+1, j, k);
    u->iy(iBy, i, jR, k) = u->uC(iBy, i, j, k); u->iy(iBy, i, jL, k) = u->uC(iBy, i, j+1, k);
    u->iz(iBz, i, j, kR) = u->uC(iBz, i, j, k); u->iz(iBz, i, j, kL) = u->uC(iBz, i, j, k+1);


    /* ----- 12. Transform the interface values to conserved variables ----- */
    if (!toConserved(u, 0, iR, j, k, gamma)) std::cout << "ERROR:::plm.h::PLM:: Transformation of intx (R) to conserved variables failed !" << std::endl;
    if (!toConserved(u, 0, iL, j, k, gamma)) std::cout << "ERROR:::plm.h::PLM:: Transformation of intx (R) to conserved variables failed !" << std::endl;
    if (!toConserved(u, 1, i, jR, k, gamma)) std::cout << "ERROR:::plm.h::PLM:: Transformation of intx (R) to conserved variables failed !" << std::endl;
    if (!toConserved(u, 1, i, jL, k, gamma)) std::cout << "ERROR:::plm.h::PLM:: Transformation of intx (R) to conserved variables failed !" << std::endl;
    if (!toConserved(u, 2, i, j, kR, gamma)) std::cout << "ERROR:::plm.h::PLM:: Transformation of intx (R) to conserved variables failed !" << std::endl;
    if (!toConserved(u, 2, i, j, kL, gamma)) std::cout << "ERROR:::plm.h::PLM:: Transformation of intx (R) to conserved variables failed !" << std::endl;
    
    return true;
}

#else
bool PLM(Arrays *u, int i, int j, int k) {

}

#endif // MHD