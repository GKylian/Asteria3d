#pragma once
#include "Arrays.h"
#include "MHD.h"




bool HDPLM(Arrays *us, int xi, int yi, Arrays *int_x, Arrays *int_y, long double dt, long double gx, long double gy, long double gz, long double rs) {
    long double dx = us->dx; long double dy = us->dy;

    //0.1) Transform the variables at (i-1, i, i+1) and (j-1, (j,) j+1)
    long double wi[8] = { 0 }; for (int i = 0; i < 5; i++) wi[i] = us->at(i, xi, yi);
    if (!toPrimitive(wi, wi)) cout << "First conversion to primitive failed." << endl;
    long double wj[8] = { 0 }; copy(wi, wi+8, wj); swap(wj[1], wj[2]);

    long double wip[8] = { 0 }; for (int i = 0; i < 5; i++) wip[i] = us->at(i, xi+1, yi);
    if (!toPrimitive(wip, wip)) cout << "First conversion to primitive failed." << endl;
    long double wim[8] = { 0 }; for (int i = 0; i < 5; i++) wim[i] = us->at(i, xi-1, yi);
    if (!toPrimitive(wim, wim)) cout << "First conversion to primitive failed." << endl;



    long double wjp[8] = { 0 }; for (int i = 0; i < 5; i++) wjp[i] = us->at(i, xi, yi+1);
    if (!toPrimitive(wjp, wjp)) cout << "First conversion to primitive failed." << endl; swap(wjp[1], wjp[2]);
    long double wjm[8] = { 0 }; for (int i = 0; i < 5; i++) wjm[i] = us->at(i, xi, yi-1);
    if (!toPrimitive(wjm, wjm)) cout << "First conversion to primitive failed." << endl; swap(wjm[1], wjm[2]);


    //0.2) Get eigensystems
    long double eigen_x[5] = { 0 }; long double Lx[5][5] = { 0 }; long double Rx[5][5] = { 0 };
    getEigenHD(wi, eigen_x, Lx, Rx);
    long double eigen_y[5] = { 0 }; long double Ly[5][5] = { 0 }; long double Ry[5][5] = { 0 };
    getEigenHD(wj, eigen_y, Ly, Ry);

    long double eigM_x = (long double)*max_element(eigen_x, eigen_x+5); long double eig0_x = (long double)*min_element(eigen_x, eigen_x+5);
    long double eigM_y = (long double)*max_element(eigen_y, eigen_y+5); long double eig0_y = (long double)*min_element(eigen_y, eigen_y+5);



    //1) Compute L/R, centered and van Leer differences.
    long double dwL_x[5] = { 0 }; long double dwR_x[5] = { 0 }; long double dwC_x[5] = { 0 }; long double dwG_x[5] = { 0 };
    long double dwL_y[5] = { 0 }; long double dwR_y[5] = { 0 }; long double dwC_y[5] = { 0 }; long double dwG_y[5] = { 0 };
    //x: (rho, vx, vy, vz, P);    y: (rho, vy, vx, vz, P)

    for (int i = 0; i < 5; i++)
    {

        dwL_x[i] = wi[i]-wim[i];   dwR_x[i] = wip[i]-wi[i];   dwC_x[i] = (wip[i]-wim[i])/2.0; if (dwL_x[i]*dwR_x[i] > 0) dwG_x[i] = 2.0*dwL_x[i]*dwR_x[i]/(dwL_x[i]+dwR_x[i]);
        dwL_y[i] = wj[i]-wjm[i];   dwR_y[i] = wjp[i]-wj[i];   dwC_y[i] = (wjp[i]-wjm[i])/2.0; if (dwL_y[i]*dwR_y[i] > 0) dwG_y[i] = 2.0*dwL_y[i]*dwR_y[i]/(dwL_y[i]+dwR_y[i]);
    }




    //2) Project the differences along characteristics
    long double daL_x[5] = { 0 }; long double daR_x[5] = { 0 }; long double daC_x[5] = { 0 }; long double daG_x[5] = { 0 };
    long double daL_y[5] = { 0 }; long double daR_y[5] = { 0 }; long double daC_y[5] = { 0 }; long double daG_y[5] = { 0 };
    for (int n = 0; n < 5; n++)
        for (int m = 0; m < 5; m++)
        {
            daL_x[n] += Lx[n][m]*dwL_x[m]; daR_x[n] += Lx[n][m]*dwR_x[m]; daC_x[n] += Lx[n][m]*dwC_x[m]; daG_x[n] += Lx[n][m]*dwG_x[m];
            daL_y[n] += Ly[n][m]*dwL_y[m]; daR_y[n] += Ly[n][m]*dwR_y[m]; daC_y[n] += Ly[n][m]*dwC_y[m]; daG_y[n] += Ly[n][m]*dwG_y[m];
        }




    //3) Apply monotonicity constraint
    long double da_x[5] = { 0 }; long double da_y[5] = { 0 };
    for (int n = 0; n < 5; n++)
    {
        /*long double lim1_x = fminl(fabsl(daL_x[n]), fabsl(daR_x[n])); long double lim2_x = fminl(fabsl(daC_x[n]), fabsl(daG_x[n]));
        long double lim1_y = fminl(fabsl(daL_y[n]), fabsl(daR_y[n])); long double lim2_y = fminl(fabsl(daC_y[n]), fabsl(daG_y[n]));

        da_x[n] = Sgn(daC_x[n]) * fminl(2.0*lim1_x, lim2_x);
        da_y[n] = Sgn(daC_y[n]) * fminl(2.0*lim1_y, lim2_y);*/

        da_x[n] = sgn(daC_x[n])*fminl(fminl(2.0*fabsl(daL_x[n]), 2.0*fabsl(daR_x[n])), fabsl(daC_x[n]));
        da_y[n] = sgn(daC_y[n])*fminl(fminl(2.0*fabsl(daL_y[n]), 2.0*fabsl(daR_y[n])), fabsl(daC_y[n]));
        //cout << dwR_x[n] << ", " << dwR_y[n] << endl;
    }


    //FIXED: --- Selon la notation dans le papier, on devrait avoir Rx[m][n] et Ry[m][n]
    //FIXED: --- [n][m] correspont à L.w puis R.a, et R.(L.w) = w, ce qui semble correcte
    //FIXED: --- [m][n] correspont à L.w puis a.R, et (L.w).R != w, ce qui semble plus étrange
    //FIXED: --- Mais d'un autre côté faire la même opération pour L et R alors que l'un a les vecteurs propres en ligne et l'autre en colonnes semble bizarre.
    //FIXED: --- --> Après des tests, mn donne des mauvais résultats.

    //4) Project back to primitive variables
    long double dw_x[5] = { 0 }; long double dw_y[5] = { 0 };
    for (int n = 0; n < 5; n++)
        for (int m = 0; m < 5; m++)
        {
            dw_x[n] += da_x[m]*Rx[n][m];
            dw_y[n] += da_y[m]*Ry[n][m];
        }
    //cout << dw_x[4] << ", " << dw_y[4] << endl;
    //cout << "x: " << dw_x[0] << ", " << dw_x[1] << ", " << dw_x[2] << ", " << dw_x[3] << ", " << dw_x[4] << endl;
    //cout << "y: " << dw_y[0] << ", " << dw_y[1] << ", " << dw_y[2] << ", " << dw_y[3] << ", " << dw_y[4] << endl;

    //5) Compute R(i-1/2) and L(i+1/2) values
    long double wL_x[5] = { 0 }; long double wR_x[5] = { 0 };
    long double wL_y[5] = { 0 }; long double wR_y[5] = { 0 };
    for (int i = 0; i < 5; i++)
    {
        //FIXED: --- The first version (in two steps, with range checking) seems to give (slightly) better results


        wL_x[i] = wi[i] + 0.5*dw_x[i];
        wR_x[i] = wi[i] - 0.5*dw_x[i];

        wL_y[i] = wj[i] + 0.5*dw_y[i];
        wR_y[i] = wj[i] - 0.5*dw_y[i];

        wL_x[i] = fmaxl(fminl(wi[i], wip[i]), wL_x[i]); wL_x[i] = fminl(fmaxl(wi[i], wip[i]), wL_x[i]);
        wR_x[i] = fmaxl(fminl(wi[i], wim[i]), wR_x[i]); wR_x[i] = fminl(fmaxl(wi[i], wim[i]), wR_x[i]);

        wL_y[i] = fmaxl(fminl(wj[i], wjp[i]), wL_y[i]); wL_y[i] = fminl(fmaxl(wj[i], wjp[i]), wL_y[i]);
        wR_y[i] = fmaxl(fminl(wj[i], wjm[i]), wR_y[i]); wR_y[i] = fminl(fmaxl(wj[i], wjm[i]), wR_y[i]);

        wL_x[i] += -dt/(2.0*dx)*fmaxl(eigM_x, 0.0)*dw_x[i];
        wR_x[i] += -dt/(2.0*dx)*fminl(eig0_x, 0.0)*dw_x[i];

        wL_y[i] += -dt/(2.0*dy)*fmaxl(eigM_y, 0.0)*dw_y[i];
        wR_y[i] += -dt/(2.0*dy)*fminl(eig0_y, 0.0)*dw_y[i];


        //FIXED: --- in wR_x: +dt or -dt ? It says -dt in --> +dt is correct

        /*wL_x[i] = wi[i] + (0.5 - dt/(2.0*dx)*fmaxl(eigM_x, 0.0))*dw_x[i];
        wR_x[i] = wi[i] - (0.5 + dt/(2.0*dx)*fminl(eig0_x, 0.0))*dw_x[i];

        wL_y[i] = wi[i] + (0.5 - dt/(2.0*dy)*fmaxl(eigM_y, 0.0))*dw_y[i];
        wR_y[i] = wi[i] - (0.5 + dt/(2.0*dy)*fminl(eig0_y, 0.0))*dw_y[i];*/
    }


    //Subtract that don't reach the interface and/or (for HLL) that move away from the interface.
    long double qa = 0.0; long double qx;
    for (int n = 0; n < 5; n++)
    {
        if (eigen_x[n] >= 0.0) {
            //Substract amount of each wave that doesn't reach the interface
            qa = 0.0; qx = 0.5*dt/dx*(eigM_x-eigen_x[n]);
            for (int m = 0; m < 5; m++)
                qa += Lx[n][m]*qx*dw_x[m];
            for (int m = 0; m < 5; m++)
                wL_x[m] += qa*Rx[m][n];

            //For HLL, substract those moving AWAY from the interface
            qa = 0.0; qx = 0.5*dt/dx*(eig0_x-eigen_x[n]);
            for (int m = 0; m < 5; m++)
                qa += Lx[n][m]*qx*dw_x[m];
            for (int m = 0; m < 5; m++)
                wR_x[m] += qa*Rx[m][n];
        }

        if (eigen_x[n] <= 0.0) {
            //Substract amount of each wave that doesn't reach the interface
            qa = 0.0; qx = 0.5*dt/dx*(eig0_x-eigen_x[n]);
            for (int m = 0; m < 5; m++)
                qa += Lx[n][m]*qx*dw_x[m];
            for (int m = 0; m < 5; m++)
                wR_x[m] += qa*Rx[m][n];

            //For HLL, substract those moving AWAY from the interface
            qa = 0.0; qx = 0.5*dt/dx*(eigM_x-eigen_x[n]);
            for (int m = 0; m < 5; m++)
                qa += Lx[n][m]*qx*dw_x[m];
            for (int m = 0; m < 5; m++)
                wL_x[m] += qa*Rx[m][n];
        }
    }

    for (int n = 0; n < 5; n++)
    {
        if (eigen_y[n] >= 0.0) {
            //Substract amount of each wave that doesn't reach the interface
            qa = 0.0; qx = 0.5*dt/dy*(eigM_y-eigen_y[n]);
            for (int m = 0; m < 5; m++)
                qa += Ly[n][m]*qx*dw_y[m];
            for (int m = 0; m < 5; m++)
                wL_y[m] += qa*Ry[m][n];

            //For HLL, substract those moving AWAY from the interface
            qa = 0.0; qx = 0.5*dt/dy*(eig0_y-eigen_y[n]);
            for (int m = 0; m < 5; m++)
                qa += Ly[n][m]*qx*dw_y[m];
            for (int m = 0; m < 5; m++)
                wR_y[m] += qa*Ry[m][n];
        }

        if (eigen_y[n] <= 0.0) {
            //Substract amount of each wave that doesn't reach the interface
            qa = 0.0; qx = 0.5*dt/dy*(eig0_y-eigen_y[n]);
            for (int m = 0; m < 5; m++)
                qa += Ly[n][m]*qx*dw_y[m];
            for (int m = 0; m < 5; m++)
                wR_y[m] += qa*Ry[m][n];

            //For HLL, substract those moving AWAY from the interface
            qa = 0.0; qx = 0.5*dt/dy*(eigM_y-eigen_y[n]);
            for (int m = 0; m < 5; m++)
                qa += Ly[n][m]*qx*dw_y[m];
            for (int m = 0; m < 5; m++)
                wL_y[m] += qa*Ry[m][n];
        }
    }


    // Add gravity-related source terms for M and E.
    long double x = us->x0 + (xi-4.0)*dx; long double y = us->y0 + (yi-4.0)*dy;

    if (rs >= 0.0) {
        long double xm = x-dx/2.0; long double xp = x+dx/2.0;   long double ym = y-dy/2.0; long double yp = y+dy/2.0;

        //S = 0, -rh*dphidx, -rh*dphidy, 0, -rho*(vx*dphidx+vy*dphidy)
        long double dpdx_m = Phix(xm, y, rs); long double dpdx_p = Phix(xp, y, rs);
        long double dpdy_m = Phiy(x, ym, rs); long double dpdy_p = Phiy(x, yp, rs);

        /*wL_x[1] -= dt/2.0*dpdx_p; wR_x[1] -= dt/2.0*dpdx_m;
        wL_y[1] -= dt/2.0*dpdy_p; wR_y[1] -= dt/2.0*dpdy_m;*/

    }
    else {
        wL_x[1] += dt/2.0*gx; wR_x[1] += dt/2.0*gx;
        wL_y[1] += dt/2.0*gy; wR_y[1] += dt/2.0*gy;
    }

    if (wL_x[0] <= 1e-6) wL_x[0] = 1e-6; if (wR_x[0] <= 1e-6) wR_x[0] = 1e-6; if (wL_y[0] <= 1e-6) wL_y[0] = 1e-6; if (wR_y[0] <= 1e-6) wR_y[0] = 1e-6;
    if (wL_x[4] <= 1e-6) wL_x[4] = 1e-6; if (wR_x[4] <= 1e-6) wR_x[4] = 1e-6; if (wL_y[4] <= 1e-6) wL_y[4] = 1e-6; if (wR_y[4] <= 1e-6) wR_y[4] = 1e-6;

    long double uxL[8] = { wL_x[0], wL_x[1], wL_x[2], wL_x[3], wL_x[4], 0, 0, 0 }; long double uxR[8] = { wR_x[0], wR_x[1], wR_x[2], wR_x[3], wR_x[4], 0, 0, 0 };
    long double uyL[8] = { wL_y[0], wL_y[2], wL_y[1], wL_y[3], wL_y[4], 0, 0, 0 }; long double uyR[8] = { wR_y[0], wR_y[2], wR_y[1], wR_y[3], wR_y[4], 0, 0, 0 };


    if (wL_x[4] <= 0.0 || wR_x[4] <= 0.0 || wL_y[4] <= 0.0 || wR_y[4] <= 0.0) {
        cout << "RANGE-ERROR:: P <= 0 at the end of the reconstruction: " << endl;
        cout << wL_x[4] << ", " << wR_x[4] << ";  " << wL_y[4] << ", " << wR_y[4] << endl;
        return false;
    }

    toConserved(uxL, uxL); toConserved(uxR, uxR);
    toConserved(uyL, uyL); toConserved(uyR, uyR);

    for (int i = 0; i < 8; i++)
    {
        int_x->at(i, xi*2+1, yi) = uxL[i]; //left of i+1/2
        int_x->at(i, xi*2, yi) = uxR[i]; //right of i-1/2
        int_y->at(i, xi, yi*2+1) = uyL[i]; //left of j+1/2
        int_y->at(i, xi, yi*2) = uyR[i]; //right of j-1/2

        if (isnan(uxL[i]) || isnan(uxR[i]) ||isnan(uyL[i]) ||isnan(uyR[i])) {
            cout << "PPM-Char:: at least one of the final values is NaN" << endl; return false;
        }
    }

    return true;
}



//TODO:-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//TODO:-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//TODO:-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//TODO:-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//TODO:-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------









bool PLM(Arrays *us, int xi, int yi, Arrays *int_x, Arrays *int_y, long double dt, long double gx, long double gy, long double gz, bool CharTrac, long double rs) {
    long double dx = us->dx; long double dy = us->dy;



    //0.1) Transform the variables at (i-1, i, i+1) and (j-1, (j,) j+1) from u = (rho, Mx, My, Mz, E, Bx, By, Bz) to w = (rho, vn, v1, v2, P, B1, B2, Bn)
    //For -x the order is (rho, vx, vy, vz, P, By, Bz, Bx)
    //For -y the order is (rho, vy, vx, vz, P, Bx, Bz, By)

    long double wi[8] = { 0 }; for (int i = 0; i < 8; i++) wi[i] = us->at(i, xi, yi);
    if (!toPrimitive(wi, wi)) cout << "First conversion to primitive failed." << endl; swap(wi[6], wi[7]);
    wi[7] = (us->at(5, xi+1, yi) + us->at(5, xi, yi))/2.0; wi[5] = (us->at(6, xi, yi+1) + us->at(6, xi, yi))/2.0; //Bx at 7 and By at 5

    long double wj[8] = { 0 }; copy(wi, wi+8, wj); swap(wj[1], wj[2]); swap(wj[5], wj[7]);

    long double wip[8] = { 0 }; for (int i = 0; i < 8; i++) wip[i] = us->at(i, xi+1, yi);
    if (!toPrimitive(wip, wip)) cout << "First conversion to primitive failed." << endl; swap(wip[6], wip[7]);
    wip[7] = (us->at(5, xi+2, yi) + us->at(5, xi+1, yi))/2.0; wip[5] = (us->at(6, xi+1, yi+1) + us->at(6, xi+1, yi))/2.0; //Bx at 7 and By at 5

    long double wim[8] = { 0 }; for (int i = 0; i < 8; i++) wim[i] = us->at(i, xi-1, yi);
    if (!toPrimitive(wim, wim)) cout << "First conversion to primitive failed." << endl; swap(wim[6], wim[7]);
    wim[7] = (us->at(5, xi, yi) + us->at(5, xi-1, yi))/2.0; wim[5] = (us->at(6, xi-1, yi+1) + us->at(6, xi-1, yi))/2.0; //Bx at 7 and By at 5


    long double wjp[8] = { 0 }; for (int i = 0; i < 8; i++) wjp[i] = us->at(i, xi, yi+1);
    if (!toPrimitive(wjp, wjp)) cout << "First conversion to primitive failed." << endl; swap(wjp[1], wjp[2]); swap(wjp[6], wjp[7]);
    wjp[5] = (us->at(5, xi+1, yi+1)+us->at(5, xi, yi+1))/2.0; wjp[7] = (us->at(6, xi, yi+2) + us->at(6, xi, yi+1))/2.0; //Bx at 5 and By at 7

    long double wjm[8] = { 0 }; for (int i = 0; i < 8; i++) wjm[i] = us->at(i, xi, yi-1);
    if (!toPrimitive(wjm, wjm)) cout << "First conversion to primitive failed." << endl; swap(wjm[1], wjm[2]); swap(wjm[6], wjm[7]);
    wjm[5] = (us->at(5, xi+1, yi-1)+us->at(5, xi, yi-1))/2.0; wjm[7] = (us->at(6, xi, yi) + us->at(6, xi, yi-1))/2.0; //Bx at 5 and By at 7




    //0.2) Compute the eigensystem
    long double eigen_x[7] = { 0 }; long double Lx[7][7] = { 0 }; long double Rx[7][7] = { 0 };
    getEigen(wi, eigen_x, Lx, Rx);
    long double eigen_y[7] = { 0 }; long double Ly[7][7] = { 0 }; long double Ry[7][7] = { 0 };
    getEigen(wj, eigen_y, Ly, Ry);

    long double eigM_x = (long double)*max_element(eigen_x, eigen_x+7); long double eig0_x = (long double)*min_element(eigen_x, eigen_x+7);
    long double eigM_y = (long double)*max_element(eigen_y, eigen_y+7); long double eig0_y = (long double)*min_element(eigen_y, eigen_y+7);





    //1) Compute left, center and right differences in primitive variables.
    long double dwL_x[7] = { 0 }; long double dwR_x[7] = { 0 }; long double dwC_x[7] = { 0 };
    long double dwL_y[7] = { 0 }; long double dwR_y[7] = { 0 }; long double dwC_y[7] = { 0 };
    //x: (rho, vx, vy, vz, P, By, Bz, Bx);    y: (rho, vy, vx, vz, P, Bx, Bz, By)

    for (int i = 0; i < 7; i++)
    {
        dwL_x[i] = wi[i]-wim[i];   dwR_x[i] = wip[i]-wi[i];   dwC_x[i] = (wip[i]-wim[i])/2.0;
        dwL_y[i] = wj[i]-wjm[i];   dwR_y[i] = wjp[i]-wj[i];   dwC_y[i] = (wjp[i]-wjm[i])/2.0;
    }





    //2) Project the differences along characteristics
    long double daL_x[7] = { 0 }; long double daR_x[7] = { 0 }; long double daC_x[7] = { 0 }; long double daG_x[7] = { 0 };
    long double daL_y[7] = { 0 }; long double daR_y[7] = { 0 }; long double daC_y[7] = { 0 }; long double daG_y[7] = { 0 };
    for (int n = 0; n < 7; n++)
    for (int m = 0; m < 7; m++)
    {
        daL_x[n] += Lx[n][m]*dwL_x[m]; daR_x[n] += Lx[n][m]*dwR_x[m]; daC_x[n] += Lx[n][m]*dwC_x[m];
        daL_y[n] += Ly[n][m]*dwL_y[m]; daR_y[n] += Ly[n][m]*dwR_y[m]; daC_y[n] += Ly[n][m]*dwC_y[m];
    }


    //3) Apply monotonicity constraint
    long double da_x[7] = { 0 }; long double da_y[7] = { 0 };
    for (int n = 0; n < 7; n++)
    {

        da_x[n] = sgn(daC_x[n])*fminl(fminl(2.0*fabsl(daL_x[n]), 2.0*fabsl(daR_x[n])), fabsl(daC_x[n]));
        da_y[n] = sgn(daC_y[n])*fminl(fminl(2.0*fabsl(daL_y[n]), 2.0*fabsl(daR_y[n])), fabsl(daC_y[n]));
        if (null(da_x[n])) da_x[n] = 0.0; if (null(da_y[n])) da_y[n] = 0.0;
    }


    //4) Project back to primitive variables
    long double dw_x[7] = { 0 }; long double dw_y[7] = { 0 };
    for (int n = 0; n < 7; n++)
    for (int m = 0; m < 7; m++)
    {
        dw_x[n] += da_x[m]*Rx[n][m];
        dw_y[n] += da_y[m]*Ry[n][m];
        if (null(dw_x[n])) dw_x[n] = 0.0; if (null(dw_y[n])) dw_y[n] = 0.0;
    }


    //5.1) Compute R(i-1/2) and L(i+1/2) values
    long double wL_x[7] = { 0 }; long double wR_x[7] = { 0 };
    long double wL_y[7] = { 0 }; long double wR_y[7] = { 0 };
    for (int i = 0; i < 7; i++)
    {
        wL_x[i] = wi[i] + 0.5*dw_x[i];
        wR_x[i] = wi[i] - 0.5*dw_x[i];

        wL_y[i] = wj[i] + 0.5*dw_y[i];
        wR_y[i] = wj[i] - 0.5*dw_y[i];

        wL_x[i] = fmaxl(fminl(wi[i], wip[i]), wL_x[i]); wL_x[i] = fminl(fmaxl(wi[i], wip[i]), wL_x[i]);

        wR_x[i] = fmaxl(fminl(wi[i], wim[i]), wR_x[i]); wR_x[i] = fminl(fmaxl(wi[i], wim[i]), wR_x[i]);

        wL_y[i] = fmaxl(fminl(wj[i], wjp[i]), wL_y[i]); wL_y[i] = fminl(fmaxl(wj[i], wjp[i]), wL_y[i]);
        wR_y[i] = fmaxl(fminl(wj[i], wjm[i]), wR_y[i]); wR_y[i] = fminl(fmaxl(wj[i], wjm[i]), wR_y[i]);


        wL_x[i] += -dt/(2.0*dx)*fmaxl(eigM_x, 0.0)*dw_x[i];
        wR_x[i] += -dt/(2.0*dx)*fminl(eig0_x, 0.0)*dw_x[i];

        wL_y[i] += -dt/(2.0*dy)*fmaxl(eigM_y, 0.0)*dw_y[i];
        wR_y[i] += -dt/(2.0*dy)*fminl(eig0_y, 0.0)*dw_y[i];

        

        //if (null(wL_x[i])) wL_x[i] = 0.0;   if (null(wR_x[i])) wR_x[i] = 0.0;
        //if (null(wL_y[i])) wL_y[i] = 0.0;   if (null(wR_y[i])) wR_y[i] = 0.0;

    }

    //if (wL_x[6] != 0.0)
    //    cout << wi[6] << "; " << wip[6] << endl;





    if(CharTrac){
    //Subtract that don't reach the interface and/or (for HLL) that move away from the interface.
    long double qa = 0.0; long double qx;
    for (int n = 0; n < 7; n++)
    {
        if (eigen_x[n] >= 0.0) {
            //Substract amount of each wave that doesn't reach the interface
            qa = 0.0; qx = 0.5*dt/dx*(eigM_x-eigen_x[n]);
            for (int m = 0; m < 7; m++)
                qa += Lx[n][m]*qx*dw_x[m];
            if (null(qa)) qa = 0.0;
            for (int m = 0; m < 7; m++)
                wL_x[m] += qa*Rx[m][n];
            

            //For HLL, substract those moving AWAY from the interface
            qa = 0.0; qx = 0.5*dt/dx*(eig0_x-eigen_x[n]);
            for (int m = 0; m < 7; m++)
                qa += Lx[n][m]*qx*dw_x[m];
            if (null(qa)) qa = 0.0;
            for (int m = 0; m < 7; m++)
                wR_x[m] += qa*Rx[m][n];
        }

        if (eigen_x[n] <= 0.0) {
            //Substract amount of each wave that doesn't reach the interface
            qa = 0.0; qx = 0.5*dt/dx*(eig0_x-eigen_x[n]);
            for (int m = 0; m < 7; m++)
                qa += Lx[n][m]*qx*dw_x[m];
            if (null(qa)) qa = 0.0;
            for (int m = 0; m < 7; m++)
                wR_x[m] += qa*Rx[m][n];

            //For HLL, substract those moving AWAY from the interface
            qa = 0.0; qx = 0.5*dt/dx*(eigM_x-eigen_x[n]);
            for (int m = 0; m < 7; m++)
                qa += Lx[n][m]*qx*dw_x[m];
            if (null(qa)) qa = 0.0;
            for (int m = 0; m < 7; m++)
                wL_x[m] += qa*Rx[m][n];
        }
    }

    for (int n = 0; n < 7; n++)
    {
        if (eigen_y[n] >= 0.0) {
            //Substract amount of each wave that doesn't reach the interface
            qa = 0.0; qx = 0.5*dt/dy*(eigM_y-eigen_y[n]);
            for (int m = 0; m < 7; m++)
                qa += Ly[n][m]*qx*dw_y[m];
            if (null(qa)) qa = 0.0;
            for (int m = 0; m < 7; m++)
                wL_y[m] += qa*Ry[m][n];

            //For HLL, substract those moving AWAY from the interface
            qa = 0.0; qx = 0.5*dt/dy*(eig0_y-eigen_y[n]);
            for (int m = 0; m < 7; m++)
                qa += Ly[n][m]*qx*dw_y[m];
            if (null(qa)) qa = 0.0;
            for (int m = 0; m < 7; m++)
                wR_y[m] += qa*Ry[m][n];
        }

        if (eigen_y[n] <= 0.0) {
            //Substract amount of each wave that doesn't reach the interface
            qa = 0.0; qx = 0.5*dt/dy*(eig0_y-eigen_y[n]);
            for (int m = 0; m < 7; m++)
                qa += Ly[n][m]*qx*dw_y[m];
            if (null(qa)) qa = 0.0;
            for (int m = 0; m < 7; m++)
                wR_y[m] += qa*Ry[m][n];

            //For HLL, substract those moving AWAY from the interface
            qa = 0.0; qx = 0.5*dt/dy*(eigM_y-eigen_y[n]);
            for (int m = 0; m < 7; m++)
                qa += Ly[n][m]*qx*dw_y[m];
            if (null(qa)) qa = 0.0;
            for (int m = 0; m < 7; m++)
                wL_y[m] += qa*Ry[m][n];
        }

        if (null(wL_x[n])) wL_x[n] = 0.0;   if (null(wR_x[n])) wR_x[n] = 0.0;
        if (null(wL_y[n])) wL_y[n] = 0.0;   if (null(wR_y[n])) wR_y[n] = 0.0;
    }
    } // if(CharTrac)





    // Adds dBy_L/R(i+-1/2) and dBx_L/R(j+-1/2)
    long double dBy = dt/2.0 * wi[2]*(us->at(5, xi+1, yi)-us->at(5, xi, yi))/dx;   long double dBx = dt/2.0 * wi[1]*(us->at(6, xi, yi+1)-us->at(6, xi, yi))/dy;

    wL_x[5] += dBy; wR_x[5] += dBy;
    wL_y[5] += dBx; wR_y[5] += dBx;



    // Add gravity-related source terms for M and E.
    long double x = us->x0 + (xi-4.0)*dx; long double y = us->y0 + (yi-4.0)*dy;

    if (rs >= 0.0) {
        long double xm = x-dx; long double xp = x+dx;   long double ym = y-dy; long double yp = y+dy;

        //S = 0, -rh*dphidx, -rh*dphidy, 0, -rho*(vx*dphidx+vy*dphidy)
        long double dpdx_m = Phix(xm, y, rs); long double dpdx_p = Phix(xp, y, rs);
        long double dpdy_m = Phiy(x, ym, rs); long double dpdy_p = Phiy(x, yp, rs);

        wL_x[1] -= dt/2.0*dpdx_p; wR_x[1] -= dt/2.0*dpdx_m;
        wL_y[1] -= dt/2.0*dpdy_p; wR_y[1] -= dt/2.0*dpdy_m;

    }
    else {
        wL_x[1] += dt/2.0*gx; wR_x[1] += dt/2.0*gx;
        wL_y[1] += dt/2.0*gy; wR_y[1] += dt/2.0*gy;
    }

    





    // Get Bx for R(i-1/2) and L(i+1/2), and By for R(j-1/2) and L(j+1/2)
    long double Bx_R = us->at(5, xi, yi); long double Bx_L = us->at(5, xi+1, yi);
    long double By_R = us->at(6, xi, yi); long double By_L = us->at(6, xi, yi+1);

    //x: (rho, vx, vy, vz, P, By, Bz, Bx);    y: (rho, vy, vx, vz, P, Bx, Bz, By)
    long double uxL[8] = { wL_x[0], wL_x[1], wL_x[2], wL_x[3], wL_x[4], Bx_L, wL_x[5], wL_x[6] }; long double uxR[8] = { wR_x[0], wR_x[1], wR_x[2], wR_x[3], wR_x[4], Bx_R, wR_x[5], wR_x[6] };
    long double uyL[8] = { wL_y[0], wL_y[2], wL_y[1], wL_y[3], wL_y[4], wL_y[5], By_L, wL_y[6] }; long double uyR[8] = { wR_y[0], wR_y[2], wR_y[1], wR_y[3], wR_y[4], wR_y[5], By_R, wR_y[6] };


    if (wL_x[4] <= 0.0 || wR_x[4] <= 0.0 || wL_y[4] <= 0.0 || wR_y[4] <= 0.0) {
        cout << "RANGE-ERROR:: P <= 0 at the end of the reconstruction: " << endl;
        cout << wL_x[4] << ", " << wR_x[4] << ";  " << wL_y[4] << ", " << wR_y[4] << endl;
        return false;
    }

    toConserved(uxL, uxL); toConserved(uxR, uxR);
    toConserved(uyL, uyL); toConserved(uyR, uyR);

    for (int i = 0; i < 8; i++)
    {
        int_x->at(i, xi*2+1, yi) = uxL[i]; //left of i+1/2
        int_x->at(i, xi*2, yi) = uxR[i]; //right of i-1/2
        int_y->at(i, xi, yi*2+1) = uyL[i]; //left of j+1/2
        int_y->at(i, xi, yi*2) = uyR[i]; //right of j-1/2

        if (isnan(uxL[i]) || isnan(uxR[i]) ||isnan(uyL[i]) ||isnan(uyR[i])) {
            cout << "PPM-Char:: at least one of the final values is NaN" << endl; return false;
        }
    }

    //if (fabsl(uxL[6]) > 1.0 || fabsl(uxR[6]) > 1.0) {
    //    cout << "-1.0 <= By_L/R <= 1.0 is now false in reconstruction." << endl;
    //}

    return true;
}



//TODO:-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//TODO:-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//TODO:-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



bool PLM_Simple(Arrays *us, int xi, int yi, Arrays *int_x, Arrays *int_y, long double dt, long double gx, long double gy, long double gz, bool CharTrac) {
    long double dx = us->dx; long double dy = us->dy;



    //0.1) Transform the variables at (i-1, i, i+1) and (j-1, (j,) j+1) from u = (rho, Mx, My, Mz, E, Bx, By, Bz) to w = (rho, vn, v1, v2, P, B1, B2, Bn)
    //For -x the order is (rho, vx, vy, vz, P, By, Bz, Bx)
    //For -y the order is (rho, vy, vx, vz, P, Bx, Bz, By)

    long double wi[8] = { 0 }; for (int i = 0; i < 8; i++) wi[i] = us->at(i, xi, yi);
    if (!toPrimitive(wi, wi)) cout << "First conversion to primitive failed." << endl; swap(wi[6], wi[7]);
    wi[7] = (us->at(5, xi+1, yi) + us->at(5, xi, yi))/2.0; wi[5] = (us->at(6, xi, yi+1) + us->at(6, xi, yi))/2.0; //Bx at 7 and By at 5

    long double wj[8] = { 0 }; copy(wi, wi+8, wj); swap(wj[1], wj[2]); swap(wj[5], wj[7]);

    long double wip[8] = { 0 }; for (int i = 0; i < 8; i++) wip[i] = us->at(i, xi+1, yi);
    if (!toPrimitive(wip, wip)) cout << "First conversion to primitive failed." << endl; swap(wip[6], wip[7]);
    wip[7] = (us->at(5, xi+2, yi) + us->at(5, xi+1, yi))/2.0; wip[5] = (us->at(6, xi+1, yi+1) + us->at(6, xi+1, yi))/2.0; //Bx at 7 and By at 5

    long double wim[8] = { 0 }; for (int i = 0; i < 8; i++) wim[i] = us->at(i, xi-1, yi);
    if (!toPrimitive(wim, wim)) cout << "First conversion to primitive failed." << endl; swap(wim[6], wim[7]);
    wim[7] = (us->at(5, xi, yi) + us->at(5, xi-1, yi))/2.0; wim[5] = (us->at(6, xi-1, yi+1) + us->at(6, xi-1, yi))/2.0; //Bx at 7 and By at 5


    long double wjp[8] = { 0 }; for (int i = 0; i < 8; i++) wjp[i] = us->at(i, xi, yi+1);
    if (!toPrimitive(wjp, wjp)) cout << "First conversion to primitive failed." << endl; swap(wjp[1], wjp[2]); swap(wjp[6], wjp[7]);
    wjp[5] = (us->at(5, xi+1, yi+1)+us->at(5, xi, yi+1))/2.0; wjp[7] = (us->at(6, xi, yi+2) + us->at(6, xi, yi+1))/2.0; //Bx at 5 and By at 7

    long double wjm[8] = { 0 }; for (int i = 0; i < 8; i++) wjm[i] = us->at(i, xi, yi-1);
    if (!toPrimitive(wjm, wjm)) cout << "First conversion to primitive failed." << endl; swap(wjm[1], wjm[2]); swap(wjm[6], wjm[7]);
    wjm[5] = (us->at(5, xi+1, yi-1)+us->at(5, xi, yi-1))/2.0; wjm[7] = (us->at(6, xi, yi) + us->at(6, xi, yi-1))/2.0; //Bx at 5 and By at 7




    //0.2) Compute the eigensystem
    long double eigen_x[7] = { 0 }; long double Lx[7][7] = { 0 }; long double Rx[7][7] = { 0 };
    getEigen(wi, eigen_x, Lx, Rx);
    long double eigen_y[7] = { 0 }; long double Ly[7][7] = { 0 }; long double Ry[7][7] = { 0 };
    getEigen(wj, eigen_y, Ly, Ry);

    long double eigM_x = (long double)*max_element(eigen_x, eigen_x+7); long double eig0_x = (long double)*min_element(eigen_x, eigen_x+7);
    long double eigM_y = (long double)*max_element(eigen_y, eigen_y+7); long double eig0_y = (long double)*min_element(eigen_y, eigen_y+7);





    //1) Compute left, center and right differences in primitive variables.
    long double dwL_x[7] = { 0 }; long double dwR_x[7] = { 0 }; long double dwC_x[7] = { 0 };
    long double dwL_y[7] = { 0 }; long double dwR_y[7] = { 0 }; long double dwC_y[7] = { 0 };
    //x: (rho, vx, vy, vz, P, By, Bz, Bx);    y: (rho, vy, vx, vz, P, Bx, Bz, By)

    for (int i = 0; i < 7; i++)
    {
        dwL_x[i] = wi[i]-wim[i];   dwR_x[i] = wip[i]-wi[i];   dwC_x[i] = (wip[i]-wim[i])/2.0;
        dwL_y[i] = wj[i]-wjm[i];   dwR_y[i] = wjp[i]-wj[i];   dwC_y[i] = (wjp[i]-wjm[i])/2.0;
    }



    //3) Apply monotonicity constraint
    long double dw_x[7] = { 0 }; long double dw_y[7] = { 0 };
    for (int n = 0; n < 7; n++)
    {

        dw_x[n] = sgn(dwC_x[n])*fminl(fminl(2.0*fabsl(dwL_x[n]), 2.0*fabsl(dwR_x[n])), fabsl(dwC_x[n]));
        dw_y[n] = sgn(dwC_y[n])*fminl(fminl(2.0*fabsl(dwL_y[n]), 2.0*fabsl(dwR_y[n])), fabsl(dwC_y[n]));
        //if (null(dw_x[n])) dw_x[n] = 0.0; if (null(dw_y[n])) dw_y[n] = 0.0;
    }




    //5.1) Compute R(i-1/2) and L(i+1/2) values
    long double wL_x[7] = { 0 }; long double wR_x[7] = { 0 };
    long double wL_y[7] = { 0 }; long double wR_y[7] = { 0 };
    for (int i = 0; i < 7; i++)
    {
        wL_x[i] = wi[i] + 0.5*dw_x[i];
        wR_x[i] = wi[i] - 0.5*dw_x[i];

        wL_y[i] = wj[i] + 0.5*dw_y[i];
        wR_y[i] = wj[i] - 0.5*dw_y[i];

        wL_x[i] = fmaxl(fminl(wi[i], wip[i]), wL_x[i]); wL_x[i] = fminl(fmaxl(wi[i], wip[i]), wL_x[i]);
        wR_x[i] = fmaxl(fminl(wi[i], wim[i]), wR_x[i]); wR_x[i] = fminl(fmaxl(wi[i], wim[i]), wR_x[i]);

        wL_y[i] = fmaxl(fminl(wj[i], wjp[i]), wL_y[i]); wL_y[i] = fminl(fmaxl(wj[i], wjp[i]), wL_y[i]);
        wR_y[i] = fmaxl(fminl(wj[i], wjm[i]), wR_y[i]); wR_y[i] = fminl(fmaxl(wj[i], wjm[i]), wR_y[i]);


        wL_x[i] += -dt/(2.0*dx)*fmaxl(eigM_x, 0.0)*dw_x[i];
        wR_x[i] += -dt/(2.0*dx)*fminl(eig0_x, 0.0)*dw_x[i];

        wL_y[i] += -dt/(2.0*dy)*fmaxl(eigM_y, 0.0)*dw_y[i];
        wR_y[i] += -dt/(2.0*dy)*fminl(eig0_y, 0.0)*dw_y[i];



        //if (null(wL_x[i])) wL_x[i] = 0.0;   if (null(wR_x[i])) wR_x[i] = 0.0;
        //if (null(wL_y[i])) wL_y[i] = 0.0;   if (null(wR_y[i])) wR_y[i] = 0.0;

    }

    //if (wL_x[6] != 0.0)
    //    cout << wi[6] << "; " << wip[6] << endl;

    



    if (CharTrac) {
        //Subtract that don't reach the interface and/or (for HLL) that move away from the interface.
        long double qa = 0.0; long double qx;
        for (int n = 0; n < 7; n++)
        {
            if (eigen_x[n] >= 0.0) {
                //Substract amount of each wave that doesn't reach the interface
                qa = 0.0; qx = 0.5*dt/dx*(eigM_x-eigen_x[n]);
                for (int m = 0; m < 7; m++)
                    qa += Lx[n][m]*qx*dw_x[m];
                if (null(qa)) qa = 0.0;
                for (int m = 0; m < 7; m++)
                    wL_x[m] += qa*Rx[m][n];


                //For HLL, substract those moving AWAY from the interface
                qa = 0.0; qx = 0.5*dt/dx*(eig0_x-eigen_x[n]);
                for (int m = 0; m < 7; m++)
                    qa += Lx[n][m]*qx*dw_x[m];
                if (null(qa)) qa = 0.0;
                for (int m = 0; m < 7; m++)
                    wR_x[m] += qa*Rx[m][n];
            }

            if (eigen_x[n] <= 0.0) {
                //Substract amount of each wave that doesn't reach the interface
                qa = 0.0; qx = 0.5*dt/dx*(eig0_x-eigen_x[n]);
                for (int m = 0; m < 7; m++)
                    qa += Lx[n][m]*qx*dw_x[m];
                if (null(qa)) qa = 0.0;
                for (int m = 0; m < 7; m++)
                    wR_x[m] += qa*Rx[m][n];

                //For HLL, substract those moving AWAY from the interface
                qa = 0.0; qx = 0.5*dt/dx*(eigM_x-eigen_x[n]);
                for (int m = 0; m < 7; m++)
                    qa += Lx[n][m]*qx*dw_x[m];
                if (null(qa)) qa = 0.0;
                for (int m = 0; m < 7; m++)
                    wL_x[m] += qa*Rx[m][n];
            }
        }

        for (int n = 0; n < 7; n++)
        {
            if (eigen_y[n] >= 0.0) {
                //Substract amount of each wave that doesn't reach the interface
                qa = 0.0; qx = 0.5*dt/dy*(eigM_y-eigen_y[n]);
                for (int m = 0; m < 7; m++)
                    qa += Ly[n][m]*qx*dw_y[m];
                if (null(qa)) qa = 0.0;
                for (int m = 0; m < 7; m++)
                    wL_y[m] += qa*Ry[m][n];

                //For HLL, substract those moving AWAY from the interface
                qa = 0.0; qx = 0.5*dt/dy*(eig0_y-eigen_y[n]);
                for (int m = 0; m < 7; m++)
                    qa += Ly[n][m]*qx*dw_y[m];
                if (null(qa)) qa = 0.0;
                for (int m = 0; m < 7; m++)
                    wR_y[m] += qa*Ry[m][n];
            }

            if (eigen_y[n] <= 0.0) {
                //Substract amount of each wave that doesn't reach the interface
                qa = 0.0; qx = 0.5*dt/dy*(eig0_y-eigen_y[n]);
                for (int m = 0; m < 7; m++)
                    qa += Ly[n][m]*qx*dw_y[m];
                if (null(qa)) qa = 0.0;
                for (int m = 0; m < 7; m++)
                    wR_y[m] += qa*Ry[m][n];

                //For HLL, substract those moving AWAY from the interface
                qa = 0.0; qx = 0.5*dt/dy*(eigM_y-eigen_y[n]);
                for (int m = 0; m < 7; m++)
                    qa += Ly[n][m]*qx*dw_y[m];
                if (null(qa)) qa = 0.0;
                for (int m = 0; m < 7; m++)
                    wL_y[m] += qa*Ry[m][n];
            }

            if (null(wL_x[n])) wL_x[n] = 0.0;   if (null(wR_x[n])) wR_x[n] = 0.0;
            if (null(wL_y[n])) wL_y[n] = 0.0;   if (null(wR_y[n])) wR_y[n] = 0.0;
        }
    } // if(CharTrac)





    // Adds dBy_L/R(i+-1/2) and dBx_L/R(j+-1/2)
    long double dBy = dt/2.0 * wi[2]*(us->at(5, xi+1, yi)-us->at(5, xi, yi))/dx;   long double dBx = dt/2.0 * wi[1]*(us->at(6, xi, yi+1)-us->at(6, xi, yi))/dy;

    wL_x[5] += dBy; wR_x[5] += dBy;
    wL_y[5] += dBx; wR_y[5] += dBx;


    // Add gravity-related source terms for M and E.
    long double x = us->x0 + (xi-4.0)*dx; long double y = us->y0 + (yi-4.0)*dy;
    long double S[5] = { 0, wi[0]*gx, wi[0]*gy, wi[0]*gz, wi[0]*(wi[1]*gx+wi[2]*gy+wi[3]*gz) };

    for (int i = 0; i < 5; i++)
    {
        wL_x[i] += dt/2*S[i]; wR_x[i] += dt/2*S[i];
        wL_y[i] += dt/2*S[i]; wR_y[i] += dt/2*S[i];

    }




    //Get Bx for R(i-1/2) and L(i+1/2), and By for R(j-1/2) and L(j+1/2)
    long double Bx_R = us->at(5, xi, yi); long double Bx_L = us->at(5, xi+1, yi);
    long double By_R = us->at(6, xi, yi); long double By_L = us->at(6, xi, yi+1);

    //x: (rho, vx, vy, vz, P, By, Bz, Bx);    y: (rho, vy, vx, vz, P, Bx, Bz, By)
    long double uxL[8] = { wL_x[0], wL_x[1], wL_x[2], wL_x[3], wL_x[4], Bx_L, wL_x[5], wL_x[6] }; long double uxR[8] = { wR_x[0], wR_x[1], wR_x[2], wR_x[3], wR_x[4], Bx_R, wR_x[5], wR_x[6] };
    long double uyL[8] = { wL_y[0], wL_y[2], wL_y[1], wL_y[3], wL_y[4], wL_y[5], By_L, wL_y[6] }; long double uyR[8] = { wR_y[0], wR_y[2], wR_y[1], wR_y[3], wR_y[4], wR_y[5], By_R, wR_y[6] };



    if (wL_x[4] <= 0.0 || wR_x[4] <= 0.0 || wL_y[4] <= 0.0 || wR_y[4] <= 0.0) {
        cout << "RANGE-ERROR:: P <= 0 at the end of the reconstruction: " << endl;
        cout << wL_x[4] << ", " << wR_x[4] << ";  " << wL_y[4] << ", " << wR_y[4] << endl;
        return false;
    }

    toConserved(uxL, uxL); toConserved(uxR, uxR);
    toConserved(uyL, uyL); toConserved(uyR, uyR);

    for (int i = 0; i < 8; i++)
    {
        int_x->at(i, xi*2+1, yi) = uxL[i]; //left of i+1/2
        int_x->at(i, xi*2, yi) = uxR[i]; //right of i-1/2
        int_y->at(i, xi, yi*2+1) = uyL[i]; //left of j+1/2
        int_y->at(i, xi, yi*2) = uyR[i]; //right of j-1/2

        if (isnan(uxL[i]) || isnan(uxR[i]) ||isnan(uyL[i]) ||isnan(uyR[i])) {
            cout << "PPM-Char:: at least one of the final values is NaN" << endl; return false;
        }
    }


    return true;
}








//TODO:------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//TODO:------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//TODO:------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//TODO:------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//TODO:------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------








//Compute interfaces at (i-1/2; j) and (i; j-1/2)
bool PLMInterfaces(Arrays *us, int xi, int yi, Arrays *int_x, Arrays *int_y) {
    //rho, Mx, My, Mz, E, Bx, By, Bz
    //No need to interpolate for Bx and By
    //cout << xi << ", " << yi << endl;

    long double r_x[8] = { 0 }; long double r_y[8] = { 0 };


    for (int i = 0; i < 8; i++)
    {
        if (i == 5 || i == 6) continue;
            

        //'r' for limiter

        if (us->at(i, xi+1, yi) != us->at(i, xi, yi))
            r_x[i] = (us->at(i, xi, yi) - us->at(i, xi-1, yi))/(us->at(i, xi+1, yi) - us->at(i, xi, yi));
        if (us->at(i, xi, yi+1) != us->at(i, xi, yi))
            r_y[i] = (us->at(i, xi, yi) - us->at(i, xi, yi-1))/(us->at(i, xi, yi+1) - us->at(i, xi, yi));


        //minus,plus interfaces for x,y
        int_x->at(i, xi*2+1, yi) = us->at(i, xi, yi) + lim_monocentral(r_x[i])/2 * (us->at(i, xi+1, yi) - us->at(i, xi, yi));
        int_x->at(i, xi*2, yi) = us->at(i, xi, yi) - lim_monocentral(r_x[i])/2 * (us->at(i, xi+1, yi) - us->at(i, xi, yi));

        int_y->at(i, xi, yi*2+1) = us->at(i, xi, yi) + lim_monocentral(r_y[i])/2 * (us->at(i, xi, yi+1) - us->at(i, xi, yi));
        int_y->at(i, xi, yi*2) = us->at(i, xi, yi) - lim_monocentral(r_y[i])/2 * (us->at(i, xi, yi+1) - us->at(i, xi, yi));

    }

    //For the x and y components of the magnetic field (stored at +1/2):
    int_x->at(5, xi*2+1, yi) = us->at(5, xi, yi);
    int_x->at(5, xi*2, yi) = us->at(5, xi-1, yi);
    int_y->at(5, xi, yi*2+1) = (   (us->at(5, xi, yi)+us->at(5, xi-1, yi))/2 + (us->at(5, xi, yi+1)+us->at(5, xi-1, yi+1))/2   )/2;
    int_y->at(5, xi, yi*2) =   (   (us->at(5, xi, yi)+us->at(5, xi-1, yi))/2 + (us->at(5, xi, yi-1)+us->at(5, xi-1, yi-1))/2   )/2;

    int_x->at(6, xi*2+1, yi) = (   (us->at(6, xi, yi)+us->at(6, xi, yi-1))/2 + (us->at(6, xi+1, yi)+us->at(6, xi+1, yi-1))/2   )/2;
    int_x->at(6, xi*2, yi) =   (   (us->at(6, xi, yi)+us->at(6, xi, yi-1))/2 + (us->at(6, xi-1, yi)+us->at(6, xi-1, yi-1))/2   )/2;
    int_y->at(6, xi, yi*2+1) = us->at(6, xi, yi);
    int_y->at(6, xi, yi*2) = us->at(6, xi-1, yi);


    for (int i = 0; i < 8; i++) {
        
        int nans = isnan(int_x->at(i, xi*2, yi)) + isnan(int_x->at(i, xi*2+1, yi)) + isnan(int_y->at(i, xi, yi*2+1)) + isnan(int_y->at(i, xi, yi*2));

        if (nans != 0) {
            cout << "PLM:: There are " << nans << " NaNs at index " << i << endl;
            return false;
        }
        
    }
    
    

    if (int_x->at(0, xi*2+1, yi) <= 0 && int_x->at(0, xi*2, yi) <= 0) {
        cout << "PLM reconstruction: the density is negative" << endl;
    }


    return true;
}

bool PLMInterfaces_Prim(Arrays *us, int xi, int yi, Arrays *int_x, Arrays *int_y, long double dt) {
    //rho, Mx, My, Mz, E, Bx, By, Bz
    //No need to interpolate for Bx and By
    //cout << xi << ", " << yi << endl;

    long double r_x[8] = { 0 }; long double r_y[8] = { 0 };
    //TODO: When converting u to w, we're using rho, M, E, Bz - which are at cell centers - and Bx, By - which are at cell INTERFACES !!!
    //TODO: Change that to have Bx, By at cell centers.
    long double wi[8] = { 0 }; for (int i = 0; i < 8; i++) wi[i] = us->at(i, xi, yi);
    /*wi[5] = (us->at(5, xi, yi) + us->at(5, xi+1, yi))/2.0; wi[6] = (us->at(6, xi, yi) + us->at(6, xi, yi+1))/2.0;*/ toPrimitive(wi, wi);
    long double wixp[8] = { 0 }; for (int i = 0; i < 8; i++) wixp[i] = us->at(i, xi+1, yi);
    /*wixp[5] = (us->at(5, xi+1, yi) + us->at(5, xi+1, yi))/2.0; wixp[6] = (us->at(6, xi+1, yi) + us->at(6, xi+1, yi+1))/2.0;*/ toPrimitive(wixp, wixp);//FIXME: should be xi+2 on the first one but out of index
    long double wixm[8] = { 0 }; for (int i = 0; i < 8; i++) wixm[i] = us->at(i, xi-1, yi);
    /*wixm[5] = (us->at(5, xi-1, yi) + us->at(5, xi, yi))/2.0; wixm[6] = (us->at(6, xi-1, yi) + us->at(6, xi-1, yi+1))/2.0;*/ toPrimitive(wixm, wixm);
    long double wiyp[8] = { 0 }; for (int i = 0; i < 8; i++) wiyp[i] = us->at(i, xi, yi+1);
    /*wiyp[5] = (us->at(5, xi, yi+1) + us->at(5, xi+1, yi+1))/2.0; wiyp[6] = (us->at(6, xi, yi+1) + us->at(6, xi, yi+1))/2.0;*/ toPrimitive(wiyp, wiyp); //FIXME: should be yi+2 on the second one but out of index
    long double wiym[8] = { 0 }; for (int i = 0; i < 8; i++) wiym[i] = us->at(i, xi, yi-1);
    /*wiym[5] = (us->at(5, xi, yi-1) + us->at(5, xi+1, yi-1))/2.0; wiym[6] = (us->at(6, xi, yi-1) + us->at(6, xi, yi))/2.0;*/ toPrimitive(wiym, wiym);

    long double wx_ip[8] = { 0 }; long double wx_im[8] = { 0 };
    long double wy_ip[8] = { 0 }; long double wy_im[8] = { 0 };
    for (int i = 0; i < 8; i++)
    {
        if (i==5 || i==6) continue;

        long double r_x = 0; long double r_y = 0;
        if (wixp[i] != wi[i])
            r_x = (wi[i] - wixm[i])/(wixp[i] - wi[i]);
        if (wiyp[i] != wi[i])
            r_y = (wi[i] - wiym[i])/(wiyp[i] - wi[i]);

        wx_ip[i] = wi[i] + lim_minmod(r_x)/2 * (wixp[i] - wi[i]);
        wx_im[i] = wi[i] - lim_minmod(r_x)/2 * (wixp[i] - wi[i]);

        wy_ip[i] = wi[i] + lim_minmod(r_y)/2 * (wiyp[i] - wi[i]);
        wy_im[i] = wi[i] - lim_minmod(r_y)/2 * (wiyp[i] - wi[i]);

    }
    if (wx_ip[4] <= 0 || wx_im[4] <= 0 || wy_ip[4] <= 0 || wy_im[4] <= 0) {
        cout << "PLM Reconstruction: the pressure is negative or null" << endl;
        return false;
    }

    //For the x and y components of the magnetic field (stored at +1/2):
    wx_ip[5] = us->at(5, xi+1, yi);
    wx_im[5] = us->at(5, xi, yi);
    wy_ip[5] = ((us->at(5, xi+1, yi)+us->at(5, xi, yi))/2 + (us->at(5, xi+1, yi+1)+us->at(5, xi, yi+1))/2)/2;
    wy_im[5] = ((us->at(5, xi+1, yi)+us->at(5, xi, yi))/2 + (us->at(5, xi+1, yi-1)+us->at(5, xi, yi-1))/2)/2;

    wx_ip[6] = ((us->at(6, xi, yi+1)+us->at(6, xi, yi))/2 + (us->at(6, xi+1, yi+1)+us->at(6, xi+1, yi))/2)/2;
    wx_im[6] = ((us->at(6, xi, yi+1)+us->at(6, xi, yi))/2 + (us->at(6, xi-1, yi+1)+us->at(6, xi-1, yi))/2)/2;
    wy_ip[6] = us->at(6, xi, yi+1);
    wy_im[6] = us->at(6, xi, yi);

    toConserved(wx_ip, wx_ip); toConserved(wx_im, wx_im); toConserved(wy_ip, wy_ip); toConserved(wy_im, wy_im);

    //Check that everything's ok.
    long double temp[8] = { 0 };
    if (!toPrimitive(wx_ip, temp)) { cout << "PLM_PRIM:: wtf ?"; return false; } if (!toPrimitive(wx_im, temp)) { cout << "PLM_PRIM:: wtf ?"; return false; }
    if (!toPrimitive(wy_ip, temp)) { cout << "PLM_PRIM:: wtf ?"; return false; } if (!toPrimitive(wy_ip, temp)) { cout << "PLM_PRIM:: wtf ?"; return false; }
    for (int i = 0; i < 8; i++)
    {
        int_x->at(i, xi*2+1, yi) = wx_ip[i];
        int_x->at(i, xi*2, yi) = wx_im[i];
        int_y->at(i, xi, yi*2+1) = wy_ip[i];
        int_y->at(i, xi, yi*2) = wy_im[i];
    }


    return true;
}





//Compute interfaces at (i-1/2; j) and (i; j-1/2)
bool PPMInterfaces(Arrays *us, int xi, int yi, Arrays *int_x, Arrays *int_y) {
    //rho, Mx, My, Mz, E, Bx, By, Bz
    //No need to interpolate for Bx and By
    //cout << xi << ", " << yi << endl;

    long double r_x[8] = { 0 }; long double r_y[8] = { 0 };

    for (int i = 0; i < 8; i++)
    {
        if (i == 5 || i == 6) continue;

        //'r' for limiter

        if (us->at(i, xi+1, yi) != us->at(i, xi, yi))
            r_x[i] = (us->at(i, xi, yi) - us->at(i, xi-1, yi))/(us->at(i, xi+1, yi) - us->at(i, xi, yi));
        if (us->at(i, xi, yi+1) != us->at(i, xi, yi))
            r_y[i] = (us->at(i, xi, yi) - us->at(i, xi, yi-1))/(us->at(i, xi, yi+1) - us->at(i, xi, yi));


        //minus,plus interfaces for x,y
        long double duxp = us->at(i, xi+1, yi)-us->at(i, xi, yi); long double duxm = us->at(i, xi, yi)-us->at(i, xi-1, yi);
        long double duyp = us->at(i, xi, yi+1)-us->at(i, xi, yi); long double duym = us->at(i, xi, yi)-us->at(i, xi, yi-1);
        long double k = 1.0/3.0;

        int_x->at(i, xi*2+1, yi) = us->at(i, xi, yi) + lim_minmod(r_x[i])/4.0 * ((1-k)*duxm + (1+k)*duxp); //Lp
        int_x->at(i, xi*2, yi) = us->at(i, xi, yi) - lim_minmod(r_x[i])/4.0 * ((1-k)*duxp + (1+k)*duxm); //Rm

        int_y->at(i, xi, yi*2+1) = us->at(i, xi, yi) + lim_minmod(r_y[i])/4.0 * ((1-k)*duym + (1+k)*duyp); //Lp
        int_y->at(i, xi, yi*2) = us->at(i, xi, yi) - lim_minmod(r_y[i])/4.0 * ((1-k)*duyp + (1+k)*duym); //Rm

    }

    if (int_x->at(0, xi*2+1, yi) <= 0 || int_x->at(0, xi*2, yi) <= 0  || int_y->at(0, xi, yi*2+1) <= 0 || int_y->at(0, xi, yi*2) <= 0) {
        cout << "PLM reconstruction: the density is negative" << endl;
        return false;
    }


    return true;
}

//Compute interfaces at (i-1/2; j) and (i; j-1/2)
bool PPMInterfaces_Prim(Arrays *us, int xi, int yi, Arrays *int_x, Arrays *int_y) {
    //rho, Mx, My, Mz, E, Bx, By, Bz
    //No need to interpolate for Bx and By
    //cout << xi << ", " << yi << endl;

    long double r_x[8] = { 0 }; long double r_y[8] = { 0 };
    long double wi[8] = { 0 }; for (int i = 0; i < 8; i++) wi[i] = us->at(i, xi, yi); toPrimitive(wi, wi);
    long double wip[8] = { 0 }; for (int i = 0; i < 8; i++) wip[i] = us->at(i, xi+1, yi); toPrimitive(wip, wip);
    long double wim[8] = { 0 }; for (int i = 0; i < 8; i++) wim[i] = us->at(i, xi-1, yi); toPrimitive(wim, wim);

    long double wx_ip[8] = { 0 }; long double wx_im[8] = { 0 };
    for (int i = 0; i < 8; i++)
    {
        if (i==5 || i==6) continue;
        long double r_x = 0;
        if (wip[i] != wi[i])
            r_x = (wi[i] - wim[i])/(wip[i] - wi[i]);

        long double dwxp = wim[i]-wi[i]; long double dwxm = wi[i]-wim[i];
        long double k = 1.0/3.0;

        wx_ip[i] = wi[i] + lim_minmod(r_x)/4.0 * ((1-k)*dwxm + (1+k)*dwxp); //Lp
        wx_im[i] = wi[i] - lim_minmod(r_x)/4.0 * ((1-k)*dwxp + (1+k)*dwxm); //Rm

    }
    if (wx_ip[4] <= 0 || wx_im[4] <= 0) {
        cout << "PLM Reconstruction: the pressure is negative" << endl;
        return false;
    }
    toConserved(wx_ip, wx_ip); toConserved(wx_im, wx_im);
    for (int i = 0; i < 8; i++)
    {
        int_x->at(i, xi*2+1, yi) = wx_ip[i];
        int_x->at(i, xi*2, yi) = wx_im[i];
    }


    for (int i = 0; i < 8; i++)
    {
        if (i == 5 || i == 6) continue;

        //'r' for limiter
        if (us->at(i, xi, yi+1) != us->at(i, xi, yi))
            r_y[i] = (us->at(i, xi, yi) - us->at(i, xi, yi-1))/(us->at(i, xi, yi+1) - us->at(i, xi, yi));


        //minus,plus interfaces for x,y
        long double duyp = us->at(i, xi, yi+1)-us->at(i, xi, yi); long double duym = us->at(i, xi, yi)-us->at(i, xi, yi-1);
        long double k = 1.0/3.0;

        int_y->at(i, xi, yi*2+1) = us->at(i, xi, yi) + lim_minmod(r_y[i])/4.0 * ((1-k)*duym + (1+k)*duyp); //Lp
        int_y->at(i, xi, yi*2) = us->at(i, xi, yi) - lim_minmod(r_y[i])/4.0 * ((1-k)*duyp + (1+k)*duym); //Rm

    }

    if (int_y->at(0, xi, yi*2+1) <= 0 || int_y->at(0, xi, yi*2) <= 0) {
        cout << "PLM reconstruction: the density is negative" << endl;
        return false;
    }


    return true;
}