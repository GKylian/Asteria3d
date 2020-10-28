#pragma once
#include <math.h>
#include <iostream>
#include "Arrays.h"
#include "MHD.h"
#include "Eigenvalues.h"


bool Vortex(Arrays *us, ld x0, ld y0, ld xn, ld yn) {
    int Nx = us->Nx; int Ny = us->Ny;
    ld dx = us->dx; ld dy = us->dy;

    for (int yi = 0; yi < Ny; yi++)
        for (int xi = 0; xi < Nx; xi++)
        {
            ld x = x0 + (xi-4.0)*dx; if (x < x0) x = x0; if (x > xn) x = xn;
            ld y = y0 + (yi-4.0)*dy; if (y < y0) y = y0; if (y > yn) y = yn;
            ld x_ix = x0 + (xi-4.5)*dx; //x at the x interface (where Bx is)
            ld y_iy = y0 + (yi-4.5)*dy; //y at the y interface (where By is)


            ld B0 = 1/M_rt4PI;
            //B0 = 0.0;

            ld a[2][2] = { 0 };

            for (int i = 0; i < 2; i++) //Computes the potential at corner (i-1/2, j-1/2)
            for (int j = 0; j < 2; j++)
            {
                int xii = i + xi; int yii = j + yi;
                ld x_ix = x0 + (xii-4.5)*dx; //Coordinates at the x and y interfaces
                ld y_iy = y0 + (yii-4.5)*dy;


                a[i][j] = B0*(cosl(4*M_PI*x_ix)/(4*M_PI)+cosl(2*M_PI*y_iy)/(2*M_PI));
            }

            us->at(5, xi, yi) = (a[0][1]-a[0][0])/dy; us->at(6, xi, yi) = -(a[1][0]-a[0][0])/dx;





            ld rho = 25.0/(36.0*M_PI);
            us->at(0, xi, yi) = rho; us->at(1, xi, yi) = -rho*sinl(2.0*M_PI*y); us->at(2, xi, yi) = rho*sinl(2.0*M_PI*x); us->at(3, xi, yi) = rho*0.0;
            //us->at(5, xi, yi) = -B0*sinl(2*M_PI*y); us->at(6, xi, yi) = B0*sinl(4*M_PI*x); us->at(7, xi, yi) = 0.0;
            

            ld P = 5.0/(12.0*M_PI);
            ld v = sqrtl(powl(us->at(1, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(2, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(3, xi, yi)/us->at(0, xi, yi), 2));
            ld B = sqrtl(powl(us->at(5, xi, yi), 2) + powl(us->at(6, xi, yi), 2) + powl(us->at(7, xi, yi), 2));

            us->at(4, xi, yi) = P;// / (gamma - 1) + us->at(0, xi, yi)*v*v/2 + B*B/2;

        }

    //Do the transformation P->E afterwards to avoid using face-centered Bx and By values as cell-centered ones.
    setBoundaries(bound::PERIODIC, bound::PERIODIC, us, Nx-8, Ny-8, 0);

    for (int yi = 0; yi < Ny; yi++)
    for (int xi = 0; xi < Nx; xi++)
    {
        ld P = us->at(4, xi, yi);
        ld Bx = (us->at(5, xi+1, yi)+us->at(5, xi, yi))/2.0;
        ld By = (us->at(6, xi, yi+1)+us->at(6, xi, yi))/2.0;

        ld v2 = powl(us->at(1, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(2, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(3, xi, yi)/us->at(0, xi, yi), 2);
        ld B2 = Bx*Bx + By*By + powl(us->at(7, xi, yi), 2);
        us->at(4, xi, yi) = P / (gamma - 1) + us->at(0, xi, yi)*v2/2 + B2/2;
            

        for (int i = 0; i < 8; i++)
        {
            if (isnan(us->at(i, xi, yi))) { cout << i << endl; return false; }
        }
    }

    return true;
}


bool RandomTest(Arrays *us, ld x0, ld y0, ld xn, ld yn) {
    int Nx = us->Nx; int Ny = us->Ny;

    for (int yi = 0; yi < Ny; yi++)
    for (int xi = 0; xi < Nx; xi++)
    {
        for (int i = 0; i < 8; i++) {
            us->at(i, xi, yi) = rand() % 10 + 1;
        }
    }
    return true;
}

bool Loops(Arrays* us, ld x0, ld y0, ld xn, ld yn) {
    int Nx = us->Nx; int Ny = us->Ny;
    ld dx = us->dx; ld dy = us->dy;

    for (int yi = 0; yi < Ny; yi++)
        for (int xi = 0; xi < Nx; xi++)
        {
            ld x = x0 + (xi-4.0)*dx;// if (x < x0) x = x0; if (x > xn) x = xn;
            ld y = y0 + (yi-4.0)*dy;// if (y < y0) y = y0; if (y > yn) y = yn;
            ld r = sqrtl(x*x+y*y);

            ld R0 = 0.3; ld A0 = 0.001;
            ld a[2][2] = { 0 };

            for (int i = 0; i < 2; i++) //Computes the potential at corner (i-1/2, j-1/2)
            for (int j = 0; j < 2; j++)
            {
                int xii = i + xi; int yii = j + yi;
                ld x_ix = x0 + (xii-4.5)*dx; //Coordinates at the x and y interfaces
                ld y_iy = y0 + (yii-4.5)*dy;

                ld rc = sqrtl(x_ix*x_ix+y_iy*y_iy);

                if (rc < R0)
                    a[i][j] = A0*(R0-rc);
                else
                    a[i][j] = 0.0;
            }

            us->at(5, xi, yi) = (a[0][1]-a[0][0])/dy; us->at(6, xi, yi) = -(a[1][0]-a[0][0])/dx;




            if (r < R0) {
                ld rho = 1.0;
                us->at(0, xi, yi) = rho; us->at(1, xi, yi) = 2.0*rho; us->at(2, xi, yi) = 1.0*rho; us->at(3, xi, yi) = 0.0*rho;
                

                ld P = 1.0;
                ld v = sqrtl(powl(us->at(1, xi, yi)/us->at(0, xi, yi),2) + powl(us->at(2, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(3, xi, yi)/us->at(0, xi, yi), 2));
                ld B = sqrtl(powl(us->at(5, xi, yi),2) + powl(us->at(6, xi, yi), 2) + powl(us->at(7, xi, yi), 2));

                us->at(4, xi, yi) = P;// / (gamma - 1) + us->at(0, xi, yi)*v*v/2 + B*B/2;
            }
            else {
                ld rho = 1.0;
                us->at(0, xi, yi) = rho; us->at(1, xi, yi) = 2.0*rho; us->at(2, xi, yi) = 1.0*rho; us->at(3, xi, yi) = 0.0;
                
                ld P = 1.0;
                ld v = sqrtl(powl(us->at(1, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(2, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(3, xi, yi)/us->at(0, xi, yi), 2));
                ld B = sqrtl(powl(us->at(5, xi, yi), 2) + powl(us->at(6, xi, yi), 2) + powl(us->at(7, xi, yi), 2));

                us->at(4, xi, yi) = P;// / (gamma - 1) + us->at(0, xi, yi)*v*v/2 + B*B/2;
            }

        }

    setBoundaries(bound::PERIODIC, bound::PERIODIC, us, Nx-8, Ny-8, 0);

    for (int yi = 0; yi < Ny; yi++)
        for (int xi = 0; xi < Nx; xi++)
        {
            ld P = us->at(4, xi, yi);
            ld Bx = (us->at(5, xi+1, yi)+us->at(5, xi, yi))/2.0;
            ld By = (us->at(6, xi, yi+1)+us->at(6, xi, yi))/2.0;

            ld v2 = powl(us->at(1, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(2, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(3, xi, yi)/us->at(0, xi, yi), 2);
            ld B2 = Bx*Bx + By*By + powl(us->at(7, xi, yi), 2);
            us->at(4, xi, yi) = P / (gamma - 1) + us->at(0, xi, yi)*v2/2 + B2/2;
            

            for (int i = 0; i < 8; i++)
            {
                if (isnan(us->at(i, xi, yi))) { cout << i << endl; return false; }
            }
        }
    return true;
}



bool Sod(Arrays *us, ld x0, ld y0, ld xn, ld yn, bool hd) {
    int Nx = us->Nx; int Ny = us->Ny;
    ld dx = us->dx; ld dy = us->dy;

    for (int yi = 0; yi < Ny; yi++)
        for (int xi = 0; xi < Nx; xi++)
        {
            ld x = x0 + (xi-4.0)*dx;// if (x < x0) x = x0; if (x > xn) x = xn;
            ld y = y0 + (yi-4.0)*dy;// if (y < y0) y = y0; if (y > yn) y = yn;


            if (x < 0.5) {
                ld rho = 1.0;
                us->at(0, xi, yi) = rho; us->at(1, xi, yi) = 0.0*rho; us->at(2, xi, yi) = 0.0*rho; us->at(3, xi, yi) = 0.0*rho;

                if (Nx > 16) {
                    us->at(5, xi, yi) = 0.75; us->at(6, xi, yi) = 1.0; us->at(7, xi, yi) = 0.0;
                }
                else {
                    us->at(5, xi, yi) = 1.0; us->at(6, xi, yi) = 0.75; us->at(7, xi, yi) = 0.0;
                }

                if (hd) {
                    us->at(5, xi, yi) = 0.0; us->at(6, xi, yi) = 0.0; us->at(7, xi, yi) = 0.0;
                }
                

                ld P = 1.0;
                ld v = sqrtl(powl(us->at(1, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(2, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(3, xi, yi)/us->at(0, xi, yi), 2));
                ld B = sqrtl(powl(us->at(5, xi, yi), 2) + powl(us->at(6, xi, yi), 2) + powl(us->at(7, xi, yi), 2));

                us->at(4, xi, yi) = P;
            }
            else {
                ld rho = 0.125;
                us->at(0, xi, yi) = rho; us->at(1, xi, yi) = 0.0*rho; us->at(2, xi, yi) = 0.0*rho; us->at(3, xi, yi) = 0.0*rho; //vz 1.0
                if (Nx > 16) {
                    us->at(5, xi, yi) = 0.75; us->at(6, xi, yi) = -1.0; us->at(7, xi, yi) = 0.0;
                }
                else {
                    us->at(5, xi, yi) = -1.0; us->at(6, xi, yi) = 0.75; us->at(7, xi, yi) = 0.0;
                }

                if (hd) {
                    us->at(5, xi, yi) = 0.0; us->at(6, xi, yi) = 0.0; us->at(7, xi, yi) = 0.0;
                }

                ld P = 0.1;
                ld v = sqrtl(powl(us->at(1, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(2, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(3, xi, yi)/us->at(0, xi, yi), 2));
                ld B = sqrtl(powl(us->at(5, xi, yi), 2) + powl(us->at(6, xi, yi), 2) + powl(us->at(7, xi, yi), 2));

                us->at(4, xi, yi) = P;
            }

        }

    setBoundaries(bound::OUTFLOW, bound::OUTFLOW, us, Nx-8, Ny-8, 0);

    for (int yi = 0; yi < Ny; yi++)
    for (int xi = 0; xi < Nx; xi++)
    {
        ld P = us->at(4, xi, yi);
        ld Bx = (us->at(5, xi+1, yi)+us->at(5, xi, yi))/2.0;
        ld By = (us->at(6, xi, yi+1)+us->at(6, xi, yi))/2.0;

        ld v2 = powl(us->at(1, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(2, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(3, xi, yi)/us->at(0, xi, yi), 2);
        ld B2 = Bx*Bx + By*By + powl(us->at(7, xi, yi), 2);
        us->at(4, xi, yi) = P / (gamma - 1) + us->at(0, xi, yi)*v2/2 + B2/2;

    }



    return true;
}



ld VP(ld R, ld r0, ld q) {
    return sqrtl(r0)/(r0-2.0)*powl(R/r0, 1.0-q);
}


bool Disk(Arrays *us, ld x0, ld y0, ld rs) {
    int Nx = us->Nx; int Ny = us->Ny;
    ld dx = us->dx; ld dy = us->dy;

    //Disk parameters:
    double q = 1.68; double r0 = 20.0; double rhomax = 300.0; double rin = 8.0; double rho0 = 0.01; double e0 = 0.0001;
    double dcut = 1.0; double beta = 100.0; double seed = 0.01;
    ld n = 1/(gamma-1); ld q1 = 2*(q-1); ld r02 = r0-rs;
    ld f = powl(r0, 2*q-1)/(q1*r02*r02); ld C = f*powl(rin, -q1)-1.0/(rin-rs); ld K = (C+1/r02-f*powl(r0, -q1))/powl(rhomax, 1.0/n);

    for (int yi = 0; yi < Ny; yi++)
        for (int xi = 0; xi < Nx; xi++)
        {
            ld R = x0 + (xi-4.0)*dx;
            ld z = y0 + (yi-4.0)*dy;
            ld R_iR = x0 + (xi-4.5)*dx; //x at the x interface (where Bx is)
            ld z_iz = y0 + (yi-4.5)*dy; //y at the y interface (where By is)
            ld r = sqrtl(R*R+z*z);

            us->at(0, xi, yi) = rho0; us->at(3, xi, yi) = VP(R, r0, q)*rho0; us->at(4, xi, yi) = e0;

            ld Ta = f*powl(R, -q1); ld T = 1.0/Ta + C;
            ld cs = 0.086;
            if (R >= rin && r <= T) {
                ld h = fmaxl(C + 1.0/(r-2.0) - f*powl(R, -2.0*(q-1.0)), 0);
                ld rho = powl(h/K, n)*(R >= rin); rho = fmaxl(rho, rho0);
                ld Eint = powl(rho, gamma)*K/((gamma-1)*(n+1)); Eint *= (1-0.01*sinl(R));
                if (Eint >= e0 && rho >= rho0) {
                    us->at(0, xi, yi) = rho; us->at(1, xi, yi) = rho*0.0; us->at(2, xi, yi) = rho*0.0; us->at(3, xi, yi) = rho*VP(R, r0, q);
                    us->at(4, xi, yi) = Eint; //a^2 = gamma*P/rho or P = cs^2 rho
                }
                

            }
        }

    setBoundaries(bound::OUTFLOW, bound::OUTFLOW, us, Nx-8, Ny-8, 0);

    //for (int yi = 0; yi < Ny; yi++)
    //    for (int xi = 0; xi < Nx; xi++)
    //    {

    //        ld a[2][2] = { 0 };

    //        for (int i = 0; i < 2; i++) //Computes the potential at corner (i-1/2, j-1/2)
    //            for (int j = 0; j < 2; j++)
    //            {
    //                int xii = i + xi; int yii = j + yi;
    //                ld x_ix = x0 + (xii-4.5)*dx; //Coordinates at the x and y interfaces
    //                ld y_iy = y0 + (yii-4.5)*dy;


    //                a[i][j] = us->at(4, xi, yi)-0.1;
    //            }

    //        us->at(5, xi, yi) = (a[0][1]-a[0][0])/dy; us->at(6, xi, yi) = -(a[1][0]-a[0][0])/dx;



    //    }

    //Do the transformation P->E afterwards to avoid using face-centered Bx and By values as cell-centered ones.
    //setBoundaries(bound::OUTFLOW, bound::PERIODIC, us, Nx-8, Ny-8, 0);

    for (int yi = 0; yi < Ny; yi++)
        for (int xi = 0; xi < Nx; xi++)
        {
            ld P = us->at(4, xi, yi);
            ld Bx = (us->at(5, xi+1, yi)+us->at(5, xi, yi))/2.0;
            ld By = (us->at(6, xi, yi+1)+us->at(6, xi, yi))/2.0;

            ld v2 = powl(us->at(1, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(2, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(3, xi, yi)/us->at(0, xi, yi), 2);
            ld B2 = Bx*Bx + By*By + powl(us->at(7, xi, yi), 2);
            us->at(4, xi, yi) = P / (gamma - 1) + us->at(0, xi, yi)*v2/2 + B2/2;


            for (int i = 0; i < 8; i++)
            {
                if (isnan(us->at(i, xi, yi))) { cout << i << endl; return false; }
            }
        }

    return true;
}




ld rotor_f(ld r, ld r0, ld r1) {
    return (r1-r)/(r1-r0);
}

bool MHDRotor(Arrays* us, ld x0, ld y0) {
    int Nx = us->Nx; int Ny = us->Ny;
    ld dx = us->dx; ld dy = us->dy;
    cout << "Initializing the MHD Rotor problem with resolution " << Nx << "x" << Ny << endl;

    for (int yi = 0; yi < Ny; yi++)
        for (int xi = 0; xi < Nx; xi++)
        {
            //cout << xi << ", " << yi << endl;
            ld x = x0 + (xi-4.0)*dx;// if (x < x0) x = x0; if (x > xn) x = xn;
            ld y = y0 + (yi-4.0)*dy;// if (y < y0) y = y0; if (y > yn) y = yn;
            ld r = sqrtl((x-0.5)*(x-0.5) + (y-0.5)*(y-0.5));
            

            ld r0 = 0.1; ld r1 = 0.115; ld v0 = 1.0; ld P = 0.5;
            if (r <= r0) {
                ld rho = 10.0;
                us->at(0, xi, yi) = rho; us->at(1, xi, yi) = -rotor_f(r, r0, r1)*v0*(y-0.5)/r0 * rho; us->at(2, xi, yi) = rotor_f(r, r0, r1)*v0*(x-0.5)/r0 * rho; us->at(3, xi, yi) = 0.0*rho;
                us->at(5, xi, yi) = 2.5/sqrtl(4*M_PI); us->at(6, xi, yi) = 0.0; us->at(7, xi, yi) = 0.0;
                
                
                ld v = sqrtl(powl(us->at(1, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(2, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(3, xi, yi)/us->at(0, xi, yi), 2));
                ld B = sqrtl(powl(us->at(5, xi, yi), 2) + powl(us->at(6, xi, yi), 2) + powl(us->at(7, xi, yi), 2));

                us->at(4, xi, yi) = P / (gamma - 1) + us->at(0, xi, yi)*v*v/2 + B*B/2;
            }
            if (r > r0 && r < r1){
                ld rho = 1.0 + 9.0*rotor_f(r, r0, r1);
                us->at(0, xi, yi) = rho; us->at(1, xi, yi) = -rotor_f(r, r0, r1)*v0*(y-0.5)/r * rho; us->at(2, xi, yi) = rotor_f(r, r0, r1)*v0*(x-0.5)/r * rho; us->at(3, xi, yi) = 0.0*rho;
                us->at(5, xi, yi) = 2.5/sqrtl(4*M_PI); us->at(6, xi, yi) = 0.0; us->at(7, xi, yi) = 0.0;

                
                ld v = sqrtl(powl(us->at(1, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(2, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(3, xi, yi)/us->at(0, xi, yi), 2));
                ld B = sqrtl(powl(us->at(5, xi, yi), 2) + powl(us->at(6, xi, yi), 2) + powl(us->at(7, xi, yi), 2));

                us->at(4, xi, yi) = P / (gamma - 1) + us->at(0, xi, yi)*v*v/2 + B*B/2;
            }
            if (r >= r1) {
                ld rho = 1.0;
                us->at(0, xi, yi) = rho; us->at(1, xi, yi) = 0.0*rho; us->at(2, xi, yi) = 0.0*rho; us->at(3, xi, yi) = 0.0;
                us->at(5, xi, yi) = 2.5/sqrtl(4*M_PI); us->at(6, xi, yi) = 0.0; us->at(7, xi, yi) = 0.0; //5.0/M_rt4PI

                
                ld v = sqrtl(powl(us->at(1, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(2, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(3, xi, yi)/us->at(0, xi, yi), 2));
                ld B = sqrtl(powl(us->at(5, xi, yi), 2) + powl(us->at(6, xi, yi), 2) + powl(us->at(7, xi, yi), 2));

                us->at(4, xi, yi) = P / (gamma - 1) + us->at(0, xi, yi)*v*v/2 + B*B/2;
            }

            ld u[8] = { 0 }; for (int i = 0; i < 8; i++) u[i] = us->at(i, xi, yi);
            ld w[8] = { 0 }; if (!toPrimitive(u, w)) { cout << "Init failed" << endl; return false; }
            /*for (int i = 0; i < 8; i++) {
                if (fabsl(u[i]) > 100) {
                    cout << "u[" << i << "] (" << xi << ", " << yi << ") > 100 !!!" << endl;
                    cout << "u = (" << u[0] << ", " << u[1] << ", " << u[2] << ", " << u[3] << ", " << u[4] << ")." << endl;
                    cout << "r = " << r << endl;
                    return false;
                    
                }
            }*/
            

        }

    return true;
}



bool Sedov(Arrays *us, ld x0, ld y0) {
    int Nx = us->Nx; int Ny = us->Ny;
    ld dx = us->dx; ld dy = us->dy;

    for (int yi = 0; yi < Ny; yi++)
        for (int xi = 0; xi < Nx; xi++)
        {
            ld x = x0 + (xi-4.0)*dx;// if (x < x0) x = x0; if (x > xn) x = xn;
            ld y = y0 + (yi-4.0)*dy;// if (y < y0) y = y0; if (y > yn) y = yn;
            ld r = sqrtl(x*x+y*y);

            ld dr = 3.5*dx; ld rho = 1.0; ld nu = 1.0;
            if (r < dr) {
                us->at(0, xi, yi) = rho; us->at(1, xi, yi) = 0.0*rho; us->at(2, xi, yi) = 0.0*rho; us->at(3, xi, yi) = 0.0*rho;
                us->at(5, xi, yi) = 0.0; us->at(6, xi, yi) = 0.0; us->at(7, xi, yi) = 0.0;


                ld P = 3*(gamma-1)*1.0/((nu+1)*M_PI*powl(dr, nu));
                ld v = sqrtl(powl(us->at(1, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(2, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(3, xi, yi)/us->at(0, xi, yi), 2));
                ld B = sqrtl(powl(us->at(5, xi, yi), 2) + powl(us->at(6, xi, yi), 2) + powl(us->at(7, xi, yi), 2));

                us->at(4, xi, yi) = P / (gamma - 1) + us->at(0, xi, yi)*v*v/2 + B*B/2;
                //us->at(4, xi, yi) = 1.0;
            }
            else {
                us->at(0, xi, yi) = rho; us->at(1, xi, yi) = 0.0*rho; us->at(2, xi, yi) = 0.0*rho; us->at(3, xi, yi) = 0.0;
                us->at(5, xi, yi) = 0.0; us->at(6, xi, yi) = 0.0; us->at(7, xi, yi) = 0.0;

                ld P = 1e-5;
                ld v = sqrtl(powl(us->at(1, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(2, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(3, xi, yi)/us->at(0, xi, yi), 2));
                ld B = sqrtl(powl(us->at(5, xi, yi), 2) + powl(us->at(6, xi, yi), 2) + powl(us->at(7, xi, yi), 2));

                us->at(4, xi, yi) = P / (gamma - 1) + us->at(0, xi, yi)*v*v/2 + B*B/2;
            }

            ld u[8] = { 0 }; for (int i = 0; i < 8; i++) u[i] = us->at(i, xi, yi);
            ld w[8] = { 0 }; if (!toPrimitive(u, w)) { cout << "Init failed" << endl; return false; }

        }

    cout << "Initialized the variables for Sedov test." << endl;
    return true;
}





bool RTInstability(Arrays *us, ld x0, ld y0, ld xn, ld yn, ld gy) {
    int Nx = us->Nx; int Ny = us->Ny;
    ld dx = us->dx; ld dy = us->dy;

    for (int yi = 0; yi < Ny; yi++)
        for (int xi = 0; xi < Nx; xi++)
        {
            ld x = x0 + (xi-4.0)*dx;// if (x < x0) x = x0; if (x > xn) x = xn;
            ld y = y0 + (yi-4.0)*dy;// if (y < y0) y = y0; if (y > yn) y = yn;


            ld kx = 2*M_PI/(xn-x0);// if (x0 == 0.0 || xn == 0.0) kx /= 2.0;
            ld ky = 2*M_PI/(yn-y0);
            ld A = 0.01;
            
            if (y > 0.0) { //0.025*cosl(2*M_PI*x)
                ld rho = 2.0; //1.08
                us->at(0, xi, yi) = rho; us->at(1, xi, yi) = 0.0*rho; us->at(2, xi, yi) = 0.01*(1+cosl(kx*x))*(1+cosl(ky*y))/4.0*rho; us->at(3, xi, yi) = 0.0*rho;
                //us->at(0, xi, yi) = rho; us->at(1, xi, yi) = 0.0*rho; us->at(2, xi, yi) = 0.0*rho; us->at(3, xi, yi) = 0.0*rho;
                us->at(5, xi, yi) = 0.0; us->at(6, xi, yi) = 0.0; us->at(7, xi, yi) = 0.0;

                //if (fabsl(y) < 2*dx) us->at(2, xi, yi) = A*sinl(kx*x)*rho;


                ld P = 1.0/gamma +gy*rho*y;

                ld v = sqrtl(powl(us->at(1, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(2, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(3, xi, yi)/us->at(0, xi, yi), 2));
                ld B = sqrtl(powl(us->at(5, xi, yi), 2) + powl(us->at(6, xi, yi), 2) + powl(us->at(7, xi, yi), 2));

                us->at(4, xi, yi) = P;
            }
            else {
                ld rho = 1.0; // 1.0
                us->at(0, xi, yi) = rho; us->at(1, xi, yi) = 0.0*rho; us->at(2, xi, yi) = 0.01*(1+cosl(kx*x))*(1+cosl(ky*y))/4.0*rho; us->at(3, xi, yi) = 0.0*rho; //vz 1.0
                //us->at(0, xi, yi) = rho; us->at(1, xi, yi) = 0.0*rho; us->at(2, xi, yi) = 0.0*rho; us->at(3, xi, yi) = 0.0*rho;
                us->at(5, xi, yi) = 0.0; us->at(6, xi, yi) = 0.0; us->at(7, xi, yi) = 0.0;

                //if (fabsl(y) < 2*dx) us->at(2, xi, yi) = A*sinl(kx*x)*rho;

                ld P = 1.0/gamma +gy*rho*y;
                
                ld v = sqrtl(powl(us->at(1, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(2, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(3, xi, yi)/us->at(0, xi, yi), 2));
                ld B = sqrtl(powl(us->at(5, xi, yi), 2) + powl(us->at(6, xi, yi), 2) + powl(us->at(7, xi, yi), 2));

                us->at(4, xi, yi) = P;
            }

        }

    setBoundaries(bound::PERIODIC, bound::PRESSURE, us, Nx-8, Ny-8, gy);

    for (int yi = 0; yi < Ny; yi++)
    for (int xi = 0; xi < Nx; xi++)
    {
        ld P = us->at(4, xi, yi);
        ld Bx = (us->at(5, xi+1, yi)+us->at(5, xi, yi))/2.0;
        ld By = (us->at(6, xi, yi+1)+us->at(6, xi, yi))/2.0;

        ld v2 = powl(us->at(1, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(2, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(3, xi, yi)/us->at(0, xi, yi), 2);
        ld B2 = Bx*Bx + By*By + powl(us->at(7, xi, yi), 2);
        us->at(4, xi, yi) = P / (gamma - 1) + us->at(0, xi, yi)*v2/2 + B2/2;

    }



    return true;
}




bool MHDRTInstability(Arrays *us, ld x0, ld y0, ld xn, ld yn, ld gy) {
    int Nx = us->Nx; int Ny = us->Ny;
    ld dx = us->dx; ld dy = us->dy;

    for (int yi = 0; yi < Ny; yi++)
        for (int xi = 0; xi < Nx; xi++)
        {
            ld x = x0 + (xi-4.0)*dx;// if (x < x0) x = x0; if (x > xn) x = xn;
            ld y = y0 + (yi-4.0)*dy;// if (y < y0) y = y0; if (y > yn) y = yn;


            ld kx = 2*M_PI/(xn-x0);// if (x0 == 0.0 || xn == 0.0) kx /= 2.0;
            ld ky = 2*M_PI/(yn-y0);
            ld A = 0.01;
            ld B0 = 0.0;
            B0 = 0.05*0.14;
            
            if (y > 0.0) { //0.025*cosl(2*M_PI*x)
                ld rho = 3.0; //1.08
                us->at(0, xi, yi) = rho; us->at(1, xi, yi) = 0.0*rho; us->at(2, xi, yi) = 0.01*(1+cosl(kx*x))*(1+cosl(ky*y))/4.0*rho; us->at(3, xi, yi) = 0.0*rho;
                //us->at(0, xi, yi) = rho; us->at(1, xi, yi) = 0.0*rho; us->at(2, xi, yi) = 0.0*rho; us->at(3, xi, yi) = 0.0*rho;
                us->at(5, xi, yi) = B0; us->at(6, xi, yi) = 0.0; us->at(7, xi, yi) = 0.0;

                //if (fabsl(y) < 2*dx) us->at(2, xi, yi) = A*sinl(kx*x)*rho;


                ld P = 1.0/gamma +gy*rho*y;

                ld v = sqrtl(powl(us->at(1, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(2, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(3, xi, yi)/us->at(0, xi, yi), 2));
                ld B = sqrtl(powl(us->at(5, xi, yi), 2) + powl(us->at(6, xi, yi), 2) + powl(us->at(7, xi, yi), 2));

                us->at(4, xi, yi) = P;
            }
            else {
                ld rho = 1.0; // 1.0
                us->at(0, xi, yi) = rho; us->at(1, xi, yi) = 0.0*rho; us->at(2, xi, yi) = 0.01*(1+cosl(kx*x))*(1+cosl(ky*y))/4.0*rho; us->at(3, xi, yi) = 0.0*rho; //vz 1.0
                //us->at(0, xi, yi) = rho; us->at(1, xi, yi) = 0.0*rho; us->at(2, xi, yi) = 0.0*rho; us->at(3, xi, yi) = 0.0*rho;
                us->at(5, xi, yi) = B0; us->at(6, xi, yi) = 0.0; us->at(7, xi, yi) = 0.0;

                //if (fabsl(y) < 2*dx) us->at(2, xi, yi) = A*sinl(kx*x)*rho;

                ld P = 1.0/gamma +gy*rho*y;
                
                ld v = sqrtl(powl(us->at(1, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(2, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(3, xi, yi)/us->at(0, xi, yi), 2));
                ld B = sqrtl(powl(us->at(5, xi, yi), 2) + powl(us->at(6, xi, yi), 2) + powl(us->at(7, xi, yi), 2));

                us->at(4, xi, yi) = P;
            }

        }

    setBoundaries(bound::PERIODIC, bound::PRESSURE, us, Nx-8, Ny-8, gy);

    for (int yi = 0; yi < Ny; yi++)
    for (int xi = 0; xi < Nx; xi++)
    {
        ld P = us->at(4, xi, yi);
        ld Bx = (us->at(5, xi+1, yi)+us->at(5, xi, yi))/2.0;
        ld By = (us->at(6, xi, yi+1)+us->at(6, xi, yi))/2.0;

        ld v2 = powl(us->at(1, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(2, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(3, xi, yi)/us->at(0, xi, yi), 2);
        ld B2 = Bx*Bx + By*By + powl(us->at(7, xi, yi), 2);
        us->at(4, xi, yi) = P / (gamma - 1) + us->at(0, xi, yi)*v2/2 + B2/2;

    }



    return true;
}






bool HDEq(Arrays *us, ld x0, ld y0, ld xn, ld yn, ld gy) {
    int Nx = us->Nx; int Ny = us->Ny;
    ld dx = us->dx; ld dy = us->dy;

    for (int yi = 0; yi < Ny; yi++)
        for (int xi = 0; xi < Nx; xi++)
        {
            ld x = x0 + (xi-4.0)*dx;// if (x < x0) x = x0; if (x > xn) x = xn;
            ld y = y0 + (yi-4.0)*dy;// if (y < y0) y = y0; if (y > yn) y = yn;

            
            ld rho = 2.0; //1.08
            //us->at(0, xi, yi) = rho; us->at(1, xi, yi) = 0.0*rho; us->at(2, xi, yi) = 0.01*(1+cosl(kx*x))*(1+cosl(ky*y))/4.0*rho; us->at(3, xi, yi) = 0.0*rho;
            us->at(0, xi, yi) = rho; us->at(1, xi, yi) = 0.0*rho; us->at(2, xi, yi) = 0.0*rho; us->at(3, xi, yi) = 0.0*rho;
            us->at(5, xi, yi) = 0.0; us->at(6, xi, yi) = 0.0; us->at(7, xi, yi) = 0.0;

            //if (fabsl(y) < 2*dx) us->at(2, xi, yi) = A*sinl(kx*x)*rho;


            ld P = 1.0/gamma;// +gy*rho*y;

            ld v = sqrtl(powl(us->at(1, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(2, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(3, xi, yi)/us->at(0, xi, yi), 2));
            ld B = sqrtl(powl(us->at(5, xi, yi), 2) + powl(us->at(6, xi, yi), 2) + powl(us->at(7, xi, yi), 2));

            us->at(4, xi, yi) = P;
            

        }

    setBoundaries(bound::PERIODIC, bound::REFLECTING, us, Nx-8, Ny-8, gy);

    for (int yi = 0; yi < Ny; yi++)
    for (int xi = 0; xi < Nx; xi++)
    {
        ld P = us->at(4, xi, yi);
        ld Bx = (us->at(5, xi+1, yi)+us->at(5, xi, yi))/2.0;
        ld By = (us->at(6, xi, yi+1)+us->at(6, xi, yi))/2.0;

        ld v2 = powl(us->at(1, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(2, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(3, xi, yi)/us->at(0, xi, yi), 2);
        ld B2 = Bx*Bx + By*By + powl(us->at(7, xi, yi), 2);
        us->at(4, xi, yi) = P / (gamma - 1) + us->at(0, xi, yi)*v2/2 + B2/2;

    }



    return true;
}




bool MultiRTInstability(Arrays *us, ld x0, ld y0, ld xn, ld yn, ld gy) {
    int Nx = us->Nx; int Ny = us->Ny;
    ld dx = us->dx; ld dy = us->dy;

    for (int yi = 0; yi < Ny; yi++)
        for (int xi = 0; xi < Nx; xi++)
        {
            ld x = x0 + (xi-4.0)*dx;// if (x < x0) x = x0; if (x > xn) x = xn;
            ld y = y0 + (yi-4.0)*dy;// if (y < y0) y = y0; if (y > yn) y = yn;

            ld kx = 2*M_PI/(xn-x0);
            ld ky = 2*M_PI/(yn-y0);
            ld A = 0.01;
            
            if (y > 0.0) { //0.025*cosl(2*M_PI*x)
                ld rho = 2.0; //1.08
                us->at(0, xi, yi) = rho; us->at(1, xi, yi) = 0.0*rho; us->at(2, xi, yi) = (rand() % 10000)/1e6*(1+cosl(ky*y/3))/2.0*rho; us->at(3, xi, yi) = 0.0*rho;
                us->at(5, xi, yi) = 0.0; us->at(6, xi, yi) = 0.0; us->at(7, xi, yi) = 0.0;


                ld P = 2.5 + gy*rho*y; // 1.0/gamma

                ld v = sqrtl(powl(us->at(1, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(2, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(3, xi, yi)/us->at(0, xi, yi), 2));
                ld B = sqrtl(powl(us->at(5, xi, yi), 2) + powl(us->at(6, xi, yi), 2) + powl(us->at(7, xi, yi), 2));

                us->at(4, xi, yi) = P;
            }
            else {
                ld rho = 1.0; // 1.0
                us->at(0, xi, yi) = rho; us->at(1, xi, yi) = 0.0*rho; us->at(2, xi, yi) = (rand() % 10000)/1e6*(1+cosl(ky*y/3))/2.0*rho; us->at(3, xi, yi) = 0.0*rho; //vz 1.0
                us->at(5, xi, yi) = 0.0; us->at(6, xi, yi) = 0.0; us->at(7, xi, yi) = 0.0;


                ld P = 2.5 + gy*rho*y;
                
                ld v = sqrtl(powl(us->at(1, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(2, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(3, xi, yi)/us->at(0, xi, yi), 2));
                ld B = sqrtl(powl(us->at(5, xi, yi), 2) + powl(us->at(6, xi, yi), 2) + powl(us->at(7, xi, yi), 2));

                us->at(4, xi, yi) = P;
            }

        }

    setBoundaries(bound::PERIODIC, bound::REFLECTING, us, Nx-8, Ny-8, gy);

    for (int yi = 0; yi < Ny; yi++)
    for (int xi = 0; xi < Nx; xi++)
    {
        ld P = us->at(4, xi, yi);
        ld Bx = (us->at(5, xi+1, yi)+us->at(5, xi, yi))/2.0;
        ld By = (us->at(6, xi, yi+1)+us->at(6, xi, yi))/2.0;

        ld v2 = powl(us->at(1, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(2, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(3, xi, yi)/us->at(0, xi, yi), 2);
        ld B2 = Bx*Bx + By*By + powl(us->at(7, xi, yi), 2);
        us->at(4, xi, yi) = P / (gamma - 1) + us->at(0, xi, yi)*v2/2 + B2/2;

    }



    return true;
}




bool KHInstability(Arrays *us, ld x0, ld y0, ld xn, ld yn, ld gy) {
    int Nx = us->Nx; int Ny = us->Ny;
    ld dx = us->dx; ld dy = us->dy;

    for (int yi = 0; yi < Ny; yi++)
        for (int xi = 0; xi < Nx; xi++)
        {
            ld x = x0 + (xi-4.0)*dx;// if (x < x0) x = x0; if (x > xn) x = xn;
            ld y = y0 + (yi-4.0)*dy;// if (y < y0) y = y0; if (y > yn) y = yn;

            ld B0 = 0.0;
            B0 = 0.5*M_rt4PI;

            if (fabsl(y) > 0.25) { //
                ld rho = 1.0; //1.08
                us->at(0, xi, yi) = rho; us->at(1, xi, yi) = (-0.5+(rand() % 10000)/1e6)*rho; us->at(2, xi, yi) = (rand() % 10000)/1e6*rho; us->at(3, xi, yi) = 0.0*rho;
                us->at(5, xi, yi) = B0; us->at(6, xi, yi) = 0.0; us->at(7, xi, yi) = 0.0;


                ld P = 2.5;

                ld v = sqrtl(powl(us->at(1, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(2, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(3, xi, yi)/us->at(0, xi, yi), 2));
                ld B = sqrtl(powl(us->at(5, xi, yi), 2) + powl(us->at(6, xi, yi), 2) + powl(us->at(7, xi, yi), 2));

                us->at(4, xi, yi) = P;
            }
            else {
                ld rho = 2.0; // 1.0
                us->at(0, xi, yi) = rho; us->at(1, xi, yi) = (0.5+(rand() % 10000)/1e6)*rho; us->at(2, xi, yi) = (rand() % 10000)/1e6*rho; us->at(3, xi, yi) = 0.0*rho; //vz 1.0
                us->at(5, xi, yi) = B0; us->at(6, xi, yi) = 0.0; us->at(7, xi, yi) = 0.0;


                ld P = 2.5;

                ld v = sqrtl(powl(us->at(1, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(2, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(3, xi, yi)/us->at(0, xi, yi), 2));
                ld B = sqrtl(powl(us->at(5, xi, yi), 2) + powl(us->at(6, xi, yi), 2) + powl(us->at(7, xi, yi), 2));

                us->at(4, xi, yi) = P;
            }

        }

    setBoundaries(bound::PERIODIC, bound::PERIODIC, us, Nx-8, Ny-8, 0);

    for (int yi = 0; yi < Ny; yi++)
        for (int xi = 0; xi < Nx; xi++)
        {
            ld P = us->at(4, xi, yi);
            ld Bx = (us->at(5, xi+1, yi)+us->at(5, xi, yi))/2.0;
            ld By = (us->at(6, xi, yi+1)+us->at(6, xi, yi))/2.0;

            ld v2 = powl(us->at(1, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(2, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(3, xi, yi)/us->at(0, xi, yi), 2);
            ld B2 = Bx*Bx + By*By + powl(us->at(7, xi, yi), 2);
            us->at(4, xi, yi) = P / (gamma - 1) + us->at(0, xi, yi)*v2/2 + B2/2;

        }



    return true;
}