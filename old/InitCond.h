#pragma once
#include <math.h>
#include <iostream>
#include "Arrays.h"
#include "MHD.h"
#include "Eigenvalues.h"


bool Vortex(Arrays *us, long double x0, long double y0, long double xn, long double yn) {
    int Nx = us->Nx; int Ny = us->Ny;
    long double dx = us->dx; long double dy = us->dy;

    for (int yi = 0; yi < Ny; yi++)
        for (int xi = 0; xi < Nx; xi++)
        {
            long double x = x0 + (xi-4.0)*dx; if (x < x0) x = x0; if (x > xn) x = xn;
            long double y = y0 + (yi-4.0)*dy; if (y < y0) y = y0; if (y > yn) y = yn;
            long double x_ix = x0 + (xi-4.5)*dx; //x at the x interface (where Bx is)
            long double y_iy = y0 + (yi-4.5)*dy; //y at the y interface (where By is)


            long double B0 = 1/M_rt4PI;
            //B0 = 0.0;

            long double a[2][2] = { 0 };

            for (int i = 0; i < 2; i++) //Computes the potential at corner (i-1/2, j-1/2)
            for (int j = 0; j < 2; j++)
            {
                int xii = i + xi; int yii = j + yi;
                long double x_ix = x0 + (xii-4.5)*dx; //Coordinates at the x and y interfaces
                long double y_iy = y0 + (yii-4.5)*dy;


                a[i][j] = B0*(cosl(4*M_PI*x_ix)/(4*M_PI)+cosl(2*M_PI*y_iy)/(2*M_PI));
            }

            us->at(5, xi, yi) = (a[0][1]-a[0][0])/dy; us->at(6, xi, yi) = -(a[1][0]-a[0][0])/dx;





            long double rho = 25.0/(36.0*M_PI);
            us->at(0, xi, yi) = rho; us->at(1, xi, yi) = -rho*sinl(2.0*M_PI*y); us->at(2, xi, yi) = rho*sinl(2.0*M_PI*x); us->at(3, xi, yi) = rho*0.0;
            //us->at(5, xi, yi) = -B0*sinl(2*M_PI*y); us->at(6, xi, yi) = B0*sinl(4*M_PI*x); us->at(7, xi, yi) = 0.0;
            

            long double P = 5.0/(12.0*M_PI);
            long double v = sqrtl(powl(us->at(1, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(2, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(3, xi, yi)/us->at(0, xi, yi), 2));
            long double B = sqrtl(powl(us->at(5, xi, yi), 2) + powl(us->at(6, xi, yi), 2) + powl(us->at(7, xi, yi), 2));

            us->at(4, xi, yi) = P;// / (gamma - 1) + us->at(0, xi, yi)*v*v/2 + B*B/2;

        }

    //Do the transformation P->E afterwards to avoid using face-centered Bx and By values as cell-centered ones.
    setBoundaries(bound::PERIODIC, bound::PERIODIC, us, Nx-8, Ny-8, 0);

    for (int yi = 0; yi < Ny; yi++)
    for (int xi = 0; xi < Nx; xi++)
    {
        long double P = us->at(4, xi, yi);
        long double Bx = (us->at(5, xi+1, yi)+us->at(5, xi, yi))/2.0;
        long double By = (us->at(6, xi, yi+1)+us->at(6, xi, yi))/2.0;

        long double v2 = powl(us->at(1, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(2, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(3, xi, yi)/us->at(0, xi, yi), 2);
        long double B2 = Bx*Bx + By*By + powl(us->at(7, xi, yi), 2);
        us->at(4, xi, yi) = P / (gamma - 1) + us->at(0, xi, yi)*v2/2 + B2/2;
            

        for (int i = 0; i < 8; i++)
        {
            if (isnan(us->at(i, xi, yi))) { cout << i << endl; return false; }
        }
    }

    return true;
}


bool RandomTest(Arrays *us, long double x0, long double y0, long double xn, long double yn) {
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

bool Loops(Arrays* us, long double x0, long double y0, long double xn, long double yn) {
    int Nx = us->Nx; int Ny = us->Ny;
    long double dx = us->dx; long double dy = us->dy;

    for (int yi = 0; yi < Ny; yi++)
        for (int xi = 0; xi < Nx; xi++)
        {
            long double x = x0 + (xi-4.0)*dx;// if (x < x0) x = x0; if (x > xn) x = xn;
            long double y = y0 + (yi-4.0)*dy;// if (y < y0) y = y0; if (y > yn) y = yn;
            long double r = sqrtl(x*x+y*y);

            long double R0 = 0.3; long double A0 = 0.001;
            long double a[2][2] = { 0 };

            for (int i = 0; i < 2; i++) //Computes the potential at corner (i-1/2, j-1/2)
            for (int j = 0; j < 2; j++)
            {
                int xii = i + xi; int yii = j + yi;
                long double x_ix = x0 + (xii-4.5)*dx; //Coordinates at the x and y interfaces
                long double y_iy = y0 + (yii-4.5)*dy;

                long double rc = sqrtl(x_ix*x_ix+y_iy*y_iy);

                if (rc < R0)
                    a[i][j] = A0*(R0-rc);
                else
                    a[i][j] = 0.0;
            }

            us->at(5, xi, yi) = (a[0][1]-a[0][0])/dy; us->at(6, xi, yi) = -(a[1][0]-a[0][0])/dx;




            if (r < R0) {
                long double rho = 1.0;
                us->at(0, xi, yi) = rho; us->at(1, xi, yi) = 2.0*rho; us->at(2, xi, yi) = 1.0*rho; us->at(3, xi, yi) = 0.0*rho;
                

                long double P = 1.0;
                long double v = sqrtl(powl(us->at(1, xi, yi)/us->at(0, xi, yi),2) + powl(us->at(2, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(3, xi, yi)/us->at(0, xi, yi), 2));
                long double B = sqrtl(powl(us->at(5, xi, yi),2) + powl(us->at(6, xi, yi), 2) + powl(us->at(7, xi, yi), 2));

                us->at(4, xi, yi) = P;// / (gamma - 1) + us->at(0, xi, yi)*v*v/2 + B*B/2;
            }
            else {
                long double rho = 1.0;
                us->at(0, xi, yi) = rho; us->at(1, xi, yi) = 2.0*rho; us->at(2, xi, yi) = 1.0*rho; us->at(3, xi, yi) = 0.0;
                
                long double P = 1.0;
                long double v = sqrtl(powl(us->at(1, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(2, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(3, xi, yi)/us->at(0, xi, yi), 2));
                long double B = sqrtl(powl(us->at(5, xi, yi), 2) + powl(us->at(6, xi, yi), 2) + powl(us->at(7, xi, yi), 2));

                us->at(4, xi, yi) = P;// / (gamma - 1) + us->at(0, xi, yi)*v*v/2 + B*B/2;
            }

        }

    setBoundaries(bound::PERIODIC, bound::PERIODIC, us, Nx-8, Ny-8, 0);

    for (int yi = 0; yi < Ny; yi++)
        for (int xi = 0; xi < Nx; xi++)
        {
            long double P = us->at(4, xi, yi);
            long double Bx = (us->at(5, xi+1, yi)+us->at(5, xi, yi))/2.0;
            long double By = (us->at(6, xi, yi+1)+us->at(6, xi, yi))/2.0;

            long double v2 = powl(us->at(1, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(2, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(3, xi, yi)/us->at(0, xi, yi), 2);
            long double B2 = Bx*Bx + By*By + powl(us->at(7, xi, yi), 2);
            us->at(4, xi, yi) = P / (gamma - 1) + us->at(0, xi, yi)*v2/2 + B2/2;
            

            for (int i = 0; i < 8; i++)
            {
                if (isnan(us->at(i, xi, yi))) { cout << i << endl; return false; }
            }
        }
    return true;
}



bool Sod(Arrays *us, long double x0, long double y0, long double xn, long double yn, bool hd) {
    int Nx = us->Nx; int Ny = us->Ny;
    long double dx = us->dx; long double dy = us->dy;

    for (int yi = 0; yi < Ny; yi++)
        for (int xi = 0; xi < Nx; xi++)
        {
            long double x = x0 + (xi-4.0)*dx;// if (x < x0) x = x0; if (x > xn) x = xn;
            long double y = y0 + (yi-4.0)*dy;// if (y < y0) y = y0; if (y > yn) y = yn;


            if (x < 0.5) {
                long double rho = 1.0;
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
                

                long double P = 1.0;
                long double v = sqrtl(powl(us->at(1, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(2, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(3, xi, yi)/us->at(0, xi, yi), 2));
                long double B = sqrtl(powl(us->at(5, xi, yi), 2) + powl(us->at(6, xi, yi), 2) + powl(us->at(7, xi, yi), 2));

                us->at(4, xi, yi) = P;
            }
            else {
                long double rho = 0.125;
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

                long double P = 0.1;
                long double v = sqrtl(powl(us->at(1, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(2, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(3, xi, yi)/us->at(0, xi, yi), 2));
                long double B = sqrtl(powl(us->at(5, xi, yi), 2) + powl(us->at(6, xi, yi), 2) + powl(us->at(7, xi, yi), 2));

                us->at(4, xi, yi) = P;
            }

        }

    setBoundaries(bound::OUTFLOW, bound::OUTFLOW, us, Nx-8, Ny-8, 0);

    for (int yi = 0; yi < Ny; yi++)
    for (int xi = 0; xi < Nx; xi++)
    {
        long double P = us->at(4, xi, yi);
        long double Bx = (us->at(5, xi+1, yi)+us->at(5, xi, yi))/2.0;
        long double By = (us->at(6, xi, yi+1)+us->at(6, xi, yi))/2.0;

        long double v2 = powl(us->at(1, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(2, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(3, xi, yi)/us->at(0, xi, yi), 2);
        long double B2 = Bx*Bx + By*By + powl(us->at(7, xi, yi), 2);
        us->at(4, xi, yi) = P / (gamma - 1) + us->at(0, xi, yi)*v2/2 + B2/2;

    }



    return true;
}



long double VP(long double R, long double r0, long double q) {
    return sqrtl(r0)/(r0-2.0)*powl(R/r0, 1.0-q);
}


bool Disk(Arrays *us, long double x0, long double y0, long double rs) {
    int Nx = us->Nx; int Ny = us->Ny;
    long double dx = us->dx; long double dy = us->dy;

    //Disk parameters:
    double q = 1.68; double r0 = 20.0; double rhomax = 300.0; double rin = 8.0; double rho0 = 0.01; double e0 = 0.0001;
    double dcut = 1.0; double beta = 100.0; double seed = 0.01;
    long double n = 1/(gamma-1); long double q1 = 2*(q-1); long double r02 = r0-rs;
    long double f = powl(r0, 2*q-1)/(q1*r02*r02); long double C = f*powl(rin, -q1)-1.0/(rin-rs); long double K = (C+1/r02-f*powl(r0, -q1))/powl(rhomax, 1.0/n);

    for (int yi = 0; yi < Ny; yi++)
        for (int xi = 0; xi < Nx; xi++)
        {
            long double R = x0 + (xi-4.0)*dx;
            long double z = y0 + (yi-4.0)*dy;
            long double R_iR = x0 + (xi-4.5)*dx; //x at the x interface (where Bx is)
            long double z_iz = y0 + (yi-4.5)*dy; //y at the y interface (where By is)
            long double r = sqrtl(R*R+z*z);

            us->at(0, xi, yi) = rho0; us->at(3, xi, yi) = VP(R, r0, q)*rho0; us->at(4, xi, yi) = e0;

            long double Ta = f*powl(R, -q1); long double T = 1.0/Ta + C;
            long double cs = 0.086;
            if (R >= rin && r <= T) {
                long double h = fmaxl(C + 1.0/(r-2.0) - f*powl(R, -2.0*(q-1.0)), 0);
                long double rho = powl(h/K, n)*(R >= rin); rho = fmaxl(rho, rho0);
                long double Eint = powl(rho, gamma)*K/((gamma-1)*(n+1)); Eint *= (1-0.01*sinl(R));
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

    //        long double a[2][2] = { 0 };

    //        for (int i = 0; i < 2; i++) //Computes the potential at corner (i-1/2, j-1/2)
    //            for (int j = 0; j < 2; j++)
    //            {
    //                int xii = i + xi; int yii = j + yi;
    //                long double x_ix = x0 + (xii-4.5)*dx; //Coordinates at the x and y interfaces
    //                long double y_iy = y0 + (yii-4.5)*dy;


    //                a[i][j] = us->at(4, xi, yi)-0.1;
    //            }

    //        us->at(5, xi, yi) = (a[0][1]-a[0][0])/dy; us->at(6, xi, yi) = -(a[1][0]-a[0][0])/dx;



    //    }

    //Do the transformation P->E afterwards to avoid using face-centered Bx and By values as cell-centered ones.
    //setBoundaries(bound::OUTFLOW, bound::PERIODIC, us, Nx-8, Ny-8, 0);

    for (int yi = 0; yi < Ny; yi++)
        for (int xi = 0; xi < Nx; xi++)
        {
            long double P = us->at(4, xi, yi);
            long double Bx = (us->at(5, xi+1, yi)+us->at(5, xi, yi))/2.0;
            long double By = (us->at(6, xi, yi+1)+us->at(6, xi, yi))/2.0;

            long double v2 = powl(us->at(1, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(2, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(3, xi, yi)/us->at(0, xi, yi), 2);
            long double B2 = Bx*Bx + By*By + powl(us->at(7, xi, yi), 2);
            us->at(4, xi, yi) = P / (gamma - 1) + us->at(0, xi, yi)*v2/2 + B2/2;


            for (int i = 0; i < 8; i++)
            {
                if (isnan(us->at(i, xi, yi))) { cout << i << endl; return false; }
            }
        }

    return true;
}




long double rotor_f(long double r, long double r0, long double r1) {
    return (r1-r)/(r1-r0);
}

bool MHDRotor(Arrays* us, long double x0, long double y0) {
    int Nx = us->Nx; int Ny = us->Ny;
    long double dx = us->dx; long double dy = us->dy;
    cout << "Initializing the MHD Rotor problem with resolution " << Nx << "x" << Ny << endl;

    for (int yi = 0; yi < Ny; yi++)
        for (int xi = 0; xi < Nx; xi++)
        {
            //cout << xi << ", " << yi << endl;
            long double x = x0 + (xi-4.0)*dx;// if (x < x0) x = x0; if (x > xn) x = xn;
            long double y = y0 + (yi-4.0)*dy;// if (y < y0) y = y0; if (y > yn) y = yn;
            long double r = sqrtl((x-0.5)*(x-0.5) + (y-0.5)*(y-0.5));
            

            long double r0 = 0.1; long double r1 = 0.115; long double v0 = 1.0; long double P = 0.5;
            if (r <= r0) {
                long double rho = 10.0;
                us->at(0, xi, yi) = rho; us->at(1, xi, yi) = -rotor_f(r, r0, r1)*v0*(y-0.5)/r0 * rho; us->at(2, xi, yi) = rotor_f(r, r0, r1)*v0*(x-0.5)/r0 * rho; us->at(3, xi, yi) = 0.0*rho;
                us->at(5, xi, yi) = 2.5/sqrtl(4*M_PI); us->at(6, xi, yi) = 0.0; us->at(7, xi, yi) = 0.0;
                
                
                long double v = sqrtl(powl(us->at(1, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(2, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(3, xi, yi)/us->at(0, xi, yi), 2));
                long double B = sqrtl(powl(us->at(5, xi, yi), 2) + powl(us->at(6, xi, yi), 2) + powl(us->at(7, xi, yi), 2));

                us->at(4, xi, yi) = P / (gamma - 1) + us->at(0, xi, yi)*v*v/2 + B*B/2;
            }
            if (r > r0 && r < r1){
                long double rho = 1.0 + 9.0*rotor_f(r, r0, r1);
                us->at(0, xi, yi) = rho; us->at(1, xi, yi) = -rotor_f(r, r0, r1)*v0*(y-0.5)/r * rho; us->at(2, xi, yi) = rotor_f(r, r0, r1)*v0*(x-0.5)/r * rho; us->at(3, xi, yi) = 0.0*rho;
                us->at(5, xi, yi) = 2.5/sqrtl(4*M_PI); us->at(6, xi, yi) = 0.0; us->at(7, xi, yi) = 0.0;

                
                long double v = sqrtl(powl(us->at(1, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(2, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(3, xi, yi)/us->at(0, xi, yi), 2));
                long double B = sqrtl(powl(us->at(5, xi, yi), 2) + powl(us->at(6, xi, yi), 2) + powl(us->at(7, xi, yi), 2));

                us->at(4, xi, yi) = P / (gamma - 1) + us->at(0, xi, yi)*v*v/2 + B*B/2;
            }
            if (r >= r1) {
                long double rho = 1.0;
                us->at(0, xi, yi) = rho; us->at(1, xi, yi) = 0.0*rho; us->at(2, xi, yi) = 0.0*rho; us->at(3, xi, yi) = 0.0;
                us->at(5, xi, yi) = 2.5/sqrtl(4*M_PI); us->at(6, xi, yi) = 0.0; us->at(7, xi, yi) = 0.0; //5.0/M_rt4PI

                
                long double v = sqrtl(powl(us->at(1, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(2, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(3, xi, yi)/us->at(0, xi, yi), 2));
                long double B = sqrtl(powl(us->at(5, xi, yi), 2) + powl(us->at(6, xi, yi), 2) + powl(us->at(7, xi, yi), 2));

                us->at(4, xi, yi) = P / (gamma - 1) + us->at(0, xi, yi)*v*v/2 + B*B/2;
            }

            long double u[8] = { 0 }; for (int i = 0; i < 8; i++) u[i] = us->at(i, xi, yi);
            long double w[8] = { 0 }; if (!toPrimitive(u, w)) { cout << "Init failed" << endl; return false; }
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



bool Sedov(Arrays *us, long double x0, long double y0) {
    int Nx = us->Nx; int Ny = us->Ny;
    long double dx = us->dx; long double dy = us->dy;

    for (int yi = 0; yi < Ny; yi++)
        for (int xi = 0; xi < Nx; xi++)
        {
            long double x = x0 + (xi-4.0)*dx;// if (x < x0) x = x0; if (x > xn) x = xn;
            long double y = y0 + (yi-4.0)*dy;// if (y < y0) y = y0; if (y > yn) y = yn;
            long double r = sqrtl(x*x+y*y);

            long double dr = 3.5*dx; long double rho = 1.0; long double nu = 1.0;
            if (r < dr) {
                us->at(0, xi, yi) = rho; us->at(1, xi, yi) = 0.0*rho; us->at(2, xi, yi) = 0.0*rho; us->at(3, xi, yi) = 0.0*rho;
                us->at(5, xi, yi) = 0.0; us->at(6, xi, yi) = 0.0; us->at(7, xi, yi) = 0.0;


                long double P = 3*(gamma-1)*1.0/((nu+1)*M_PI*powl(dr, nu));
                long double v = sqrtl(powl(us->at(1, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(2, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(3, xi, yi)/us->at(0, xi, yi), 2));
                long double B = sqrtl(powl(us->at(5, xi, yi), 2) + powl(us->at(6, xi, yi), 2) + powl(us->at(7, xi, yi), 2));

                us->at(4, xi, yi) = P / (gamma - 1) + us->at(0, xi, yi)*v*v/2 + B*B/2;
                //us->at(4, xi, yi) = 1.0;
            }
            else {
                us->at(0, xi, yi) = rho; us->at(1, xi, yi) = 0.0*rho; us->at(2, xi, yi) = 0.0*rho; us->at(3, xi, yi) = 0.0;
                us->at(5, xi, yi) = 0.0; us->at(6, xi, yi) = 0.0; us->at(7, xi, yi) = 0.0;

                long double P = 1e-5;
                long double v = sqrtl(powl(us->at(1, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(2, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(3, xi, yi)/us->at(0, xi, yi), 2));
                long double B = sqrtl(powl(us->at(5, xi, yi), 2) + powl(us->at(6, xi, yi), 2) + powl(us->at(7, xi, yi), 2));

                us->at(4, xi, yi) = P / (gamma - 1) + us->at(0, xi, yi)*v*v/2 + B*B/2;
            }

            long double u[8] = { 0 }; for (int i = 0; i < 8; i++) u[i] = us->at(i, xi, yi);
            long double w[8] = { 0 }; if (!toPrimitive(u, w)) { cout << "Init failed" << endl; return false; }

        }

    cout << "Initialized the variables for Sedov test." << endl;
    return true;
}





bool RTInstability(Arrays *us, long double x0, long double y0, long double xn, long double yn, long double gy) {
    int Nx = us->Nx; int Ny = us->Ny;
    long double dx = us->dx; long double dy = us->dy;

    for (int yi = 0; yi < Ny; yi++)
        for (int xi = 0; xi < Nx; xi++)
        {
            long double x = x0 + (xi-4.0)*dx;// if (x < x0) x = x0; if (x > xn) x = xn;
            long double y = y0 + (yi-4.0)*dy;// if (y < y0) y = y0; if (y > yn) y = yn;


            long double kx = 2*M_PI/(xn-x0);// if (x0 == 0.0 || xn == 0.0) kx /= 2.0;
            long double ky = 2*M_PI/(yn-y0);
            long double A = 0.01;
            
            if (y > 0.0) { //0.025*cosl(2*M_PI*x)
                long double rho = 2.0; //1.08
                us->at(0, xi, yi) = rho; us->at(1, xi, yi) = 0.0*rho; us->at(2, xi, yi) = 0.01*(1+cosl(kx*x))*(1+cosl(ky*y))/4.0*rho; us->at(3, xi, yi) = 0.0*rho;
                //us->at(0, xi, yi) = rho; us->at(1, xi, yi) = 0.0*rho; us->at(2, xi, yi) = 0.0*rho; us->at(3, xi, yi) = 0.0*rho;
                us->at(5, xi, yi) = 0.0; us->at(6, xi, yi) = 0.0; us->at(7, xi, yi) = 0.0;

                //if (fabsl(y) < 2*dx) us->at(2, xi, yi) = A*sinl(kx*x)*rho;


                long double P = 1.0/gamma +gy*rho*y;

                long double v = sqrtl(powl(us->at(1, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(2, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(3, xi, yi)/us->at(0, xi, yi), 2));
                long double B = sqrtl(powl(us->at(5, xi, yi), 2) + powl(us->at(6, xi, yi), 2) + powl(us->at(7, xi, yi), 2));

                us->at(4, xi, yi) = P;
            }
            else {
                long double rho = 1.0; // 1.0
                us->at(0, xi, yi) = rho; us->at(1, xi, yi) = 0.0*rho; us->at(2, xi, yi) = 0.01*(1+cosl(kx*x))*(1+cosl(ky*y))/4.0*rho; us->at(3, xi, yi) = 0.0*rho; //vz 1.0
                //us->at(0, xi, yi) = rho; us->at(1, xi, yi) = 0.0*rho; us->at(2, xi, yi) = 0.0*rho; us->at(3, xi, yi) = 0.0*rho;
                us->at(5, xi, yi) = 0.0; us->at(6, xi, yi) = 0.0; us->at(7, xi, yi) = 0.0;

                //if (fabsl(y) < 2*dx) us->at(2, xi, yi) = A*sinl(kx*x)*rho;

                long double P = 1.0/gamma +gy*rho*y;
                
                long double v = sqrtl(powl(us->at(1, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(2, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(3, xi, yi)/us->at(0, xi, yi), 2));
                long double B = sqrtl(powl(us->at(5, xi, yi), 2) + powl(us->at(6, xi, yi), 2) + powl(us->at(7, xi, yi), 2));

                us->at(4, xi, yi) = P;
            }

        }

    setBoundaries(bound::PERIODIC, bound::PRESSURE, us, Nx-8, Ny-8, gy);

    for (int yi = 0; yi < Ny; yi++)
    for (int xi = 0; xi < Nx; xi++)
    {
        long double P = us->at(4, xi, yi);
        long double Bx = (us->at(5, xi+1, yi)+us->at(5, xi, yi))/2.0;
        long double By = (us->at(6, xi, yi+1)+us->at(6, xi, yi))/2.0;

        long double v2 = powl(us->at(1, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(2, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(3, xi, yi)/us->at(0, xi, yi), 2);
        long double B2 = Bx*Bx + By*By + powl(us->at(7, xi, yi), 2);
        us->at(4, xi, yi) = P / (gamma - 1) + us->at(0, xi, yi)*v2/2 + B2/2;

    }



    return true;
}




bool MHDRTInstability(Arrays *us, long double x0, long double y0, long double xn, long double yn, long double gy) {
    int Nx = us->Nx; int Ny = us->Ny;
    long double dx = us->dx; long double dy = us->dy;

    for (int yi = 0; yi < Ny; yi++)
        for (int xi = 0; xi < Nx; xi++)
        {
            long double x = x0 + (xi-4.0)*dx;// if (x < x0) x = x0; if (x > xn) x = xn;
            long double y = y0 + (yi-4.0)*dy;// if (y < y0) y = y0; if (y > yn) y = yn;


            long double kx = 2*M_PI/(xn-x0);// if (x0 == 0.0 || xn == 0.0) kx /= 2.0;
            long double ky = 2*M_PI/(yn-y0);
            long double A = 0.01;
            long double B0 = 0.0;
            B0 = 0.05*0.14;
            
            if (y > 0.0) { //0.025*cosl(2*M_PI*x)
                long double rho = 3.0; //1.08
                us->at(0, xi, yi) = rho; us->at(1, xi, yi) = 0.0*rho; us->at(2, xi, yi) = 0.01*(1+cosl(kx*x))*(1+cosl(ky*y))/4.0*rho; us->at(3, xi, yi) = 0.0*rho;
                //us->at(0, xi, yi) = rho; us->at(1, xi, yi) = 0.0*rho; us->at(2, xi, yi) = 0.0*rho; us->at(3, xi, yi) = 0.0*rho;
                us->at(5, xi, yi) = B0; us->at(6, xi, yi) = 0.0; us->at(7, xi, yi) = 0.0;

                //if (fabsl(y) < 2*dx) us->at(2, xi, yi) = A*sinl(kx*x)*rho;


                long double P = 1.0/gamma +gy*rho*y;

                long double v = sqrtl(powl(us->at(1, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(2, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(3, xi, yi)/us->at(0, xi, yi), 2));
                long double B = sqrtl(powl(us->at(5, xi, yi), 2) + powl(us->at(6, xi, yi), 2) + powl(us->at(7, xi, yi), 2));

                us->at(4, xi, yi) = P;
            }
            else {
                long double rho = 1.0; // 1.0
                us->at(0, xi, yi) = rho; us->at(1, xi, yi) = 0.0*rho; us->at(2, xi, yi) = 0.01*(1+cosl(kx*x))*(1+cosl(ky*y))/4.0*rho; us->at(3, xi, yi) = 0.0*rho; //vz 1.0
                //us->at(0, xi, yi) = rho; us->at(1, xi, yi) = 0.0*rho; us->at(2, xi, yi) = 0.0*rho; us->at(3, xi, yi) = 0.0*rho;
                us->at(5, xi, yi) = B0; us->at(6, xi, yi) = 0.0; us->at(7, xi, yi) = 0.0;

                //if (fabsl(y) < 2*dx) us->at(2, xi, yi) = A*sinl(kx*x)*rho;

                long double P = 1.0/gamma +gy*rho*y;
                
                long double v = sqrtl(powl(us->at(1, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(2, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(3, xi, yi)/us->at(0, xi, yi), 2));
                long double B = sqrtl(powl(us->at(5, xi, yi), 2) + powl(us->at(6, xi, yi), 2) + powl(us->at(7, xi, yi), 2));

                us->at(4, xi, yi) = P;
            }

        }

    setBoundaries(bound::PERIODIC, bound::PRESSURE, us, Nx-8, Ny-8, gy);

    for (int yi = 0; yi < Ny; yi++)
    for (int xi = 0; xi < Nx; xi++)
    {
        long double P = us->at(4, xi, yi);
        long double Bx = (us->at(5, xi+1, yi)+us->at(5, xi, yi))/2.0;
        long double By = (us->at(6, xi, yi+1)+us->at(6, xi, yi))/2.0;

        long double v2 = powl(us->at(1, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(2, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(3, xi, yi)/us->at(0, xi, yi), 2);
        long double B2 = Bx*Bx + By*By + powl(us->at(7, xi, yi), 2);
        us->at(4, xi, yi) = P / (gamma - 1) + us->at(0, xi, yi)*v2/2 + B2/2;

    }



    return true;
}






bool HDEq(Arrays *us, long double x0, long double y0, long double xn, long double yn, long double gy) {
    int Nx = us->Nx; int Ny = us->Ny;
    long double dx = us->dx; long double dy = us->dy;

    for (int yi = 0; yi < Ny; yi++)
        for (int xi = 0; xi < Nx; xi++)
        {
            long double x = x0 + (xi-4.0)*dx;// if (x < x0) x = x0; if (x > xn) x = xn;
            long double y = y0 + (yi-4.0)*dy;// if (y < y0) y = y0; if (y > yn) y = yn;

            
            long double rho = 2.0; //1.08
            //us->at(0, xi, yi) = rho; us->at(1, xi, yi) = 0.0*rho; us->at(2, xi, yi) = 0.01*(1+cosl(kx*x))*(1+cosl(ky*y))/4.0*rho; us->at(3, xi, yi) = 0.0*rho;
            us->at(0, xi, yi) = rho; us->at(1, xi, yi) = 0.0*rho; us->at(2, xi, yi) = 0.0*rho; us->at(3, xi, yi) = 0.0*rho;
            us->at(5, xi, yi) = 0.0; us->at(6, xi, yi) = 0.0; us->at(7, xi, yi) = 0.0;

            //if (fabsl(y) < 2*dx) us->at(2, xi, yi) = A*sinl(kx*x)*rho;


            long double P = 1.0/gamma;// +gy*rho*y;

            long double v = sqrtl(powl(us->at(1, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(2, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(3, xi, yi)/us->at(0, xi, yi), 2));
            long double B = sqrtl(powl(us->at(5, xi, yi), 2) + powl(us->at(6, xi, yi), 2) + powl(us->at(7, xi, yi), 2));

            us->at(4, xi, yi) = P;
            

        }

    setBoundaries(bound::PERIODIC, bound::REFLECTING, us, Nx-8, Ny-8, gy);

    for (int yi = 0; yi < Ny; yi++)
    for (int xi = 0; xi < Nx; xi++)
    {
        long double P = us->at(4, xi, yi);
        long double Bx = (us->at(5, xi+1, yi)+us->at(5, xi, yi))/2.0;
        long double By = (us->at(6, xi, yi+1)+us->at(6, xi, yi))/2.0;

        long double v2 = powl(us->at(1, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(2, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(3, xi, yi)/us->at(0, xi, yi), 2);
        long double B2 = Bx*Bx + By*By + powl(us->at(7, xi, yi), 2);
        us->at(4, xi, yi) = P / (gamma - 1) + us->at(0, xi, yi)*v2/2 + B2/2;

    }



    return true;
}




bool MultiRTInstability(Arrays *us, long double x0, long double y0, long double xn, long double yn, long double gy) {
    int Nx = us->Nx; int Ny = us->Ny;
    long double dx = us->dx; long double dy = us->dy;

    for (int yi = 0; yi < Ny; yi++)
        for (int xi = 0; xi < Nx; xi++)
        {
            long double x = x0 + (xi-4.0)*dx;// if (x < x0) x = x0; if (x > xn) x = xn;
            long double y = y0 + (yi-4.0)*dy;// if (y < y0) y = y0; if (y > yn) y = yn;

            long double kx = 2*M_PI/(xn-x0);
            long double ky = 2*M_PI/(yn-y0);
            long double A = 0.01;
            
            if (y > 0.0) { //0.025*cosl(2*M_PI*x)
                long double rho = 2.0; //1.08
                us->at(0, xi, yi) = rho; us->at(1, xi, yi) = 0.0*rho; us->at(2, xi, yi) = (rand() % 10000)/1e6*(1+cosl(ky*y/3))/2.0*rho; us->at(3, xi, yi) = 0.0*rho;
                us->at(5, xi, yi) = 0.0; us->at(6, xi, yi) = 0.0; us->at(7, xi, yi) = 0.0;


                long double P = 2.5 + gy*rho*y; // 1.0/gamma

                long double v = sqrtl(powl(us->at(1, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(2, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(3, xi, yi)/us->at(0, xi, yi), 2));
                long double B = sqrtl(powl(us->at(5, xi, yi), 2) + powl(us->at(6, xi, yi), 2) + powl(us->at(7, xi, yi), 2));

                us->at(4, xi, yi) = P;
            }
            else {
                long double rho = 1.0; // 1.0
                us->at(0, xi, yi) = rho; us->at(1, xi, yi) = 0.0*rho; us->at(2, xi, yi) = (rand() % 10000)/1e6*(1+cosl(ky*y/3))/2.0*rho; us->at(3, xi, yi) = 0.0*rho; //vz 1.0
                us->at(5, xi, yi) = 0.0; us->at(6, xi, yi) = 0.0; us->at(7, xi, yi) = 0.0;


                long double P = 2.5 + gy*rho*y;
                
                long double v = sqrtl(powl(us->at(1, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(2, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(3, xi, yi)/us->at(0, xi, yi), 2));
                long double B = sqrtl(powl(us->at(5, xi, yi), 2) + powl(us->at(6, xi, yi), 2) + powl(us->at(7, xi, yi), 2));

                us->at(4, xi, yi) = P;
            }

        }

    setBoundaries(bound::PERIODIC, bound::REFLECTING, us, Nx-8, Ny-8, gy);

    for (int yi = 0; yi < Ny; yi++)
    for (int xi = 0; xi < Nx; xi++)
    {
        long double P = us->at(4, xi, yi);
        long double Bx = (us->at(5, xi+1, yi)+us->at(5, xi, yi))/2.0;
        long double By = (us->at(6, xi, yi+1)+us->at(6, xi, yi))/2.0;

        long double v2 = powl(us->at(1, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(2, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(3, xi, yi)/us->at(0, xi, yi), 2);
        long double B2 = Bx*Bx + By*By + powl(us->at(7, xi, yi), 2);
        us->at(4, xi, yi) = P / (gamma - 1) + us->at(0, xi, yi)*v2/2 + B2/2;

    }



    return true;
}




bool KHInstability(Arrays *us, long double x0, long double y0, long double xn, long double yn, long double gy) {
    int Nx = us->Nx; int Ny = us->Ny;
    long double dx = us->dx; long double dy = us->dy;

    for (int yi = 0; yi < Ny; yi++)
        for (int xi = 0; xi < Nx; xi++)
        {
            long double x = x0 + (xi-4.0)*dx;// if (x < x0) x = x0; if (x > xn) x = xn;
            long double y = y0 + (yi-4.0)*dy;// if (y < y0) y = y0; if (y > yn) y = yn;

            long double B0 = 0.0;
            B0 = 0.5*M_rt4PI;

            if (fabsl(y) > 0.25) { //
                long double rho = 1.0; //1.08
                us->at(0, xi, yi) = rho; us->at(1, xi, yi) = (-0.5+(rand() % 10000)/1e6)*rho; us->at(2, xi, yi) = (rand() % 10000)/1e6*rho; us->at(3, xi, yi) = 0.0*rho;
                us->at(5, xi, yi) = B0; us->at(6, xi, yi) = 0.0; us->at(7, xi, yi) = 0.0;


                long double P = 2.5;

                long double v = sqrtl(powl(us->at(1, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(2, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(3, xi, yi)/us->at(0, xi, yi), 2));
                long double B = sqrtl(powl(us->at(5, xi, yi), 2) + powl(us->at(6, xi, yi), 2) + powl(us->at(7, xi, yi), 2));

                us->at(4, xi, yi) = P;
            }
            else {
                long double rho = 2.0; // 1.0
                us->at(0, xi, yi) = rho; us->at(1, xi, yi) = (0.5+(rand() % 10000)/1e6)*rho; us->at(2, xi, yi) = (rand() % 10000)/1e6*rho; us->at(3, xi, yi) = 0.0*rho; //vz 1.0
                us->at(5, xi, yi) = B0; us->at(6, xi, yi) = 0.0; us->at(7, xi, yi) = 0.0;


                long double P = 2.5;

                long double v = sqrtl(powl(us->at(1, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(2, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(3, xi, yi)/us->at(0, xi, yi), 2));
                long double B = sqrtl(powl(us->at(5, xi, yi), 2) + powl(us->at(6, xi, yi), 2) + powl(us->at(7, xi, yi), 2));

                us->at(4, xi, yi) = P;
            }

        }

    setBoundaries(bound::PERIODIC, bound::PERIODIC, us, Nx-8, Ny-8, 0);

    for (int yi = 0; yi < Ny; yi++)
        for (int xi = 0; xi < Nx; xi++)
        {
            long double P = us->at(4, xi, yi);
            long double Bx = (us->at(5, xi+1, yi)+us->at(5, xi, yi))/2.0;
            long double By = (us->at(6, xi, yi+1)+us->at(6, xi, yi))/2.0;

            long double v2 = powl(us->at(1, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(2, xi, yi)/us->at(0, xi, yi), 2) + powl(us->at(3, xi, yi)/us->at(0, xi, yi), 2);
            long double B2 = Bx*Bx + By*By + powl(us->at(7, xi, yi), 2);
            us->at(4, xi, yi) = P / (gamma - 1) + us->at(0, xi, yi)*v2/2 + B2/2;

        }



    return true;
}