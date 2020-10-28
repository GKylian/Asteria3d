#pragma once
#include <math.h>
#include <iostream>
#include "Arrays.h"


//const double gamma = 1.4;



//Returns dPhi/dx
double Phix(ld x, ld y, ld rs) {
    return 0;
    ld r = sqrtl(x*x+y*y);
    if (r <= rs) return 0;
    return x/(r*(r-rs)*(r-rs));
}

//Returns dPhi/dy
double Phiy(ld x, ld y, ld rs) {
    //return 0.0;
    ld r = sqrtl(x*x+y*y);
    if (r == rs) return 0;
    return y/(r*(r-rs)*(r-rs));
}

double Phi_Newt(ld x, ld y) {
    ld r = sqrtl(x*x+y*y);
    if (r == 0) return 0;
    return -1.0/r;
}

double Phi_Rel(ld x, ld y, ld rs) {
    ld r = sqrtl(x*x+y*y);
    if (r <= rs) return 0;
    return -1.0/(r-rs);
}


bool toPrimitive(ld *u, ld *w) {
    w[0] = u[0]; w[1] = u[1]/u[0]; w[2] = u[2]/u[0]; w[3] = u[3]/u[0];
    w[5] = u[5]; w[6] = u[6]; w[7] = u[7];

    ld v2 = w[1]*w[1]+w[2]*w[2]+w[3]*w[3];
    ld B2 = w[5]*w[5]+w[6]*w[6]+w[7]*w[7];
    w[4] = (u[4] - w[0]*v2/2 - B2/2) * (gamma-1);
    if (w[4] <= 0) {
        std::cout << "toPrimitive::The pressure is negative or null !!!" << std::endl; return false;
    }

    for (int i = 0; i < 8; i++){
        if (isnan(w[i])) { cout << "toPrimitive:: One of the values is NaN" << endl; return false; }
    }

    return true;
}

bool toConserved(ld *w, ld *u) {
    ld vx = w[1]; ld vy = w[2]; ld vz = w[3];
    u[0] = w[0]; u[1] = w[0]*vx; u[2] = w[0]*vy; u[3] = w[0]*vz;
    u[5] = w[5]; u[6] = w[6]; u[7] = w[7];

    ld v2 = vx*vx+vy*vy+vz*vz;
    ld B2 = w[5]*w[5]+w[6]*w[6]+w[7]*w[7];

    u[4] = w[4]/(gamma-1) + w[0]*v2/2 + B2/2;

    for (int i = 0; i < 8; i++) {
        if (isnan(u[i])) { cout << "toConserved:: One of the values is NaN" << endl; return false; }
    }


    return true;
}

bool F(ld u[8], ld *flux) {
    //rho, Mx, My, Mz, E, Bx, By, Bz
    if (u[0] <= 0)
        return false;
    ld vx = u[1]/u[0]; ld vy = u[2]/u[0]; ld vz = u[3]/u[0];
    ld v = sqrtl(vx*vx + vy*vy + vz*vz);
    ld B = sqrtl(u[5]*u[5] + u[6]*u[6] + u[7]*u[7]);
    ld P = (u[4] - u[0]*v*v/2 - B*B/2) * (gamma - 1);
    if (P <= 0) {
        std::cout << "FUNCTION::ERROR::The pressure is negative or null: " << P << std::endl;
        return false;
    }

    ld Pstar = P + B*B/2;

    flux[0] = u[1];
    flux[1] = u[1]*vx + Pstar - u[5]*u[5];
    flux[2] = u[1]*vy - u[5]*u[6];
    flux[3] = u[1]*vz - u[5]*u[7];
    flux[4] = (u[4]+Pstar)*vx - (u[5]*vx+u[6]*vy+u[7]*vz)*u[5];
    flux[5] = 0;
    flux[6] = u[6]*vx - u[5]*vy;
    flux[7] = u[7]*vx - u[5]*vz;

    return true;
}

bool G(ld u[8], ld *flux) {
    //0=rho, 1=Mx, 2=My, 3=Mz, 4=E, 5=Bx, 6=By, 7=Bz
    if (u[0] <= 0)
        return false;
    ld vx = u[1]/u[0]; ld vy = u[2]/u[0]; ld vz = u[3]/u[0];
    ld v = sqrtl(vx*vx + vy*vy + vz*vz);
    ld B = sqrtl(u[5]*u[5] + u[6]*u[6] + u[7]*u[7]);
    ld P = (u[4] - u[0]*v*v/2 - B*B/2) * (gamma - 1);
    if (P <= 0) {
        std::cout << "FUNCTION::ERROR::The pressure is negative or null: " << P << std::endl;
        return false;
    }

    ld Pstar = P + B*B/2;

    flux[0] = u[2];
    flux[1] = u[2]*vx - u[6]*u[5];
    flux[2] = u[2]*vy + Pstar - u[6]*u[6];
    flux[3] = u[2]*vz - u[6]*u[7];
    flux[4] = (u[4]+Pstar)*vy - (u[5]*vx+u[6]*vy+u[7]*vz)*u[6];
    flux[5] = u[5]*vy - u[6]*vx;
    flux[6] = 0;
    flux[7] = u[7]*vy - u[6]*vz;

    return true;
}

enum bound
{
    PERIODIC,
    OUTFLOW,
    REFLECTING,
    PRESSURE,
    NONE
};


void setBoundaries(bound xbounds, bound ybounds, Arrays *us, int Nx, int Ny, ld gy) {
    ld y0 = us->y0; ld yn = us->yn; ld dy = us->dy;

    if (xbounds == bound::OUTFLOW) {
        for (int yi = 0; yi < Ny+8; yi++) {
        for (int i = 0; i < 8; i++) {
            
            us->at(i, 0, yi) = us->at(i, 4, yi); us->at(i, 1, yi) = us->at(i, 4, yi); us->at(i, 2, yi) = us->at(i, 4, yi); us->at(i, 3, yi) = us->at(i, 4, yi);
            us->at(i, Nx+7, yi) = us->at(i, Nx+3, yi); us->at(i, Nx+6, yi) = us->at(i, Nx+3, yi); us->at(i, Nx+5, yi) = us->at(i, Nx+3, yi); us->at(i, Nx+4, yi) = us->at(i, Nx+3, yi);
        }
        us->at(5, Nx+8, yi) = us->at(5, Nx+3, yi);
        }
    }
    if (xbounds == bound::PERIODIC) {
        for (int yi = 0; yi < Ny+8; yi++) {
        for (int i = 0; i < 8; i++) {
            us->at(i, 0, yi) = us->at(i, Nx, yi); us->at(i, 1, yi) = us->at(i, Nx+1, yi); us->at(i, 2, yi) = us->at(i, Nx+2, yi); us->at(i, 3, yi) = us->at(i, Nx+3, yi);
            us->at(i, Nx+7, yi) = us->at(i, 7, yi);  us->at(i, Nx+6, yi) = us->at(i, 6, yi); us->at(i, Nx+5, yi) = us->at(i, 5, yi); us->at(i, Nx+4, yi) = us->at(i, 4, yi);
        }

        us->at(5, Nx+8, yi) = us->at(5, 8, yi);
        }
    }
    if (xbounds == bound::REFLECTING) {
        for (int yi = 0; yi < Ny+8; yi++) {
        for (int i = 0; i < 8; i++) {
            if (i == 5) continue;
            us->at(i, 0, yi) = us->at(i, 7, yi); us->at(i, 1, yi) = us->at(i, 6, yi); us->at(i, 2, yi) = us->at(i, 5, yi); us->at(i, 3, yi) = us->at(i, 4, yi);
            us->at(i, Nx+7, yi) = us->at(i, Nx, yi); us->at(i, Nx+6, yi) = us->at(i, Nx+1, yi); us->at(i, Nx+5, yi) = us->at(i, Nx+2, yi); us->at(i, Nx+4, yi) = us->at(i, Nx+3, yi);
        }
        for (int gh = 0; gh < 4; gh++)
        {
            us->at(1, gh, yi) = -us->at(1, gh, yi);
            us->at(1, Nx+4+gh, yi) = -us->at(1, Nx+4+gh, yi);
        }
        
        }
        for (int yi = 0; yi < Ny+8; yi++) {
            us->at(5, 0, yi) = us->at(5, 8, yi); us->at(5, 1, yi) = us->at(5, 7, yi); us->at(5, 2, yi) = us->at(5, 6, yi); us->at(5, 3, yi) = us->at(5, 5, yi);
            us->at(5, Nx+8, yi) = us->at(5, Nx, yi); us->at(5, Nx+7, yi) = us->at(5, Nx+1, yi); us->at(5, Nx+6, yi) = us->at(5, Nx+2, yi); us->at(5, Nx+5, yi) = us->at(5, Nx+3, yi);
        }

    }

    if (ybounds == bound::OUTFLOW) {
        for (int xi = 0; xi < Nx+8; xi++){
        for (int i = 0; i < 8; i++) {
            us->at(i, xi, 3) = us->at(i, xi, 4); us->at(i, xi, 2) = us->at(i, xi, 4); us->at(i, xi, 1) = us->at(i, xi, 4); us->at(i, xi, 0) = us->at(i, xi, 4);
            us->at(i, xi, Ny+7) = us->at(i, xi, Ny+3); us->at(i, xi, Ny+6) = us->at(i, xi, Ny+3); us->at(i, xi, Ny+5) = us->at(i, xi, Ny+3); us->at(i, xi, Ny+4) = us->at(i, xi, Ny+3);
        }
            
        us->at(6, xi, Ny+8) = us->at(6, xi, Ny+3);
        }
    }
    if (ybounds == bound::PERIODIC) {
        for (int xi = 0; xi < Nx+8; xi++){
        for (int i = 0; i < 8; i++) {
            us->at(i, xi, 3) = us->at(i, xi, Ny+3); us->at(i, xi, 2) = us->at(i, xi, Ny+2); us->at(i, xi, 1) = us->at(i, xi, Ny+1); us->at(i, xi, 0) = us->at(i, xi, Ny);
            us->at(i, xi, Ny+7) = us->at(i, xi, 7); us->at(i, xi, Ny+6) = us->at(i, xi, 6); us->at(i, xi, Ny+5) = us->at(i, xi, 5); us->at(i, xi, Ny+4) = us->at(i, xi, 4);
        }
        
        us->at(6, xi, Ny+8) = us->at(6, xi, 8);
        }
    }
    if (ybounds == bound::REFLECTING) {
        for (int xi = 0; xi < Nx+8; xi++){
        for (int i = 0; i < 8; i++) {
            if (i == 6) continue;
            us->at(i, xi, 0) = us->at(i, xi, 7); us->at(i, xi, 1) = us->at(i, xi, 6); us->at(i, xi, 2) = us->at(i, xi, 5); us->at(i, xi, 3) = us->at(i, xi, 4);
            us->at(i, xi, Ny+7) = us->at(i, xi, Ny); us->at(i, xi, Ny+6) = us->at(i, xi, Ny+1); us->at(i, xi, Ny+5) = us->at(i, xi, Ny+2); us->at(i, xi, Ny+4) = us->at(i, xi, Ny+3);
        }
        for (int gh = 0; gh < 4; gh++)
        {
            us->at(2, xi, gh) = -us->at(2, xi, gh);
            us->at(2, xi, Ny+4+gh) = -us->at(2, xi, Ny+4+gh);
        }
        
        
        }

        for (int xi = 0; xi < Nx+8; xi++){
            us->at(6, xi, 0) = us->at(6, xi, 8); us->at(6, xi, 1) = us->at(6, xi, 7); us->at(6, xi, 2) = us->at(6, xi, 6); us->at(6, xi, 3) = us->at(6, xi, 5);
            us->at(6, xi, Ny+8) = us->at(6, xi, Ny); us->at(6, xi, Ny+7) = us->at(6, xi, Ny+1); us->at(6, xi, Ny+6) = us->at(6, xi, Ny+2); us->at(6, xi, Ny+5) = us->at(6, xi, Ny+3);
        }
    }


    if (ybounds == bound::PRESSURE) {
        for (int xi = 0; xi < Nx+8; xi++){
        
            for (int gh = 1; gh <= 4; gh++)
            {
                /*  At the y0 boundary  */
                ld wout[8] = { 0 }; /*  Ghost cell  */ 
                ld win[8]  = { 0 }; for (int i = 0; i < 8; i++) win[i]  = us->at(i, xi, 4+gh-1); /*  Corresponding cell*/
                win[5] = (us->at(5, xi, 4+gh-1)+us->at(5, xi+1, 4+gh-1))/2.0;   win[6] = (us->at(6, xi, 4+gh-1)+us->at(6, xi, 4+gh))/2.0; /*  Compute the cell-centered magnetic fields  */

                toPrimitive(win, win); 
                for (int i = 0; i < 8; i++)
                    wout[i] = win[i];

                wout[2] = -win[2]; /*  Reflect the y-velocity  */
                wout[6] = -win[6]; /*  Reflect the cell-centered By  */
                wout[4] = win[4] + win[0]*gy*(2.0*gh-1.0)*dy; /*  Project the pressure to the ghost cells  */
                toConserved(wout, wout);

                wout[5] = us->at(5, xi, 4+gh-1);    wout[6] = -us->at(6, xi, 4+gh); /*  Pick the correct Bx/By and reflect By.  */

                for (int i = 0; i < 8; i++) us->at(i, xi, 4-gh) = wout[i];


                /*  At the yn boundary  */
                for (int i = 0; i < 8; i++) win[i] = us->at(i, xi, Ny+3-gh+1); /*  Corresponding cell  */
                win[5] = (us->at(5, xi, Ny+3-gh+1)+us->at(5, xi+1, Ny+3-gh+1))/2.0;   win[6] = (us->at(6, xi, Ny+3-gh+1)+us->at(6, xi, Ny+3-gh+2))/2.0; /*  Compute cell-centered Bx & By  */

                toPrimitive(win, win);
                for (int i = 0; i < 8; i++)
                    wout[i] = win[i];

                wout[2] = -win[2]; /*  Reflect the y-velocity  */
                wout[6] = -win[6]; /*  Reflect the cell-centered By  */
                wout[4] = win[4] - win[0]*gy*(2.0*gh-1.0)*dy;
                toConserved(wout, wout);

                wout[5] = us->at(5, xi, Ny+3-gh+1);

                for (int i = 0; i < 8; i++) us->at(i, xi, Ny+3+gh) = wout[i];

                us->at(6, xi, Ny+3+gh+1) = -us->at(6, xi, Ny+3-gh+1);

            }
        
        
        }
    }

    
}