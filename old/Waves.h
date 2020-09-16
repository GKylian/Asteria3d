#pragma once
#include <iostream>
#define _USE_MATH_DEFINES
#include <math.h>
using namespace std;


#define M_rt4PI 3.54490770181103205460
//const double gamma = 1.4;

//TODO://////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//TODO://////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//TODO://////////////////---  WAVESPEEDS  ---////////////////////////////////////////////////////////////////////////////////
//TODO://////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool getXWaveSpeeds(long double u[8], long double *waves) { //--OK--

    for (int i = 0; i < 8; i++) {
        if (isnan(u[i])) {
            cout << "XWAVES:: The value at index " << i << " is NaN !" << endl; return false;
        }
    }
    if (u[0] <= 0) {
        cout << "X-WAVES::ERROR::The density is null or negative" << endl;
        return false;
    }
        

    /// 0. Variables used later
    long double v = sqrtl(u[1]/u[0]* u[1]/u[0] + u[2]/u[0]* u[2]/u[0] + u[3]/u[0]* u[3]/u[0]);
    long double B = sqrtl(u[5]*u[5] + u[6]*u[6] + u[7]*u[7]);// long double b = B/M_rt4PI;

    long double P = (u[4] - u[0]*v*v/2 - B*B/2) * (gamma - 1);
    
    

    if (P <= 0) {
        cout << "X-WAVES::WARNING::The pressure is negative or null: " << P << endl;
        P = 0;
        return false;
    } 


    /// 1. Sound wave speed:
    long double a = sqrtl(gamma*P/u[0]);
    waves[0] = a;

    /// 3. Alfvén wave speed:
    long double CA = sqrtl(B*B/u[0]);
    waves[2] = CA;

    /// 2/4. Slow/Fast magnetosonic wave speed
    long double CAx = sqrtl(u[5]*u[5]/u[0]);
    waves[1] = sqrtl(0.5*((a*a+CA*CA) - sqrtl(powl(a*a+CA*CA, 2) - 4*a*a*CAx*CAx)));
    waves[3] = sqrtl(0.5*((a*a+CA*CA) + sqrtl(powl(a*a+CA*CA, 2) - 4*a*a*CAx*CAx)));

    int nans = isnan(waves[0])+isnan(waves[1])+isnan(waves[2])+isnan(waves[3]);
    if (nans != 0) {
        cout << "XWAVES::" << nans << " of the wavespeeds is/are NaN !!! " << endl;
        cout << "     --> " << waves[0] << ",  " << waves[1] << ",  " << waves[2] << ",  " << waves[3] << endl;
        cout << "a^2: " << gamma*P/u[0] << endl << endl;
        return false;
    }

    return true;
}

bool getYWaveSpeeds(long double u[8], long double *waves) {

    if (u[0] <= 0) {
        cout << "Y-WAVES::ERROR::The density is null or negative" << endl;
        return false;
    }
        

    /// 0. Variables used later
    long double v = sqrtl(u[1]/u[0]* u[1]/u[0] + u[2]/u[0]* u[2]/u[0] + u[3]/u[0]* u[3]/u[0]);
    long double B = sqrtl(u[5]*u[5] + u[6]*u[6] + u[7]*u[7]); //long double b = B/M_rt4PI;
    long double P = (u[4] - u[0]*v*v/2 - B*B/2) * (gamma - 1);
    if (P <= 0) {
        cout << "Y-WAVES::WARNING::The pressure is negative or null: " << P << endl;
        P = 0;
        return false;
    } 


    /// 1. Sound wave speed:
    long double a = sqrtl(gamma*P/u[0]);
    waves[0] = a;

    /// 3. Alfvén wave speed:
    long double CA = sqrtl(B*B/u[0]);
    waves[2] = CA;

    /// 2/4. Slow/Fast magnetosonic wave speed
    long double CAy = sqrtl(u[6]*u[6]/u[0]);
    waves[1] = sqrtl(0.5*((a*a+CA*CA) - sqrtl(powl(a*a+CA*CA, 2) - 4*a*a*CAy*CAy)));
    waves[3] = sqrtl(0.5*((a*a+CA*CA) + sqrtl(powl(a*a+CA*CA, 2) - 4*a*a*CAy*CAy)));

    if (isnan(waves[0]) || isnan(waves[1]) || isnan(waves[2]) || isnan(waves[3])) {
        cout << "YWAVES::One of the wavespeeds is NaN !!!" << endl;
        cout << "     --> " << waves[0] << ",  " << waves[1] << ",  " << waves[2] << ",  " << waves[3] << endl;
        return false;
    }

    return true;
}
