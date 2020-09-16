// 2dMHD.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <string>
#include "Arrays.h"
#include <vector>
#include <algorithm>
#include <chrono>
#include <fstream>
#include "phymaths.h"
#include "Eigenvalues.h"
#include "Waves.h"
#include "Reconstruction.h"
#include "Riemann.h"
#include "InitCond.h"


using namespace std;
using namespace std::chrono;

//const double gamma = 1.4;









//Computes the fluxes at i-1/2 and j-1/2
bool getFluxes(Arrays *u, Arrays *ix_u, Arrays *iy_u, int xi, int yi, long double *F_m, long double *G_m, int lvl) { // --OK-- (?!?!?!)


    lvl = 2;

    /////0. A few things
    long double ui[8] = { 0 };
    long double u_xL[8] = { 0 }; long double u_xR[8] = { 0 };
    long double u_yL[8] = { 0 }; long double u_yR[8] = { 0 };
    

    for (int i = 0; i < 8; i++)
    {
        u_xL[i] = ix_u->at(i, xi*2-1, yi); //L_m of xi = L_p of xi-1
        u_xR[i] = ix_u->at(i, xi*2, yi); //R_m of xi

        u_yL[i] = iy_u->at(i, xi, yi*2-1);
        u_yR[i] = iy_u->at(i, xi, yi*2);

        ui[i] =  u->at(i, xi, yi);
    }
    
    
    
    //!!! Bx and By are stored at (i-1/2, j) and (i, j-1/2) respectively
    
    //Bx and By at their next respective interface (Bx at i+1/2,j and By at i,j+1/2)
    long double Bx_ix = u->at(5, xi+1, yi); long double By_iy = u->at(6, xi, yi+1); //!!! They are not cell-centered but face-centered 


    //Make ui cell-centered Bx and By:
    ui[5] = (ui[5]+Bx_ix)/2; ui[6] = (ui[6]+By_iy)/2;

    //Bx and By at other interface (Bx at i,j-1/2 and By at i-1/2,j)
    /*long double By_ix = ( ui[6] + (u->at(6, xi-1, yi+1)+u->at(6, xi-1, yi))/2 )/2;
    long double Bx_iy = ( ui[5] + (u->at(5, xi+1, yi-1)+u->at(5, xi, yi-1))/2 )/2;*/

    //Values at interfaces L/R for Bx and By are the same (B_i-1/2 L/R = B_i-1/2):
    //u_xL[5] = ui[5]; u_xL[6] = By_ix; u_xR[5] = ui[5]; u_xR[6] = By_ix;
    //u_yL[5] = Bx_iy; u_yL[6] = ui[6]; u_yR[5] = Bx_iy; u_yR[6] = ui[6];

    //u_xL[5] = Bx_ix; u_xR[5] = Bx_ix;
    //u_yL[6] = By_iy; u_yR[6] = By_iy;
    /*u_xL[5] = ui[5]; u_xR[5] = ui[5];
    u_yL[6] = ui[6]; u_yR[6] = ui[6];*/



    if (isnan(u_xL[5]) || isnan(u_yL[5]) || isnan(u_xL[6]) || isnan(u_yL[6])) {
        cout << "X-getFluxes:: one of the interface magnetic fields is NaN !" << endl; return false;
    }



    long double F_x[8] = { 0 }; long double G_y[8] = { 0 };
    //1. Compute F_x at i-1/2,j, so Bx=Bx_ix and By=By_ix
    if (true) {

        if (u_xL[5] != u_xR[5])
            cout << "ERROR::2DMHD.CPP::getFluxes:: Bx is different at left and right interfaces." << endl;
        //Find max wavespeeds:
        long double waves_L[4] = { 0 }; if (!getXWaveSpeeds(u_xL, waves_L)) { cout << "X-getFluxes::XWaves left returned false" << endl; return false; }
        long double waves_R[4] = { 0 }; if(!getXWaveSpeeds(u_xR, waves_R)) { cout << "X-getFluxes::XWaves right returned false" << endl; return false; }
        long double c_L = *max_element(waves_L, waves_L+4); long double c_R = *max_element(waves_R, waves_R+4);


        //Find eigenvalues (from smallest to largest):
        long double eigen[7] = { 0 }; if(!getXEigenvalues(u_xL, u_xR, eigen)) return false;
        long double lambda_m = *max_element(eigen, eigen+7); long double lambda_0 = *min_element(eigen, eigen+7);
        


        //long double SR = fmaxl(lambda_m, u_xR[1]/u_xR[0] + c_R); //b_plus
        //long double SL = fminl(lambda_0, u_xL[1]/u_xL[0] - c_L); //b_min


        long double SL = fminl(u_xL[1]/u_xL[0], u_xR[1]/u_xR[0]) - fmaxl(c_L, c_R);
        long double SR = fmaxl(u_xL[1]/u_xL[0], u_xR[1]/u_xR[0]) + fmaxl(c_L, c_R);



        long double F_L[8] = { 0 }; if(!F(u_xL, F_L)) return false;
        long double F_R[8] = { 0 }; if(!F(u_xR, F_R)) return false;


        //We include Bx and By, but it's not important, we just won't use them

        if (lvl == 0) { //Preferred choice: HLLD
            if (!X_HLLDFlux(u_xL, u_xR, SL, SR, F_L, F_R, F_m)) return false;
        }
        else if (lvl == 1) { //Intermediate: HLLC
            if (!X_HLLCFlux(u_xL, u_xR, SL, SR, F_L, F_R, F_m)) return false;
        }
        else { //Should always work: HLLE
            if (!X_HLLFlux(u_xL, u_xR, SL, SR, F_L, F_R, F_m)) return false;
            //for (int i = 0; i < 8; i++)
            //{
            //    long double Fi = (SR*F_L[i] - SL*F_R[i] + SR*SL*(u_xR[i]-u_xL[i])) / (SR-SL);
            //    if (SL > 0)
            //        F_m[i] = F_L[i];//u_xL[i];//u->at(i, xi-1, yi);
            //    else if (SR < 0)
            //        F_m[i] = F_R[i];//u_xR[i];//ui[i];
            //    else
            //        F_m[i] = Fi;
            //}
        }


        
    }

    //2. Compute G_x at i,j-1/2 , so Bx=Bx_iy and By=By_iy
    if (true) {
        if (u_yL[6] != u_yR[6])
            cout << "ERROR::2DMHD.CPP::getFluxes:: By is different at left and right interfaces." << endl;

        //Find max wavespeeds:
        long double waves_L[4] = { 0 }; if(!getYWaveSpeeds(u_yL, waves_L)) return false;
        long double waves_R[4] = { 0 }; if(!getYWaveSpeeds(u_yR, waves_R)) return false;
        long double c_L = *max_element(waves_L, waves_L+4); long double c_R = *max_element(waves_R, waves_R+4);

        //Find eigenvalues (from smallest to largest):
        long double eigen[7] = { 0 }; if(!getYEigenvalues(u_yL, u_yR, eigen)) return false;
        long double lambda_m = *max_element(eigen, eigen+7); long double lambda_0 = *min_element(eigen, eigen+7);

        //long double SR = fmaxl(lambda_m, u_yR[2]/u_yR[0] + c_R);
        //long double SL = fminl(lambda_0, u_yL[2]/u_yL[0] - c_L);

        long double SL = fminl(u_yL[2]/u_yL[0], u_yR[2]/u_yR[0]) - fmaxl(c_L, c_R);
        long double SR = fmaxl(u_yL[2]/u_yL[0], u_yR[2]/u_yR[0]) + fmaxl(c_L, c_R);


        long double G_L[8] = { 0 }; if(!G(u_yL, G_L)) return false;
        long double G_R[8] = { 0 }; if(!G(u_yR, G_R)) return false;



        //We include Bx and By, but it's not important, we just won't use them


        if (lvl == 0) { //Preferred choice: HLLD
            if (!Y_HLLDFlux(u_yL, u_yR, SL, SR, G_L, G_R, G_m)) return false;
        }
        else if (lvl == 1) { //Intermediate: HLLC
            if (!Y_HLLCFlux(u_yL, u_yR, SL, SR, G_L, G_R, G_m)) return false;
        }
        else { //Should always work: HLLE
            if (!Y_HLLFlux(u_yL, u_yR, SL, SR, G_L, G_R, G_m)) return false;
            /*for (int i = 0; i < 8; i++)
            {
                long double Gi = (SR*G_L[i] - SL*G_R[i] + SR*SL*(u_yR[i]-u_yL[i])) / (SR-SL);
                if (SL > 0)
                    G_m[i] = G_L[i];
                else if (SR < 0)
                    G_m[i] = G_R[i];
                else
                    G_m[i] = Gi;
            }*/
        }

    }


    //if (F_m[7] != 0.0 || G_m[7] != 0.0)
    //    cout << "F[7] or G[7] is nonzero: " << F_m[7] << "; " << G_m[7] << endl;

    return true;
}


long double refEMF(Arrays *u, int xi, int yi) {
    long double vx = u->at(1, xi, yi)/u->at(0, xi, yi); long double vy = u->at(2, xi, yi)/u->at(0, xi, yi);
    long double Bx = (u->at(5, xi+1, yi)+u->at(5, xi, yi))/2; long double By = (u->at(6, xi, yi+1)+u->at(6, xi, yi))/2;
    
    long double val = -(vx*By - vy*Bx);
    if (isnan(val))
        cout << "refEMF::Value is NaN" << endl;
    
    return (!null(val))*val;
}


//Returns the bottom left corner EMF: Ez at (i-1/2, j-1/2)
long double getEMFCorner(Arrays *u, Arrays *fs, Arrays *gs, int xi, int yi, long double dx, long double dy) {

    /////1. Get reference EMFs at cell centers i(-1), j(-1)
    long double Er_i_j = refEMF(u, xi, yi);
    long double Er_im_j = refEMF(u, xi-1, yi); long double Er_i_jm = refEMF(u, xi, yi-1);
    long double Er_im_jm = refEMF(u, xi-1, yi-1);


    
    /////2. Get face-centered EMF from fluxes at i-1/2, j(-1) from -F and at i(-1), j-1/2 from G
    long double E_jp = -fs->at(6, xi, yi); long double E_jm = -fs->at(6, xi, yi-1); //On the y interface, above and under the corner
    long double E_ip = gs->at(5, xi, yi); long double E_im = gs->at(5, xi-1, yi); //On the x interface, left and right to the corner
    if (null(E_jp)) E_jp = 0.0; if (null(E_jm)) E_jm = 0.0; if (null(E_ip)) E_ip = 0.0; if (null(E_im)) E_im = 0.0;
    /*if (E_jp != 0 || E_jm != 0 || E_ip != 0 || E_im != 0)
        cout << E_jp << ", " << E_jm << ", " << E_ip << ", " << E_im << endl;*/

    /*if (Er_i_j != Er_i_jm || Er_im_j != Er_im_jm) {
        cout << Er_i_j << ", " << Er_i_jm << ",,   " << Er_im_j << ", " << Er_im_jm << endl;
        cout << "Bx: " << (u->at(5, xi+1, yi)+u->at(5, xi, yi))/2 << " vs " << (u->at(5, xi+1, yi+1)+u->at(5, xi, yi+1))/2 << endl;
        cout << "By: " << (u->at(6, xi, yi+1)+u->at(6, xi, yi))/2 << " vs " << (u->at(6, xi, yi+2)+u->at(6, xi, yi+1))/2 << endl;
        cout << "Bx-1: " << (u->at(5, xi, yi)+u->at(5, xi-1, yi))/2 << " vs " << (u->at(5, xi, yi+1)+u->at(5, xi-1, yi+1))/2 << endl;
        cout << "By-1: " << (u->at(6, xi-1, yi+1)+u->at(6, xi-1, yi))/2 << " vs " << (u->at(6, xi-1, yi+2)+u->at(6, xi-1, yi+1))/2 << endl;
        cout << endl;
    }

    if (E_jp != E_jm)
        cout << "E_jp != E_jm: " << E_jp-E_jm << endl;*/
    

    /////3. Compute the four velocities: v_x(i-1/2, j(-1)) and v_y(i(-1), j-1/2)
    //TODO: Comment trouver les valeurs des vitesses aux interfaces ?
    long double vx_jm = (u->at(1, xi, yi-1)/u->at(0, xi, yi-1) + u->at(1, xi-1, yi-1)/u->at(0, xi-1, yi-1))/2; long double vx_jp = (u->at(1, xi, yi)/u->at(0, xi, yi) + u->at(1, xi-1, yi)/u->at(0, xi-1, yi))/2;
    long double vy_im = (u->at(2, xi-1, yi)/u->at(0, xi-1, yi) + u->at(2, xi-1, yi-1)/u->at(0, xi-1, yi-1))/2; long double vy_ip = (u->at(2, xi, yi)/u->at(0, xi, yi) + u->at(2, xi, yi-1)/u->at(0, xi, yi-1))/2;
    long double rvx_jm = fs->at(0, xi, yi-1); long double rvx_jp = fs->at(0, xi, yi);
    long double rvy_im = gs->at(0, xi-1, yi); long double rvy_ip = gs->at(0, xi, yi);


    //3.1) wrt to y at i-1/2, j-1/4 (so jp):
    long double dEdy_jp = 0;
    if (rvx_jp > 0)
        dEdy_jp = (Er_im_j - E_im)*2/dy;
    else if (rvx_jp < 0)
        dEdy_jp = (Er_i_j - E_ip)*2/dy;
    else
        dEdy_jp = 0.5*( (Er_im_j - E_im)*2/dy + (Er_i_j - E_ip)*2/dy );
    if (null(dEdy_jp)) dEdy_jp = 0.0;

    //3.2) wrt to y at i-1/2,j-3/4 (so jm):
    long double dEdy_jm = 0;
    if (rvx_jm > 0) // dE/dy at (i-1, j-3/4) --> 2/dy * ( E(i-1, j-1/2) - E(i-1, j-1) )
        dEdy_jm = (E_im - Er_im_jm)*2/dy;
    else if (rvx_jm < 0)
        dEdy_jm = (E_ip - Er_i_jm)*2/dy;
    else
        dEdy_jm = 0.5*( (E_im - Er_im_jm)*2/dy + (E_ip - Er_i_jm)*2/dy );
    if (null(dEdy_jm)) dEdy_jm = 0.0;



    //3.3) wrt to x at i-1/4 (so ip)
    long double dEdx_ip = 0;
    if (rvy_ip > 0)
        dEdx_ip = (Er_i_jm - E_jm)/(dx/2);
    else if (rvy_ip < 0)
        dEdx_ip = (Er_i_j - E_jp)/(dx/2);
    else
        dEdx_ip = 0.5*( (Er_i_jm - E_jm)/(dx/2) + (Er_i_j - E_jp)/(dx/2) );
    if (null(dEdx_ip)) dEdx_ip = 0.0;

    //3.4) wrt to x at i-3/4 (so im)
    long double dEdx_im = 0;
    if (rvy_im > 0)
        dEdx_im = (E_jm - Er_im_jm)/(dx/2);
    else if (rvy_im < 0)
        dEdx_im = (E_jp - Er_im_j)/(dx/2);
    else
        dEdx_im = 0.5*( (E_jm - Er_im_jm)/(dx/2) + (E_jp - Er_im_j)/(dx/2) );
    if (null(dEdx_im)) dEdx_im = 0.0;



    /////4. Compute the EMF at corner i-1/2, j-1/2
    long double value = (long double) 0.25*(E_jp + E_jm + E_ip + E_im) + dy/8.0 * (dEdy_jm - dEdy_jp) + dx/8 * (dEdx_im - dEdx_ip);
    /*if (value != 0)
        cout << "Only a sith deals in absolutes " << value << endl;*/
    if (fabsl(value) < 1e-12) value = 0;

    return value;

}

bool getSourceTerms(Arrays *u, int xi, int yi, long double *sx, long double *sy, long double dx, long double dy) {
    long double dBx = (u->at(5, xi+1, yi) - u->at(5, xi, yi))/dx;
    long double dBy = (u->at(6, xi, yi+1) - u->at(6, xi, yi))/dy;

    long double Bxi = (u->at(5, xi+1, yi) + u->at(5, xi, yi))/2;
    long double Byi = (u->at(6, xi, yi+1) + u->at(6, xi, yi))/2;
    long double Bzi = u->at(7, xi, yi);

    long double vz = u->at(3, xi, yi)/u->at(0, xi, yi);

    sx[0] = 0;
    sx[1] = Bxi*dBx; sx[2] = Byi*dBx; sx[3] = Bzi*dBx;
    sx[4] = Bzi*vz*dBx; sx[5] = 0;
    sx[6] = 0; sx[7] = vz*dBx;

    sy[0] = 0;
    sy[1] = Bxi*dBy; sy[2] = Byi*dBy; sy[3] = Bzi*dBy;
    sy[4] = Bzi*vz*dBy; sy[5] = 0;
    sy[6] = 0; sy[7] = vz*dBy;

    return true;
}






//TODO:-----------------------------------------------------------------------------------------------------------------------------
//TODO:-----------------------------------------------------------------------------------------------------------------------------
//TODO:-----------------------------------------------------------------------------------------------------------------------------
//TODO:-----------------------------------------------------------------------------------------------------------------------------
//TODO:-----------------------------------------------------------------------------------------------------------------------------
//TODO:-----------------------------------------------------------------------------------------------------------------------------
//TODO:-----------------------------------------------------------------------------------------------------------------------------
//TODO:-----------------------------------------------------------------------------------------------------------------------------
//TODO:-----------------------------------------------------------------------------------------------------------------------------
//TODO:-----------------------------------------------------------------------------------------------------------------------------
//TODO:-----------------------------------------------------------------------------------------------------------------------------
//TODO:-----------------------------------------------------------------------------------------------------------------------------
//TODO:-----------------------------------------------------------------------------------------------------------------------------
//TODO:-----------------------------------------------------------------------------------------------------------------------------
//TODO:-----------------------------------------------------------------------------------------------------------------------------
//TODO:-----------------------------------------------------------------------------------------------------------------------------




int main()
{
    std::cout << "Hello World!\n";

    long double gx = 0.0; long double gy = 0.0; long double gz = 0.0; long double rs = 2.0; //-1: constant g,  0: newtonian potential,  1+: pseudo-newtonian potential
    int k = 100; long double CFL = 0.3; int smax = 50000;

    double x0 = 4.0; double xn = 44.0; int Nx = 128; double dx = (xn-x0)/(Nx-1.0);
    double y0 = -20.0; double yn = 20.0; int Ny = 128; double dy = (yn-y0)/(Ny-1.0);
    double t0 = 0.0; double tn = 1000.0; double dt0 = 0.000000001; double dt = dt0;


    bool HD = true; bool CT = true; bool CharTracing = true; bool simpleReconstr = false;
    bound xbound = bound::OUTFLOW;
    bound ybound = bound::OUTFLOW;

    //Number of cells: Nx ; Ny
    //Number of interfaces: (Nx+8)*2 ; (Ny+8)*2

    //Bx and By are stored at (i+1/2, j) and (i, j+1/2) respectively
    //  The arrays for the values, interfaces, fluxes, ... cover the whole grid (0 to Nx+5), but not everything is computed (not possible
    //  to compute interfaces for first and last cells for example). That way, we don't have to include offsets in every part of the algorithm.
    //FIXME: They are now at (i-1/2, j) and (i, j-1/2)
    
    //TODO: Change the boundary conditions for Bx and By

    cout << "gamma: " << gamma << endl;

    Arrays us; us.initAll(Nx+8, Ny+8); us.seth(dx, dy);
    if (!Disk(&us, x0, y0, rs)) {
        cout << "Initialization failed." << endl; return -3;
    }
    setBoundaries(xbound, ybound, &us, Nx, Ny, gy);

    Arrays _us; _us = us;



    
    
    //The interfaces are stored at i-1/2 and j-1/2 respectively
    Arrays ix_us; ix_us.initAll((Nx+8)*2, Ny+8);
    Arrays iy_us; iy_us.initAll(Nx+8, (Ny+8)*2);

    //1d fluxes of the conserved variables ("predictor") at time n at i-1/2 and j-1/2 respectively
    Arrays fs; fs.initAll(Nx+8, Ny+8); Arrays gs; gs.initAll(Nx+8, Ny+8);

    ofstream times("out/times.csv");// times << to_string(t0);
    long double ts[6] = { 0 };

    int s = 0; int level = 0;
    for (double t = t0; t < tn; t+=dt)
    {
        /*for (int yi = 0; yi < Ny+8; yi++)
            for (int xi = 0; xi < Nx+8; xi++)
            {
                for (int i = 0; i < 8; i++)
                    if (us.at(i, xi, yi) != _us.at(i, xi, yi)) cout << "fiusfiues" << endl;
            }*/
        
        if(level == 0)
            cout << s << ": Computing from time " << t << " with time-step " << dt << endl;

        if (s == smax) {
            cout << "Reaching the end of the iteration domain s_max = " << smax << endl;
            goto stop;
        }
        s++;

        double nextdt = 10.0;
        long double maxDiv = 0;
        long double maxu[8] = { 0 };
        long double maxw[8] = { 0 };
        
        if (level > 0) {
            if (level > 2) {
                cout << "We've tried everything..." << endl; goto stop;
            }

            //We've already tried this iteration but failed
            cout << "Trying again with level " << level << endl;
        }
        auto start = high_resolution_clock::now();
        auto stop = high_resolution_clock::now();

        auto globstart = high_resolution_clock::now();
        auto globstop = high_resolution_clock::now();



        if(t == t0)
            goto div; //Compute the inital divergence, check the coherence of the initial conditions and compute the initial timestep

        /// 0. Reconstruct values at left and right interfaces from cell-centered values.

        start = high_resolution_clock::now();
        for (int yi = 1; yi < Ny+7; yi++)
        for (int xi = 1; xi < Nx+7; xi++)
        {
            //if (!PPMInterfaces_Prim(&us, xi, yi, &ix_us, &iy_us)) {
            if (HD) {
                //if (!PLMInterfaces_Prim(&us, xi, yi, &ix_us, &iy_us, dt))
                if (!HDPLM(&us, xi, yi, &ix_us, &iy_us, dt, gx, gy, gz, rs))
                    goto stop;
            }
            else {
                if (simpleReconstr) {
                    if (!PLM_Simple(&us, xi, yi, &ix_us, &iy_us, dt, gx, gy, gz, CharTracing)) {
                        cout << "PLM_Simple at (" << xi << ", " << yi << ") failed!" << endl;
                        goto stop;
                    }
                } else {
                    if (!PLM(&us, xi, yi, &ix_us, &iy_us, dt, gx, gy, gz, CharTracing, rs)) {
                        cout << "PLM at (" << xi << ", " << yi << ") failed!" << endl;
                        goto stop;
                    }
                }

                
            }
                
            
            
            
        }
        stop = high_resolution_clock::now();
        ts[0] += duration<double>(stop-start).count();
        cout << "LOOP:: Reconstructed interfaces from cell centers" << endl;

        

        /// 1. Compute and store 1d fluxes at interfaces
        

        start = high_resolution_clock::now();
        for (int yi = 2; yi < Ny+7; yi++)
        for (int xi = 2; xi < Nx+7; xi++)
        {
            
            long double Fm[8] = { 0 }; long double Gm[8] = { 0 };
            if(!getFluxes(&us, &ix_us, &iy_us, xi, yi, Fm, Gm, 0))
                if (!getFluxes(&us, &ix_us, &iy_us, xi, yi, Fm, Gm, 1))
                    if (!getFluxes(&us, &ix_us, &iy_us, xi, yi, Fm, Gm, 2)) {
                        level += 1; goto retry;
                    }
            


            for (int i = 0; i < 8; i++)
            {
                fs.at(i, xi, yi) = Fm[i]; gs.at(i, xi, yi) = Gm[i];
                if (isnan(fs.at(i, xi, yi)) || isnan(gs.at(i, xi, yi)))
                    cout << "One of the fluxes at n is NaN !!!" << endl;
            }
             
        }
        stop = high_resolution_clock::now();
        ts[1] += duration<double>(stop-start).count();
        
        cout << "LOOP:: Computed the fluxes" << endl;

        







        start = high_resolution_clock::now();

        if (!HD && CT){
        //6. Compute EMFs at cell corners, update rho, M, E and Bz to n+1, update Bx and By to n+1
        for (int yi = 4; yi < Ny+5; yi++)
        for (int xi = 4; xi < Nx+5; xi++)
        {
            
            //5.1) Compute EMFs at cell corners
            long double E_im_jm = getEMFCorner(&us, &fs, &gs, xi, yi, dx, dy);
            long double E_ip_jm = getEMFCorner(&us, &fs, &gs, xi+1, yi, dx, dy);
            long double E_im_jp = getEMFCorner(&us, &fs, &gs, xi, yi+1, dx, dy);
            
            
            //5.3) Update Bx and By:
            _us.at(5, xi, yi) += -dt/dy * (E_im_jp - E_im_jm);
            _us.at(6, xi, yi) += dt/dx * (E_ip_jm - E_im_jm);
            
        }
        } /*  if(!MHD)  */

        for (int yi = 4; yi < Ny+4; yi++)
        for (int xi = 4; xi < Nx+4; xi++)
        {
            
            //5.2) Update the conserved variables from time n to time n+1
            long double x = x0 + (xi-4.0)*dx; long double y = y0 + (yi-4.0)*dy;
            
            long double hrho = us.at(0, xi, yi)  -  dt/(2.0*dx) * (fs.at(0, xi+1, yi) - fs.at(0, xi, yi))  -dt/(2.0*dy) * (gs.at(0, xi, yi+1) - gs.at(0, xi, yi));

            //Gravitational source terms:
            if (rs < 0) {
                _us.at(1, xi, yi) += dt*hrho*gx;
                _us.at(2, xi, yi) += dt*hrho*gy;
                _us.at(3, xi, yi) += dt*hrho*gz;
                _us.at(4, xi, yi) += dt*gy*(gs.at(0, xi, yi+1)+gs.at(0, xi, yi))/2.0;
            } else
            {
                //TODO: Or (Phi(x+dx/2, y, rs)-Phi(x-dx/2, y, rs))/dx
                long double phic = Phi_Rel(x, y, rs); long double phil_x = Phi_Rel(x-dx/2.0, y, rs); long double phir_x = Phi_Rel(x+dx/2.0, y, rs);
                long double phil_y = Phi_Rel(x, y-dy/2.0, rs); long double phir_y = Phi_Rel(x, y+dy/2.0, rs);
                //long double dpdx = (Phix(x+dx/2.0, y, rs)-Phix(x-dx/2.0, y, rs))/dx; long double dpdy = (Phiy(x, y+dy/2.0, rs)-Phiy(x, y-dy/2.0, rs))/dy;
                long double Mx = (fs.at(0, xi+1, yi)+fs.at(0, xi, yi))/2.0;
                long double My = (gs.at(0, xi, yi+1)+gs.at(0, xi, yi))/2.0;

                _us.at(1, xi, yi) -= dt/dx*hrho*(phir_x-phil_x);
                _us.at(2, xi, yi) -= dt/dx*hrho*(phir_y-phil_y);
                _us.at(4, xi, yi) -= dt/dx*(fs.at(0, xi, yi)*(phic-phil_x) + fs.at(0, xi+1, yi)*(phir_x-phic));
                _us.at(4, xi, yi) -= dt/dy*(gs.at(0, xi, yi)*(phic-phil_y) + gs.at(0, xi, yi+1)*(phir_y-phic));
            }
            
              

            for (int i = 0; i < 8; i++) {
                if (i == 5 || i == 6 && CT) continue;
                _us.at(i, xi, yi) += -dt/dx * (fs.at(i, xi+1, yi) - fs.at(i, xi, yi)) -dt/dy * (gs.at(i, xi, yi+1) - gs.at(i, xi, yi));// *dx *dx;
                
            }

            if (HD && (_us.at(5, xi, yi) != 0.0 || _us.at(6, xi, yi) != 0.0 || _us.at(7, xi, yi) != 0.0)) {
                cout << "Selected hydrodynamics, but B != 0" << endl; goto stop;
            }
            
            long double u[8] = { 0 }; for (int i = 0; i < 8; i++) u[i] = _us.at(i, xi, yi);
            if (!toPrimitive(u, u)) {
                cout << "Transformation to primitive after update at (" << xi << ", " << yi << ") failed !" << endl;
                for (int i = 0; i < 8; i++) cout << u[i] << ", "; cout << "..." << endl;
                goto stop;
            }
                
            //5.4) Check if anyting is NaN
            for (int i = 0; i < 8; i++)
            {
                if (isnan(us.at(i, xi, yi)))
                    cout << "Component " << i << " of the conserved variables at n+1 is NaN !!!" << endl;
            }
            if (us.at(0, xi, yi) <= 0)
                cout << "Rho at n+1 is either null or negative !!!" << endl;

            /*if (us.at(5, xi, yi) != 0 || us.at(6, xi, yi) != 0 || us.at(7, xi, yi) != 0)
                cout << "One of the magnetic fields is non-zero" << endl;*/
        }
        stop = high_resolution_clock::now();
        ts[2] += duration<double>(stop-start).count();

        cout << "LOOP:: Updated to n+1" << endl;

        
        //Evolution finished, transfer to other array

        ///6. Apply boundary conditions


        start = high_resolution_clock::now();
        setBoundaries(xbound, ybound, &_us, Nx, Ny, gy);
        us = _us;
        stop = high_resolution_clock::now();
        ts[3] += duration<double>(stop-start).count();
        


        

        div:
        ///7. Check that the divergence of the magnetic field is zero 


        for (int yi = 4; yi < Ny+2; yi++)
        for (int xi = 4; xi < Nx+2; xi++)
        {
            
            long double DivB = (us.at(5, xi+1, yi)-us.at(5, xi, yi))/dx + (us.at(6, xi, yi+1)-us.at(6, xi, yi))/dy;
            if (fabsl(DivB) > fabsl(maxDiv)) maxDiv = DivB;
        }
        cout << "Maximum divergence: " << maxDiv << endl;


        start = high_resolution_clock::now();
        ///8. Compute the new time-steps from variables at n+1
        for (int yi = 2; yi < Ny+7; yi++)
        for (int xi = 2; xi < Nx+7; xi++)
        {
            /*if (xi > 60 && xi < 75) {
                cout << xi << ": " << us.at(0, xi, yi) << ", " << us.at(1, xi, yi) << ", " << us.at(4, xi, yi) << endl;
            }
            if (yi > 60 && yi < 75) {
                cout << xi << ": " << us.at(0, xi, yi) << ", " << us.at(2, xi, yi) << ", " << us.at(4, xi, yi) << endl;
            }*/

            //For HD check that the magnetic fields stay null:
            /*if (us.at(5, xi, yi) != 0 || us.at(6, xi, yi) != 0 || us.at(7, xi, yi) != 0) {
                cout << "The magnetic field is non-zero !!!" << endl; goto stop;
            }*/
            //if (us.at(2, xi, yi) != 0 || us.at(3, xi, yi)) cout << "The y or z component of the momentum is non-zero !!!" << endl;
            
            

            long double ui[8]; for (int i = 0; i < 8; i++) ui[i] = us.at(i, xi, yi);
            ui[5] = (us.at(5, xi+1, yi)+us.at(5, xi, yi))/2.0; ui[6] = (us.at(6, xi, yi+1)+us.at(6, xi, yi))/2.0;
            long double wi[8]; toPrimitive(ui, wi);

            for (int i = 0; i < 8; i++)
            {
                if (fabsl(ui[i]) > fabsl(maxu[i])) maxu[i] = ui[i];
                if (fabsl(wi[i]) > fabsl(maxw[i])) maxw[i] = wi[i];
            }


            long double Cxs[4]; if (!getXWaveSpeeds(ui, Cxs)) cout << "dt xwaves" << endl; long double cx = *max_element(Cxs, Cxs+4);
            long double Cys[4]; if (!getYWaveSpeeds(ui, Cys)) cout << "dt ywaves" << endl; long double cy = *max_element(Cys, Cys+4);

            long double vx = fabsl(us.at(1, xi, yi)/us.at(0, xi, yi)); long double vy = fabsl(us.at(2, xi, yi)/us.at(0, xi, yi));
            long double thdt = CFL * fminl(dx/(vx+cx), dy/(vy+cy));
            if (thdt < nextdt)
                nextdt = thdt;
        }


        cout << "Maximum values (conserved): (" << maxu[0];
        for (int i = 1; i < 8; i++)
            cout << ", " << maxu[i];
        cout << ")." << endl;
        
        cout << "Maximum values (primitives): (" << maxw[0];
        for (int i = 1; i < 8; i++)
            cout << ", " << maxw[i];
        cout << ")." << endl;

        //cout << "dt for next time-step: " << nextdt << endl;
        dt = nextdt;
        
        stop = high_resolution_clock::now();
        ts[4] += duration<double>(stop-start).count();

        
        level = 0;//We've completed this iteration without a need to retry --> reset level to 0
        
        
        if (s == 1 || s % k == 0 || t+dt >= tn) {
            
            cout << "Saving data..." << endl;
            string data2out = "out/2data_" + to_string(t) + ".csv";
            
            if(Nx > 30 || Ny > 30)
                times << "," << to_string(t); times.flush();

            if (Nx > 32 && Ny > 32) {
                ofstream out2(data2out); out2 << Nx << "," << x0 << "," << xn << "," << Ny << "," << y0 << "," << yn << endl;
                for (int yi = 0; yi <= Ny+7; yi++)
                {
                    //out2 << us.at(0, 0, yi);
                    for (int xi = 0; xi <= Nx+7; xi++)
                    {
                        long double Bx = (us.at(5, xi+1, yi)+us.at(5, xi, yi))/2.0; long double By = (us.at(6, xi, yi+1)+us.at(6, xi, yi))/2.0; long double Bz = us.at(7, xi, yi);
                        long double B2 = Bx*Bx+By*By+Bz*Bz;
                        long double v2 = powl(us.at(1, xi, yi)/us.at(0, xi, yi), 2)+powl(us.at(2, xi, yi)/us.at(0, xi, yi), 2)+powl(us.at(3, xi, yi)/us.at(0, xi, yi), 2);
                        long double PB = B2/2.0;
                        long double P = (us.at(4, xi, yi) - us.at(0, xi, yi)*v2/2 - B2/2) * (gamma-1);
                        long double DivB = (us.at(5, xi+1, yi)-us.at(5, xi, yi))/dx + (us.at(6, xi, yi+1)-us.at(6, xi, yi))/dy;
                        //out2 << sqrtl(powl(us.at(1, xi, yi)/us.at(0, xi, yi),2)+powl(us.at(2, xi, yi)/us.at(0, xi, yi), 2)) << ",";
                        //out2 << PB << ",";
                        out2 << us.at(0, xi, yi) << ","; //us.at(2, xi, yi)/us.at(0, xi, yi)
                    }
                }
                out2 << endl; out2.close();
            }

            

            string dataout = "out/data_" + to_string(t) + ".csv";
            
            

            if (Nx > 32 && Ny <= 32) {
                ofstream out(dataout);
                out << Nx << "," << x0 << "," << xn << endl;
                int yii = (int)round(Ny/2.0+4.0);   //31(Vortex)
                //int yii = 32;
                out << us.at(0, 0, yii);// us.at(2, 0, yii)/us.at(0, 0, yii)
                for (int xi = 1; xi < Nx+8; xi++)
                {
                    //cout << xi << "; " << yii << endl;
                    out << "," << us.at(0, xi, yii);// us.at(2, xi, yii)/us.at(0, xi, yii)
                }
                out << endl; out.close();
            }
            if (Ny > 32 && Nx <= 32) {
                ofstream out(dataout);
                out << Ny << "," << y0 << "," << yn << endl;
                int xii = (int)round(Nx/2.0+4.0);   //31(Vortex)

                long double B2 = us.at(5, xii, 0)*us.at(5, xii, 0) + us.at(6, xii, 0)*us.at(6, xii, 0) + us.at(7, xii, 0)*us.at(7, xii, 0);
                long double v2 = powl(us.at(1, xii, 0)/us.at(0, xii, 0), 2)+powl(us.at(2, xii, 0)/us.at(0, xii, 0), 2)+powl(us.at(3, xii, 0)/us.at(0, xii, 0), 2);
                long double PB = B2/2.0;
                long double P = (us.at(4, xii, 0) - us.at(0, xii, 0)*v2/2 - B2/2) * (gamma-1);

                out << P;
                //out << us.at(2, xii, 0)/us.at(0, xii, 0);// us.at(1, xii, 0)/us.at(0, xii, 0)
                for (int yi = 1; yi < Ny+8; yi++)
                {
                    long double B2 = us.at(5, xii, yi)*us.at(5, xii, yi) + us.at(6, xii, yi)*us.at(6, xii, yi) + us.at(7, xii, yi)*us.at(7, xii, yi);
                    long double v2 = powl(us.at(1, xii, yi)/us.at(0, xii, yi), 2)+powl(us.at(2, xii, yi)/us.at(0, xii, yi), 2)+powl(us.at(3, xii, yi)/us.at(0, xii, yi), 2);
                    long double PB = B2/2.0;
                    long double P = (us.at(4, xii, yi) - us.at(0, xii, yi)*v2/2 - B2/2) * (gamma-1);
                    out << "," << us.at(0, xii, yi);

                    //out << "," << us.at(2, xii, yi)/us.at(0, xii, yi);// us.at(1, xii, yi)/us.at(0, xii, yi)
                }
                out << endl; out.close();
            }
            
        }

        cout << endl; cout << endl;
        globstop = high_resolution_clock::now();
        ts[5] += duration<double>(globstop-globstart).count();


        retry:
        int col = 3;
        


        //return 3003;
        //break;
    }
    //dt for 12h (tn=3000) at 128x64: 0.02 --> 216'000 time-steps (dt if v+c=1 and dx=60/128: 0.23)
    stop:

    //0 = reconstruction,   1 = fluxes,   2 = update,   3 = boundary conditions,   4 = next dt,   5 = total
    cout << endl;
    cout << "Average time for each step: reconstruction = " << ts[0]/s << ", fluxes = " << ts[1]/s << ", update = " << ts[2]/s << ", boundary conditions = " << ts[3]/s << ", next dt = " << ts[4]/s << endl;
    cout << "Average time per timestep: " << ts[5]/s << endl;
    cout << "Total time (estimated): " << ts[5] << endl;

    int iueuf = 12;
    //Save to file
    //cout << "Bx: ";

    //ofstream out("data.csv");

    //if (Ny > 50 && Nx < 50) {
        //out << us.at(0, (int)round(Nx/2.0+4.0), 0);
        ////out << iy_us.at(0, (int)round(Nx/2.0+4.0), 0);
        //for (int yi = 1; yi < Ny+8; yi++)
        //{
        //    out << "," << us.at(0, (int)round(Nx/2.0+4.0), yi);
        //    //out << "," << iy_us.at(0, (int)round(Nx/2.0+4.0), yi+2);
        //}
    //}
    //if (Nx > 50 && Ny < 50) {
        
    //}

    

    //1d along diagonal:
    //ofstream out("data.csv");
    //out << us.at(0, 0, Nx+7); //0,Nx+5  1,Nx+4   2,Nx+3   3,Nx+2   4,Nx+1   ...
    //for (int xi = 0; xi < Nx+7; xi++)
    //{
    //    out << "," << us.at(0, xi, Nx+7-xi);
    //}

    //1d along other diagonal (y=1-x):
    //ofstream out("data.csv");
    //out << us.at(0, 0, 0); //0,Nx+5  1,Nx+4   2,Nx+3   3,Nx+2   4,Nx+1   ...
    //for (int yi = 1; yi < Ny+8; yi++)
    //{
    //    out << "," << us.at(0, yi, yi);
    //    //cout << Ny+7-yi << ", " << yi << endl;
    //}

}