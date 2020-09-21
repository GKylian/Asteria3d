// Asteria3d.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include "Defs/Definitions.h"
#include "Defs/Arrays.h"
#include "prob.h"
#include "Export/vtk.h"
#include "Export/csv.h"

#ifdef PRIMBOUNDS
  #include "primBounds.h"
#else
  #include "consBounds.h"
#endif // PRIMBOUNDS

#include "Reconstruction/plm.h"
#include "Reconstruction/pcm.h"

#ifdef HLLE
  #include "Riemann/hlle.h"
#endif // HLLE
#ifdef HLLC
  #include "Riemann/hllc.h"
#endif // HLLC
#ifdef HLLD
  #include "Riemann/hlld.h"
#endif // HLLD


Arrays u;
Arrays _u;

using namespace std;

bool init(); bool Loop();
long double getdt(long double CFL);
long double cfl = 0.3;
vector<Export> out;


int main()
{

    cout << "Hello World!\n";


    /*  Initialize the arrays and set initial conditions  */
    if (!init()) {
        cout << "Asteria3d::main:: The initialization of the simulation failed !" << endl;
        return -1;
    }


    /* Big loop */
    if (!Loop()) {
        cout << "Asteria3d::main:: The computation failed !" << endl;
    }
}

bool init() {
    /* Coordinates:
    <-->
    ---> Cartesian coordinates: (x,y,z) = (x,y,z)
    ---> Cylindrical coordinates: (x,y,z) = (r,phi,z)
    <-->
    */
    //For cartesian coordinates: (x,y,z) = (x,y,z)




    /* Initialize the arrays with the problem chosen */
    if (!Problem(&u, &out)) {
        cout << "Asteria3D::init:: The initialization of the values from Problem() failed !" << endl; return false;
    }; cout << "Initialized the problem-specific intial values with Problem()." << endl;


    /* Apply the boundary conditions */
    if (!allBounds(&u)) {
        cout << "Asteria3d::init:: Could not apply the boundary conditions !" << endl; return false;
    } 
    _u = u;
    cout << "Applied the boundary conditions to the inital values." << endl;

    /* Temporary: exports the initial data */
    /*exports exp[2] = { exports::PRIM, exports::DIVB };
    vtk_3d(&u, false, exp, 2);*/
    vtk_3d(&u, out);
    csv(&u, out);

    //DoInLoop(&u); 
    if (!check(&u)) {
        cout << "ERROR:::INIT:: Call to check() after initialization returned false -> aborting the simulation..." << endl;
        return false;
    }

    

    /* Compute the initial dt */
    //cout << test(0.23) << endl;

    long double dt0 = getdt(cfl);
    cout << "INIT:: Initial dt: " << dt0 << endl;
    u.dt = dt0;

    
    
    return true;
}


bool Loop() {
    /*
    ---> To do in loop:
    --->   - Reconstruct values at interfaces from cell averages (don't forget source terms)
    --->   - Compute fluxes at interfaces
    --->   - If MHD: compute EMFs at cell corners and use them to update Bx, By and Bz
    --->   - Update rho to n+1/2
    --->   - Add source terms
    --->   - Update all the variables to n+1
    --->   - Apply boundary conditions
    --->   - Compute a new dt
    --->   - When necessary, export data (.csv, .vtk, ...)
    --->   - Update t
    */
    long double t0 = 0.0; long double tn = 0.5;
    bool nx = u.Nx>1; bool ny = u.Ny>1; bool nz = u.Nz>1; long double hrho = 0.0;

    for (u.t = t0; u.t < tn; u.t=u.t+u.dt) {

        //if (u.s == 1) break;
        u.s = u.s+1;
        cout << endl << endl;
        std::cout << "LOOP:: Computing at time " << u.t << " with time-step " << u.dt << endl;

        /* ----- 1. Reconstruct interfaces from cell averages ----- */
        cout << "LOOP:: Starting the reconstruction..." << endl;
        for(int i = u.i_dl+1*nx; i <= u.i_dr-1*nx; i++)
        for(int j = u.j_dl+1*ny; j <= u.j_dr-1*ny; j++)
        for(int k = u.k_dl+1*nz; k <= u.k_dr-1*nz; k++)
        {
            
            if (!PLM(&u, i, j, k)) {
                cout << "ERROR:::LOOP:: Reconstruction failed !" << endl;
                cout << "---> Cell (" << i << ", " << j << ", " << k << ")" << endl;
                return false;
            }
        }
        cout << "LOOP:: Reconstructed interfaces from cell averages." << endl;


        /* ----- 2. Compute fluxes from newly computed interface values ----- */
        cout << "LOOP:: Starting to compute the fluxes..." << endl;
        for(int i = u.i_dl+2*nx; i <= u.i_dr-1*nx; i++)
        for(int j = u.j_dl+2*ny; j <= u.j_dr-1*ny; j++)
        for(int k = u.k_dl+2*nz; k <= u.k_dr-1*nz; k++)
        {
            if (!getFluxes(&u, i, j, k)) {
                cout << "ERROR:::LOOP:: getFluxes failed !" << endl;
                cout << "---> Cell (" << i << ", " << j << ", " << k << ")" << endl;
                return false;
            }
            for (int n = 0; n < NVAL; n++) {
                if (isnan(u.Fx(n, i, j, k))) {
                    cout << "ERROR:::LOOP:: Fx[" << n << "] is NaN !" << endl; return false;
                }
                if (isnan(u.Fy(n, i, j, k))) {
                    cout << "ERROR:::LOOP:: Fy[" << n << "] is NaN !" << endl; return false;
                }
                if (isnan(u.Fz(n, i, j, k))) {
                    cout << "ERROR:::LOOP:: Fz[" << n << "] is NaN !" << endl; return false;
                }
                    
            }
        }
        cout << "LOOP:: Computed fluxes." << endl;
        

#ifdef CT
        /* ----- 3. Compute EM fields at cell corners and update B ----- */
#endif // CT

        /* ----- 4. Update to n+1 ----- */
        for(int i = u.i_dl+3*nx; i <= u.i_dr-2*nx; i++)
        for(int j = u.j_dl+3*ny; j <= u.j_dr-2*ny; j++)
        for(int k = u.k_dl+3*nz; k <= u.k_dr-2*nz; k++)
        {
            /* --- 4.1) Compute rho at t(n+1/2) --- */
            hrho = u.uC(0, i, j, k);
            if (u.Nx > 1) hrho += -u.dt/(2.0*u.dx) * (u.Fx(0, i+1, j, k)-u.Fx(0, i, j, k));
            if (u.Ny > 1) hrho += -u.dt/(2.0*u.dy) * (u.Fy(0, i, j+1, k)-u.Fy(0, i, j, k));
            if (u.Nz > 1) hrho += -u.dt/(2.0*u.dz) * (u.Fz(0, i, j, k+1)-u.Fz(0, i, j, k));


            /* --- 4.2) Add the gravitational source terms to _u --- */
            //TODO: Add them
            /* Gravity with acceleration \vec{g} */
            if (gx != 0.0 || gy != 0.0 || gz != 0.0) {
                _u.uC(1, i, j, k) += u.dt*hrho*gx;
                _u.uC(2, i, j, k) += u.dt*hrho*gy;
                _u.uC(3, i, j, k) += u.dt*hrho*gz;
#ifndef ISO
                if (u.Nx > 1) _u.uC(4, i, j, k) += u.dt*gx*(u.Fx(0, i+1, j, k)+u.Fx(0, i, j, k))/2.0;
                if (u.Ny > 1) _u.uC(4, i, j, k) += u.dt*gy*(u.Fy(0, i, j+1, k)+u.Fy(0, i, j, k))/2.0;
                if (u.Nz > 1) _u.uC(4, i, j, k) += u.dt*gz*(u.Fz(0, i, j, k+1)+u.Fz(0, i, j, k))/2.0;
#endif // !ISO

            }
            /* Gravity with static potential defined in prob.h */
            else if (Phii(&u, i, j, k) != 0.0){
                long double x[3] = { 0 }; u.pos(x, i, j, k); long double Phic = Phii(&u, i, j, k); long double Phir, Phil = 0.0;
                if (u.Nx>1) {
                    Phir = Phi(x[0]+0.5*u.dx, x[1], x[2]); Phil = Phi(x[0]-0.5*u.dx, x[1], x[2]);
                    _u.uC(1, i, j, k) -= u.dt/u.dx*hrho*(Phir-Phil);
#ifndef ISO
                    _u.uC(4, i, j, k) -= u.dt/u.dx*(u.Fx(0, i, j, k)*(Phic-Phil) + u.Fx(0, i+1, j, k)*(Phir-Phic));
#endif // !ISO
                }
                if (u.Ny>1) {
                    Phir = Phi(x[0], x[1]+0.5*u.dy, x[2]); Phil = Phi(x[0], x[1]-0.5*u.dy, x[2]);
                    _u.uC(2, i, j, k) -= u.dt/u.dy*hrho*(Phir-Phil);
#ifndef ISO
                    _u.uC(4, i, j, k) -= u.dt/u.dy*(u.Fy(0, i, j, k)*(Phic-Phil) + u.Fy(0, i, j+1, k)*(Phir-Phic));
#endif // !ISO
                }
                if (u.Nz>1) {
                    Phir = Phi(x[0], x[1], x[2]+0.5*u.dz); Phil = Phi(x[0], x[1], x[2]-0.5*u.dz);
                    _u.uC(3, i, j, k) -= u.dt/u.dz*hrho*(Phir-Phil);
#ifndef ISO
                    _u.uC(4, i, j, k) -= u.dt/u.dz*(u.Fz(0, i, j, k)*(Phic-Phil) + u.Fz(0, i, j, k+1)*(Phir-Phic));
#endif // !ISO
                }

            }




            /* --- 4.3) Update u to n+1 in _u using the fluxes previously computed --- */
            int Bstop = 0;
#ifdef CT
            if (MHD) Bstop = 3;
#endif // CT

            for (int n = 0; n < NWAVE-Bstop; n++) {
                long double dFxdx = 0.0; if (u.Nx>1) dFxdx = u.dt/u.dx * (u.Fx(n, i+1, j, k)-u.Fx(n, i, j, k));
                long double dFydy = 0.0; if (u.Ny>1) dFydy = u.dt/u.dy * (u.Fy(n, i, j+1, k)-u.Fy(n, i, j, k));
                long double dFzdz = 0.0; if (u.Nz>1) dFzdz = u.dt/u.dz * (u.Fz(n, i, j, k+1)-u.Fz(n, i, j, k));
                _u.uC(n, i, j, k) += -dFxdx - dFydy - dFzdz;

                /*if (null(dFxdx) && null(dFydy) && null(dFydy)) {
                    cout << "_u.uC[" << n << "] doesn't change !" << endl;
                    cout << "---> dFxdx: " << dFxdx << ", dFydy: " << dFydy << ", dFzdz: " << dFzdz << endl;
                    cout << "---> F(i+1) = " << u.Fx(n, i+1, j, k) << ", f(i) = " << u.Fx(n, i, j, k) << endl;
                }*/

                if (isnan(_u.uC(n, i, j, k))) {
                    cout << "ERROR:::LOOP:: After evolution to n+1, uC[" << n << "] is NaN !" << endl;
                    cout << "---> Cell (" << i << ", " << j << ", " << k << ")" << endl;
                    cout << "---> dFndn = " << dFxdx << ", " << dFydy << ", " << dFzdz << endl;
                    cout << "---> dFn = " << (u.Fx(n, i+1, j, k)-u.Fx(n, i, j, k)) << ", " << (u.Fy(n, i, j+1, k)-u.Fy(n, i, j, k)) << ", " << (u.Fz(n, i, j, k+1)-u.Fz(n, i, j, k)) << endl;
                    cout << "---> Fn(n+1) = " << u.Fx(n, i+1, j, k) << ", " << u.Fy(n, i, j+1, k) << ", " << u.Fz(n, i, j, k+1) << endl;
                    return false;
                }
            }
        }
        for(int i = u.i_dl; i <= u.i_dr; i++)
        for(int j = u.j_dl; j <= u.j_dr; j++)
        for(int k = u.k_dl; k <= u.k_dr; k++)
        {
            for (int n = 0; n < NVAL; n++) {
                u.uC(n, i, j, k) = _u.uC(n, i, j, k);
                u.uP(n, i, j, k) = _u.uP(n, i, j, k);
            }
        }
        cout << "LOOP:: Updated to n+1." << endl;
        


        /* ----- 5. Get primitive variables and B center ----- */
        for(int i = u.i_dl; i <= u.i_dr; i++)
        for(int j = u.j_dl; j <= u.j_dr; j++)
        for(int k = u.k_dl; k <= u.k_dr; k++)
        {
#ifdef MHD
            BCenter(&u, i, j, k);
#endif

            toPrimitive(&u, i, j, k, gamma);

            if (u.uC(0, i, j, k) <= 0.0 || null(u.uC(0, i, j, k))) {
                cout << "ERROR:::LOOP:: rho in uC is either NaN or null/negative: " << u.uC(0, i, j, k) << endl;
                cout << "---> Cell (" << i << ", " << j << ", " << k << ")" << endl;
                return false;
            }
            if (u.uP(0, i, j, k) <= 0.0 || null(u.uP(0, i, j, k))) {
                cout << "ERROR:::LOOP:: rho in uP is either NaN or null/negative: " << u.uP(0, i, j, k) << endl;
                cout << "---> Cell (" << i << ", " << j << ", " << k << ")" << endl;
                return false;
            }

            if (u.uC(4, i, j, k) <= 0.0 || null(u.uC(4, i, j, k))) {
                cout << "ERROR:::LOOP:: E in uC is either NaN or null/negative: " << u.uC(4, i, j, k) << endl;
                cout << "---> Cell (" << i << ", " << j << ", " << k << ")" << endl;
                return false;
            }
            if (u.uP(4, i, j, k) <= 0.0 || null(u.uP(4, i, j, k))) {
                cout << "ERROR:::LOOP:: P in uP is either NaN or null/negative: " << u.uP(4, i, j, k) << endl;
                cout << "---> Cell (" << i << ", " << j << ", " << k << ")" << endl;
                return false;
            }
        }



        /* ----- 6. Apply boundary conditions ----- */
        allBounds(&u); //TODO: Transfer the handling of the BCs to the problem header.


        /* ----- 7. Call DoInLoop in the problem header ----- */
        DoInLoop(&u);


        /* ----- 8. Re-apply boundary conditions in case we changed values in DoInLoop ----- */
        allBounds(&u); //TODO: Transfer the handling of the BCs to the problem header.
        //Also check max/min values of uP and uC
        long double maxC[NVAL] = { -1000 }; long double minC[NVAL] = { 1000 };
        long double maxP[NVAL] = { -1000 }; long double minP[NVAL] = { 1000 };
        for(int i = u.i_dl; i <= u.i_dr; i++)
        for(int j = u.j_dl; j <= u.j_dr; j++)
        for(int k = u.k_dl; k <= u.k_dr; k++)
        {
            for (int n = 0; n < NVAL; n++) {
                _u.uC(n, i, j, k) = u.uC(n, i, j, k);
                _u.uP(n, i, j, k) = u.uP(n, i, j, k);

                maxC[n] = fmaxl(maxC[n], u.uC(n, i, j, k));
                minC[n] = fminl(minC[n], u.uC(n, i, j, k));
                maxP[n] = fmaxl(maxP[n], u.uP(n, i, j, k));
                minP[n] = fminl(minP[n], u.uP(n, i, j, k));
            }
        }
        cout << "LOOP:: Maximum values of uC: [" << maxC[0]; for (int n = 1; n < NVAL; n++) cout << ", " << maxC[n]; cout << "]." << endl;
        cout << "LOOP:: Minimum values of uC: [" << minC[0]; for (int n = 1; n < NVAL; n++) cout << ", " << minC[n]; cout << "]." << endl;
        cout << "LOOP:: Maximum values of uP: [" << maxP[0]; for (int n = 1; n < NVAL; n++) cout << ", " << maxP[n]; cout << "]." << endl;
        cout << "LOOP:: Minimum values of uP: [" << minP[0]; for (int n = 1; n < NVAL; n++) cout << ", " << minP[n]; cout << "]." << endl;

        /* ----- 9. Call the Check() function in the problem file ----- */
        if (!check(&u)) {
            cout << "ERROR:::LOOP:: Call to check() returned false -> stoping the loop..." << endl;
            return false;
        }

        /* ----- 10. Compute a new dt ----- */
        u.dt = getdt(cfl);

        
        vtk_3d(&u, out);
        csv(&u, out);
    }


    return true;
}



long double getdt(long double CFL) {
    long double dt = 10000.0;
    
    long double cx[WV] = { 0 }; long double cy[WV] = { 0 }; long double cz[WV] = { 0 };
    long double dts[3] = { 0 };
    long double Cx = 0.0, Cy = 0.0, Cz = 0.0;

    // Loop through the whole domain (including ghost cells)
    for(int i = u.i_dl; i <= u.i_dr; i++)
    for(int j = u.j_dl; j <= u.j_dr; j++)
    for(int k = u.k_dl; k <= u.k_dr; k++)
    {
        
        getWavespeeds(&u, i, j, k, gamma, cx, cy, cz); 
        Cx = (long double)*std::max_element(cx, cx+WV); Cy = (long double)*std::max_element(cy, cy+WV); Cz = (long double)*std::max_element(cz, cz+WV);
        

        dts[0] = u.dx/(u.uP(1, i, j, k)+Cx); if (u.Nx<=1) dts[0] = 10000.0; if (isnan(dts[0])) dts[0] = 10000.0;
        dts[1] = u.dy/(u.uP(2, i, j, k)+Cy); if (u.Ny<=1) dts[1] = 10000.0; if (isnan(dts[1])) dts[1] = 10000.0;
        dts[2] = u.dz/(u.uP(3, i, j, k)+Cz); if (u.Nz<=1) dts[2] = 10000.0; if (isnan(dts[2])) dts[2] = 10000.0;
        long double thdt = CFL * (*std::min_element(dts, dts+3));

        if (thdt <= 0) {
            thdt = 100.0;
        }

        dt = fminl(dt, thdt);

    }
    cout << "Computed dt: " << dt << endl;

    return dt;

}






