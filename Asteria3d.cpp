// Asteria3d.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include "Defs/Definitions.h"
#include "Defs/Arrays.h"
#include "prob.h"
#include "Export/vtk.h"

#ifdef PRIMBOUNDS
  #include "primBounds.h"
#else
  #include "consBounds.h"
#endif // PRIMBOUNDS

#include "Reconstruction/plm.h" 


Arrays u;
Arrays _u;

using namespace std;

bool init(); bool Loop();
long double getdt(long double CFL);
long double cfl = 0.3;


int main()
{
    std::cout << "Hello World!\n";


    /*  Initialize the arrays and set initial conditions  */
    if (!init()) {
        cout << "Asteria3d::main:: The initialization of the simulation failed !" << endl;
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



    /* Set array parameters */
    long double x0 = 0.0; long double xn = 1.0; int Nx = 32; long double dx = (xn-x0)/(Nx-1.0); if (Nx == 1) dx = 1.0;
    long double y0 = 0.0; long double yn = 1.0; int Ny = 32; long double dy = (yn-y0)/(Ny-1.0); if (Ny == 1) dy = 1.0;
    long double z0 = 0.0; long double zn = 1.0; int Nz = 1; long double dz = (zn-z0)/(Nz-1.0); if (Nz == 1) dz = 1.0;

    u.boundaries[0] = bounds::PERIODIC; u.boundaries[1] = bounds::PERIODIC; /* x */
    u.boundaries[2] = bounds::PERIODIC; u.boundaries[3] = bounds::PERIODIC;

    /* Create the arrays */
    u.initAll(Nx, Ny, Nz); u.setRange(x0, xn, y0, yn, z0, zn); u.seth(dx, dy, dz);
    cout << "INIT:: Domain size: (" << u.Nx << ", " << u.Ny << ", " << u.Nz << ")." << endl;

    /* Initialize the arrays with the problem choosen */
    if (!Problem(&u)) {
        cout << "Asteria3D::init:: The initialization of the values from Problem() failed !" << endl; return false;
    }; cout << "Initialized the problem-specific intial values with Problem()." << endl;
    
    /* Apply the boundary conditions */
    if (!allBounds(&u)) {
        cout << "Asteria3d::init:: Could not apply the boundary conditions !" << endl; return false;
    } 
    _u = u;
    cout << "Applied the boundary conditions to the inital values." << endl;
    

    /* Temporary: exports the initial data */
    exports exp[2] = { exports::PRIM, exports::DIVB };
    //vtk_3d(&u, false, exp, 2);

    //DoInLoop(&u);


    /* Compute the initial dt */
    u.dt = getdt(cfl);
    cout << "INIT:: Initial dt: " << u.dt << endl;

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
    long double t0 = 0.0; long double tn = 0.1;
    bool nx = u.Nx>1; bool ny = u.Ny>1; bool nz = u.Nz>1;

    for (u.t = t0; u.t < tn; u.t+=u.dt) {
        std::cout << "LOOP:: Computing at time " << u.t << "with time-step " << u.dt << endl;

        /* ----- 1. Reconstruct interfaces from cell averages ----- */
        for(int i = u.i_dl+1*nx; i <= u.i_dr-1*nx; i++)
        for(int j = u.j_dl+1*ny; j <= u.j_dr-1*ny; j++)
        for(int k = u.k_dl+1*nz; k <= u.k_dr-1*nz; k++)
        {
            cout << "Reconstructing at " << i << ", " << j << ", " << k << endl;
            if(!PLM(&u, i, j, k)) return false;
        }
    }


    return true;
}








long double getdt(long double CFL) {
    long double dt = 10000.0;

    long double cx[WV] = { 0 }; long double cy[WV] = { 0 }; long double cz[WV] = { 0 };
    long double dts[3] = { 0 };
    long double Cx = 0.0, Cy = 0.0, Cz = 0.0;

    /* Loop through the whole domain (including ghost cells) */
    for(int i = u.i_dl; i <= u.i_dr; i++)
    for(int j = u.j_dl; j <= u.j_dr; j++)
    for(int k = u.k_dl; k <= u.k_dr; k++)
    {
        getWavespeeds(&u, i, j, k, gamma, cx, cy, cz);
        Cx = *std::max_element(cx, cx+WV); Cy = *std::max_element(cy, cy+WV); Cz = *std::max_element(cz, cz+WV);

        dts[0] = u.dx/(u.uP(1, i, j, k)+Cx); if (u.Nx<=1) dts[0] = 10000.0; if (isnan(dts[0])) dts[0] = 10000.0;
        dts[1] = u.dy/(u.uP(2, i, j, k)+Cy); if (u.Ny<=1) dts[1] = 10000.0; if (isnan(dts[1])) dts[1] = 10000.0;
        dts[2] = u.dz/(u.uP(3, i, j, k)+Cz); if (u.Nz<=1) dts[2] = 10000.0; if (isnan(dts[2])) dts[2] = 10000.0;
        long double thdt = CFL * (*std::min_element(dts, dts+3));

        if (thdt <= 0) {
            thdt = 100.0;
        }

        dt = fminl(dt, thdt);

    }

    return dt;

}