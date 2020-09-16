// Asteria3d.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include "Definitions.h"
#include "Arrays.h"
#include "otvortex.h"
#include "vtk.h"

#ifdef PRIMBOUNDS
  #include "primBounds.h"
#else
  #include "consBounds.h"
#endif // PRIMBOUNDS

#include "Arrays.h"


Arrays u;
Arrays _u;

using namespace std;

bool init(); bool Loop();



int main()
{
    std::cout << "Hello World!\n";

    /*  Check that the definitions are good:  */
#ifdef HLLD
  #ifdef HLLC
    #error Both the HLLD and HLLC solvers are defined !
  #endif // HLLC
  #ifndef MHD
    #error Using the HLLD solver for hydrodynamics !
  #endif // !MHD

  #ifdef HLL
    #error Both the HLLD and HLL solvers are defined !
  #endif // HLL

#endif // HLLD


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
    vtk_3d(&u, false, exp, 2);

    DoInLoop(&u);
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

    return true;
}