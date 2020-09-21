#pragma once
#include "Defs/Definitions.h"

#ifdef MHD
  #ifndef ISO
    #include "Eq/AdiabMHD.h"
  #else
    #include "Eq/IsoMHD.h"
  #endif // !ISO
#else
  #ifndef ISO
    #include "Eq/AdiabHD.h" 
  #endif // !ISO
#endif // MHD


/*
<->
---> Problem: ----- Orszag-Tang Vortex (MHD and HD versions)
---> Source: ------ Orszag & Tang, 1998
---> EOS: --------- Adiabatic
---> Gamma: ------- 5/3
---> Gas: --------- MHD or HD
---> Domain: ------ x in (0;1), y in (0;1)
---> Resolution: -- (Nx, Ny) = (N, N)
---> CFL: --------- 0.4
<->
<->
---> Parameters:
--->  - B0: The maximum intensity of the magnetic field
*/

#define gamma 1.6666666667
#ifdef Iso
  #undef Iso
#endif // Iso
#ifdef MHD
  #define B0 0.282094792
#endif // MHD

#define gx 0.0
#define gy 0.0
#define gz 0.0

long double Phi(long double x, long double y, long double z) {
	return 0;
}

long double Phii(Arrays *u, int i, int j, int k) {
	long double x[3] = { 0 }; 
	u->pos(x, i, j, k);
	return Phi(x[0], x[1], x[2]);
}

/* 
<-------------------- Functions -------------------->
<->
---> bool Problem(Arrays *u):
--->   - Initialized all the values in the domain (w/o ghost cells) in the arrays u
--->   - Param::Arrays *u: pointer to the data on the whole grid
--->   - Returns bool: Has the initialization worked ?
<->
---> bool DoInLoop(Arrays u):
--->   - Do something (if anything) once the evolution n->n+1 is complete
--->   - Returns bool: Has the job worked ?

*/

bool Problem(Arrays *u, std::vector<Export> *out);
bool DoInLoop(Arrays *u);


bool Problem(Arrays *u, std::vector<Export> *out) {

	/* Set array parameters */
	long double x0 = 0.0; long double xn = 1.0; int Nx = 64; long double dx = (xn-x0)/(Nx-1.0); if (Nx <= 1) dx = 1.0;
	long double y0 = 0.0; long double yn = 1.0; int Ny = 64; long double dy = (yn-y0)/(Ny-1.0); if (Ny <= 1) dy = 1.0;
	long double z0 = 0.0; long double zn = 1.0; int Nz = 1; long double dz = (zn-z0)/(Nz-1.0); if (Nz <= 1) dz = 1.0;

	u->boundaries[0] = bounds::PERIODIC; u->boundaries[1] = bounds::PERIODIC; /* x */
	u->boundaries[2] = bounds::PERIODIC; u->boundaries[3] = bounds::PERIODIC;

	/* Create the arrays */
	u->initAll(Nx, Ny, Nz); u->setRange(x0, xn, y0, yn, z0, zn); u->seth(dx, dy, dz);
	cout << "INIT:: Domain size: (" << u->Nx << ", " << u->Ny << ", " << u->Nz << ")." << endl;

	/* Cycle through the domain (with one ghost cell) */
	for(int k = u->k_gl; k <= u->k_gr; k++)
	for(int j = u->j_gl; j <= u->j_gr; j++)
	for(int i = u->i_gl; i <= u->i_gr; i++)
	{
		/* Cell center coordinates */
		long double x = u->x0 + u->dx*(i-NGHOST);
		long double y = u->y0 + u->dy*(j-NGHOST);
		//  long double z = u->z0 + u->dz*(k-NGHOST);  --> 2d, not needed

		/* Cell interface coordinates */
		long double ix = x - 0.5*u->dx; long double iy = y - 0.5*u->dy;

		/* Set HD variables as primitive variables */
		u->uP(0, i, j, k) = 25.0/(36*M_PI); u->uP(4, i, j, k) = 5.0/(12.0*M_PI);
		u->uP(1, i, j, k) = -sinl(2*M_PI*y); u->uP(2, i, j, k) = sinl(2*M_PI*x); u->uP(3, i, j, k) = 0.0;
		

#ifdef MHD
		/* Compute potential at cell corners */
		long double Az[2][2] = { 0 };
		for(int jj = 0; jj < 2; jj++)
		for(int ii = 0; ii < 2; ii++)
		{
			/* Coordinates at cell corner of (i+ii, j+jj) */
			long double xc = u->x0 + u->dx*(i+ii-4.5); long double yc = u->y0 + u->dy*(j+jj-4.5);
			Az[ii][jj] = B0 * ( cosl(4*M_PI*xc)/(4*M_PI) + cosl(2*M_PI*yc)/(2*M_PI) );
		}

		/* Compute Bx and By at interface from potentials at cell corners */
		u->uC(5, i, j, k) = (Az[0][1]-Az[0][0])/u->dy; u->uC(6, i, j, k) = -(Az[1][0]-Az[0][0])/u->dx;
		u->uC(7, i, j, k) = 0.0;

#endif // MHD
	}

	/* Finish setting up both arrays */
	for(int k = u->k_cl; k <= u->k_cr; k++)
	for(int j = u->j_cl; j <= u->j_cr; j++)
	for(int i = u->i_cl; i <= u->i_cr; i++)
	{
#ifdef MHD
		/* Compute cell-centered magnetic fields (in uP) from interfaces (in uC) */
		BCenter(u, i, j, k);
#endif // MHD

		
		/* Compute primitive variables */
		toConserved(u, i, j, k, gamma);


	}
#ifdef MHD
	long double maxDiv = 0.0;
	for (int k = u->k_cl; k <= u->k_cr; k++)
		for (int j = u->j_cl; j <= u->j_cr; j++)
			for (int i = u->i_cl; i <= u->i_cr; i++)
			{
				maxDiv = fmaxl(fabsl(maxDiv), fabsl(getDiv(u, i, j, k)));

			}
	std::cout << "Maximum divergence after problem initialization: " << maxDiv << std::endl;


#endif // MHD

	/* Set up the output */
	std::vector<exports> vtkvar{ exports::PRIM };
	Export vtkout{ 0, vtkvar, 0, 0, 0, 0, 0, 0, false, "prims" };
	out->push_back(vtkout);

	return true;
}


bool DoInLoop(Arrays *u) {

#ifdef MHD
	//Check the divergence
	long double maxDiv = 0.0; 
	for(int k = u->k_cl; k <= u->k_cr; k++)
	for(int j = u->j_cl; j <= u->j_cr; j++)
	for(int i = u->i_cl; i <= u->i_cr; i++)
	{
		maxDiv = fmaxl(fabsl(maxDiv), fabsl(getDiv(u, i, j, k)));
		
	}
	std::cout << "Maximum divergence in DoInLoop: " << maxDiv << std::endl;

	

#endif // MHD
	return true;
}

bool check(Arrays *u){
	/* Check that nothing NaN, that rho, P and E are positive */
	for(int k = u->k_dl; k <= u->k_dr; k++)
	for(int j = u->j_dl; j <= u->j_dr; j++)
	for(int i = u->i_dl; i <= u->i_dr; i++)
	{
		for(int n = 0; n < NVAL; n++){
			if (std::isnan(u->uP(n, i, j, k))) { std::cout << "check:: uP[" << n << "](" << i << ", " << j << ", " << k << ") is NaN." << std::endl; return false; }
			if (std::isnan(u->uC(n, i, j, k))) { std::cout << "check:: uC[" << n << "](" << i << ", " << j << ", " << k << ") is NaN." << std::endl; return false; }
        }
		if (u->uC(0, i, j, k) <= 0) { std::cout << "check:: rho at (" << i << ", " << j << ", " << k << ") is null or negative." << std::endl; return false; }
		if (u->uP(0, i, j, k) <= 0) { std::cout << "check:: rho at (" << i << ", " << j << ", " << k << ") is null or negative." << std::endl; return false; }

		if (u->uC(4, i, j, k) <= 0) { std::cout << "check:: E at (" << i << ", " << j << ", " << k << ") is null or negative." << std::endl; return false; }
		if (u->uP(4, i, j, k) <= 0) { std::cout << "check:: P at (" << i << ", " << j << ", " << k << ") is null or negative." << std::endl; return false; }
	}

	return true;
}