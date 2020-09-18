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

bool Problem(Arrays *u);
bool DoInLoop(Arrays *u);


bool Problem(Arrays *u) {

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
		/* Compute cell-centered magnetic fields (in uP) from interfaces (in uC) */
		BCenter(u, i, j, k);
		
		/* Compute primitive variables */
		toConserved(u, i, j, k, gamma);


	}

	long double maxDiv = 0.0; 
	for(int k = u->k_cl; k <= u->k_cr; k++)
	for(int j = u->j_cl; j <= u->j_cr; j++)
	for(int i = u->i_cl; i <= u->i_cr; i++)
	{
		maxDiv = fmaxl(fabsl(maxDiv), fabsl(getDiv(u, i, j, k)));
		
	}
	std::cout << "Maximum divergence after problem initialization: " << maxDiv << std::endl;

	return true;
}


bool DoInLoop(Arrays *u) {

	//Check the divergence
	long double maxDiv = 0.0; 
	for(int k = u->k_cl; k <= u->k_cr; k++)
	for(int j = u->j_cl; j <= u->j_cr; j++)
	for(int i = u->i_cl; i <= u->i_cr; i++)
	{
		maxDiv = fmaxl(fabsl(maxDiv), fabsl(getDiv(u, i, j, k)));
		if (u->uP(0, i, j, k) != 25.0/(36*M_PI) || u->uC(0, i, j, k) != 25.0/(36*M_PI)) {
			std::cout << "Density: " << u->uP(0, i, j, k) << ", or C: " << u->uC(0, i, j, k) << std::endl;
		}
		
	}
	std::cout << "Maximum divergence in DoInLoop: " << maxDiv << std::endl;

	return true;
}