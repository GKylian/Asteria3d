#pragma once
#include "Defs/Definitions.h"
#include "Misc/readParams.h"

#include <string>
#include <map>

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
#include "Bounds/primBounds.h"

#define PNAME "rtinstability"

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

#ifdef Iso
  #undef Iso
#endif // Iso

#define gx 0.0
#define gy 0.1
#define gz 0.0

ld Phi(ld x, ld y, ld z) {
	return 0;
}

ld Phii(Arrays *u, int i, int j, int k) {
	ld x[3] = { 0 }; 
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

bool Problem(Arrays *u, std::vector<Export> *out, map2d *params, std::string fname);
bool DoInLoop(Arrays *u);
bool ApplyBounds(Arrays *u);

void RTreflect_yin(Arrays *u);
void RTreflect_yout(Arrays *u);




bool Problem(Arrays *u, std::vector<Export> *out, map2d *params, std::string fname) {
	cout << "problem():: Loading the problem" << endl;

	/* Load map with parameters from the prob.params file (copy-pasted in the launch bash script) */
	readParams(fname, params);

	u->Nx = stod((*params)["domain"]["Nx"]); u->Ny = stod((*params)["domain"]["Ny"]); u->Nz = stod((*params)["domain"]["Nz"]);
	u->x0 = stod((*params)["domain"]["x0"]); u->y0 = stod((*params)["domain"]["y0"]); u->z0 = stod((*params)["domain"]["z0"]);
	u->xn = stod((*params)["domain"]["xn"]); u->yn = stod((*params)["domain"]["yn"]); u->zn = stod((*params)["domain"]["zn"]);

	u->boundaries[0] = strBC[(*params)["domain"]["BCxin"]]; u->boundaries[1] = strBC[(*params)["domain"]["BCxout"]];
	u->boundaries[2] = strBC[(*params)["domain"]["BCyin"]]; u->boundaries[3] = strBC[(*params)["domain"]["BCyout"]];
	u->boundaries[4] = strBC[(*params)["domain"]["BCzin"]]; u->boundaries[5] = strBC[(*params)["domain"]["BCzout"]];

	u->gamma = stod((*params)["problem"]["gamma"]); ld gamma = u->gamma;

	ld rhoUp = stod((*params)["problem"]["rhoUp"]);
	ld rhoDown = stod((*params)["problem"]["rhoDown"]);
	ld P0 = stod((*params)["problem"]["P0"]);
	

	/* Create the arrays */
	u->seth();
	u->initAll();
	cout << "problem():: Domain size: (" << u->Nx << ", " << u->Ny << ", " << u->Nz << ")." << endl;
	cout << "problem():: Steps: (" << u->dx << ", " << u->dy << ", " << u->dz << ")." << endl;

	/* Cycle through the domain (with one ghost cell) */
	for(int k = u->k_gl; k <= u->k_gr; k++)
	for(int j = u->j_gl; j <= u->j_gr; j++)
	for(int i = u->i_gl; i <= u->i_gr; i++)
	{
		/* Cell center coordinates */
		ld x = u->x0 + u->dx*(i-NGHOST);
		ld y = u->y0 + u->dy*(j-NGHOST);
		ld z = u->z0 + u->dz*(k-NGHOST);

		/* Cell interface coordinates */
		ld ix = x - 0.5*u->dx; ld iy = y - 0.5*u->dy; ld iz = z - 0.5*u->dz;

		/* Set HD variables as primitive variables */
		
		if (y > 0.0)
			u->uP(0, i, j, k) = rhoUp;
		else
			u->uP(0, i, j, k) = rhoDown;

		u->uP(4, i, j, k) = 1.0/gamma - gy*u->uP(0, i, j, k)*y;
		u->uP(1, i, j, k) = 0.0; u->uP(2, i, j, k) = 0.0; u->uP(3, i, j, k) = 0.0;
		
		if(fabsl(y) < 0.3)
			u->uP(2, i, j, k) = 0.0025*(1+cosl(4*M_PI*x))*(1+cosl(3*M_PI*y));

#ifdef MHD
		u->uC(5, i, j, k) = 0.0; u->uC(6, i, j, k) = 0.0; u->uC(7, i, j, k) = 0.0;
#endif // MHD


	}
	cout << "problem():: Initialized the variables" << endl;

		


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
	cout << "problem():: Transformed to conserved variables and computed cell centered magnetic fields" << endl;

	DoInLoop(u);

	return true;
}


bool DoInLoop(Arrays *u) {

#ifdef MHD
	//Check the divergence
	ld maxDiv = 0.0; int mi = 0, mj = 0, mk = 0;
	ld div = 0.0;
	for(int k = u->k_cl; k <= u->k_cr; k++)
	for(int j = u->j_cl; j <= u->j_cr; j++)
	for(int i = u->i_cl; i <= u->i_cr; i++)
	{
		div = getDiv(u, i, j, k);
		if (fabsl(div) > fabsl(maxDiv)) {
			maxDiv = div;
			mi = i; mj = j; mk = k;
		}
		
	}
	std::cout << "Maximum divergence in DoInLoop: " << maxDiv << " at (" << mi << ", " << mj << ", " << mk << ")." << std::endl;
	if (u->dim == 1 && maxDiv != 0.0) {
		cout << "\tThere is one dimension, but the divergence is not zero !" << endl;
		return false;
	}
	

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
			if(std::isnan(u->uP(n, i, j, k))){ std::cout << "check:: uP[" << n << "](" << i << ", " << j << ", " << k << ") is NaN." << std::endl; return false; }
			if(std::isnan(u->uC(n, i, j, k))){ std::cout << "check:: uC[" << n << "](" << i << ", " << j << ", " << k << ") is NaN." << std::endl; return false; }
        }
		if(u->uP(0, i, j, k) <= 0) { std::cout << "check:: rho in uP at (" << i << ", " << j << ", " << k << ") is null or negative: " << u->uP(0, i, j, k) << std::endl; return false; }
		if(u->uC(0, i, j, k) <= 0) { std::cout << "check:: rho in uC at (" << i << ", " << j << ", " << k << ") is null or negative: " << u->uC(0, i, j, k) << std::endl; return false; }

		if(u->uC(4, i, j, k) <= 0) { std::cout << "check:: E at (" << i << ", " << j << ", " << k << ") is null or negative: " << u->uC(4, i, j, k) << std::endl; return false; }
		if(u->uP(4, i, j, k) <= 0) { std::cout << "check:: P at (" << i << ", " << j << ", " << k << ") is null or negative: " << u->uP(4, i, j, k) << std::endl; return false; }
	}

	//TODO: Remove this when I figured out what's causing the difference...
	/* Check that the cell-centered average of B in uP is correct */
	/*for(int k = u->k_cl; k <= u->k_cr; k++)
	for(int j = u->j_cl; j <= u->j_cr; j++)
	for(int i = u->i_cl; i <= u->i_cr; i++)
	{
		if (0.5*(u->uC(5, i, j, k)+u->uC(5, i+1, j, k)) != u->uP(5, i, j, k)) {
			std::cout << "\tcheck:: The cell-centered value of Bx in uP is not correct: " << u->uP(5, i, j, k) << " but should be " << 0.5*(u->uC(5, i, j, k)+u->uC(5, i+1, j, k)) << std::endl;
			return false;
		}
		if (0.5*(u->uC(6, i, j, k)+u->uC(6, i, j+1, k)) != u->uP(6, i, j, k)) {
			std::cout << "\tcheck:: The cell-centered value of By in uP is not correct: " << u->uP(6, i, j, k) << " but should be " << 0.5*(u->uC(6, i, j, k)+u->uC(6, i, j+1, k)) << std::endl;
			return false;
		}
	}*/
	return true;
}


bool ApplyBounds(Arrays *u) {
	if (u->Nx > 1 && u->boundaries[0] != bounds::USER) bounds_xin(u);
	if (u->Nx > 1 && u->boundaries[1] != bounds::USER) bounds_xout(u);

	if (u->Ny > 1 && u->boundaries[2] != bounds::USER) bounds_yin(u);
	else if (u->Ny > 1 && u->boundaries[2] == bounds::USER) RTreflect_yin(u);
	if (u->Ny > 1 && u->boundaries[3] != bounds::USER) bounds_yout(u);
	else if (u->Ny > 1 && u->boundaries[3] == bounds::USER) RTreflect_yout(u);

	if (u->Nz > 1 && u->boundaries[4] != bounds::USER) bounds_zin(u);
	if (u->Nz > 1 && u->boundaries[5] != bounds::USER) bounds_zout(u);

	/* For the whole grid (including ghost cells), compute the cell-centered magnetic fields (if MHD) and the conserved variables */
	for(int k = u->k_dl; k <= u->k_dr; k++)
	for(int j = u->j_dl; j <= u->j_dr; j++)
	for(int i = u->i_dl; i <= u->i_dr; i++)
	{
#ifdef MHD
		BCenter(u, i, j, k);
#endif // MHD
		toConserved(u, i, j, k, u->gamma);

	}
	cout << "prob.h:::ApplyBounds:: Applied the boundary conditions" << endl;

	return true;
}






/* ---------- Boundary functions (if any) ---------- */

void RTreflect_yin(Arrays *u) {
	/* Loop through the transverse dimensions (x,z) */
	for(int k = u->k_dl; k <= u->k_dr; k++)
	for(int i = u->i_dl; i <= u->i_dr; i++)
	{
		/* Loop through the ghost cells */
		for (int j = 0; j < NGHOST; j++) {
			u->uP(0, i, u->j_gl-j, k) = u->uP(0, i, u->j_cl+j, k);
			u->uP(1, i, u->j_gl-j, k) = u->uP(1, i, u->j_cl+j, k); u->uP(3, i, u->j_gl-j, k) = u->uP(3, i, u->j_cl+j, k);
			u->uP(2, i, u->j_gl-j, k) = -u->uP(2, i, u->j_cl+j, k); // we reflect so '-'

			u->uP(4, i, u->j_gl-j, k) = u->uP(4, i, u->j_cl+j, k);
			u->uP(4, i, u->j_gl-j, k) += u->uP(0, i, u->j_gl-j, k)*gy*(2*j+1)*u->dy; // += rho*g*(2j+1)dy
		}
	}
	
}

void RTreflect_yout(Arrays *u) {
	for(int k = u->k_dl; k <= u->k_dr; k++)
	for(int i = u->i_dl; i <= u->i_dr; i++)
	{
		/* Loop through the ghost cells */
		for (int j = 0; j < NGHOST; j++) {
			u->uP(0, i, u->j_gr+j, k) = u->uP(0, i, u->j_cr-j, k);
			u->uP(1, i, u->j_gr+j, k) = u->uP(1, i, u->j_cr-j, k); u->uP(3, i, u->j_gr+j, k) = u->uP(3, i, u->j_cr-j, k);
			u->uP(2, i, u->j_gr+j, k) = -u->uP(2, i, u->j_cr-j, k); // we reflect so '-'

			u->uP(4, i, u->j_gr+j, k) = u->uP(4, i, u->j_cr-j, k);
			u->uP(4, i, u->j_gr+j, k) -= u->uP(0, i, u->j_gr+j, k)*gy*(2*j+1)*u->dy; // -= rho*g*(2j+1)dy
		}
	}
}