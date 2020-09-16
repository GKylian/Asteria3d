#pragma once

#include <iostream>
#include <fstream>
using namespace std;

#include "Arrays.h"

bool vtk_add(Arrays *u, bool ghosts, exports exp)
{
	
	ofstream vtk("exports/out.vtk", ios::out | ios::app);
	if (exp == exports::Bvec || exp == exports::Vvec || exp == exports::Mvec) {
		/* Write header for  */
		vtk << "VECTORS " << exp_names[(int)exp] << " float" << endl;
	}
	else {
		/* Write header for scalar */
		vtk << "SCALARS " << exp_names[(int)exp] << " float" << endl;
		vtk << "LOOKUP_TABLE default" << endl;
	}

	/* Store the 'for' loops' limits to avoid computing them every iteration */
	int imax = u->Nx+(!ghosts)*NGHOST; int jmax = u->Ny+(!ghosts)*NGHOST; int kmax = u->Nz+(!ghosts)*NGHOST;





	/* We put the ifs outside the loops instead of inside a single big
	   one for better performance (1 'if' instead of Nx*Ny*Nz 'if's)   */

	/* -------------------- VECTOR EXPORT -------------------- */

	if (exp == exports::Bvec) {
		/* Write the three components each line */
		for(int k = (!ghosts)*NGHOST; k < kmax; k++)
		for(int j = (!ghosts)*NGHOST; j < jmax; j++)
		for(int i = (!ghosts)*NGHOST; i < imax; i++)
		{
			vtk << u->uP(5, i, j, k) << " " << u->uP(6, i, j, k) << " " << u->uP(7, i, j, k) << endl;
		}
	}
	if (exp == exports::Vvec) {
		/* Write the three components each line */
		for(int k = (!ghosts)*NGHOST; k < kmax; k++)
		for(int j = (!ghosts)*NGHOST; j < jmax; j++)
		for(int i = (!ghosts)*NGHOST; i < imax; i++)
		{
			vtk << u->uP(1, i, j, k) << " " << u->uP(2, i, j, k) << " " << u->uP(3, i, j, k) << endl;
		}
	}
	if (exp == exports::Mvec) {
		/* Write the three components each line */
		for(int k = (!ghosts)*NGHOST; k < kmax; k++)
		for(int j = (!ghosts)*NGHOST; j < jmax; j++)
		for(int i = (!ghosts)*NGHOST; i < imax; i++)
		{
			vtk << u->uC(1, i, j, k) << " " << u->uC(2, i, j, k) << " " << u->uC(3, i, j, k) << endl;
		}
	}
	




	/* -------------------- SCALAR EXPORT -------------------- */

	if (exp == exports::Rho) {
		/* Write the scalar each line */
		for(int k = (!ghosts)*NGHOST; k < kmax; k++)
		for(int j = (!ghosts)*NGHOST; j < jmax; j++)
		for(int i = (!ghosts)*NGHOST; i < imax; i++)
		{
			//cout << i << ", " << j << endl;
			//cout << u->uC(0, i, j, k) << endl;
			vtk << u->uC(0, i, j, k) << endl;
		}
	}
	if (exp == exports::P) {
		/* Write the scalar each line */
		for(int k = (!ghosts)*NGHOST; k < kmax; k++)
		for(int j = (!ghosts)*NGHOST; j < jmax; j++)
		for(int i = (!ghosts)*NGHOST; i < imax; i++)
		{
			vtk << u->uP(4, i, j, k) << endl;
		}
	}
	if (exp == exports::E) {
		/* Write the scalar each line */
		for(int k = (!ghosts)*NGHOST; k < kmax; k++)
		for(int j = (!ghosts)*NGHOST; j < jmax; j++)
		for(int i = (!ghosts)*NGHOST; i < imax; i++)
		{
			vtk << u->uC(4, i, j, k) << endl;
		}
	}
	long double maxDiv = 0.0;
	if (exp == exports::DIVB) {
		/* Write the scalar each line */
		for(int k = (!ghosts)*NGHOST; k < kmax; k++)
		for(int j = (!ghosts)*NGHOST; j < jmax; j++)
		for(int i = (!ghosts)*NGHOST; i < imax; i++)
		{
			long double DivB = getDiv(u, i, j, k);
			if(fabsl(DivB) > 1e-12) cout << i << ", " << j << ", " << k << ": " << DivB << endl;
			vtk << DivB << endl;
			if (fabsl(DivB)>fabsl(maxDiv))
				maxDiv = DivB;
		}
		cout << "Maximum divergence: " << maxDiv << endl;
	}
	



	vtk << endl; vtk.close();
	return true;
}

bool vtk_3d(Arrays *u, bool ghosts, exports *exp, int Nexp)
{

	int Nx = u->Nx+ghosts*2*NGHOST*(u->Nx>1); int Ny = u->Ny+ghosts*2*NGHOST*(u->Ny>1); int Nz = u->Nz+ghosts*2*NGHOST*(u->Nz>1);

	ofstream vtk("exports/out.vtk", ios::out | ios::trunc);

	/* Headers of the .vtk file */
	vtk << "# vtk DataFile Version 2.0" << endl;
	vtk << "testing" << endl;
	vtk << "ASCII" << endl << endl;

	vtk << "DATASET STRUCTURED_POINTS" << endl;
	vtk << "DIMENSIONS " << Nx << " " << Ny << " " << Nz << endl;
	vtk << "ORIGIN " << u->x0-ghosts*NGHOST*u->dx*(u->Nx>1) << " " << u->y0-ghosts*NGHOST*u->dy*(u->Ny>1) << " " << u->z0-ghosts*NGHOST*u->dz*(u->Nz>1) << endl;
	vtk << "SPACING " << u->dx << " " << u->dy << " " << u->dz << endl << endl;

	vtk << "POINT_DATA " << Nx*Ny*Nz << endl;
	vtk.close();


	
	for (int e = 0; e < Nexp; e++)
	{
		exports expi = exp[e];
		if (expi == exports::CONS) {

			vtk_add(u, ghosts, exports::Rho);
			vtk_add(u, ghosts, exports::Mvec);
			vtk_add(u, ghosts, exports::E);
#ifdef MHD
			vtk_add(u, ghosts, exports::Bvec);
#endif // MHD

		}
		else if (expi == exports::PRIM) {

			vtk_add(u, ghosts, exports::Rho);
			vtk_add(u, ghosts, exports::Vvec);
			vtk_add(u, ghosts, exports::P);
#ifdef MHD
			vtk_add(u, ghosts, exports::Bvec);
#endif // MHD

		}
		else {
			vtk_add(u, ghosts, expi);
		}
	}

	




	vtk.close();
	return true;
}