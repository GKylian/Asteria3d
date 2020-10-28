#pragma once

#include <iostream>
#include <fstream>
#include <string>
using namespace std;

#include "../Defs/Arrays.h"

bool vtk_add(Arrays *u, bool ghosts, exports out, string fname) 
{
	
	ofstream vtk(fname, ios::out | ios::app);
	if (out == exports::Bvec || out == exports::Vvec || out == exports::Mvec) {
		/* Write header for  */
		vtk << "VECTORS " << exp_names[(int)out] << " float" << endl;
	}
	else {
		/* Write header for scalar */
		vtk << "SCALARS " << exp_names[(int)out] << " float" << endl;
		vtk << "LOOKUP_TABLE default" << endl;
	}

	/* Store the 'for' loops' limits to avoid computing them every iteration */
	//int imax = u->Nx+(!ghosts)*NGHOST; int jmax = u->Ny+(!ghosts)*NGHOST; int kmax = u->Nz+(!ghosts)*NGHOST;
	int imin = ghosts*u->i_dl+(!ghosts)*u->i_cl; int imax = ghosts*u->i_dr+(!ghosts)*u->i_cr;
	int jmin = ghosts*u->j_dl+(!ghosts)*u->j_cl; int jmax = ghosts*u->j_dr+(!ghosts)*u->j_cr;
	int kmin = ghosts*u->k_dl+(!ghosts)*u->k_cl; int kmax = ghosts*u->k_dr+(!ghosts)*u->k_cr;




	/* We put the ifs outside the loops instead of inside a single big
	   one for better performance (1 'if' instead of Nx*Ny*Nz 'if's)   */

	/* -------------------- VECTOR EXPORT -------------------- */

	if (out == exports::Bvec) {
		/* Write the three components each line */
		for(int k = kmin; k <= kmax; k++)
		for(int j = jmin; j <= jmax; j++)
		for(int i = imin; i <= imax; i++)
		{
			vtk << u->uP(5, i, j, k) << " " << u->uP(6, i, j, k) << " " << u->uP(7, i, j, k) << endl;
		}
	}
	if (out == exports::Vvec) {
		/* Write the three components each line */
		for(int k = kmin; k <= kmax; k++)
		for(int j = jmin; j <= jmax; j++)
		for(int i = imin; i <= imax; i++)
		{
			vtk << u->uP(1, i, j, k) << " " << u->uP(2, i, j, k) << " " << u->uP(3, i, j, k) << endl;
		}
	}
	if (out == exports::Mvec) {
		/* Write the three components each line */
		for(int k = kmin; k <= kmax; k++)
		for(int j = jmin; j <= jmax; j++)
		for(int i = imin; i <= imax; i++)
		{
			vtk << u->uC(1, i, j, k) << " " << u->uC(2, i, j, k) << " " << u->uC(3, i, j, k) << endl;
		}
	}
	




	/* -------------------- SCALAR EXPORT -------------------- */

	if (out == exports::Rho) {
		/* Write the scalar each line */
		for(int k = imin; k <= kmax; k++)
		for(int j = jmin; j <= jmax; j++)
		for(int i = kmin; i <= imax; i++)
		{
			//cout << i << ", " << j << endl;
			//cout << u->uC(0, i, j, k) << endl;
			vtk << u->uC(0, i, j, k) << endl;
		}
	}
	if (out == exports::P) {
		/* Write the scalar each line */
		for(int k = kmin; k <= kmax; k++)
		for(int j = jmin; j <= jmax; j++)
		for(int i = imin; i <= imax; i++)
		{
			vtk << u->uP(4, i, j, k) << endl;
		}
	}
	if (out == exports::E) {
		/* Write the scalar each line */
		for(int k = kmin; k <= kmax; k++)
		for(int j = jmin; j <= jmax; j++)
		for(int i = imin; i <= imax; i++)
		{
			vtk << u->uC(4, i, j, k) << endl;
		}
	}
#ifdef MHD
	ld maxDiv = 0.0;
	if (out == exports::DIVB) {
		/* Write the scalar each line */
		for(int k = kmin; k <= kmax; k++)
		for(int j = jmin; j <= jmax; j++)
		for(int i = imin; i <= imax; i++)
		{
			ld DivB = getDiv(u, i, j, k);
			vtk << DivB << endl;
			if (fabsl(DivB)>fabsl(maxDiv))
				maxDiv = DivB;
		}
		cout << "Maximum divergence: " << maxDiv << endl;
	}
#endif



	vtk << endl; vtk.close();
	return true;
}

bool vtk_3d(Arrays *u, vector<Export> out)
{

	
	for (Export E : out) {
		if (E.type != 0) continue; /* if not vtk */
		string fname = "exports/"+E.name+"_"+to_string(u->s)+".vtk";
		cout << "EXPORT:: Exporting to vtk file " << fname << endl;
		ofstream vtk(fname, ios::out | ios::trunc);

		int Nx = u->Nx+E.ghosts*2*NGHOST*(u->Nx>1); int Ny = u->Ny+E.ghosts*2*NGHOST*(u->Ny>1); int Nz = u->Nz+E.ghosts*2*NGHOST*(u->Nz>1);

		/* Headers of the .vtk file */
		vtk << "# vtk DataFile Version 2.0" << endl;
		vtk << "testing" << endl;
		vtk << "ASCII" << endl << endl;

		vtk << "DATASET STRUCTURED_POINTS" << endl;
		vtk << "DIMENSIONS " << Nx << " " << Ny << " " << Nz << endl;
		vtk << "ORIGIN " << u->x0-E.ghosts*NGHOST*u->dx*(u->Nx>1) << " " << u->y0-E.ghosts*NGHOST*u->dy*(u->Ny>1) << " " << u->z0-E.ghosts*NGHOST*u->dz*(u->Nz>1) << endl;
		vtk << "SPACING " << u->dx << " " << u->dy << " " << u->dz << endl << endl;

		vtk << "POINT_DATA " << Nx*Ny*Nz << endl;
		vtk.close();

		for (exports e : E.exp) {

			if (e == exports::CONS) {

				vtk_add(u, E.ghosts, exports::Rho, fname);
				vtk_add(u, E.ghosts, exports::Mvec, fname);
				vtk_add(u, E.ghosts, exports::E, fname);
#ifdef MHD
				vtk_add(u, E.ghosts, exports::Bvec, fname);
#endif // MHD

			}
			else if (e == exports::PRIM) {

				vtk_add(u, E.ghosts, exports::Rho, fname);
				vtk_add(u, E.ghosts, exports::Vvec, fname);
				vtk_add(u, E.ghosts, exports::P, fname);
#ifdef MHD
				vtk_add(u, E.ghosts, exports::Bvec, fname);
#endif // MHD

			}
			else {
				vtk_add(u, E.ghosts, e, fname);
			}
		}
	}

	
	return true;
}