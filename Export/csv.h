#pragma once

#include <iostream>
#include <fstream>
#include <string>
using namespace std;

#include "../Defs/Arrays.h"

bool csv_add(Arrays *u, bool ghosts, exports out, string fname)
{
	
	ofstream csv(fname, ios::out | ios::app);
	//if (exp == exports::Bvec || exp == exports::Vvec || exp == exports::Mvec) {
	//	/* Write header for  */
	//	csv << "VECTORS " << exp_names[(int)exp] << " float" << endl;
	//}
	//else {
	//	/* Write header for scalar */
	//	csv << "SCALARS " << exp_names[(int)exp] << " float" << endl;
	//	csv << "LOOKUP_TABLE default" << endl;
	//}

	/* Store the 'for' loops' limits to avoid computing them every iteration */
	//int imax = u->Nx+(!ghosts)*NGHOST; int jmax = u->Ny+(!ghosts)*NGHOST; int kmax = u->Nz+(!ghosts)*NGHOST;
	int imin = ghosts*u->i_dl+(!ghosts)*u->i_cl; int imax = ghosts*u->i_dr+(!ghosts)*u->i_cr;
	int jmin = ghosts*u->j_dl+(!ghosts)*u->j_cl; int jmax = ghosts*u->j_dr+(!ghosts)*u->j_cr;
	int kmin = ghosts*u->k_dl+(!ghosts)*u->k_cl; int kmax = ghosts*u->k_dr+(!ghosts)*u->k_cr;




	/* We put the ifs outside the loops instead of inside a single big
	   one for better performance (1 'if' instead of Nx*Ny*Nz 'if's)   */

	/* -------------------- VECTOR EXPORT -------------------- */

	if (out == exports::Bvec) {
		int x = (out==exports::Bvec)*5+(out!=exports::Bvec)*1;
		int y = (out==exports::Bvec)*6+(out!=exports::Bvec)*2;
		int z = (out==exports::Bvec)*7+(out!=exports::Bvec)*3;

		/* Write one line for each component */
		for(int k = kmin; k <= kmax; k++)
		for(int j = jmin; j <= jmax; j++)
		for(int i = imin; i <= imax; i++)
		{
			csv << u->uP(x, i, j, k);
			if (k != kmax || j != jmax || i != imax) csv << ",";
		} csv << endl;
		for(int k = kmin; k <= kmax; k++)
		for(int j = jmin; j <= jmax; j++)
		for(int i = imin; i <= imax; i++)
		{
			csv << u->uP(y, i, j, k);
			if (k != kmax || j != jmax || i != imax) csv << ",";
		} csv << endl;
		for(int k = kmin; k <= kmax; k++)
		for(int j = jmin; j <= jmax; j++)
		for(int i = imin; i <= imax; i++)
		{
			csv << u->uP(z, i, j, k);
			if (k != kmax || j != jmax || i != imax) csv << ",";
		} csv << endl;
	}


	if (out == exports::Mvec) {
		int x = (out==exports::Bvec)*5+(out!=exports::Bvec)*1;
		int y = (out==exports::Bvec)*6+(out!=exports::Bvec)*2;
		int z = (out==exports::Bvec)*7+(out!=exports::Bvec)*3;

		/* Write one line for each component */
		for(int k = kmin; k <= kmax; k++)
		for(int j = jmin; j <= jmax; j++)
		for(int i = imin; i <= imax; i++)
		{
			csv << u->uC(x, i, j, k);
			if (k != kmax || j != jmax || i != imax) csv << ",";
		} csv << endl;
		for(int k = kmin; k <= kmax; k++)
		for(int j = jmin; j <= jmax; j++)
		for(int i = imin; i <= imax; i++)
		{
			csv << u->uC(y, i, j, k);
			if (k != kmax || j != jmax || i != imax) csv << ",";
		} csv << endl;
		for(int k = kmin; k <= kmax; k++)
		for(int j = jmin; j <= jmax; j++)
		for(int i = imin; i <= imax; i++)
		{
			csv << u->uC(z, i, j, k);
			if (k != kmax || j != jmax || i != imax) csv << ",";
		} csv << endl;
	}
	




	/* -------------------- SCALAR EXPORT -------------------- */

	if (out == exports::Rho) {
		cout << "EXPORT::csv.h:: Exporting the density" << endl;
		
		/* Write the scalar each line */
		for(int k = kmin; k <= kmax; k++)
		for(int j = jmin; j <= jmax; j++)
		for(int i = imin; i <= imax; i++)
		{
			//cout << i << ", " << j << endl;
			//cout << u->uC(0, i, j, k) << endl;
			csv << u->uC(0, i, j, k);
			if (k != kmax || j != jmax || i != imax) csv << ",";
		}
		csv << endl;
	}
	if (out == exports::P) {
		/* Write the scalar each line */
		for(int k = kmin; k <= kmax; k++)
		for(int j = jmin; j <= jmax; j++)
		for(int i = imin; i <= imax; i++)
		{
			csv << u->uP(4, i, j, k);
			if (k != kmax || j != jmax || i != imax) csv << ",";
		}
		csv << endl;
	}
	if (out == exports::E) {
		/* Write the scalar each line */
		for(int k = kmin; k <= kmax; k++)
		for(int j = jmin; j <= jmax; j++)
		for(int i = imin; i <= imax; i++)
		{
			csv << u->uC(4, i, j, k);
			if (k != kmax || j != jmax || i != imax) csv << ",";
		}
		csv << endl;
	}

#ifdef MHD
	long double maxDiv = 0.0;
	if (out == exports::DIVB) {
		/* Write the scalar each line */
		for(int k = kmin; k <= kmax; k++)
		for(int j = jmin; j <= jmax; j++)
		for(int i = imin; i <= imax; i++)
		{
			long double DivB = getDiv(u, i, j, k);
			csv << DivB;
			if (k != kmax || j != jmax || i != imax) csv << ",";
			if (fabsl(DivB)>fabsl(maxDiv))
				maxDiv = DivB;
		}
		csv << endl;
	}
#endif



	csv << endl; csv.close();
	return true;
}

bool csv(Arrays *u, vector<Export> out)
{

	
	for (Export E : out) {
		if (E.type != 1) continue; /* if not csv */
		string fname = "exports/"+E.name+"_"+to_string(u->s)+".csv";
		cout << "EXPORT:: Exporting to csv file " << fname << endl;
		ofstream csv(fname, ios::out | ios::trunc);

		int Nx = u->Nx+E.ghosts*2*NGHOST*(u->Nx>1); int Ny = u->Ny+E.ghosts*2*NGHOST*(u->Ny>1); int Nz = u->Nz+E.ghosts*2*NGHOST*(u->Nz>1);

		/* Headers of the .csv file */
		csv << Nx << "," << Ny << "," << Nz << endl;
		csv << u->dx << "," << u->dy << "," << u->dz << endl;
		csv << u->x0-E.ghosts*NGHOST*u->dx*(u->Nx>1) << "," << u->y0-E.ghosts*NGHOST*u->dy*(u->Ny>1) << "," << u->z0-E.ghosts*NGHOST*u->dz*(u->Nz>1) << endl;
		csv.close();

		for (exports e : E.exp) {

			if (e == exports::CONS) {

				csv_add(u, E.ghosts, exports::Rho, fname);
				csv_add(u, E.ghosts, exports::Mvec, fname);
				csv_add(u, E.ghosts, exports::E, fname);
#ifdef MHD
				csv_add(u, E.ghosts, exports::Bvec, fname);
#endif // MHD

			}
			else if (e == exports::PRIM) {

				csv_add(u, E.ghosts, exports::Rho, fname);
				csv_add(u, E.ghosts, exports::Vvec, fname);
				csv_add(u, E.ghosts, exports::P, fname);
#ifdef MHD
				csv_add(u, E.ghosts, exports::Bvec, fname);
#endif // MHD

			}
			else {
				csv_add(u, E.ghosts, e, fname);
			}
		}
	}

	
	return true;
}