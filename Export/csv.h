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
		/* Write the scalar each line */
		for(int k = kmin; k <= kmax; k++)
		for(int j = jmin; j <= jmax; j++)
		for(int i = imin; i <= imax; i++)
		{
			csv << u->uC(0, i, j, k);
			if (k != kmax || j != jmax || i != imax) csv << ",";
		}
		csv << endl;
	}


	if (out == exports::V) {
		/* Write the scalar each line */
		for(int k = kmin; k <= kmax; k++)
		for(int j = jmin; j <= jmax; j++)
		for(int i = imin; i <= imax; i++)
		{
			csv << sqrtl( SQ(u->uP(1, i, j, k)) + SQ(u->uP(2, i, j, k)) + SQ(u->uP(3, i, j, k)) );
			if (k != kmax || j != jmax || i != imax) csv << ",";
		}
		csv << endl;
	}
	if (out == exports::VX) {
		/* Write the scalar each line */
		for(int k = kmin; k <= kmax; k++)
		for(int j = jmin; j <= jmax; j++)
		for(int i = imin; i <= imax; i++)
		{
			csv << u->uP(1, i, j, k);
			if (k != kmax || j != jmax || i != imax) csv << ",";
		}
		csv << endl;
	}
	if (out == exports::VY) {
		/* Write the scalar each line */
		for(int k = kmin; k <= kmax; k++)
		for(int j = jmin; j <= jmax; j++)
		for(int i = imin; i <= imax; i++)
		{
			csv << u->uP(2, i, j, k);
			if (k != kmax || j != jmax || i != imax) csv << ",";
		}
		csv << endl;
	}
	if (out == exports::VZ) {
		/* Write the scalar each line */
		for(int k = kmin; k <= kmax; k++)
		for(int j = jmin; j <= jmax; j++)
		for(int i = imin; i <= imax; i++)
		{
			csv << u->uP(3, i, j, k);
			if (k != kmax || j != jmax || i != imax) csv << ",";
		}
		csv << endl;
	}


	if (out == exports::M) {
		/* Write the scalar each line */
		for(int k = kmin; k <= kmax; k++)
		for(int j = jmin; j <= jmax; j++)
		for(int i = imin; i <= imax; i++)
		{
			csv << sqrtl( SQ(u->uC(1, i, j, k)) + SQ(u->uC(2, i, j, k)) + SQ(u->uC(3, i, j, k)) );
			if (k != kmax || j != jmax || i != imax) csv << ",";
		}
		csv << endl;
	}
	if (out == exports::MX) {
		/* Write the scalar each line */
		for(int k = kmin; k <= kmax; k++)
		for(int j = jmin; j <= jmax; j++)
		for(int i = imin; i <= imax; i++)
		{
			csv << u->uC(1, i, j, k);
			if (k != kmax || j != jmax || i != imax) csv << ",";
		}
		csv << endl;
	}
	if (out == exports::MY) {
		/* Write the scalar each line */
		for(int k = kmin; k <= kmax; k++)
		for(int j = jmin; j <= jmax; j++)
		for(int i = imin; i <= imax; i++)
		{
			csv << u->uC(2, i, j, k);
			if (k != kmax || j != jmax || i != imax) csv << ",";
		}
		csv << endl;
	}
	if (out == exports::MZ) {
		/* Write the scalar each line */
		for(int k = kmin; k <= kmax; k++)
		for(int j = jmin; j <= jmax; j++)
		for(int i = imin; i <= imax; i++)
		{
			csv << u->uC(3, i, j, k);
			if (k != kmax || j != jmax || i != imax) csv << ",";
		}
		csv << endl;
	}


	if (out == exports::B) {
		/* Write the scalar each line */
		for(int k = kmin; k <= kmax; k++)
		for(int j = jmin; j <= jmax; j++)
		for(int i = imin; i <= imax; i++)
		{
			csv << sqrtl( SQ(u->uP(5, i, j, k)) + SQ(u->uP(6, i, j, k)) + SQ(u->uP(7, i, j, k)) );
			if (k != kmax || j != jmax || i != imax) csv << ",";
		}
		csv << endl;
	}
	if (out == exports::BX) {
		/* Write the scalar each line */
		for(int k = kmin; k <= kmax; k++)
		for(int j = jmin; j <= jmax; j++)
		for(int i = imin; i <= imax; i++)
		{
			csv << u->uP(5, i, j, k);
			if (k != kmax || j != jmax || i != imax) csv << ",";
		}
		csv << endl;
	}
	if (out == exports::BY) {
		/* Write the scalar each line */
		for(int k = kmin; k <= kmax; k++)
		for(int j = jmin; j <= jmax; j++)
		for(int i = imin; i <= imax; i++)
		{
			csv << u->uP(6, i, j, k);
			if (k != kmax || j != jmax || i != imax) csv << ",";
		}
		csv << endl;
	}
	if (out == exports::BZ) {
		/* Write the scalar each line */
		for(int k = kmin; k <= kmax; k++)
		for(int j = jmin; j <= jmax; j++)
		for(int i = imin; i <= imax; i++)
		{
			csv << u->uP(7, i, j, k);
			if (k != kmax || j != jmax || i != imax) csv << ",";
		}
		csv << endl;
	}
	

	if (out == exports::J) { // Magnitude of the Current density
		/* Write the scalar each line */
		for(int k = kmin; k <= kmax; k++)
		for(int j = jmin; j <= jmax; j++)
		for(int i = imin; i <= imax; i++)
		{
			ld dBxdy = 0.0; ld dBxdz = 0.0;   ld dBydx = 0.0; ld dBydz = 0.0;   ld dBzdx = 0.0; ld dBzdy = 0.0;
			if (u->Nx > 1) {
				dBydx = (u->uC(6, i+1, j, k)-u->uC(6, i, j, k))/u->dx;   dBzdx = (u->uC(7, i+1, j, k)-u->uC(7, i, j, k))/u->dx;
			}
			if (u->Ny > 1) {
				dBxdy = (u->uC(5, i, j+1, k)-u->uC(5, i, j, k))/u->dy;   dBzdy = (u->uC(7, i, j+1, k)-u->uC(7, i, j, k))/u->dy;
			}
			if (u->Nz > 1) {
				dBxdz = (u->uC(5, i, j, k+1)-u->uC(5, i, j, k))/u->dz;   dBydz = (u->uC(6, i, j, k)-u->uC(6, i, j, k+1))/u->dz;
			}
			csv << sqrtl( SQ(dBzdy - dBydz) + SQ(dBxdz - dBzdx) + SQ(dBydx - dBxdy) );
			if (k != kmax || j != jmax || i != imax) csv << ",";
		}
		csv << endl;
	}
	if (out == exports::Jx) { // Magnitude of the Current density
		/* Write the scalar each line */
		for(int k = kmin; k <= kmax; k++)
		for(int j = jmin; j <= jmax; j++)
		for(int i = imin; i <= imax; i++)
		{
			ld dBydz = 0.0; ld dBzdy = 0.0;

			if (u->Ny > 1)
				dBzdy = (u->uC(7, i, j+1, k)-u->uC(7, i, j, k))/u->dy;
			if (u->Nz > 1)
				dBydz = (u->uC(6, i, j, k)-u->uC(6, i, j, k+1))/u->dz;
			
			csv << dBzdy - dBydz;
			if (k != kmax || j != jmax || i != imax) csv << ",";
		}
		csv << endl;
	}
	if (out == exports::Jy) { // Magnitude of the Current density
		/* Write the scalar each line */
		for(int k = kmin; k <= kmax; k++)
		for(int j = jmin; j <= jmax; j++)
		for(int i = imin; i <= imax; i++)
		{
			ld dBxdz = 0.0; ld dBzdx = 0.0;
			if (u->Nx > 1)
				dBzdx = (u->uC(7, i+1, j, k)-u->uC(7, i, j, k))/u->dx;
			if (u->Nz > 1)
				dBxdz = (u->uC(5, i, j, k+1)-u->uC(5, i, j, k))/u->dz;
			
			csv << dBxdz - dBzdx;
			if (k != kmax || j != jmax || i != imax) csv << ",";
		}
		csv << endl;
	}
	if (out == exports::Jz) { // Magnitude of the Current density
		/* Write the scalar each line */
		for(int k = kmin; k <= kmax; k++)
		for(int j = jmin; j <= jmax; j++)
		for(int i = imin; i <= imax; i++)
		{
			ld dBxdy = 0.0; ld dBydx = 0.0;
			if (u->Nx > 1)
				dBydx = (u->uC(6, i+1, j, k)-u->uC(6, i, j, k))/u->dx;
			if (u->Ny > 1)
				dBxdy = (u->uC(5, i, j+1, k)-u->uC(5, i, j, k))/u->dy;
			csv << dBydx - dBxdy;
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
	ld maxDiv = 0.0;
	if (out == exports::DIVB) {
		/* Write the scalar each line */
		for(int k = kmin; k <= kmax; k++)
		for(int j = jmin; j <= jmax; j++)
		for(int i = imin; i <= imax; i++)
		{
			ld DivB = getDiv(u, i, j, k);
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

bool csv(Arrays *u, vector<Export> *out)
{

	for (int i = 0; i < out->size(); i++) {
		Export E = out->at(i);
	//for (Export E : out) {
		if (E.type != 1) continue; /* if not csv */

		/* Check if we need to export based on dt */
		bool toexp = false;
		for (double t = 0; t <= u->t+2*u->dt; t+=E.dt)
		{
			/* I'm going to have to find a new way of doing this. */
			if ((u->t < t && u->t+u->dt >= t) || u->t == 0.0) {
				/* If the multiple of the export Dt is between t and t+dt, then export now. Else continue to next variable. */
				toexp = true; break;
			}
		}
		if (toexp == false) continue;
		

		string fname = "exports/"+E.name+"_"+to_string(E.id)+".csv";
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


		out->at(i).id += 1;
	}

	
	return true;
}