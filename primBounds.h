#pragma once

#include "Arrays.h"

using namespace std;
/*
<-------------------- Functions -------------------->
<->
---> bool allBounds(Arrays *u)
	 <-> Applies all the boundary conditions (BCs) by calling the corresponding functions.
	 - Arrays *u: the arrays containing the primitive values we're applying the BCs to.
<->
<->
---> bool bounds_xin(Arrays *u)
	 <-> Apply the BC at the inner boundary of x (x0).
<->
---> bool bounds_xout(Arrays *u)
	 <-> Apply the BC at the outer boundary of x (xn).
<->
*/
bool allBounds(Arrays *u);
bool bounds_xin(Arrays *u);
bool bounds_xout(Arrays *u);
bool bounds_yin(Arrays *u);
bool bounds_yout(Arrays *u);
bool bounds_zin(Arrays *u);
bool bounds_zout(Arrays *u);







bool allBounds(Arrays *u) {

	/* Call boundary functions if needed (N > 1) */
	if (u->Nx > 1) {
		if (!bounds_xin(u))  { cout << "primBounds::allBounds:: Could not set the boundaries for x_in !"  << endl; return false; }
		if (!bounds_xout(u)) { cout << "primBounds::allBounds:: Could not set the boundaries for x_out !" << endl; return false; }
	}
	if (u->Ny > 1) {
		if (!bounds_yin(u))  { cout << "primBounds::allBounds:: Could not set the boundaries for y_in !"  << endl; return false; }
		if (!bounds_yout(u)) { cout << "primBounds::allBounds:: Could not set the boundaries for y_out !" << endl; return false; }
	}
	if (u->Nz > 1) {
		if (!bounds_zin(u))  { cout << "primBounds::allBounds:: Could not set the boundaries for z_in !"  << endl; return false; }
		if (!bounds_zout(u)) { cout << "primBounds::allBounds:: Could not set the boundaries for z_out !" << endl; return false; }
	}


	/* For the whole grid (including ghost cells), compute the cell-centered magnetic fields (if MHD) and the conserved variables */
	for(int k = 0; k < u->Nz+2*NGHOST; k++)
	for(int j = 0; j < u->Ny+2*NGHOST; j++)
	for(int i = 0; i < u->Nx+2*NGHOST; i++)
	{
#ifdef MHD
		BCenter(u, i, j, k);
#endif // MHD
		toConserved(u, i, j, k, gamma);

	}


	return true;
}





bool bounds_xin(Arrays *u) {
	
	if (u->boundaries[0] == bounds::PERIODIC) {

		/* Loop through the transverse dimensions (y, z) */
		for(int k = 0; k < u->Nz+2*NGHOST*(u->Nz>1); k++)
		for(int j = 0; j < u->Ny+2*NGHOST*(u->Ny>1); j++)
		{
			/* Loop though the ghost cells for (M)HD variables */
			for (int i = 0; i < NGHOST; i++) {
				u->uP(0, i, j, k) = u->uP(0, u->Nx+i, j, k);

				u->uP(1, i, j, k) = u->uP(1, u->Nx+i, j, k); u->uP(2, i, j, k) = u->uP(2, u->Nx+i, j, k); u->uP(3, i, j, k) = u->uP(3, u->Nx+i, j, k);
				u->uP(4, i, j, k) = u->uP(4, u->Nx+i, j, k);
#ifdef MHD
				u->uC(6, i, j, k) = u->uC(6, u->Nx+i, j, k); u->uC(7, i, j, k) = u->uC(7, u->Nx+i, j, k);
				u->uC(5, i, j, k) = u->uC(5, u->Nx+i, j, k);
#endif // MHD

			}

		}
	}

	return true;
}


bool bounds_xout(Arrays *u) {


	if (u->boundaries[1] == bounds::PERIODIC) {

		/* Loop through the transverse dimensions (y, z) */
		for(int k = u->k_dl; k <= u->k_dr; k++)
		for(int j = u->j_dl; j <= u->j_dr; j++)
		{
			/* Loop though the ghost cells for (M)HD variables */
			for (int i = 0; i <= u->i_gl; i++) {
				u->uP(0, u->i_gr+i, j, k) = u->uP(0, u->i_cl+i, j, k);
				u->uP(1, u->i_gr+i, j, k) = u->uP(1, u->i_cl+i, j, k); u->uP(2, u->i_gr+i, j, k) = u->uP(2, u->i_cl+i, j, k); u->uP(3, u->i_gr+i, j, k) = u->uP(3, u->i_cl+i, j, k);
				u->uP(4, u->i_gr+i, j, k) = u->uP(4, u->i_cl+i, j, k);
#ifdef MHD
				u->uC(6, u->i_gr+i, j, k) = u->uC(6, u->i_cl+i, j, k); u->uC(7, u->i_gr+i, j, k) = u->uC(7, u->i_cl+i, j, k);
				u->uC(5, u->i_gr+i, j, k) = u->uC(5, u->i_cl+i, j, k);
				
#endif // MHD

			}
#ifdef MHD
			u->uC(5, u->i_dr+1, j, k) = u->uC(5, u->i_cl+u->i_cl, j, k);
#endif // MHD


		}
	}

	return true;
}





bool bounds_yin(Arrays *u) {

	if (u->boundaries[2] == bounds::PERIODIC) {

		/* Loop through the transverse dimensions (x,z) */
		for(int k = 0; k < u->Nz+2*NGHOST*(u->Nz>1); k++)
		for(int i = 0; i < u->Nx+2*NGHOST*(u->Nx>1); i++)
		{
			/* Loop through the ghost cells */
			for (int j = 0; j < NGHOST; j++) {
				u->uP(0, i, j, k) = u->uP(0, i, u->Ny+j, k);
				u->uP(1, i, j, k) = u->uP(1, i, u->Ny+j, k); u->uP(2, i, j, k) = u->uP(2, i, u->Ny+j, k); u->uP(3, i, j, k) = u->uP(3, i, u->Ny+j, k);
				u->uP(4, i, j, k) = u->uP(4, i, u->Ny+j, k);
#ifdef MHD
				u->uC(5, i, j, k) = u->uC(5, i, u->Ny+j, k); u->uC(7, i, j, k) = u->uC(7, i, u->Ny+j, k);
				u->uC(6, i, j, k) = u->uC(6, i, u->Ny+j, k);
#endif // MHD

			}
		}


	}

	return true;
}


bool bounds_yout(Arrays *u) {

	if (u->boundaries[3] == bounds::PERIODIC) {

		/* Loop through the transverse dimensions (y, z) */

		for(int k = u->k_dl; k <= u->k_dr; k++)
		for(int i = u->i_dl; i <= u->i_dr; i++)
		{
			/* Loop though the ghost cells for (M)HD variables */
			for (int j = 0; j <= u->j_gl; j++) {
				u->uP(0, i, u->j_gr+j, k) = u->uP(0, i, u->j_cl+j, k);

				u->uP(1, i, u->j_gr+j, k) = u->uP(1, i, u->j_cl+j, k); u->uP(2, i, u->j_gr+j, k) = u->uP(2, i, u->j_cl+j, k); u->uP(3, i, u->j_gr+j, k) = u->uP(3, i, u->j_cl+j, k);
				u->uP(4, i, u->j_gr+j, k) = u->uP(4, i, u->j_cl+j, k);
#ifdef MHD
				u->uC(5, i, u->j_gr+j, k) = u->uC(5, i, u->j_cl+j, k); u->uC(7, i, u->j_gr+j, k) = u->uC(7, i, u->j_cl+j, k);
				//u->uC(5, u->i_gr+i, j, k) = u->uC(5, u->i_cl+i, j, k);
				u->uC(6, i, u->j_gr+j, k) = u->uC(6, i, u->j_cl+j, k);
				//cout << u->uC(6, i, u->j_gr+j, k) << " --> " << u->uC(6, i, u->j_cl+j, k) << endl;
#endif // MHD

			}
#ifdef MHD
			//u->uC(5, u->i_dr+1, j, k) = u->uC(5, u->i_cl, j, k); 
			u->uC(6, i, u->j_dr+1, k) = u->uC(6, i, u->j_cl+u->j_cl, k);
#endif // MHD


		}
	}

	return true;
}





bool bounds_zin(Arrays *u) {


	return true;
}


bool bounds_zout(Arrays *u) {


	return true;
}