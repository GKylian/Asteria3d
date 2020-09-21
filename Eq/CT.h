#pragma once
#include "../Defs/Arrays.h"


/* -------------------------------------- */
/* ---------- Reference fields ---------- */
/* -------------------------------------- */

long double Erx(Arrays *u, int i, int j, int k) {
	//Erx = -(v x B)_x = -(vy*Bz - vz*By)
	return - ( u->uP(2, i, j, k)*u->uP(NVAL-1, i, j, k) - u->uP(3, i, j, k)*u->uP(NVAL-2, i, j, k) );
}

long double Ery(Arrays *u, int i, int j, int k) {
	//Erx = -(v x B)_y = (vx*Bz - vz*Bx)
	return - ( u->uP(1, i, j, k)*u->uP(NVAL-1, i, j, k) - u->uP(3, i, j, k)*u->uP(NVAL-3, i, j, k) );
}

long double Erz(Arrays *u, int i, int j, int k) {
	//Erx = -(v x B)_z = -(vx*By - vy*Bx)
	return - ( u->uP(1, i, j, k)*u->uP(NVAL-2, i, j, k) - u->uP(2, i, j, k)*u->uP(NVAL-3, i, j, k) );
}


/* ----- Get Ex at (i, j-1/2, k-1/2) ----- */
long double getEx(Arrays *u, int i, int j, int k) {

}


/* ----- Get Ey at (i-1/2, j, k-1/2) ----- */
long double getEy(Arrays *u, int i, int j, int k) {

}


/* ----- Get Ey at (i-1/2, j-1/2, k) ----- */
long double getEz(Arrays *u, int i, int j, int k) {

}




bool getEMFs(Arrays *u, long double *E, int i, int j, int k) {

	if (u->Ny > 1 && u->Nz > 1) E[0] = getEx(u);
	if (u->Nx > 1 && u->Nz > 1) E[1] = getEy(u);
	if (u->Nx > 1 && u->Ny > 1) E[2] = getEz(u);
}

