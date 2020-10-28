#pragma once
#include "../Defs/Arrays.h"


/* -------------------------------------- */
/* ---------- Reference fields ---------- */
/* -------------------------------------- */

ld Erx(Arrays *u, int i, int j, int k) {
	//Erx = -(v x B)_x = -(vy*Bz - vz*By)
	return - ( u->uP(2, i, j, k)*u->uP(NVAL-1, i, j, k) - u->uP(3, i, j, k)*u->uP(NVAL-2, i, j, k) );
}

ld Ery(Arrays *u, int i, int j, int k) {
	//Ery = -(v x B)_y = (vx*Bz - vz*Bx)
	return  ( u->uP(1, i, j, k)*u->uP(NVAL-1, i, j, k) - u->uP(3, i, j, k)*u->uP(NVAL-3, i, j, k) );
}

ld Erz(Arrays *u, int i, int j, int k) {
	//Erz = -(v x B)_z = -(vx*By - vy*Bx)
	return - ( u->uP(1, i, j, k)*u->uP(NVAL-2, i, j, k) - u->uP(2, i, j, k)*u->uP(NVAL-3, i, j, k) );
}


/* ----- Get Ex at (i, j-1/2, k-1/2) ----- */
//TODO: Changed the initial four, the 'ifs', the dEydz and dEydx (except averages), the averages   ->   SHOULD BE OK
ld getEx(Arrays *u, int i, int j, int k) {
	int iBx = NVAL-3; int iBy = NVAL-2; int iBz = NVAL-1;
	ld Ex = 0.0;
	Ex += 0.25*(  u->Fz(iBy, i, j-1, k) + u->Fz(iBy, i, j, k) + (-1.0)*u->Fy(iBz, i, j, k-1) + (-1.0)*u->Fy(iBz, i, j, k)  );

	//TODO: Here we're using -Fy[iBz] and Fz[iBy]

	/* Add correct term for dEx/dz (i-1/2, j-3/4, k) */
	if (u->Fy(0, i, j, k-1) > 0) {
		Ex += u->dz/8.0 * (u->Fz(iBy, i, j-1, k) - Erx(u, i, j-1, k-1))/(u->dz/2.0);
	}
	else if (u->Fy(0, i, j, k-1) < 0) {
		Ex += u->dz/8.0 * (u->Fz(iBy, i, j, k) - Erx(u, i, j, k-1))/(u->dz/2.0);
	}
	else {
		/* Average of the two */
		Ex += 0.5*u->dz/8.0*(   (u->Fz(iBy, i, j-1, k) - Erx(u, i, j-1, k-1))/(u->dz/2.0)   +   (u->Fz(iBy, i, j, k) - Erx(u, i, j, k-1))/(u->dz/2.0)   );
	}

	/* Add correct term for dEx/dz (i-1/2, j-1/4, k) */
	if (u->Fy(0, i, j, k) > 0) {
		Ex += -u->dz/8.0 * (Erx(u, i, j-1, k) - u->Fz(iBy, i, j-1, k))/(u->dz/2.0);
	}
	else if (u->Fy(0, i, j, k) < 0) {
		Ex += -u->dz/8.0 * (Erx(u, i, j, k) - u->Fz(iBy, i, j, k))/(u->dz/2.0);
	}
	else {
		/* Average of the two */
		Ex += -0.5*u->dz/8.0*(   (Erx(u, i, j-1, k) - u->Fz(iBy, i, j-1, k))/(u->dz/2.0)   +   (Erx(u, i, j, k) - u->Fz(iBy, i, j, k))/(u->dz/2.0)   );
	}



	/* Add correct term for dEx/dy (i-3/4, j-1/2, k) */
	if (u->Fz(0, i, j-1, k) > 0) {
		Ex += u->dy/8.0 * ((-1.0)*u->Fy(iBz, i, j, k-1) - Erx(u, i, j-1, k-1))/(u->dy/2.0);
	}
	else if (u->Fz(0, i, j-1, k) < 0) {
		Ex += u->dy/8.0 * ((-1.0)*u->Fy(iBz, i, j, k) - Erx(u, i, j-1, k))/(u->dy/2.0);
	}
	else {
		/* Average of the two */
		Ex += 0.5*u->dy/8.0*(   ((-1.0)*u->Fy(iBz, i, j, k-1) - Erx(u, i, j-1, k-1))/(u->dy/2.0)   +   ((-1.0)*u->Fy(iBz, i, j, k) - Erx(u, i, j-1, k))/(u->dy/2.0)   );
	}

	/* Add correct term for dEx/dy (i-1/4, j-1/2, k) */
	if (u->Fz(0, i, j, k) > 0) {
		Ex += -u->dy/8.0 * (Erx(u, i, j, k-1) - (-1.0)*u->Fy(iBz, i, j, k-1))/(u->dy/2.0);
	}
	else if (u->Fz(0, i, j, k) < 0) {
		Ex += -u->dy/8.0 * (Erx(u, i, j, k) - (-1.0)*u->Fy(iBz, i, j, k))/(u->dy/2.0);
	}
	else {
		/* Average of the two */
		Ex += -0.5*u->dy/8.0 *(   (Erx(u, i, j, k-1) - (-1.0)*u->Fy(iBz, i, j, k-1))/(u->dy/2.0)   +   (Erx(u, i, j, k) - (-1.0)*u->Fy(iBz, i, j, k))/(u->dy/2.0)   );
	}

	return Ex;
}


/* ----- Get Ey at (i-1/2, j, k-1/2) ----- */
//TODO: Changed the initial 4, the 'ifs', the dEydz and dEydx (except averages), the averages   -> SHOULD BE OK
ld getEy(Arrays *u, int i, int j, int k) {
	int iBx = NVAL-3; int iBy = NVAL-2; int iBz = NVAL-1;
	ld Ey = 0.0;
	Ey += 0.25*(  (-1.0)*u->Fz(iBx, i-1, j, k) + (-1.0)*u->Fz(iBx, i, j, k) +u->Fx(iBz, i, j, k-1) + u->Fx(iBz, i, j, k)  );

	//return Ey;

	//TODO: Here we're using Fx[iBz] and -Fz[iBx]

	/* Add correct term for dEy/dy (i-1/2, j-3/4, k) */
	if (u->Fx(0, i, j, k-1) > 0) {
		Ey += u->dz/8.0 * ((-1.0)*u->Fz(iBx, i-1, j, k) - Ery(u, i-1, j, k-1))/(u->dz/2.0);
	}
	else if (u->Fx(0, i, j, k-1) < 0) {
		Ey += u->dz/8.0 * ((-1.0)*u->Fz(iBx, i, j, k) - Ery(u, i, j, k-1))/(u->dz/2.0);
	}
	else {
		/* Average of the two */
		Ey += 0.5*u->dz/8.0*(   ((-1.0)*u->Fz(iBx, i-1, j, k) - Ery(u, i-1, j, k-1))/(u->dz/2.0)   +   ((-1.0)*u->Fz(iBx, i, j, k) - Ery(u, i, j, k-1))/(u->dz/2.0)   );
	}

	/* Add correct term for dEy/dy (i-1/2, j-1/4, k) */
	if (u->Fx(0, i, j, k) > 0) {
		Ey += -u->dz/8.0 * (Ery(u, i-1, j, k) - (-1.0)*u->Fz(iBx, i-1, j, k))/(u->dz/2.0);
	}
	else if (u->Fx(0, i, j, k) < 0) {
		Ey += -u->dz/8.0 * (Ery(u, i, j, k) - (-1.0)*u->Fz(iBx, i, j, k))/(u->dz/2.0);
	}
	else {
		/* Average of the two */
		Ey += -0.5*u->dz/8.0*(   (Ery(u, i-1, j, k) - (-1.0)*u->Fz(iBx, i-1, j, k))/(u->dz/2.0)   +   (Ery(u, i, j, k) - (-1.0)*u->Fz(iBx, i, j, k))/(u->dz/2.0)   );
	}



	/* Add correct term for dEy/dx (i-3/4, j-1/2, k) */
	if (u->Fz(0, i-1, j, k) > 0) {
		Ey += u->dx/8.0 * (u->Fx(iBz, i, j, k-1) - Ery(u, i-1, j, k-1))/(u->dx/2.0);
	}
	else if (u->Fz(0, i-1, j, k) < 0) {
		Ey += u->dx/8.0 * (u->Fx(iBz, i, j, k) - Ery(u, i-1, j, k))/(u->dx/2.0);
	}
	else {
		/* Average of the two */
		Ey += 0.5*u->dx/8.0*(   (u->Fx(iBz, i, j, k-1) - Ery(u, i-1, j, k-1))/(u->dx/2.0)   +   (u->Fx(iBz, i, j, k) - Ery(u, i-1, j, k))/(u->dx/2.0)   );
	}

	/* Add correct term for dEy/dx (i-1/4, j-1/2, k) */
	if (u->Fz(0, i, j, k) > 0) {
		Ey += -u->dx/8.0 * (Ery(u, i, j, k-1) - u->Fx(iBz, i, j, k-1))/(u->dx/2.0);
	}
	else if (u->Fz(0, i, j, k) < 0) {
		Ey += -u->dx/8.0 * (Ery(u, i, j, k) - u->Fx(iBz, i, j, k))/(u->dx/2.0);
	}
	else {
		/* Average of the two */
		Ey += -0.5*u->dx/8.0 *(   (Ery(u, i, j, k-1) - u->Fx(iBz, i, j, k-1))/(u->dx/2.0)   +   (Ery(u, i, j, k) - u->Fx(iBz, i, j, k))/(u->dx/2.0)   );
	}

	return Ey;
}


/* ----- Get Ey at (i-1/2, j-1/2, k) ----- */ 
ld getEz(Arrays *u, int i, int j, int k) {
	int iBx = NVAL-3; int iBy = NVAL-2; int iBz = NVAL-1;
	ld Ez = 0.0;
	Ez += 0.25*(  (-1.0)*u->Fx(iBy, i, j-1, k) + (-1.0)*u->Fx(iBy, i, j, k) +u->Fy(iBx, i-1, j, k) + u->Fy(iBx, i, j, k)  );

	//return Ez;

	/* Add correct term for dEz/dy (i-1/2, j-3/4, k) */
	if (u->Fx(0, i, j-1, k) > 0) { //TODO: Once and for all verify if it's j or j-1
		Ez += u->dy/8.0 * (u->Fy(iBx, i-1, j, k) - Erz(u, i-1, j-1, k))/(u->dy/2.0);
	}
	else if (u->Fx(0, i, j-1, k) < 0) {
		Ez += u->dy/8.0 * (u->Fy(iBx, i, j, k) - Erz(u, i, j-1, k))/(u->dy/2.0);
	}
	else {
		/* Average of the two */
		Ez += 0.5*u->dy/8.0*(   (u->Fy(iBx, i-1, j, k) - Erz(u, i-1, j-1, k))/(u->dy/2.0)   +   (u->Fy(iBx, i, j, k) - Erz(u, i, j-1, k))/(u->dy/2.0)   );
	}

	/* Add correct term for dEz/dy (i-1/2, j-1/4, k) */
	if (u->Fx(0, i, j, k) > 0) {
		Ez += -u->dy/8.0 * (Erz(u, i-1, j, k) - u->Fy(iBx, i-1, j, k))/(u->dy/2.0);
	}
	else if (u->Fx(0, i, j, k) < 0) {
		Ez += -u->dy/8.0 * (Erz(u, i, j, k) - u->Fy(iBx, i, j, k))/(u->dy/2.0);
	}
	else {
		/* Average of the two */
		Ez += -0.5*u->dy/8.0*(   (Erz(u, i-1, j, k) - u->Fy(iBx, i-1, j, k))/(u->dy/2.0)   +   (Erz(u, i, j, k) - u->Fy(iBx, i, j, k))/(u->dy/2.0)   );
	}



	/* Add correct term for dEz/dx (i-3/4, j-1/2, k) */
	if (u->Fy(0, i-1, j, k) > 0) {
		Ez += u->dx/8.0 * ((-1.0)*u->Fx(iBy, i, j-1, k) - Erz(u, i-1, j-1, k))/(u->dx/2.0);
	}
	else if (u->Fy(0, i-1, j, k) < 0) {
		Ez += u->dx/8.0 * ((-1.0)*u->Fx(iBy, i, j, k) - Erz(u, i-1, j, k))/(u->dx/2.0);
	}
	else {
		/* Average of the two */
		Ez += 0.5*u->dx/8.0*(   ((-1.0)*u->Fx(iBy, i, j-1, k) - Erz(u, i-1, j-1, k))/(u->dx/2.0)   +   ((-1.0)*u->Fx(iBy, i, j, k) - Erz(u, i-1, j, k))/(u->dx/2.0)   );
	}

	/* Add correct term for dEz/dx (i-1/4, j-1/2, k) */
	if (u->Fy(0, i, j, k) > 0) {
		Ez += -u->dx/8.0 * (Erz(u, i, j-1, k) - (-1.0)*u->Fx(iBy, i, j-1, k))/(u->dx/2.0);
	}
	else if (u->Fy(0, i, j, k) < 0) {
		Ez += -u->dx/8.0 * (Erz(u, i, j, k) - (-1.0)*u->Fx(iBy, i, j, k))/(u->dx/2.0);
	}
	else {
		/* Average of the two */
		Ez += -0.5*u->dx/8.0 *(   (Erz(u, i, j-1, k) - (-1.0)*u->Fx(iBy, i, j-1, k))/(u->dx/2.0)   +   (Erz(u, i, j, k) - (-1.0)*u->Fx(iBy, i, j, k))/(u->dx/2.0)   );
	}

	return Ez;
}




bool getEMFs(Arrays *u, ld *E, int i, int j, int k) {

	if (u->Ny > 1 && u->Nz > 1) E[0] = getEx(u, i, j, k);
	if (u->Nx > 1 && u->Nz > 1) E[1] = getEy(u, i, j, k);
	if (u->Nx > 1 && u->Ny > 1) E[2] = getEz(u, i, j, k);

	return true;
}

