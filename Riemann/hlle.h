#pragma once
#include "../Defs/Arrays.h"
#include "../prob.h"

bool hlle(ld *uL, ld *wL, ld *uR, ld *wR, ld *FL, ld *FR, ld *F, ld SL, ld SR);

bool getFluxes(Arrays *u, int i, int j, int k) {
	int iBx = NVAL-3; int iBy = NVAL-2; int iBz = NVAL-1;
	ld gamma = u->gamma;

	ld cL[WV] = { 0 }; ld cR[WV] = { 0 };
	int iR = i*2; int iL = i*2-1; int jR = j*2; int jL = j*2-1; int kR = k*2; int kL = k*2-1;
	ld FL[NVAL] = { 0 }; ld FR[NVAL] = { 0 }; ld Fs[NVAL] = { 0 };

	if (u->Nx > 1) {
		 
		/* ----- 0. Load values into arrays ----- */
		ld uL[NVAL] = { 0 }; ld uR[NVAL] = { 0 }; ld wL[NVAL] = { 0 }; ld wR[NVAL] = { 0 };
		for (int n = 0; n < NVAL; n++) {    uL[n] = u->ix(n, iL, j, k); uR[n] = u->ix(n, iR, j, k);    }
		/* Order: (rho, Mn, M1, M2, E, Bn, B1, B2) */

		if (!toPrimitive(uL, wL, gamma)) { std::cout << "ERROR::: hlle.h::getFluxes:: Could not transform uL to wL." << std::endl; return false; }
		if (!toPrimitive(uR, wR, gamma)) { std::cout << "ERROR::: hlle.h::getFluxes:: Could not transform uR to wR." << std::endl; return false; }
		//std::cout << "Bx(i, j, k) = " << u->uC(5, i, j, k) << "y " << "By(i, j, k) = " << u->uC(6, i, j, k) << "Bz(i, j, k) = " << u->uC(7, i, j, k) << std::endl;



		/* ----- 1. Get wavespeeds and max/min Roe eigenvalues ----- */
		ld CL = maxWavespeed(wL, gamma);   ld CR = maxWavespeed(wR, gamma);
		ld l[2] = { 0 }; minmaxRoeEigenvalue(uL, wL, uR, wR, gamma, l);


		/* ----- 2. Compute left and right signal speeds ----- */
		/*ld SL = fminl(wL[1], wR[1]) - fmaxl(CL, CR);
		ld SR = fmaxl(wL[1], wR[1]) + fmaxl(CL, CR);*/
		ld SL = fminl(l[0], wL[1]-CL);   ld SR = fmaxl(l[1], wR[1]+CR);


		//ld bmin = fminl(SL, 0); ld bplus = fmaxl(SR, 0);

		/* ----- 3. Compute left/right interfaces fluxes ----- */
		F(wL, uL, FL); F(wR, uR, FR);
		if (null(*std::max_element(FL, FL+NVAL)) && null(*std::min_element(FL, FL+NVAL))) {
			std::cout << "ERROR:::hlle.h::getFluxes(x): All the components of FL are null !" << std::endl; return false;
		}
		if (null(*std::max_element(FR, FR+NVAL)) && null(*std::min_element(FR, FR+NVAL))) {
			std::cout << "ERROR:::hlle.h::getFluxes(x): All the components of FR are null !" << std::endl; return false;
		}

		/* ----- 4. Compute the final fluxes with the riemann solver ----- */
		if (!hlle(uL, wL, uR, wR, FL, FR, Fs, SL, SR)) {
			std::cout << "ERROR:::hlle.h::getFluxes(x): HLLE did not work" << std::endl;
			return false;
		}



		if (null(*std::max_element(Fs, Fs+NVAL)) && null(*std::min_element(Fs, Fs+NVAL))) {
			std::cout << "ERROR:::hlle.h::getFluxes(x): All the components of Fs are null !" << std::endl;
			std::cout << "---> SL = " << SL << ", SR = " << SR << std::endl;
			std::cout << "---> CL = " << CL << ", CR = " << CR << std::endl;
			std::cout << "---> wL[1] = " << wL[1] << ", wR[1] = " << wR[1] << std::endl;
			return false;
		}

		/* ----- 5. Transfer the fluxes to the arrays ----- */
		for (int n = 0; n < NVAL; n++)
			u->Fx(n, i, j, k) = Fs[n];
	}


	if (u->Ny > 1) {

		/* ----- 0. Load values into arrays ----- */
		ld uL[NVAL] = { 0 }; ld uR[NVAL] = { 0 }; ld wL[NVAL] = { 0 }; ld wR[NVAL] = { 0 };
		for (int n = 0; n < NVAL; n++) { uL[n] = u->iy(n, i, jL, k); uR[n] = u->iy(n, i, jR, k); }

		/* Order: (vx, vy, vz) -> (vy, vx, vz) */
		swap(uL[1], uL[2]); swap(uR[1], uR[2]);
#ifdef MHD
		swap(uL[iBx], uL[iBy]); swap(uR[iBx], uR[iBy]);
#endif // MHD


		if (!toPrimitive(uL, wL, gamma)) { std::cout << "ERROR:::hlle.h::getFluxes:: Could not transform uL to wL." << std::endl; return false; }
		if (!toPrimitive(uR, wR, gamma)) { std::cout << "ERROR:::hlle.h::getFluxes:: Could not transform uR to wR." << std::endl; return false; }


		/* ----- 1. Get wavespeeds and max/min Roe eigenvalues ----- */
		ld CL = maxWavespeed(wL, gamma);   ld CR = maxWavespeed(wR, gamma);
		ld l[2] = { 0 }; minmaxRoeEigenvalue(uL, wL, uR, wR, gamma, l);


		/* ----- 2. Compute left and right signal speeds ----- */
		/*ld SL = fminl(wL[1], wR[1]) - fmaxl(CL, CR);
		ld SR = fmaxl(wL[1], wR[1]) + fmaxl(CL, CR);*/
		ld SL = fminl(l[0], wL[1]-CL);   ld SR = fmaxl(l[1], wR[1]+CR);


		/* ----- 3. Compute left/right interfaces fluxes ----- */
		F(wL, uL, FL); F(wR, uR, FR);
		
		/* ----- 4. Compute the final fluxes with the riemann solver and swap them back ----- */
		if (!hlle(uL, wL, uR, wR, FL, FR, Fs, SL, SR)) {
			std::cout << "ERROR:::hlle.h::getFluxes(y): HLLE did not work" << std::endl;
			return false;
		}
		/* Order: (vy, vx, vz) -> (vx, vy, vz) */
		swap(Fs[1], Fs[2]);
#ifdef MHD
		swap(Fs[iBx], Fs[iBy]);
#endif // MHD

		//TODO:cout << "\t" << Fs[iBy] << endl;
		/* ----- 5. Transfer the fluxes to the arrays ----- */
		for (int n = 0; n < NVAL; n++)
			u->Fy(n, i, j, k) = Fs[n];

	}


	if (u->Nz > 1) {

		/* ----- 0. Load values into arrays ----- */
		ld uL[NVAL] = { 0 }; ld uR[NVAL] = { 0 }; ld wL[NVAL] = { 0 }; ld wR[NVAL] = { 0 };
		for (int n = 0; n < NVAL; n++) { uL[n] = u->iz(n, i, j, kL); uR[n] = u->iz(n, i, j, kR); }

		/* Order: (vx, vy, vz) -> (vz, vx, vy) */
		swap(uL[1], uL[3]); swap(uL[2], uL[3]);   swap(uR[1], uR[3]); swap(uR[2], uR[3]);
#ifdef MHD
		swap(uL[iBx], uL[iBz]); swap(uL[iBy], uL[iBz]);   swap(uR[iBx], uR[iBz]); swap(uR[iBy], uR[iBz]);
#endif // MHD


		if (!toPrimitive(uL, wL, gamma)) { std::cout << "ERROR:::hlle.h::getFluxes:: Could not transform uL to wL." << std::endl; return false; }
		if (!toPrimitive(uR, wR, gamma)) { std::cout << "ERROR:::hlle.h::getFluxes:: Could not transform uR to wR." << std::endl; return false; }


		/* ----- 1. Get wavespeeds and max/min Roe eigenvalues ----- */
		ld CL = maxWavespeed(wL, gamma);   ld CR = maxWavespeed(wR, gamma);
		ld l[2] = { 0 }; minmaxRoeEigenvalue(uL, wL, uR, wR, gamma, l);


		/* ----- 2. Compute left and right signal speeds ----- */
		/*ld SL = fminl(wL[1], wR[1]) - fmaxl(CL, CR);
		ld SR = fmaxl(wL[1], wR[1]) + fmaxl(CL, CR);*/
		ld SL = fminl(l[0], wL[1]-CL);   ld SR = fmaxl(l[1], wR[1]+CR);


		/* ----- 3. Compute left/right interfaces fluxes ----- */
		F(wL, uL, FL); F(wR, uR, FR);

		/* ----- 4. Compute the final fluxes with the riemann solver and swap them back ----- */
		if (!hlle(uL, wL, uR, wR, FL, FR, Fs, SL, SR)) {
			std::cout << "ERROR:::hlle.h::getFluxes(z): HLLE did not work" << std::endl;
			return false;
		}
		/* Order: (vz, vx, vy) -> (vx, vy, vz) */
		swap(Fs[2], Fs[3]); swap(Fs[1], Fs[3]);
#ifdef MHD
		swap(Fs[iBy], Fs[iBz]); swap(Fs[iBx], Fs[iBz]);
#endif // MHD


		/* ----- 5. Transfer the fluxes to the arrays ----- */
		for (int n = 0; n < NVAL; n++)
			u->Fz(n, i, j, k) = Fs[n];

	}

	return true;

}



bool hlle(ld *uL, ld *wL, ld *uR, ld *wR, ld *FL, ld *FR, ld *F, ld SL, ld SR) {

	if (SL > 0.0) {
		for (int n = 0; n < NVAL; n++) 
			F[n] = FL[n];
	}
	if (SL <= 0.0 && SR >= 0.0) {
		for (int n = 0; n < NVAL; n++)
			F[n] = (SR*FL[n] - SL*FR[n] + SR*SL*(uR[n]-uL[n])) / (SR-SL);
	}
	if (SR < 0.0) {
		for (int n = 0; n < NVAL; n++)
			F[n] = FR[n];
	}
	
	return true;
}
