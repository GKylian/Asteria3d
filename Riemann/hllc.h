#pragma once
#include "../Defs/Arrays.h"
#include "../prob.h"

bool checkS(ld *wL, ld *wR, ld SL, ld SR, ld gamma);

bool hllcG(ld *uL, ld *wL, ld *uR, ld *wR, ld *FL, ld *FR, ld *F, ld SL, ld SR);
bool hllcL(ld *uL, ld *wL, ld *uR, ld *wR, ld *FL, ld *FR, ld *F, ld SL, ld SR);



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
		if (!checkS(wL, wR, SL, SR, gamma)) return false;
		F(wL, uL, FL); F(wR, uR, FR);
		if (null(*std::max_element(FL, FL+NVAL)) && null(*std::min_element(FL, FL+NVAL))) {
			std::cout << "ERROR:::hllc.h::getFluxes(x): All the components of FL are null !" << std::endl; return false;
		}
		if (null(*std::max_element(FR, FR+NVAL)) && null(*std::min_element(FR, FR+NVAL))) {
			std::cout << "ERROR:::hllc.h::getFluxes(x): All the components of FR are null !" << std::endl; return false;
		}

		/* ----- 4. Compute the final fluxes with the riemann solver ----- */
		if (!hllcL(uL, wL, uR, wR, FL, FR, Fs, SL, SR)) {
			std::cout << "ERROR:::hllc.h::getFluxes(x): HLLE did not work" << std::endl;
			return false;
		}



		if (null(*std::max_element(Fs, Fs+NVAL)) && null(*std::min_element(Fs, Fs+NVAL))) {
			std::cout << "ERROR:::hllc.h::getFluxes(x): All the components of Fs are null !" << std::endl;
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
		if (!checkS(wL, wR, SL, SR, gamma)) return false;
		F(wL, uL, FL); F(wR, uR, FR);
		
		/* ----- 4. Compute the final fluxes with the riemann solver and swap them back ----- */
		if (!hllcL(uL, wL, uR, wR, FL, FR, Fs, SL, SR)) {
			std::cout << "ERROR:::hllc.h::getFluxes(y): HLLE did not work" << std::endl;
			return false;
		}
		/* Order: (vy, vx, vz) -> (vx, vy, vz) */
		swap(Fs[1], Fs[2]);
#ifdef MHD
		swap(Fs[iBx], Fs[iBy]);
#endif // MHD
		
		//TODO:cout << "\t" << FR[iBy] << endl;
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
		if (!checkS(wL, wR, SL, SR, gamma)) return false;
		F(wL, uL, FL); F(wR, uR, FR);

		/* ----- 4. Compute the final fluxes with the riemann solver and swap them back ----- */
		if (!hllcL(uL, wL, uR, wR, FL, FR, Fs, SL, SR)) {
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



bool hllcG(ld *uL, ld *wL, ld *uR, ld *wR, ld *FL, ld *FR, ld *F, ld SL, ld SR) {

	ld wLs[NVAL] = { 0 }; ld wRs[NVAL] = { 0 };
#ifdef MHD
	if (wL[5] == 0.0) {
		ld pL = wL[4] + 0.5*SQ(wL[5])+SQ(wL[6])+SQ(wL[7]);
		ld pR = wR[4] + 0.5*SQ(wR[5])+SQ(wR[6])+SQ(wR[7]);

#else
		ld pL = wL[4];
		ld pR = wR[4];
#endif // MHD

		ld SM = ( (SR-wR[1])*wR[0]*wR[1] - (SL-wL[1])*wL[0]*wL[1] - pR + pL ) / ((SR-wR[1])*wR[0] - (SL-wL[1])*wL[0]);

		ld pLs = pL + wL[0]*(SL-wL[1])*(SM-wL[1]);   ld pRs = pR + wR[0]*(SR-wR[1])*(SM-wR[1]);
		ld ps = (pLs+pRs)/2.0;

		wLs[0] = wL[0]*(SL-wL[1])/(SL-SM);   wRs[0] = wR[0]*(SR-wR[1])/(SR-SM);
		wLs[1] = SM;   wRs[1] = SM;
		wLs[2] = wL[2];   wRs[2] = wR[2];
		wLs[3] = wL[3];   wRs[3] = wR[3];
		wLs[4] = ((SL-wL[1])*wL[4] - pL*wL[1] + ps*SM) / (SL-SM);   wRs[4] = ((SR-wR[1])*wR[4] - pR*wR[1] + ps*SM) / (SR-SM);

#ifdef MHD
		wLs[6] = wL[6]*(SL-wL[1])/(SL-SM);   wRs[6] = wR[6]*(SR-wR[1])/(SR-SM);
		wLs[7] = wL[7]*(SL-wL[1])/(SL-SM);   wRs[7] = wR[7]*(SR-wR[1])/(SR-SM);
#endif // MHD

#ifdef MHD
	}
#endif // MHD
#ifdef MHD
	else {
		ld pL = wL[4] + SQ(wL[5])+SQ(wL[6])+SQ(wL[7]);
		ld pR = wR[4] + SQ(wR[5])+SQ(wR[6])+SQ(wR[7]);

		ld SM = ((SR-wR[1])*wR[0]*wR[1] - (SL-wL[1])*wL[0]*wL[1] - pR + pL) / ((SR-wR[1])*wR[0] - (SL-wL[1])*wL[0]);

		ld pLs = pL + wL[0]*(SL-wL[1])*(SM-wL[1]);   ld pRs = pR + wR[0]*(SR-wR[1])*(SM-wR[1]);
		ld ps = (pLs+pRs)/2.0;

		wLs[0] = wL[0]*(SL-wL[1])/(SL-SM);   wRs[0] = wR[0]*(SR-wR[1])/(SR-SM);
		wLs[1] = SM;   wRs[1] = SM;
		wLs[6] = wL[6];   wRs[6] = wR[6];
		wLs[7] = wL[7];   wRs[7] = wR[7];

		/*wLs[2] = wL[2]*wL[5]*()*/

		wLs[4] = ((SL-wL[1])*wL[4] - pL*wL[1] + ps*SM) / (SL-SM);   wRs[4] = ((SR-wR[1])*wR[4] - pR*wR[1] + ps*SM) / (SR-SM);



	}	
#endif // MHD



	if (SL > 0.0) {
		for (int n = 0; n < NVAL; n++) 
			F[n] = FL[n];
	}
	if (SL <= 0.0 && SR >= 0.0) {
		for (int n = 0; n < NVAL; n++)
			F[n] = (SR*FL[n] - SL*FR[n] + SR*SL*(uR[n]-uL[n])) / (SR-SL); //F[i] = (SR*FL[i] - SL*FR[i] + SR*SL*(uR[i]-uL[i])) / (SR-SL);
	}
	if (SR < 0.0) {
		for (int n = 0; n < NVAL; n++)
			F[n] = FR[n];
	}
	
	return true;
}




bool checkS(ld *wL, ld *wR, ld SL, ld SR, ld gamma) {
	ld aL = sqrtl(gamma*wL[4]/wL[0]);   ld aR = sqrtl(gamma*wR[4]/wR[0]);

	if (SL >= wL[1]-sqrtl((gamma-1)/(2*gamma))*aL) {
		cout << "\thllc.h:::checkS:: The positivity condition for SL is not satisfied: " << SL << " >= " << wL[1]-sqrtl((gamma-1)/(2*gamma))*aL << endl;
		return false;
	}
	if (SR <= wR[1]+sqrtl((gamma-1)/(2*gamma))*aR) {
		cout << "\thllc.h:::checkS:: The positivity condition for SL is not satisfied: " << SL << " >= " << wL[1]-sqrtl((gamma-1)/(2*gamma))*aL << endl;
		return false;
	}

	return true;
}


bool hllcL(ld *uL, ld *wL, ld *uR, ld *wR, ld *FL, ld *FR, ld *F, ld SL, ld SR) {

	ld uLs[NVAL] = { 0 }; ld uRs[NVAL] = { 0 };

	if (SL > 0.0) {
		for (int n = 0; n < NVAL; n++)
			F[n] = FL[n];
		return true;
	}
	if (SR < 0.0) {
		for (int n = 0; n < NVAL; n++)
			F[n] = FR[n];
		return true;
	}
	
#ifdef MHD
	ld pL = wL[4] + 0.5*(SQ(wL[5])+SQ(wL[6])+SQ(wL[7]));
	ld pR = wR[4] + 0.5*(SQ(wR[5])+SQ(wR[6])+SQ(wR[7]));

	ld SM = ((SR-wR[1])*wR[0]*wR[1] - (SL-wL[1])*wL[0]*wL[1] - pR + pL) / ((SR-wR[1])*wR[0] - (SL-wL[1])*wL[0]);

	uLs[0] = uL[0]*(SL-wL[1])/(SL-SM);   uRs[0] = uR[0]*(SR-wR[1])/(SR-SM);
	uLs[1] = uLs[0]*SM;   uRs[1] = uRs[0]*SM;

	// Compute all the HLL averages:
	
	ld Bx = wR[5]; /* (SR*wR[5] - SL*wL[5])/(SR-SL); -> wR[5]=wL[5] */     uLs[5] = Bx;   uRs[5] = Bx;
	ld By = (SR*wR[6] - SL*wL[6] - (FR[6]-FL[6]))/(SR-SL); uLs[6] = By; uRs[6] = By;
	ld Bz = (SR*wR[7] - SL*wL[7] - (FR[7]-FL[7]))/(SR-SL); uLs[7] = Bz; uRs[7] = Bz;
	if (wR[5] != wL[5]) cout << "\thllc.h::hllcL:: BnL != BnL" << endl;

	ld rho = (SR*wR[0] - SL*wL[0] - (FR[0]-FL[0]))/(SR-SL);
	ld vx = (SR*uR[1] - SL*uL[1] - (FR[1]-FL[1]))/(SR-SL) / rho;
	ld vy = (SR*uR[2] - SL*uL[2] - (FR[2]-FL[2]))/(SR-SL) / rho;
	ld vz = (SR*uR[3] - SL*uL[3] - (FR[3]-FL[3]))/(SR-SL) / rho;


	ld pLs = wL[0]*(SL-wL[1])*(SM-wL[1])+pL-SQ(wL[5])+SQ(Bx);
	ld pRs = wR[0]*(SR-wR[1])*(SM-wR[1])+pR-SQ(wR[5])+SQ(Bx);
	ld ps = (pLs+pRs)/2.0;


	//TODO: What values do My and Mz get then ?
	uLs[2] = uL[2]*(SL-wL[1])/(SL-SM) - (Bx*uLs[6]-wL[5]*wL[6])/(SL-SM);   uRs[2] = uR[2]*(SR-wR[1])/(SR-SM) - (Bx*uRs[6]-wR[5]*wR[6])/(SR-SM);
	uLs[3] = uL[3]*(SL-wL[1])/(SL-SM) - (Bx*uLs[7]-wL[5]*wL[7])/(SL-SM);   uRs[3] = uR[3]*(SR-wR[1])/(SR-SM) - (Bx*uRs[7]-wR[5]*wR[7])/(SR-SM);
	if (isnan(uLs[2])) {
		cout << "\tBx = " << Bx << ", uLs[6] = " << uLs[6] << ", uLs[7] = " << uLs[7] << endl;
	}

	
	ld BdotuHLL = Bx*vx+By*vy+Bz*vz;

	ld BdotuL = wL[1]*wL[5] + wL[2]*wL[6] + wL[3]*wL[7];
	uLs[4] = uL[4]*(SL-wL[1])/(SL-SM) + ( (ps*SM-pL*wL[1]) - (Bx*BdotuHLL-wL[5]*BdotuL) ) / (SL-SM);

	ld BdotuR = wR[1]*wR[5] + wR[2]*wR[6] + wR[3]*wR[7];
	uRs[4] = uR[4]*(SR-wR[1])/(SR-SM) + ( (ps*SM-pR*wR[1]) - (Bx*BdotuHLL-wR[5]*BdotuR) ) / (SR-SM);

#else

	ld pL = wL[4];
	ld pR = wR[4];

	ld SM = ((SR-wR[1])*wR[0]*wR[1] - (SL-wL[1])*wL[0]*wL[1] - pR + pL) / ((SR-wR[1])*wR[0] - (SL-wL[1])*wL[0]);

	uLs[0] = uL[0]*(SL-wL[1])/(SL-SM);   uRs[0] = uR[0]*(SR-wR[1])/(SR-SM);
	uLs[1] = uLs[0]*SM;   uRs[1] = uRs[0]*SM;

	ld pLs = pL+wL[0]*(SL-wL[1])*(SM-wL[1]);
	ld pRs = pR+wR[0]*(SR-wR[1])*(SM-wR[1]);
	ld ps = (pLs+pRs)/2.0;
	/*if (fabsl(pLs-pRs) > pLs/10.0) {
		cout << "\tpLs and pRs: " << pLs << ", " << pRs << endl;
	}*/

	uLs[2] = wL[2]*uLs[0];   uRs[2] = wR[2]*uRs[0];
	uLs[3] = wL[3]*uLs[0];   uRs[3] = wR[3]*uRs[0];

	uLs[4] = uL[4]*(SL-wL[1])/(SL-SM) + (ps*SM-pL*wL[1])/(SL-SM);
	uRs[4] = uR[4]*(SR-wR[1])/(SR-SM) + (ps*SM-pR*wR[1])/(SR-SM);

	/*for (int n = 0; n < NVAL; n++) {
		if (!(uLs[n] <= fminl(uL[n], uR[n]) && uLs[n] >= fmaxl(uL[n], uR[n]))) cout << "\tuLs[" << n << "] is not between uL and uR !" << endl;
		if (!(uRs[n] <= fminl(uL[n], uR[n]) && uRs[n] >= fmaxl(uL[n], uR[n]))) cout << "\tuRs[" << n << "] is not between uL and uR !" << endl;
	}*/

#endif // MHD


	
	if (SL <= 0 && SM >= 0) {
		for (int n = 0; n < NVAL; n++)
			F[n] = FL[n] + SL*(uLs[n]-uL[n]);
		return true;
	}
	if (SR >= 0 && SM <= 0) {
		for (int n = 0; n < NVAL; n++)
			F[n] = FR[n] + SR*(uRs[n]-uR[n]);
		return true;
	}
	
	else {
		cout << "\thllc.h::hllcL:: No SL/SR condition was met !" << endl;
	}

	for (int n = 0; n < NVAL; n++) {
		if (isnan(F[n])) {
			cout << "\tF[" << n << "] is NaN !!!" << endl;
			cout << "\tSL = " << SL << ", SM = " << SM << ", SR = " << SR << endl;
			cout << "\tFL = " << FL[n] << ", FR = " << FR[n] << ", uLs = " << uLs[n] << ", " << uRs[n] << endl;
		}
	}
	
	return true;
}





