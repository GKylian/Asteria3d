#pragma once
#include "../Defs/Arrays.h"
#include "../prob.h"


bool checkS(ld *wL, ld *wR, ld SL, ld SR, ld gamma);

bool hlld(ld *uL, ld *wL, ld *uR, ld *wR, ld *FL, ld *FR, ld *F, ld SL, ld SR);



bool getFluxes(Arrays *u, int i, int j, int k) {
	int iBx = NVAL-3; int iBy = NVAL-2; int iBz = NVAL-1;
	ld gamma = u->gamma;

	ld cL[WV] = { 0 }; ld cR[WV] = { 0 };
	int iR = i*2; int iL = i*2-1; int jR = j*2; int jL = j*2-1; int kR = k*2; int kL = k*2-1;
	ld FL[NVAL] = { 0 }; ld FR[NVAL] = { 0 }; ld Fs[NVAL] = { 0 };

	if (u->Nx > 1) {

		/* ----- 0. Load values into arrays ----- */
		ld uL[NVAL] = { 0 }; ld uR[NVAL] = { 0 }; ld wL[NVAL] = { 0 }; ld wR[NVAL] = { 0 };
		for (int n = 0; n < NVAL; n++) { uL[n] = u->ix(n, iL, j, k); uR[n] = u->ix(n, iR, j, k); }
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
		if (FL[iBx] != 0.0 || FR[iBx] != 0.0) {
			std::cout << "\thlld.h:: FL[Bx] or FR[Bx] is not null !" << std::endl;
		}

		/* ----- 4. Compute the final fluxes with the riemann solver ----- */
		if (!hlld(uL, wL, uR, wR, FL, FR, Fs, SL, SR)) {
			std::cout << "ERROR:::hllc.h::getFluxes(x): HLLE did not work" << std::endl;
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
		if (!hlld(uL, wL, uR, wR, FL, FR, Fs, SL, SR)) {
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
		if (!hlld(uL, wL, uR, wR, FL, FR, Fs, SL, SR)) {
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








bool checkS(ld *wL, ld *wR, ld SL, ld SR, ld gamma) {
	//return true;

	ld CfL = maxWavespeed(wL, gamma); ld CfR = maxWavespeed(wR, gamma);

	//ld aL = sqrtl(gamma*wL[4]/wL[0]);   ld aR = sqrtl(gamma*wR[4]/wR[0]);

	if (SL >= wL[1]-sqrtl((gamma-1)/(2*gamma))*CfL) {
		cout << "\thllc.h:::checkS:: The positivity condition for SL is not satisfied: " << SL << " >= " << wL[1]-sqrtl((gamma-1)/(2*gamma))*CfL << endl;
		return false;
	}
	if (SR <= wR[1]+sqrtl((gamma-1)/(2*gamma))*CfR) {
		cout << "\thllc.h:::checkS:: The positivity condition for SL is not satisfied: " << SL << " >= " << wL[1]-sqrtl((gamma-1)/(2*gamma))*CfR << endl;
		return false;
	}

	return true;
}




















bool hlld(ld *uL, ld *wL, ld *uR, ld *wR, ld *FL, ld *FR, ld *F, ld SL, ld SR) {

	/* First set of conditions (no need to compute F* and F**) */
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

	ld uLs[NVAL] = { 0 }; ld uRs[NVAL] = { 0 };
	ld uLss[NVAL] = { 0 }; ld uRss[NVAL] = { 0 };

#ifdef MHD
	ld Bx = wL[5]; ld Bx2 = Bx*Bx;
	ld pL = wL[4] + 0.5*(Bx2+SQ(wL[6])+SQ(wL[7]));
	ld pR = wR[4] + 0.5*(Bx2+SQ(wR[6])+SQ(wR[7]));
	

	ld SM = ((SR-wR[1])*wR[0]*wR[1] - (SL-wL[1])*wL[0]*wL[1] - pR + pL) / ((SR-wR[1])*wR[0] - (SL-wL[1])*wL[0]);

	ld ps = (  (SR-wR[1])*wR[0]*pL - (SL-wL[1])*wL[0]*pR + wL[0]*wR[0]*(SR-wR[1])*(SL-wL[1])*(wR[1]-wL[1])  ) / (  (SR-wR[1])*wR[0] - (SL-wL[1])*wL[0]  );

	uLs[0] = wL[0]*(SL-wL[1])/(SL-SM);   uRs[0] = wR[0]*(SR-wR[1])/(SR-SM);
	uLs[1] = SM;   uRs[1] = SM;
	uLs[5] = Bx;   uRs[5] = Bx;


	ld divL = wL[0]*(SL-wL[1])*(SL-SM) - Bx2;   ld divR = wR[0]*(SR-wR[1])*(SR-SM) - Bx2;
	uLs[2] = wL[2] - Bx*wL[6]*(SM-wL[1])/divL;   uRs[2] = wR[2] - Bx*wR[6]*(SM-wR[1])/divR;
	uLs[3] = wL[3] - Bx*wL[7]*(SM-wL[1])/divL;   uRs[3] = wR[3] - Bx*wR[7]*(SM-wR[1])/divR;


	uLs[6] = wL[6]*( wL[0]*SQ(SL-wL[1]) - Bx2 ) / divL;   uRs[6] = wR[6]*( wR[0]*SQ(SR-wR[1]) - Bx2 ) / divR;
	uLs[7] = wL[7]*( wL[0]*SQ(SL-wL[1]) - Bx2 ) / divL;   uRs[7] = wR[7]*( wR[0]*SQ(SR-wR[1]) - Bx2 ) / divR;

	/* If vy*, vz*, By* or Bz* contains a 0/0, set other values */
	if (divL == 0.0) {
		uLs[2] = wL[2];   uLs[3] = wL[3];   
		uLs[6] = 0.0;   uLs[7] = 0.0;
	}
	if (divR == 0.0) {
		uRs[2] = wR[2];   uRs[3] = wR[3];
		uRs[6] = 0.0;   uRs[7] = 0.0;
	}

	uLs[4] = (  (SL-wL[1])*uL[4] - pL*wL[1] + ps*SM + Bx*( (wL[1]*wL[5]+wL[2]*wL[6]+wL[3]*wL[7]) - (uLs[1]*uLs[5]+uLs[2]*uLs[6]+uLs[3]*uLs[7]) )  ) / (SL-SM);
	uRs[4] = (  (SR-wR[1])*uR[4] - pR*wR[1] + ps*SM + Bx*( (wR[1]*wR[5]+wR[2]*wR[6]+wR[3]*wR[7]) - (uRs[1]*uRs[5]+uRs[2]*uRs[6]+uRs[3]*uRs[7]) )  ) / (SR-SM);


	ld SLs = SM - fabsl(Bx)/sqrt(uLs[0]);   ld SRs = SM + fabsl(Bx)/sqrtl(uRs[0]);

	/* HLLC-like conditions for intermediate states */
	if (SL <= 0 && SLs >= 0) { // SLs vs SM
		/* v to M */
		uLs[1] *= uLs[0];   uLs[2] *= uLs[0];   uLs[3] *= uLs[0];

		for (int n = 0; n < NVAL; n++)
			F[n] = FL[n] + SL*(uLs[n]-uL[n]);

		if (F[5] != 0.0) std::cout << "\t1. F[5] != 0.0" << std::endl;
		return true;
	}
	if (SR >= 0 && SRs <= 0) { // SRs vs SM
		/* v to M */
		uRs[1] *= uRs[0];   uRs[2] *= uRs[0];   uRs[3] *= uRs[0];

		for (int n = 0; n < NVAL; n++)
			F[n] = FR[n] + SR*(uRs[n]-uR[n]);
		if (F[5] != 0.0) std::cout << "\t2. F[5] != 0.0" << std::endl;
		return true;
	}


	//std::cout << "nope" << std::endl;
	if (Bx == 0.0) {
		std::cout << "\thlld.d: Bx == 0.0, but we went past the hllc conditions !" << std::endl;
	}


	uLss[0] = uLs[0];   uRss[0] = uRs[0];
	uLss[1] = SM;   uRss[1] = SM;
	uLss[5] = Bx;   uRss[5] = Bx;
	ld rtL = sqrtl(uLs[0]);   ld rtR = sqrtl(uRs[0]);
	
	uLss[2] = (rtL*uLs[2] + rtR*uRs[2] + (uRs[6]-uLs[6])*sgn(Bx)) / (rtL+rtR); uRss[2] = uLss[2];
	uLss[3] = (rtL*uLs[3] + rtR*uRs[3] + (uRs[7]-uLs[7])*sgn(Bx)) / (rtL+rtR); uRss[3] = uLss[3];

	uLss[6] = (rtL*uLs[6] + rtR*uRs[6] + rtL*rtR*(uRs[2]-uLs[2])*sgn(Bx)) / (rtL+rtR); uRss[6] = uLss[6];
	uLss[7] = (rtL*uLs[7] + rtR*uRs[7] + rtL*rtR*(uRs[3]-uLs[3])*sgn(Bx)) / (rtL+rtR); uRss[7] = uLss[7];


	uLss[4] = uLs[4] - rtL*( (uLs[1]*uLs[5]+uLs[2]*uLs[6]+uLs[3]*uLs[7]) - (uLss[1]*uLss[5]+uLss[2]*uLss[6]+uLss[3]*uLss[7]) )*sgn(Bx);
	uRss[4] = uRs[4] + rtR*( (uRs[1]*uRs[5]+uRs[2]*uRs[6]+uRs[3]*uRs[7]) - (uRss[1]*uRss[5]+uRss[2]*uRss[6]+uRss[3]*uRss[7]) )*sgn(Bx);


	/* v to M */
	uLs[1] *= uLs[0]; uRs[1] *= uRs[0];   uLs[2] *= uLs[0]; uRs[2] *= uRs[0];   uLs[3] *= uLs[0]; uRs[3] *= uRs[0];
	uLss[1] *= uLss[0]; uRss[1] *= uRss[0];   uLss[2] *= uLss[0]; uRss[2] *= uRss[0];   uLss[3] *= uLss[0]; uRss[3] *= uRss[0];


	if (SM >= 0 && SLs <= 0) {

		for (int n = 0; n < NVAL; n++)
			F[n] = FL[n] + SL*(uLs[n]-uL[n]) + SLs*(uLss[n]-uLs[n]);
		if (F[5] != 0.0) std::cout << "\t3. F[5] != 0.0" << std::endl;
		return true;
	}
	if (SM <= 0 && SRs >= 0) {

		for (int n = 0; n < NVAL; n++)
			F[n] = FR[n] + SR*(uRs[n]-uR[n]) + SRs*(uRss[n]-uRs[n]);
		if (F[5] != 0.0)
			std::cout << "\t4. F[5] != 0.0" << std::endl;
		
			
		return true;
	}
	

	std::cout << "\thlld.h::No condition was met ?" << std::endl;
	return false;



#else

	ld pL = wL[4];
	ld pR = wR[4];
	

	ld SM = ((SR-wR[1])*wR[0]*wR[1] - (SL-wL[1])*wL[0]*wL[1] - pR + pL) / ((SR-wR[1])*wR[0] - (SL-wL[1])*wL[0]);
	ld ps = (  (SR-wR[1])*wR[0]*pL - (SL-wL[1])*wL[0]*pR + wL[0]*wR[0]*(SR-wR[1])*(SL-wL[1])*(wR[1]-wL[1])  ) / (  (SR-wR[1])*wR[0] - (SL-wL[1])*wL[0]  );

	uLs[0] = wL[0]*(SL-wL[1])/(SL-SM);   uRs[0] = wR[0]*(SR-wR[1])/(SR-SM);
	uLs[1] = SM;   uRs[1] = SM;


	ld divL = (wL[0]*(SL-wL[1])*(SL-SM));   ld divR = (wR[0]*(SR-wR[1])*(SR-SM));
	uLs[2] = wL[2];   uRs[2] = wR[2];
	uLs[3] = wL[3];   uRs[3] = wR[3];


	/* If vy*, vz*, By* or Bz* contains a 0/0, set other values */
	if (divL == 0.0) {
		uLs[2] = wL[2];   uLs[3] = wL[3];
	}
	if (divR == 0.0) {
		uRs[2] = wR[2];   uRs[3] = wR[3];
	}

	uLs[4] = (  (SL-wL[1])*uL[4] - pL*wL[1] + ps*SM ) / (SL-SM);
	uRs[4] = (  (SR-wR[1])*uR[4] - pR*wR[1] + ps*SM ) / (SR-SM);

	/* HLLC-like conditions for intermediate states */
	if (SL <= 0 && SM >= 0) { // SLs vs SM
		/* v to M */
		uLs[1] *= uLs[0]; uRs[1] *= uRs[0];   uLs[2] *= uLs[0]; uRs[2] *= uRs[0];   uLs[3] *= uLs[0]; uRs[3] *= uRs[0];

		for (int n = 0; n < NVAL; n++)
			F[n] = FL[n] + SL*(uLs[n]-uL[n]);
		return true;
	}
	if (SR >= 0 && SM <= 0) { // SRs vs SM
		/* v to M */
		uLs[1] *= uLs[0]; uRs[1] *= uRs[0];   uLs[2] *= uLs[0]; uRs[2] *= uRs[0];   uLs[3] *= uLs[0]; uRs[3] *= uRs[0];

		for (int n = 0; n < NVAL; n++)
			F[n] = FR[n] + SR*(uRs[n]-uR[n]);
		return true;
	}

	std::cout << "\thlld.h::No condition was met ?" << std::endl;
	return false;
#endif // MHD


}