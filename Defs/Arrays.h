#pragma once

#include <vector>
#include <algorithm>
#include <iostream>

#include "Definitions.h"
using namespace std;

template<typename T>
class Array
{

public:
	T &get(int x, int y, int z) {
		return m_array[z*Nx*Ny + y*Nx + x];
	}

	T &get(int x, int y) {
		return m_array[y * Nx + x];
	}


	const T &get(int x, int y, int z) const {
		return m_array[z*Nx*Ny + y*Nx + x];
	}

	const T &get(int x, int y) const {
		return m_array[y * Nx + x];
	}


	T &getMaxValue() {
		return *std::max_element(m_array.begin(), m_array.end());
	}

	void setSize(int _Nx, int _Ny, int _Nz) {
		Nx = _Nx; int Ny = _Ny; int Nz = _Nz;
		m_array.resize(Nx * Ny * Nz);
		//if (_Nz == 1) dim = 2; else dim = 3;
		dim = 3;

		for (size_t i = 0; i < m_array.size(); ++i) {
			m_array[i] = 0.0;
		}
	}

	void setSize(int _Nx, int _Ny) {
		Nx = _Nx; int Ny = _Ny; int Nz = 1;
		m_array.resize(Nx * Ny * Nz);
		dim = 2;
	}

	std::vector<T> &getVector() {
		return m_array;
	}

	void deleteData() {
		std::vector<T>().swap(m_array);
	}

private:
	std::vector<T> m_array;
	int Nx = 1, Ny = 1, Nz = 1;
	int dim = 1;
};



struct Arrays
{
	//rho, Mx, My, Mz, E, Bx, By, Bz
	Array<long double> cons[NVAL]; /*  rho, Mx, My, Mz, E, Bx, By, Bz   OR   rho, Mx, My, Mz, cs, Bx, By, Bz  */
								   /*  rho, Mx, My, Mz, E               OR   rho, Mx, My, Mz, cs              */

	Array<long double> prim[NVAL];
	//tex: $\vec{v} \to \vec{M}, E \to P$, and $\vec{B}$ are stored at cell centers

	//The three interface arrays
	Array<long double> intx[NVAL];
	Array<long double> inty[NVAL];
	Array<long double> intz[NVAL];

	
	
	int Nx, Ny, Nz = 0;
	long double dx, dy, dz = 0;
	long double x0, xn, y0, yn, z0, zn = 0.0;
	bounds boundaries[6]; /* xin, xout, yin, yout, zin, zout */
	long double t; long double dt; 


	//dl/dr: first/last cell of the whole domain	gl/gr: last/first ghost cell	cl/cr: first/last cell of the computation domain
	int i_dl, i_gl, i_cl, i_dr, i_gr, i_cr;
	int j_dl, j_gl, j_cl, j_dr, j_gr, j_cr;
	int k_dl, k_gl, k_cl, k_dr, k_gr, k_cr;

	void initAll(int nx, int ny, int nz) {
		Nx = nx; Ny = ny; Nz = nz;
		for (int i = 0; i < NVAL; i++)
		{
			//if (i == 5 || i == 6) continue;
			cons[i].setSize(nx+2*NGHOST+(i==5), ny+2*NGHOST+(i==6), nz+2*NGHOST+(i==7));
			prim[i].setSize(nx+2*NGHOST, ny+2*NGHOST, nz+2*NGHOST);

			intx[i].setSize(2*(nx+2*NGHOST), ny+2*NGHOST, nz+2*NGHOST);
			inty[i].setSize(nx+2*NGHOST, 2*(ny+2*NGHOST), nz+2*NGHOST);
			intz[i].setSize(nx+2*NGHOST, ny+2*NGHOST, 2*(nz+2*NGHOST));
		}
		setIndexes();
	}

	void seth(long double _dx, long double _dy) {
		dx = _dx; dy = _dy;
	}

	void setRange(long double _x0, long double _xn, long double _y0, long double _yn) {
		x0 = _x0; xn = _xn; y0 = _y0; yn = _yn;
	}

	void seth(long double _dx, long double _dy, long double _dz) {
		dx = _dx; dy = _dy; dz = _dz;
	}

	void setRange(long double _x0, long double _xn, long double _y0, long double _yn, long double _z0, long double _zn) {
		x0 = _x0; xn = _xn; y0 = _y0; yn = _yn; z0 = _z0; zn = _zn;
	}

	long double &uC(int i, int xi, int yi) {
		return cons[i].get(xi, yi);
	}
	long double &uC(int i, int xi, int yi, int zi) {
		return cons[i].get(xi, yi, zi);
	}
	Array<long double> &getC(int i) {
		return cons[i];
	}
	long double &uP(int i, int xi, int yi) {
		return prim[i].get(xi, yi);
	}
	long double &uP(int i, int xi, int yi, int zi) {
		return prim[i].get(xi, yi, zi);
	}
	Array<long double> &getP(int i) {
		return prim[i];
	}
	long double &ix(int i, int xi, int yi, int zi) {
		return intx[i].get(xi, yi, zi);
	}
	long double &iy(int i, int xi, int yi, int zi) {
		return inty[i].get(xi, yi, zi);
	}
	long double &iz(int i, int xi, int yi, int zi) {
		return intz[i].get(xi, yi, zi);
	}

	void setIndexes() {
		//Array: first of all; last of the ghost; first of domain; last of domain (except Bx); first of ghost; last of all;
		if (Nx > 1) {
			i_dl = 0;				i_gl = NGHOST-1;	i_cl = NGHOST;
			i_cr = NGHOST+Nx-1;		i_gr = NGHOST+Nx;	i_dr = 2*NGHOST+Nx-1;
		}
		else {
			i_dl = 0; i_gl = 0; i_cl = 0; i_dr = 1; i_gr = 1; i_cr = 1;
		}

		if (Ny > 1) {
			j_dl = 0;				j_gl = NGHOST-1;	j_cl = NGHOST;
			j_cr = NGHOST+Ny-1;		j_gr = NGHOST+Ny;	j_dr = 2*NGHOST+Ny-1;
		}
		else {
			j_dl = 0; j_gl = 0; j_cl = 0; j_dr = 1; j_gr = 1; j_cr = 1;
		}

		if (Nz > 1) {
			k_dl = 0;				k_gl = NGHOST-1;	k_cl = NGHOST;
			k_cr = NGHOST+Nz-1;		k_gr = NGHOST+Nz;	k_dr = 2*NGHOST+Nz-1;
		}
		else {
			k_dl = 0; k_gl = 0; k_cl = 0; k_dr = 1; k_gr = 1; k_cr = 1;
		}

	}
};