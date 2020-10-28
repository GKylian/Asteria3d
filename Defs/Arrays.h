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
		if (x >= Nx || y >= Ny || z >= Nz)
			cout << "ERROR:::Arrays.h::Array::get() Trying to access an out-of-bound element !" << endl;
		return m_array[z*Nx*Ny + y*Nx + x];
	}

	const T &get(int x, int y) const {
		if (x >= Nx || y >= Ny)
			cout << "ERROR:::Arrays.h::Array::get() Trying to access an out-of-bound element !" << endl;
		return m_array[y * Nx + x];
	}


	T &getMaxValue() {
		return *std::max_element(m_array.begin(), m_array.end());
	}

	void setSize(int _Nx, int _Ny, int _Nz) {
		Nx = _Nx; Ny = _Ny; Nz = _Nz;
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

	void printSize() {
		cout << "(" << Nx << ", " << Ny << ", " << Nz << ")";
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
	Array<ld> cons[NVAL]; /*  rho, Mx, My, Mz, E, Bx, By, Bz   OR   rho, Mx, My, Mz, cs, Bx, By, Bz  */
								   /*  rho, Mx, My, Mz, E               OR   rho, Mx, My, Mz, cs              */

	Array<ld> prim[NVAL];
	//tex: $\vec{v} \to \vec{M}, E \to P$, and $\vec{B}$ are stored at cell centers

	//The three interface arrays
	Array<ld> intx[NVAL];
	Array<ld> inty[NVAL];
	Array<ld> intz[NVAL];

	//The three fluxes arrays.
	Array<ld> F_x[NVAL];
	Array<ld> F_y[NVAL];
	Array<ld> F_z[NVAL];

	
	
	int Nx, Ny, Nz = 0; int dim = 0;
	ld dx, dy, dz = 0;
	ld x0, xn, y0, yn, z0, zn = 0.0;
	bounds boundaries[6]; /* xin, xout, yin, yout, zin, zout */
	ld t; ld dt; int s; ld tn = 1.0;
	ld gamma = 1.4;


	//dl/dr: first/last cell of the whole domain	gl/gr: last/first ghost cell	cl/cr: first/last cell of the computation domain
	int i_dl, i_gl, i_cl, i_dr, i_gr, i_cr;
	int j_dl, j_gl, j_cl, j_dr, j_gr, j_cr;
	int k_dl, k_gl, k_cl, k_dr, k_gr, k_cr;

	void initAll(int nx, int ny, int nz) {
		Nx = nx; Ny = ny; Nz = nz; s = 0;
		for (int i = 0; i < NVAL; i++)
		{
			//if (i == 5 || i == 6) continue;
			cons[i].setSize(Nx+(Nx>1)*(2*NGHOST+(i==5)), Ny+(Ny>1)*(2*NGHOST+(i==6)), Nz+(Nz>1)*(2*NGHOST+(i==7)));
			//cout << "cons[" << i << "] array size: "; cons[i].printSize(); cout << endl;
			prim[i].setSize(Nx+(Nx>1)*2*NGHOST, Ny+(Ny>1)*2*NGHOST, Nz+(Nz>1)*2*NGHOST);
			//cout << "prim[" << i << "] array size: "; prim[i].printSize(); cout << endl;

			intx[i].setSize((Nx>1)*2*(Nx+2*NGHOST)+(Ny<=1), Ny+(Ny>1)*2*NGHOST, Nz+(Nz>1)*2*NGHOST);
			//cout << "intx[" << i << "] array size: "; intx[i].printSize(); cout << endl;
			inty[i].setSize(Nx+(Nx>1)*2*NGHOST, (Ny>1)*2*(Ny+2*NGHOST)+(Ny<=1), Nz+(Nz>1)*2*NGHOST);
			//cout << "inty[" << i << "] array size: "; inty[i].printSize(); cout << endl;
			intz[i].setSize(Nx+(Nx>1)*2*NGHOST, Ny+(Ny>1)*2*NGHOST, (Nz>1)*2*(Nz+2*NGHOST)+(Nz<=1));
			//cout << "intz[" << i << "] array size: "; intz[i].printSize(); cout << endl;

			F_x[i].setSize((Nx>1)*(Nx+2*NGHOST)+1, Ny+(Ny>1)*2*NGHOST, Nz+(Nz>1)*2*NGHOST);
			//cout << "F_x[" << i << "] array size: "; F_x[i].printSize(); cout << endl;
			//cout << "--> should be " << (Nx>1)*(Nx+2*NGHOST)+1 << ", " << Ny+(Ny>1)*2*NGHOST << ", " << Nz+(Nz>1)*2*NGHOST << endl;
			F_y[i].setSize(Nx+(Nx>1)*2*NGHOST, (Ny>1)*(Ny+2*NGHOST)+1, Nz+(Nz>1)*2*NGHOST);
			//cout << "F_y[" << i << "] array size: "; F_y[i].printSize(); cout << endl;
			F_z[i].setSize(Nx+(Nx>1)*2*NGHOST, Ny+(Ny>1)*2*NGHOST, (Nz>1)*(Nz+2*NGHOST)+1);
			//cout << "F_z[" << i << "] array size: "; F_z[i].printSize(); cout << endl;
		}
		setIndexes();
		dim = (nx > 1)+(ny > 1)+(nz > 1);
	}

	void initAll() {
		initAll(Nx, Ny, Nz);
	}

	void seth(ld _dx, ld _dy) {
		dx = _dx; dy = _dy;
	}

	void setRange(ld _x0, ld _xn, ld _y0, ld _yn) {
		x0 = _x0; xn = _xn; y0 = _y0; yn = _yn;
	}

	void seth(ld _dx, ld _dy, ld _dz) {
		dx = _dx; dy = _dy; dz = _dz;
	}

	void seth() {
		dx = (xn-x0)/(Nx-1.0); if (Nx <= 1) dx = 1.0;
		dy = (yn-y0)/(Ny-1.0); if (Ny <= 1) dy = 1.0;
		dz = (zn-z0)/(Nz-1.0); if (Nz <= 1) dz = 1.0;
	}

	void setRange(ld _x0, ld _xn, ld _y0, ld _yn, ld _z0, ld _zn) {
		x0 = _x0; xn = _xn; y0 = _y0; yn = _yn; z0 = _z0; zn = _zn;
	}

	ld &uC(int i, int xi, int yi) {
		return cons[i].get(xi, yi);
	}
	ld &uC(int i, int xi, int yi, int zi) {
		return cons[i].get(xi, yi, zi);
	}
	Array<ld> &getC(int i) {
		return cons[i];
	}
	ld &uP(int i, int xi, int yi) {
		return prim[i].get(xi, yi);
	}
	ld &uP(int i, int xi, int yi, int zi) {
		return prim[i].get(xi, yi, zi);
	}
	Array<ld> &getP(int i) {
		return prim[i];
	}

	ld &ix(int i, int xi, int yi, int zi) {
		return intx[i].get(xi, yi, zi);
	}
	ld &iy(int i, int xi, int yi, int zi) {
		return inty[i].get(xi, yi, zi);
	}
	ld &iz(int i, int xi, int yi, int zi) {
		return intz[i].get(xi, yi, zi);
	}

	ld &Fx(int i, int xi, int yi, int zi) {
		if (xi > i_dr || yi > j_dr || zi > k_dr) cout << "ERROR:::Arrays.h::Fx:: Out of bound: " << xi << ", " << yi << ", " << zi << endl;
		return F_x[i].get(xi, yi, zi);
	}
	ld &Fy(int i, int xi, int yi, int zi) {
		if (xi > i_dr || yi > j_dr || zi > k_dr) cout << "ERROR:::Arrays.h::Fy:: Out of bound: " << xi << ", " << yi << ", " << zi << endl;
		return F_y[i].get(xi, yi, zi);
	}
	ld &Fz(int i, int xi, int yi, int zi) {
		if (xi > i_dr || yi > j_dr || zi > k_dr) cout << "ERROR:::Arrays.h::Fz:: Out of bound: " << xi << ", " << yi << ", " << zi << endl;
		return F_z[i].get(xi, yi, zi);
	}


	void pos(ld *x, int i, int j, int k) {
		x[0] = x0+(i-i_cl)*dx;
		x[1] = y0+(j-j_cl)*dy;
		x[2] = z0+(k-k_cl)*dz;
	}

	void setIndexes() {
		//Array: first of all; last of the ghost; first of domain; last of domain (except Bx); first of ghost; last of all;
		if (Nx > 1) {
			i_dl = 0;				i_gl = NGHOST-1;	i_cl = NGHOST;
			i_cr = NGHOST+Nx-1;		i_gr = NGHOST+Nx;	i_dr = 2*NGHOST+Nx-1;
		}
		else {
			i_dl = 0; i_gl = 0; i_cl = 0; i_dr = 0; i_gr = 0; i_cr = 0;
		}

		if (Ny > 1) {
			j_dl = 0;				j_gl = NGHOST-1;	j_cl = NGHOST;
			j_cr = NGHOST+Ny-1;		j_gr = NGHOST+Ny;	j_dr = 2*NGHOST+Ny-1;
		}
		else {
			j_dl = 0; j_gl = 0; j_cl = 0; j_dr = 0; j_gr = 0; j_cr = 0;
		}

		if (Nz > 1) {
			k_dl = 0;				k_gl = NGHOST-1;	k_cl = NGHOST;
			k_cr = NGHOST+Nz-1;		k_gr = NGHOST+Nz;	k_dr = 2*NGHOST+Nz-1;
		}
		else {
			k_dl = 0; k_gl = 0; k_cl = 0; k_dr = 0; k_gr = 0; k_cr = 0;
		}

	}
};