#pragma once

#include <vector>
#include <algorithm>



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
	Array<ld> arr[8];
	int Nx, Ny = 0;
	ld dx, dy = 0;
	ld x0, xn, y0, yn = 0.0;

	void initAll(int nx, int ny) {
		Nx = nx; Ny = ny;
		for (int i = 0; i < 8; i++)
		{
			//if (i == 5 || i == 6) continue;
			arr[i].setSize(nx+(i==5), ny+(i==6));
		}
		
	}

	void seth(ld _dx, ld _dy) {
		dx = _dx; dy = _dy;
	}

	void setRange(ld _x0, ld _xn, ld _y0, ld _yn) {
		x0 = _x0; xn = _xn; y0 = _y0; yn = _yn;
	}

	ld &at(int i, int xi, int yi) {
		return arr[i].get(xi, yi);
	}
	
	

	Array<ld> &get(int i) {
		return arr[i];
	}
};