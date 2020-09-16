#pragma once
#include "MHD.h"
#include <string>
#include <fstream>
#include <sstream>

using namespace std;

bool LoadValues(string path, Arrays *us, int Nx, int Ny) {
	ifstream csv; csv.open(path);
	string line;

	if (!csv.good()) {
		cout << "ERROR::InOut.h:: Could not load csv file at " << path << endl;
		return false;
	}

	
	getline(csv, line); stringstream iss(line);
	//1. Check that we have the same Nx and Ny
	string val; getline(iss, line, ",");
	if (stoi(val) != Nx) {
		cout << "ERROR::InOut.h:: Nx defined in .csv file is not the one expected: " << stoi(val) << endl;
		return false;
	} getline(iss, line, ","); getline(iss, line, ","); getline(iss, line, ",");
	if (stoi(val) != Ny) {
		cout << "ERROR::InOut.h:: Ny defined in .csv file is not the one expected: " << stoi(val) << endl;
		return false;
	}


	getline(csv, line); stringstream data(line);

	for (int yi = 0; yi <= Ny+7; yi++)
	{
		for (int xi = 0; xi <= Nx+7; xi++)
		{
			int i = 1;
		}
	}


	return true;
}