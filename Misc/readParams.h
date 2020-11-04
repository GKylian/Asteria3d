#pragma once

#include <iostream>
#include <string>
#include <map>
#include <fstream>

using namespace std;

/* Reads the blocks in the input .params file and puts them (also by block) in a 2d c++ map. 
/  They can be set/accessed through params["blockname"]["parameter name"].
/  --> They are stored as strings for a greater flexibility, so a call to stoi or stod is necessary when dealing with numbers
/  --> A map is probably not the least computationally expensive choice, but since most of (if not all of) the
/      parameters are only read at the beginning of the simulation, it doesn't change much and is way simpler to implement.
*/


bool readParams(string fname, map<string, map<string, string>> *params) {
	ifstream input(fname);
	cout << "Reading simulation parameters from prob.params file..." << endl;

	string block = ""; /* What block are we in (if any) ? */
	for (string line; getline(input, line);) {
		if (!(line.find('_') != string::npos) && !(line.find('<') != string::npos)) continue; /* Line without any parameter or block name -> skip */
		if (line.find('_') != string::npos) {
			unsigned first = line.find('_');
			unsigned last = line.find_last_of('_');
			block = line.substr(first+1, last-first-1);
			cout << "Entering block named " << block << endl;
			continue;
		}
		/* Split the line in two (before and after the equal */
		string left = line.substr(0, line.find('='));
		string right = line.substr(line.find('=')+1, line.length()-1);

		string key = ""; string val = "";
		if (left != "") {
			unsigned first = left.find('<'); unsigned last = left.find('>');
			key = left.substr(first+1, last-first-1);
		}
		if (right != "") {
			unsigned first = right.find('<'); unsigned last = right.find('>');
			val = right.substr(first+1, last-first-1);
		}
		


		(*params)[block][key] = val;
	}

	return true;
}