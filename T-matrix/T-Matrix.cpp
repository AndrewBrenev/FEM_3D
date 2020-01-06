#include "T-Matrix.h"


void T_MATRIX::readFromFiles(const char* config, const char* _ig, const char* _jg, const char* _gg) {

	// read file config
	ifstream fin(config, ios::in);
	fin >> nc;
	fin >> nt;
	fin.close();

	ig.resize(nt + 1);

	fin = ifstream(_ig, ios::in);
	for (int i = 0; i < nt + 1; ++i)
		fin >> ig[i];
	fin.close();

	jg.resize(ig[nt]);
	gg.resize(ig[nt]);
	fin = ifstream(_jg, ios::in); 
	ifstream gg_in(_gg, ios::in);
	for (int i = 0; i < ig[nt]; ++i) {
		fin >> jg[i];
		gg_in >> gg[i];
	}
	gg_in.close();
	fin.close();


};

vector<uint32_t> T_MATRIX::getNotZeroRowsOfColum(const uint32_t i) 
{
	int columId = i - nc;
	vector<uint32_t> notZeroRows;

	for (int j = ig[columId]; j < ig[columId + 1]; ++j) 
		notZeroRows.push_back(jg[j]);

	return notZeroRows;
}

double T_MATRIX::getElement(const uint32_t i, const uint32_t j) {
	int columId = j - nc;

	for (int _j = ig[columId]; _j < ig[columId + 1]; ++_j)
		if (jg[_j] == i)
			return gg[_j];
	return 0;
}