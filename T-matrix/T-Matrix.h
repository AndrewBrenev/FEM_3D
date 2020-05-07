#pragma once

#include "stdafx.h"
#include <fstream>

// matrix in CSC format
class T_MATRIX :public MATRIX
{
public:
	int nc;
	int nt;

	T_MATRIX() {};
	T_MATRIX(vector<double> _gg, vector<uint32_t> _ig, vector<uint32_t> _jg) {
		gg = _gg;
		ig = _ig;
		jg = _jg;
	};
	~T_MATRIX() {};

	void readFromFiles(const char* config, const char* ig, const char* jg, const char* gg);
	vector<uint32_t> getNotZeroRowsOfColum(const uint32_t i);
	double getElement(const uint32_t i, const uint32_t j);
};
