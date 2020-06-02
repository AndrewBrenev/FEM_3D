#pragma once

#include "stdafx.h"

class Grid {
public:
	vector<ELEM> elems;
	vector<NODE> nodes;
	vector<uint32_t> firstBoundary;
	vector<PLANE> secondBoundary;
	size_t n;

	bool readGridFromFiles(const char* inf, const char* xyz, const char* nver);
	bool buildTestInconsistentGrid();

};
