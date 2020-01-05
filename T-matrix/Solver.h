#pragma once

#pragma comment(lib, "mkl_intel_lp64.lib")
#pragma comment(lib, "mkl_core.lib")

#if 0
#pragma comment(lib, "mkl_sequential.lib")
#else
#pragma comment(lib, "mkl_intel_thread.lib")
#pragma comment(lib, "libiomp5md.lib")
#endif

#include "mkl.h"
#include <vector>
#include <cstdint>

using std::vector;

class SymmetricSolver
{
private:
	void *pt[64];
	vector<double> aelem;
	vector<MKL_INT64> iptr, jptr;
	MKL_INT64 iparm[64];
public:
	void PutMatrix(const vector<double> &gg, const vector<double> &di, const vector<uint32_t> &ig, const vector<uint32_t> &jg);
	int Solve(vector<double> &gg, vector<double> &di, vector<uint32_t> &ig, vector<uint32_t> &jg, vector<double> &b, vector<double> &x);
	int Solve(vector<double> &ag, vector<MKL_INT64> &ig, vector<MKL_INT64> &jg, vector<double> &b, vector<double> &x);
	SymmetricSolver();
	~SymmetricSolver();
};

