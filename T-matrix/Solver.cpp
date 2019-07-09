#include "Solver.h"

#include <cstdint>

void SymmetricSolver::PutMatrix(const vector<double> &gg, const vector<double> &di, const vector<uint32_t> &ig, const vector<uint32_t> &jg)
{
	vector<MKL_INT64> col;
	col.resize(di.size(), 0);

	for (uint32_t i = 0; i < ig.size() - 1; i++)
	{
		for (uint32_t k = ig[i]; k < ig[i + 1]; k++)
		{
			uint32_t j = jg[k];
			col[j]++;
		}
		col[i]++;
	}

	iptr.resize(col.size() + 1);
	iptr[0] = 0;

	for (uint32_t i = 0; i < col.size(); i++)
	{
		iptr[i + 1] = iptr[i] + col[i];
		col[i] = iptr[i];
	}

	aelem.resize((size_t) iptr.back());
	jptr.resize((size_t) iptr.back());
	for (uint32_t i = 0; i < di.size(); i++)
	{
		jptr[col[i]] = i;
		aelem[col[i]] = di[i];
		col[i]++;
	}

	for (uint32_t i = 0; i < ig.size() - 1; i++)
	{
		for (uint32_t k = ig[i]; k < ig[i + 1]; k++)
		{
			int j = jg[k];
			jptr[col[j]] = i;
			aelem[col[j]] = gg[k];
			col[j]++;
		}
	}


}

int SymmetricSolver::Solve(vector<double> &gg, vector<double> &di, vector<uint32_t> &ig, vector<uint32_t> &jg, vector<double>& b, vector<double>& x)
{
	PutMatrix(gg, di, ig, jg);
	return Solve(aelem, iptr, jptr, b, x);
}

int SymmetricSolver::Solve(vector<double> &ag, vector<MKL_INT64> &ig, vector<MKL_INT64> &jg, vector<double>& b, vector<double>& x)
{
	x.resize(ig.size() - 1, 0);
	MKL_INT64 maxfct, mnum, phase, error, msglvl, mtype, nrhs;
	maxfct = 1;
	mnum = 1;
	msglvl = 0;
	error = 0;

	MKL_INT64 N = ig.size() - 1;
	//double ddum;          /* Double dummy */
	MKL_INT64 idum;         /* Integer dummy. */

	nrhs = 1;
	mtype = -2;
	phase = 13;
	
	PARDISO_64(pt, &maxfct, &mnum, &mtype, &phase,
		&N, ag.data(), ig.data(), jg.data(), &idum, &nrhs, iparm, &msglvl, b.data(), x.data(), &error);

	MKL_INT64 save_error = error;

	return (int)save_error;
}

SymmetricSolver::SymmetricSolver()
{
	memset(pt, 0, sizeof(pt));
	memset(iparm, 0, sizeof(iparm));
	iparm[0] = 1;         /* No solver default */
	iparm[1] = 2;         /* Fill-in reordering from METIS */
	iparm[3] = 0;         /* No iterative-direct algorithm */
	iparm[4] = 0;         /* No user fill-in reducing permutation */
	iparm[5] = 0;         /* Write solution into x */
	iparm[6] = 0;         /* Not in use */
	iparm[7] = 2;         /* Max numbers of iterative refinement steps */
	iparm[8] = 0;         /* Not in use */
	iparm[9] = 13;        /* Perturb the pivot elements with 1E-13 */
	iparm[10] = 1;        /* Use nonsymmetric permutation and scaling MPS */
	iparm[11] = 0;        /* Conjugate transposed/transpose solve */
	iparm[12] = 1;        /* Maximum weighted matching algorithm is switched-on (default for non-symmetric) */
	iparm[13] = 0;        /* Output: Number of perturbed pivots */
	iparm[14] = 0;        /* Not in use */
	iparm[15] = 0;        /* Not in use */
	iparm[16] = 0;        /* Not in use */
	iparm[17] = -1;       /* Output: Number of nonzeros in the factor LU */
	iparm[18] = -1;       /* Output: Mflops for LU factorization */
	iparm[19] = 0;        /* Output: Numbers of CG Iterations */
	iparm[34] = 1;
}

SymmetricSolver::~SymmetricSolver()
{
	MKL_INT64 phase = -1;
	double ddum;          /* Double dummy */
	MKL_INT64 idum, msglvl, error;         /* Integer dummy. */
	msglvl = 0;
	PARDISO_64(pt, &idum, &idum, &idum, &phase,
		&idum, &ddum, &idum, &idum, &idum, &idum,
		iparm, &msglvl, &ddum, &ddum, &error);
}
