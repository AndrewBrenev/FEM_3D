#ifndef PROBLEM_H_

#define PROBLEM_H_

#define EL_SIZE 8

//----------------------------
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <vector>
#include <stdexcept>

#include "Solver.h"

using namespace std;

const double LAMBDA = 1;

const double GAMMA1 = 1;
const double GAMMA2 = 1;


// Локальные элементы
struct ELEM
{
	int id;
	int value[EL_SIZE];
	double lambda;
	double gamma;
};

//Координаты узлов
struct NODE
{
	// 0 - x ,1 - y , 2 - z
	double value[3];
	int id;
};

struct MATRIX
{
	int n;
	vector<double> di;
	vector<double> gg;
	vector<uint32_t> ig;
	vector<uint32_t> jg;
};

struct REGION {
	double X0, Xn, Xh;
	double Y0, Yn, Yh;
	double Z0, Zn, Zh;
};

class PROBLEM_3D
{
private:

	MATRIX slae, M;
	vector<ELEM> elems;
	vector<NODE> nodes;

	SymmetricSolver solver;

	const REGION outter, inner;

	vector<double> f;
	vector<double> x;

	vector<uint32_t> firstBoundary;

	void setSlaeExactValue(uint32_t, double, MATRIX&);

	vector<double> solveM(vector<double>&);
public:

	PROBLEM_3D();

	PROBLEM_3D(const REGION& _outter, const REGION& _inner) :inner(_inner), outter(_outter)
	{	};

	~PROBLEM_3D();


	bool readFromFiles(const char* inf, const char* xyz, const char* nver);

	void buildPortrait();

	void fillTheMatrix();

	void add(int, int, double);
	void addM(int, int, double);

	void print_result();

	void buildGrid();

	void solveMatrix();

};

#endif
