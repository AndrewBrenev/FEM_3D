#ifndef PROBLEM_H_
#define PROBLEM_H_
#define EL_SIZE 8
//----------------------------
#include "Solver.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

const double lambda = 1;
const double gamma = 1;

// Локальные элементы
struct ELEM
{
	static int count;
	int mas[EL_SIZE];
	double lambda;
	double gamma;
};

//Координаты узлов
struct NODE
{
	double mas[3];
	static int count;
};


class MATRIX
{
public:
	int n;
	vector<double> di;
	vector<double> gg;
	vector<uint32_t> ig;
	vector<uint32_t> jg;
	~MATRIX();
};



class PROBLEM_3D: public MATRIX
{
private:
	ELEM *elems;
	NODE *nodes;
	SymmetricSolver solver;

	vector<double> f;
	
	vector<double> vect_v;
	vector<double> x;

public:

	const char *FILE_ELEM_NUMBER = "local.txt";
	const char *FILE_NODE = "cross.txt";
	const char *FILE_OUT;
	const char *FILE_CRAEV_1 = "boundary.txt";
	const char *FILE_CRAEV_2;
	const char *FILE_CRAEV_3;
	double POR = 1e30;

	PROBLEM_3D();
	~PROBLEM_3D();
	void read_ELEM();
	void read_NODE();
	void buildPortrait();
	void fillTheMatrix();
	void add(int,int,double);
	
	void kraev_1();
	void kraev_2();
	void kraev_3();
	void print_result();

	bool readFromFiles(const char* inf,  const char *xyz, const char *nver) {
		FILE* file = fopen(inf, "r");
		if (file == NULL) return false;
		int nodes_size, el_size;
		fscanf(file, " ISLAU=       0 INDKU1=       0 INDFPO=       0\nKUZLOV=%8d   KPAR=%8d    KT1=       0   KTR2=       0   KTR3=       0\nKISRS1=       0 KISRS2=       0 KISRS3=       0   KBRS=       0\n   KT7=       0   KT10=       0   KTR4=       0  KTSIM=       0\n   KT6=       0\n", &nodes_size, &el_size);
		fclose(file);

		
		NODE::count = nodes_size ;
		ELEM::count = el_size;
		n = NODE::count;


		nodes = new NODE[nodes_size];
		elems = new ELEM[ELEM::count];

		int *tmp = new int[el_size];
		size_t t;
		FILE *fout;
	
		double *tmp1 = new double[3 * nodes_size];
		if (fopen_s(&fout, xyz, "rb")) return false;
		t = fread(tmp1, sizeof(double), 3 * nodes_size, fout);
		int k = feof(fout);
		fclose(fout);

		if (fopen_s(&fout, nver, "rb")) return false;
		int  *tmp2 = new int[el_size * 14];
		t = fread(tmp2, sizeof(int), el_size * 14, fout);

		//Преобразуем считанные данные в массивы
		double a[3];
		for (int i = 0, k = 0; i < 3 * nodes_size; i += 3, k++) {
			nodes[k].mas[0] = tmp1[i];
			nodes[k].mas[1] = tmp1[i + 1];
			nodes[k].mas[2] = tmp1[i + 2];
		}

		for (int i = 0, k = 0; k < el_size; k++, i = i + 14) {
			for (int j = 0; j < 8; j++) {
				elems[k].mas[j] = tmp2[i + j];
				elems[k].lambda = lambda;
				elems[k].gamma = gamma;
			}
		}

		di.resize(n);
		f.resize(n);
		x.resize(n);
		ig.resize(n + 1);

		return true;
	};

	double U(double nx, double y, double z);
	double get_psi(int num_fe, int num_basis, double x,double y, double z);

	void solveMatrix() {
		for (int i = 0; i < ig.size(); i++)
			ig[i]--;
		for (int i = 0; i < jg.size(); i++)
			jg[i]--;
		solver.Solve(gg, di, ig, jg, f, x);
	}

};

#endif
