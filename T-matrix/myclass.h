#ifndef MYCLASS_H_
#define MYCLASS_H_
#define M 8
//----------------------------
//#include "local_matrix.h"
//#include <alloc.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

//----------------------------
double right(double *x);
void mul_c_vector(double*,double*);

const double lambda = 1;
const double gamma = 1;

// Локальные области
struct local
{
	static int count;
	int mas[M];
	double lambda;
	double gamma;
};

//Координаты узлов
struct cross
{
	double mas[3];
	static int count;
};



class MATRIX
{
public:
	int n;
	double *d;
	double *gg;
	int *ig;
	int *jg;


	void SST(MATRIX *);
	void mul_matrix_vector(double *vector,double *result);
	void solution_x(double*,double*);
	void solution_x_l(double*,double*);
	void solution_x_u(double*,double*);
	~MATRIX();
};



class GLOBAL_MATRIX : public MATRIX
{
private:
	local *matr;
	cross *set;
	double *f;
	double *x;
	double *vect_v;


public:
	GLOBAL_MATRIX();
	~GLOBAL_MATRIX();
	void read_local();
	void read_cross();
	void formier_profil();
	void formier_matrix();
	void add(int,int,double);
	void MSG();
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

		
		cross::count = nodes_size ;
		local::count = el_size;
		n = cross::count;


		set = new cross[nodes_size];
		matr = new local[local::count];

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
			set[k].mas[0] = tmp1[i];
			set[k].mas[1] = tmp1[i + 1];
			set[k].mas[2] = tmp1[i + 2];
		}

		for (int i = 0, k = 0; k < el_size; k++, i = i + 14) {
			for (int j = 0; j < 8; j++) {
				matr[k].mas[j] = tmp2[i + j];
				matr[k].lambda = lambda;
				matr[k].gamma = gamma;
			}
		}

		return true;
	};

	double U(double nx, double y, double z);
	double get_psi(int num_fe, int num_basis, double x,double y, double z);
	const char *FILE_LOCAL_NUMBER;
	const char *FILE_CROSS;
	const char *FILE_OUT;
	const char *FILE_CRAEV_1;
	const char *FILE_CRAEV_2;
	const char *FILE_CRAEV_3;
	double POR;
};

#endif
