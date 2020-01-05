//---------------------------------------------------------------------------
#include "problem.h"
//---------------------------------------------------------------------------


int G1[2][2] = { {1,-1},{-1,1} };
double M1[2][2] = { {1. / 3,1. / 6} , {1. / 6,1. / 3} };


// -div( lambda grad(u)) + gamma*u  = right

double right(double* x)
{
	//return x[0] + x[1] + x[2];
	return 5;
}
double u_real(double* x)
{
	//return x[0] + x[1] + x[2];
	return 5;
}

//u с хвостиком
int u(int i) {
	return i / 4;
}
int my(int i) {
	return i % 2;
}
int v(int i) {
	return (i / 2) % 2;
}


//--------------------------------
void mul_c_vector(double* vector, double* h)
{
	double result[EL_SIZE];

	for (int i = 0; i < EL_SIZE; i++) {
		result[i] = 0;
		for (int j = 0; j < EL_SIZE; j++) {
			double value = h[0] * M1[my(i)][my(j)] * h[1] * M1[v(i)][v(j)] * h[2] * M1[u(i)][u(j)];
			result[i] += value * vector[j];
		}
	}

	for (int i = 0; i < EL_SIZE; i++)
		vector[i] = result[i];
}


double PROBLEM_3D::Gauss(double (*f)(double), double begin, double end)
{
	double iLength = (end - begin) / GAUSS_FRAGMENTATION;
	double sum = 0.0;
	for (int i = 0; i < GAUSS_FRAGMENTATION; ++i)
	{
		double a = begin + i * iLength;
		double b = begin + (i + 1) * iLength;

		const double Xi[GAUSS_DEGREE] = { -0.9061798,-0.5384693,0,0.5384693,0.9061798 };
		const double Ci[GAUSS_DEGREE] = { 0.4786287,0.2369269,0.5688888 ,0.2369269,0.4786287 };

		double ra = (b - a) / 2;
		double su = (a + b) / 2;
		double Q, S = 0.0;
		for (int i = 0; i < GAUSS_DEGREE; i++)
		{
			Q = su + ra * Xi[i];
			S += Ci[i] * f(Q);
		}
		sum += ra * S;
	}
	return sum;
}

//------------------------------------------------------
void PROBLEM_3D::setSlaeExactValue(uint32_t i, double value, MATRIX& slae)
{
	if (i >= slae.di.size())
		throw std::invalid_argument("Out of matrix range");

	f[i] = value;
	slae.di[i] = 1;

	for (uint32_t k = slae.ig[i]; k < slae.ig[i + 1]; k++)
	{
		uint32_t j = slae.jg[k];
		f[j] -= value * slae.gg[k];
		slae.gg[k] = 0.0;
	}

	for (uint32_t ii = i + 1; ii < slae.di.size(); ii++)
	{
		for (uint32_t k = slae.ig[ii]; k < slae.ig[ii + 1]; k++)
		{
			uint32_t j = slae.jg[k];
			if (j == i)
			{
				f[ii] -= value * slae.gg[k];
				slae.gg[k] = 0.0;
				break;
			}
		}
	}


}
//------------------------------------------------------
vector<double> PROBLEM_3D::solveM(vector<double>& f) {

	vector<double> x_res(f.size());

	M = slae;

	// 0 - dx, 1 - dy, 3 - dz 
	double h[3];


	for (int k = 0; k < elems.size(); k++)
	{
		//¬ычисление шага
		//h----------------------------------------------------------
		h[0] = abs(nodes[elems[k].value[0]].value[0] - nodes[elems[k].value[1]].value[0]);
		h[1] = abs(nodes[elems[k].value[0]].value[1] - nodes[elems[k].value[2]].value[1]);
		h[2] = abs(nodes[elems[k].value[0]].value[2] - nodes[elems[k].value[4]].value[2]);

		//------------------------------------------------------------
		//формирование элементов матрицы масс


		for (int i = 0; i < EL_SIZE; i++)
		{
			for (int j = 0; j <= i; j++)
			{
				double Mij = elems[k].gamma * h[0] * M1[my(i)][my(j)] * h[1] * M1[v(i)][v(j)] * h[2] * M1[u(i)][u(j)];

				if (i != j)
				{
					//добавка в gg
					(elems[k].value[i] < elems[k].value[j]) ?
						addM(elems[k].value[j], elems[k].value[i], Mij) :
						addM(elems[k].value[i], elems[k].value[j], Mij);
				}
				else
					//добавка в диагональ
					M.di[elems[k].value[i]] += Mij;
			}


		}
	}

	solver.Solve(M.gg, M.di, M.ig, M.jg, f, x_res);

	return x_res;
}
//------------------------------------------------------

void PROBLEM_3D::buildPortrait()
{
	int** prof;
	prof = new int* [slae.n - 1];
	for (int i = 0; i < slae.n - 1; i++)
		prof[i] = new int[i + 1];

	for (int i = 0; i < slae.n - 1; i++)
		for (int j = 0; j <= i; j++)
			prof[i][j] = 0;

	int s = 0;// ол-во ненулевых эл-тов
	for (int k = 0; k < elems.size(); k++)
		for (int i = 0; i < EL_SIZE; i++)
			for (int j = i + 1; j < EL_SIZE; j++) {
				int i1 = elems[k].value[i], j1 = elems[k].value[j], k_b;
				if (i1 < j1) {
					k_b = i1;
					i1 = j1;
					j1 = k_b;
				}
				if (!prof[i1 - 1][j1]) {
					s++;
					prof[i1 - 1][j1] = 1;
				}
			}

	//‘ормирование массива ig и jg	
	slae.jg.resize(s);
	slae.gg.resize(s);
	slae.ig[0] = 0;
	slae.ig[1] = 0;
	for (int i = 0, d = 0; i < slae.n - 1; i++) {
		int k = 0;
		for (int j = 0; j <= i; j++)
			if (prof[i][j] == 1) {
				slae.jg[d] = j;
				d++;
				k++;
			}
		slae.ig[i + 2] = slae.ig[i + 1] + k;
	}

	for (int i = 0; i < slae.n - 1; i++)
		delete[] prof[i];
	delete[] prof;
}
//------------------------------------------------------
void PROBLEM_3D::add(int i, int j, double x)
{
	int k;
	for (k = slae.ig[i]; k < slae.ig[i + 1]; k++)
		if (slae.jg[k] == j)
			break;
	slae.gg[k] += x;
}
//------------------------------------------------------
void PROBLEM_3D::addM(int i, int j, double x)
{
	int k;
	for (k = M.ig[i]; k < M.ig[i + 1]; k++)
		if (M.jg[k] == j)
			break;
	M.gg[k] += x;
}
//------------------------------------------------------
void PROBLEM_3D::fillTheMatrix()
{
	double right_vector[EL_SIZE];

	// 0 - dx, 1 - dy, 3 - dz 
	double h[3];


	for (int k = 0; k < elems.size(); k++)
	{
		//¬ычисление шага
		//h----------------------------------------------------------
		h[0] = abs(nodes[elems[k].value[0]].value[0] - nodes[elems[k].value[1]].value[0]);
		h[1] = abs(nodes[elems[k].value[0]].value[1] - nodes[elems[k].value[2]].value[1]);
		h[2] = abs(nodes[elems[k].value[0]].value[2] - nodes[elems[k].value[4]].value[2]);

		//------------------------------------------------------------
		//формирование элементов матрицы

		//¬ычисление вектоора правой части
		for (int i = 0; i < EL_SIZE; i++)
			right_vector[i] = right(nodes[elems[k].value[i]].value);

		mul_c_vector(right_vector, h);

		double gX = 1. / h[0], gY = 1. / h[1], gZ = 1. / h[2];

		for (int i = 0; i < EL_SIZE; i++)
		{
			for (int j = 0; j <= i; j++)
			{
				double Gij =
					gX * G1[my(i)][my(j)] * h[1] * M1[v(i)][v(j)] * h[2] * M1[u(i)][u(j)] +
					h[0] * M1[my(i)][my(j)] * gY * G1[v(i)][v(j)] * h[2] * M1[u(i)][u(j)] +
					h[0] * M1[my(i)][my(j)] * h[1] * M1[v(i)][v(j)] * gZ * G1[u(i)][u(j)];

				double Mij = h[0] * M1[my(i)][my(j)] * h[1] * M1[v(i)][v(j)] * h[2] * M1[u(i)][u(j)];

				double value = elems[k].lambda * Gij + elems[k].gamma * Mij;

				if (i != j)
				{
					//добавка в gg
					(elems[k].value[i] < elems[k].value[j]) ?
						add(elems[k].value[j], elems[k].value[i], value) :
						add(elems[k].value[i], elems[k].value[j], value);
				}
				else
					//добавка в диагональ
					slae.di[elems[k].value[i]] += value;
			}

			//добавка в правую часть
			f[elems[k].value[i]] += right_vector[i];
		}
	}

	for (auto it = firstBoundary.begin(); it < firstBoundary.end(); it++)
		setSlaeExactValue(*it, u_real(nodes[*it].value), slae);


}
//------------------------------------------------------

PROBLEM_3D::~PROBLEM_3D()
{
	elems.clear();
	nodes.clear();
	f.clear();
	x.clear();
}

bool PROBLEM_3D::readFromFiles(const char* inf, const char* xyz, const char* nver)
{
	FILE* file = fopen(inf, "r");
	if (file == NULL) return false;
	int nodes_size, el_size;
	fscanf(file, " ISLAU=       0 INDKU1=       0 INDFPO=       0\nKUZLOV=%8d   KPAR=%8d    KT1=       0   KTR2=       0   KTR3=       0\nKISRS1=       0 KISRS2=       0 KISRS3=       0   KBRS=       0\n   KT7=       0   KT10=       0   KTR4=       0  KTSIM=       0\n   KT6=       0\n", &nodes_size, &el_size);
	fclose(file);

	slae.n = nodes_size;


	nodes.resize(nodes_size);
	nodes.reserve(nodes_size);

	elems.reserve(el_size);
	elems.resize(el_size);

	int* tmp = new int[el_size];
	size_t t;
	FILE* fout;

	double* tmp1 = new double[3 * nodes_size];
	if (fopen_s(&fout, xyz, "rb")) return false;
	t = fread(tmp1, sizeof(double), 3 * nodes_size, fout);
	int k = feof(fout);
	fclose(fout);

	if (fopen_s(&fout, nver, "rb")) return false;
	int* tmp2 = new int[el_size * 14];
	t = fread(tmp2, sizeof(int), el_size * 14, fout);

	//ѕреобразуем считанные данные в массивы
	double a[3];
	for (int i = 0, k = 0; i < 3 * nodes_size; i += 3, k++) {
		nodes[k].value[0] = tmp1[i];
		nodes[k].value[1] = tmp1[i + 1];
		nodes[k].value[2] = tmp1[i + 2];
	}

	for (int i = 0, k = 0; k < el_size; k++, i = i + 14) {
		for (int j = 0; j < 8; j++) {
			elems[k].value[j] = tmp2[i + j];
			elems[k].lambda = LAMBDA;
			elems[k].gamma = GAMMA1;
		}
	}

	slae.di.resize(slae.n);
	f.resize(slae.n);
	x.resize(slae.n);
	slae.ig.resize(slae.n + 1);

	return true;
};

void PROBLEM_3D::solveMatrix()
{
	solver.Solve(slae.gg, slae.di, slae.ig, slae.jg, f, x);
}

//--------------------------------------
void PROBLEM_3D::print_result()
{
	FILE* f = fopen("resault.txt", "w");
	for (int i = 0; i < slae.n; i++)
	{
		double real_solution = u_real(nodes[i].value);

		fprintf(f, "%.3lf %.3lf %.3lf %.15lf %.15lf %.15lf\n", nodes[i].value[0], nodes[i].value[1], nodes[i].value[2], real_solution, x[i], abs(real_solution - x[i]));
	}
	fclose(f);
}