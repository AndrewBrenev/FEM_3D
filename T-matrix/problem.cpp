//---------------------------------------------------------------------------
#include "problem.h"
//---------------------------------------------------------------------------

int G1[2][2] = { {1,-1},{-1,1} };
int M1[2][2] = { {2,1},{1,2} };


// -div( lambda grad(u)) + gamma*u  = right

double right(double *x)
{
	return 5;

}
double u_real(double *x)
{
	return 5;
}

int my(int i) {
	return i % 2;
}
int v(int i) {
	return (i/2) % 2;
}
//u с хвостиком
int u(int i) {
	return i / 4;
}

//--------------------------------
void mul_c_vector(double *vector,double *h)
{
	double result[EL_SIZE];

	for (int i = 0; i < EL_SIZE; i++) {
		result[i] = 0;
		for (int j = 0; j < EL_SIZE; j++) {
			double value = (h[0] / 6) * M1[my(i)][my(j)] * (h[1] / 6) * M1[v(i)][v(j)] * (h[2] / 6) * M1[u(i)][u(j)];
			result[i] += value * vector[j];
		}
	}

	for (int i = 0; i < EL_SIZE; i++)
		vector[i] = result[i];
}

int ELEM::count=0;
int NODE::count=0;

MATRIX::~MATRIX()
{
	di.clear();
	ig.clear();
	gg.clear();
	jg.clear();
}

//------------------------------------------------------
void PROBLEM_3D::buildPortrait()
{
	int** prof;
	prof = new int* [n - 1];
	for (int i = 0; i < n - 1; i++)
		prof[i] = new int[i + 1];

	for (int i = 0; i < n - 1; i++)
		for (int j = 0; j <= i; j++)
			prof[i][j] = 0;

	int s = 0;// ол-во ненулевых эл-тов
	for (int k = 0; k < ELEM::count; k++)
		for (int i = 0; i < EL_SIZE; i++)
			for (int j = i + 1; j < EL_SIZE; j++) {
				int i1 = elems[k].value[i], j1 = elems[k].value[j], k_b;
				if (i1 < j1) {
					k_b = i1;
					i1 = j1;
					j1 = k_b;
				}
				if (prof[i1 - 2][j1 - 1] == 0) {
					s++;
					prof[i1 - 2][j1 - 1] = 1;
				}
			}

	//‘ормирование массива ig и jg	
	jg.resize(s);
	gg.resize(s);
	ig[0] = 1;
	ig[1] = 1;
	for (int i = 0, d = 0; i < n - 1; i++) {
		int k = 0;
		for (int j = 0; j <= i; j++)
			if (prof[i][j] == 1) {
				jg[d] = j + 1;
				d++;
				k++;
			}
		ig[i + 2] = ig[i + 1] + k;
	}
	/*
	printf("\n");
	for(int i =0; i <= n; i++)
		printf(" %d",ig[i]);
	printf("\n");
	for(int i =0; i < s; i++)
		printf(" %d",jg[i]);
		*/
	for (int i = 0; i < n - 1; i++)
		delete[] prof[i];
	delete[] prof;
}
//------------------------------------------------------
void PROBLEM_3D::add(int i, int j, double x)
{
	int k;
	for (k = ig[i] - 1; k < ig[i + 1] - 1; k++)
		if (jg[k] == j + 1)
			break;
	gg[k] += x;
}
//------------------------------------------------------
void PROBLEM_3D::fillTheMatrix()
{
	double right_vector[EL_SIZE];
	double h[3];
	for (int k = 0; k < ELEM::count; k++)
	{
		//¬ычисление шага
		//h----------------------------------------------------------
		for (int i = 0; i < 3; i++) { //i-по x,y,z    
			int flag = 1;
			int j;
			for (j = 1; j < EL_SIZE && flag; j++) {//1 узел фиксируем,пробегаем по остальным    
				flag = 0;
				for (int l = 0; l < 3 && !flag; l++)//провер€ем,лежат ли точки на нужном ребре
					if (i != l && nodes[elems[k].value[0] - 1].value[l] != nodes[elems[k].value[j] - 1].value[l])
						flag = 1;
			}
			if (!flag)
				h[i] = fabs(nodes[elems[k].value[0] - 1].value[i] - nodes[elems[k].value[j - 1] - 1].value[i]);
		}
		//------------------------------------------------------------
		//формирование элементов матрицы

		//¬ычисление вектоора правой части
		for (int i = 0; i < EL_SIZE; i++)
			right_vector[i] = right(nodes[elems[k].value[i] - 1].value);
		mul_c_vector(right_vector, h);

		double gX = 1 / h[0], gY = 1 / h[1], gZ = 1 / h[2];
		for (int i = 0; i < EL_SIZE; i++)
		{
			for (int j = 0; j <= i; j++)
			{
				double Gij =
					gX * G1[my(i)][my(j)] * (h[1] / 6) * M1[v(i)][v(j)] * (h[2] / 6) * M1[u(i)][u(j)] +
					(h[0] / 6) * M1[my(i)][my(j)] * gY * G1[v(i)][v(j)] * (h[2] / 6) * M1[u(i)][u(j)] +
					(h[0] / 6) * M1[my(i)][my(j)] * (h[1] / 6) * M1[v(i)][v(j)] * gZ * G1[u(i)][u(j)];

				double Mij = (h[0] / 6) * M1[my(i)][my(j)] * (h[1] / 6) * M1[v(i)][v(j)] * (h[2] / 6) * M1[u(i)][u(j)];

				double value = elems[k].lambda * Gij + elems[k].gamma * Mij;

				if (i != j)
				{
					//добавка в gg
					(elems[k].value[i] < elems[k].value[j]) ?
						add(elems[k].value[j] - 1, elems[k].value[i] - 1, value) :
						add(elems[k].value[i] - 1, elems[k].value[j] - 1, value);
				}
				else
					//добавка в диагональ
					di[elems[k].value[i] - 1] = value;
			}

			//добавка в правую часть
			f[elems[k].value[i] - 1] += right_vector[i];
		}
	}
}
//------------------------------------------------------
PROBLEM_3D:: PROBLEM_3D(){} 

PROBLEM_3D::~PROBLEM_3D()
{
	delete [] elems;
	delete [] nodes;
	f.clear();
	x.clear();
}

//--------------------------------------
void PROBLEM_3D::print_result()
{
	FILE* f = fopen("resault.txt", "w");
	for (int i = 0; i < n; i++)
	{
		double real_solution = u_real(nodes[i].value);

		fprintf(f, "%.3lf %.3lf %.3lf %.15lf %.15lf %.15lf\n", nodes[i].value[0], nodes[i].value[1], nodes[i].value[2], real_solution, x[i], abs(real_solution - x[i]));
	}
	fclose(f);
}