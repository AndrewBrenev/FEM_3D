#include "problem.h"

int G1[2][2] = { {1,-1},{-1,1} };
double M1[2][2] = { { 1. / 3., 1. / 6. }, { 1. / 6., 1. / 3. } };

//u с хвостиком
int u(int i) { return i / 4; };
int my(int i) { return i % 2; };
int v(int i) { return (i / 2) % 2; };

// -div( lambda grad(u)) + gamma*u  = right

double right(double* x)
{
	return 4 * (x[0] * x[0] + x[1] * x[1] + x[2] * x[2]) + 18;
	//return 5;
}
double u_real(double* x)
{
	return x[0] * x[0] + x[1] * x[1] + x[2] * x[2];
	//return 5;
}



void mul_c_vector(vector <double>& _vector, double* h) {
	vector <double> result(EL_SIZE);

	for (int i = 0; i < EL_SIZE; ++i)
		for (int j = 0; j < EL_SIZE; ++j)
			result[i] += h[0] * M1[my(i)][my(j)] * h[1] * M1[v(i)][v(j)] * h[2] * M1[u(i)][u(j)] * _vector[j];

	_vector = result;
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

vector<double> PROBLEM_3D::solveM(vector<double>& f) {

	vector<double> x_res(f.size());

	M = slae;

	// 0 - dx, 1 - dy, 3 - dz 
	double h[3];


	for (int k = 0; k < grid.elems.size(); k++)
	{
		//Вычисление шага
		//h----------------------------------------------------------
		h[0] = abs(grid.nodes[grid.elems[k].value[0]].value[0] - grid.nodes[grid.elems[k].value[1]].value[0]);
		h[1] = abs(grid.nodes[grid.elems[k].value[0]].value[1] - grid.nodes[grid.elems[k].value[2]].value[1]);
		h[2] = abs(grid.nodes[grid.elems[k].value[0]].value[2] - grid.nodes[grid.elems[k].value[4]].value[2]);

		//------------------------------------------------------------
		//формирование элементов матрицы масс
		for (int i = 0; i < EL_SIZE; i++)
			for (int j = 0; j <= i; j++)
			{
				double Mij = grid.elems[k].gamma * h[0] * M1[my(i)][my(j)] * h[1] * M1[v(i)][v(j)] * h[2] * M1[u(i)][u(j)];

				M.addValue(grid.elems[k].value[j], grid.elems[k].value[i], Mij);
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

	int s = 0;//Кол-во ненулевых эл-тов
	for (int k = 0; k < grid.elems.size(); k++)
		for (int i = 0; i < EL_SIZE; i++)
			for (int j = i + 1; j < EL_SIZE; j++) {
				int i1 = grid.elems[k].value[i], j1 = grid.elems[k].value[j], k_b;
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

	//Формирование массива ig и jg	
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

void PROBLEM_3D::destruct()
{
	slae.gg.clear();
	slae.di.clear();
	slae.ig.clear();
	slae.jg.clear();
	f.clear();
}

void PROBLEM_3D::buildSetBasedPortrait()
{
	cout << "Building matrix portrait...";
	vector<set<uint32_t>> Portrait;

	destruct();

	uint32_t N = (uint32_t)grid.n;

	Portrait.resize(N);

	for (auto elem : grid.elems)
	{
		vector<uint32_t> sorted_points;

		for (uint8_t i = 0; i < EL_SIZE; i++)
		{
			uint32_t node_id = elem.value[i];
			// find proper position in descending order
			auto it = std::lower_bound(sorted_points.begin(), sorted_points.end(), node_id, std::greater<int>());
			sorted_points.insert(it, node_id); // insert before iterator it
		}

		for (int i = 0; i < sorted_points.size(); i++)
		{
			for (int j = i + 1; j < sorted_points.size(); j++)
				Portrait[sorted_points[i]].insert(Portrait[sorted_points[i]].begin(), sorted_points[j]);
		}
	}

	slae.ig.resize(N + 1);
	slae.ig[0] = 0;

	for (int i = 0; i < slae.ig.size() - 1; i++)
	{
		slae.ig[i + 1] = slae.ig[i] + (uint32_t)Portrait[i].size();
	}

	uint32_t elements = slae.ig[N];

	slae.di.resize(N, 0);
	slae.gg.resize(elements, 0);

	slae.jg.reserve(elements);
	slae.jg.clear();
	for (uint32_t i = 0; i < N; i++)
	{
		for (int col : Portrait[i])
		{
			slae.jg.push_back(col);
		}
	}

	/*
	slae.helper.resize(N);

	for (int i = 0; i < N; i++)
	{
		for (int k = slae.ig[i]; k < slae.ig[i + 1]; k++)
		{
			int j = slae.jg[k];
			slae.helper[i].push_back(PortraitHelperElem(i, j, k));
			slae.helper[j].push_back(PortraitHelperElem(i, j, k));
		}
	}
	*/

	f.resize(N, 0);

	slae.gg.shrink_to_fit();
	slae.di.shrink_to_fit();
	slae.ig.shrink_to_fit();
	slae.jg.shrink_to_fit();
	f.shrink_to_fit();

	cout << "done!" << endl;
}

void PROBLEM_3D::processTerminalNode(int i, double value)
{
	if (i < t_matix.nc) {
		if (i == 37)
			int r = 0;
		f[i] += value;
	}
	else
	{
		auto notZeroRows = t_matix.getNotZeroRowsOfColum(i);
		for (int row : notZeroRows) {
			auto tRowValue = t_matix.getElement(row, i);
			f[row] += tRowValue * value;
		}
	}
}
void PROBLEM_3D::processTerminalNode(int i, int j, double value)
{
	if (i < t_matix.nc && j < t_matix.nc)
		slae.addValue(i, j, value);

	if (i >= t_matix.nc && j < t_matix.nc)
	{
		auto notZeroRows = t_matix.getNotZeroRowsOfColum(i);
		for (int row : notZeroRows) {
			auto tRowValue = t_matix.getElement(row, i);
			slae.addValue(row, j, tRowValue *value);
		}
	}

	if (i < t_matix.nc && j >= t_matix.nc)
	{
		auto notZeroRows = t_matix.getNotZeroRowsOfColum(j);
		for (int row : notZeroRows) {
			auto tRowValue = t_matix.getElement( row,j);
			slae.addValue( i,row, tRowValue * value);
		}
	}

	if (i >= t_matix.nc && j >= t_matix.nc)
	{
		auto notZeroRows = t_matix.getNotZeroRowsOfColum(i);
		auto notZeroColums = t_matix.getNotZeroRowsOfColum(j);

		for (int row : notZeroRows) 
		for(int colum: notZeroColums)
		{
			auto tRowValue = t_matix.getElement(row, i);
			auto tColumValue = t_matix.getElement(colum, j);
			slae.addValue(row, colum, tColumValue *tRowValue * value);
		}

	}
}
//------------------------------------------------------
void PROBLEM_3D::fillTheMatrix()
{
	// 0 - dx, 1 - dy, 3 - dz 
	double h[3];

	for (int k = 0; k < grid.elems.size(); k++) {
		//Вычисление шага
		//h----------------------------------------------------------
		h[0] = abs(grid.nodes[grid.elems[k].value[0]].value[0] - grid.nodes[grid.elems[k].value[1]].value[0]);
		h[1] = abs(grid.nodes[grid.elems[k].value[0]].value[1] - grid.nodes[grid.elems[k].value[4]].value[1]);
		h[2] = abs(grid.nodes[grid.elems[k].value[0]].value[2] - grid.nodes[grid.elems[k].value[2]].value[2]);


		double gX = 1. / h[0], gY = 1. / h[1], gZ = 1. / h[2];

		for (int i = 0; i < EL_SIZE; ++i)
			for (int j = 0; j <= i; ++j) {
				double Gij =
					gX * G1[my(i)][my(j)] * h[1] * M1[v(i)][v(j)] * h[2] * M1[u(i)][u(j)] +
					h[0] * M1[my(i)][my(j)] * gY * G1[v(i)][v(j)] * h[2] * M1[u(i)][u(j)] +
					h[0] * M1[my(i)][my(j)] * h[1] * M1[v(i)][v(j)] * gZ * G1[u(i)][u(j)];

				double Mij = h[0] * M1[my(i)][my(j)] * h[1] * M1[v(i)][v(j)] * h[2] * M1[u(i)][u(j)];
				double value = grid.elems[k].lambda * Gij + grid.elems[k].gamma * Mij;

				if (inconsistentProblem)
					processTerminalNode(grid.elems[k].value[j], grid.elems[k].value[i], value);
				else
					slae.addValue(grid.elems[k].value[j], grid.elems[k].value[i], value);
			}
	}
}

void PROBLEM_3D::calculateRightPart() {
	double h[3];
	vector<double> elementRightPart(EL_SIZE);

	for (int k = 0; k < grid.elems.size(); k++) {

		h[0] = abs(grid.nodes[grid.elems[k].value[0]].value[0] - grid.nodes[grid.elems[k].value[1]].value[0]);
		h[1] = abs(grid.nodes[grid.elems[k].value[0]].value[1] - grid.nodes[grid.elems[k].value[4]].value[1]);
		h[2] = abs(grid.nodes[grid.elems[k].value[0]].value[2] - grid.nodes[grid.elems[k].value[2]].value[2]);

		for (int i = 0; i < EL_SIZE; i++)
			elementRightPart[i] = right(grid.nodes[grid.elems[k].value[i]].value);

		mul_c_vector(elementRightPart, h);

		for (int i = 0; i < EL_SIZE; i++)
			if (inconsistentProblem)
				processTerminalNode(grid.elems[k].value[i], elementRightPart[i]);
			else
				f[grid.elems[k].value[i]] += elementRightPart[i];
	}
}

vector<uint32_t> PROBLEM_3D::readFirstBoundaryFromFile(const char* input) {

	// read file config
	ifstream fin(input, ios::in);

	int edges_count;
	fin >> edges_count;
	
	vector<uint32_t> nodes(edges_count);
	nodes.resize(edges_count);

	for (int i = 0; i < edges_count; ++i)
		fin >> nodes[i];
	fin.close();
	return nodes;
};

void PROBLEM_3D::applyFirstEdge() 
{
	grid.firstBoundary = readFirstBoundaryFromFile(boundary_path.c_str());
	for (auto it = grid.firstBoundary.begin(); it < grid.firstBoundary.end(); it++)
		setSlaeExactValue(*it, u_real(grid.nodes[*it].value), slae);
}
PROBLEM_3D::~PROBLEM_3D()
{
	grid.elems.clear();
	grid.nodes.clear();
	f.clear();
	x.clear();
}

void PROBLEM_3D::readGridFromFiles(const char* inf, const char* xyz, const char* nver)
{

	if (inconsistentProblem)
		grid.buildTestInconsistentGrid();
	else
		grid.readGridFromFiles(inf, xyz, nver);
	
	setMatrixSize(grid.n);
	
}
void PROBLEM_3D::setMatrixSize(const size_t size) {
	slae.n = size;
	slae.di.resize(slae.n);
	f.resize(slae.n);
	x.resize(slae.n);
	slae.ig.resize(slae.n + 1);
}


void PROBLEM_3D::solveMatrix() {

	solver.Solve(slae.gg, slae.di, slae.ig, slae.jg, f, x);
}

void PROBLEM_3D::showResault()
{
	double real_solution;
	for (int i = 0; i < slae.n; i++)
	{
		real_solution = u_real(grid.nodes[i].value);
		printf("%.3lf %.3lf %.3lf %.15lf %.15lf %.15lf\n", grid.nodes[i].value[0], grid.nodes[i].value[1], grid.nodes[i].value[2], real_solution, x[i], abs(real_solution - x[i]));
	}
}
//--------------------------------------
void PROBLEM_3D::printResult(const char* fname)
{
	FILE* f = fopen(fname, "w");
	for (int i = 0; i < slae.n; i++)
	{
		double real_solution = u_real(grid.nodes[i].value);

		fprintf(f, "%.3lf %.3lf %.3lf %.15lf %.15lf %.15lf\n", grid.nodes[i].value[0], grid.nodes[i].value[1], grid.nodes[i].value[2], real_solution, x[i], abs(real_solution - x[i]));
	}
	fclose(f);
}



// Jacobian && local matrix

double PROBLEM_3D::getDeterminant()
{
	double result = 0;
	result = jacobianMatrix[0][0] * jacobianMatrix[1][1] * jacobianMatrix[2][2] -
		jacobianMatrix[0][0] * jacobianMatrix[1][2] * jacobianMatrix[2][1] -
		jacobianMatrix[0][1] * jacobianMatrix[1][0] * jacobianMatrix[2][2] +
		jacobianMatrix[0][1] * jacobianMatrix[2][0] * jacobianMatrix[1][2] +
		jacobianMatrix[1][0] * jacobianMatrix[0][2] * jacobianMatrix[2][1] -
		jacobianMatrix[0][2] * jacobianMatrix[1][1] * jacobianMatrix[2][0];

	return result;
}

void PROBLEM_3D::calculateInverseMatrix()
{
	double determinant = getDeterminant();

	inverseMatrix[0][0] = (jacobianMatrix[1][1] * jacobianMatrix[2][2] - jacobianMatrix[1][2] * jacobianMatrix[2][1]) / fabs(determinant);
	inverseMatrix[0][1] = (jacobianMatrix[0][2] * jacobianMatrix[2][1] - jacobianMatrix[0][1] * jacobianMatrix[2][2]) / fabs(determinant);
	inverseMatrix[0][2] = (jacobianMatrix[0][1] * jacobianMatrix[1][2] - jacobianMatrix[0][2] * jacobianMatrix[1][1]) / fabs(determinant);
	inverseMatrix[1][0] = (jacobianMatrix[2][0] * jacobianMatrix[1][2] - jacobianMatrix[1][0] * jacobianMatrix[2][2]) / fabs(determinant);
	inverseMatrix[1][1] = (jacobianMatrix[0][0] * jacobianMatrix[2][2] - jacobianMatrix[0][2] * jacobianMatrix[2][0]) / fabs(determinant);
	inverseMatrix[1][2] = (jacobianMatrix[1][0] * jacobianMatrix[0][2] - jacobianMatrix[0][0] * jacobianMatrix[1][2]) / fabs(determinant);
	inverseMatrix[2][0] = (jacobianMatrix[1][0] * jacobianMatrix[2][1] - jacobianMatrix[1][1] * jacobianMatrix[2][0]) / fabs(determinant);
	inverseMatrix[2][1] = (jacobianMatrix[0][1] * jacobianMatrix[2][0] - jacobianMatrix[0][0] * jacobianMatrix[2][1]) / fabs(determinant);
	inverseMatrix[2][2] = (jacobianMatrix[0][0] * jacobianMatrix[1][1] - jacobianMatrix[0][1] * jacobianMatrix[1][0]) / fabs(determinant);
}
/*
void PROBLEM_3D::calculateJacobianMatrix(double ksi, double eta, double teta, int numElement)
{
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			jacobianMatrix[i][j] = 0;
		}
	}

	auto elem = nver[numElement];

	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < EL_SIZE; j++)
		{
			jacobianMatrix[i][0] += XYZ[elem.grid.nodes[j]].x * gradPhi[i][j](ksi, eta, teta);
			jacobianMatrix[i][1] += XYZ[elem.grid.nodes[j]].y * gradPhi[i][j](ksi, eta, teta);
			jacobianMatrix[i][2] += XYZ[elem.grid.nodes[j]].z * gradPhi[i][j](ksi, eta, teta);
		}
	}
}

void PROBLEM_3D::calculateTemplateJacobianMatrix(int index, int numElement)
{
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			jacobianMatrix[i][j] = 0;
		}
	}

	auto elem = nver[numElement];

	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < EL_SIZE; j++)
		{
			jacobianMatrix[i][0] += XYZ[elem.grid.nodes[j]].x * gradPhiRes[index][i][j];
			jacobianMatrix[i][1] += XYZ[elem.grid.nodes[j]].y * gradPhiRes[index][i][j];
			jacobianMatrix[i][2] += XYZ[elem.grid.nodes[j]].z * gradPhiRes[index][i][j];
		}
	}
}

double GetIntegrandRightPart(int index, int i, int numElement)
{
	double result = 0;
	double part[3];
	vector<double> HcLocal(3, 0);

	for (int w = 0; w < 8; w++)
	{
		for (int j = 0; j < 3; j++)
		{
			HcLocal[j] += phiRes[index][w] * Hc_Ofgrid.nodes[nver[numElement].grid.nodes[w]][j];
		}
	}

	for (int k = 0; k < 3; k++)
	{
		part[k] = 0;
		for (int l = 0; l < 3; l++)
		{
			part[k] += inverseMatrix[k][l] * gradPhiRes[index][l][i];
		}

		result += part[k] * HcLocal[k];
	}

	result *= fabsl(determinant);
	return result;
}

*/