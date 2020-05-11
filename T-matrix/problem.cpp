#include "problem.h"

// -div( lambda grad(u)) + gamma*u  = right



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
void PROBLEM_3D::addLocalMatrix(int local_id) {

	for (int i = 0; i < EL_SIZE; ++i) {

		for (int j = 0; j <= i; ++j) {
			double value = mesh.getMatrixElement(i,j);
			if (inconsistentProblem)
				processTerminalNode(grid.elems[local_id].value[j], grid.elems[local_id].value[i], value);
			else
				slae.addValue(grid.elems[local_id].value[j], grid.elems[local_id].value[i], value);

		}

		/*
		double r_value  = mesh.getRightPartRow(i)
		if (inconsistentProblem)
			//Add elem to matrix
			processTerminalNode(grid.elems[k].value[i], r_value);
		else
			// add right part
			f[grid.elems[k].value[i]] += r_value;
			*/

	}
}

void PROBLEM_3D::buildMatrixAndRightPart() {
	for (int i = 0; i < grid.elems.size(); i++)
	{
		mesh.setLocalMatrixAndRightPart();
		mesh.calculateLocal(i);
		addLocalMatrix(i);
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

void PROBLEM_3D::applyFirstBoundaryConditions()
{
	grid.firstBoundary = readFirstBoundaryFromFile(boundary_path.c_str());
	for (auto it = grid.firstBoundary.begin(); it < grid.firstBoundary.end(); it++)
		setSlaeExactValue(*it, u_real(grid.nodes[*it]), slae);
}

void PROBLEM_3D::applySecondBoundaryConditions(){}

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
		real_solution = u_real(grid.nodes[i]);
		printf("%.3lf %.3lf %.3lf %.15lf %.15lf %.15lf\n", grid.nodes[i].x, grid.nodes[i].y, grid.nodes[i].z, real_solution, x[i], abs(real_solution - x[i]));
	}
}
//--------------------------------------
void PROBLEM_3D::printResult(const char* fname)
{
	FILE* f = fopen(fname, "w");
	for (int i = 0; i < slae.n; i++)
	{
		double real_solution = u_real(grid.nodes[i]);

		fprintf(f, "%.3lf %.3lf %.3lf %.15lf %.15lf %.15lf\n", grid.nodes[i].x, grid.nodes[i].y, grid.nodes[i].z, real_solution, x[i], abs(real_solution - x[i]));
	}
	fclose(f);
}
