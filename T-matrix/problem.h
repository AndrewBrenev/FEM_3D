 #pragma once

#include "stdafx.h"
#include "Solver.h"
#include "Grid.h"
#include "hexagon-elements.h"
#include "second-boundaries.h"
#include "T-Matrix.h"

const string t_matrix_config = "../input/Pipe/t-matrix/params.txt";
const string t_matrix_ig = "../input/Pipe/t-matrix/ig.txt";
const string t_matrix_jg = "../input/Pipe/t-matrix/jg.txt";
const string t_matrix_gg = "../input/Pipe/t-matrix/gg.txt";

const string mesh_inftry_path = "../input/Pipe/mesh/inftry.dat";
const string mesh_xyz_path = "../input/Pipe/mesh/xyz.dat";
const string mesh_nver_path = "../input/Pipe/mesh/nver.dat";

const string first_boundary_path = "../input/Pipe/first_boundary.txt";
const string second_boundary_path = "../input/Test/second_boundary.txt";



class PROBLEM_3D
{
private:

	Grid grid;
	HexagonElements mesh;

	//in CSR format
	MATRIX slae, M;
	SymmetricSolver solver;

	const bool inconsistentProblem;
	// in CSC format
	T_MATRIX t_matix;
	
	vector<double> f;
	vector<double> x;
	
	void setSlaeExactValue(uint32_t, double, MATRIX&);
	void addLocalMatrix(int);
	void destruct();
	void setMatrixSize(const size_t size);
	void processTerminalNode(int, int, double);
	void processTerminalNode(int, double);

	vector<uint32_t>  readFirstBoundaryFromFile(const char*);
	vector<PLANE> readSecondBoundaryFromFile(const char* input);
public:

	PROBLEM_3D(bool _inconsistentProblem) :inconsistentProblem(_inconsistentProblem), mesh(grid){

		if (inconsistentProblem)
			t_matix.readFromFiles(t_matrix_config.c_str(), t_matrix_ig.c_str(), t_matrix_jg.c_str(), t_matrix_gg.c_str());

	}; 
	PROBLEM_3D(bool _inconsistentProblem,T_MATRIX _t_matrix) :inconsistentProblem(_inconsistentProblem),t_matix(_t_matrix) ,mesh(grid){};

	~PROBLEM_3D();

	void buildPortrait();
	void buildSetBasedPortrait();
	void buildMatrixAndRightPart();

	void readGridFromFiles(const char* inf, const char* xyz, const char* nver);

	void applyFirstBoundaryConditions();
	
	void applySecondBoundaryConditions();

	void printResult(const char* fname);

	void showResault();
	void solveMatrix();

};
