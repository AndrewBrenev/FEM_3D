#pragma once
#define  GAUSS_DEGREE  4
#define JACOBIAN_SIZE  3

#define _USE_MATH_DEFINES
#include "stdafx.h"
#include "grid.h"
#include <cmath>
#include "functional"


class HexagonElements {
private:

	const Grid& grid;

	double determinant;
	double gradPhiRes[64][JACOBIAN_SIZE][EL_SIZE];
	double phiRes[64][EL_SIZE];
	double tauProduct[64];
	double jacobianMatrix[JACOBIAN_SIZE][JACOBIAN_SIZE];
	double inverseMatrix[JACOBIAN_SIZE][JACOBIAN_SIZE];

	vector<vector<double>> Hc_OfNodes;
	vector<CONDUCTOR> conductors;
	vector<double> volumesOfElements;
	vector<double> volumesOfNodes;

	double A[EL_SIZE][EL_SIZE], b[EL_SIZE];
	double t[GAUSS_DEGREE], tau[GAUSS_DEGREE];
	
	double Gauss(double (*f)(double), double begin, double end);

	// Jacobian && local matrix
	double getDeterminant();
	double getIntegrand(int index, int i, int j, int numElement);
	double getIntegrandRightPart(int index, int i, int numElement);
	void calculateJacobianMatrix(double ksi, double eta, double teta, int numElement);
	void calculateTemplateJacobianMatrix(int index, int numElement);
	void calculateInverseMatrix();
	void calculateGradPhiAndPhiAtGaussPoints();

	vector<vector<double>> getHcFromOneConductor(CONDUCTOR conductor, int hn);
	vector<double> getHcUsingGauss4(NODE* fArgs, double* Jtok, int hn, NODE X);

public:
	HexagonElements(const Grid& _grid) :grid(_grid) {
		t[0] = -0.8611363115940525752;
		t[1] = 0.3399810435848562648;
		t[2] = -0.3399810435848562648;
		t[3] = 0.8611363115940525752;
		tau[0] = 0.3478548451374538573;
		tau[1] = 0.6521451548625461426;
		tau[2] = 0.6521451548625461426;
		tau[3] = 0.3478548451374538573;
	};

	void ñalculateVolumes();
	void calculateHcFromConductors();
	void calculateLocal(int elem_id);
	void setLocalMatrixAndRightPart();

};
