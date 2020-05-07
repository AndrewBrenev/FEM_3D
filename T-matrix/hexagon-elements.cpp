#include "hexagon-elements.h"


inline double w1(double a) { return (1 - a) / 2; }
inline double w2(double a) { return (1 + a) / 2; }
inline double phi1(double ksi, double eta, double teta) { return w1(ksi) * w1(eta) * w1(teta); }
inline double phi2(double ksi, double eta, double teta) { return w2(ksi) * w1(eta) * w1(teta); }
inline double phi3(double ksi, double eta, double teta) { return w1(ksi) * w2(eta) * w1(teta); }
inline double phi4(double ksi, double eta, double teta) { return w2(ksi) * w2(eta) * w1(teta); }
inline double phi5(double ksi, double eta, double teta) { return w1(ksi) * w1(eta) * w2(teta); }
inline double phi6(double ksi, double eta, double teta) { return w2(ksi) * w1(eta) * w2(teta); }
inline double phi7(double ksi, double eta, double teta) { return w1(ksi) * w2(eta) * w2(teta); }
inline double phi8(double ksi, double eta, double teta) { return w2(ksi) * w2(eta) * w2(teta); }

static vector<function<double(double, double, double)>> phi = { phi1, phi2, phi3, phi4, phi5, phi6,phi7, phi8 };

inline double PhiKsi1(double ksi, double eta, double teta) { return -w1(eta) * w1(teta) / 2; }
inline double PhiKsi2(double ksi, double eta, double teta) { return w1(eta) * w1(teta) / 2; }
inline double PhiKsi3(double ksi, double eta, double teta) { return -w2(eta) * w1(teta) / 2; }
inline double PhiKsi4(double ksi, double eta, double teta) { return w2(eta) * w1(teta) / 2; }
inline double PhiKsi5(double ksi, double eta, double teta) { return -w1(eta) * w2(teta) / 2; }
inline double PhiKsi6(double ksi, double eta, double teta) { return w1(eta) * w2(teta) / 2; }
inline double PhiKsi7(double ksi, double eta, double teta) { return -w2(eta) * w2(teta) / 2; }
inline double PhiKsi8(double ksi, double eta, double teta) { return w2(eta) * w2(teta) / 2; }
inline double PhiEta1(double ksi, double eta, double teta) { return -w1(ksi) * w1(teta) / 2; }
inline double PhiEta2(double ksi, double eta, double teta) { return -w2(ksi) * w1(teta) / 2; }
inline double PhiEta3(double ksi, double eta, double teta) { return w1(ksi) * w1(teta) / 2; }
inline double PhiEta4(double ksi, double eta, double teta) { return w2(ksi) * w1(teta) / 2; }
inline double PhiEta5(double ksi, double eta, double teta) { return -w1(ksi) * w2(teta) / 2; }
inline double PhiEta6(double ksi, double eta, double teta) { return -w2(ksi) * w2(teta) / 2; }
inline double PhiEta7(double ksi, double eta, double teta) { return w1(ksi) * w2(teta) / 2; }
inline double PhiEta8(double ksi, double eta, double teta) { return w2(ksi) * w2(teta) / 2; }
inline double PhiTeta1(double ksi, double eta, double teta) { return -w1(eta) * w1(teta) / 2; }
inline double PhiTeta2(double ksi, double eta, double teta) { return -w2(eta) * w1(teta) / 2; }
inline double PhiTeta3(double ksi, double eta, double teta) { return -w1(eta) * w2(teta) / 2; }
inline double PhiTeta4(double ksi, double eta, double teta) { return -w2(eta) * w2(teta) / 2; }
inline double PhiTeta5(double ksi, double eta, double teta) { return w1(eta) * w1(teta) / 2; }
inline double PhiTeta6(double ksi, double eta, double teta) { return w2(eta) * w1(teta) / 2; }
inline double PhiTeta7(double ksi, double eta, double teta) { return w1(eta) * w2(teta) / 2; }
inline double PhiTeta8(double ksi, double eta, double teta) { return w2(eta) * w2(teta) / 2; }

static vector<vector<function<double(double, double, double)>>> gradPhi = {
{ PhiKsi1, PhiKsi2, PhiKsi3, PhiKsi4, PhiKsi5, PhiKsi6, PhiKsi7, PhiKsi8 },
{ PhiEta1, PhiEta2, PhiEta3, PhiEta4, PhiEta5, PhiEta6, PhiEta7, PhiEta8 },
{ PhiTeta1, PhiTeta2, PhiTeta3, PhiTeta4, PhiTeta5, PhiTeta6, PhiTeta7, PhiTeta8 }
};


// Number of interval splits in Gauss Integration
static const uint32_t GAUSS_FRAGMENTATION = 30;
static const int GAUSS_DEGREE_2 = 5;
double HexagonElements::Gauss(double (*f)(double), double begin, double end)
{
	const double Xi[GAUSS_DEGREE_2] = { -0.9061798,-0.5384693,0,0.5384693,0.9061798 };
	const double Ci[GAUSS_DEGREE_2] = { 0.4786287,0.2369269,0.5688888 ,0.2369269,0.4786287 };

	double iLength = (end - begin) / GAUSS_FRAGMENTATION;
	double sum = 0.0;
	for (int i = 0; i < GAUSS_FRAGMENTATION; ++i)
	{
		double a = begin + i * iLength;
		double b = begin + (i + 1) * iLength;

		double ra = (b - a) / 2;
		double su = (a + b) / 2;
		double Q, S = 0.0;
		for (int i = 0; i < GAUSS_DEGREE_2; i++)
		{
			Q = su + ra * Xi[i];
			S += Ci[i] * f(Q);
		}
		sum += ra * S;
	}
	return sum;
};

double HexagonElements::getDeterminant()
{
	double result = 0;
	result = jacobianMatrix[0][0] * jacobianMatrix[1][1] * jacobianMatrix[2][2] -
		jacobianMatrix[0][0] * jacobianMatrix[1][2] * jacobianMatrix[2][1] -
		jacobianMatrix[0][1] * jacobianMatrix[1][0] * jacobianMatrix[2][2] +
		jacobianMatrix[0][1] * jacobianMatrix[2][0] * jacobianMatrix[1][2] +
		jacobianMatrix[1][0] * jacobianMatrix[0][2] * jacobianMatrix[2][1] -
		jacobianMatrix[0][2] * jacobianMatrix[1][1] * jacobianMatrix[2][0];

	return result;
};

void HexagonElements::calculateInverseMatrix()
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
};

void HexagonElements::calculateGradPhiAndPhiAtGaussPoints()
{
	for (int k = 0; k < 4; k++)
		for (int l = 0; l < 4; l++)
			for (int r = 0; r < 4; r++)
			{
				int index = k * 16 + l * 4 + r;
				for (int m = 0; m < 8; m++)
				{
					phiRes[index][m] = phi[m](t[k], t[l], t[r]);
					tauProduct[index] = tau[k] * tau[l] * tau[r];
					for (int j = 0; j < 3; j++)
					{
						gradPhiRes[index][j][m] = gradPhi[j][m](t[k], t[l], t[r]);
					}
				}
			}
};
double HexagonElements::getIntegrand(int index, int i, int j, int numElement)
{
	double result = 0;
	double part1[3], part2[3];
	fill_n(part1, 3, 0);
	fill_n(part2, 3, 0);
	for (int k = 0; k < 3; k++)
	{
		for (int l = 0; l < 3; l++)
		{
			part1[k] += inverseMatrix[k][l] * gradPhiRes[index][l][i];
			part2[k] += inverseMatrix[k][l] * gradPhiRes[index][l][j];
		}
		result += part1[k] * part2[k];
	}
	result *= fabs(determinant);
	return result;
}


void HexagonElements::calculateJacobianMatrix(double ksi, double eta, double teta, int numElement)
{
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			jacobianMatrix[i][j] = 0;

	auto elem = grid.elems[numElement];

	for (int i = 0; i < 3; i++)
		for (int j = 0; j < EL_SIZE; j++)
		{
			jacobianMatrix[i][0] += grid.nodes[elem.value[j]].x * gradPhi[i][j](ksi, eta, teta);
			jacobianMatrix[i][1] += grid.nodes[elem.value[j]].y * gradPhi[i][j](ksi, eta, teta);
			jacobianMatrix[i][2] += grid.nodes[elem.value[j]].z * gradPhi[i][j](ksi, eta, teta);
		}

};

void HexagonElements::calculateTemplateJacobianMatrix(int index, int numElement)
{
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			jacobianMatrix[i][j] = 0;

	auto elem = grid.elems[numElement];
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 8; j++)
		{
			jacobianMatrix[i][0] += grid.nodes[elem.value[j]].x * gradPhiRes[index][i][j];
			jacobianMatrix[i][1] += grid.nodes[elem.value[j]].y * gradPhiRes[index][i][j];
			jacobianMatrix[i][2] += grid.nodes[elem.value[j]].z * gradPhiRes[index][i][j];
		}

};
double HexagonElements::getIntegrandRightPart(int index, int i, int numElement)
{
	double result = 0;
	double part[3];
	vector<double> HcLocal(3, 0);
	for (int w = 0; w < 8; w++)
	{
		for (int j = 0; j < 3; j++)
		{
			HcLocal[j] += phiRes[index][w] * Hc_OfNodes[grid.elems[numElement].value[w]][j];
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
};


vector<double> HexagonElements::getHcUsingGauss4(NODE* fArgs, double* Jtok, int hn, NODE X)
{
	vector<double> Hc;
	Hc.resize(3, 0.0);
	NODE h;
	h = (fArgs[1] - fArgs[0]);
	double Lh = h.getLengthVector();
	double SumI[3];
	fill_n(SumI, 3, 0);
	for (int i = 0; i < 4; i++)
	{
		double SumJ[3];
		fill_n(SumJ, 3, 0);
		for (int j = 1; j < hn; j++)
		{
			NODE curX, delta;
			curX = (fArgs[j - 1] + fArgs[j] + h * t[i]) / 2.0;
			delta = X - curX;
			double nominator = sqrtl(delta.x * delta.x + delta.y * delta.y + delta.z * delta.z);
			nominator = powl(nominator, 3);
			SumJ[0] += (Jtok[1] * delta.z - Jtok[2] * delta.y) / nominator;
			SumJ[1] += (Jtok[2] * delta.x - Jtok[0] * delta.z) / nominator;
			SumJ[2] += (Jtok[0] * delta.y - Jtok[1] * delta.x) / nominator;
		}
		for (int k = 0; k < 3; k++)
			SumI[k] += SumJ[k] * tau[i];
	}
	for (int k = 0; k < 3; k++)
		Hc[k] = SumI[k] * Lh / 2.0 / (4 * M_PI);

	return Hc;
};

vector<vector<double>> HexagonElements::getHcFromOneConductor(CONDUCTOR conductor, int hn)
{
	vector<vector<double>> result;
	result.resize(grid.nodes.size());
	NODE* fArgs = new NODE[hn];
	NODE str = conductor.str;
	NODE end = conductor.end;
	NODE h;
	h = (end - str) / (hn - 1);
	for (int i = 0; i < hn; i++)
		fArgs[i] = str + h * i;

	for (int i = 0; i < grid.nodes.size(); i++)
		result[i] = getHcUsingGauss4(fArgs, conductor.Jtok, hn, grid.nodes[i]);

	return result;
};

void HexagonElements::calculateHcFromConductors()
{
	Hc_OfNodes.resize(grid.nodes.size());
	int hn = 5;
	for (int i = 0; i < grid.nodes.size(); i++)
	{
		Hc_OfNodes[i].resize(3, 0);
	}
	for (unsigned int i = 0; i < conductors.size(); i++)
	{
		vector<vector<double>> HcFromOneConductor = getHcFromOneConductor(conductors[i], hn);
		for (int j = 0; j < grid.nodes.size(); j++)
		{
			for (int k = 0; k < 3; k++)
			{
				Hc_OfNodes[j][k] += HcFromOneConductor[j][k];
			}
		}
		printf_s("%d %d\n", i, conductors.size());
	}
}

void HexagonElements::ñalculateVolumes()
{
	volumesOfNodes.resize(grid.nodes.size(), 0);
	volumesOfElements.resize(grid.elems.size(), 0);
	for (int i = 0; i < grid.elems.size(); i++)
	{
		auto element = grid.elems[i];
		double volume = 0;
		for (int index = 0; index < 64; index++)
		{
			calculateTemplateJacobianMatrix(index, i);
			volume += tauProduct[index] * fabsl(getDeterminant());
		}
		volumesOfElements[i] = volume;
		for (int j = 0; j < 8; j++)
		{
			volumesOfNodes[element.value[j]] += volumesOfElements[i];
		}
	}
}

void HexagonElements::setLocalMatrixAndRightPart()
{
	for (int i = 0; i < EL_SIZE; i++)
	{
		b[i] = 0;
		for (int j = 0; j < EL_SIZE; j++)
		{
			A[i][j] = 0;
		}
	}
}

void HexagonElements::calculateLocal(int numElement) {
	double lam = 1;
	for (int index = 0; index < 64; index++)
	{
		calculateTemplateJacobianMatrix(index, numElement);
		calculateInverseMatrix();
		determinant = getDeterminant();
		for (int i = 0; i < EL_SIZE; i++)
		{
			b[i] += tauProduct[index] * getIntegrandRightPart(index, i, numElement) * lam;
			for (int j = i; j < EL_SIZE; j++)
			{
				double sum = tauProduct[index] * getIntegrand(index, i, j, numElement) * lam;
				A[i][j] += sum;
				if (j != i)
					A[j][i] += sum;

			}
		}
	}
};