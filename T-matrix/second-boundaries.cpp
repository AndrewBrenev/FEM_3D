#include "second-boundaries.h"

inline double w1(double a) { return (1 - a) / 2; }
inline double w2(double a) { return (1 + a) / 2; }
inline double phi1(double ksi, double eta) { return w1(ksi) * w1(eta); }
inline double phi2(double ksi, double eta) { return w2(ksi) * w1(eta); }
inline double phi3(double ksi, double eta) { return w1(ksi) * w2(eta); }
inline double phi4(double ksi, double eta) { return w2(ksi) * w2(eta); }

static vector<function<double(double, double)>> phi = { phi1, phi2, phi3, phi4 };

double test_function(double x, double y) {
	return x + y;
}

// Number of interval splits in Gauss Integration
static const uint32_t GAUSS_FRAGMENTATION = 30;
static const int GAUSS_DEGREE = 5;
double Gauss(double (*f)(double), double begin, double end)
{
	const double Xi[GAUSS_DEGREE] = { -0.9061798,-0.5384693,0,0.5384693,0.9061798 };
	const double Ci[GAUSS_DEGREE] = { 0.4786287,0.2369269,0.5688888 ,0.2369269,0.4786287 };

	double iLength = (end - begin) / GAUSS_FRAGMENTATION;
	double ra = iLength / 2;

	double sum = 0.0;
	for (int i = 0; i < GAUSS_FRAGMENTATION; ++i)
	{
		
		double su = begin + i * iLength + ra;
		double Q, S = 0.0;
		for (int i = 0; i < GAUSS_DEGREE; i++)
		{
			Q = su + ra * Xi[i];
			S += Ci[i] * f(Q);
		}
		sum += ra * S;
	}
	return sum;
};


double SecondBoundary::calculateGaussIntegral( int _i)
{
	const double Xi[GAUSS_DEGREE] = { -0.9061798,-0.5384693,0,0.5384693,0.9061798 };
	const double Ci[GAUSS_DEGREE] = { 0.4786287,0.2369269,0.5688888 ,0.2369269,0.4786287 };

	double begin(0), end(1);


	double iLength = (end - begin) / GAUSS_FRAGMENTATION;
	double ra = iLength / 2;

	double sum = 0.0;
	for (int i = 0; i < GAUSS_FRAGMENTATION; ++i)
	{
		
		double su = begin + i * iLength + ra;

		double ksi, S = 0.0;

		for (int j = 0; j < GAUSS_DEGREE; j++)
		{
			ksi = su + ra * Xi[j];

			// Внутренний интеграл
			double inner_sum = 0.0;
			for (int k = 0; k < GAUSS_FRAGMENTATION; ++k)
			{
		
				double i_su = begin + k * iLength + ra;

				double eta, i_S = 0.0;
				for (int m = 0; m < GAUSS_DEGREE; m++)
				{
					eta = i_su + ra * Xi[m];
					// i_S += Ci[m] * test_function(ksi, eta);
					i_S += Ci[m] * tetta  * phi[i](ksi, eta) * JacobianGran2D(ksi, eta);
				}
				inner_sum += ra * i_S;
			}
			// Конец внутреннего интеграла
		
			S += Ci[j] * inner_sum;
		}
		sum += ra * S;
	}
	return sum;
};

double SecondBoundary::JacobianGran2D( double ksi, double eta)
{
    double ru[3], rv[3], g11(0.), g22(0.), g12(0.);
    ru[0] = dr_du_t(x, eta);
    ru[1] = dr_du_t(y, eta);
    ru[2] = dr_du_t(z, eta);

    rv[0] = dr_dv_t(x, ksi);
    rv[1] = dr_dv_t(y, ksi);
    rv[2] = dr_dv_t(z, ksi);

    for (int i = 0; i < 3; i++)
    {
        g11 += ru[i] * ru[i];
        g22 += rv[i] * rv[i];
        g12 += ru[i] * rv[i];
    }
    return sqrt(fabs(g11 * g22 - g12 * g12));
}

double SecondBoundary:: dr_du_t(const double* t, const double v)
{
    return (-t[0] * (1 - v) + t[1] * (1 - v) - t[2] * (1 + v) + t[3] * (1 + v)) / 4.0;
}

double SecondBoundary:: dr_dv_t(const double* t, const double u)
{
    return (-t[0] * (1 - u) - t[1] * (1 + u) + t[2] * (1 - u) + t[3] * (1 + u)) / 4.0;
}

 SecondBoundary::SecondBoundary(const Grid& _grid): grid(_grid){ }

 SecondBoundary::~SecondBoundary(){}

 double SecondBoundary::getB(int _i){
	 return (_i >= 0 && _i <= 3) ? b[_i] : 0;
 }

 void SecondBoundary::applySecondBoundaryConditions(const PLANE& _plane) {

	// double g = calculateGaussIntegral(4);

	 for (int i = 0; i < 4; i++)
	 {
		 x[i] = grid.nodes[_plane.node_id[i]].x;
		 y[i] = grid.nodes[_plane.node_id[i]].y;
		 z[i] = grid.nodes[_plane.node_id[i]].z;
	 }

	 for (int i = 0; i < 4; i++)
		 b[i] = calculateGaussIntegral(i);
 }