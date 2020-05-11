
#define EL_SIZE 8
#pragma once

//----------------------------
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <vector>

#include <stdexcept>
#include <iostream>
#include <set>

using namespace std;

const double LAMBDA = 3;
const double GAMMA = 4;



//u с хвостиком
inline int u(int i) { return i / 4; };
inline int my(int i) { return i % 2; };
inline int v(int i) { return (i / 2) % 2; };

// Локальные элементы
struct ELEM
{
	int id;
	int value[EL_SIZE];
	double lambda;
	double gamma;

};

struct NODE
{
	int id;
	double x, y, z;
	NODE() {
		x = 0;
		y = 0;
		z = 0;
		id = -1;
	}
	NODE(double _x, double _y, double _z) {
		x = _x;
		y = _y;
		z = _z;
		id = -1;
	}
	NODE(double _x, double _y, double _z, int _id) {
		x = _x;
		y = _y;
		z = _z;
		id = _id;
	}
	NODE operator= (NODE B)
	{
		x = B.x;
		y = B.y;
		z = B.z;
		return *this;
	}
	NODE operator+ (NODE B)
	{
		return NODE(x + B.x, y + B.y, z + B.z);
	}
	NODE operator- (NODE B)
	{
		return NODE(x - B.x, y - B.y, z - B.z);
	}
	NODE operator/ (double scale)
	{
		return NODE(x / scale, y / scale, z / scale);
	}
	NODE operator* (double scale)
	{
		return NODE(x * scale, y * scale, z * scale);
	}
	double DotProduct(NODE B)
	{
		return x * B.x + y * B.y + z * B.z;
	}
	double getLengthVector() {
		return sqrtl(x * x + y * y + z * z);
	}
};

inline double right(NODE x) {
	return 4 * (x.x * x.x + x.y * x.y + x.z * x.z) + 18;
}
inline double u_real(NODE x) {
	return x.x * x.x + x.y * x.y + x.z * x.z;
}

class MATRIX
{
public:
	int n;
	vector<double> di;
	vector<double> gg;
	vector<uint32_t> ig;
	vector<uint32_t> jg;

	virtual ~MATRIX() {
		di.clear();
		gg.clear();
		ig.clear();
		jg.clear();
	}

	void addValue(int i, int j, double value)
	{
		if (i != j) {
			if (i < j) {
				int tmp = j;
				j = i;
				i = tmp;
			}

			int k;
			for (k = ig[i]; k < ig[i + 1]; k++)
				if (jg[k] == j)
					break;
			gg[k] += value;
		}
		else
			di[i] += value;
	};

	vector<double> multOnVector(const vector<double>& x)const
	{
		vector<double> y(x.size());

		for (int i = 0; i < x.size(); ++i) {
			y[i] = di[i] * x[i];
			for (size_t k = ig[i]; k < ig[i + 1]; k++)
			{
				size_t j = jg[k];
				y[i] += gg[k] * x[j];
				y[j] += gg[k] * x[i];
			}
		}
		return y;
	}
};

struct CONDUCTOR
{
	NODE str;
	NODE end;
	double Jtok[3];
};
