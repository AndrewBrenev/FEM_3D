#pragma once

#define EL_SIZE 8

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

// Number of interval splits in Gauss Integration
const uint32_t GAUSS_FRAGMENTATION = 30;
const int GAUSS_DEGREE = 5;

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
	// 0 - x ,1 - y , 2 - z
	double value[3];
	int id;

	NODE() {
		value[0] = 0;
		value[1] = 0;
		value[2] = 0;
		id = -1;
	}
	NODE(double _x, double _y, double _z, int _id) {
		value[0] = _x;
		value[1] = _y;
		value[2] = _z;
		id = _id;
	}
};

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

