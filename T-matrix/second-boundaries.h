#pragma once

#include "stdafx.h"
#include "grid.h"

const double tetta = 10;

class SecondBoundary {
private:

    const Grid& grid;

    double b[4];

    double x[4], y[4], z[4];

    double JacobianGran2D(double ksi, double eta);
    double dr_du_t(const double* t, const double v);
    double dr_dv_t(const double* t, const double u);

    double calculateGaussIntegral(int i);

public:
    SecondBoundary(const Grid& _grid);
    ~SecondBoundary();

    void applySecondBoundaryConditions(const PLANE &);

    double getB(int);
};
