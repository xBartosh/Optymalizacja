#pragma once

#include"ode_solver.h"

#include <cmath>

#define M_PI 3.14159265359


matrix ff0T(matrix, matrix = NAN, matrix = NAN);
matrix ff0R(matrix, matrix = NAN, matrix = NAN);

matrix df0(double, matrix, matrix = NAN, matrix = NAN);
matrix df1(double t, matrix Y, matrix ud1 = NAN, matrix ud2 = NAN);
matrix df2(double t, matrix Y, matrix ud1 = NAN, matrix ud2 = NAN);

matrix ff1T(matrix x, matrix ud1 = NAN, matrix ud2 = NAN);
matrix* getSimulationData1R( matrix x, matrix ud1 = NAN, matrix ud2 = NAN);
matrix ff1R(matrix x, matrix ud1 = NAN, matrix ud2 = NAN);

double ff2T(matrix x, matrix ud1 = NAN, matrix ud2 = NAN);
double boundary(int i, double x, double a);
matrix ff2Tz(matrix x, matrix ud1 = NAN, matrix ud2 = NAN);
matrix ff2Tw(matrix x, matrix ud1 = NAN, matrix ud2 = NAN);
matrix ff2R(matrix x, matrix ud1 = NAN, matrix ud2 = NAN);


