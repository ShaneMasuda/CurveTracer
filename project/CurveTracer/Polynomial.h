#pragma once
#include <vector>
#include "../Model/alglib/solvers.h"

class Polynomial
{
public:
  Polynomial();
  Polynomial(double* coefficients, int size);
  std::vector<double> coefficients;
};

std::vector<double> getRealRoots(Polynomial a);

double evaluate(Polynomial &a, double value);

double minimize(Polynomial &a, double value);

Polynomial operator +(Polynomial a, Polynomial b);

Polynomial operator *(double c, Polynomial a);

Polynomial operator *(Polynomial a, Polynomial b);

class Curve3d
{
public:
  Curve3d();
  Polynomial components[3];
};
