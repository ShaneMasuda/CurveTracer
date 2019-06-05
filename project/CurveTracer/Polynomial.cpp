#include <algorithm>
#include <cmath>
#include "Polynomial.h"

Polynomial::Polynomial()
{
  coefficients = std::vector<double>();
}

Polynomial::Polynomial(double* coefficients, int size)
{
  this->coefficients = std::vector<double>(coefficients, coefficients + size);
}

std::vector<double> getRealRoots(Polynomial a)
{
  alglib::real_1d_array curve;
  curve.attach_to_ptr(a.coefficients.size(), a.coefficients.data());
  alglib::complex_1d_array roots;
  alglib::polynomialsolverreport rep;
  try
  {
    alglib::polynomialsolve(curve, a.coefficients.size() - 1, roots, rep);
  }
  catch(alglib::ap_error)
  {
    return std::vector<double>();
  }
  std::vector<double> realRoots;
  for (int i = 0; i < a.coefficients.size() - 1; i++)
  {
    if (roots[i].y == 0 && std::find(realRoots.begin(), realRoots.end(), roots[i].x) == realRoots.end())
    {
      realRoots.push_back(roots[i].x);
    }
  }
  std::sort(realRoots.begin(), realRoots.end(), std::greater<double>());
  return realRoots;
}

double evaluate(Polynomial &a, double value)
{
  if (a.coefficients.size() == 0)
  {
    return 0;
  }
  double result = a.coefficients[0];
  for (int i = 1; i < a.coefficients.size(); i++)
  {
    result += a.coefficients[i] * pow(value, i);
  }
  return result;
}

Polynomial derivative(Polynomial &a)
{
  Polynomial d;
  for (int i = 1; i < a.coefficients.size(); i++)
  {
    d.coefficients.push_back(a.coefficients[i] * (double)i);
  }
  return d;
}

double minimize(Polynomial &a, double value)
{
  if (a.coefficients.size() < 2)
  {
    return 0;
  }
  Polynomial derivative;
  for (int i = 1; i < a.coefficients.size(); i++)
  {
    derivative.coefficients.push_back(a.coefficients[i] * (double)i);
  }
  return evaluate(derivative, value);
}

Polynomial operator +(Polynomial a, Polynomial b)
{
  Polynomial result;
  if (a.coefficients.size() > b.coefficients.size())
  {
    result = a;
    for (int i = 0; i < b.coefficients.size(); i++)
    {
      result.coefficients[i] += b.coefficients[i];
    }
  }
  else
  {
    result = b;
    for (int i = 0; i < a.coefficients.size(); i++)
    {
      result.coefficients[i] += a.coefficients[i];
    }
  }
  return result;
}

Polynomial operator *(double c, Polynomial a)
{
  Polynomial result = a;
  for (int i = 0; i < result.coefficients.size(); i++)
  {
    result.coefficients[i] *= c;
  }
  return result;
}

Polynomial operator *(Polynomial a, Polynomial b)
{
  Polynomial result;
  int size = a.coefficients.size() + b.coefficients.size() - 1;
  for (int i = 0; i < size; i++)
  {
    result.coefficients.push_back(0);
  }
  for (int i = 0; i < a.coefficients.size(); i++)
  {
    for (int j = 0; j < b.coefficients.size(); j++)
    {
      result.coefficients[i + j] += a.coefficients[i] * b.coefficients[j];
    }
  }
  return result;
}

Curve3d::Curve3d()
{
  for (int i = 0; i < 3; i++)
  {
    components[i] = Polynomial();
  }
}
