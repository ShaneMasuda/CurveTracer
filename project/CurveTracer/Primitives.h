#pragma once
#include <vector>
#include "Polynomial.h"
#include "../Model/Vector.h"
#include "../Simplex/SimplexNoise.h"

class Primitives
{
public:
  double origin[3] = {0.0, 0.0, 0.0};
  virtual std::vector<double> intersects(Curve3d c) = 0;
};

class Sphere : public Primitives
{
public:
  Sphere();
  Sphere(double radius, double* origin);
  double radius;
  std::vector<double> intersects(Curve3d c);
};

class Cylinder : public Primitives
{
public:
  Cylinder();
  Cylinder(double innerRadius, double radius, double height, double* origin);
  double innerRadius;
  double radius;
  double height;
  std::vector<double> intersects(Curve3d c);
};

class Torus : public Primitives
{
public:
  Torus();
  Torus(double R, double r, double* origin);
  double R;
  double r;
  std::vector<double> intersects(Curve3d c);
};

class Disk : public Primitives
{
public:
  Disk();
  Disk(double r, double* origin);
  double r;
  std::vector<double> intersects(Curve3d c);
};

int findFirstPrimitiveHit(std::vector<std::vector<double>> rootList, Curve3d c, Ray ray, double &t);

int findIntersectedPrimitive(std::vector<Primitives*> primitives, Curve3d c, double &t, double t0, double tf);
