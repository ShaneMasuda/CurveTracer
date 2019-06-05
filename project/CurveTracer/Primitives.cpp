#include "Primitives.h"
#include <float.h>

Sphere::Sphere()
{
  radius = 0;
}

Sphere::Sphere(double radius, double* origin)
{
  this->radius = radius;
  this->origin[0] = origin[0];
  this->origin[1] = origin[1];
  this->origin[2] = origin[2];
}

std::vector<double> Sphere::intersects(Curve3d c)
{
  double constant[1] = {1.0};
  Polynomial constantTerm(constant, 1);
  return getRealRoots((origin[0] * constantTerm + -1.0 * c.components[0]) * (origin[0] * constantTerm + -1.0 * c.components[0])
                    + (origin[1] * constantTerm + -1.0 * c.components[1]) * (origin[1] * constantTerm + -1.0 * c.components[1])
                    + (origin[2] * constantTerm + -1.0 * c.components[2]) * (origin[2] * constantTerm + -1.0 * c.components[2])
                    + -1.0 * radius * radius * constantTerm);
}

Cylinder::Cylinder()
{
  radius = 0;
  height = 0;
}

Cylinder::Cylinder(double innerRadius, double radius, double height, double* origin)
{
  this->innerRadius = innerRadius;
  this->radius = radius;
  this->height = height;
  this->origin[0] = origin[0];
  this->origin[1] = origin[1];
  this->origin[2] = origin[2];
}

std::vector<double> Cylinder::intersects(Curve3d c)
{
  double constant[1] = {1.0};
  Polynomial constantTerm(constant, 1);
  std::vector<double> potentialRoots;
  std::vector<double> outerEdge = getRealRoots((origin[0] * constantTerm + -1.0 * c.components[0]) * (origin[0] * constantTerm + -1.0 * c.components[0])
                                                  + (origin[2] * constantTerm + -1.0 * c.components[2]) * (origin[2] * constantTerm + -1.0 * c.components[2])
                                                  + -1.0 * radius * radius * constantTerm);
  for (double root: outerEdge)
  {
    if (evaluate(c.components[1], root) <= height / 2.0 && evaluate(c.components[1], root) >= -height / 2.0)
    {
      potentialRoots.push_back(root);
    }
  }
  std::vector<double> top = getRealRoots(c.components[1] + (-1.0 * height / 2.0) * constantTerm);
  for (double root : top)
  {
    double x = evaluate(c.components[0], root);
    double z = evaluate(c.components[2], root);
    double r = sqrt(x*x + z*z);
    if (r >= innerRadius && r <= radius)
    {
      potentialRoots.push_back(root);
    }
  }
  std::vector<double> bottom = getRealRoots(c.components[1] + (height / 2.0) * constantTerm);
  for (double root : bottom)
  {
    double x = evaluate(c.components[0], root);
    double z = evaluate(c.components[2], root);
    double r = sqrt(x*x + z*z);
    if (r >= innerRadius && r <= radius)
    {
      potentialRoots.push_back(root);
    }
  }
  std::sort(potentialRoots.begin(), potentialRoots.end(), std::greater<double>());
  return potentialRoots;
}

Torus::Torus()
{
  R = 0;
  r = 0;
}

Torus::Torus(double R, double r, double* origin)
{
  this->R = R;
  this->r = r;
  this->origin[0] = origin[0];
  this->origin[1] = origin[1];
  this->origin[2] = origin[2];
}

std::vector<double> Torus::intersects(Curve3d c)
{
  double constant[1] = {1.0};
  Polynomial constantTerm(constant, 1);
  return getRealRoots((c.components[0] * c.components[0] + c.components[1] * c.components[1] + c.components[2] * c.components[2] + R * R * constantTerm + -1.0 * r * r * constantTerm) *
                      (c.components[0] * c.components[0] + c.components[1] * c.components[1] + c.components[2] * c.components[2] + R * R * constantTerm + -1.0 * r * r * constantTerm) +
                      -4.0 * R * R * (c.components[0] * c.components[0] + c.components[2] * c.components[2]));
}

Disk::Disk()
{
  r = 0;
}

Disk::Disk(double r, double* origin)
{
  this->r = r;
  this->origin[0] = origin[0];
  this->origin[1] = origin[1];
  this->origin[2] = origin[2];
}

std::vector<double> Disk::intersects(Curve3d c)
{
  std::vector<double> planeRoots = getRealRoots(c.components[1]);
  std::vector<double> roots;
  for (double root : planeRoots)
  {

    double x = evaluate(c.components[0], root);
    double z = evaluate(c.components[2], root);
    if (sqrt(x*x + z*z) < r)
    {
      roots.push_back(root);
    }
  }
  return roots;
}

int findFirstPrimitiveHit(std::vector<std::vector<double>> rootList, Curve3d c, Ray ray, double &t)
{
  int index = -1;
  double smallestPositiveDistance = DBL_MAX;
  for (int i = 0; i < rootList.size(); i++)
  {
    for (double root : rootList[i])
    {
      Point3d p(evaluate(c.components[0], root), evaluate(c.components[1], root), evaluate(c.components[2], root));
      double distance = ray.D * (p - ray.O);
      if (distance > 0 && distance < smallestPositiveDistance)
      {
          t = root;
          index = i;
          smallestPositiveDistance = distance;
      }
    }
  }
  return index;
}

int findIntersectedPrimitive(std::vector<Primitives*> primitives, Curve3d c, double &t, double t0, double tf)
{
  std::vector<std::vector<double>> rootList;
  for (Primitives* primitive : primitives)
  {
    rootList.push_back(primitive->intersects(c));
  }
  int index = -1;
  for (int i = 0; i < rootList.size(); i++)
  {
    for (double root : rootList[i])
    {
      if (root >= t0 && root <= tf)
      {
        t = root;
        index = i;
      }
    }
  }
  return index;
}
