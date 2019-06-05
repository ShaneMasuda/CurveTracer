#include <math.h>
#include <iostream>
#include <assert.h>
#include <fstream>
#include "LightPath.h"

const double dist[] = {
            30.0, 29.0, 28.0, 27.0, 26.0, 25.0, 24.0, 23.0, 22.0, 21.0,
            20.0, 19.0, 18.0, 17.0, 16.0, 15.0, 14.0, 13.0, 12.0, 11.0,
            10.0, 9.5, 9.0, 8.5, 8.0, 7.5, 7.0, 6.75, 6.5, 6.25,
            6.0, 5.75, 5.5, 5.25, 5.0, 4.8, 4.6, 4.4, 4.2, 4.0,
            3.9, 3.8, 3.7, 3.6, 3.5, 3.4, 3.3, 3.2, 3.1, 3.0,
            2.9, 2.8, 2.7, 2.6, 2.5, 2.4, 2.3, 2.2, 2.1, 2.0,
            1.95, 1.9, 1.85, 1.8, 1.75, 1.7, 1.65, 1.6, 1.55, 1.5,
            1.45, 1.4, 1.35, 1.3, 1.25, 1.2, 1.15, 1.1, 1.05};

const char *model_path = "../Model/Data/model.dat";

double LightPath::buffer[MAX_DIST_STEPS][MAX_ALPHA_STEPS][3] = {{{0}}};

int LightPath::get_dist_steps(const double d) {
  int steps = 0;

  // binary search
  int min = 0;
  int max = MAX_DIST_STEPS - 1;
  do {
    steps = (min + max) / 2;
    if (dist[steps+1] > d)
      min = steps;
    else if (dist[steps] < d)
      max = steps;
    else
      break;

  } while (true);

  return steps;
}

void LightPath::LoadModel() {
  std::ifstream data_file(model_path, std::ios::binary);
  data_file.read(reinterpret_cast<char*>(buffer), 8 * MAX_DIST_STEPS * MAX_ALPHA_STEPS * 3);
  data_file.close();
}

Curve3d LightPath::ModelQuery(Ray &ray, const Metric &metric) {

  const double r = ray.O.dist(metric.origin);
  const Vector3d X = (ray.O - metric.origin) * (1.0 / r);

  const double d = r / metric.rs;
  const int d_step = get_dist_steps(d);
  const double alpha = RAD_2_DEG(acos(ray.D * X));
  const int alpha_step = (int)(alpha / D_ALPHA);

  assert(d_step >= 0 && d_step < MAX_DIST_STEPS - 1);
  if (alpha_step < 0 || alpha_step >= MAX_ALPHA_STEPS - 1)
  {
    std::cout << alpha_step << ' ' << MAX_ALPHA_STEPS << '\n';
  }
  assert(alpha_step >= 0 && alpha_step < MAX_ALPHA_STEPS - 1);

  const double alpha_weight = 1.0 - (alpha - (alpha_step * D_ALPHA)) / D_ALPHA;
  const double d_weight = 1.0 - (dist[d_step] - d) / (dist[d_step] - dist[d_step+1]);

  const double w1 = d_weight * alpha_weight;
  const double w2 = (1.0 - d_weight) * alpha_weight;
  const double w3 = d_weight * (1.0 - alpha_weight);
  const double w4 = (1.0 - d_weight) * (1.0 - alpha_weight);

  Polynomial lerp = w1 * Polynomial(buffer[d_step][alpha_step], 3)
        +  w2 * Polynomial(buffer[d_step+1][alpha_step], 3)
        +  w3 * Polynomial(buffer[d_step][alpha_step+1], 3)
        +  w4 * Polynomial(buffer[d_step+1][alpha_step+1], 3);

  const Vector3d Y = (ray.D - X * (X * ray.D)).normalize();

  double c[1] = {1.0};
  double t[2] = {0.0, 1.0};
  Polynomial constantTerm(c, 1);
  Polynomial linearTerm(t, 2);

  Curve3d result;
  result.components[0] = X.x * linearTerm + Y.x * lerp;
  result.components[1] = X.y * linearTerm + Y.y * lerp;
  result.components[2] = X.z * linearTerm + Y.z * lerp;

  return result;
}

/* int LightPath::March(Ray ray, const Metric &metric, std::vector<Primitives*> primitives, double marchDistance, double &t, Curve3d &c)
{
  int index = -1;
  int counter = 0;
  double r;
  Vector3d R;
  while (counter < 30)
  {
    r = ray.O.dist(metric.origin);
    if (r < 1.1 * metric.rs)
    {
      return 0;
    }
    R = (ray.O - metric.origin) * (1.0 / r);
    c = ModelQuery(ray, metric);
    std::cout << "Direction: " << ray.O.x << '+'<< ray.D.x << "t," << ray.O.y << '+'<< ray.D.y << "t," << ray.O.z << '+'<< ray.D.z << "t\n";
    std::cout << "x: " << c.components[0].coefficients[0] << '+' << c.components[0].coefficients[1] << "t+" << c.components[0].coefficients[2] << "t^2\n";
    std::cout << "y: " << c.components[1].coefficients[0] << '+' << c.components[1].coefficients[1] << "t+" << c.components[1].coefficients[2] << "t^2\n";
    std::cout << "z: " << c.components[2].coefficients[0] << '+' << c.components[2].coefficients[1] << "t+" << c.components[2].coefficients[2] << "t^2\n";
    std::vector<std::vector<double>> rootList;
    for (Primitives* primitive : primitives)
    {
      rootList.push_back(primitive->intersects(c));
    }

    // Find the t value where the curve is highest
    Vector3d Y = (ray.D - R * (R * ray.D)).normalize();
    std::vector<double> minima = getRealRoots(Y.x * derivative(c.components[0]) + Y.y * derivative(c.components[1]) + Y.z * derivative(c.components[2]));
    assert(minima.size() == 1);
    double minimum = minima[0];

    double scale = ray.D * R;
    double sign = scale > 0 ? 1.0 : -1.0;
    if (r < 6.0 * metric.rs)
    {
      sign = -1.0;
    }
    index = findFirstPrimitiveHit(rootList, c, Ray(ray.O, R * sign), t);
    if (sign > 0)
    {
      std::cout << "t: " << r << ' ' << 2*r << '\n';
      if (t >= r && t <= 2.0 * r)
      {
        return index;
      }
      Point3d O(evaluate(c.components[0], 2.0 * r), evaluate(c.components[1], 2.0 * r), evaluate(c.components[2], 2.0 * r));
      Point3d D(minimize(c.components[0], 2.0 * r), minimize(c.components[1], 2.0 * r), minimize(c.components[2], 2.0 * r));
      ray = Ray(O, (D * (minimum > 2.0 * r ? 1.0 : -1.0)).normalize());
      counter++;
    }
    else
    {
      std::cout << "t: " << r << ' ' << 100 << '\n';
      if (t >= 100 && t <= r)
      {
        return index;
      }
      Point3d O(evaluate(c.components[0], 100), evaluate(c.components[1], 100), evaluate(c.components[2], 100));
      Point3d D(minimize(c.components[0], 100), minimize(c.components[1], 100), minimize(c.components[2], 100));
      ray = Ray(O, (D * (minimum < 100? 1.0 : -1.0)).normalize());
      counter++;
    }
  }
  // std::cout << "Too many iterations\n";
  return index;
} */

int LightPath::March(Ray ray, const Metric &metric, std::vector<Primitives*> primitives, double marchDistance, double &t, Curve3d &c)
{
  int index = -1;
  int counter = 0;
  double r;
  Vector3d R;
  while (counter < 30)
  {
    r = ray.O.dist(metric.origin);
    if (r < 1.1 * metric.rs)
    {
      return 0;
    }
    R = (ray.O - metric.origin) * (1.0 / r);
    c = ModelQuery(ray, metric);
    std::cout << "Direction: " << ray.O.x << '+'<< ray.D.x << "t," << ray.O.y << '+'<< ray.D.y << "t," << ray.O.z << '+'<< ray.D.z << "t\n";
    std::cout << "x: " << c.components[0].coefficients[0] << '+' << c.components[0].coefficients[1] << "t+" << c.components[0].coefficients[2] << "t^2\n";
    std::cout << "y: " << c.components[1].coefficients[0] << '+' << c.components[1].coefficients[1] << "t+" << c.components[1].coefficients[2] << "t^2\n";
    std::cout << "z: " << c.components[2].coefficients[0] << '+' << c.components[2].coefficients[1] << "t+" << c.components[2].coefficients[2] << "t^2\n";

    // Find the t value where the curve is highest
    Vector3d Y = (ray.D - R * (R * ray.D)).normalize();
    std::vector<double> minima = getRealRoots(Y.x * derivative(c.components[0]) + Y.y * derivative(c.components[1]) + Y.z * derivative(c.components[2]));
    assert(minima.size() == 1);
    double minimum = minima[0];

    double emissionAngle = acos(ray.D * R);
    std::cout << emissionAngle << '\n';

    if (emissionAngle < M_PI / 2.0)
    {
      std::cout << "t: " << r << ' ' << 2*r << '\n';
      index = findIntersectedPrimitive(primitives, c, t, r, 2.0 * r);
      if (index != -1)
      {
        return index;
      }
      Point3d O(evaluate(c.components[0], 2.0 * r), evaluate(c.components[1], 2.0 * r), evaluate(c.components[2], 2.0 * r));
      Point3d D(minimize(c.components[0], 2.0 * r), minimize(c.components[1], 2.0 * r), minimize(c.components[2], 2.0 * r));
      ray = Ray(O, (D * (minimum > 2.0 * r ? -1.0 : 1.0)).normalize());
      counter++;
    }
    else
    {
      std::cout << "t: " << 100 << ' ' << r << '\n';
      index = findIntersectedPrimitive(primitives, c, t, 100, r);
      if (index != -1)
      {
        return index;
      }
      Point3d O(evaluate(c.components[0], 100), evaluate(c.components[1], 100), evaluate(c.components[2], 100));
      Point3d D(minimize(c.components[0], 100), minimize(c.components[1], 100), minimize(c.components[2], 100));
      ray = Ray(O, (D * (minimum < 100 ? -1.0 : 1.0)).normalize());
      counter++;
    }
  }
  return index;
}

void LightPath::BackgroundTextureLookup(Ray &ray, const Metric &metric, double &theta, double &phi)
{
  /* const double r = ray.O.dist(metric.origin);
  const Vector3d X = (ray.O - metric.origin) * (1.0 / r);

  const double d = r / metric.rs;
  const int d_step = get_dist_steps(d);
  const double alpha = RAD_2_DEG(acos(ray.D * X));
  const int alpha_step = (int)(alpha / D_ALPHA);

  const double alpha_weight = 1.0 - (alpha - (alpha_step * D_ALPHA)) / D_ALPHA;
  const double d_weight = 1.0 - (dist[d_step] - d) / (dist[d_step] - dist[d_step+1]);

  const double w1 = d_weight * alpha_weight;
  const double w2 = (1.0 - d_weight) * alpha_weight;
  const double w3 = d_weight * (1.0 - alpha_weight);
  const double w4 = (1.0 - d_weight) * (1.0 - alpha_weight);

  Polynomial lerp = w1 * Polynomial(buffer[d_step][alpha_step], 3)
        +  w2 * Polynomial(buffer[d_step+1][alpha_step], 3)
        +  w3 * Polynomial(buffer[d_step][alpha_step+1], 3)
        +  w4 * Polynomial(buffer[d_step+1][alpha_step+1], 3);

  const Vector3d Y = (ray.D - X * (X * ray.D)).normalize();

  double T[2] = {0.0, 1.0};
  Polynomial linearTerm(T, 2);

  Polynomial rSquared = (linearTerm * linearTerm) + (lerp * lerp);
  double r3 = minimize(rSquared);

  if (r3 <= metric.rs)
  {
    theta = 0;
    phi = 0;
    return;
  }

  double angularDeflection = 2.0 * metric.rs / (r3 * sqrt(r3 / (r3 - metric.rs)));
  double xo = r3;
  double yo = evaluate(lerp, r3);
  Vector3d direction(xo - r, yo, 0);
  direction.x = direction.x * cos(angularDeflection) - direction.y * sin(angularDeflection);
  direction.y = direction.x * sin(angularDeflection) + direction.y * cos(angularDeflection);
  direction.normalize();
  double x = direction.x;
  double y = direction.y;
  double t = (-1.0 * (x * xo + y * yo) + sqrt(2.0 * x * xo * y * yo - x * x * yo * yo - y * y * xo * xo + (x * x + y * y) * 1000000000 * 1000000000)) / (x * x + y * y);
  double xi = xo + t * x;
  double yi = yo + t * y;
  Point3d wsPoint(X.x * xi + Y.x * yi, X.y * xi + Y.y * yi, X.z * xi + Y.z * yi);
  theta = acos(wsPoint.z / 1000000000);
  phi = atan(wsPoint.y / wsPoint.x); */
}

bool LightPath::MetricQuery(Ray &ray, const double dt, const double tfinal, const Metric &metric) {
  double curr_time = 0.0;
  for (; curr_time <= tfinal - dt; curr_time += dt) {
    if (!metric.step(dt, ray))
      return false;
  }

  if (curr_time < tfinal)
    return metric.step(tfinal - curr_time, ray);

  return true;
}
