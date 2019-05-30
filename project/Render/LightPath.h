#pragma once
#include "Vector.h"
#include "Metric.h"
#include "../CurveTracer/Polynomial.h"
#include "../CurveTracer/Primitives.h"

#define D_ALPHA 1.0

// the step in normalized time
#define STEPS_PER_RADIUS 20.0

#define MAX_ALPHA_STEPS 181
#define MAX_DIST_STEPS 79
#define MAX_T_STEPS ((int)(30 * STEPS_PER_RADIUS))

extern const double dist[];
extern const char *model_path;

class LightPath {
  static double tstep;
  static double rs;
  static int get_dist_steps(const double d);
public:
  static double buffer[MAX_DIST_STEPS][MAX_ALPHA_STEPS][3];
  static void LoadModel();
  static bool MetricQuery(Ray &ray, const double dt, const double tfinal, const Metric &metric);
  static int March(Ray ray, const Metric &metric, std::vector<Primitives*> primitives, double marchDistance, double &t, Curve3d &c);
  static void BackgroundTextureLookup(Ray &ray, const Metric &metric, double &theta, double &phi);
  static Curve3d ModelQuery(Ray &ray, const Metric &metric);


};
