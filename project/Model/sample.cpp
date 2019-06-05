#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include "Metric.h"
#include "Vector.h"
#include "LightPath.h"
#include "Timer.h"
#include "alglib/interpolation.h"

void sample();

const double dt = 5.0;
const double rs = 100.0;
const double jump_dist = 0.05;
const int INC_PER_STEP = 100;

static const Metric metric = Metric(100.0, Point3d());

struct Curve
{
  Curve()
  {
    numberOfPoints = 0;
  }
  Curve(int size)
  {
    numberOfPoints = size;
    xPoints = new double[size];
    yPoints = new double[size];
  }

  int numberOfPoints;
  double* xPoints;
  double* yPoints;
};

std::string doubleToString(double d)
{
  std::ostringstream oss;
  oss << std::fixed;
  oss << d;
  return oss.str();
}

Curve lookupTable[MAX_DIST_STEPS][MAX_ALPHA_STEPS];

double interpolatedCurves[MAX_DIST_STEPS][MAX_ALPHA_STEPS][3];

int main(int argc, char **argv) {

  sample();

  if (argc == 3)
  {
    int r = atoi(argv[1]);
    int alpha = atoi(argv[2]);
    std::ofstream plot;
    plot.open("./Data/plot.dat");
    for (int i = 0; i < lookupTable[r][alpha].numberOfPoints; i++)
    {
      plot << lookupTable[r][alpha].xPoints[i] << ' ' << lookupTable[r][alpha].yPoints[i] << std::endl;
    }
    std::string gnuplotFunction = "f(x)=" + doubleToString(interpolatedCurves[r][alpha][0]);

    for (int i = 1; i < 3; i++)
    {
      gnuplotFunction += "+" + doubleToString(interpolatedCurves[r][alpha][i]) + "*x**" + std::to_string(i);
    }
    gnuplotFunction += ';';
    plot.close();

    std::string gnuplotCommand = "gnuplot -e \"" + gnuplotFunction + " plot [-700:700] f(x), 'Data/plot.dat'; pause -1\"";
    std::cout << gnuplotCommand << std::endl;
    system(gnuplotCommand.c_str());
  }

  return 0;
}

void sample() {

  std::ofstream data_file(model_path, std::ios::binary);

  Ray ray;

  unsigned short input[3];
  double output[3];

  Timer t;

  for (int d_step = 0; d_step < MAX_DIST_STEPS; d_step++) {
    for (int alpha_step = 0; alpha_step < MAX_ALPHA_STEPS; alpha_step++) {
      double alpha = alpha_step * D_ALPHA;
      double cos_alpha = cos(DEG_2_RAD(alpha));
      double sin_alpha = sin(DEG_2_RAD(alpha));
      double d = dist[d_step];
      ray.D = Vector3d(cos_alpha, sin_alpha, 0.0);
      ray.O = Vector3d(rs * d, 0.0, 0.0);

      int tmax = MIN(MAX_T_STEPS, (int)((d - 1.0) * STEPS_PER_RADIUS));
      lookupTable[d_step][alpha_step] = Curve(tmax + 1);
      double new_d, beta, gamma;
      for (int t_step = 0; t_step <= tmax; t_step++) {

        new_d = ray.O.dist(metric.origin) / rs;

        beta = atan2(ray.O.y, ray.O.x);
        beta += beta < 0.0 ? 2.0 * M_PI : 0.0;

        gamma = atan2(ray.D.y, ray.D.x);
        gamma += gamma < 0.0 ? 2.0 * M_PI : 0.0;

        input[0] = d_step;
        input[1] = alpha_step;
        input[2] = t_step;
        output[0] = new_d;
        output[1] = beta;
        output[2] = gamma;

        lookupTable[d_step][alpha_step].xPoints[t_step] = ray.O.x;
        lookupTable[d_step][alpha_step].yPoints[t_step] = ray.O.y;

        for (int i = 0; i < INC_PER_STEP; i++) {
          metric.step(jump_dist, ray);
        }
      }
    }
  }

  for (int d_step = 0; d_step < MAX_DIST_STEPS; d_step++) {
    for (int alpha_step = 0; alpha_step < MAX_ALPHA_STEPS; alpha_step++) {
      alglib::real_1d_array xPoints;
      alglib::real_1d_array yPoints;
      xPoints.attach_to_ptr(lookupTable[d_step][alpha_step].numberOfPoints, lookupTable[d_step][alpha_step].xPoints);
      yPoints.attach_to_ptr(lookupTable[d_step][alpha_step].numberOfPoints, lookupTable[d_step][alpha_step].yPoints);
      alglib::ae_int_t info;
      alglib::barycentricinterpolant b;
      alglib::polynomialfitreport report;

      alglib::polynomialfit(xPoints, yPoints, lookupTable[d_step][alpha_step].numberOfPoints, alglib::ae_int_t(3), info, b, report);
      alglib::real_1d_array curve;
      alglib::xparams degree;
      alglib::polynomialbar2pow(b, curve, degree);

      double* interpolatedCurve = curve.getcontent();

      interpolatedCurves[d_step][alpha_step][0] = interpolatedCurve[0];
      interpolatedCurves[d_step][alpha_step][1] = interpolatedCurve[1];
      interpolatedCurves[d_step][alpha_step][2] = interpolatedCurve[2];
    }
  }

  data_file.write(reinterpret_cast<char*>(interpolatedCurves), sizeof(interpolatedCurves));

  t.stop();
  data_file.close();
  std::cout << "Generated light paths in ";
  t.print_elapsed_time();
}
