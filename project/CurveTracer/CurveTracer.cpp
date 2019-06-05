#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include "CurveTracer.h"
#include "Polynomial.h"
#include "Primitives.h"
#include "Camera.h"
#include "../Model/alglib/solvers.h"
#include "../Render/LightPath.h"
#include "../Render/Timer.h"
#include "CurveTracer.h"
#include "../png++/png.hpp"

#define pi 3.141592653589793238462643383279502884

int main(int argc, char* argv[])
{
  if (argc < 2)
  {
    std::cout << "Usage: ./CurveTracer (path to scene file)" << std::endl;
    return 0;
  }
  LightPath::LoadModel();

  Camera camera;
  Metric metric;
  std::vector<Primitives*> primitives;

  // Initialize the scene
  parseScene(argv[1], camera, metric, primitives);
  SimplexNoise sn;

  png::image<png::rgb_pixel> image(camera.xres, camera.yres);
  png::image<png::rgb_pixel> background("../Render/Textures/stars2.png");
  if (argc == 3 && argv[2] == std::string("background"))
  {
    std::cout << "Generating only the background...\n";
    for (size_t y = 0; y < image.get_height(); ++y)
    {
      for (size_t x = 0; x < image.get_width(); ++x)
      {
        double dx = (double)x / image.get_width();
        double dy = (double)y / image.get_height();
        dx = 2.0 * dx - 1.0;
        dy = 2.0 * dy - 1.0;
        Point3d p = camera.center + (camera.sx * dx) + (camera.sy * dy);
        Ray ray(camera.pos, (p - camera.pos).normalize());

        double C[1] = {1.0};
        double T[2] = {0.0, 1.0};
        Polynomial constantTerm(C, 1);
        Polynomial linearTerm(T, 2);
        Curve3d c;
        c.components[0] = ray.O.x * constantTerm + ray.D.x * linearTerm;
        c.components[1] = ray.O.y * constantTerm + ray.D.y * linearTerm;
        c.components[2] = ray.O.z * constantTerm + ray.D.z * linearTerm;
        std::vector<double> bhroots = primitives[0]->intersects(c);
        if (bhroots.size() > 0)
        {
          image[y][x] = png::rgb_pixel(0, 0, 0);
          continue;
        }
        std::vector<double> roots = primitives[1]->intersects(c);
        double theta;
        double phi;
        double xi, yi, zi;

        double ti = roots[0];

        xi = evaluate(c.components[0], ti);
        yi = evaluate(c.components[1], ti);
        zi = evaluate(c.components[2], ti);

        theta = acos(zi / 3000);
        phi = atan2(yi, xi);
        double u = phi / (2.0 * pi) + 0.5;
        double v = theta / pi;

        int tx = (int)((u) * (background.get_width() - 1));
        int ty = (int)((v) * (background.get_height() - 1));
        png::rgb_pixel pixel = background[ty][tx];
        image[y][x] = pixel;
      }
    }
    image.write("straightrays.png");
    return 0;
  }

  Timer timer;
  for (size_t y = 0; y < image.get_height(); ++y)
  {
    for (size_t x = 0; x < image.get_width(); ++x)
    {
      double dx = (double)x / image.get_width();
      double dy = (double)y / image.get_height();
      dx = 2.0 * dx - 1.0;
      dy = 2.0 * dy - 1.0;
      Point3d p = camera.center + (camera.sx * dx) + (camera.sy * dy);
      Ray ray(camera.pos, (p - camera.pos).normalize());

      double r = ray.O.dist(metric.origin);
      Vector3d R = (ray.O - metric.origin) * (1.0 / r);
      double alphaMax = pi - asin((3 / 2) * sqrt(3 * (1 - metric.rs / r)) * metric.rs / r);
      double t;
      Curve3d c;
      int index;
      index = LightPath::March(ray, metric, primitives, 100, t, c);
      std::cout << "DONE\n";
      /*if (acos(ray.D * R) >= alphaMax)
      {
        index = LightPath::March(ray, metric, primitives, 100, t, c);
      }
      else
      {
        c = LightPath::ModelQuery(ray, metric);
        std::vector<std::vector<double>> rootList;
        for (Primitives* primitive : primitives)
        {
          rootList.push_back(primitive->intersects(c));
        }
        index = findFirstPrimitiveHit(rootList, c, Ray(ray.O, R * -1.0), t);
      }*/
      /* if (acos(ray.D * R) >= alphaMax && index == 1)
      {
        index = 0;
      } */
      if (index == 0)
      {
        image[y][x] = png::rgb_pixel(0, 0, 0);
      }
      else if (index == 1)
      {
        double xi, yi, zi, theta, phi;

        xi = evaluate(c.components[0], t);
        yi = evaluate(c.components[1], t);
        zi = evaluate(c.components[2], t);
        theta = acos(zi / 3000);
        phi = atan2(yi, xi);

        double u = phi / (2.0 * pi) + 0.5;
        double v = theta / pi;

        int tx = (int)((u) * (background.get_width() - 1));
        int ty = (int)((v) * (background.get_height() - 1));
        png::rgb_pixel pixel = background[ty][tx];
        image[y][x] = pixel;
      }
      else if (index == 2)
      {
        Point3d intersection(evaluate(c.components[0], t), evaluate(c.components[1], t), evaluate(c.components[2], t));
        double r = intersection.dist(Point3d(0, 0, 0));
        // double radiance = (1.0 - sqrt(metric.rs / r)) / pow(r, 3);
        image[y][x] = png::rgb_pixel(0, 255, 0);
      }
      else
      {
        image[y][x] = png::rgb_pixel(0, 0, 255);
      }
    }
  }
  timer.stop();
  image.write("test.png");
  std::cout << "Rendered image in ";
  timer.print_elapsed_time();

  return 0;
}
