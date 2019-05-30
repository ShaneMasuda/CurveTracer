#pragma once
#include <vector>
#include <fstream>
#include "../Model/Metric.h"
#include "Primitives.h"


struct Camera {
  Point3d pos, center;
  Vector3d dir, up, right, sx, sy;
  int xres, yres;
  int xsamples = 1;
  int ysamples = 1;
  double focalDist, aspect, verticleFov;

  Camera()
  {
  }

  Camera(const Point3d &pos, const Vector3d &dir, const Vector3d &up,
    const int xres, const int yres, const double focalDist, const double verticleFov)
  {
      this->pos = pos;
      this->dir = dir;
      this->up = up;
      this->xres = xres;
      this->yres = yres;
      this->focalDist = focalDist;
      this->aspect = (double)xres / yres;
      this->verticleFov = verticleFov;

      this->normalize();
  }

  void normalize()
  {
    this->dir = this->dir.normalize();
    this->right = this->dir.cross(this->up).normalize();
    this->up = this->right.cross(this->dir).normalize();

    this->center = this->pos + this->dir * focalDist;
    this->sy = this->up * tan(DEG_2_RAD(this->verticleFov) / 2.0) * focalDist;
    this->sx = this->right * tan(DEG_2_RAD(this->verticleFov) / 2.0) * focalDist * this->aspect;
  }
};

void parseScene(const char* sceneFile, Camera &camera, Metric &metric, std::vector<Primitives*> &primitives)
{
  std::ifstream scene(sceneFile);
  if (scene.is_open())
  {
    std::string sceneElement;
    while (scene >> sceneElement)
    {
      if (sceneElement == "blackhole")
      {
        double schwarzchildRadius;
        Vector3d blackHoleOrigin;
        scene >> schwarzchildRadius >> blackHoleOrigin.x >> blackHoleOrigin.y >> blackHoleOrigin.z;
        double o[3] = {blackHoleOrigin.x, blackHoleOrigin.y, blackHoleOrigin.z};
        primitives.push_back(new Sphere(schwarzchildRadius, o));
        metric = Metric(schwarzchildRadius, blackHoleOrigin);
      }
      else if (sceneElement == "camera")
      {
        Vector3d origin;
        Vector3d direction;
        Vector3d up;
        int xRes;
        int yRes;
        double focalDist;
        double verticalFov;
        scene >> origin.x >> origin.y >> origin.z >> direction.x >> direction.y >> direction.z
              >> up.x >> up.y >> up.z >> xRes >> yRes >> focalDist >> verticalFov;
        camera = Camera(origin, direction, up, xRes, yRes, focalDist, verticalFov);
      }
      else if (sceneElement == "sphere")
      {
        double radius;
        double origin[3];
        scene >> radius >> origin[0] >> origin[1] >> origin[2];
        primitives.push_back(new Sphere(radius, origin));
      }
      else if (sceneElement == "cylinder")
      {
        double innerRadius;
        double radius;
        double height;
        double origin[3];
        scene >> innerRadius >> radius >> height >> origin[0] >> origin[1] >> origin[2];
        primitives.push_back(new Cylinder(innerRadius, radius, height, origin));
      }
      else if (sceneElement == "torus")
      {
        double R;
        double r;
        double origin[3];
        scene >> R >> r >> origin[0] >> origin[1] >> origin[2];
        primitives.push_back(new Torus(R, r, origin));
      }
      else if (sceneElement == "disk")
      {
        double r;
        double origin[3];
        scene >> r >> origin[0] >> origin[1] >> origin[2];
        primitives.push_back(new Disk(r, origin));
      }
    }
  }
}
