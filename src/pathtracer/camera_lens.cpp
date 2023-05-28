#include "camera.h"

#include <iostream>
#include <sstream>
#include <fstream>

#include "CGL/misc.h"
#include "CGL/vector2D.h"
#include "CGL/vector3D.h"

using std::cout;
using std::endl;
using std::max;
using std::min;
using std::ifstream;
using std::ofstream;

namespace CGL {

using Collada::CameraInfo;

Ray Camera::generate_ray_for_thin_lens(double x, double y, double rndR, double rndTheta) const {

  // TODO Assignment 7: Part 4
  // compute position and direction of ray from the input sensor sample coordinate.
  // Note: use rndR and rndTheta to uniformly sample a unit disk.
  double nwidth = 2 * tan(radians(hFov/2)), nheight = 2 * tan(radians(vFov/2));
  double nx = x * nwidth, ny = y * nheight;
  nx -= tan(radians(hFov/2)); ny -= tan(radians(vFov/2));

  Vector3D pLens = Vector3D(lensRadius * sqrt(rndR) * cos(rndTheta),
                            lensRadius * sqrt(rndR) * sin(rndTheta),
                            0);
  Vector3D pFocus = Vector3D(focalDistance * nx, focalDistance * ny, -focalDistance);
  Ray res = Ray(pos + c2w * pLens,
                (c2w * (pFocus-pLens)).unit());
  res.min_t = nClip;
  res.max_t = fClip;
  return res;
}


} // namespace CGL
