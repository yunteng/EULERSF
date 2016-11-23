/*
This file is part of ESFC (Eulerian Solid-FLuid Coupling).

ESFC is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

ESFC is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with ESFC.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef TRIANGLE_H
#define TRIANGLE_H

#include <SETTINGS.h>
#include <iostream>

#if _WIN32
#include <gl/glut.h>
#elif USING_OSX
#include <GLUT/glut.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#endif

using namespace std;

class TRIANGLE {
 public:
  TRIANGLE(VEC3F& v0, VEC3F& v1, VEC3F& v2, VEC3F color = VEC3F(1, 1, 1));
  TRIANGLE(const TRIANGLE& triangle);

  bool operator==(const TRIANGLE& RHS) const;

  void draw();
  void draw(const vector<VEC3F>& colors);

  VEC3F center();
  Real area() {
    _area =
        0.5 *
        (*vertices[1] - *vertices[0]).cross(*vertices[2] - *vertices[0]).norm();
    return _area;
  };
  const VEC3F& normal() {
    _normal = (*vertices[1] - *vertices[0]).cross(*vertices[2] - *vertices[0]);
    _normal.normalize();
    return _normal;
  };
  VEC3F restNormal() {
    VEC3F n = (*restPose[1] - *restPose[0]).cross(*restPose[2] - *restPose[0]);
    n.normalize();
    return n;
  }
  bool baryCenter(const VEC3F& point, VEC3F& lambda, bool rest = false);
  VEC3F interpPosition(VEC3F& lambda);

  void boundingBox(VEC3F& mins, VEC3F& maxs);
  // intersect with line segment
  bool intersects(const VEC3F& start, const VEC3F& end);
  VEC3F projection(const VEC3F& point);

  Real maxEdgeLength() {
    Real maxl = (*vertices[0] - *vertices[1]).norm();
    maxl = max(maxl, (*vertices[0] - *vertices[2]).norm());
    maxl = max(maxl, (*vertices[2] - *vertices[1]).norm());
    return maxl;
  }
  VEC3F* vertex(int i) { return vertices[i]; };

 public:
  VEC3F* vertices[3];
  VEC3F* restPose[3];

 private:
  VEC3F _color;
  VEC3F _normal;
  Real _area;
};
#endif