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
#ifndef TRIANGLE_BVH_H
#define TRIANGLE_BVH_H

#include <SETTINGS.h>
#include <deformCD/aabb.h>
#include <geometry/OBJ.h>

#if _WIN32
#include <gl/glut.h>
#elif USING_OSX
#include <GLUT/glut.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#endif

using namespace ::std;

class TRIANGLE_BVH_NODE {
 public:
  TRIANGLE_BVH_NODE();
  TRIANGLE_BVH_NODE(unsigned int);
  TRIANGLE_BVH_NODE(unsigned int* lst, unsigned int lst_num, aabb* triBoxes,
                    VEC3F* triCenters);

  ~TRIANGLE_BVH_NODE();

  aabb& refit(OBJ* mesh);

  bool collide(TRIANGLE_BVH_NODE* other);

  friend class TRIANGLE_BVH;

 private:
  aabb _box;
  unsigned int _id;
  TRIANGLE_BVH_NODE* _left;
  TRIANGLE_BVH_NODE* _right;
};
class TRIANGLE_BVH {
 public:
  TRIANGLE_BVH(OBJ* mesh);
  ~TRIANGLE_BVH();

  bool collide(TRIANGLE_BVH* other);
  void refit();

  bool inside(const aabb& wrapBox);

 private:
  void init(vector<TRIANGLE>& triangles);

 private:
  TRIANGLE_BVH_NODE* _root;
  OBJ* _mesh;
};

#endif
