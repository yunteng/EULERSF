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
#include <assert.h>
#include <limits.h>
#include <stdlib.h>

#include <deformCD/aap.h>
#include <geometry/TRIANGLE_BVH.h>

extern float middle_xyz(char xyz, const VEC3F& p1, const VEC3F& p2,
                        const VEC3F& p3);

TRIANGLE_BVH::TRIANGLE_BVH(OBJ* mesh) : _mesh(mesh) {
  vector<TRIANGLE>& triangles = _mesh->triangles;

  init(triangles);
}

void TRIANGLE_BVH::init(vector<TRIANGLE>& triangles) {
  aabb total;

  for (unsigned int x = 0; x < triangles.size(); x++) {
    for (int y = 0; y < 3; y++) {
      total += *(triangles[x].vertices[y]);
    }
  }

  unsigned int totalTris = triangles.size();

  aabb* triBoxes = new aabb[totalTris];
  VEC3F* triCenters = new VEC3F[totalTris];

  aap pln(total);
  unsigned int* idx_buffer = new unsigned int[totalTris];
  unsigned int left_idx = 0, right_idx = totalTris;

  for (unsigned int triID = 0; triID < totalTris; triID++) {
    triCenters[triID] = triangles[triID].center();
    if (pln.inside(triCenters[triID]))
      idx_buffer[left_idx++] = triID;
    else
      idx_buffer[--right_idx] = triID;

    for (int y = 0; y < 3; y++) {
      triBoxes[triID] += *(triangles[triID].vertices[y]);
    }
  }

  _root = new TRIANGLE_BVH_NODE();
  _root->_box = total;
  if (left_idx == 0 || left_idx == totalTris) {
    int hf = totalTris / 2;
    if (hf > 0) {
      _root->_left =
          new TRIANGLE_BVH_NODE(idx_buffer, hf, triBoxes, triCenters);
      _root->_right = new TRIANGLE_BVH_NODE(idx_buffer + hf, totalTris - hf,
                                            triBoxes, triCenters);
    } else {
      _root->_left =
          new TRIANGLE_BVH_NODE(idx_buffer, totalTris, triBoxes, triCenters);
      _root->_right = NULL;
    }
  } else {
    _root->_left =
        new TRIANGLE_BVH_NODE(idx_buffer, left_idx, triBoxes, triCenters);
    _root->_right = new TRIANGLE_BVH_NODE(
        idx_buffer + left_idx, totalTris - left_idx, triBoxes, triCenters);
  }

  delete[] triBoxes;
  delete[] triCenters;
  delete[] idx_buffer;
}
TRIANGLE_BVH::~TRIANGLE_BVH() { delete _root; }

void TRIANGLE_BVH::refit() { _root->refit(_mesh); }

bool TRIANGLE_BVH::inside(const aabb& wrapBox) {
  return wrapBox.inside(_root->_box._min) && wrapBox.inside(_root->_box._max);
}

bool TRIANGLE_BVH::collide(TRIANGLE_BVH* other) {
  return _root->collide(other->_root);
}
//#################################################################

TRIANGLE_BVH_NODE::TRIANGLE_BVH_NODE() {
  _id = UINT_MAX;
  _left = _right = NULL;
}

TRIANGLE_BVH_NODE::TRIANGLE_BVH_NODE(unsigned int id) {
  _left = _right = NULL;
  _id = id;
}

TRIANGLE_BVH_NODE::TRIANGLE_BVH_NODE(unsigned int* lst, unsigned int lst_num,
                                     aabb* triBoxes, VEC3F* triCenters) {
  assert(lst_num > 0);
  _left = _right = NULL;
  _id = UINT_MAX;

  if (lst_num == 1) {
    _id = lst[0];
    _box = triBoxes[_id];
  } else {  // try to split them
    for (unsigned int t = 0; t < lst_num; t++) {
      int i = lst[t];
      _box += triBoxes[i];
    }

    if (lst_num == 2) {  // must split it!
      _left = new TRIANGLE_BVH_NODE(lst[0]);
      _right = new TRIANGLE_BVH_NODE(lst[1]);
    } else {
      aap pln(_box);
      unsigned int left_idx = 0, right_idx = lst_num - 1;

      for (unsigned int t = 0; t < lst_num; t++) {
        int i = lst[left_idx];
        if (pln.inside(triCenters[i]))
          left_idx++;
        else {  // swap it
          unsigned int tmp = i;
          lst[left_idx] = lst[right_idx];
          lst[right_idx--] = tmp;
        }
      }

      int hal = lst_num / 2;
      if (left_idx == 0 || left_idx == lst_num) {
        _left = new TRIANGLE_BVH_NODE(lst, hal, triBoxes, triCenters);
        _right = new TRIANGLE_BVH_NODE(lst + hal, lst_num - hal, triBoxes,
                                       triCenters);
      } else {
        _left = new TRIANGLE_BVH_NODE(lst, left_idx, triBoxes, triCenters);
        _right = new TRIANGLE_BVH_NODE(lst + left_idx, lst_num - left_idx,
                                       triBoxes, triCenters);
      }
    }
  }
}

TRIANGLE_BVH_NODE::~TRIANGLE_BVH_NODE() {
  if (_left) delete _left;
  if (_right) delete _right;
}

aabb& TRIANGLE_BVH_NODE::refit(OBJ* mesh) {
  if (_left == NULL) {
    _box.empty();
    // for(int i = 0; i < _ids.size(); i++) {
    TRIANGLE& tri = mesh->triangles[_id];
    for (int y = 0; y < 3; y++) {
      _box += *(tri.vertices[y]);
    }
    // }
  } else {
    _box = _left->refit(mesh);
    _box += _right->refit(mesh);
  }
  return _box;
}

bool TRIANGLE_BVH_NODE::collide(TRIANGLE_BVH_NODE* other) {
  if (!_box.overlaps(other->_box)) return false;
  if (_left == NULL && other->_left == NULL) {
    return _box.overlaps(other->_box);
  }
  if (_left == NULL) {
    return collide(other->_left) || collide(other->_right);
  } else if (other->_left != NULL) {
    return _left->collide(other->_left) || _left->collide(other->_right) ||
           _right->collide(other->_left) || _right->collide(other->_right);
  } else {
    return _left->collide(other) || _right->collide(other);
  }
}
