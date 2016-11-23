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
#ifndef MATRIX_UTIL_H
#define MATRIX_UTIL_H

#include <SETTINGS.h>
#include <linearalgebra/COO_MATRIX.h>
#include <Eigen/SVD>

class MATRIX_UTIL {
 public:
  static MATRIX3 cross(const VEC3F& vec) {
    MATRIX3 final;
    final(0, 0) = 0.0;
    final(1, 0) = vec[2];
    final(2, 0) = -vec[1];

    final(0, 1) = -vec[2];
    final(1, 1) = 0.0;
    final(2, 1) = vec[0];

    final(0, 2) = vec[1];
    final(1, 2) = -vec[0];
    final(2, 2) = 0.0;

    return final;
  }
};

inline std::istream& operator>>(std::istream& in, VEC3F& v) {
  return in >> v[0] >> v[1] >> v[2];
}

inline std::istream& operator>>(std::istream& in, VEC2F& v) {
  return in >> v[0] >> v[1];
}

inline VEC3F cross(const VEC3F& a, const VEC3F& b) { return a.cross(b); }
inline void unitize(VEC3F a) { a.normalize(); }
inline Real norm(const VEC3F& a) { return a.norm(); }
inline Real norm2(const VEC3F& a) { return a.dot(a); }

inline VEC3F operator^(const VEC3F& a, const VEC3F& b) { return a.cross(b); }
inline MATRIX3 operator^(const MATRIX3& n, const MATRIX3& m) {
  return n.transpose() * m;
}
inline Real trace(const MATRIX3& mat) { return mat.trace(); }

#endif