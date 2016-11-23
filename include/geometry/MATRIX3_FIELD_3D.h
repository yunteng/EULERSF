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
#ifndef MATRIX3_FIELD_3D_H
#define MATRIX3_FIELD_3D_H

#include <SETTINGS.h>
#include <assert.h>
#include <util/MATH_UTILS.h>
#include <zlib.h>
#include <cstdlib>
#include <iostream>

using namespace std;

class MATRIX3_FIELD_3D {
 public:
  // construction and destruction
  MATRIX3_FIELD_3D()
      : _xRes(0),
        _yRes(0),
        _zRes(0),
        _origin(0, 0, 0),
        _dh(0),
        _slabSize(0),
        _totalCells(0){};
  MATRIX3_FIELD_3D(int xRes, int yRes, int zRes, Real dh, const VEC3F& origin)
      : _xRes(xRes), _yRes(yRes), _zRes(zRes), _origin(origin), _dh(dh) {
    _slabSize = _xRes * _yRes;
    _totalCells = _xRes * _yRes * _zRes;
    _data.resize(_totalCells, MATRIX3::Identity());
  };

  virtual ~MATRIX3_FIELD_3D(){};

  // accessors
  int xRes() const { return _xRes; };
  int yRes() const { return _yRes; };
  int zRes() const { return _zRes; };
  Real dh() const { return _dh; };
  const VEC3F& origin() const { return _origin; };
  int totalCells() const { return _totalCells; };

  // (x,y,z) accessor
  inline MATRIX3& operator()(int x, int y, int z) {
    assert(x >= 0);
    assert(x < _xRes);
    assert(y >= 0);
    assert(y < _yRes);
    assert(z >= 0);
    assert(z < _zRes);
    return _data[x + y * _xRes + z * _slabSize];
  };

  inline MATRIX3 operator()(int x, int y, int z) const {
    assert(x >= 0);
    assert(x < _xRes);
    assert(y >= 0);
    assert(y < _yRes);
    assert(z >= 0);
    assert(z < _zRes);
    return _data[x + y * _xRes + z * _slabSize];
  };

  // direct indexing
  inline MATRIX3& operator[](int index) {
    assert(index >= 0);
    assert(index < _totalCells);
    return _data[index];
  }

  inline MATRIX3 operator[](int index) const {
    assert(index >= 0);
    assert(index < _totalCells);
    return _data[index];
  }
  // set one field to another
  MATRIX3_FIELD_3D& operator=(const MATRIX3_FIELD_3D& field) {
    _xRes = field.xRes();
    _yRes = field.yRes();
    _zRes = field.zRes();
    _origin = field.origin();
    _dh = field.dh();
    _slabSize = _xRes * _yRes;
    _totalCells = _xRes * _yRes * _zRes;
    _data = field._data;
    return *this;
  }

  // swap the data between this and another field
  inline void swap(MATRIX3_FIELD_3D& other) {
    assert(other.xRes() == _xRes && other.yRes() == _yRes &&
           other.zRes() == _zRes);
    assert(abs(other.dh() - _dh) < 1e-6);
    assert((other.origin() - _origin).norm() < 1e-6);
    _data.swap(other._data);
  }

  // set the size of the array
  void resize(int xRes, int yRes, int zRes, Real dh, const VEC3F& origin) {
    _xRes = xRes;
    _yRes = yRes;
    _zRes = zRes;
    _origin = origin;
    _dh = dh;
    _totalCells = _xRes * _yRes * _zRes;
    _slabSize = _xRes * _yRes;
    _data.resize(_totalCells, MATRIX3::Identity());
  }

  // set to zero
  inline void clear() { _data.resize(_totalCells, MATRIX3::Identity()); }

  MATRIX3 interpolate_value(const VEC3F& point, int& index) {
    int i, j, k;
    get_barycentric(point[0] - _origin[0], i, 0, _xRes - 1);
    get_barycentric(point[1] - _origin[1], j, 0, _yRes - 1);
    get_barycentric(point[2] - _origin[2], k, 0, _zRes - 1);
    index = k * _slabSize + j * _xRes + i;

    return _data[index];
  }

  void writeGz(const string& filename) const {
    gzFile file = gzopen(filename.c_str(), "wb1");
    if (file == NULL) {
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : "
           << endl;
      cout << " FIELD_3D write failed! " << endl;
      cout << " Could not open file " << filename.c_str() << endl;
      exit(0);
    }
    writeGz(file);
    gzclose(file);
  }

  void writeGz(gzFile& file) const {
    // write dimensions
    gzwrite(file, (void*)&_xRes, sizeof(int));
    gzwrite(file, (void*)&_yRes, sizeof(int));
    gzwrite(file, (void*)&_zRes, sizeof(int));
    gzwrite(file, (void*)&_dh, sizeof(Real));
    gzwrite(file, (void*)&_origin[0], sizeof(Real));
    gzwrite(file, (void*)&_origin[1], sizeof(Real));
    gzwrite(file, (void*)&_origin[2], sizeof(Real));

    for (int x = 0; x < _totalCells; x++) {
      double mat[9];
      for (int j = 0; j < 3; j++)
        for (int i = 0; i < 3; i++) mat[j * 3 + i] = _data[x](i, j);
      gzwrite(file, (void*)mat, 9 * sizeof(double));
    }
  }

  void readGz(const string& filename) {
    gzFile file = gzopen(filename.c_str(), "rb");
    if (file == NULL) {
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : "
           << endl;
      cout << " FIELD_3D read failed! " << endl;
      cout << " Could not open file " << filename << endl;
      exit(0);
    }
    readGz(file);
    gzclose(file);
  }

  void readGz(gzFile& file) {
    // read dimensions
    gzread(file, (void*)&_xRes, sizeof(int));
    gzread(file, (void*)&_yRes, sizeof(int));
    gzread(file, (void*)&_zRes, sizeof(int));
    gzread(file, (void*)&_dh, sizeof(Real));
    gzread(file, (void*)&_origin[0], sizeof(Real));
    gzread(file, (void*)&_origin[1], sizeof(Real));
    gzread(file, (void*)&_origin[2], sizeof(Real));

    _slabSize = _xRes * _yRes;

    _totalCells = _xRes * _yRes * _zRes;

    _data.resize(_totalCells);

    for (int x = 0; x < _totalCells; x++) {
      double mat[9];
      gzread(file, (void*)mat, 9 * sizeof(double));
      for (int j = 0; j < 3; j++)
        for (int i = 0; i < 3; i++) _data[x](i, j) = mat[j * 3 + i];
    }
  }

 private:
  inline void get_barycentric(Real x, int& i, int i_low, int i_high) const {
    x /= _dh;
    Real s = std::floor(x);
    i = int(s);
    if (i <= i_low) {
      i = i_low;
    } else if (i >= i_high) {
      i = i_high;
    }
  }

 private:
  int _xRes;
  int _yRes;
  int _zRes;
  VEC3F _origin;
  Real _dh;
  int _slabSize;
  int _totalCells;

  vector<MATRIX3> _data;
};

#endif
