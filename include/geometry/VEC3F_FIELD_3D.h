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
#ifndef VEC3F_FIELD_3D_H
#define VEC3F_FIELD_3D_H

#include <SETTINGS.h>
#include <util/MATH_UTILS.h>
#include <zlib.h>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <vector>

using namespace std;

class VEC3F_FIELD_3D {
 public:
  // construction and destruction
  VEC3F_FIELD_3D()
      : _xRes(0),
        _yRes(0),
        _zRes(0),
        _origin(0, 0, 0),
        _dh(0),
        _slabSize(0),
        _totalCells(0){};

  VEC3F_FIELD_3D(int xRes, int yRes, int zRes, Real dh, const VEC3F& origin)
      : _xRes(xRes), _yRes(yRes), _zRes(zRes), _origin(origin), _dh(dh) {
    _slabSize = _xRes * _yRes;
    _totalCells = _xRes * _yRes * _zRes;
    _data.resize(xRes * yRes * zRes, VEC3F(0, 0, 0));
  };

  virtual ~VEC3F_FIELD_3D(){};

  // accessors
  int xRes() const { return _xRes; };
  int yRes() const { return _yRes; };
  int zRes() const { return _zRes; };
  Real dh() const { return _dh; };
  int totalCells() const { return _totalCells; };
  const VEC3F& origin() const { return _origin; };

  vector<VEC3F>& data() { return _data; };

  const vector<VEC3F>& data() const { return _data; };

  // (x,y,z) accessor
  inline VEC3F& operator()(int x, int y, int z) {
    assert(x >= 0);
    assert(x < _xRes);
    assert(y >= 0);
    assert(y < _yRes);
    assert(z >= 0);
    assert(z < _zRes);
    return _data[x + y * _xRes + z * _slabSize];
  };

  inline const VEC3F operator()(int x, int y, int z) const {
    assert(x >= 0);
    assert(x < _xRes);
    assert(y >= 0);
    assert(y < _yRes);
    assert(z >= 0);
    assert(z < _zRes);
    return _data[x + y * _xRes + z * _slabSize];
  };

  // direct indexing
  inline VEC3F& operator[](int index) {
    assert(index >= 0);
    assert(index < _totalCells);
    return _data[index];
  }

  inline const VEC3F& operator[](int index) const {
    assert(index >= 0);
    assert(index < _totalCells);
    return _data[index];
  }

  // scale the field by a constant scalar
  VEC3F_FIELD_3D& operator*=(const Real& scalar) {
    for (unsigned int x = 0; x < _data.size(); x++) _data[x] *= scalar;
    return *this;
  }
  // scale the field by a constant VEC3F
  VEC3F_FIELD_3D& operator*=(const VEC3F& scale) {
    for (unsigned int x = 0; x < _data.size(); x++) {
      _data[x][0] *= scale[0];
      _data[x][1] *= scale[1];
      _data[x][1] *= scale[2];
    }
    return *this;
  }

  // add another field to this field
  VEC3F_FIELD_3D& operator+=(const VEC3F_FIELD_3D& field) {
    assert(field.xRes() == _xRes && field.yRes() == _yRes &&
           field.zRes() == _zRes);
    assert(abs(field.dh() - _dh) < 1e-6);
    assert((field.origin() - _origin).norm() < 1e-6);

    const vector<VEC3F>& other = field.data();
    for (unsigned int x = 0; x < _data.size(); x++) _data[x] += other[x];
    return *this;
  }

  // subtract another field from this field
  VEC3F_FIELD_3D& operator-=(const VEC3F_FIELD_3D& field) {
    assert(field.xRes() == _xRes && field.yRes() == _yRes &&
           field.zRes() == _zRes);
    assert(abs(field.dh() - _dh) < 1e-6);
    assert((field.origin() - _origin).norm() < 1e-6);
    const vector<VEC3F>& other = field.data();
    for (unsigned int x = 0; x < _data.size(); x++) _data[x] -= other[x];
    return *this;
  }

  // squared sum of entries
  Real squaredNorm() const {
    Real sum = 0;
    for (unsigned int x = 0; x < _data.size(); x++)
      sum += _data[x][0] * _data[x][0] + _data[x][1] * _data[x][1] +
             _data[x][2] * _data[x][2];
    return sum / _totalCells;
  };

  // swap the data between this and another field
  inline void swap(VEC3F_FIELD_3D& other) {
    assert(other.xRes() == _xRes && other.yRes() == _yRes &&
           other.zRes() == _zRes);
    assert(abs(other.dh() - _dh) < 1e-6);
    assert((other.origin() - _origin).norm() < 1e-6);
    _data.swap(other.data());
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
    _data.resize(_totalCells, VEC3F(0, 0, 0));
  }

  void resizeLike(const VEC3F_FIELD_3D& field) {
    _xRes = field.xRes();
    _yRes = field.yRes();
    _zRes = field.zRes();
    _origin = field.origin();
    _dh = field.dh();
    _totalCells = _xRes * _yRes * _zRes;
    _slabSize = _xRes * _yRes;
    _data.resize(_totalCells, VEC3F(0, 0, 0));
  }

  // set to zero
  inline void clear() { _data.resize(_totalCells, VEC3F(0, 0, 0)); };

  VEC3F interpolate_value(const VEC3F& point) {
    int i, j, k;
    Real fx, fy, fz;

    get_barycentric(point[0] - _origin[0], i, fx, 0, _xRes - 2);
    get_barycentric(point[1] - _origin[1], j, fy, 0, _yRes - 2);
    get_barycentric(point[2] - _origin[2], k, fz, 0, _zRes - 2);

    const VEC3F& v000 = (*this)(i, j, k);
    const VEC3F& v001 = (*this)(i, j, k + 1);
    const VEC3F& v010 = (*this)(i, j + 1, k);
    const VEC3F& v011 = (*this)(i, j + 1, k + 1);
    const VEC3F& v100 = (*this)(i + 1, j, k);
    const VEC3F& v101 = (*this)(i + 1, j, k + 1);
    const VEC3F& v110 = (*this)(i + 1, j + 1, k);
    const VEC3F& v111 = (*this)(i + 1, j + 1, k + 1);

    return trilerp(v000, v100, v010, v110, v001, v101, v011, v111, fx, fy, fz);
  }

  void writeGz(const string& filename) const {
    gzFile file = gzopen(filename.c_str(), "wb1");
    if (file == NULL) {
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : "
           << endl;
      cout << " VEC3F_FIELD_3D write failed! " << endl;
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

    double* dataDouble = new double[_totalCells * 3];
    for (int i = 0; i < _totalCells; i++) {
      dataDouble[i * 3] = _data[i][0];
      dataDouble[i * 3 + 1] = _data[i][1];
      dataDouble[i * 3 + 2] = _data[i][2];
    }
    gzwrite(file, (void*)dataDouble, _totalCells * 3 * sizeof(double));
    delete[] dataDouble;
  }

  void readGz(const string& filename) {
    gzFile file = gzopen(filename.c_str(), "rb");
    if (file == NULL) {
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : "
           << endl;
      cout << " VEC3F_FIELD_3D read failed! " << endl;
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

    _totalCells = _xRes * _yRes * _zRes;
    _slabSize = _xRes * _yRes;

    _data.resize(_totalCells);

    double* dataDouble = new double[_totalCells * 3];
    gzread(file, (void*)dataDouble, _totalCells * 3 * sizeof(double));

    for (int i = 0; i < _totalCells; i++) {
      _data[i][0] = dataDouble[i * 3];
      _data[i][1] = dataDouble[i * 3 + 1];
      _data[i][2] = dataDouble[i * 3 + 2];
    }
    delete[] dataDouble;
  }

 private:
  inline void get_barycentric(Real x, int& i, Real& f, int i_low,
                              int i_high) const {
    x /= _dh;
    Real s = std::floor(x);
    i = int(s);
    if (i < i_low) {
      i = i_low;
      f = 0.5;
    } else if (i > i_high) {
      i = i_high;
      f = 0.5;
    } else {
      f = x - s;
    }
  }

 private:
  int _xRes;
  int _yRes;
  int _zRes;
  // lower left corner
  VEC3F _origin;
  Real _dh;
  int _slabSize;
  int _totalCells;

  vector<VEC3F> _data;
};

#endif
