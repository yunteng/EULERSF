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
#ifndef FIELD_3D_H
#define FIELD_3D_H

#include <SETTINGS.h>
#include <assert.h>
#include <util/MATH_UTILS.h>
#include <zlib.h>
#include <cstdlib>
#include <fstream>
#include <iostream>

using namespace std;

template <class T>
class FIELD_3D {
 public:
  // construction and destruction
  FIELD_3D()
      : _xRes(0),
        _yRes(0),
        _zRes(0),
        _origin(0, 0, 0),
        _dh(0),
        _slabSize(0),
        _totalCells(0){};
  FIELD_3D(int xRes, int yRes, int zRes, Real dh, const VEC3F& origin)
      : _xRes(xRes), _yRes(yRes), _zRes(zRes), _origin(origin), _dh(dh) {
    _slabSize = _xRes * _yRes;
    _totalCells = _xRes * _yRes * _zRes;
    _data = Eigen::Array<T, Eigen::Dynamic, 1>::Zero(xRes * yRes * zRes);
    cout << "data size " << _data.size() << " " << xRes * yRes * zRes << endl;
  };

  virtual ~FIELD_3D(){};

  // accessors
  int xRes() const { return _xRes; };
  int yRes() const { return _yRes; };
  int zRes() const { return _zRes; };
  Real dh() const { return _dh; };
  const VEC3F& origin() const { return _origin; };
  int totalCells() const { return _totalCells; };

  Eigen::Array<T, Eigen::Dynamic, 1>& data() { return _data; };

  const Eigen::Array<T, Eigen::Dynamic, 1>& data() const { return _data; };

  // (x,y,z) accessor
  inline T& operator()(int x, int y, int z) {
    assert(x >= 0);
    assert(x < _xRes);
    assert(y >= 0);
    assert(y < _yRes);
    assert(z >= 0);
    assert(z < _zRes);
    return _data[x + y * _xRes + z * _slabSize];
  };

  inline T operator()(int x, int y, int z) const {
    assert(x >= 0);
    assert(x < _xRes);
    assert(y >= 0);
    assert(y < _yRes);
    assert(z >= 0);
    assert(z < _zRes);
    return _data[x + y * _xRes + z * _slabSize];
  };

  // direct indexing
  inline T& operator[](int index) {
    assert(index >= 0);
    assert(index < _totalCells);
    return _data[index];
  }

  inline const T operator[](int index) const {
    assert(index >= 0);
    assert(index < _totalCells);
    return _data[index];
  }
  // set one field to another
  FIELD_3D<T>& operator=(const FIELD_3D<T>& field) {
    _xRes = field.xRes();
    _yRes = field.yRes();
    _zRes = field.zRes();
    _origin = field.origin();
    _dh = field.dh();
    _slabSize = _xRes * _yRes;
    _totalCells = _xRes * _yRes * _zRes;
    _data = field.data();
    return *this;
  }
  FIELD_3D<T>& operator=(T alpha) {
    _data = Eigen::Array<T, Eigen::Dynamic, 1>::Constant(_totalCells, alpha);
    return *this;
  }
  // scale the field by a constant
  FIELD_3D<T> operator*=(const Real& scalar) {
    _data *= scalar;
    return *this;
  }

  // add another field to this field
  FIELD_3D& operator+=(const FIELD_3D<T>& field) {
    assert(field.xRes() == _xRes && field.yRes() == _yRes &&
           field.zRes() == _zRes);
    assert(abs(field.dh() - _dh) < 1e-6);
    assert((field.origin() - _origin).norm() < 1e-6);
    _data += field.data();
    return *this;
  }

  // swap the data between this and another field
  inline void swap(FIELD_3D<T>& other) {
    assert(other.xRes() == _xRes && other.yRes() == _yRes &&
           other.zRes() == _zRes);
    assert(abs(other.dh() - _dh) < 1e-6);
    assert((other.origin() - _origin).norm() < 1e-6);
    _data.swap(other.data());
  }

  // subtract another field from this field
  FIELD_3D& operator-=(const FIELD_3D<T>& field) {
    assert(field.xRes() == _xRes && field.yRes() == _yRes &&
           field.zRes() == _zRes);
    assert(abs(field.dh() - _dh) < 1e-6);
    assert((field.origin() - _origin).norm() < 1e-6);
    _data -= field.data();
    return *this;
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
    _data = Eigen::Array<T, Eigen::Dynamic, 1>::Zero(xRes * yRes * zRes);
  }

  void resizeLike(const FIELD_3D<T>& field) {
    _xRes = field.xRes();
    _yRes = field.yRes();
    _zRes = field.zRes();
    _origin = field.origin();
    _dh = field.dh();
    _totalCells = _xRes * _yRes * _zRes;
    _slabSize = _xRes * _yRes;
    _data = Eigen::Array<T, Eigen::Dynamic, 1>::Zero(_totalCells);
  }

  // set to zero
  inline void clear() { _data.setZero(); }

  // squared sum of entries
  Real squaredNorm() const { return _data.matrix().squaredNorm(); };

  Real norm() const { return _data.matrix().norm(); };

  // dot product with another field
  Real dot(const FIELD_3D<T>& field) const {
    assert(field.xRes() == _xRes && field.yRes() == _yRes &&
           field.zRes() == _zRes);
    assert(abs(field.dh() - _dh) < 1e-6);
    assert((field.origin() - _origin).norm() < 1e-6);
    return (_data * field.data()).sum();
  }

  // axpy operation with another field
  // this = this + scalar * field
  void axpy(const Real& scalar, const FIELD_3D<T>& field) {
    assert(field.xRes() == _xRes && field.yRes() == _yRes &&
           field.zRes() == _zRes);
    assert(abs(field.dh() - _dh) < 1e-6);
    assert((field.origin() - _origin).norm() < 1e-6);
    _data += field.data() * scalar;
  }

  T interpolate_value(const VEC3F& point) const {
    int i, j, k;
    Real fx, fy, fz;

    get_barycentric(point[0] - _origin[0], i, fx, 0, _xRes - 2);
    get_barycentric(point[1] - _origin[1], j, fy, 0, _yRes - 2);
    get_barycentric(point[2] - _origin[2], k, fz, 0, _zRes - 2);

    T v000 = (*this)(i, j, k);
    T v001 = (*this)(i, j, k + 1);
    T v010 = (*this)(i, j + 1, k);
    T v011 = (*this)(i, j + 1, k + 1);
    T v100 = (*this)(i + 1, j, k);
    T v101 = (*this)(i + 1, j, k + 1);
    T v110 = (*this)(i + 1, j + 1, k);
    T v111 = (*this)(i + 1, j + 1, k + 1);

    return trilerp(v000, v100, v010, v110, v001, v101, v011, v111, fx, fy, fz);
  }

  T interpolate_gradient(VEC3F& gradient, const VEC3F& point) const {
    int i, j, k;
    Real fx, fy, fz;

    get_barycentric(point[0] - _origin[0], i, fx, 0, _xRes - 2);
    get_barycentric(point[1] - _origin[1], j, fy, 0, _yRes - 2);
    get_barycentric(point[2] - _origin[2], k, fz, 0, _zRes - 2);

    T v000 = (*this)(i, j, k);
    T v001 = (*this)(i, j, k + 1);
    T v010 = (*this)(i, j + 1, k);
    T v011 = (*this)(i, j + 1, k + 1);
    T v100 = (*this)(i + 1, j, k);
    T v101 = (*this)(i + 1, j, k + 1);
    T v110 = (*this)(i + 1, j + 1, k);
    T v111 = (*this)(i + 1, j + 1, k + 1);

    T ddx00 = (v100 - v000);
    T ddx10 = (v110 - v010);
    T ddx01 = (v101 - v001);
    T ddx11 = (v111 - v011);
    T dv_dx = bilerp(ddx00, ddx10, ddx01, ddx11, fy, fz);

    T ddy00 = (v010 - v000);
    T ddy10 = (v110 - v100);
    T ddy01 = (v011 - v001);
    T ddy11 = (v111 - v101);
    T dv_dy = bilerp(ddy00, ddy10, ddy01, ddy11, fx, fz);

    T ddz00 = (v001 - v000);
    T ddz10 = (v101 - v100);
    T ddz01 = (v011 - v010);
    T ddz11 = (v111 - v110);
    T dv_dz = bilerp(ddz00, ddz10, ddz01, ddz11, fx, fy);

    gradient[0] = dv_dx;
    gradient[1] = dv_dy;
    gradient[2] = dv_dz;

    Real mag = gradient.norm();
    if (mag > 0) gradient /= mag;

    // return value for good measure.
    return trilerp(v000, v100, v010, v110, v001, v101, v011, v111, fx, fy, fz);
  }
  void writePbrt(const string& filename) const {
    ofstream out(filename.c_str());
    if (out.fail()) {
      cerr << "Failed to open " << filename << " to save PBRT" << endl;
    }
    out << "MakeNamedMedium \"smoke\" ";
    out << "\"string type\" \"heterogeneous\" ";
    out << "\"integer nx\" " << _xRes << " ";
    out << "\"integer ny\" " << _yRes << " ";
    out << "\"integer nz\" " << _zRes << endl;
    out << "\"point p0\" [" << _origin[0] << " " << _origin[1] << " "
        << _origin[2] << "]\n";
    out << "\"point p1\" [" << _origin[0] + (_xRes - 1) * _dh << " "
        << _origin[1] + (_yRes - 1) * _dh << " "
        << _origin[2] + (_zRes - 1) * _dh << "]";
    out << "\"float density\" [\n";
    for (int x = 0; x < _totalCells; x++) {
      out << _data[x] << " ";
    }
    out << "]" << endl;
    out.close();
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

    T test = 0;
    // always write out as a double
    if (checkType(test) == FLOAT) {
      double* dataDouble = new double[_totalCells];
      for (int x = 0; x < _totalCells; x++) dataDouble[x] = _data[x];

      gzwrite(file, (void*)dataDouble, _totalCells * sizeof(double));
      delete[] dataDouble;
    } else {
      const T* addr = &(_data[0]);
      gzwrite(file, (void*)addr, _totalCells * sizeof(T));
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

    _data.conservativeResize(_totalCells);

    T test = 0;
    // always read in as a double
    if (checkType(test) == FLOAT) {
      double* dataDouble = new double[_totalCells];
      gzread(file, (void*)dataDouble, _totalCells * sizeof(double));

      for (int x = 0; x < _totalCells; x++) _data[x] = dataDouble[x];

      delete[] dataDouble;
    } else {
      T* addr = &(_data[0]);
      gzread(file, (void*)addr, _totalCells * sizeof(T));
    }
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
  VEC3F _origin;
  Real _dh;
  int _slabSize;
  int _totalCells;

  Eigen::Array<T, Eigen::Dynamic, 1> _data;

  enum DataType { CHAR, INT, UNSIGNED_INT, FLOAT, DOUBLE };

  DataType checkType(char f) const { return CHAR; };
  DataType checkType(int f) const { return INT; };
  DataType checkType(unsigned int f) const { return UNSIGNED_INT; };
  DataType checkType(float f) const { return FLOAT; };
  DataType checkType(double f) const { return DOUBLE; };
};

typedef FIELD_3D<Real> FIELD_3Df;
typedef FIELD_3D<unsigned int> FIELD_3Dui;
typedef FIELD_3D<char> FIELD_3Dc;

#endif
