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
#ifndef COO_MATRIX_H
#define COO_MATRIX_H

#if USING_LINUX
#include <parallel/algorithm>
#else
#include <algorithm>  // std::sort
#endif

#include <SETTINGS.h>
#include <iostream>
#include <vector>  // std::vector

#if USING_OPENMP
#include <omp.h>
#endif

using namespace ::std;

class COO_MATRIX {
 public:
  COO_MATRIX() : _rows(0), _cols(0){};
  COO_MATRIX(int rows, int cols) : _rows(rows), _cols(cols){};

  COO_MATRIX(const MATRIX& mat) {
    _rows = mat.rows();
    _cols = mat.cols();

    _diag.resize(min(_rows, _cols));
    _diag.setZero();

    for (int x = 0; x < mat.rows(); x++)
      for (int y = 0; y < mat.cols(); y++) {
        if (abs(mat(x, y)) > 1e-16) _matrix.push_back(TRIPLET(x, y, mat(x, y)));
        if (x == y) _diag[x] = mat(x, y);
      }
  }
  COO_MATRIX(const COO_MATRIX& other) {
    _rows = other.rows();
    _cols = other.cols();
    _matrix = other.matrix();
    _diag = other.diag();
  }
  void setFromSpMat(const SpMat& mat) {
    _rows = mat.rows();
    _cols = mat.cols();
    _matrix.resize(mat.nonZeros());
    int index = 0;
    for (int k = 0; k < mat.outerSize(); ++k)
      for (SpMat::InnerIterator it(mat, k); it; ++it) {
        _matrix[index] = TRIPLET(it.row(), it.col(), it.value());
        index++;
      }
    // assert(index == _matrix.size());
    order();
    aggregate();
  }

  inline int nnZ() const { return _matrix.size(); };
  inline int lowerNNz() const {
    int cnt = 0;
    for (unsigned int x = 0; x < _matrix.size(); x++)
      if (_matrix[x].row() >= _matrix[x].col()) cnt++;
    return cnt;
  };
  int rows() const { return _rows; };
  int cols() const { return _cols; };

  void resize(int rows, int cols) {
    _rows = rows;
    _cols = cols;
    _matrix.clear();
    _diag.resize(min(_rows, _cols));
    _diag.setZero();
  };
  void conservativeResize(int rows, int cols) {
    _rows = rows;
    _cols = cols;
  }

  vector<TRIPLET>& matrix() { return _matrix; };
  VECTOR& diag() { return _diag; };
  const vector<TRIPLET>& matrix() const { return _matrix; };
  const VECTOR& diag() const { return _diag; };

  void order() {
#if USING_LINUX
    __gnu_parallel::sort(_matrix.begin(), _matrix.end(), sortFunc);
#else
    std::sort(_matrix.begin(), _matrix.end(), sortFunc);
#endif
  };

  inline void clear() {
    _matrix.clear();
    _diag.setZero();
  };

  int maxRowElements() {
    int ret = 0;

    int previousRow = -1;
    int rowEleCount = 0;
    for (unsigned int x = 0; x < _matrix.size(); x++) {
      int row = _matrix[x].row();
      if (row != previousRow) {
        previousRow = row;
        if (previousRow != -1) {
          ret = ret > rowEleCount ? ret : rowEleCount;
        }
        rowEleCount = 1;
      } else {
        rowEleCount++;
      }
    }
    ret = ret > rowEleCount ? ret : rowEleCount;
    return ret;
  }
  void toSpMat(SpMat& spMat) const {
    spMat.resize(_rows, _cols);
    spMat.setFromTriplets(_matrix.begin(), _matrix.end());
  }

  void entries(vector<int>& rows, vector<int>& cols,
               vector<Real>& values) const {
    // wipe previous entries
    rows.clear();
    cols.clear();
    values.clear();

    const vector<TRIPLET>& triplets = this->matrix();

    for (unsigned int x = 0; x != triplets.size(); x++) {
      // int i = triplet[x].first.first;
      // int k = triplet[x].first.second;
      rows.push_back(triplets[x].row());
      cols.push_back(triplets[x].col());
      values.push_back(triplets[x].value());
    }
  }
  void aggregate(bool chopUpperTriangle = false) {
    if (_matrix.empty()) return;
    _diag.conservativeResize(min(_rows, _cols));
    _diag.setZero();

    int cleanedPtr = 0;
    if (!chopUpperTriangle) {
      for (unsigned int arrayPtr = 1; arrayPtr < _matrix.size(); arrayPtr++) {
        if (_matrix[arrayPtr].row() == _matrix[cleanedPtr].row() &&
            _matrix[arrayPtr].col() == _matrix[cleanedPtr].col()) {
          _matrix[cleanedPtr].value() += _matrix[arrayPtr].value();
        } else {
          cleanedPtr++;
          _matrix[cleanedPtr] = _matrix[arrayPtr];
        }
        if (_matrix[arrayPtr].row() == _matrix[arrayPtr].col())
          _diag[_matrix[arrayPtr].row()] += _matrix[arrayPtr].value();
      }
    } else {
      for (unsigned int arrayPtr = 1; arrayPtr < _matrix.size(); arrayPtr++) {
        // row < col, skip
        if (_matrix[arrayPtr].row() < _matrix[arrayPtr].col()) continue;
        if (_matrix[arrayPtr].row() == _matrix[cleanedPtr].row() &&
            _matrix[arrayPtr].col() == _matrix[cleanedPtr].col())
          _matrix[cleanedPtr].value() += _matrix[arrayPtr].value();
        else {
          cleanedPtr++;
          _matrix[cleanedPtr] = _matrix[arrayPtr];
        }
        if (_matrix[arrayPtr].row() == _matrix[arrayPtr].col())
          _diag[_matrix[arrayPtr].row()] += _matrix[arrayPtr].value();
      }
    }

    _matrix.resize(cleanedPtr + 1);
  }
  inline void add(const COO_MATRIX& A, int row, int col) {
    const vector<TRIPLET>& otherMat = A.matrix();
    for (unsigned int x = 0; x < otherMat.size(); x++) {
      const TRIPLET& ele = otherMat[x];
      _matrix.push_back(TRIPLET(row + ele.row(), col + ele.col(), ele.value()));
    }
    _rows = _rows > row + A.rows() ? _rows : row + A.rows();
    _cols = _cols > col + A.cols() ? _cols : col + A.cols();
  }
  inline void add(const MATRIX& A, int row, int col) {
    for (int y = 0; y < A.rows(); y++)
      for (int x = 0; x < A.cols(); x++)
        _matrix.push_back(TRIPLET(y + row, x + col, A(y, x)));

    _rows = _rows > row + A.rows() ? _rows : row + A.rows();
    _cols = _cols > col + A.cols() ? _cols : col + A.cols();
  }
  inline void add3x3(const MATRIX3& A, int row, int col) {
    assert(_rows >= row + 3 && _cols >= col + 3);
    for (int x = 0; x < 3; x++)
      for (int y = 0; y < 3; y++)
        _matrix.push_back(TRIPLET(x + row, y + col, A(x, y)));
  }

  inline void add(const Real& entry, int row, int col) {
    // if(row < _rows && col < _cols)
    _matrix.push_back(TRIPLET(row, col, entry));
  }

  inline void subtract(const MATRIX& A, int row, int col) {
    for (int y = 0; y < A.rows(); y++)
      for (int x = 0; x < A.cols(); x++)
        _matrix.push_back(TRIPLET(y + row, x + col, -A(y, x)));

    _rows = _rows > A.rows() ? _rows : A.rows();
    _cols = _cols > A.cols() ? _cols : A.cols();
  }

  inline void subtract(Real& entry, int row, int col) {
    _matrix.push_back(TRIPLET(row, col, -entry));
  }

  VECTOR operator*(const VECTOR& v) {
    assert(this->cols() == v.size());
    VECTOR res(this->rows());
    res.setZero();
    const vector<TRIPLET>& mat = this->matrix();

    if (mat.size() == 0) return res;

#if USING_OPENMP
#pragma omp parallel

    {
      // VECTOR localRes(this->rows());
      // localRes.setZero();
      const int id = omp_get_thread_num();
      _mulVecCopies[id].setZero();
#pragma omp for schedule(static)
      for (unsigned int x = 0; x < mat.size(); x++) {
        // localRes[mat[x].row()] += mat[x].value() * v[mat[x].col()];
        _mulVecCopies[id][mat[x].row()] += mat[x].value() * v[mat[x].col()];
        // res[mat[x].row()] += mat[x].value() * v[mat[x].col()];
      }
      // #pragma omp barrier // ===========================
      // #pragma omp critical
      // {
      //   res += localRes;
      // }
    }
    for (int i = 0; i < _mulVecCopies.size(); i++) res += _mulVecCopies[i];
#else
    for (unsigned int x = 0; x < mat.size(); x++) {
      res[mat[x].row()] += mat[x].value() * v[mat[x].col()];
    }
#endif
    return res;
  }

  // this.transpose() * v
  VECTOR transposeTimes(const VECTOR& v) const {
    assert(this->rows() == v.size());
    VECTOR res(this->cols());
    res.setZero();
    const vector<TRIPLET>& mat = this->matrix();

    if (mat.size() == 0) return res;

#if USING_OPENMP
#pragma omp parallel

    {
      VECTOR localRes(this->cols());
      localRes.setZero();

#pragma omp for schedule(static)
      for (unsigned int x = 0; x < mat.size(); x++) {
        localRes[mat[x].col()] += mat[x].value() * v[mat[x].row()];
        // res[mat[x].row()] += mat[x].value() * v[mat[x].col()];
      }
#pragma omp barrier  // ===========================
#pragma omp critical
      { res += localRes; }
    }
#else
    for (unsigned int x = 0; x < mat.size(); x++) {
      res[mat[x].col()] += mat[x].value() * v[mat[x].row()];
    }
#endif
    return res;
  }

  COO_MATRIX& operator*=(const Real& alpha) {
#if USING_OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (unsigned int x = 0; x < _matrix.size(); x++)
      _matrix[x].value() *= alpha;
    return *this;
  }
  COO_MATRIX& operator+=(const COO_MATRIX& A) {
    // const
    copy(A.matrix().begin(), A.matrix().end(), back_inserter(_matrix));
    _rows = _rows > A.rows() ? _rows : A.rows();
    _cols = _cols > A.cols() ? _cols : A.cols();
    return *this;
  }
  COO_MATRIX& operator+=(const MATRIX& A) {
    add(A, 0, 0);
    return *this;
  }
  COO_MATRIX& operator-=(const COO_MATRIX& A) {
    const vector<TRIPLET>& other = A.matrix();
    for (unsigned int x = 0; x < other.size(); x++)
      _matrix.push_back(
          TRIPLET(other[x].row(), other[x].col(), -other[x].value()));
    _rows = _rows > A.rows() ? _rows : A.rows();
    _cols = _cols > A.cols() ? _cols : A.cols();
    return *this;
  }

  COO_MATRIX transpose() {
    COO_MATRIX trans(_cols, _rows);
    for (unsigned int x = 0; x < _matrix.size(); x++)
      trans.add(_matrix[x].value(), _matrix[x].col(), _matrix[x].row());
    return trans;
  }

  void offset(int row, int col) {
    for (unsigned int x = 0; x < _matrix.size(); x++) {
      _matrix[x].row() += row;
      _matrix[x].col() += col;
    }
  }

  Real sum2() {
    Real sum = 0;
    for (unsigned int x = 0; x < _matrix.size(); x++) {
      sum += _matrix[x].value() * _matrix[x].value();
    }
    return sum;
  }
  Real sum() {
    Real sum = 0;
    for (unsigned int x = 0; x < _matrix.size(); x++) {
      sum += _matrix[x].value();
    }
    return sum;
  }

  // denseMat * this
  MATRIX leftMult(const MATRIX& denseMat) const {
    assert(denseMat.cols() == _rows);
    MATRIX ret(denseMat.rows(), _cols);
    ret.setZero();

    const vector<TRIPLET>& triplet = this->matrix();

    for (unsigned int x = 0; x != triplet.size(); x++) {
      int k = triplet[x].row();
      int j = triplet[x].col();
      ret.col(j) += denseMat.col(k) * triplet[x].value();
    }
    return ret;
  };

  // this * denseMat
  MATRIX rightMult(const MATRIX& denseMat) const {
    assert(_cols == denseMat.rows());
    MATRIX ret(_rows, denseMat.cols());
    ret.setZero();

    const vector<TRIPLET>& triplet = this->matrix();

    for (unsigned int x = 0; x != triplet.size(); x++) {
      int i = triplet[x].row();
      int k = triplet[x].col();
      ret.row(i) += triplet[x].value() * denseMat.row(k);
    }
    return ret;
  }

 protected:
  vector<TRIPLET> _matrix;
  VECTOR _diag;

  vector<VECTOR> _mulVecCopies;

  int _rows;
  int _cols;

  struct orderByIndexRowFirst {
    bool operator()(TRIPLET const& a, TRIPLET const& b) const {
      return a.row() < b.row();
    }
  } sortFunc;
};

#endif
