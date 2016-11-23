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
#include <util/SIMPLE_PARSER.h>
#if USING_OPENMP
#include <omp.h>
#endif

template <class PRECONDITIONER>
CONJUGATE_RESIDUAL<PRECONDITIONER>::CONJUGATE_RESIDUAL(
    PRECONDITIONER* preconditioner)
    : _preconditioner(preconditioner), _totalFrames(0), _totalIterations(0) {
  _maxIteration = SIMPLE_PARSER::getInt("cr max iterations", 10);
  _eps = SIMPLE_PARSER::getFloat("cr eps", 1e-8);
}
template <class PRECONDITIONER>
CONJUGATE_RESIDUAL<PRECONDITIONER>::~CONJUGATE_RESIDUAL() {}
template <class PRECONDITIONER>
void CONJUGATE_RESIDUAL<PRECONDITIONER>::parallelMatVecMult(SpMat& A,
                                                            const VECTOR& b,
                                                            VECTOR& x) {
  x.conservativeResize(A.rows());
  x.setZero();
#if USING_OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (int k = 0; k < A.rows(); ++k)
    for (SpMat::InnerIterator it(A, k); it; ++it) {
      x[k] += it.value() * b[it.col()];
    }
}

template <class PRECONDITIONER>
bool CONJUGATE_RESIDUAL<PRECONDITIONER>::solveCR(
    const VECTOR& b, VECTOR& x,
    void (*getMatVecMult)(void*, const VECTOR&, VECTOR&), void* app) {
  x.conservativeResize(b.size());
  x.setZero();

  (*getMatVecMult)(app, x, _q);
  _direction = _residual = b - _q;

  (*getMatVecMult)(app, _residual, _s);

  _Ap = _s;

  Real deltaNew = _residual.dot(_s);

  Real delta0 = deltaNew;

  bool converged = false;
  int i = 0;
  for (; i < _maxIteration && !converged; i++) {
    // (*getMatVecMult)(app, _direction, _q);

    Real alpha = _Ap.squaredNorm();
    if (fabs(alpha) > 0.0) alpha = deltaNew / alpha;

    x += alpha * _direction;

    if (i > 0 && i % 25 == 0) {
      (*getMatVecMult)(app, x, _residual);
      _residual *= -1;
      _residual += b;
    } else
      _residual -= alpha * _Ap;

    Real deltaOld = deltaNew;

    (*getMatVecMult)(app, _residual, _s);
    deltaNew = _residual.dot(_s);

    Real beta = deltaNew / deltaOld;

    _direction *= beta;
    _direction += _residual;

    _Ap *= beta;
    _Ap += _s;

    if (abs(deltaNew) <= abs(_eps * delta0)) converged = true;
  }

  if (converged) {
    cout << "Conjugate Residual solver converged in " << i << " iterations"
         << endl;
    _totalIterations++;
    _totalFrames += i;
    _numberOfIterations.push_back(_totalFrames);
  } else
    cout << "Conjugate Residual solver didn't converge in " << _maxIteration
         << " iterations, residual " << deltaNew << endl;
  return converged;
}

template <class PRECONDITIONER>
bool CONJUGATE_RESIDUAL<PRECONDITIONER>::solvePCR(
    const VECTOR& b, VECTOR& x,
    void (*getMatVecMult)(void*, const VECTOR&, VECTOR&), void* app) {
  if (x.size() != b.size()) {
    x.conservativeResize(b.size());
    x.setZero();
  }

  (*getMatVecMult)(app, x, _q);
  _residual = b - _q;

  _preconditioner->solve(_residual, _direction);
  _residual = _direction;

  (*getMatVecMult)(app, _residual, _s);
  _Ap = _s;

  Real deltaNew = _residual.dot(_s);

  // Real initNorm = _residual.norm();

  Real delta0 = deltaNew;
  Real minDelta = delta0;

  Real bNorm = b.norm();

  bool converged = false;
  int i = 0;
  for (; i < _maxIteration && !converged; i++) {
    _preconditioner->solve(_Ap, _q);

    Real alpha = _Ap.dot(_q);
    if (fabs(alpha) > 0.0) alpha = deltaNew / alpha;

    x += alpha * _direction;

    _residual -= alpha * _q;

    (*getMatVecMult)(app, _residual, _s);

    Real deltaOld = deltaNew;

    deltaNew = _residual.dot(_s);

    Real beta = deltaNew / deltaOld;

    _direction *= beta;
    _direction += _residual;

    _Ap *= beta;
    _Ap += _s;

    if (i > 0 && i % 5 == 0) {
      (*getMatVecMult)(app, x, _tmp);
      _tmp -= b;
      // initNorm = _tmp.norm();
      Real tmpNorm = _tmp.norm();
      if (tmpNorm < _eps * bNorm) {
        converged = true;
        break;
      }
    }
  }

  // if(converged) {
  //   cout << "absolute residual (Ax-b).norm " << initNorm << endl;
  //   cout << "relative residual (Ax-b).norm / b.norm " << initNorm / bNorm <<
  //   endl;
  // }

  if (converged) {
    _totalIterations++;
    _totalFrames += i;
    _numberOfIterations.push_back(i);
    cout << "Preconditioned Conjugate Residual solver converged in " << i
         << " iterations " << endl;
  } else
    cout << "Preconditioned Conjugate Residual solver didn't converge in "
         << _maxIteration << " iterations, relative residual " << deltaNew
         << " absolute residual " << _residual.norm() << endl;

  return converged;
}

template <class PRECONDITIONER>
bool CONJUGATE_RESIDUAL<PRECONDITIONER>::solvePCR(SpMat& A, const VECTOR& b,
                                                  VECTOR& x) {
  x.conservativeResize(b.size());
  x.setZero();

  _q.conservativeResize(b.size());
  _q.setZero();
  _residual = b;

  _preconditioner->solve(_residual, _direction);
  _residual = _direction;

  parallelMatVecMult(A, _residual, _s);

  _Ap = _s;

  Real deltaNew = _residual.dot(_s);

  // Real initNorm = _residual.norm();

  Real delta0 = deltaNew;
  Real minDelta = delta0;

  Real bNorm = b.norm();

  bool converged = false;
  int i = 0;
  for (; i < _maxIteration && !converged; i++) {
    _preconditioner->solve(_Ap, _q);

    Real alpha = _Ap.dot(_q);
    if (fabs(alpha) > 0.0) alpha = deltaNew / alpha;

    x += alpha * _direction;

    _residual -= alpha * _q;

    parallelMatVecMult(A, _residual, _s);

    Real deltaOld = deltaNew;

    deltaNew = _residual.dot(_s);

    Real beta = deltaNew / deltaOld;

    _direction *= beta;
    _direction += _residual;

    _Ap *= beta;
    _Ap += _s;

    if (i > 0 && i % 5 == 0) {
      parallelMatVecMult(A, x, _tmp);
      _tmp -= b;
      // initNorm = _tmp.norm();
      Real tmpNorm = _tmp.norm();
      if (tmpNorm < _eps * bNorm) {
        converged = true;
        break;
      }
    }
  }

  // if(converged) {
  //   cout << "absolute residual (Ax-b).norm " << initNorm << endl;
  //   cout << "relative residual (Ax-b).norm / b.norm " << initNorm / bNorm <<
  //   endl;
  // }

  if (converged) {
    _totalIterations++;
    _totalFrames += i;
    _numberOfIterations.push_back(i);
    cout << "Preconditioned Conjugate Residual solver converged in " << i
         << " iterations " << endl;
  } else
    cout << "Preconditioned Conjugate Residual solver didn't converge in "
         << _maxIteration << " iterations, relative residual " << deltaNew
         << " absolute residual " << _residual.norm() << endl;

  return converged;
}
