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
#ifndef CONJUGATE_RESIDUAL_H
#define CONJUGATE_RESIDUAL_H

#include <SETTINGS.h>
#include <iostream>

using namespace ::std;

template <class PRECONDITIONER>
class CONJUGATE_RESIDUAL {
 public:
  CONJUGATE_RESIDUAL(PRECONDITIONER* preconditioner = NULL);
  ~CONJUGATE_RESIDUAL();

  int& maxIteration() { return _maxIteration; };
  Real& eps() { return _eps; };

  bool solvePCR(SpMat& A, const VECTOR& b, VECTOR& x);

  bool solveCR(const VECTOR& b, VECTOR& x,
               void (*getMatVecMult)(void*, const VECTOR&, VECTOR&), void* app);
  bool solvePCR(const VECTOR& b, VECTOR& x,
                void (*getMatVecMult)(void*, const VECTOR&, VECTOR&),
                void* app);

  vector<int>& numberOfIterations() { return _numberOfIterations; };

 protected:
  void parallelMatVecMult(SpMat& A, const VECTOR& b, VECTOR& x);

  PRECONDITIONER* _preconditioner;
  int _maxIteration;
  Real _eps;

  int _totalFrames;
  int _totalIterations;
  vector<int> _numberOfIterations;

  VECTOR _residual;
  VECTOR _direction;
  VECTOR _q;
  VECTOR _s;
  VECTOR _Ap;
  VECTOR _tmp;
};

#include "CONJUGATE_RESIDUAL.inl"
#endif