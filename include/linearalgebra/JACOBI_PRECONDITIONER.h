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
#ifndef JACOBI_PRECONDITIONER_H
#define JACOBI_PRECONDITIONER_H
#include <SETTINGS.h>
#include <iostream>
#if USING_OPENMP
#include <omp.h>
#endif

using namespace ::std;
class JACOBI_PRECONDITIONER {
 public:
  JACOBI_PRECONDITIONER(VECTOR& diagonal) : _diag(diagonal) {}
  ~JACOBI_PRECONDITIONER(){};

  void init() {}
  void inverseSolve(const VECTOR& b, VECTOR& x) {
    x = (b.array() * _diag.array()).matrix();
  }
  void solve(const VECTOR& b, VECTOR& x) {
    x = (b.array() / _diag.array()).matrix();
  }

 private:
  const VECTOR& _diag;
};
#endif