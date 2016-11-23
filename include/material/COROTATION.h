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
#ifndef COROTATION_H
#define COROTATION_H

#include <SETTINGS.h>
#include <iostream>

using namespace std;

//////////////////////////////////////////////////////////////////////
// Co-rotational material class
//////////////////////////////////////////////////////////////////////
class COROTATION {
 public:
  COROTATION(Real lambda, Real mu);
  ~COROTATION(){};

  Real lambda() { return _lambda; }
  Real mu() { return _mu; }

  COROTATION* copy();

  void init(const MATRIX3& F);

  Real strainEnergy();

  MATRIX3 firstPiolaKirchhoff();

  MATRIX3 firstPKDifferential(const MATRIX3& dF);

  MATRIX3& U() { return _U; };
  MATRIX3& V() { return _V; };
  MATRIX3& Fhat() { return _Fhat; };

 private:
  Real _lambda;
  Real _mu;
  static const MATRIX3 _I3;
  MATRIX3 _U, _Fhat, _V, _R, _S, _L;
  string _materialName;

  void decomposeF(const MATRIX3& F, MATRIX3& U, MATRIX3& Fhat, MATRIX3& V,
                  MATRIX3& R, MATRIX3& S, MATRIX3& L);

  void svd(const MATRIX3& F, MATRIX3& U, MATRIX3& Fhat, MATRIX3& V);
  void removeUReflection(MATRIX3& U, MATRIX3& Fhat);
  void orthogonalizeU(MATRIX3& U, MATRIX3& Fhat);
};

#endif
