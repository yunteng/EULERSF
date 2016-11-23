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
#include <linearalgebra/MATRIX_UTIL.h>
#include <material/COROTATION.h>
#include <Eigen/SVD>

const MATRIX3 COROTATION::_I3 = MATRIX3::Identity();
//////////////////////////////////////////////////////////////////////
// Constructor for COROTATION
//////////////////////////////////////////////////////////////////////
COROTATION::COROTATION(Real lambda, Real mu) : _lambda(lambda), _mu(mu) {
  _materialName.assign("COROTATION");
}

//////////////////////////////////////////////////////////////////////
// make a copy
//////////////////////////////////////////////////////////////////////
COROTATION* COROTATION::copy() {
  COROTATION* material = new COROTATION(_lambda, _mu);
  return material;
}

void COROTATION::init(const MATRIX3& F) {
  decomposeF(F, _U, _Fhat, _V, _R, _S, _L);
}

// E = mu * ||F - R||_F + lambda / 2 * (S - I)
Real COROTATION::strainEnergy() {
  MATRIX3 E = _S - _I3;

  Real energy = _mu * (E.transpose() * E).trace() +
                0.5 * _lambda * (E.trace() * E.trace());

  return energy;
}
//////////////////////////////////////////////////////////////////////
// implementation of first PK stress tensor w.r.t. the singular values
//////////////////////////////////////////////////////////////////////
MATRIX3 COROTATION::firstPiolaKirchhoff() {
  MATRIX3 E = _S - _I3;
  MATRIX3 tmp = 2.0 * _mu * E;
  Real tr = _lambda * trace(E);
  tmp(0, 0) += tr;
  tmp(1, 1) += tr;
  tmp(2, 2) += tr;
  // tmp = (tmp.array() + _lambda * trace(E)).matrix();
  return _R * tmp;
}

MATRIX3 COROTATION::firstPKDifferential(const MATRIX3& dF) {
  MATRIX3 dFhat = _R.transpose() * dF;
  MATRIX3 dFsym = 0.5 * (dFhat + dFhat.transpose());
  MATRIX3 dFskew = 0.5 * (dFhat - dFhat.transpose());

  MATRIX3 dPsym =
      (2 * _mu * dFsym) + (_lambda * trace(dFsym) * MATRIX3::Identity());

  VEC3F f(-dFskew(1, 2) + dFskew(2, 1), -dFskew(2, 0) + dFskew(0, 2),
          -dFskew(0, 1) + dFskew(1, 0));

  Real tr = _lambda * dFhat.trace();
  dPsym(0, 0) += tr;
  dPsym(1, 1) += tr;
  dPsym(2, 2) += tr;

  MATRIX3 dPskew = MATRIX_UTIL::cross(_L * f);
  MATRIX3 deltaP = dPsym + dPskew;
  return _R * deltaP;
  ;
}

//////////////////////////////////////////////////////////////////////////////
// F = U * Fhat * VT
// R = U * VT
// S = V * Fhat * VT
//////////////////////////////////////////////////////////////////////////////
void COROTATION::decomposeF(const MATRIX3& F, MATRIX3& U, MATRIX3& Fhat,
                            MATRIX3& V, MATRIX3& R, MATRIX3& S, MATRIX3& L) {
  svd(F, U, Fhat, V);

  R = U * V.transpose();
  S = V * Fhat * V.transpose();

  MATRIX3 Ld = (_mu * _I3) +
               (_lambda * (Fhat - _I3).trace() - 2 * _mu) *
                   ((Fhat.trace() * _I3) - Fhat).inverse();

  Ld(0, 0) = Ld(0, 0) > 0 ? Ld(0, 0) : 0;
  Ld(1, 1) = Ld(1, 1) > 0 ? Ld(1, 1) : 0;
  Ld(2, 2) = Ld(2, 2) > 0 ? Ld(2, 2) : 0;

  L = V * Ld * V.transpose();
}

//////////////////////////////////////////////////////////////////////////////
// Correct both Fhat and U if U contains a reflection
//////////////////////////////////////////////////////////////////////////////
void COROTATION::removeUReflection(MATRIX3& U, MATRIX3& Fhat) {
  // find the smallest singular value
  int smallest = (Fhat(0, 0) < Fhat(1, 1)) ? 0 : 1;
  smallest = (Fhat(smallest, smallest) < Fhat(2, 2)) ? smallest : 2;

  // negate it, and the corresponding column in U
  Fhat(smallest, smallest) *= -1.0;
  U(0, smallest) *= -1.0;
  U(1, smallest) *= -1.0;
  U(2, smallest) *= -1.0;
}

//////////////////////////////////////////////////////////////////////////////
// Orthogonalize U if one of the singular values is near zero
//////////////////////////////////////////////////////////////////////////////
void COROTATION::orthogonalizeU(MATRIX3& U, MATRIX3& Fhat) {
  // record the sign of U (ie whether it is a reflection) so that
  // the sign can be preserved.
  //
  // If it contains a reflection, it will be corrected in
  // removeUReflection, so make sure we don't double-handle it here.
  Real signU = (U.determinant() >= 0.0) ? 1.0 : -1.0;

  Real smallEigenvalue = 1e-4;

  // see if an eigenvalue is really small
  bool invalid[] = {false, false, false};
  int totalInvalid = 0;
  for (int x = 0; x < 3; x++)
    if (Fhat(x, x) < smallEigenvalue) {
      invalid[x] = true;
      totalInvalid++;
    }

  // if none are bad, do nothing
  if (totalInvalid == 0) return;

  // correct U according to the number of small eigenvalues
  if (totalInvalid == 1) {
    // figure out which column is bad
    int whichInvalid = 0;
    if (invalid[1]) whichInvalid = 1;
    if (invalid[2]) whichInvalid = 2;

    // get the other two columns
    VEC3F col0 = U.col((whichInvalid + 1) % 3);
    VEC3F col1 = U.col((whichInvalid + 2) % 3);

    // compute an orthogonal vector
    VEC3F col2 = col0.cross(col1);
    col2.normalize();

    // store it in U
    U(0, whichInvalid) = col2[0];
    U(1, whichInvalid) = col2[1];
    U(2, whichInvalid) = col2[2];

    // see if we changed any reflections
    if (signU * U.determinant() < 0.0) {
      U(0, whichInvalid) = -col2[0];
      U(1, whichInvalid) = -col2[1];
      U(2, whichInvalid) = -col2[2];
    }

    return;
  } else if (totalInvalid == 2) {
    // find the good vector
    int whichValid = 0;
    if (!invalid[1]) whichValid = 1;
    if (!invalid[2]) whichValid = 2;
    VEC3F validVec = U.col(whichValid);

    // find the smallest entry in the good vector
    int smallest = (fabs(validVec[0]) < fabs(validVec[1])) ? 0 : 1;
    smallest = (fabs(validVec[2]) < fabs(validVec[smallest])) ? 2 : smallest;

    // compute something arbitrarily orthogonal to
    // the good vector
    VEC3F secondValid;
    int next = (smallest + 1) % 3;
    int nextnext = (smallest + 2) % 3;
    secondValid[smallest] = 0.0;
    secondValid[next] = -validVec[nextnext];
    secondValid[nextnext] = validVec[next];

    // compute a third good vector
    VEC3F thirdValid = validVec ^ secondValid;
    thirdValid.normalize();

    // copy the good vectors into U
    next = (whichValid + 1) % 3;
    nextnext = (whichValid + 2) % 3;
    U(0, next) = secondValid[0];
    U(1, next) = secondValid[1];
    U(2, next) = secondValid[2];
    U(0, nextnext) = thirdValid[0];
    U(1, nextnext) = thirdValid[1];
    U(2, nextnext) = thirdValid[2];

    // see if we changed any reflections
    if (signU * U.determinant() < 0.0) {
      // negate the one with the smallest F, which is the one
      // removeUReflection will be negating again
      if (Fhat(next, next) < Fhat(nextnext, nextnext)) {
        U(0, next) = -secondValid[0];
        U(1, next) = -secondValid[1];
        U(2, next) = -secondValid[2];
      } else {
        U(0, nextnext) = -thirdValid[0];
        U(1, nextnext) = -thirdValid[1];
        U(2, nextnext) = -thirdValid[2];
      }
    }
    return;
  }

  // all three are bad, just use identity
  U = _I3;
}

//////////////////////////////////////////////////////////////////////////////
// Get the SVD of F
//////////////////////////////////////////////////////////////////////////////
void COROTATION::svd(const MATRIX3& F, MATRIX3& U, MATRIX3& Fhat, MATRIX3& V) {
  // compute the SVD
  MATRIX3 Fnormal3 = F.transpose() * F;

  Eigen::JacobiSVD<MATRIX3> eigenSystem(Fnormal3, Eigen::ComputeFullV);
  VEC3F eigenvalues = eigenSystem.singularValues();
  V = eigenSystem.matrixV();

  VEC3F oEig = eigenvalues;
  MATRIX3 oV = V;

  // populate the diagonal
  Fhat.setZero();
  Fhat(0, 0) = sqrt(eigenvalues(0));
  Fhat(1, 1) = sqrt(eigenvalues(1));
  Fhat(2, 2) = sqrt(eigenvalues(2));

#ifdef _WIN32
  if (_isnan(Fhat(0, 0))) Fhat(0, 0) = 0.0;
  if (_isnan(Fhat(1, 1))) Fhat(1, 1) = 0.0;
  if (_isnan(Fhat(2, 2))) Fhat(2, 2) = 0.0;
#else
  if (std::isnan(Fhat(0, 0))) Fhat(0, 0) = 0.0;
  if (std::isnan(Fhat(1, 1))) Fhat(1, 1) = 0.0;
  if (std::isnan(Fhat(2, 2))) Fhat(2, 2) = 0.0;
#endif

  // if V is a reflection, multiply a column by -1
  // to ensure it is a rotation
  Real detV = V.determinant();
  if (detV < 0.0) {
    V(0, 0) *= -1.0;
    V(1, 0) *= -1.0;
    V(2, 0) *= -1.0;
  }

  // compute U
  U.setZero();
  U(0, 0) = (Fhat(0, 0) > 0.0f) ? 1.0f / Fhat(0, 0) : 0.0f;
  U(1, 1) = (Fhat(1, 1) > 0.0f) ? 1.0f / Fhat(1, 1) : 0.0f;
  U(2, 2) = (Fhat(2, 2) > 0.0f) ? 1.0f / Fhat(2, 2) : 0.0f;
  U = F * V * U;
  orthogonalizeU(U, Fhat);

  // correct any reflections in U to ensure it is a rotation
  if (F.determinant() < 0.0) removeUReflection(U, Fhat);

///////////////////////////////////////////////////////////////
// All error checking starting here.
///////////////////////////////////////////////////////////////
#if MATERIAL_DEBUG
  // basic error checking

  if (fabs(U.determinant() - 1.0) > 1e-4 ||
      fabs(V.determinant() - 1.0) > 1e-4) {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " U and V of SVD have drifted from pure rotation!" << endl;
    cout << " det(U): " << U.determinant() << endl;
    cout << " det(V): " << V.determinant() << endl;
  }

  // sanity check that this decomposition does hand back F when combined
  MATRIX3 sanity = U * Fhat * V.transpose();
  for (int x = 0; x < 3; x++)
    for (int y = 0; y < 3; y++) {
      Real diff = sanity(x, y) - F(x, y);
      // if (F(x,y) > 0.0)
      // diff = diff / F(x,y);
      if (diff > 1e-4) {
        cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : "
             << endl;
        cout << "SVD and original matrix differ!" << endl;
        cout << " relative diff: " << diff << " sanity: " << sanity(x, y)
             << " original: " << F(x, y) << endl;
        cout << " SVD product:\n" << sanity << endl;
        cout << " original:\n" << F << endl;
        cout << " normal F:\n" << Fnormal3 << endl;
        cout << " eigenvalues: " << eigenvalues.transpose() << endl;
        cout << " eigenvectors: \n" << oV << endl;
        cout << " eigenproduct:\n " << oV * oEig.asDiagonal() * oV.transpose();
        // exit(0);
      }
    }
#endif
}
