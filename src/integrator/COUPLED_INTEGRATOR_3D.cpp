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
#include <integrator/COUPLED_INTEGRATOR_3D.h>
#include <util/IO.h>
#include <util/SIMPLE_PARSER.h>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseLU>
#if USING_OPENMP
#include <omp.h>
#endif

#if _WIN32
#include <gl/glut.h>
#elif USING_OSX
#include <GLUT/glut.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#endif

COUPLED_INTEGRATOR_3D::COUPLED_INTEGRATOR_3D(int xRes, int yRes, int zRes,
                                             Real dh, VEC3F origin)
    : _xRes(xRes),
      _yRes(yRes),
      _zRes(zRes),
      _dh(dh),
      _origin(origin),
      _defaultDt(SIMPLE_PARSER::getFloat("dt", 0.01)),
      _dt(_defaultDt),
      _fluidDensity(1),
      _fluidViscocity(0),
      _preconditioner(_systemMatrixDiag),
      _solver(&_preconditioner),
      _divfreeSolverPreconditioner(_JMInvJTDiag),
      _divFreeSolver(&_divfreeSolverPreconditioner) {
  _slabSize = _xRes * _yRes;
  _totalCells = _xRes * _yRes * _zRes;
  _totalActiveCells = (_xRes - 2) * (_yRes - 2) * (_zRes - 2);
  _dhInv = 1.0 / _dh;

  _solidMass.resize(_xRes, _yRes, _zRes, _dh, _origin);
  _intensity.resizeLike(_solidMass);
  _workspace.resizeLike(_solidMass);

  _nodeSolidVolumeFraction.resizeLike(_solidMass);
  _quarterSolidVolumeFractions.resize(_totalCells * 8);

  _velocity.resize(_xRes, _yRes, _zRes, _dh, _origin);
  _vecWorkspace.resizeLike(_velocity);
  _vorticity.resizeLike(_velocity);

  _vorticityEps = SIMPLE_PARSER::getFloat("vorticity eps", 2.0);
  int maxRes = max(max(_xRes, _yRes), _zRes);
  Real scaling = 64.0f / maxRes;
  scaling = (scaling < 1.0f) ? 1.0f : scaling;
  _vorticityEps /= scaling;

  _solver.eps() = SIMPLE_PARSER::getFloat("implicit solver eps", 1e-8);
  _divFreeSolver.eps() = SIMPLE_PARSER::getFloat("div free solver eps", 1e-8);

  cout << "implicit solver eps " << _solver.eps() << endl;
  cout << "div free solver eps " << _divFreeSolver.eps() << endl;
}
COUPLED_INTEGRATOR_3D::~COUPLED_INTEGRATOR_3D() {
  for (DEFOR_SOLID_3D* solid : _solids) {
    if (solid) {
      delete solid;
      solid = NULL;
    }
  }
}
// Compute the time step size given a CFL condition
Real COUPLED_INTEGRATOR_3D::CFL(Real alpha) {
  Real dt = 1e9;
  int index = 0;
  for (int z = 1; z < _zRes - 1; z++)
    for (int y = 1; y < _yRes - 1; y++)
      for (int x = 1; x < _xRes - 1; x++, index++) {
        int index = z * _slabSize + y * _xRes + x;
        if (_solidMass[index] > 0) {
          Real vmax = _velocity[index].array().abs().maxCoeff() * _dhInv;
          dt = min(dt, alpha / vmax);
        }
      }
  dt = min(dt, _defaultDt);
  return dt;
}

void COUPLED_INTEGRATOR_3D::init() {
  for (int i = 0; i < _solids.size(); i++) {
    cout << "initializing solid " << i << endl;
    _solids[i]->init();
  }

  for (int i = 0; i < _solids.size(); i++) {
    for (int j = 0; j < _solids.size(); j++) {
      if (i != j) _solids[i]->addNeighbor(_solids[j]);
    }
  }

  gatherSolidMasses();
#if USING_OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (int i = 0; i < _solids.size(); i++)
    _solids[i]->recoverVelocityField(_velocity);

  _gridIdxToSolutionIdx.resize(_totalCells, -1);

  int idx = 0;
  for (int z = 1; z < _zRes - 1; z++)
    for (int y = 1; y < _yRes - 1; y++)
      for (int x = 1; x < _xRes - 1; x++) {
        int index = z * _slabSize + y * _xRes + x;
        _gridIdxToSolutionIdx[index] = idx++;
      }
#if !HYBRID_DIV_CONSTRAINT
  computeDivergenceConstraintJacobian();
#endif
}

void COUPLED_INTEGRATOR_3D::computeDivergenceConstraintJacobian() {
  static int offsets[8] = {0,
                           1,
                           _xRes,
                           _xRes + 1,
                           _slabSize,
                           _slabSize + 1,
                           _slabSize + _xRes,
                           _slabSize + _xRes + 1};
  static int faceCentersOffsets[3] = {1, _xRes, _slabSize};
  static int edgeCentersOffsets[12] = {_xRes, -_xRes, _slabSize, -_slabSize,
                                       1,     -1,     _slabSize, -_slabSize,
                                       1,     -1,     _xRes,     -_xRes};
  static int vertexOffsets[12] = {
      _xRes + _slabSize,  _xRes - _slabSize, -_xRes + _slabSize,
      -_xRes - _slabSize,  // y+-1, z+-1
      1 + _slabSize,      1 - _slabSize,     -1 + _slabSize,
      -1 - _slabSize,  // x+-1,z+-1
      1 + _xRes,          1 - _xRes,         -1 + _xRes,
      -1 - _xRes  // x+-1,y+-1
  };

  vector<TRIPLET> jacobiTripletList;

  int cnt = 0;
#if HYBRID_DIV_CONSTRAINT
  _isConstrained.resize(_totalCells, false);
  _divergenceConstraintsRHS.conservativeResize(_totalCells);
  _divergenceConstraintsRHS.setZero();
  for (int z = 2; z < _zRes - 2; z++)
    for (int y = 2; y < _yRes - 2; y++)
      for (int x = 2; x < _xRes - 2; x++) {
        int index = z * _slabSize + y * _xRes + x;
        if (_solidMass[index] == 0) continue;
        for (int i = 0; i < 8; i++) _isConstrained[index - offsets[i]] = true;

        // for each direction
        for (int i = 0; i < 3; i++) {
          int idx = _gridIdxToSolutionIdx[index + faceCentersOffsets[i]];
          if (idx >= 0)
            jacobiTripletList.push_back(TRIPLET(cnt, idx * 3 + i, 4));
          else
            _divergenceConstraintsRHS[cnt] -=
                4 * _velocity[index + faceCentersOffsets[i]][i];
          idx = _gridIdxToSolutionIdx[index - faceCentersOffsets[i]];
          if (idx >= 0)
            jacobiTripletList.push_back(TRIPLET(cnt, idx * 3 + i, -4));
          else
            _divergenceConstraintsRHS[cnt] +=
                4 * _velocity[index - faceCentersOffsets[i]][i];

          for (int j = i * 4; j < (i + 1) * 4; j++) {
            idx = _gridIdxToSolutionIdx[index + edgeCentersOffsets[j] +
                                        faceCentersOffsets[i]];
            if (idx >= 0)
              jacobiTripletList.push_back(TRIPLET(cnt, idx * 3 + i, 2));
            else
              _divergenceConstraintsRHS[cnt] -=
                  2 * _velocity[index + edgeCentersOffsets[j] +
                                faceCentersOffsets[i]][i];

            idx = _gridIdxToSolutionIdx[index + edgeCentersOffsets[j] -
                                        faceCentersOffsets[i]];
            if (idx >= 0)
              jacobiTripletList.push_back(TRIPLET(cnt, idx * 3 + i, -2));
            else
              _divergenceConstraintsRHS[cnt] +=
                  2 * _velocity[index + edgeCentersOffsets[j] -
                                faceCentersOffsets[i]][i];

            idx = _gridIdxToSolutionIdx[index + vertexOffsets[j] +
                                        faceCentersOffsets[i]];
            if (idx >= 0)
              jacobiTripletList.push_back(TRIPLET(cnt, idx * 3 + i, 1));
            else
              _divergenceConstraintsRHS[cnt] -=
                  _velocity[index + vertexOffsets[j] + faceCentersOffsets[i]]
                           [i];

            idx = _gridIdxToSolutionIdx[index + vertexOffsets[j] -
                                        faceCentersOffsets[i]];
            if (idx >= 0)
              jacobiTripletList.push_back(TRIPLET(cnt, idx * 3 + i, -1));
            else
              _divergenceConstraintsRHS[cnt] +=
                  _velocity[index + vertexOffsets[j] - faceCentersOffsets[i]]
                           [i];
          }
        }
        cnt++;
      }
#endif

  FIELD_3Dc& solidCell = _solids[0]->_valid;
  int index = 0;
  for (int z = 0; z < _zRes - 1; z++, index += _xRes)
    for (int y = 0; y < _yRes - 1; y++, index++)
      for (int x = 0; x < _xRes - 1; x++, index++) {
        assert(index == z * _slabSize + y * _xRes + x);
#if HYBRID_DIV_CONSTRAINT
        if (_isConstrained[index]) continue;
#endif
        for (int i = 0; i < 8; i++) {
          int idx = _gridIdxToSolutionIdx[index + offsets[i]];
          int signX = i % 2 == 0 ? -1 : 1;
          int signY = i % 4 <= 1 ? -1 : 1;
          int signZ = i < 4 ? -1 : 1;
          if (idx >= 0) {
            jacobiTripletList.push_back(TRIPLET(cnt, idx * 3, signX));
            jacobiTripletList.push_back(TRIPLET(cnt, idx * 3 + 1, signY));
            jacobiTripletList.push_back(TRIPLET(cnt, idx * 3 + 2, signZ));
          }
        }
#if HYBRID_DIV_CONSTRAINT
        else {
          _divergenceConstraintsRHS[cnt] -=
              signX * _velocity[index + offsets[i]][0];
          _divergenceConstraintsRHS[cnt] -=
              signY * _velocity[index + offsets[i]][1];
          _divergenceConstraintsRHS[cnt] -=
              signZ * _velocity[index + offsets[i]][2];
        }
#endif
        cnt++;
      }

  _divergenceConstraintJacobian.resize(cnt, _totalActiveCells * 3);
  _divergenceConstraintJacobian.setFromTriplets(jacobiTripletList.begin(),
                                                jacobiTripletList.end());
  _divergenceConstraintJacobian.makeCompressed();

#if HYBRID_DIV_CONSTRAINT
  _divergenceConstraintsRHS.conservativeResize(cnt);
#endif

#if !HYBRID_DIV_CONSTRAINT
  _JJTDiagEntries.clear();
  for (int k = 0; k < _divergenceConstraintJacobian.rows(); ++k) {
    vector<pair<int, Real> > entries;
    for (SpMat::InnerIterator it(_divergenceConstraintJacobian, k); it; ++it) {
      entries.push_back(make_pair(it.col(), it.value() * it.value()));
    }
    _JJTDiagEntries.push_back(entries);
  }
  _JT = _divergenceConstraintJacobian.transpose();
  _JT.makeCompressed();
#endif
}

void COUPLED_INTEGRATOR_3D::computeDivergenceConstraintRHS() {
  _divergenceConstraintsRHS.conservativeResize(
      _divergenceConstraintJacobian.rows());

  static int offsets[8] = {0,
                           1,
                           _xRes,
                           _xRes + 1,
                           _slabSize,
                           _slabSize + 1,
                           _slabSize + _xRes,
                           _slabSize + _xRes + 1};

#if USING_OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (int z = 0; z < _zRes - 1; z++) {
    int cnt = (_yRes - 1) * (_xRes - 1) * z;
    int index = _slabSize * z;
    for (int y = 0; y < _yRes - 1; y++, index++)
      for (int x = 0; x < _xRes - 1; x++, index++) {
        _divergenceConstraintsRHS[cnt] = 0;
        for (int i = 0; i < 8; i++) {
          int idx = _gridIdxToSolutionIdx[index + offsets[i]];
          int signX = i % 2 == 0 ? -1 : 1;
          int signY = i % 4 <= 1 ? -1 : 1;
          int signZ = i < 4 ? -1 : 1;
          if (idx < 0) {
            _divergenceConstraintsRHS[cnt] -=
                signX * _velocity[index + offsets[i]][0];
            _divergenceConstraintsRHS[cnt] -=
                signY * _velocity[index + offsets[i]][1];
            _divergenceConstraintsRHS[cnt] -=
                signZ * _velocity[index + offsets[i]][2];
          }
        }
        cnt++;
      }
  }
}

void COUPLED_INTEGRATOR_3D::gatherSolidMasses() {
  static int offsets[8] = {0,
                           1,
                           _xRes,
                           _xRes + 1,
                           _slabSize,
                           _slabSize + 1,
                           _slabSize + _xRes,
                           _slabSize + _xRes + 1};

  TIMING_BREAKDOWN::tic();
  _solidMass.clear();
  _quarterSolidVolumeFractions.setZero();
  for (DEFOR_SOLID_3D* solid : _solids) {
    _solidMass += solid->_mass;
    _quarterSolidVolumeFractions += solid->_quarterVolumeFractions;
  }
  _nodeSolidVolumeFraction.clear();

  for (int i = 0; i < 8; i++) {
#if USING_OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int z = 1; z < _zRes - 1; z++) {
      int index = z * _slabSize + _xRes + 1;
      for (int y = 1; y < _yRes - 1; y++, index += 2)
        for (int x = 1; x < _xRes - 1; x++, index++) {
          _nodeSolidVolumeFraction[index + offsets[i]] +=
              _quarterSolidVolumeFractions[index * 8 + i];
        }
    }
  }
  _nodeSolidVolumeFraction *= 0.125;

  _gridIdxToSolidIdx.resize(_totalCells, -1);
  int idx = 0;
  int index = _slabSize + _xRes + 1;
  for (int z = 1; z < _zRes - 1; z++, index += 2 * _xRes)
    for (int y = 1; y < _yRes - 1; y++, index += 2)
      for (int x = 1; x < _xRes - 1; x++, index++) {
        assert(index == z * _slabSize + y * _xRes + x);
        if (_solidMass[index] > 0)
          _gridIdxToSolidIdx[index] = idx++;
        else
          _gridIdxToSolidIdx[index] = -1;
      }
  _totalSolidCells = idx;
  TIMING_BREAKDOWN::toc("gatherSolidMasses");
}

Real COUPLED_INTEGRATOR_3D::stepImplicit() {
  cout << "Check for collisions" << endl;
  TIMING_BREAKDOWN::tic();
#if USING_OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (int i = 0; i < _solids.size(); i++) _solids[i]->_bvh->refit();
  TIMING_BREAKDOWN::toc("Refit BVH");

  TIMING_BREAKDOWN::tic();
  bool proneToContact = false;
  for (int i = 0; i < _solids.size() && !proneToContact; i++) {
    if (_solids[i]->closeToWall()) {
      cout << "solid " << i << " is close to wall" << endl;
      proneToContact = true;
    }
  }
  if (!proneToContact) {
    for (int i = 0; i < _solids.size() - 1; i++)
      for (int j = i + 1; j < _solids.size(); j++) {
        if (_solids[i]->overlap(_solids[j])) {
          proneToContact = true;
          cout << "solid " << i << " is close to solid " << j << endl;
          i = _solids.size();
          break;
        }
      }
  }
  TIMING_BREAKDOWN::toc("Broad phase collision detection");
  if (proneToContact) {
    _dt = CFL(0.6);
  } else {
    _dt = CFL(1.8);
  }
  cout << "Use timestep " << _dt << endl;

  _mv.conservativeResize(_totalActiveCells * 3);
  _unitForce.conservativeResize(_totalActiveCells * 3);
  _unitForce.setZero();
  _nzMassVec.conservativeResize(_totalActiveCells * 3);
  _nzMassVecInv.conservativeResize(_totalActiveCells * 3);

  _solidUnitForce.conservativeResize(_totalSolidCells * 3);
  _solidUnitForce.setZero();
  _solidMv.conservativeResize(_totalSolidCells * 3);
  _solidMv.setZero();
  _dampedSolidMassVec.conservativeResize(_totalSolidCells * 3);
  _dampedSolidMassVec.setZero();

  vector<TRIPLET> tripletList;

  Real energy = 0;
  cout << "Computing solid stiffness matrix and material force ...";
  flush(cout);
  TIMING_BREAKDOWN::tic();
  for (int i = 0; i < _solids.size(); i++) {
    // cout << "solid " << i << endl;
    energy += _solids[i]->computeStiffnessMatrixAndMaterialForce(
        _solidUnitForce, tripletList, _gridIdxToSolidIdx);
    energy += _solids[i]->computeMV(_solidMv, _gridIdxToSolidIdx);
  }

  TIMING_BREAKDOWN::toc("Compute Solid Material Forces and Stiffness Matrix");
  cout << "Done." << endl;

  Real gridSize = 0.5f / _dh;
  _workspace.clear();

  static Real rayleighAlpha = SIMPLE_PARSER::getFloat("rayleigh alpha", 0);
  static Real rayleighBeta = SIMPLE_PARSER::getFloat("rayleigh beta", 0);

  Real solidMassDampCoeff = 1 + _dt * rayleighAlpha;

  // calculate vorticity
  if (_vorticityEps > 0) {
#if USING_OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int z = 1; z < _zRes - 1; z++)
      for (int y = 1; y < _yRes - 1; y++)
        for (int x = 1; x < _xRes - 1; x++) {
          int index = z * _slabSize + y * _xRes + x;
          if (_nodeSolidVolumeFraction[index] >= 1) continue;
          int up =
              _nodeSolidVolumeFraction[index + _xRes] >= 1 || y == _yRes - 2
                  ? index
                  : index + _xRes;
          int down = _nodeSolidVolumeFraction[index - _xRes] >= 1 || y == 1
                         ? index
                         : index - _xRes;
          Real dy = (up == index || down == index) ? _dhInv : gridSize;

          int out =
              _nodeSolidVolumeFraction[index + _slabSize] >= 1 || z == _zRes - 2
                  ? index
                  : index + _slabSize;
          int in = _nodeSolidVolumeFraction[index - _slabSize] >= 1 || z == 1
                       ? index
                       : index - _slabSize;

          Real dz = (out == index || in == index) ? _dhInv : gridSize;
          int right = _nodeSolidVolumeFraction[index + 1] >= 1 || x == _xRes - 2
                          ? index
                          : index + 1;
          int left = _nodeSolidVolumeFraction[index - 1] >= 1 || x == 1
                         ? index
                         : index - 1;
          Real dx = (right == index || left == index) ? _dhInv : gridSize;

          _vorticity[index][0] = ((_velocity[up][2] - _velocity[down][2]) +
                                  (-_velocity[out][1] + _velocity[in][1])) *
                                 dz;
          _vorticity[index][1] = ((_velocity[out][0] - _velocity[in][0]) +
                                  (-_velocity[right][2] + _velocity[left][2])) *
                                 dx;
          _vorticity[index][2] = ((_velocity[right][1] - _velocity[left][1]) +
                                  (-_velocity[up][0] + _velocity[down][0])) *
                                 dy;
          _workspace[index] = _vorticity[index].norm();
        }
  }

  static Real gravity = SIMPLE_PARSER::getFloat("gravity", 0);

  static Real bouyancy = SIMPLE_PARSER::getFloat("buoyancy", 0.1);

  TIMING_BREAKDOWN::tic();

#if USING_OPENMP
#pragma omp parallel for schedule(static) reduction(+ : energy)
#endif
  for (int z = 1; z < _zRes - 1; z++)
    for (int y = 1; y < _yRes - 1; y++)
      for (int x = 1; x < _xRes - 1; x++) {
        int index = z * _slabSize + y * _xRes + x;
        int idx = _gridIdxToSolutionIdx[index];
        int solidIdx = _gridIdxToSolidIdx[index];

        Real fluidMass = max(1.0 - _nodeSolidVolumeFraction[index], 0.0) *
                         _fluidDensity * _dh * _dh * _dh;
        Real mass = fluidMass + _solidMass[index];

        Real h = _dh * y + _origin[1];
        energy += mass * gravity * h;

        if (solidIdx < 0) {
          // gravity
          _unitForce.segment<3>(idx * 3) += mass * gravity * VEC3F(0, -1, 0);

          // buoyancy
          _unitForce[idx * 3 + 1] += bouyancy * fluidMass * _intensity[index];

          _nzMassVec[idx * 3] = _nzMassVec[idx * 3 + 1] =
              _nzMassVec[idx * 3 + 2] = mass;

          _mv.segment<3>(idx * 3) = fluidMass * _velocity[index];
          // add vorticity confinement force
          if (_vorticityEps > 0) {
            Real N[3];
            if (_nodeSolidVolumeFraction[index + _xRes] == 0 &&
                _nodeSolidVolumeFraction[index - _xRes] == 0 &&
                _nodeSolidVolumeFraction[index + 1] == 0 &&
                _nodeSolidVolumeFraction[index - 1] == 0 &&
                _nodeSolidVolumeFraction[index + _slabSize] == 0 &&
                _nodeSolidVolumeFraction[index - _slabSize] == 0) {
              int up = y == _yRes - 2 ? index : index + _xRes;
              int down = _nodeSolidVolumeFraction[index - _xRes] >= 1 || y == 1
                             ? index
                             : index - _xRes;
              Real dy = (up == index || down == index) ? _dhInv : gridSize;

              int out = _nodeSolidVolumeFraction[index + _slabSize] >= 1 ||
                                z == _zRes - 2
                            ? index
                            : index + _slabSize;
              int in =
                  _nodeSolidVolumeFraction[index - _slabSize] >= 1 || z == 1
                      ? index
                      : index - _slabSize;

              Real dz = (out == index || in == index) ? _dhInv : gridSize;
              int right =
                  _nodeSolidVolumeFraction[index + 1] >= 1 || x == _xRes - 2
                      ? index
                      : index + 1;
              int left = _nodeSolidVolumeFraction[index - 1] >= 1 || x == 1
                             ? index
                             : index - 1;
              Real dx = (right == index || left == index) ? _dhInv : gridSize;

              N[0] = (_workspace[right] - _workspace[left]) * dx;
              N[1] = (_workspace[up] - _workspace[down]) * dy;
              N[2] = (_workspace[out] - _workspace[in]) * dz;

              Real magnitude = sqrt(N[0] * N[0] + N[1] * N[1] + N[2] * N[2]);
              if (magnitude > 0.0f) {
                magnitude = 1.0 / magnitude;
                N[0] *= magnitude;
                N[1] *= magnitude;
                N[2] *= magnitude;

                _unitForce[idx * 3] += (N[1] * _vorticity[index][2] -
                                        N[2] * _vorticity[index][1]) *
                                       _dh * _vorticityEps * fluidMass;
                _unitForce[idx * 3 + 1] -= (N[0] * _vorticity[index][2] -
                                            N[2] * _vorticity[index][0]) *
                                           _dh * _vorticityEps * fluidMass;
                _unitForce[idx * 3 + 2] += (N[0] * _vorticity[index][1] -
                                            N[1] * _vorticity[index][0]) *
                                           _dh * _vorticityEps * fluidMass;
              }
            }
          }
        } else {
          _dampedSolidMassVec[solidIdx * 3] =
              _dampedSolidMassVec[solidIdx * 3 + 1] =
                  _dampedSolidMassVec[solidIdx * 3 + 2] =
                      _solidMass[index] * solidMassDampCoeff + fluidMass;
          _solidMv.segment<3>(solidIdx * 3) += fluidMass * _velocity[index];
          _solidUnitForce.segment<3>(solidIdx * 3) +=
              mass * gravity * VEC3F(0, -1, 0);
          // buoyancy
          _solidUnitForce[solidIdx * 3 + 1] +=
              bouyancy * fluidMass * _intensity[index];
        }
        _nzMassVecInv[idx * 3] = _nzMassVecInv[idx * 3 + 1] =
            _nzMassVecInv[idx * 3 + 2] = 1.0 / mass;

        energy += 0.5 * fluidMass * _velocity[index].dot(_velocity[index]);
      }
  TIMING_BREAKDOWN::toc("Compute Body forces");

  TIMING_BREAKDOWN::tic();
  _systemMatrix.resize(_totalSolidCells * 3, _totalSolidCells * 3);
  _systemMatrix.setFromTriplets(tripletList.begin(), tripletList.end());

  _systemMatrix *= _dt * _dt + rayleighBeta * _dt;

  for (unsigned int i = 0; i < _dampedSolidMassVec.size(); i++) {
    _systemMatrix.coeffRef(i, i) += _dampedSolidMassVec[i];
  }

  _systemMatrix.makeCompressed();
  TIMING_BREAKDOWN::toc("construct global system matrix");

#if HYBRID_DIV_CONSTRAINT
  TIMING_BREAKDOWN::tic();
  computeDivergenceConstraintJacobian();
  TIMING_BREAKDOWN::toc("computeDivergenceConstraintJacobian");
#else
  TIMING_BREAKDOWN::tic();
  computeDivergenceConstraintRHS();
  TIMING_BREAKDOWN::toc("computeDivergenceConstraintRHS");
#endif

  pcrImplicitSolve();
  finalizeSolution();

  return energy;
}

void COUPLED_INTEGRATOR_3D::getMatVecMult(const VECTOR& input, VECTOR& output) {
  static VECTOR midResult(_JT.rows());
#if USING_OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (int k = 0; k < _JT.rows(); ++k) {
    midResult[k] = 0;
    for (SpMat::InnerIterator it(_JT, k); it; ++it) {
      midResult[k] += it.value() * input[it.col()];
    }
  }

  output.conservativeResize(_divergenceConstraintJacobian.rows());
#if USING_OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (int k = 0; k < _divergenceConstraintJacobian.rows(); ++k) {
    output[k] = 0;
    for (SpMat::InnerIterator it(_divergenceConstraintJacobian, k); it; ++it) {
      output[k] += it.value() * midResult[it.col()] * _nzMassVecInv[it.col()];
    }
  }
}

void COUPLED_INTEGRATOR_3D::pcrImplicitSolve() {
  _systemMatrixDiag = _systemMatrix.diagonal();
  _solidUnitForce *= _dt;
  _solidUnitForce += _solidMv;

  // cout << "solid force " << _solidUnitForce.norm() << endl;

  _unitForce *= _dt;
  _unitForce += _mv;

  cout << "Implicit update" << endl;
  TIMING_BREAKDOWN::tic();
  // Eqn (8) in the paper
  _solver.solvePCR(_systemMatrix, _solidUnitForce, _solidVelocityAtDt0);
  TIMING_BREAKDOWN::toc("Solid velocity solve 1");

  _velocityAtDt0.conservativeResize(_totalActiveCells * 3);
#if USING_OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (int index = 0; index < _totalCells; index++) {
    int idx = _gridIdxToSolutionIdx[index];
    if (idx < 0) continue;
    int solidIdx = _gridIdxToSolidIdx[index];
    if (solidIdx < 0)
      _velocityAtDt0.segment<3>(idx * 3) =
          _unitForce.segment<3>(idx * 3) * _nzMassVecInv[idx * 3];
    else
      _velocityAtDt0.segment<3>(idx * 3) =
          _solidVelocityAtDt0.segment<3>(solidIdx * 3);
  }

  _JMInvJTDiag.conservativeResize(_divergenceConstraintJacobian.rows());
  TIMING_BREAKDOWN::tic();
#if USING_OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (int i = 0; i < _JJTDiagEntries.size(); i++) {
    Real sum = 0;
    for (pair<int, Real>& entry : _JJTDiagEntries[i]) {
      sum += _nzMassVecInv[entry.first] * entry.second;
    }
    _JMInvJTDiag[i] = sum;
  }
  TIMING_BREAKDOWN::toc("Compute JMInvJT diagonal");

  VECTOR b1 = _divergenceConstraintJacobian * _velocityAtDt0 -
              _divergenceConstraintsRHS;

  TIMING_BREAKDOWN::tic();
  static VECTOR lambda1;
  cout << "Solve for pressure" << endl;
  // Eqn (12) in the paper
  _divFreeSolver.solvePCR(b1, lambda1,
                          &COUPLED_INTEGRATOR_3D::staticGetMatVecMult, this);
  TIMING_BREAKDOWN::toc("Pressure Solve");

  b1 = -_divergenceConstraintJacobian.transpose() * lambda1;
  VECTOR solidB1(_totalSolidCells * 3);
#if USING_OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (int x = 0; x < _totalCells; x++) {
    int idx = _gridIdxToSolutionIdx[x];
    int solidIdx = _gridIdxToSolidIdx[x];
    if (solidIdx >= 0)
      solidB1.segment<3>(solidIdx * 3) = b1.segment<3>(idx * 3);
  }

  cout << "Remove divergent part" << endl;
  TIMING_BREAKDOWN::tic();
  // Eqn (13) in the paper
  _solver.solvePCR(_systemMatrix, solidB1, _solidPressureCorrection);
  TIMING_BREAKDOWN::toc("Solid velocity solve 2");

  _pressureCorrection.conservativeResize(_totalActiveCells * 3);
#if USING_OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (int index = 0; index < _totalCells; index++) {
    int idx = _gridIdxToSolutionIdx[index];
    if (idx < 0) continue;
    int solidIdx = _gridIdxToSolidIdx[index];
    if (solidIdx < 0)
      _pressureCorrection.segment<3>(idx * 3) =
          b1.segment<3>(idx * 3) * _nzMassVecInv[idx * 3];
    else
      _pressureCorrection.segment<3>(idx * 3) =
          _solidPressureCorrection.segment<3>(solidIdx * 3);
  }

  _velocityAtDt1 = _velocityAtDt0 + _pressureCorrection;

  _velocity.clear();

#if USING_OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (unsigned int i = 0; i < _totalCells; i++) {
    int idx = _gridIdxToSolutionIdx[i];
    if (idx >= 0) _velocity[i] = _velocityAtDt1.segment<3>(idx * 3);
  }
}
// Semi-Lagrangian advection
void COUPLED_INTEGRATOR_3D::advectField(FIELD_3Df& field, FIELD_3Df& workspace,
                                        Real dt) {
  workspace.clear();
#if USING_OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (int z = 0; z < _zRes - 0; z++)
    for (int y = 0; y < _yRes - 0; y++)
      for (int x = 0; x < _xRes - 0; x++) {
        int index = z * _slabSize + y * _xRes + x;
        if (_solidMass[index] > 0) continue;

        VEC3F node(x * _dh + _origin[0], y * _dh + _origin[1],
                   z * _dh + _origin[2]);
        // VEC3F v = _velocity.interpolate_value(node - 0.5 * _velocity[index] *
        // dt);
        VEC3F v = _velocity[index];
        VEC3F newPos = node - v * dt;

        Real mass = _solidMass.interpolate_value(newPos);
        Real newdt = _dt / 2;
        while (mass > 0 && newdt > dt / 16) {
          newPos += v * newdt;
          newdt /= 2;
          mass = _solidMass.interpolate_value(newPos);
        }

        workspace[index] = field.interpolate_value(newPos);
      }
  field.swap(workspace);
}
void COUPLED_INTEGRATOR_3D::advectFluidVelocityField(VEC3F_FIELD_3D& workspace,
                                                     Real dt) {
  TIMING_BREAKDOWN::tic();
#if USING_OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (int z = 1; z < _zRes - 1; z++)
    for (int y = 1; y < _yRes - 1; y++)
      for (int x = 1; x < _xRes - 1; x++) {
        int index = z * _slabSize + y * _xRes + x;
        if (_solidMass[index] > 0) continue;

        VEC3F node(x * _dh + _origin[0], y * _dh + _origin[1],
                   z * _dh + _origin[2]);
        VEC3F v =
            _velocity.interpolate_value(node - 0.5 * _velocity[index] * dt);
        VEC3F newPos = node - v * dt;

        Real mass = _solidMass.interpolate_value(newPos);
        Real newdt = dt / 2;
        while (mass > 0 && newdt > dt / 16) {
          newPos += v * newdt;
          newdt /= 2;
          mass = _solidMass.interpolate_value(newPos);
        }
        workspace[index] = _velocity.interpolate_value(newPos);
      }
  _velocity.swap(workspace);
  TIMING_BREAKDOWN::toc("advectFluidVelocityField");
}

void COUPLED_INTEGRATOR_3D::finalizeSolution() {
  cout << "Copy global velocity to solids" << endl;
  for (int i = 0; i < _solids.size(); i++)
    _solids[i]->finalizeSolution(_velocity, _dt);

  cout << "FLIP::grid to particle" << endl;
  for (int i = 0; i < _solids.size(); i++) _solids[i]->gridToParticle();

  cout << "FLIP::particle to grid" << endl;
  TIMING_BREAKDOWN::tic();
  for (int i = 0; i < _solids.size(); i++) {
    _solids[i]->particleToGrid();
  }
  TIMING_BREAKDOWN::toc("Particle to Grid");

  TIMING_BREAKDOWN::tic();
#if USING_OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (int i = 0; i < _solids.size(); i++) {
    _solids[i]->extendVectorField(_solids[i]->_velocity);
  }
  TIMING_BREAKDOWN::toc("ExtendVectorField");

  cout << "Update Solid Displacement" << endl;
  TIMING_BREAKDOWN::tic();
#if USING_OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (int i = 0; i < _solids.size(); i++) _solids[i]->finalizeSolution();
  TIMING_BREAKDOWN::toc("Update Solid Displacement");

  cout << "Compute Solid Masses" << endl;
  TIMING_BREAKDOWN::tic();
  for (int i = 0; i < _solids.size(); i++) _solids[i]->computeSolidMass(true);
  TIMING_BREAKDOWN::toc("Compute Solid Masses");

  cout << "Advect surface meshes" << endl;
  TIMING_BREAKDOWN::tic();
  for (int i = 0; i < _solids.size(); i++) _solids[i]->advectSurfaceMesh();
  TIMING_BREAKDOWN::toc("Advect surface meshes");

  cout << "updater solid sdf & init particles" << endl;
  TIMING_BREAKDOWN::tic();
#if USING_OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (int i = 0; i < _solids.size(); i++) {
    _solids[i]->computeSDF();
    _solids[i]->initMaterialParticles();
  }
  TIMING_BREAKDOWN::toc("solid sdf & init particles");

  // cout << "gatherSolidMasses" << endl;
  gatherSolidMasses();

  cout << "Advect fluid velocity field" << endl;
  TIMING_BREAKDOWN::tic();
  advectFluidVelocityField(_vecWorkspace, _dt);
  TIMING_BREAKDOWN::toc("Advect fluid velocity field");

  cout << "Advect fluid intensity field" << endl;
  TIMING_BREAKDOWN::tic();
  advectField(_intensity, _workspace, _dt);
  TIMING_BREAKDOWN::toc("Advect fluid particles");

#if USING_OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (int i = 0; i < _solids.size(); i++)
    _solids[i]->recoverVelocityField(_velocity);

  // Set Dirichlet Boundary condition
  for (int x = 0; x < _xRes; x++) {
    for (int y = 0; y < _yRes; y++) {
      _velocity(x, y, 0).setZero();
      _velocity(x, y, _zRes - 1).setZero();
    }
    for (int z = 0; z < _zRes; z++) {
      _velocity(x, 0, z).setZero();
      _intensity(x, 0, z) = 0;
      _velocity(x, _yRes - 1, z).setZero();
    }
  }
  for (int y = 0; y < _yRes; y++)
    for (int z = 0; z < _zRes; z++) {
      _velocity(0, y, z).setZero();
      _velocity(_xRes - 1, y, z).setZero();
    }
}

void COUPLED_INTEGRATOR_3D::drawSolidMasses() {
  for (DEFOR_SOLID_3D* solid : _solids) solid->drawMass();
  // _solids[2]->drawMass();
}
void COUPLED_INTEGRATOR_3D::drawVelocity() {
  glBegin(GL_LINES);
  for (int z = 1; z < _zRes - 1; z++)
    for (int y = 1; y < _yRes - 1; y++)
      for (int x = 1; x < _xRes - 1; x++) {
        if (_solidMass(x, y, z) == 0)
          glColor3f(0, 1, 0);
        else
          glColor3f(0, 0, 1);

        VEC3F node(x * _dh + _origin[0], y * _dh + _origin[1],
                   z * _dh + _origin[2]);
        VEC3F end = node + _velocity(x, y, z) / _velocity(x, y, z).norm() * _dh;

        glVertex3f(node[0], node[1], node[2]);
        glVertex3f(end[0], end[1], end[2]);
      }
  glEnd();
}

void COUPLED_INTEGRATOR_3D::drawSolidVelocities() {
  for (DEFOR_SOLID_3D* solid : _solids) solid->drawVelocity();
}
void COUPLED_INTEGRATOR_3D::drawSolidParticles() {
  for (DEFOR_SOLID_3D* solid : _solids) solid->drawMaterialParticles();
}
void COUPLED_INTEGRATOR_3D::drawSolidSurfaces() {
  glColor4f(1, 0, 0, 1);
  for (DEFOR_SOLID_3D* solid : _solids) solid->drawSurfaceMesh();
}

void COUPLED_INTEGRATOR_3D::drawIntensity() {
  Real strideDx = _dh * 0.5;
  Real strideDy = _dh * 0.5;
  Real strideDz = _dh * 0.5;

  glPointSize(10.f);
  glBegin(GL_POINTS);
  int index = 0;
  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++, index++) {
        VEC3F center = VEC3F(x + 0.5, y + 0.5, z + 0.5) * _dh + _origin;
        if (fabs(_intensity[index]) < 1e-3) {
          continue;
        }
        Real c = fabs(_intensity[index]);
        glColor4f(1, 1, 1, c);
        glColor4f(1, 1, 1, c);
        glVertex3f(center[0], center[1], center[2]);
      }
  glEnd();
}

void COUPLED_INTEGRATOR_3D::readFrame(int frame, bool debug) {
#if USING_OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (int i = 0; i < _solids.size(); i++) {
    string filename = SIMPLE_PARSER::getString("data path", "") + "frame." +
                      IO::itoPaddedString(frame) + ".solid." + to_string(i);
    _solids[i]->readFrame(filename, frame, debug);
  }
  if (debug || (frame % 10 == 0)) {
    _velocity.readGz(SIMPLE_PARSER::getString("data path", "") + "frame." +
                     IO::itoPaddedString(frame) + ".velocity");
  }
  _intensity.readGz(SIMPLE_PARSER::getString("data path", "") + "frame." +
                    IO::itoPaddedString(frame) + ".intensity");
  gatherSolidMasses();
  for (DEFOR_SOLID_3D* solid : _solids) {
    solid->recoverVelocityField(_velocity);
  }
}
void COUPLED_INTEGRATOR_3D::writeFrame(int frame, bool debug) {
#if USING_OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (int i = 0; i < _solids.size(); i++) {
    string filename = SIMPLE_PARSER::getString("data path", "") + "frame." +
                      IO::itoPaddedString(frame) + ".solid." + to_string(i);
    _solids[i]->writeFrame(filename, frame, debug);
  }
  if (debug || (frame % 10 == 0)) {
    _velocity.writeGz(SIMPLE_PARSER::getString("data path", "") + "frame." +
                      IO::itoPaddedString(frame) + ".velocity");
  }
  _intensity.writeGz(SIMPLE_PARSER::getString("data path", "") + "frame." +
                     IO::itoPaddedString(frame) + ".intensity");
}