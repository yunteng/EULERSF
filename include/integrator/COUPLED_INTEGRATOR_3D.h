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
#ifndef COUPLED_INTEGRATOR_3D_H
#define COUPLED_INTEGRATOR_3D_H

#include <SETTINGS.h>
#include <integrator/DEFOR_SOLID_3D.h>
#include <linearalgebra/CONJUGATE_RESIDUAL.h>
#include <linearalgebra/JACOBI_PRECONDITIONER.h>
#include <util/TIMING_BREAKDOWN.h>
#if USING_OPENMP
#include <omp.h>
#endif

class COUPLED_INTEGRATOR_3D {
 public:
  COUPLED_INTEGRATOR_3D(int xRes, int yRes, int zRes, Real dh, VEC3F origin);
  ~COUPLED_INTEGRATOR_3D();
  int xRes() const { return _xRes; };
  int yRes() const { return _yRes; };
  int zRes() const { return _zRes; };
  Real dh() const { return _dh; };
  const VEC3F& origin() const { return _origin; };

  FIELD_3Df& intensity() { return _intensity; };
  FIELD_3Df& solidMass() { return _solidMass; };

  void addSolid(DEFOR_SOLID_3D* solid) { _solids.push_back(solid); };
  void setFluidDensity(Real d) { _fluidDensity = d; };
  void setFluidViscocity(Real v) { _fluidViscocity = v; };
  void setMaxDt(Real dt) { _defaultDt = dt; };
  Real dt() { return _dt; };
  // Compute the time step size given a CFL condition
  Real CFL(Real alpha);

  void drawSolidMasses();
  void drawSolidVelocities();
  void drawVelocity();
  void drawSolidParticles();
  void drawSolidSurfaces();
  void drawIntensity();

  void init();
  Real stepImplicit();

  void pcrImplicitSolve();
  void finalizeSolution();

  void computeDivergenceConstraintJacobian();
  void computeDivergenceConstraintRHS();

  void gatherSolidMasses();

  void getMatVecMult(const VECTOR& input, VECTOR& output);
  static void staticGetMatVecMult(void* integrator, const VECTOR& input,
                                  VECTOR& output) {
    static_cast<COUPLED_INTEGRATOR_3D*>(integrator)
        ->getMatVecMult(input, output);
  }

  void writeFrame(int frame, bool debug = false);
  void readFrame(int frame, bool debug = false);

  void reAdvectIntensity() { advectField(_intensity, _workspace, _dt); }

 private:
  void advectFluidVelocityField(VEC3F_FIELD_3D& workspace, Real dt);
  void advectField(FIELD_3Df& field, FIELD_3Df& workspace, Real dt);

 private:
  int _xRes;
  int _yRes;
  int _zRes;
  int _slabSize;
  int _totalCells;
  int _totalActiveCells;
  int _totalSolidCells;
  VEC3F _origin;
  Real _dh;
  Real _dhInv;

  vector<DEFOR_SOLID_3D*> _solids;

  VEC3F_FIELD_3D _velocity;
  VEC3F_FIELD_3D _vecWorkspace;
  VEC3F_FIELD_3D _vorticity;

  Real _defaultDt;
  Real _dt;

  VECTOR _solidUnitForce;
  VECTOR _solidMv;

  VECTOR _unitForce;
  VECTOR _mv;
  VECTOR _nzMassVec;
  VECTOR _dampedSolidMassVec;
  VECTOR _nzMassVecInv;

  FIELD_3Df _solidMass;
  FIELD_3Df _nodeSolidVolumeFraction;
  VECTOR _quarterSolidVolumeFractions;

  vector<int> _gridIdxToSolutionIdx;
  vector<int> _gridIdxToSolidIdx;

  VECTOR _velocityAtDt0;
  VECTOR _velocityAtDt1;
  VECTOR _pressureCorrection;

  VECTOR _solidVelocityAtDt0;
  VECTOR _solidVelocityAtDt1;
  VECTOR _solidPressureCorrection;

  SpMat _systemMatrix;
  VECTOR _systemMatrixDiag;
  SpMat _JT;
  SpMat _JMInvJT;
  VECTOR _JMInvJTDiag;
  JACOBI_PRECONDITIONER _preconditioner;
  CONJUGATE_RESIDUAL<JACOBI_PRECONDITIONER> _solver;
  JACOBI_PRECONDITIONER _divfreeSolverPreconditioner;
  CONJUGATE_RESIDUAL<JACOBI_PRECONDITIONER> _divFreeSolver;

  vector<vector<pair<int, Real> > > _JJTDiagEntries;

  Real _fluidDensity;
  Real _fluidViscocity;
  SpMat _divergenceConstraintJacobian;
  VECTOR _divergenceConstraintsRHS;

  FIELD_3Df _intensity;
  FIELD_3Df _workspace;

  Real _vorticityEps;
  vector<bool> _isConstrained;

  vector<VEC3F> _fluidParticles;
};

#endif