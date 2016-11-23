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
#ifndef DEFOR_SOLID_3D_H
#define DEFOR_SOLID_3D_H

#include <SETTINGS.h>
#include <deformCD/aabb.h>
#include <dtgrid/SPARSE_SDF.h>
#include <geometry/FIELD_3D.h>
#include <geometry/MATRIX3_FIELD_3D.h>
#include <geometry/OBJ.h>
#include <geometry/TRIANGLE_BVH.h>
#include <geometry/VEC3F_FIELD_3D.h>
#include <material/COROTATION.h>
#include <set>

struct PARTICLE {
  VEC3F velocity;
  VEC3F restPos;
  VEC3F pos;
  Real mass;
  bool pinned;
};

class DEFOR_SOLID_3D {
  friend class COUPLED_INTEGRATOR_3D;

 public:
  DEFOR_SOLID_3D(int xRes, int yRes, int zRes, Real dh, VEC3F origin,
                 const FIELD_3Df& restMass);
  ~DEFOR_SOLID_3D();
  int xRes() const { return _xRes; };
  int yRes() const { return _yRes; };
  int zRes() const { return _zRes; };
  Real dh() const { return _dh; };
  const VEC3F& origin() const { return _origin; };

  void setInitialState(MATRIX3 rotation, VEC3F translation, VEC3F scaling);
  void setInitialVelocity(VEC3F v);
  void setMaterial(Real lambda, Real mu);
  void setPlasticYield(Real v) {
    _isPlastic = true;
    _plasticYield = v;
  }

  void setId(int id) { _id = id; };

  void loadRestSDF(const string& filename);
  OBJ* surfaceMesh() { return _surfaceMesh; };
  void setMassDensity(Real density) {
    _quarterCellSolidMass = density * _dh * _dh * _dh * 0.125;
  }
  void loadSurfaceMesh(const string& filename);

  FIELD_3Df& phi() { return _phi; };

  void init();
  void finalizeSolution();
  void finalizeSolution(VEC3F_FIELD_3D& velocity, Real dt);

  void recoverVelocityField(VEC3F_FIELD_3D& velocity);

  void drawMaterialParticles();
  void drawMass();
  void drawVelocity();
  void drawSurfaceMesh() {
    if (_surfaceMesh) _surfaceMesh->draw();
  };

  bool overlap(DEFOR_SOLID_3D* other) { return _bvh->collide(other->_bvh); }
  bool closeToWall() { return !_bvh->inside(_safeFromWallBox); }

  void addNeighbor(DEFOR_SOLID_3D* neighbor) {
    _neighbors.push_back(neighbor);
  };

  void readFrame(const string& filename, int frame, bool debug);
  void writeFrame(const string& filename, int frame, bool debug);

 protected:
  Real computeF(const MATRIX& coord, const MATRIX& pNpF, MATRIX3& F);
  void computePFPxhat(const MATRIX3& F, const MATRIX& pNpx, MATRIX& pFpxhat);

  void computeSolidMass(bool usesdf = false);
  Real computeStiffnessMatrixAndMaterialForce(VECTOR& force,
                                              vector<TRIPLET>& tripletList,
                                              const vector<int>& indexMap);

  void initMaterialParticles();
  void gridToParticle();
  void particleToGrid();

  void computeSDF();
  void sweepPhi();
  void extendVectorField(VEC3F_FIELD_3D& field);
  void sweepVectorField(VEC3F_FIELD_3D& field);
  void sweepOnePass(VEC3F_FIELD_3D& field, int i0, int i1, int j0, int j1,
                    int k0, int k1);

  void advectSurfaceMesh();
  void advectField(VEC3F_FIELD_3D& field, VEC3F_FIELD_3D& workspace, Real dt);
  Real computeMV(VECTOR& mv, const vector<int>& indexMap);

 protected:
  VEC3F_FIELD_3D _X;
  VEC3F_FIELD_3D _u;
  VEC3F_FIELD_3D _vecWorkspace;

  aabb _safeFromWallBox;
  OBJ* _surfaceMesh;
  TRIANGLE_BVH* _bvh;

  COROTATION* _solidMaterial;
  vector<COROTATION*> _materialCopies;

  const FIELD_3Df _restMass;
  SPARSE_SDF* _restSDF;

  vector<DEFOR_SOLID_3D*> _neighbors;

  // Plastic deformation gradient
  MATRIX3_FIELD_3D _FpInv;

  int _xmin, _ymin, _zmin, _xmax, _ymax, _zmax;
  VEC3I _tightMin;
  VEC3I _tightMax;

  vector<int> _activeCells;
  vector<int> _gridIndexToActiveCellIndex;

  SpMat _systemMatrix;

  int _xRes;
  int _yRes;
  int _zRes;
  int _slabSize;
  int _totalCells;

  VEC3F _origin;
  Real _dh;
  Real _dhInv;
  Real _subsampleDh;

  Real _defaultDt;
  Real _dt;
  Real _quarterCellSolidMass;

  FIELD_3Dc _valid;
  FIELD_3Df _phi;
  FIELD_3Df _mass;
  VECTOR _quarterVolumeFractions;

  FIELD_3Df _workspace;

  VEC3F_FIELD_3D _velocity;
  VECTOR _unitForce;
  VECTOR _mv;

  bool _isPlastic;
  Real _plasticYield;

  bool _solveVelocityAtDt0;
  VECTOR _velocityAtDt0;
  VECTOR _velocityAtDt1;
  VECTOR _gradient;
  VECTOR _constraintObjective;
  Real _objective;

  FIELD_3Df _detFInv;

  MATRIX _G;
  vector<VEC3F> _quadratureParams;
  vector<MATRIX> _quadraturePNpx;
  vector<VECTOR> _quadratureBaryCentric;
  vector<PARTICLE> _materialParticles;

  set<int> _moniters;

  int _id;

  VEC3F _initialTranslation, _initalScaling;
  MATRIX3 _initialRotation;
};

#endif