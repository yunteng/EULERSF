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
#include <integrator/DEFOR_SOLID_3D.h>
#include <unordered_map>
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

#include <linearalgebra/CONJUGATE_RESIDUAL.h>
#include <linearalgebra/COO_MATRIX.h>
#include <linearalgebra/JACOBI_PRECONDITIONER.h>

#include <util/MERSENNETWISTER.h>
#include <util/SIMPLE_PARSER.h>
#include <util/TIMING_BREAKDOWN.h>

#define mymap unordered_map

DEFOR_SOLID_3D::DEFOR_SOLID_3D(int xRes, int yRes, int zRes, Real dh,
                               VEC3F origin, const FIELD_3Df& restMass)
    : _xRes(xRes),
      _yRes(yRes),
      _zRes(zRes),
      _dh(dh),
      _origin(origin),
      _restMass(restMass),
      _restSDF(NULL),
      _quarterCellSolidMass(1),
      _surfaceMesh(NULL),
      _bvh(NULL),
      _defaultDt(SIMPLE_PARSER::getFloat("dt", 0.01)),
      _solveVelocityAtDt0(false),
      _solidMaterial(NULL),
      _isPlastic(false) {
  _slabSize = _xRes * _yRes;
  _totalCells = _xRes * _yRes * _zRes;
  _dhInv = 1.0 / _dh;
  _subsampleDh = SIMPLE_PARSER::getFloat("mass sample dh", _dh / 10);

  _safeFromWallBox._min = _origin + VEC3F(2, 2, 2) * _dh;
  _safeFromWallBox._max =
      _origin + VEC3F(_xRes - 3, _yRes - 3, _zRes - 3) * _dh;

  _X.resize(_xRes, _yRes, _zRes, _dh, _origin);
  _velocity.resizeLike(_X);
  _u.resizeLike(_X);
  _vecWorkspace.resizeLike(_X);

  _mass.resize(_xRes, _yRes, _zRes, _dh, _origin);
  _workspace.resizeLike(_mass);
  _phi.resizeLike(_mass);

  _quarterVolumeFractions.resize(_totalCells * 8);
  _quarterVolumeFractions.setZero();

  _valid.resize(_xRes - 1, _yRes - 1, _zRes - 1, _dh, _origin);
  _detFInv.resize(_xRes, _yRes, _zRes, _dh, _origin);

  _dt = _defaultDt;
}

DEFOR_SOLID_3D::~DEFOR_SOLID_3D() {
  if (_bvh != NULL) {
    delete _bvh;
    _bvh = NULL;
  }
  if (_surfaceMesh != NULL) {
    delete _surfaceMesh;
    _surfaceMesh = NULL;
  }
  if (_solidMaterial != NULL) {
    delete _solidMaterial;
    _solidMaterial = NULL;
  }

  for (int i = 0; i < _materialCopies.size(); i++) {
    if (_materialCopies[i]) {
      delete _materialCopies[i];
      _materialCopies[i] = NULL;
    }
  }

  if (_restSDF != NULL) {
    delete _restSDF;
    _restSDF = NULL;
  }
}

void DEFOR_SOLID_3D::setMaterial(Real lambda, Real mu) {
  if (_solidMaterial) delete _solidMaterial;
  _solidMaterial = new COROTATION(lambda, mu);
#if USING_OPENMP
  int totalCores = omp_get_max_threads();
#else
  int totalCores = 1;
#endif
  _materialCopies.resize(totalCores);
  for (int i = 0; i < totalCores; i++)
    _materialCopies[i] = new COROTATION(lambda, mu);
}

void DEFOR_SOLID_3D::setInitialState(MATRIX3 rotation, VEC3F translation,
                                     VEC3F scaling) {
  _initialRotation = rotation;
  _initialTranslation = translation;
  _initalScaling = scaling;
  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++) {
        VEC3F point(x, y, z);
        point *= _dh;
        point += _origin;
        point -= translation;
        point = rotation.transpose() * point;
        point[0] /= scaling[0];
        point[1] /= scaling[1];
        point[2] /= scaling[2];
        _X(x, y, z) = point;
      }
  int index = 0;
  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++, index++) {
        VEC3F pos(x * _dh, y * _dh, z * _dh);
        pos += _origin;
        _u[index] = pos - _X[index];
      }
}
void DEFOR_SOLID_3D::setInitialVelocity(VEC3F v) {
  for (int i = 0; i < _totalCells; i++) _velocity[i] = v;
}

void DEFOR_SOLID_3D::loadSurfaceMesh(const string& filename) {
  _surfaceMesh = new OBJ();
  _surfaceMesh->Load(filename);

  _bvh = new TRIANGLE_BVH(_surfaceMesh);

  vector<VEC3F>& vertices = _surfaceMesh->restPose;
  VEC3F boxMin(1e9, 1e9, 1e9);
  VEC3F boxMax(-1e9, -1e9, -1e9);
  for (int x = 0; x < vertices.size(); x++) {
    boxMin[0] = min(boxMin[0], vertices[x][0]);
    boxMin[1] = min(boxMin[1], vertices[x][1]);
    boxMin[2] = min(boxMin[2], vertices[x][2]);

    boxMax[0] = max(boxMax[0], vertices[x][0]);
    boxMax[1] = max(boxMax[1], vertices[x][1]);
    boxMax[2] = max(boxMax[2], vertices[x][2]);
  }

  boxMin -= VEC3F(_dh, _dh, _dh);
  boxMax += VEC3F(_dh, _dh, _dh);

  // initialize the inverse plastic deformation gradient field
  if (_isPlastic) {
    Real materialDh = _dh * .5;
    int x = (boxMax[0] - boxMin[0]) / materialDh + 1;
    int y = (boxMax[1] - boxMin[1]) / materialDh + 1;
    int z = (boxMax[2] - boxMin[2]) / materialDh + 1;

    _FpInv.resize(x, y, z, materialDh, boxMin);
  }

  for (int i = 0; i < 3; i++) {
    boxMin[i] *= _initalScaling[i];
    boxMax[i] *= _initalScaling[i];
  }
  boxMin = _initialRotation * boxMin + _initialTranslation;
  boxMax = _initialRotation * boxMax + _initialTranslation;

  boxMin -= _origin;
  boxMin *= _dhInv;

  boxMax -= _origin;
  boxMax *= _dhInv;

  VEC3F res(_xRes, _yRes, _zRes);
  for (int i = 0; i < 3; i++) {
    _tightMin[i] = max(boxMin[i], 0.0);
    _tightMax[i] = min(ceil(boxMax[i]), res[i]);
  }
}
// Load narrow-banded sdf generated by DT-Grid
// (https://code.google.com/archive/p/dt-grid/)
void DEFOR_SOLID_3D::loadRestSDF(const string& filename) {
  _restSDF = new SPARSE_SDF();
  _restSDF->load(filename);
  _restSDF->initHashSurface(_surfaceMesh->triangles);
}

void DEFOR_SOLID_3D::computeSDF() {
  int index = 0;
  Real large_distance = _xRes + _yRes + _zRes + 2;

  _xmax = _ymax = _zmax = 0;
  _xmin = _xRes;
  _ymin = _yRes;
  _zmin = _zRes;

  _phi = large_distance;

  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++, index++) {
        if (_mass[index] > 0) {
          _phi[index] = -0.5;
          if (x < _xmin) _xmin = x;
          if (x > _xmax) _xmax = x;

          if (y < _ymin) _ymin = y;
          if (y > _ymax) _ymax = y;

          if (z < _zmin) _zmin = z;
          if (z > _zmax) _zmax = z;
        }
      }

  static int tightBorder = SIMPLE_PARSER::getInt("tight border", 3);
  static int extrap = SIMPLE_PARSER::getInt("extrapolation cells", 20);

  cout << "solid center velocity "
       << _velocity((_xmin + _xmax) / 2, (_ymin + _ymax) / 2,
                    (_zmin + _zmax) / 2)
              .transpose()
       << endl;

  _tightMin = VEC3I(max(_xmin - tightBorder, 0), max(_ymin - tightBorder, 0),
                    max(_zmin - tightBorder, 0));

  _tightMax =
      VEC3I(min(_xmax + tightBorder, _xRes), min(_ymax + tightBorder, _yRes),
            min(_zmax + tightBorder, _zRes));

  _xmin = max(_xmin - extrap, 0);
  _ymin = max(_ymin - extrap, 0);
  _zmin = max(_zmin - extrap, 0);

  _xmax = min(_xmax + extrap, _xRes);
  _ymax = min(_ymax + extrap, _yRes);
  _zmax = min(_zmax + extrap, _zRes);

  for (int i = 0; i < 2; ++i) sweepPhi();
}
void DEFOR_SOLID_3D::init() {
  // set up 8-point Newtwon quadrature for evaluating the material force
  _quadratureParams.resize(8);
  Real oneoversqrt3 = 1.0 / sqrt(3);

  for (int i = 0; i < 8; i++) {
    int signX = i % 2 == 0 ? -1 : 1;
    int signY = i % 4 <= 1 ? -1 : 1;
    int signZ = i < 4 ? -1 : 1;
    _quadratureParams[i] =
        VEC3F(signX * oneoversqrt3, signY * oneoversqrt3, signZ * oneoversqrt3);
  }

  _quadraturePNpx.clear();
  _quadratureBaryCentric.clear();
  for (int x = 0; x < 8; x++) {
    const VEC3F& param = _quadratureParams[x];
    MATRIX pNpx(8, 3);
    VECTOR n(8);
    for (int i = 0; i < 8; i++) {
      int signX = i % 2 == 0 ? -1 : 1;
      int signY = i % 4 <= 1 ? -1 : 1;
      int signZ = i < 4 ? -1 : 1;
      pNpx(i, 0) = signX * (1 + signY * param[1]) * (1 + signZ * param[2]);
      pNpx(i, 1) = signY * (1 + signX * param[0]) * (1 + signZ * param[2]);
      pNpx(i, 2) = signZ * (1 + signX * param[0]) * (1 + signY * param[1]);
      n[i] = (1 + signX * param[0]) * (1 + signY * param[1]) *
             (1 + signZ * param[2]);
    }
    pNpx *= 0.25 * _dhInv;
    n *= 0.125;
    _quadraturePNpx.push_back(pNpx);
    _quadratureBaryCentric.push_back(n);
  }

  _G.resize(8, 3);
  for (int i = 0; i < 8; i++) {
    _G(i, 0) = i % 2 == 0 ? -1 : 1;
    _G(i, 1) = i % 4 <= 1 ? -1 : 1;
    _G(i, 2) = i < 4 ? -1 : 1;
  }
  _G *= 0.25 * _dhInv;

  _xmin = 0;
  _ymin = 0;
  _zmin = 0;

  _xmax = _xRes;
  _ymax = _yRes;
  _zmax = _zRes;

  _tightMin = VEC3I(0, 0, 0);
  _tightMax = VEC3I(_xRes, _yRes, _zRes);

  cout << "compute solid mass " << endl;
  computeSolidMass();
  cout << "compute sdf" << endl;
  computeSDF();
  cout << "init material particles" << endl;
  initMaterialParticles();
  cout << "advect surface mesh" << endl;
  advectSurfaceMesh();
}

void DEFOR_SOLID_3D::initMaterialParticles() {
  _materialParticles.clear();

  static int numberOfParticlesPerCell =
      SIMPLE_PARSER::getInt("number of particles per material cell", 8);
  static Real density = _quarterCellSolidMass / (_dh * _dh * _dh * 0.125);

  static Real volume = _dh * _dh * _dh / numberOfParticlesPerCell;

  MERSENNETWISTER twister;

  for (int z = _tightMin[2]; z < _tightMax[2] - 1; z++)
    for (int y = _tightMin[1]; y < _tightMax[1] - 1; y++)
      for (int x = _tightMin[0]; x < _tightMax[0] - 1; x++) {
        if (_valid(x, y, z) != 1) continue;

        VEC3F grid(_origin[0] + x * _dh, _origin[1] + y * _dh,
                   _origin[2] + z * _dh);
        bool hasMass = false;
        for (int i = 0; i < numberOfParticlesPerCell; i++) {
          VEC3F pos = grid + VEC3F(twister.rand(_dh), twister.rand(_dh),
                                   twister.rand(_dh));
          VEC3F rPos = _X.interpolate_value(pos);
          Real d = _restMass.interpolate_value(rPos);
          PARTICLE p;
          p.pinned = false;
          p.restPos = rPos;
          p.pos = pos;
          p.velocity = _velocity.interpolate_value(pos);
          p.mass = density * volume;
          _materialParticles.push_back(p);
          hasMass = true;
        }
      }
}

static inline void solveDistance(Real p, Real q, Real r, Real& v) {
  Real distances[3];
  distances[0] = min(p, min(q, r));
  distances[2] = max(p, max(q, r));
  distances[1] = p + q + r - distances[0] - distances[2];
  Real b[3];
  Real c[3];
  b[0] = distances[0];
  c[0] = distances[0] * distances[0] - 1.0f;
  for (int x = 1; x < 3; x++) {
    b[x] = distances[x] + b[x - 1];
    c[x] = distances[x] * distances[x] + c[x - 1];
  }

  int i = 2;
  Real discrim = b[i] * b[i] - (Real)(i + 1) * c[i];
  Real newDist = (b[i] + sqrt(discrim)) / (Real)(i + 1);
  while ((discrim < 0.0f || newDist < distances[i]) && i != -1) {
    i--;
    discrim = b[i] * b[i] - (Real)(i + 1) * c[i];

    if (discrim > 0.0f) newDist = (b[i] + sqrt(discrim)) / (Real)(i + 1);
  }

  if (newDist < v) v = newDist;
}

// fast sweeping outside the solid in all four sweep directions
void DEFOR_SOLID_3D::sweepPhi() {
  int i, j, k;

  for (k = _zmin + 1; k < _zmax; k++)
    for (j = _ymin + 1; j < _ymax; j++)
      for (i = _xmin + 1; i < _xmax; i++)
        if (_mass(i, j, k) == 0)
          solveDistance(_phi(i - 1, j, k), _phi(i, j - 1, k), _phi(i, j, k - 1),
                        _phi(i, j, k));

  for (k = _zmin + 1; k < _zmax; k++)
    for (j = _ymax - 2; j >= _ymin; j--)
      for (i = _xmin + 1; i < _xmax; i++)
        if (_mass(i, j, k) == 0)
          solveDistance(_phi(i - 1, j, k), _phi(i, j + 1, k), _phi(i, j, k - 1),
                        _phi(i, j, k));

  for (k = _zmin + 1; k < _zmax; k++)
    for (j = _ymin + 1; j < _ymax; j++)
      for (i = _xmax - 2; i >= _xmin; i--)
        if (_mass(i, j, k) == 0)
          solveDistance(_phi(i + 1, j, k), _phi(i, j - 1, k), _phi(i, j, k - 1),
                        _phi(i, j, k));

  for (k = _zmin + 1; k < _zmax; k++)
    for (j = _ymax - 2; j >= _ymin; j--)
      for (i = _xmax - 2; i >= _xmin; i--)
        if (_mass(i, j, k) == 0)
          solveDistance(_phi(i + 1, j, k), _phi(i, j + 1, k), _phi(i, j, k - 1),
                        _phi(i, j, k));

  for (k = _zmax - 2; k >= _zmin; k--)
    for (j = _ymin + 1; j < _ymax; j++)
      for (i = _xmin + 1; i < _xmax; i++)
        if (_mass(i, j, k) == 0)
          solveDistance(_phi(i - 1, j, k), _phi(i, j - 1, k), _phi(i, j, k + 1),
                        _phi(i, j, k));

  for (k = _zmax - 2; k >= _zmin; k--)
    for (j = _ymax - 2; j >= _ymin; j--)
      for (i = _xmin + 1; i < _xmax; i++)
        if (_mass(i, j, k) == 0)
          solveDistance(_phi(i - 1, j, k), _phi(i, j + 1, k), _phi(i, j, k + 1),
                        _phi(i, j, k));

  for (k = _zmax - 2; k >= _zmin; k--)
    for (j = _ymin + 1; j < _ymax; j++)
      for (i = _xmax - 2; i >= _xmin; i--)
        if (_mass(i, j, k) == 0)
          solveDistance(_phi(i + 1, j, k), _phi(i, j - 1, k), _phi(i, j, k + 1),
                        _phi(i, j, k));

  for (k = _zmax - 2; k >= _zmin; k--)
    for (j = _ymax - 2; j >= _ymin; j--)
      for (i = _xmax - 2; i >= _xmin; i--)
        if (_mass(i, j, k) == 0)
          solveDistance(_phi(i + 1, j, k), _phi(i, j + 1, k), _phi(i, j, k + 1),
                        _phi(i, j, k));
}

void DEFOR_SOLID_3D::sweepOnePass(VEC3F_FIELD_3D& field, int i0, int i1, int j0,
                                  int j1, int k0, int k1) {
  int di = (i0 < i1) ? 1 : -1;
  int dj = (j0 < j1) ? 1 : -1;
  int dk = (k0 < k1) ? 1 : -1;
  Real dp, dq, dr, alpha, beta;
  Real dpn, dqn, drn;
  bool solvableX = true;
  bool solvableY = true;
  bool solvableZ = true;
  for (int k = k0; k != k1; k += dk)
    for (int j = j0; j != j1; j += dj)
      for (int i = i0; i != i1; i += di) {
        // if(!_dialatedValid(i, j))
        // continue;
        if (_mass(i, j, k) == 0) {
          dp = _phi(i, j, k) - _phi(i - di, j, k);
          dpn = _phi(i, j, k) - _phi(i + di, j, k);
          if (dp < 0) {
            if (dpn >= 0)
              continue;  // not useful on this sweep direction
            else
              solvableX = false;
          }
          dq = _phi(i, j, k) - _phi(i, j - dj, k);
          dqn = _phi(i, j, k) - _phi(i, j + dj, k);
          if (dq < 0) {
            if (dqn >= 0)
              continue;  // not useful on this sweep direction
            else
              solvableY = false;
          }
          dr = _phi(i, j, k) - _phi(i, j, k - dk);
          drn = _phi(i, j, k) - _phi(i, j, k + dk);
          if (dr < 0) {
            if (drn >= 0)
              continue;  // not useful on this sweep direction
            else
              solvableZ = false;
          }
          /*if(!solvableX && !solvableY && !solvableZ) {
            if(_moniters.find(i + j * _xRes + k * _slabSize) != _moniters.end())
              continue;
            _moniters.insert(i + j * _xRes + k * _slabSize);
          }*/
          // assert(solvableX || solvableY);
          if (!solvableX) dp = 0;
          if (!solvableY) dq = 0;
          if (!solvableZ) dr = 0;

          if (dp + dq + dr < 1e-3) {
            alpha = 0.333;
            beta = 0.333;
          } else {
            alpha = dp / (dp + dq + dr);
            beta = dq / (dp + dq + dr);
          }
          field(i, j, k) = alpha * field(i - di, j, k) +
                           beta * field(i, j - dj, k) +
                           (1 - alpha - beta) * field(i, j, k - dk);
        }
      }
}

void DEFOR_SOLID_3D::sweepVectorField(VEC3F_FIELD_3D& field) {
  int i, j, k;
  // sweep into the air
  sweepOnePass(field, _xmin + 1, _xmax - 1, _ymin + 1, _ymax - 1, _zmin + 1,
               _zmax - 1);

  sweepOnePass(field, _xmin + 1, _xmax - 1, _ymax - 2, _ymin, _zmin + 1,
               _zmax - 1);

  sweepOnePass(field, _xmax - 2, _xmin, _ymin + 1, _ymax - 1, _zmin + 1,
               _zmax - 1);

  sweepOnePass(field, _xmax - 2, _xmin, _ymax - 2, _ymin, _zmin + 1, _zmax - 1);

  sweepOnePass(field, _xmin + 1, _xmax - 1, _ymin + 1, _ymax - 1, _zmax - 2,
               _zmin);

  sweepOnePass(field, _xmin + 1, _xmax - 1, _ymax - 2, _ymin, _zmax - 2, _zmin);

  sweepOnePass(field, _xmax - 2, _xmin, _ymin + 1, _ymax - 1, _zmax - 2, _zmin);

  sweepOnePass(field, _xmax - 2, _xmin, _ymax - 2, _ymin, _zmax - 2, _zmin);

  for (i = _xmin; i < _xmax; i++)
    for (j = _ymin; j < _ymax; j++) {
      if (_mass(i, j, _zmin) == 0) field(i, j, _zmin) = field(i, j, _zmin + 1);
      if (_mass(i, j, _zmax - 1) == 0)
        field(i, j, _zmax - 1) = field(i, j, _zmax - 2);
    }

  for (i = _xmin; i < _xmax; i++)
    for (k = _zmin; k < _zmax; k++) {
      if (_mass(i, _ymin, k) == 0) field(i, _ymin, k) = field(i, _ymin + 1, k);
      if (_mass(i, _ymax - 1, k) == 0)
        field(i, _ymax - 1, k) = field(i, _ymax - 2, k);
    }

  for (j = _ymin; j < _ymax; j++)
    for (k = _zmin; k < _zmax; k++) {
      if (_mass(_xmin, j, k) == 0) field(_xmin, j, k) = field(_xmin + 1, j, k);
      if (_mass(_xmax - 1, j, k) == 0)
        field(_xmax - 1, j, k) = field(_xmax - 2, j, k);
    }
}

void DEFOR_SOLID_3D::extendVectorField(VEC3F_FIELD_3D& field) {
  for (int i = 0; i < 4; ++i) {
    // _vecWorkspace = _velocity;
    _moniters.clear();
    sweepVectorField(field);
    // _vecWorkspace -= _velocity;
    // cout << _vecWorkspace.squaredNorm() << endl;
  }
}

void DEFOR_SOLID_3D::computeSolidMass(bool usesdf) {
  _mass.clear();
  _valid.clear();
  _detFInv.clear();
  _quarterVolumeFractions.setZero();

  static Real subsampleVolume = _subsampleDh * _subsampleDh * _subsampleDh;
  static int massSampleRes = _dh / _subsampleDh;
  static int halefMassSampleRes = massSampleRes >> 1;
  static Real fraction = 1.0 / (massSampleRes * massSampleRes * massSampleRes);

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
  for (int z = _tightMin[2]; z < _tightMax[2] - 1; z++) {
    for (int y = _tightMin[1]; y < _tightMax[1] - 1; y++)
      for (int x = _tightMin[0]; x < _tightMax[0] - 1; x++) {
        int index = z * _slabSize + y * _xRes + x;
        MATRIX coord(3, 8);
        for (int i = 0; i < 8; i++) coord.col(i) = _u[index + offsets[i]];
        MATRIX3 F;
        _detFInv[index] = abs(1.0 / computeF(coord, _G, F));
      }
  }
#if USING_OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (int z = _tightMin[2]; z < _tightMax[2] - 1; z++) {
    for (int y = _tightMin[1]; y < _tightMax[1] - 1; y++)
      for (int x = _tightMin[0]; x < _tightMax[0] - 1; x++) {
        int index = z * _slabSize + y * _xRes + x;

        VEC3F node(x * _dh + _origin[0], y * _dh + _origin[1],
                   z * _dh + _origin[2]);

        if (usesdf && _phi[index] > 2) continue;

        for (int p = 0; p < massSampleRes; p++)
          for (int q = 0; q < massSampleRes; q++)
            for (int r = 0; r < massSampleRes; r++) {
              VEC3F pos = node + VEC3F((p + 0.5) * _subsampleDh,
                                       (q + 0.5) * _subsampleDh,
                                       (r + 0.5) * _subsampleDh);
              VEC3F rPos = _X.interpolate_value(pos);
              Real d = _restMass.interpolate_value(rPos) * subsampleVolume;
              if (d < 1e-6) continue;

              int offset = p < halefMassSampleRes ? 0 : 1;
              offset += q < halefMassSampleRes ? 0 : 2;
              offset += r < halefMassSampleRes ? 0 : 4;
              _quarterVolumeFractions[index * 8 + offset] += fraction;
            }
        Real sum = 0;
        for (int i = 0; i < 8; i++) {
          sum += _quarterVolumeFractions[index * 8 + i];
        }
        Real threadhold = 0.5;
        if (sum > threadhold) {
          for (int i = 0; i < 8; i++) {
            _quarterVolumeFractions[index * 8 + i] = 1;
          }
          _valid(x, y, z) = 1;
        } else {
          if (sum > 0) _valid(x, y, z) = 2;
          for (int i = 0; i < 8; i++) {
            _quarterVolumeFractions[index * 8 + i] = 0;
          }
        }
      }
  }
  for (int i = 0; i < 8; i++) {
#if USING_OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int z = _tightMin[2]; z < _tightMax[2] - 1; z++)
      for (int y = _tightMin[1]; y < _tightMax[1] - 1; y++)
        for (int x = _tightMin[0]; x < _tightMax[0] - 1; x++) {
          if (!_valid(x, y, z)) continue;
          int index = z * _slabSize + y * _xRes + x;
          _mass[index + offsets[i]] += _quarterVolumeFractions[index * 8 + i] *
                                       _detFInv[index] * _quarterCellSolidMass;
        }
  }
  cout << "mass " << _mass.data().sum() << endl;
}

Real DEFOR_SOLID_3D::computeF(const MATRIX& coord, const MATRIX& pNpF,
                              MATRIX3& F) {
  MATRIX3 Finv = MATRIX3::Identity() - coord * pNpF;
  Real Jinv = Finv.determinant();
  if (abs(Jinv) > 1e-4) {
    F = Finv.inverse();
  } else {
    F = MATRIX3::Identity();
    Jinv = 1;
  }
  return 1.0 / Jinv;
}

void DEFOR_SOLID_3D::computePFPxhat(const MATRIX3& F, const MATRIX& pNpx,
                                    MATRIX& pFpxhat) {
  MATRIX GF = pNpx * F;

  for (int i = 0; i < 8; i++) {
    for (int j = 0; j < 3; j++) {
      pFpxhat(j * 3, i * 3 + j) = GF(i, 0);
      pFpxhat(j * 3 + 1, i * 3 + j) = GF(i, 1);
      pFpxhat(j * 3 + 2, i * 3 + j) = GF(i, 2);
    }
  }
}

Real DEFOR_SOLID_3D::computeStiffnessMatrixAndMaterialForce(
    VECTOR& force, vector<TRIPLET>& tripletList, const vector<int>& indexMap)

{
  if (!_solidMaterial) {
    cout << "!!!material is not set!!!" << endl;
    return 0;
  }
  static Real flowRate = SIMPLE_PARSER::getFloat("flow rate", 0);
  static int offsets[8] = {0,
                           1,
                           _xRes,
                           _xRes + 1,
                           _slabSize,
                           _slabSize + 1,
                           _slabSize + _xRes,
                           _slabSize + _xRes + 1};

  vector<int> matOffsets(_materialCopies.size() + 1, 0);
  int startSize = tripletList.size();

  Real energy = 0;
#if USING_OPENMP
#pragma omp parallel
#endif
  {
#if USING_OPENMP
    const int id = omp_get_thread_num();
#else
    const int id = 0;
#endif
    VECTOR localForce(force.size());
    localForce.setZero();
    vector<TRIPLET> localMat;
    mymap<int, MATRIX3> updatedFps;
#if USING_OPENMP
#pragma omp for schedule(static) reduction(+ : energy)
#endif
    for (int z = _tightMin[2]; z < _tightMax[2] - 1; z++)
      for (int y = _tightMin[1]; y < _tightMax[1] - 1; y++)
        for (int x = _tightMin[0]; x < _tightMax[0] - 1; x++) {
          if (_valid(x, y, z) != 1) continue;
          int index = z * _slabSize + y * _xRes + x;

          MATRIX nodalForcesFEM = MATRIX::Constant(3, 8, 0);
          MATRIX stiffness = MATRIX::Constant(24, 24, 0);

          MATRIX coord(3, 8);
          for (int i = 0; i < 8; i++) {
            coord.col(i) = _u[index + offsets[i]];
          }
          VEC3F restPos(0, 0, 0);
          MATRIX X(3, 8);
          for (int i = 0; i < 8; i++) {
            X.col(i) = _X[index + offsets[i]];
            restPos += _X[index + offsets[i]];
          }
          restPos /= 8;

          COROTATION* material = _materialCopies[id];

          // One point quadrature
          MATRIX3 F;
          Real Jinv = 1.0 / computeF(coord, _G, F);
          material->init(F);
          MATRIX3 P = material->firstPiolaKirchhoff();
          MATRIX3 cauchyStress = Jinv * P * F.transpose();

          if (_isPlastic) {
            if (cauchyStress.norm() > _plasticYield) {
              int materialGridIndex = 0;
              MATRIX3 FpInv =
                  _FpInv.interpolate_value(restPos, materialGridIndex);
              F *= FpInv;
              material->init(F);
              P = material->firstPiolaKirchhoff();
              cauchyStress = Jinv * P * F.transpose();
              Real cNorm = cauchyStress.norm();
              if (cNorm > _plasticYield) {
                Real gamma = flowRate * (cNorm - _plasticYield) / cNorm;
                MATRIX3& Fhat = material->Fhat();
                MATRIX3& V = material->V();

                Real detFhatCubicRoot = cbrt(1.0 / max(abs(Jinv), 1e-3));
                MATRIX3 oldFhat = Fhat;
                Fhat(0, 0) = pow(Fhat(0, 0) / detFhatCubicRoot, -gamma);
                Fhat(1, 1) = pow(Fhat(1, 1) / detFhatCubicRoot, -gamma);
                Fhat(2, 2) = pow(Fhat(2, 2) / detFhatCubicRoot, -gamma);

                MATRIX3 deltaFpInv = V * Fhat * V.transpose();
                MATRIX3 previousFpInv = FpInv;
                FpInv = deltaFpInv * FpInv;
                if ((FpInv.array() != FpInv.array()).any()) {
                  cout << "Bad update!!! FpInv has nan!!!" << endl;
                } else
                  updatedFps[materialGridIndex] = FpInv;
              }
            }
          }

          MATRIX pFpxhat = MATRIX::Constant(9, 24, 0);
          computePFPxhat(F, _G, pFpxhat);
          MATRIX GF = _G * F;

          nodalForcesFEM -= cauchyStress * _G.transpose() * _dh * _dh * _dh;

          energy += material->strainEnergy() * _dh * _dh * _dh;

          for (int k = 0; k < 24; k++) {
            int col = index + offsets[k / 3];
            if (_mass[col] == 0) continue;
            MATRIX3 dF;
            dF(0, 0) = pFpxhat(0, k);
            dF(0, 1) = pFpxhat(1, k);
            dF(0, 2) = pFpxhat(2, k);
            dF(1, 0) = pFpxhat(3, k);
            dF(1, 1) = pFpxhat(4, k);
            dF(1, 2) = pFpxhat(5, k);
            dF(2, 0) = pFpxhat(6, k);
            dF(2, 1) = pFpxhat(7, k);
            dF(2, 2) = pFpxhat(8, k);

            MATRIX3 deltaP = material->firstPKDifferential(dF);
            MATRIX deltaForce =
                deltaP * GF.transpose() * Jinv * _dh * _dh * _dh;
            for (int m = 0; m < 8; m++) {
              stiffness(m * 3, k) += deltaForce(0, m);
              stiffness(m * 3 + 1, k) += deltaForce(1, m);
              stiffness(m * 3 + 2, k) += deltaForce(2, m);
            }
          }

          for (int i = 0; i < 8; i++) {
            int idx = indexMap[index + offsets[i]];
            if (idx >= 0)
              localForce.segment<3>(idx * 3) += nodalForcesFEM.col(i);
          }
          for (int i = 0; i < 24; i++) {
            int col = index + offsets[i / 3];
            if (_mass[col] == 0) continue;
            if (indexMap[col] < 0) continue;
            int newCol = indexMap[col] * 3 + i % 3;
            for (int k = 0; k < 24; k++) {
              int row = index + offsets[k / 3];
              if (_mass[row] == 0) continue;
              if (indexMap[row] < 0) continue;
              int newRow = indexMap[row] * 3 + k % 3;
              localMat.push_back(TRIPLET(newRow, newCol, stiffness(k, i)));
            }
          }
        }
    matOffsets[id + 1] = localMat.size();
#if USING_OPENMP
#pragma omp barrier
#pragma omp single
#endif
    {
      for (int t = 1; t < matOffsets.size(); t++)
        matOffsets[t] += matOffsets[t - 1];
      tripletList.resize(matOffsets.back() + startSize);
    }
    vector<TRIPLET>::iterator dest = tripletList.begin() + startSize;
    std::advance(dest, matOffsets[id]);
    std::copy(localMat.begin(), localMat.end(), dest);
#if USING_OPENMP
#pragma omp critical
#endif
    {
      force += localForce;
      for (const pair<int, MATRIX3>& ele : updatedFps)
        _FpInv[ele.first] = ele.second;
    }
  }
  return energy;
}

Real DEFOR_SOLID_3D::computeMV(VECTOR& mv, const vector<int>& indexMap) {
  Real energy = 0;
#if USING_OPENMP
#pragma omp parallel for schedule(static) reduction(+ : energy)
#endif
  for (int z = _zmin; z < _zmax - 1; z++)
    for (int y = _ymin; y < _ymax - 1; y++)
      for (int x = _xmin; x < _xmax - 1; x++) {
        int index = z * _slabSize + y * _xRes + x;
        int idx = indexMap[index];
        if (idx >= 0) {
          mv.segment<3>(idx * 3) += _mass[index] * _velocity[index];
          energy += 0.5 * _mass[index] * _velocity[index].dot(_velocity[index]);
        }
      }
  return energy;
}

void DEFOR_SOLID_3D::gridToParticle() {
  static VEC3F wallNormals[6] = {VEC3F(1, 0, 0),  VEC3F(0, 1, 0),
                                 VEC3F(0, 0, 1),  VEC3F(-1, 0, 0),
                                 VEC3F(0, -1, 0), VEC3F(0, 0, -1)};
  static Real boundaries[6] = {_origin[0],
                               _origin[1],
                               _origin[2],
                               _origin[0] + _dh * (_xRes - 1),
                               _origin[1] + _dh * (_yRes - 1),
                               _origin[2] + _dh * (_zRes - 1)};
  static Real alpha = SIMPLE_PARSER::getFloat("flip alpha", 0.95);
  TIMING_BREAKDOWN::tic();
#if USING_OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (int c = 0; c < _materialParticles.size(); c++) {
    PARTICLE& p = _materialParticles[c];
    Real x = (p.pos[0] - _origin[0]) / _dh;
    Real y = (p.pos[1] - _origin[1]) / _dh;
    Real z = (p.pos[2] - _origin[2]) / _dh;
    int i = std::floor(x) - 1;
    int j = std::floor(y) - 1;
    int k = std::floor(z) - 1;
    VEC3F pic(0, 0, 0);
    VEC3F flip(0, 0, 0);
    Real weights = 0;
    for (int r = k; r < k + 4; r++)
      for (int t = j; t < j + 4; t++)
        for (int s = i; s < i + 4; s++) {
          if (s < 0 || t < 0 || r < 0 || s > _xRes - 1 || t > _yRes - 1 ||
              r > _zRes - 1)
            continue;

          if (_mass(s, t, r) == 0) continue;
          VEC3F diff = _dhInv * (p.pos - _origin);
          diff -= VEC3F(s, t, r);
          VEC3F n(0, 0, 0);
          for (int k = 0; k < 3; k++) {
            if (abs(diff[k]) < 1) {
              n[k] = 0.5 * abs(pow(diff[k], 3)) - diff[k] * diff[k] + 2.0 / 3;
            } else if (abs(diff[k]) < 2) {
              n[k] = -1.0 / 6 * abs(pow(diff[k], 3)) + diff[k] * diff[k] -
                     2 * abs(diff[k]) + 4.0 / 3;
            }
          }
          Real w = n[0] * n[1] * n[2];
          weights += w;
          pic += _velocity(s, t, r) * w;
          flip += (_velocity(s, t, r) - _vecWorkspace(s, t, r)) * w;
        }
    flip += p.velocity;

    p.velocity = (1 - alpha) * pic + alpha * flip;
  }
  TIMING_BREAKDOWN::toc("gridToParticle");

  TIMING_BREAKDOWN::tic();
  Real thickness = SIMPLE_PARSER::getFloat("contact safe distance", _dh * 0.5);
  Real wallThickness = SIMPLE_PARSER::getFloat("wall thickness", _dh * 1.5);
  Real springConstant =
      SIMPLE_PARSER::getFloat("contact spring stiffness", 10000);
  Real friction = SIMPLE_PARSER::getFloat("friction", 0.5);
#if USING_OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
  for (int x = 0; x < _materialParticles.size(); x++) {
    PARTICLE& p = _materialParticles[x];
    VEC3F futurePos = p.pos;
    Real fx = (futurePos[0] - _origin[0]) / _dh;
    Real fy = (futurePos[1] - _origin[1]) / _dh;
    Real fz = (futurePos[2] - _origin[2]) / _dh;
    int i = std::floor(fx);
    int j = std::floor(fy);
    int k = std::floor(fz);
    i = min(max(i, 0), _xRes - 2);
    j = min(max(j, 0), _yRes - 2);
    k = min(max(k, 0), _zRes - 2);
    // handle collision with other solids
    for (int y = 0; y < _neighbors.size(); y++) {
      DEFOR_SOLID_3D& neighbor = *(_neighbors[y]);
      if (neighbor._phi.interpolate_value(futurePos) > 2) continue;
      VEC3F otherRestPos = neighbor._X.interpolate_value(futurePos);
      VEC3F normal;
      Real sdf = neighbor._restSDF->signedDistAndCurrentNormal(
          otherRestPos, normal, thickness);

      if (sdf > thickness) {
        continue;
      }
      Real overlap = thickness - sdf;
      VEC3F v_other = neighbor._velocity.interpolate_value(p.pos);
      VEC3F vrel = p.velocity - v_other;
      Real vn = vrel.dot(normal);

      if (vn <= 0 || sdf <= 0) {
        Real impulse = -min(_dt * springConstant * overlap,
                            p.mass * (0.1 * overlap / _dt - vn));
        Real delta_vn = -impulse / p.mass;
        VEC3F vt = vrel - vn * normal;
        VEC3F new_vt = max(1 - friction * delta_vn / vt.norm(), 0.0) * vt;
        VEC3F new_vrel = (vn + delta_vn) * normal + new_vt;
        p.velocity = new_vrel + v_other;
        break;
      }
    }
    // handle collision with walls
    for (int b = 0; b < 6; b++) {
      Real dist = abs(p.pos[b % 3] - boundaries[b]);
      if (dist < wallThickness) {
        Real overlap = wallThickness - dist;
        int xIndex = i;
        if (b == 0)
          xIndex = 0;
        else if (b == 3)
          xIndex = _xRes - 1;
        int yIndex = j;
        if (b == 1)
          yIndex = 0;
        else if (b == 4)
          yIndex = _yRes - 1;
        int zIndex = k;
        if (b == 2)
          zIndex = 0;
        else if (b == 5)
          zIndex = _zRes - 1;

        VEC3F v_other = _velocity(xIndex, yIndex, zIndex);
        VEC3F vrel = p.velocity - v_other;
        Real vn = vrel.dot(wallNormals[b]);
        if (vn <= 0) {
          Real impulse = -min(_dt * springConstant * overlap,
                              p.mass * (0.1 * overlap / _dt - vn));
          Real delta_vn = -impulse / p.mass;
          VEC3F vt = vrel - vn * wallNormals[b];
          VEC3F new_vt = max(0.1 - friction * delta_vn / vt.norm(), 0.0) * vt;
          p.velocity = (vn + delta_vn) * wallNormals[b] + new_vt + v_other;
          break;
        }
      }
    }
    p.pos += p.velocity * _dt;
  }
  TIMING_BREAKDOWN::toc("collision hanlding");
}

void DEFOR_SOLID_3D::particleToGrid() {
  mymap<int, pair<VEC3F, Real> > indexToWeight;
#if USING_OPENMP
#pragma omp parallel
#endif
  {
    mymap<int, pair<VEC3F, Real> > localMap;
#if USING_OPENMP
#pragma omp for schedule(static)
#endif
    for (int iterator = 0; iterator < _materialParticles.size(); iterator++) {
      PARTICLE& p = _materialParticles[iterator];
      Real x = (p.pos[0] - _origin[0]) / _dh;
      Real y = (p.pos[1] - _origin[1]) / _dh;
      Real z = (p.pos[2] - _origin[2]) / _dh;
      int i = std::floor(x) - 1;
      int j = std::floor(y) - 1;
      int k = std::floor(z) - 1;

      for (int r = k; r < k + 4; r++)
        for (int t = j; t < j + 4; t++)
          for (int s = i; s < i + 4; s++) {
            if (s < 0 || t < 0 || r < 0 || s > _xRes - 1 || t > _yRes - 1 ||
                r > _zRes - 1)
              continue;

            VEC3F diff = _dhInv * (p.pos - _origin);
            diff -= VEC3F(s, t, r);
            VEC3F n(0, 0, 0);
            for (int k = 0; k < 3; k++) {
              if (abs(diff[k]) < 1) {
                n[k] = 0.5 * abs(pow(diff[k], 3)) - diff[k] * diff[k] + 2.0 / 3;
              } else if (abs(diff[k]) < 2) {
                n[k] = -1.0 / 6 * abs(pow(diff[k], 3)) + diff[k] * diff[k] -
                       2 * abs(diff[k]) + 4.0 / 3;
              }
            }
            Real w = n[0] * n[1] * n[2];
            int index = r * _slabSize + t * _xRes + s;

            if (localMap.find(index) == localMap.end()) {
              localMap[index].first = w * p.velocity;
              localMap[index].second = w;
            } else {
              localMap[index].first += w * p.velocity;
              localMap[index].second += w;
            }
          }
    }
#if USING_OPENMP
#pragma omp barrier
#pragma omp critical
#endif
    {
      for (auto iToW : localMap) {
        if (indexToWeight.find(iToW.first) == indexToWeight.end())
          indexToWeight[iToW.first] = iToW.second;
        else {
          indexToWeight[iToW.first].first += iToW.second.first;
          indexToWeight[iToW.first].second += iToW.second.second;
        }
      }
    }
  }
  _mass.clear();
  static int cutoff = SIMPLE_PARSER::getInt("flip cutoff", 3);
  for (auto iToW : indexToWeight) {
    if (iToW.second.second < cutoff) continue;
    _velocity[iToW.first] = iToW.second.first / iToW.second.second;
    _mass[iToW.first] = 1;
  }
}
void DEFOR_SOLID_3D::advectField(VEC3F_FIELD_3D& field,
                                 VEC3F_FIELD_3D& workspace, Real dt) {
  for (int z = 0; z < field.zRes(); z++)
    for (int y = 0; y < field.yRes(); y++)
      for (int x = 0; x < field.xRes(); x++) {
        int index = z * _slabSize + y * _xRes + x;
        VEC3F node(x, y, z);
        node = node * field.dh() + field.origin();
        // VEC3F v = _velocity.interpolate_value(node - 0.5 * _velocity[index] *
        // dt);
        VEC3F newPos = node - _velocity[index] * dt;
        // VEC3F newPos = node - v * dt;
        workspace[index] = field.interpolate_value(newPos);
      }
  field.swap(workspace);
}

void DEFOR_SOLID_3D::advectSurfaceMesh() {
  if (!_surfaceMesh) return;
  static Real epsilon = 1e-4;
  static int offsets[8] = {0,
                           1,
                           _xRes,
                           _xRes + 1,
                           _slabSize,
                           _slabSize + 1,
                           _slabSize + _xRes,
                           _slabSize + _xRes + 1};

  vector<VEC3F>& restPose = _surfaceMesh->restPose;
  vector<VEC3F>& vertices = _surfaceMesh->vertices;
#if USING_OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (unsigned int x = 0; x < vertices.size(); x++) {
    VEC3F copy = vertices[x];
    VEC3F& vertex = vertices[x];
    VEC3F v = _velocity.interpolate_value(vertex);
    vertex += v * _dt;
    Real residual = 1;
    int count = 0;
    while (residual > epsilon && count++ < 10) {
      VEC3F u = _u.interpolate_value(vertex);
      VEC3F prest = vertex - u;
      VEC3F diff = prest - restPose[x];
      residual = diff.norm();

      int i = (vertex[0] - _origin[0]) / _dh;
      if (i < 0) i = 0;
      if (i > _xRes - 2) i = _xRes - 2;
      int j = (vertex[1] - _origin[1]) / _dh;
      if (j < 0) j = 0;
      if (j > _yRes - 2) j = _yRes - 2;
      int k = (vertex[2] - _origin[2]) / _dh;
      if (k < 0) k = 0;
      if (k > _zRes - 2) k = _zRes - 2;

      int index = k * _slabSize + j * _xRes + i;
      MATRIX coord(3, 8);
      for (int i = 0; i < 8; i++) coord.col(i) = _u[index + offsets[i]];
      MATRIX3 F;
      computeF(coord, _G, F);
      vertex -= F * diff;
    }
    vertex[0] = max(_origin[0], vertex[0]);
    vertex[1] = max(_origin[1], vertex[1]);
    vertex[2] = max(_origin[2], vertex[2]);

    vertex[0] = min(_origin[0] + _dh * (_xRes - 1), vertex[0]);
    vertex[1] = min(_origin[1] + _dh * (_yRes - 1), vertex[1]);
    vertex[2] = min(_origin[2] + _dh * (_zRes - 1), vertex[2]);
  }
}

void DEFOR_SOLID_3D::finalizeSolution(VEC3F_FIELD_3D& velocity, Real dt) {
  _velocity.clear();
#if USING_OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (int z = _zmin; z < _zmax; z++)
    for (int y = _ymin; y < _ymax; y++)
      for (int x = _xmin; x < _xmax; x++) {
        // if(_mass(x, y, z) == 0)
        // continue;
        int index = _slabSize * z + y * _xRes + x;
        _velocity[index] = velocity[index];
      }
  _dt = dt;
  // finalizeSolution();
}
void DEFOR_SOLID_3D::finalizeSolution() {
  // Advect displacement
  for (int i = 0; i < _totalCells; i++) {
    _u[i] += _velocity[i] * _dt;
  }
  advectField(_u, _vecWorkspace, _dt);

  // Adjust boundary displacement
  for (int x = _xmin; x < _xmax; x++)
    for (int y = _ymin; y < _ymax; y++) {
      _u(x, y, _zmin) = 2 * _u(x, y, _zmin + 1) - _u(x, y, _zmin + 2);
      _u(x, y, _zmax - 1) = 2 * _u(x, y, _zmax - 2) - _u(x, y, _zmax - 3);
    }
  for (int x = _xmin; x < _xmax; x++)
    for (int z = _zmin; z < _zmax; z++) {
      _u(x, _ymin, z) = 2 * _u(x, _ymin + 1, z) - _u(x, _ymin + 2, z);
      _u(x, _ymax - 1, z) = 2 * _u(x, _ymax - 2, z) - _u(x, _ymax - 3, z);
    }
  for (int y = _ymin; y < _ymax; y++)
    for (int z = _zmin; z < _zmax; z++) {
      _u(_xmin, y, z) = 2 * _u(_xmin + 1, y, z) - _u(_xmin + 2, y, z);
      _u(_xmax - 1, y, z) = 2 * _u(_xmax - 2, y, z) - _u(_xmax - 3, y, z);
    }

  // Recover material position field
  int index = 0;
  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++, index++) {
        VEC3F pos(x * _dh, y * _dh, z * _dh);
        pos += _origin - _u[index];
        _X[index] = pos;
      }

  _vecWorkspace = _velocity;
}

void DEFOR_SOLID_3D::recoverVelocityField(VEC3F_FIELD_3D& v) {
  for (int z = _zmin; z < _zmax - 1; z++)
    for (int y = _ymin; y < _ymax - 1; y++)
      for (int x = _xmin; x < _xmax - 1; x++) {
        if (_mass(x, y, z) > 0) v(x, y, z) = _velocity(x, y, z);
      }
}

void DEFOR_SOLID_3D::drawMaterialParticles() {
  glColor3f(0, 0, 1);
  glPointSize(1.f);
  glBegin(GL_POINTS);
  for (const PARTICLE& p : _materialParticles) {
    glVertex3f(p.pos[0], p.pos[1], p.pos[2]);
  }
  glEnd();
}
void DEFOR_SOLID_3D::drawMass() {
  glColor3f(1, 0, 0);
  glPointSize(5.f);
  glBegin(GL_POINTS);
  for (int z = _zmin; z < _zmax - 1; z++)
    for (int y = _ymin; y < _ymax - 1; y++) {
      for (int x = _xmin; x < _xmax - 1; x++) {
        VEC3F node(x * _dh + _origin[0], y * _dh + _origin[1],
                   z * _dh + _origin[2]);
        if (_mass(x, y, z) != 0) {
          glColor3f(1, 0, 0);
          glVertex3f(node[0], node[1], node[2]);
        }
      }
    }
  glPointSize(10.f);
  for (int index : _moniters) {
    int z = index / _slabSize;
    int y = (index - z * _slabSize) / _xRes;
    int x = (index - z * _slabSize) % _xRes;
    VEC3F node(x * _dh + _origin[0], y * _dh + _origin[1],
               z * _dh + _origin[2]);
    glColor3f(0, 1, 0);
    glVertex3f(node[0], node[1], node[2]);
  }
  glEnd();
}
void DEFOR_SOLID_3D::drawVelocity() {
  glBegin(GL_LINES);
  for (int z = _zmin; z < _zmax; z++)
    for (int y = _ymin; y < _ymax; y++)
      for (int x = _xmin; x < _xmax; x++) {
        if (_mass(x, y, z) == 0)
          continue;
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

void DEFOR_SOLID_3D::readFrame(const string& filename, int frame, bool debug) {
  if (debug || (frame % 10 == 0)) {
    _u.readGz(filename + ".displacement");
    _velocity.readGz(filename + ".velocity");
    _phi.readGz(filename + ".phi");
    if (_isPlastic) _FpInv.readGz(filename + ".plasticdefor");
  }
  _vecWorkspace = _velocity;
  int index = 0;
  for (int z = 0; z < _zRes; z++)
    for (int y = 0; y < _yRes; y++)
      for (int x = 0; x < _xRes; x++, index++) {
        VEC3F pos(x * _dh, y * _dh, z * _dh);
        pos += _origin;
        _X[index] = pos - _u[index];
      }
  OBJ tmpMesh;
  tmpMesh.Load(filename + ".obj");
  _surfaceMesh->vertices = tmpMesh.vertices;

  vector<VEC3F>& vertices = _surfaceMesh->vertices;
  VEC3F boxMin(1e9, 1e9, 1e9);
  VEC3F boxMax(-1e9, -1e9, -1e9);
  for (int x = 0; x < vertices.size(); x++) {
    boxMin[0] = min(boxMin[0], vertices[x][0]);
    boxMin[1] = min(boxMin[1], vertices[x][1]);
    boxMin[2] = min(boxMin[2], vertices[x][2]);

    boxMax[0] = max(boxMax[0], vertices[x][0]);
    boxMax[1] = max(boxMax[1], vertices[x][1]);
    boxMax[2] = max(boxMax[2], vertices[x][2]);
  }

  _xmin = (boxMin[0] - _origin[0]) / _dh;
  _ymin = (boxMin[1] - _origin[1]) / _dh;
  _zmin = (boxMin[2] - _origin[2]) / _dh;

  _xmax = (boxMax[0] - _origin[0]) / _dh;
  _ymax = (boxMax[1] - _origin[1]) / _dh;
  _zmax = (boxMax[2] - _origin[2]) / _dh;

  static int tightBorder = SIMPLE_PARSER::getInt("tight border", 3);
  static int extrap = SIMPLE_PARSER::getInt("extrapolation cells", 20);

  _tightMin = VEC3I(max(_xmin - tightBorder, 0), max(_ymin - tightBorder, 0),
                    max(_zmin - tightBorder, 0));

  _tightMax =
      VEC3I(min(_xmax + tightBorder, _xRes), min(_ymax + tightBorder, _yRes),
            min(_zmax + tightBorder, _zRes));

  _xmin = max(_xmin - extrap, 0);
  _ymin = max(_ymin - extrap, 0);
  _zmin = max(_zmin - extrap, 0);

  _xmax = min(_xmax + extrap, _xRes);
  _ymax = min(_ymax + extrap, _yRes);
  _zmax = min(_zmax + extrap, _zRes);

  computeSolidMass(true);
  initMaterialParticles();
}
void DEFOR_SOLID_3D::writeFrame(const string& filename, int frame, bool debug) {
  if (debug || (frame % 10 == 0)) {
    _u.writeGz(filename + ".displacement");
    _velocity.writeGz(filename + ".velocity");
    _phi.writeGz(filename + ".phi");
    if (_isPlastic) _FpInv.writeGz(filename + ".plasticdefor");
  }
  _surfaceMesh->Save(filename + ".obj");
  _surfaceMesh->ComputeVertexNormals();
  _surfaceMesh->SaveToPBRT(filename + ".pbrt");
}