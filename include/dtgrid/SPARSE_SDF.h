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
#ifndef SPARSE_SDF_H
#define SPARSE_SDF_H

#include <SETTINGS.h>

#include <string>
#include <vector>

#include <Algorithms/LevelSet/HeaderFiles/LevelSet3D.h>
#include <Algorithms/LevelSet/ScalarFields/HeaderFiles/MeanCurvatureFlowScalarField3D.h>
#include <Algorithms/Math/Interpolators/HeaderFiles/TrilinearInterpolator.h>
#include <DataStructures/Grids/DTGrid/HeaderFiles/DTGrid.h>

#include <geometry/TRIANGLE.h>
using namespace std;

// typedef double Real;
typedef Real Data;
typedef short Index;
typedef unsigned int UInt;

/** The same as Grids::DTGridTraitsDefault except that this Trait allows
 * processing without a safe band (includeSafeBand=false).
  * Note that this slows down stencil iteration somewhat and should not be
 * enabled if not needed.
  * Processing without a safe band is required if rebuild() has not been called
 * on the grid before iterating over it using stencil iterators.
  */
template <typename Data, typename Index, typename Real, typename UInt>
struct DTGridTraitsClosedLevelSet {
  typedef Data DataType;
  typedef Index IndexType;
  typedef Real RealType;
  typedef Data* DataPtr;
  typedef UInt UIntType;

  /** The width of dilation used when rebuilding the narrow band (corresponds to
   * the max number of voxels the interface can move between rebuilds) */
  static const UInt rebuildDilationWidth = 1;
  /** If the region includes an extra band of voxels around the gamma tube
   * (added by a narrowBandRebuild() operation), some optimizations can be
   * applied for iterators that rely on the DTGrid storing a signed distance
   * field, and includeSafeBand should be defined. Note that when using the
   * rebuild method supplied with the DT-Grid, the narrow band will always
   * include a safeband! Furthermore note that when using FMM-based
   * reinitialization, includeSafeBand must be disabled because FMM will not add
   * a safeband. */
  static const bool includeSafeBand = false;
  /** If this is defined, the calls to increment operations in the
   * updateStencil() method are inlined directly in the updateStencilMethod.
   * Improves performance slightly. */
  static const bool inlineUpdateStencilCalls = true;
  /** The maximum number of grid points allowed in a dilation operation. Used to
   * allocate static arrays. */
  static const int maxDilationWidth = 10;
  static const Grids::DTGridSearchType randomAccessType =
      Grids::USE_LINEAR_SEARCH;
  // static const Grids::DTGridSearchType randomAccessType =
  // Grids::USE_BINARY_SEARCH;
  /** If this is false, linear array indices in off-center stencil iterators are
   * not maintained (and hence cannot be used). Speeds up iteration */
  static const bool maintainArrayIndexInStencil = false;
  /** The width in grid cells of the zero crossing band */
  static const int zeroCrossingWidth = 1;
  /** If this boolean is true, open level sets are supported in random access.
   *  Open level sets are always supported for iteration with a stencil
   * iterator. */
  static const bool openLevelSetsSupport = false;
};

typedef DTGridTraitsClosedLevelSet<Data, Index, Real, UInt>
    MyTraitsClosedLevelSet;

typedef Grids::DTGrid<MyTraitsClosedLevelSet> MyGrid;
typedef MyGrid::InitParams MyInitParams;
typedef Matrix::Matrix4x4<Real> MyMatrix;
typedef Matrix::Vector3<Real> MyVec3;
typedef MyGrid::Locator MyLocator;

class SPARSE_SDF {
 public:
  SPARSE_SDF() { _grid = new MyGrid(_initParams); }
  ~SPARSE_SDF();
  void load(string inputFileName) { _grid->loadSparseVolume(inputFileName); }

  void initHashSurface(vector<TRIANGLE>& surfaces);
  Real signedDist(const VEC3F& position, VEC3F& normal);
  Real nearestSurfacePoint(const VEC3F& position, VEC3F& surface,
                           Real safeDist);
  Real signedDistAndCurrentNormal(const VEC3F& point, VEC3F& normal,
                                  Real safeDist);
  Real signedDistAndRestNormal(const VEC3F& point, VEC3F& normal,
                               Real safeDist);

  int trianglesAt(const VEC3F& position, TRIANGLE**& triangles);

  bool getLocator(const Real pos[3], MyLocator* loc) {
    return _grid->getLocator(pos, loc);
  };

 private:
  MyGrid* _grid;
  MyInitParams _initParams;
  vector<TRIANGLE**> _hashSurface;
  vector<Index> _gridTriangleNum;

  // compute the bounding box of triangle and push it to the grid cells that it
  // intersects with
  void hashTriangle(TRIANGLE& triangle, map<int, vector<TRIANGLE*> >& tempHash);
};
#endif
