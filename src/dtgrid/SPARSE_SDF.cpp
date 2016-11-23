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
#include <dtgrid/SPARSE_SDF.h>

int SPARSE_SDF::trianglesAt(const VEC3F& position, TRIANGLE**& triangles) {
  MyLocator loc;
  MyVec3 pos(position[0], position[1], position[2]);
  if (_grid->getLocator(pos, &loc)) {
    triangles = _hashSurface[loc.iv3D];
    return _gridTriangleNum[loc.iv3D];
  } else
    return 0;
  // int n = _gridTriangleNum[loc.iv3D];
  // return n;
}
SPARSE_SDF::~SPARSE_SDF() {
  for (unsigned int x = 0; x < _hashSurface.size(); x++) {
    if (_hashSurface[x] != NULL) {
      delete[] _hashSurface[x];
    }
  }
  delete _grid;
}
void SPARSE_SDF::hashTriangle(TRIANGLE& triangle,
                              map<int, vector<TRIANGLE*> >& tempHash) {
  VEC3F mins, maxs;
  triangle.boundingBox(mins, maxs);
  MyVec3 dtMin(mins[0], mins[1], mins[2]);
  MyVec3 dtMax(maxs[0], maxs[1], maxs[2]);

  Index iMins[3];
  Index iMaxs[3];
  _grid->boundingBoxIndices(dtMin, dtMax, iMins, iMaxs);

  for (Index z = iMins[2]; z <= iMaxs[2]; z++)
    for (Index y = iMins[1]; y <= iMaxs[1]; y++)
      for (Index x = iMins[0]; x <= iMaxs[0]; x++) {
        MyLocator loc;
        bool gridExist = _grid->getLocator(x, y, z, &loc);
        if (gridExist) tempHash[loc.iv3D].push_back(&triangle);
      }
}

void SPARSE_SDF::initHashSurface(vector<TRIANGLE>& surfaces) {
  _hashSurface.resize(_grid->getNumVa3D(), NULL);
  _gridTriangleNum.resize(_grid->getNumVa3D(), 0);

  map<int, vector<TRIANGLE*> > tempHash;
  for (int i = 0; i < surfaces.size(); i++) {
    hashTriangle(surfaces[i], tempHash);
  }
  for (map<int, vector<TRIANGLE*> >::iterator iter = tempHash.begin();
       iter != tempHash.end(); iter++) {
    int index = iter->first;
    vector<TRIANGLE*>& triangles = iter->second;
    // cout << "grid has " << triangles.size() << " triangles " << endl;
    _hashSurface[index] = new TRIANGLE*[triangles.size()];
    _gridTriangleNum[index] = triangles.size();
    for (unsigned int x = 0; x < triangles.size(); x++) {
      _hashSurface[index][x] = triangles[x];
    }
  }
}

Real SPARSE_SDF::signedDistAndRestNormal(const VEC3F& point, VEC3F& normal,
                                         Real safeDist) {
  VEC3F restSurfacePos;
  Real sd = nearestSurfacePoint(point, restSurfacePos, safeDist);
  if (sd >= safeDist) return sd;
  TRIANGLE** surfaceTriangles = NULL;
  int numTris = trianglesAt(restSurfacePos, surfaceTriangles);
  // cout << numTris << endl;
  for (int tri = 0; tri < numTris; tri++) {
    TRIANGLE*& triangle = surfaceTriangles[tri];
    VEC3F lambda;
    if (triangle->baryCenter(restSurfacePos, lambda, true)) {
      normal = triangle->restNormal();
      return sd;
    }
  }
  return 1;
}

Real SPARSE_SDF::signedDistAndCurrentNormal(const VEC3F& point, VEC3F& normal,
                                            Real safeDist) {
  VEC3F restSurfacePos;
  Real sd = nearestSurfacePoint(point, restSurfacePos, safeDist);
  if (sd > safeDist) {
    return 1;
  }
  TRIANGLE** surfaceTriangles = NULL;
  int numTris = trianglesAt(restSurfacePos, surfaceTriangles);
  // cout << numTris << endl;
  for (int tri = 0; tri < numTris; tri++) {
    TRIANGLE*& triangle = surfaceTriangles[tri];
    VEC3F lambda;
    if (triangle->baryCenter(restSurfacePos, lambda, true)) {
      normal = triangle->normal();
      return sd;
    }
  }
  return 1;
}

Real SPARSE_SDF::nearestSurfacePoint(const VEC3F& position, VEC3F& surface,
                                     Real safeDist) {
  MyVec3 vecNormal;
  Real sdf;

  MyVec3 positionVec(position[0], position[1], position[2]);
  MyVec3 tmpSurface = positionVec;

  bool inNarrowBand =
      _grid->lookupDistanceAndNormalFast(positionVec, sdf, vecNormal);

  if (sdf >= safeDist) return sdf;

  Real initSdf = sdf;
  int iter = 0;

  while (abs(sdf) > 1e-9 && iter < 100) {
    tmpSurface += (-sdf) * vecNormal;
    sdf = (*_grid)(tmpSurface);
    iter++;
  }

  surface[0] = tmpSurface[0];
  surface[1] = tmpSurface[1];
  surface[2] = tmpSurface[2];
  return sdf;
}
Real SPARSE_SDF::signedDist(const VEC3F& position, VEC3F& normal) {
  MyVec3 vecNormal;
  Real sdf;
  MyVec3 positionVec(position[0], position[1], position[2]);

  bool inNarrowBand =
      _grid->lookupDistanceAndNormalFast(positionVec, sdf, vecNormal);
  normal[0] = vecNormal[0];
  normal[1] = vecNormal[1];
  normal[2] = vecNormal[2];
  return sdf;
}
