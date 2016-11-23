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
#include <SETTINGS.h>
#include <dtgrid/SPARSE_SDF.h>
#include <geometry/FIELD_3D.h>
#include <geometry/OBJ.h>
#include <util/SIMPLE_PARSER.h>
#include <iostream>

using namespace std;

int main(int argc, char* argv[]) {
  if (argc < 2) {
    cout << "./bin/rasterizeSolid *.cfg" << endl;
    exit(0);
  }
  string configName(argv[1]);
  if (!SIMPLE_PARSER::parse(configName)) {
    cout << "Cannot open config file " << configName << endl;
    exit(0);
  }
  string meshfilename = SIMPLE_PARSER::getString("mesh file name", "");
  string sdffilename = SIMPLE_PARSER::getString("sdf file name", "");

  OBJ mesh;
  mesh.Load(meshfilename);
  SPARSE_SDF sdf;
  sdf.load(sdffilename);

  VEC3F boxmin, boxmax;
  mesh.BoundingBox(boxmin, boxmax);

  Real mass_dh = SIMPLE_PARSER::getFloat("mass dh", 0.01);
  VEC3F expand(mass_dh * 10, mass_dh * 10, mass_dh * 10);
  boxmin -= expand;
  boxmax += expand;

  int xRes = (boxmax[0] - boxmin[0]) / mass_dh;
  int yRes = (boxmax[1] - boxmin[1]) / mass_dh;
  int zRes = (boxmax[2] - boxmin[2]) / mass_dh;

  Real density = SIMPLE_PARSER::getFloat("solid mass density", 1.0);
  FIELD_3Df mass(xRes, yRes, zRes, mass_dh, boxmin);
  int index = 0;
  Real sumMass = 0;
  VEC3F dummyNormal;
  for (int z = 0; z < zRes; z++)
    for (int y = 0; y < yRes; y++)
      for (int x = 0; x < xRes; x++, index++) {
        VEC3F pos(x, y, z);
        pos *= mass_dh;
        pos += boxmin;
        Real sd = sdf.signedDist(pos, dummyNormal);
        if (sd <= mass_dh / 2) {
          mass[index] = density;
          sumMass += density * mass_dh * mass_dh * mass_dh;
        }
      }

  cout << "sum mass " << sumMass << endl;

  string massfilename =
      SIMPLE_PARSER::getString("mass file name", meshfilename + ".mass");
  mass.writeGz(massfilename);

  return 0;
}