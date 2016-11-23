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
#include <integrator/COUPLED_INTEGRATOR_3D.h>
#include <integrator/DEFOR_SOLID_3D.h>
#include <util/IO.h>
#include <util/VIEWER.h>
#include <Eigen/Geometry>
#include <iostream>
#if USING_OPENMP
#include <omp.h>
#endif

using namespace std;

class APPLICATION {
 public:
  APPLICATION(string config) : configName(config), integrator(NULL) {
    if (!SIMPLE_PARSER::parse(configName)) exit(0);

    renderPath = SIMPLE_PARSER::getString("render path", "");
    dataPath = SIMPLE_PARSER::getString("data path", "");

    // create working directories
    string mkdirRender = string("mkdir -p ") + renderPath;
    system(mkdirRender.c_str());
    string mkdirData = string("mkdir -p ") + dataPath;
    system(mkdirData.c_str());

    startFrame = SIMPLE_PARSER::getInt("start frame", 0);
    endFrame = SIMPLE_PARSER::getInt("end frame", 100);
    currentFrame = startFrame;

    drawParticles = drawMass = false;
    drawSurface = true;
    drawSmoke = true;

    debug = SIMPLE_PARSER::getBool("debug", false);
  };
  ~APPLICATION() {
    if (integrator) delete integrator;
  }
  int xRes() { return integrator->xRes(); };
  int yRes() { return integrator->yRes(); };
  int zRes() { return integrator->zRes(); };
  Real dh() { return integrator->dh(); };
  const VEC3F& origin() const { return integrator->origin(); };

  void init() {
    int xRes = SIMPLE_PARSER::getInt("xRes", 100);
    int yRes = SIMPLE_PARSER::getInt("yRes", 100);
    int zRes = SIMPLE_PARSER::getInt("zRes", 100);
    Real dh = SIMPLE_PARSER::getFloat("dh", 0.01);
    VEC3F origin = SIMPLE_PARSER::getVEC3F("origin", VEC3F(0, 0, 0));
    integrator = new COUPLED_INTEGRATOR_3D(xRes, yRes, zRes, dh, origin);

    addSolids();

    Real fluidDensity = SIMPLE_PARSER::getFloat("fluid density", 1);
    Real fluidViscocity = SIMPLE_PARSER::getFloat("fluid viscocity", 0);
    integrator->setFluidDensity(fluidDensity);
    integrator->setFluidViscocity(fluidViscocity);

    integrator->init();

    Real carpetThickness =
        SIMPLE_PARSER::getFloat("smoke carpet thickness", 0.06);
    int carpetLayers = carpetThickness / dh;
    // Generate a carpet of smoke
    FIELD_3Df& intensity = integrator->intensity();
    for (int z = 1; z < zRes - 1; z++)
      for (int y = 1; y < carpetLayers; y++)
        for (int x = 1; x < xRes - 1; x++) {
          VEC2F pos(x * dh + origin[0], z * dh + origin[2]);
          intensity(x, y, z) = 0.8f;
        }

    if (startFrame > 0) integrator->readFrame(startFrame - 1, debug);

#if USING_OPENMP
    cout << "Using openmp, number of threads " << omp_get_max_threads() << endl;
#endif
  }

  void addSolids() {
    // initialize the solid objects
    int numberOfSolids = SIMPLE_PARSER::getInt("number of solids", 1);
    for (int i = 0; i < numberOfSolids; i++) {
      string index = to_string(i);
      string meshPath = SIMPLE_PARSER::getString("solid path " + index, "");
      string meshName = SIMPLE_PARSER::getString("solid name " + index, "");
      meshPath += meshName;
      // load mass
      FIELD_3Df restMass;
      restMass.readGz(meshPath + ".mass");
      DEFOR_SOLID_3D* solid =
          new DEFOR_SOLID_3D(xRes(), yRes(), zRes(), dh(), origin(), restMass);

      // assign a corotational material
      Real lambda = SIMPLE_PARSER::getFloat("material lambda " + index, 100);
      Real mu = SIMPLE_PARSER::getFloat("material mu " + index, 100);
      solid->setMaterial(lambda, mu);
      Real density = SIMPLE_PARSER::getFloat("solid mass density " + index, 1);
      solid->setMassDensity(density);

      Real plasticYield =
          SIMPLE_PARSER::getFloat("solid plastic yield " + index, -1);
      if (plasticYield > 0) {
        solid->setPlasticYield(plasticYield);
      }

      // set world transformation
      VEC3F translation =
          SIMPLE_PARSER::getVEC3F("solid translation " + index, VEC3F(0, 0, 0));
      VEC3F scaling =
          SIMPLE_PARSER::getVEC3F("solid scaling " + index, VEC3F(1, 1, 1));
      VEC3F rotation =
          SIMPLE_PARSER::getVEC3F("solid rotation " + index, VEC3F(0, 0, 0));
      MATRIX3 rot;
      rot = Eigen::AngleAxisd(rotation[0] * M_PI, VEC3F::UnitX()) *
            Eigen::AngleAxisd(rotation[1] * M_PI, VEC3F::UnitY()) *
            Eigen::AngleAxisd(rotation[2] * M_PI, VEC3F::UnitZ());

      solid->setInitialState(rot, translation, scaling);

      VEC3F velocity = SIMPLE_PARSER::getVEC3F(
          "solid initial velocity " + index, VEC3F(0, 0, 0));
      solid->setInitialVelocity(velocity);

      // transform the surface mesh
      solid->loadSurfaceMesh(meshPath + ".obj");

      string sdffilename(meshPath + ".svol");
      solid->loadRestSDF(sdffilename);

      OBJ* surfaceMesh = solid->surfaceMesh();
      surfaceMesh->scale(scaling);
      surfaceMesh->transformMesh(rot);
      surfaceMesh->translateMesh(translation);
      surfaceMesh->ComputeVertexNormals();

      integrator->addSolid(solid);

      mySolid = solid;
    }
  }

  void display() {
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

    if (drawParticles) integrator->drawSolidParticles();

    if (drawMass) integrator->drawSolidMasses();

    if (drawSurface) integrator->drawSolidSurfaces();

    if (drawSmoke) integrator->drawIntensity();
  }

  void step() {
    if (currentFrame >= endFrame) exit(0);

    cout << "================= Frame " << currentFrame
         << " =================\n";

    Real frameTimeStep = SIMPLE_PARSER::getFloat("dt", 0.005);

    Real simTimeStep = 0;
    Real energy;
    int subStep = 0;
    while (simTimeStep < frameTimeStep) {
      cout << "========== Sub step " << subStep << " ==========\n";
      Real maxDt = frameTimeStep - simTimeStep;
      integrator->setMaxDt(maxDt);

      TIMING_BREAKDOWN::startFrame();

      energy = integrator->stepImplicit();

      TIMING_BREAKDOWN::endFrame();

      simTimeStep += integrator->dt();
      subStep++;
    }
    TIMING_BREAKDOWN::printTimingBreakdown();

    static bool saveData = SIMPLE_PARSER::getBool("save data", false);
    if (saveData) integrator->writeFrame(currentFrame, debug);

    currentFrame++;
  }

  void keyboardFunc(unsigned char key) {
    switch (key) {
      case 'p':
        TIMING_BREAKDOWN::printTimingBreakdown();
        break;
      case 'v':
        drawParticles = !drawParticles;
        break;
      case 'o':
        drawSmoke = !drawSmoke;
        break;
      case 'i':
        drawSurface = !drawSurface;
        break;
      case 'm':
        drawMass = !drawMass;
        break;
    }
  }

  void click(VEC3F point) {}

 public:
  string configName;
  COUPLED_INTEGRATOR_3D* integrator;
  DEFOR_SOLID_3D* mySolid;
  bool drawParticles, drawSmoke, drawSurface, drawMass;

  bool debug;

  string renderPath;
  string dataPath;

  int currentFrame;
  int startFrame;
  int endFrame;
};

template <class T>
GLVU VIEWER<T>::glvu;

template <class T>
T* VIEWER<T>::simulator = NULL;

template <class T>
bool VIEWER<T>::animate = false;

template <class T>
bool VIEWER<T>::step = true;

template <class T>
bool VIEWER<T>::showGrid = true;

int main(int argc, char* argv[]) {
  if (argc < 2) {
    cout << "./bin/ball *.cfg" << endl;
    exit(0);
  }
  Eigen::initParallel();
  string configName(argv[1]);
  VIEWER<APPLICATION>::simulator = new APPLICATION(configName);
  VIEWER<APPLICATION>::simulator->init();
  VIEWER<APPLICATION>::init();
  return 0;
}
