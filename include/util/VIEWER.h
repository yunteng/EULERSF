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
#ifndef VIEWER_H
#define VIEWER_H

#include <SETTINGS.h>
#include <util/SIMPLE_PARSER.h>
#include <util/TIMING_BREAKDOWN.h>
#include <glvu.hpp>
#include <iostream>

using namespace std;
template <class T>
class VIEWER {
 public:
  static void init();
  static void screenshot(string renderPath, int frame);
  ~VIEWER() { delete simulator; };

 public:
  static T* simulator;
  static GLVU glvu;
  static bool showGrid;
  static bool animate;
  static bool step;

 private:
  static void displayFunc();
  static void keyboardFunc(unsigned char Key, int x, int y);
  static void idleFunc();
  static void mouseFunc(int button, int state, int x, int y);
  static void motionFunc(int x, int y);
  static int glvuWindow(glvuVec3f bboxCenter, Real eyeDistanceScale = 1.0);
  static void drawGrid();
  static VEC3F unproject(float x, float y, float z);

 private:
  static const int windowStartX = 0;
  static const int windowStartY = 0;
  static const int windowWidth = 700;
  static const int windowHeight = 700;
};
#include "VIEWER.inl"

#endif
