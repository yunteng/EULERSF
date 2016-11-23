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
//
// Simple config file parser
//
// Original source code was provided courtesy of
// Nils Theurey (www.ntoken.com)
//
// Minor modifications started on 6/27/2008 by Ted Kim
//
//////////////////////////////////////////////////////////////////////

#ifndef SIMPLE_PARSER_H
#define SIMPLE_PARSER_H

#include <SETTINGS.h>
#include <util/IO.h>
#include <iostream>
#include <map>
#include <string>
#include <vector>
// #include <VEC3.h>

using namespace std;

// read the whole test file of pairs
// a = b
// and look up values with conversion
class SIMPLE_PARSER {
 public:
  static bool parse(std::string file);
  // get values from file
  static int getInt(string name, int defaultValue, bool needed = false);
  static bool getBool(string name, bool defaultValue, bool needed = false);
  static double getFloat(string name, double defaultValue, bool needed = false);
  static string getString(string name, string defaultValue,
                          bool needed = false);
  static VEC3F getVEC3F(string name, VEC3F defaultValue, bool needed = false);
  static VEC2F getVEC2F(string name, VEC2F defaultValue, bool needed = false);
  // check if there were any unused pairs
  static bool haveUnusedValues();
  static string printAllUnused();

  // check if the a parameters was specified in the config file
  static bool defined(string name);

  // static MATERIAL *READ_MATERIAL(string filename);

 protected:
  // value pair storage
  static map<string, string> mVals;
  // value pair check if used...
  static map<string, bool> mUsed;

  template <class T>
  T static getScalarValue(string name, T defaultValue, bool needed = false);
};

#endif
