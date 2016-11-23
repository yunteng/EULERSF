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
#ifndef IO_H
#define IO_H

#include <string>

class IO {
 public:
  ////////////////////////////////////////////////////////////////
  // Print integer to a zero-padded string
  //////////////////////////////////////////////////////////////////
  static std::string itoPaddedString(int frame) {
    char buffer[256];
    sprintf(buffer, "%i", frame);

    std::string number = std::string(buffer);
    if (frame < 10) number = std::string("0") + number;
    if (frame < 100) number = std::string("0") + number;
    if (frame < 1000) number = std::string("0") + number;

    return number;
  }
};

#endif
