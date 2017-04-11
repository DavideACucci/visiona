/*******************************************************************************
**
** Copyright (C) 2015-2019 EPFL (Swiss Federal Institute of Technology)
**
** Contact:
**   Dr. Davide A. Cucci, post-doctoral researcher
**   E-mail: davide.cucci@epfl.ch
**
**   Geodetic Engineering Laboratory,
**   1015 Lausanne, Switzerland (www.topo.epfl.ch).
**
**
**
**   This file is part of visiona.
**
**   visiona is free software: you can redistribute it and/or modify
**   it under the terms of the GNU General Public License as published by
**   the Free Software Foundation, either version 3 of the License, or
**   (at your option) any later version.
**
**   visiona is distributed in the hope that it will be useful,
**   but WITHOUT ANY WARRANTY; without even the implied warranty of
**   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
**   GNU General Public License for more details.
**
**   You should have received a copy of the GNU General Public License
**   along with visiona.  If not, see <http://www.gnu.org/licenses/>.
**
*******************************************************************************/

#include "Timer.h"

#include <ctime>
#include <vector>

#include <iostream>

namespace visiona {

std::vector<clock_t> timers;

void tic() {
//  std::cerr << "TIC level " << timers.size() << std::endl;

  timers.push_back(clock());
}

double toc() {
//  std::cerr << "TOC level " << timers.size() << std::endl;

  clock_t last = timers.back();
  timers.pop_back();
  return (double)(clock()-last)/CLOCKS_PER_SEC;
}

} /* namespace orbitlib */
