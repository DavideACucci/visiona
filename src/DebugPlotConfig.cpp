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

/*
 * DebugPlotConfig.cpp
 *
 *  Created on: Mar 20, 2015
 *      Author: Davide A. Cucci (davide.cucci@epfl.ch)
 */

#include "DebugPlotConfig.h"

using namespace std;

namespace visiona {

TimeStats::TimeStats() {
  memset(this, 0, sizeof(struct TimeStats));
}

DebugPlotConfig::DebugPlotConfig() :
    enabled(false),
        enableEdgeDetection(false),
        enableContours(false),
        enableFilteredContours(false),
        enableCirclesClusters(false),
        enableSelectedTargets(false),
        enableEsposure(false),
        enableRoughMeasure(false),
        enableSubPixelEllipses(false),
        blitSubRegion(false),
        frameNumber(0), writeContours(false), writeSubpixelContours(
            false) {
}

string writeTimeWithPercent(double time, double wrt) {
  stringstream s;

  s << fixed << setprecision(3) << time * 1e3
      << " (" << fixed << setprecision(2)
      << time * 100.0 / wrt << "%)";

  return s.str();
}

ostream &operator<<(ostream &s, const TimeStats &ts) {

  s << ts.tEdgeDetection << ", " << ts.tContours << ", " << ts.tFilteredContours
      << ", " << ts.tClusterCircles << ", " << ts.tSelectMarker << ", "
      << ts.tTotalDetection << ", " << ts.tExposure << ", " << ts.tMeasureRough
      << ", " << ts.tMeasure << endl;

  /* for output on the screen
   s << "ED:  " << writeTimeWithPercent(ts.tEdgeDetection, ts.tTotal) << endl;
   s << "CT:  " << writeTimeWithPercent(ts.tContours, ts.tTotal) << endl;
   s << "FCT: " << writeTimeWithPercent(ts.tFilteredContours, ts.tTotal) << endl;
   s << "CC:  " << writeTimeWithPercent(ts.tClusterCircles, ts.tTotal) << endl;
   s << "SM:  " << writeTimeWithPercent(ts.tSelectMarker, ts.tTotal) << endl;
   s << "DET: " << writeTimeWithPercent(ts.tTotalDetection, ts.tTotal) << endl;
   s << "EXP: " << writeTimeWithPercent(ts.tExposure, ts.tTotal) << endl;
   s << "MR:  " << writeTimeWithPercent(ts.tMeasureRough, ts.tTotal) << endl;
   s << "M:   " << writeTimeWithPercent(ts.tMeasure, ts.tTotal) << endl;

   s << "TOT: " << fixed << setprecision(3) << ts.tTotal * 1e3
   << endl;
   */

  return s;
}

}
