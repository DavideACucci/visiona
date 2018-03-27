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
 * DebugPlotConfig.h
 *
 *  Created on: Mar 20, 2015
 *      Author: Davide A. Cucci (davide.cucci@epfl.ch)
 */

#ifndef DEBUGPLOTCONFIG_H_
#define DEBUGPLOTCONFIG_H_

#include <iomanip>
#include <string>
#include <opencv/cv.hpp>

namespace visiona {

struct TimeStats {
    double tEdgeDetection;
    double tContours;
    double tFilteredContours;
    double tClusterCircles;
    double tSelectMarker;
    double tTotalDetection;

    double tExposure;

    double tMeasureRough;

    double tMeasure;

    double tTotal;

    TimeStats();
};

std::ostream &operator<<(std::ostream &s, const TimeStats &ts);

struct DebugPlotConfig {

  DebugPlotConfig();

  /**
   * enable specific debug plot
   */
  bool enabled;
  bool enableEdgeDetection; // plots the raw edge pixels
  bool enableContours; // plots the clustered contours
  bool enableFilteredContours; // plots the filtered contours
  bool enableCirclesClusters; // plots the couples of concentric circles
  bool enableSelectedTargets; // plots the detected target or where the signal was searched in case of failure
  bool enableEsposure; // plots pixels sampled for exposure determination
  bool enableRoughMeasure; // plot a cube on the target in case of success, seed points in case of failure
  bool enableSubPixelEllipses; // plots the subpixel ellipses locations

  /**
   * the name of the window in which image should be plot
   */
  std::string windowName;

  /**
   * target debug plots can plot a zoomed regiond around the target
   */
  bool blitSubRegion;

  /**
   * 1 selects just the target, extra area can be plot increasing it
   */
  double blitRegionWidthMultiplier;

  /**
   * input image
   */
  cv::Mat rawImage;

  /**
   * temporary image on which the debug plot is performed
   */
  cv::Mat dbgImage;

  /**
   * path in which debug output files should be saved
   */
  std::string debugFilesPath;

  /**
   * image counter, useful for debug file name generation
   */
  long int frameNumber;

  /**
   * if ellipse contours coordinates file should be produced
   */
  bool writeContours;
  bool writeSubpixelContours;

  TimeStats tstats;
};

}

#endif /* DEBUGPLOTCONFIG_H_ */
