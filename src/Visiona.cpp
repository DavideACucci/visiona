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
 * Visiona.cpp
 *
 *  Created on: Apr 21, 2016
 *      Author: Davide A. Cucci (davide.cucci@epfl.ch)
 */

#include "Visiona.h"

#include <limits>

#include "MarkerDetector_impl.h"

namespace visiona {

Circle::Circle(const Contour& cnt_in, const cv::Point2f &center_in, float r_in) :
    cnt(cnt_in), center(center_in), r(r_in) {
}

Circle::Circle() {
}

Target::Target() :
    detected(false), roughlyMeasured(false), measured(false), white(255), black(
        0), meanReprojectionError(std::numeric_limits<float>::infinity()) {
}

MarkerDetector *MarkerDetectorFactory::makeMarkerDetector(
    const MarkerDetectorConfig &cfg) {
  return new MarkerDetector_impl(cfg);
}

}
