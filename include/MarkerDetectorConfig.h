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
 * MarkerDetectorConfig.h
 *
 *  Created on: Mar 20, 2015
 *      Author: Davide A. Cucci (davide.cucci@epfl.ch)
 */

#ifndef MARKERDETECTORCONFIG_H_
#define MARKERDETECTORCONFIG_H_

#include <iostream>
#include <opencv/cv.hpp>
#include <libconfig.h++>

namespace visiona {

struct MarkerDetectorConfig {

    int CannyBlurKernelSize;
    int CannyLowerThreshold;
    int CannyHigherThreshold;

    float contourFilterMinSize;

    float markerDiameter;
    float markerInnerDiameter;

    float markerSignalRadiusPercentage;

    float markerSignalStartsWith;
    std::vector<float> markerSignalModel; // percentage of the flips

    float markerxCorrThreshold;

    cv::Mat K;
    cv::Mat distortion;

    MarkerDetectorConfig();

    bool loadConfig(const std::string &file);

  protected:

    template <typename T>
    void fromSettingToVector(libconfig::Setting &s, std::vector<T> &v) {
      v.resize(s.getLength());

      for (int i = 0; i < s.getLength(); i++) {
        v[i] = s[i];
      }
    }
};

}

#endif /* MARKERDETECTORCONFIG_H_ */
