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
 * Visiona.h
 *
 *  Created on: Apr 21, 2016
 *      Author: Davide A. Cucci (davide.cucci@epfl.ch)
 */

#ifndef VISIONA_H_
#define VISIONA_H_

#include <vector>
#include <memory>

#include <opencv/cv.h>

#include "MarkerDetectorConfig.h"
#include "DebugPlotConfig.h"

namespace visiona {

typedef std::vector<cv::Point> Contour;

/**
 * Circle approximation of a contour
 */

struct Circle {
    Contour cnt;
    cv::Point2f center;
    float r;

    Circle(const Contour& cnt_in, const cv::Point2f &center_in, float r_in);
    Circle();
};

/**
 * A target
 * includes temporary and final results of the detect and measure process
 */

class Target {
  public:
    // detect stage
    bool detected; /**< if the target has been detected in the image **/
    Circle inner, outer; /**< inner and outer circle circle approximations **/
    float heading; /**< heading of the target on the image plane **/

    // evaluate exposure stage;
    float black, white; /**< estimated grey level for 'black' and 'white' pixels **/

    // measure-rough stage;
    bool roughlyMeasured; /**< if the rough measurement succeded **/
    std::vector<cv::Point2f> codePoints; /**< image measurements for the white dots used in PnP **/
    cv::Mat roughR, rought; /**< rotation and translation from camera to target **/
    double meanReprojectionError; /**< mean reprojection error of codePoints after PnP **/

    // measure stage
    bool measured; /**< if the accurate measurement step succeded **/

    // temporary, should become a rotation and a translation
    double cx, cy, scale;
    double phi, kappa;

    Target();
};

class MarkerDetector {
  protected:
    MarkerDetectorConfig _cfg;

  public:
    /**
     * \brief search the image for the target
     *
     * param image the input image, must be of type CV_8UC1
     * param dbg the debug configuration, default, NO
     *
     * return A vector containing the targets found (currently only one is returned)
     */
    virtual std::vector<std::shared_ptr<Target>> detect(const cv::Mat &image, DebugPlotConfig *dbg = NULL) = 0;

    /**
     * \brief given an already detected target, determines value for the black and white levels
     *
     * param tg the detected target
     * param dbg the debug configuration, default, NO
     */
    virtual void evaluateExposure(const cv::Mat &image,
        std::shared_ptr<Target> tg, DebugPlotConfig *dbg = NULL) = 0;

    /**
     * \brief roughly measures the target based on the center of the white circles forming the code
     *
     * param tg the detected target
     * param dbg the debug configuration, default, NO
     */
    virtual bool measureRough(const cv::Mat &image, std::shared_ptr<Target> tg, DebugPlotConfig *dbg = NULL) = 0;

    /**
     * \brief measures accurately the target based on the concentric circle algorithm
     *
     * param tg the detected target
     * param dbg the debug configuration, default, NO
     */
    virtual bool measure(const cv::Mat &image, std::shared_ptr<Target> tg, DebugPlotConfig *dbg = NULL) = 0;


    MarkerDetector(const MarkerDetectorConfig &cfg) :
        _cfg(cfg) {
    }

    virtual ~MarkerDetector() {
    }
};

class MarkerDetectorFactory {
  public:
    static MarkerDetector *makeMarkerDetector(const MarkerDetectorConfig &cfg);
};

}

#endif /* VISIONA_H_ */
