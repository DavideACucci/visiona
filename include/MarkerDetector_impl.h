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
 * MarkerDetector.h
 *
 *  Created on: Mar 20, 2015
 *      Author: Davide A. Cucci (davide.cucci@epfl.ch)
 */

#ifndef MARKERDETECTOR_H_
#define MARKERDETECTOR_H_

#include <iostream>
#include <fstream>

#include <Eigen/Dense>

#include "Visiona.h"

namespace visiona {

class MarkerDetector_impl: public MarkerDetector {

  public:

    typedef std::vector<Contour> Contours;

    typedef cv::RotatedRect Ellipse;
    typedef Eigen::Matrix<double, 6, 1> EllipsePoly;

    typedef std::vector<Circle> Circles;

    struct CircleCluster {
        cv::Point2f center;

        int rep; // the representative circle
        std::vector<int> circleIds;
    };

    cv::Mat _floodfillMask;
    std::vector<cv::Point3f> _worldPoints;

    MarkerDetector_impl(const MarkerDetectorConfig &cfg);

    // TODO: remove this wrapper method
    inline std::vector<std::shared_ptr<Target>> detect(const cv::Mat &raw,
        DebugPlotConfig *dbg = NULL) {
      std::shared_ptr<Target> tg(new Target);
      std::vector<std::shared_ptr<Target>> ret;
      ret.push_back(tg);

      tg->detected = detect(raw, tg->outer, tg->inner, tg->heading, dbg);

      return ret;
    }

    bool detect(const cv::Mat &raw, Circle &outer, Circle &inner,
        float &heading, DebugPlotConfig *dbg = NULL);

    bool measure(const cv::Mat &image, std::shared_ptr<Target> tg,
        DebugPlotConfig *dbg = NULL);

    inline bool measureRough(const cv::Mat &image, std::shared_ptr<Target> tg,
        std::vector<cv::Point2f> *seedPoints = NULL,
        DebugPlotConfig *dbg = NULL);

    // TODO: remove this wrapper method
    void evaluateExposure(const cv::Mat &image, std::shared_ptr<Target> tg,
        DebugPlotConfig *dbg = NULL) {
      evaluateExposure(image, tg->outer, tg->heading, tg->black, tg->white,
          dbg);
    }

    void evaluateExposure(const cv::Mat &raw, const Circle &outer,
        float heading, float &black, float &white, DebugPlotConfig *dbg = NULL);

    void detectEdges(const cv::Mat &raw, cv::Mat &edges, DebugPlotConfig *dbg =
    NULL);

    void parallelThreshold(const cv::Mat &raw, cv::Mat &edges, int nThreads,
        DebugPlotConfig *dbg = NULL);

    void detectContours(const cv::Mat &edges, Contours &out,
        DebugPlotConfig *dbg = NULL);

    void filterContours(const Contours &in, Circles &out, DebugPlotConfig *dbg =
    NULL);

    void clusterCircles(const Circles &in, std::vector<CircleCluster> &out,
        DebugPlotConfig *dbg = NULL);

    bool selectMarker(const cv::Mat &image, const Circles &in,
        const std::vector<int> &ids, int &mId, float &theta,
        DebugPlotConfig *dbg = NULL);

    bool findConcentricCircles(const Circles &circles,
        const CircleCluster &cluster, int &innerCircle, DebugPlotConfig *dbg =
        NULL);

    void fitEllipse(const Contour &cnt, Ellipse &e,
        DebugPlotConfig *dbg = NULL);

    void getEllipseMatrix(const Ellipse &elps, Eigen::Matrix3d &Q,
        DebugPlotConfig *dbg = NULL);

    void getEllipsePolynomialCoeff(const Ellipse &in, EllipsePoly &out);

    void getEllipseLineIntersections(const EllipsePoly &in, double x0,
        double y0, double theta, cv::Point2d &p1, cv::Point2d &p2);

    void getDistanceGivenCenter(const EllipsePoly& elps, const cv::Point2d& c,
        double r, double &mu, double &std, int N);

    void getPoseGivenCenter(const EllipsePoly& elps, const cv::Point2d& c,
        double r, double &d, double &phi, double &kappa, DebugPlotConfig *dbg =
        NULL);

    void getDistanceWithGradientDescent(const EllipsePoly &outer,
        const EllipsePoly &inner, const cv::Point2d x0, double step,
        double lambda, cv::Point2d &xend, double tolX, double tolFun);

    double evalDistanceF(const EllipsePoly &outer, const EllipsePoly &inner,
        const cv::Point2d &x, const cv::Point2d &x0);

    void getProjectedCircleCenter(const Circles &in, cv::Point2f &projCenter,
        DebugPlotConfig *dbg = NULL);

    void buildCircleHolesImage(const EllipsePoly &elps, const cv::Point2f &c,
        double d, const std::string &fname, int N, float dx);

  protected:
    void subpixelEdgeWithLeastSquares(const cv::Mat &image, const Ellipse &elps,
        const EllipsePoly &poly, float theta, float a, float b,
        cv::Point2f &edge, int N = 2, DebugPlotConfig *dbg =
        NULL);

    void refineEllipseCntWithSubpixelEdges(const cv::Mat &image,
        const Target &tg, const Ellipse &elps,
        bool ignoreSignalAreas, int N, std::vector<cv::Point2f> &cnt,
        std::vector<double> &angles,
        DebugPlotConfig *dbg = NULL);

    void getSignalInsideCircle(const cv::Mat &image, const Circle &circle,
        float radiusPercentage, std::vector<float> &signal);

    void getSignalInsideEllipse(const cv::Mat &image, const Ellipse &ellipse,
        float radiusPercentage, std::vector<float> &signal, float &increment,
        float thetasmall = -M_PI, float thetabig = M_PI,
        std::vector<cv::Point2f> *pts = NULL);

    void normalizeSignal(std::vector<float> &signal);

    void computeNormalizedxCorr(const std::vector<float> &signal, cv::Mat &out);

    cv::Point2f evalEllipse(float at, const cv::Point2f &c, float a, float b,
        float phi);

    void initColorPlotSurfacte(cv::Mat &in, const DebugPlotConfig *dbg) const;

    void initZoomedSubregionSurface(cv::Point2i center, float r, cv::Mat &in,
        cv::Mat &out, int size, float &ratio, cv::Point2i &basept);

    inline cv::Point2i transformPoint(const cv::Point2f &in,
        const cv::Point2i &basept, float ratio) {
      return (in + cv::Point2f(-basept.x + 0.5, -basept.y + 0.5)) * ratio;
    }

    inline Ellipse transofrmEllipse(const Ellipse &in,
        cv::Point2i &basept, float ratio) {
      Ellipse out;
      out.angle = in.angle;
      out.size = in.size * ratio;
      out.center = transformPoint(in.center, basept, ratio);

      return out;
    }

    inline float getSubpix(const cv::Mat &img, const cv::Point2f &pt) {
      cv::Mat patch;
      cv::getRectSubPix(img, cv::Size(1, 1), pt, patch, CV_32FC1);
      return patch.at<float>(0, 0);
    }

    cv::Point2f distort(const cv::Point2f &p);

    void drawCube(cv::Mat &image, const cv::Mat &R, const cv::Mat &t,
        const float l, const cv::Scalar &color = cv::Scalar(0, 0, 255),
        const cv::Point2i &basept = cv::Point2i(0, 0), float ratio = 1.0);

    void drawPyramid(cv::Mat &image, const cv::Mat &R, const cv::Mat &t,
        const float l, const cv::Scalar &color = cv::Scalar(0, 0, 255));

    inline void wrapTo2PiInPlace(float &a) {
      bool was_neg = a < 0;
      a = std::fmod(a, 2.0 * M_PI);
      if (was_neg) {
        a += 2.0 * M_PI;
      }
    }

    inline void wrapToSize(unsigned int &towrap, int size) {
      towrap = (towrap + size) % size;
    }

    template<class T>
    inline void printVector(const std::vector<T> &in) {
      for (unsigned int i = 0; i < in.size(); ++i) {
        std::cerr << in[i] << " ";
      }
      std::cerr << std::endl;
    }

    template<class T>
    bool outputContour(const std::vector<T> &cnt, const std::string &path,
        const std::string &prefix, long int frameNumber,
        const std::string &suffix) {

      std::stringstream fname;
      fname << path << "/" << prefix << frameNumber << suffix;

      std::ofstream f(fname.str(), std::ios_base::out);

      if (!f.is_open()) {
        return false;
      }

      f << std::fixed;
      f.precision(12);

      for_each(cnt.begin(), cnt.end(),
          [&f](const T &p) {f << p.x << " " << p.y << std::endl;});

      f.close();

      return true;
    }

    template<class T>
    bool outputContour(const std::vector<T> &cnt,
        const std::vector<double> &angles, const std::string &path,
        const std::string &prefix, long int frameNumber,
        const std::string &suffix) {

      assert(cnt.size() == angles.size());

      std::stringstream fname;
      fname << path << "/" << prefix << frameNumber << suffix;

      std::ofstream f(fname.str(), std::ios_base::out);

      if (!f.is_open()) {
        return false;
      }

      f << std::fixed;
      f.precision(12);

      for (int i = 0; i < cnt.size(); ++i) {
        f << cnt[i].x << " " << cnt[i].y << " " << angles[i] << std::endl;
      }

      f.close();

      return true;
    }
};

}

#endif /* MARKERDETECTOR_H_ */
