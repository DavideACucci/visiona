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
 * MarkerDetector.cpp
 *
 *  Created on: Mar 20, 2015
 *      Author: Davide A. Cucci (davide.cucci@epfl.ch)
 */

#include "MarkerDetector_impl.h"
#include "Timer.h"

#include <thread>

#include <Eigen/Eigenvalues>

#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/calib3d/calib3d.hpp>

using namespace cv;
using namespace std;

namespace visiona {

MarkerDetector_impl::MarkerDetector_impl(const MarkerDetectorConfig &cfg) :
    MarkerDetector(cfg) {

  for (int cnt = 0; cnt < _cfg.markerSignalModel.size() / 2; ++cnt) {
    int i = (_cfg.markerSignalStartsWith == 1.0 ? 0 : 1) + 2 * cnt;

    float maxAngle, minAngle, angle;

    if (i == 0) {
      minAngle = 2 * M_PI
          * (_cfg.markerSignalModel[_cfg.markerSignalModel.size() - 1] - 1);
    } else {
      minAngle = 2 * M_PI * _cfg.markerSignalModel[i - 1];
    }

    maxAngle = 2 * M_PI * _cfg.markerSignalModel[i];

    angle = 0.5 * (maxAngle + minAngle);

    Point3f wp;
    wp.x = cos(angle);
    wp.y = sin(angle);
    wp.z = 0.0;
    wp *= _cfg.markerDiameter * _cfg.markerSignalRadiusPercentage / 2.0;

    _worldPoints.push_back(wp);
  }

}

bool MarkerDetector_impl::detect(const cv::Mat &raw, Circle &outer,
    Circle &inner, float &heading, DebugPlotConfig *dbg) {

  if (dbg != NULL) {
    tic();
  }

  Mat edges;
  detectEdges(raw, edges, dbg);

  Contours ctr;
  detectContours(edges, ctr, dbg);

  Circles circles;
  filterContours(ctr, circles, dbg);

  std::vector<CircleCluster> clusters;
  clusterCircles(circles, clusters, dbg);

  bool found = false;

  if (circles.size() > 0) {
    // Select marker among the cluster representatives
    std::vector<int> representatives;
    for (auto it = clusters.begin(); it != clusters.end(); ++it) {
      representatives.push_back(it->rep);
    }

    int selectedCluster;
    found = selectMarker(raw, circles, representatives, selectedCluster,
        heading, dbg);

    if (found) {

      // in which we have exactly two circles per cluster
      assert(clusters[selectedCluster].circleIds.size() == 2);

      outer = circles[clusters[selectedCluster].circleIds[0]];
      inner = circles[clusters[selectedCluster].circleIds[1]];
      //*/

      /* in which in each cluster there can be many circles
       int selectedOuterId = clusters[selectedCluster].rep;
       outer = circles[selectedOuterId];

       // for the selected cluster find the concentric circle
       int selectedInnerId;
       bool innerCircleFound = findConcentricCircles(circles,
       clusters[selectedCluster], selectedInnerId);

       if (innerCircleFound) {
       inner = circles[selectedInnerId];
       return true;
       }
       return false;
       */
    }
  }

  if (dbg != NULL) {
    dbg->tstats.tTotalDetection = toc();
  }

  return found;
}

bool MarkerDetector_impl::measure(const cv::Mat &image,
    std::shared_ptr<Target> tg,
    DebugPlotConfig *dbg) {

//bool MarkerDetector_impl::measure(const Mat &image, const Circle& outer,
//    const Circle &inner, float heading, double &cx, double &cy, double &scale,
//    double &phi, double &kappa, float black, float white,
//    DebugPlotConfig *dbg) {

  if (dbg != NULL) {
    tic();
  }

  // fit ellipses on the markers
  Ellipse outerElps, innerElps;
  fitEllipse(tg->outer.cnt, outerElps);
  fitEllipse(tg->inner.cnt, innerElps);

  // refine with subpixel
  vector<Point2f> outerElpsCntSubpx;
  vector<double> outerAngles;
  refineEllipseCntWithSubpixelEdges(image, *tg, outerElps, true, 2,
      outerElpsCntSubpx, outerAngles);

  vector<Point2f> innerElpsCntSubpx;
  vector<double> innerAngles;
  refineEllipseCntWithSubpixelEdges(image, *tg, innerElps, true, 2,
      innerElpsCntSubpx, innerAngles);

  // undistort ellipse points
  vector<Point2f> outerElpsCntSubpx_und;
  undistortPoints(outerElpsCntSubpx, outerElpsCntSubpx_und, _cfg.K,
      _cfg.distortion, Mat(), _cfg.K);
  vector<Point2f> innerElpsCntSubpx_und;
  undistortPoints(innerElpsCntSubpx, innerElpsCntSubpx_und, _cfg.K,
      _cfg.distortion, Mat(), _cfg.K);
  //*/

  /* or not
   vector<Point2f> outerElpsCntSubpx_und = outerElpsCntSubpx;
   vector<Point2f> innerElpsCntSubpx_und = innerElpsCntSubpx;
   //*/

  // fit the ellipse
  Ellipse outerElpsSubpx = cv::fitEllipse(outerElpsCntSubpx_und);
  outerElps = outerElpsSubpx;
  Ellipse innerElpsSubpx = cv::fitEllipse(innerElpsCntSubpx_und);
  innerElps = innerElpsSubpx;

  // get ellipse polynomials
  EllipsePoly outerPoly, innerPoly;
  getEllipsePolynomialCoeff(outerElpsSubpx, outerPoly);
  getEllipsePolynomialCoeff(innerElpsSubpx, innerPoly);

  // compute center
  Point2d center;
  getDistanceWithGradientDescent(outerPoly, innerPoly, outerElps.center, 1e-6,
      1e-2, center, 1e-8, 0);

  // compute distance
  getPoseGivenCenter(outerPoly, center, _cfg.markerDiameter / 2.0, tg->distance,
      tg->phi, tg->kappa, dbg);

  // get the center in the distorted image
  center = distort(center);
  tg->cx = center.x;
  tg->cy = center.y;

  tg->measured = true;

  if (dbg != NULL) {

    dbg->tstats.tMeasure = toc();

    if (dbg->enabled && dbg->enableSubPixelEllipses) {
      float ratio;
      Point2i basept;

      initZoomedSubregionSurface(tg->outer.center,
          tg->outer.r * dbg->blitRegionWidthMultiplier,
          dbg->rawImage, dbg->dbgImage, 1280, ratio, basept);

//    circle(dbg->dbgImage, transformPoint(Point2f(cx, cy), basept, ratio),
//        0.5 * ratio, Scalar(255, 0, 0), -1);
//    ellipse(dbg->dbgImage, transofrmEllipse(outerElpsSubpx, basept, ratio),
//        Scalar(255, 0, 0), 0.25 * ratio);
//    ellipse(dbg->dbgImage, transofrmEllipse(innerElpsSubpx, basept, ratio),
//        Scalar(255, 0, 0), 0.25 * ratio);

      for (int i = 0; i < outerElpsCntSubpx.size(); ++i) {
        circle(dbg->dbgImage,
            transformPoint(outerElpsCntSubpx[i], basept, ratio),
            0.25 * ratio, Scalar(0, 255, 0), -1);
      }

      for (int i = 0; i < innerElpsCntSubpx.size(); ++i) {
        circle(dbg->dbgImage,
            transformPoint(innerElpsCntSubpx[i], basept, ratio),
            0.25 * ratio, Scalar(0, 0, 255), -1);
      }

      imshow(dbg->windowName, dbg->dbgImage);
      waitKey(1);

//      printf(" Center: (%.2f, %.2f), Distance: %.5f, Phi: %.5f, Kappa: %.5f\n",
//          tg->cx, tg->cy, tg->scale, tg->phi / M_PI * 180.0, tg->kappa / M_PI * 180.0);

      if (dbg->writeSubpixelContours) {
        outputContour(outerElpsCntSubpx, outerAngles, dbg->debugFilesPath,
            "contour_sbpx_",
            dbg->frameNumber, "_out.txt");
        outputContour(innerElpsCntSubpx, innerAngles, dbg->debugFilesPath,
            "contour_sbpx_",
            dbg->frameNumber, "_in.txt");
      }
    }
  }

  return true;
}

void MarkerDetector_impl::evaluateExposure(const cv::Mat &raw,
    const Circle &outer, float heading, float &black, float &white,
    DebugPlotConfig *dbg) {

  if (dbg != NULL) {
    tic();
  }

  Ellipse outerElps;
  fitEllipse(outer.cnt, outerElps);

  float outMargin = 10;
  float inMargin = 2;

  float outerElpsSize = (outerElps.size.width + outerElps.size.height) / 4.0;

  const float outScale = (outerElpsSize + outMargin) / outerElpsSize;
  const float inScale = (_cfg.markerDiameter + _cfg.markerInnerDiameter) / 2.0
      / _cfg.markerDiameter;

  // TODO: should be at the correct percentage
  float safetyAngle = atan(inMargin / (outerElpsSize * inScale)); // to avoid getting the white dots

  float inc;
  vector<float> signal;

  getSignalInsideEllipse(raw, outerElps, outScale, signal, inc);

  white = 0;
  for (const float &d : signal) {
    white += d;
  }
  white /= signal.size();

  signal.clear();

  // TODO: check minimum and maximun angles, they are somewhat wrong
  vector<Point2f> pts;
  getSignalInsideEllipse(raw, outerElps, inScale, signal, inc,
      heading + 2.0 * M_PI * _cfg.markerSignalModel[4] + safetyAngle,
      heading + 2.0 * M_PI * _cfg.markerSignalModel[5] - safetyAngle, &pts);

  black = 0;
  for (const float &d : signal) {
    black += d;
  }
  black /= signal.size();

  if (dbg != NULL) {

    dbg->tstats.tExposure = toc();

    if (dbg->enabled && dbg->enableEsposure) {
      cerr << "black: " << black << " white: " << white << endl;

      initColorPlotSurfacte(dbg->dbgImage, dbg);

      Ellipse tmpOutElps = outerElps;
      tmpOutElps.size.height *= outScale;
      tmpOutElps.size.width *= outScale;

      ellipse(dbg->dbgImage, tmpOutElps, Scalar(0, 0, 255));

//    Ellipse tmpInElps = outerElps;
//    tmpInElps.size.height *= inScale;
//    tmpInElps.size.width *= inScale;
//
//    ellipse(dbg->dbgImage, tmpInElps, Scalar(0, 255, 0));

      for_each(pts.begin(), pts.end(), [&dbg](const Point2f &p) {
        circle(dbg->dbgImage, p, 1, Scalar(0,255,0),-1);
      });

      imshow(dbg->windowName, dbg->dbgImage);
      waitKey(1);
    }
  }
}

void MarkerDetector_impl::detectEdges(const Mat& raw, Mat& edges,
    DebugPlotConfig *dbg) {

  if (dbg != NULL) {
    tic();
  }

  Mat tmp;

  // with canny
  if (_cfg.CannyBlurKernelSize > 0) {
    blur(raw, tmp, Size(_cfg.CannyBlurKernelSize, _cfg.CannyBlurKernelSize));

    Canny(tmp, edges, _cfg.CannyLowerThreshold, _cfg.CannyHigherThreshold, 3,
        true);
  } else {
    Canny(raw, edges, _cfg.CannyLowerThreshold, _cfg.CannyHigherThreshold, 3,
        true);
  }
  //*/

  /* with thresholding (single thread)
   adaptiveThreshold(raw, edges, 255, ADAPTIVE_THRESH_MEAN_C,
   THRESH_BINARY_INV, 15, 7);
   //*/

  /* parallel version, multithread
   parallelThreshold(raw, edges, 4);
   //*/

  /* debug parallel threasholding
   Mat test;
   tStart = clock();
   adaptiveThreshold(raw, test, 255, ADAPTIVE_THRESH_MEAN_C,
   THRESH_BINARY_INV, 7, 7);
   printf("Single: %.6fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);

   double minVal, maxVal;
   Point minLoc, maxLoc;
   minMaxLoc(test-edges, &minVal, &maxVal, &minLoc, &maxLoc );

   cout << "min val : " << minVal << endl;
   cout << "max val: " << maxVal << endl;

   //*/

  if (dbg != NULL) {

    dbg->tstats.tEdgeDetection = toc();

    if (dbg->enabled && dbg->enableEdgeDetection) {

      initColorPlotSurfacte(dbg->dbgImage, dbg);

      for (int y = 0; y < dbg->dbgImage.rows; ++y) {
        for (int x = 0; x < dbg->dbgImage.cols; ++x) {
          if (edges.at<uchar>(y, x)) {
            dbg->dbgImage.at<Vec3b>(y, x).val[2] = 255;
          }
        }
      }

      imshow(dbg->windowName, dbg->dbgImage);
      waitKey(1);
    }
  }
}

void MarkerDetector_impl::parallelThreshold(const cv::Mat& raw, cv::Mat& edges,
    int nThreads, DebugPlotConfig* dbg)
    {

  edges.create(raw.rows, raw.cols, CV_8UC1);

  // divide the image in N horizontal slices and process them
  int sliceH = floor((float) raw.rows / nThreads);

  thread t[nThreads];

  for (int i = 0; i < nThreads; ++i) {
    Rect roi(0, i * sliceH, raw.cols, sliceH);

    t[i] = thread(adaptiveThreshold, raw(roi), edges(roi), 255,
        ADAPTIVE_THRESH_MEAN_C,
        THRESH_BINARY_INV, 7, 7);
  }

  // in the meanwhile compute the correct junction bands

  Mat tmpM[nThreads - 1];

  for (int i = 0; i < nThreads - 1; ++i) {
    Rect roi(0, (i + 1) * sliceH, raw.cols, 12);
    tmpM[i].create(12, raw.cols, CV_8UC1);

    adaptiveThreshold(raw(roi), tmpM[i], 255, ADAPTIVE_THRESH_MEAN_C,
        THRESH_BINARY_INV, 7, 7);
  }

  // wait for threads to finish
  for (int i = 0; i < nThreads; ++i) {
    t[i].join();
  }

  // overwrite
  for (int i = 0; i < nThreads - 1; ++i) {
    Rect roi(0, (i + 1) * sliceH, raw.cols, 12);
    edges(roi) = tmpM[i];
  }
}

void MarkerDetector_impl::detectContours(const Mat &edges, Contours &ctrs,
    DebugPlotConfig *dbg) {

  if (dbg != NULL) {
    tic();
  }

  ctrs.clear();

  findContours(edges, ctrs, CV_RETR_LIST, CV_CHAIN_APPROX_NONE, Point(0, 0));

  if (dbg != NULL) {

    dbg->tstats.tContours = toc();

    if (dbg->enabled && dbg->enableContours) {

      initColorPlotSurfacte(dbg->dbgImage, dbg);

      for (unsigned int i = 0; i < ctrs.size(); i++) {
        Scalar color(0, 0, 255);

        drawContours(dbg->dbgImage, ctrs, i, color, 1, 8);
      }
      imshow(dbg->windowName, dbg->dbgImage);
      waitKey(1);
    }
  }
}

void MarkerDetector_impl::filterContours(const Contours& in, Circles& out,
    DebugPlotConfig *dbg) {

  if (dbg != NULL) {
    tic();
  }

  out.clear();

  for (auto cit = in.begin(); cit != in.end(); ++cit) {
    const Contour &c = *cit;

    if (c.size() < _cfg.contourFilterMinSize) {
      continue;
    }

    // compute center
    float cx = 0, cy = 0;
    for (auto it = c.begin(); it != c.end(); ++it) {
      cx += it->x;
      cy += it->y;
    }

    cx /= c.size();
    cy /= c.size();

    // compute average and std of distance

    float sumd = 0, sumd2 = 0, d2;
    for (auto it = c.begin(); it != c.end(); ++it) {
      d2 = pow(it->x - cx, 2) + pow(it->y - cy, 2);
      sumd2 += d2;
      sumd += sqrt(d2);
    }

    float rMean = sumd / c.size();
    float rStd = sqrt(
        c.size() * (sumd2 / c.size() - pow(rMean, 2)) / (c.size() - 1));

    // Attempt to have a metric that adapts to depth
    if ((rStd / rMean) < 0.075) {
      Point2f center(cx, cy);
      out.emplace_back(c, center, rMean);
    }

  }

  if (dbg != NULL) {

    dbg->tstats.tFilteredContours = toc();

    if (dbg->enabled && dbg->enableFilteredContours) {
      initColorPlotSurfacte(dbg->dbgImage, dbg);

      for (unsigned int i = 0; i < out.size(); i++) {
        Scalar color(0, 255, 0);

        Contours tmp;
        tmp.push_back(out[i].cnt);

        drawContours(dbg->dbgImage, tmp, 0, color, 2, 2);
      }
      imshow(dbg->windowName, dbg->dbgImage);
      waitKey(1);
    }
  }
}

void MarkerDetector_impl::clusterCircles(const Circles &in,
    std::vector<CircleCluster> &out, DebugPlotConfig *dbg) {

  if (dbg != NULL) {
    tic();
  }

  // for each circle, do a cluster with the best matching inner circle, if found

  float radiusRatio = _cfg.markerInnerDiameter / _cfg.markerDiameter;

  for (int i = 0; i < in.size(); ++i) {

    int bestMatch = -1;
    float bestDiff = numeric_limits<float>::infinity();

    for (int j = 0; j < in.size(); ++j) {
      if (in[i].r > in[j].r
          && norm(in[i].center - in[j].center) < in[i].r * 0.5) {
        float curDiff = fabs(in[i].r * radiusRatio - in[j].r);

        if (curDiff < bestDiff) {
          bestMatch = j;
          bestDiff = curDiff;
        }
      }
    }

    if (bestMatch != -1) {
      if (fabs(in[bestMatch].r / in[i].r / radiusRatio - 1) < 0.25) {
        out.resize(out.size() + 1);
        out.back().center = in[i].center;
        out.back().rep = i;
        out.back().circleIds.push_back(i);
        out.back().circleIds.push_back(bestMatch);
      }
    }
  }

  if (dbg != NULL) {

    dbg->tstats.tClusterCircles = toc();

    if (dbg->enabled && dbg->enableCirclesClusters) {

      initColorPlotSurfacte(dbg->dbgImage, dbg);

      std::vector<Scalar> colors;
      colors.push_back(Scalar(255, 0, 0));
      colors.push_back(Scalar(128, 128, 0));
      colors.push_back(Scalar(0, 255, 0));
      colors.push_back(Scalar(0, 128, 128));
      colors.push_back(Scalar(0, 0, 255));
      colors.push_back(Scalar(0, 128, 255));
      colors.push_back(Scalar(0, 255, 255));
      colors.push_back(Scalar(128, 128, 255));
      colors.push_back(Scalar(255, 0, 255));
      colors.push_back(Scalar(255, 128, 128));
      colors.push_back(Scalar(255, 255, 0));

      for (unsigned int i = 0; i < out.size(); i++) {

        Scalar &color = colors[i % colors.size()];

        for (int j = 0; j < out[i].circleIds.size(); ++j) {

          Contours tmpcnt;
          tmpcnt.push_back(in[out[i].circleIds[j]].cnt);

          drawContours(dbg->dbgImage, tmpcnt, 0, color);
        }
      }
      imshow(dbg->windowName, dbg->dbgImage);
      waitKey(1);
    }
  }
}

bool MarkerDetector_impl::selectMarker(const Mat& image,
    const Circles &candidates, const vector<int> &representativesIds,
    int &selectedCluster, float & theta, DebugPlotConfig *dbg) {

  if (dbg != NULL) {
    tic();
  }

  bool found = false;

  float maxCorr = _cfg.markerxCorrThreshold;

  for (auto it = 0; it < representativesIds.size(); ++it) {

    const Circle &c = candidates[representativesIds[it]];

    vector<float> signal;
    getSignalInsideCircle(image, c, _cfg.markerSignalRadiusPercentage, signal);

    normalizeSignal(signal);

    Mat corr;
    computeNormalizedxCorr(signal, corr);

    double m, M;
    Point2i mLoc, MLoc;

    minMaxLoc(corr, &m, &M, &mLoc, &MLoc);

    if (M > maxCorr) {
      maxCorr = M;
      // compute theta
//      theta = -1.0 / c.r * MLoc.x;
      theta = +M_PI - (float) MLoc.x / signal.size() * 2.0 * M_PI;

      found = true;
      selectedCluster = it;
    }

//    cerr << M << endl;
  }

  if (dbg != NULL) {
    dbg->tstats.tSelectMarker = toc();

    if (dbg->enabled && dbg->enableSelectedTargets) {
      if (found) {
        // copy it so I can modify
        const Circle marker = candidates[representativesIds[selectedCluster]];

        Ellipse e;
        fitEllipse(marker.cnt, e);

        stringstream s;
        s << fixed << setprecision(2) << "w: " << e.size.width << " h: "
            << e.size.height << endl;

        float ratio = 1.0;
        Point2i basept(0, 0);

        if (dbg->blitSubRegion == true) {
          initZoomedSubregionSurface(marker.center,
              marker.r * dbg->blitRegionWidthMultiplier, dbg->rawImage,
              dbg->dbgImage, 1280, ratio, basept);
        } else {
          initColorPlotSurfacte(dbg->dbgImage, dbg);
        }

        Scalar color(0, 0, 255);

        e.center.x = (e.center.x - basept.x + 0.5) * ratio;
        e.center.y = (e.center.y - basept.y + 0.5) * ratio;
        e.size.width *= ratio;
        e.size.height *= ratio;

        ellipse(dbg->dbgImage, e, color, 3);

        Point2f ls, le;

        ls.x = e.center.x + cos(theta) * marker.r * 2.0 * ratio;
        ls.y = e.center.y + sin(theta) * marker.r * 2.0 * ratio;

        le.x = e.center.x - cos(theta) * marker.r * 2.0 * ratio;
        le.y = e.center.y - sin(theta) * marker.r * 2.0 * ratio;

        line(dbg->dbgImage, ls, le, color, 3);

        putText(dbg->dbgImage, s.str(), Point(30, 30), FONT_HERSHEY_SIMPLEX, 1,
            Scalar(0, 0, 255), 3);

      } else {
        initColorPlotSurfacte(dbg->dbgImage, dbg);

        std::vector<Scalar> colors;
        colors.push_back(Scalar(255, 0, 0));
        colors.push_back(Scalar(0, 255, 0));
        colors.push_back(Scalar(0, 0, 255));
        colors.push_back(Scalar(0, 255, 255));
        colors.push_back(Scalar(255, 0, 255));
        colors.push_back(Scalar(255, 255, 0));

        for (auto it = 0; it < representativesIds.size() && it < 6; ++it) {
          const Circle th = candidates[representativesIds[it]];

          Ellipse e;
          fitEllipse(th.cnt, e);

          e.size.width *= _cfg.markerSignalRadiusPercentage;
          e.size.height *= _cfg.markerSignalRadiusPercentage;

          ellipse(dbg->dbgImage, e, colors[it], 3);
        }
      }
      imshow(dbg->windowName, dbg->dbgImage);
      waitKey(1);
    }
  }

  return found;
}

bool MarkerDetector_impl::findConcentricCircles(const Circles &circles,
    const CircleCluster &cluster, int &innerCircle, DebugPlotConfig *dbg) {

  float repRadius = circles[cluster.rep].r;
  float searchedRadius = repRadius * _cfg.markerInnerDiameter
      / _cfg.markerDiameter;

  innerCircle = 0;
  float diff = numeric_limits<float>::infinity();

  // search the circle in the cluster for which the radius is nearest wrt searchedRadius
  for (auto it = cluster.circleIds.begin(); it != cluster.circleIds.end();
      ++it) {

    float d = fabs(circles[*it].r - searchedRadius);
    if (d < diff) {
      innerCircle = *it;
      diff = d;
    }
  }

  bool found = true;
  if (fabs((circles[innerCircle].r - searchedRadius) / searchedRadius) > 0.25) {
    cerr << "WARNING, inner circle possibly wrong" << endl;
    found = false;
  }

  if (dbg != NULL) {

    Scalar color(0, 0, 255);

    const Circle &outer = circles[cluster.rep];
    const Circle &inner = circles[innerCircle];

    if (dbg->writeContours) {
      outputContour(outer.cnt, dbg->debugFilesPath, "contour_",
          dbg->frameNumber, "_out.txt");
      outputContour(inner.cnt, dbg->debugFilesPath, "contour_",
          dbg->frameNumber, "_in.txt");
    }

    if (dbg->blitSubRegion) {

      float ratio;
      Point2i basept;

      initZoomedSubregionSurface(outer.center,
          outer.r * dbg->blitRegionWidthMultiplier, dbg->rawImage,
          dbg->dbgImage, 1280, ratio, basept);

      auto transformfun = [&basept, &ratio](Point r) {
        return (r - basept)*ratio + Point2i(0.5*ratio,0.5*ratio);
      };

      Contour outerCntAdapted(outer.cnt.size()), innerCntAdapted(
          inner.cnt.size());

      std::transform(outer.cnt.begin(), outer.cnt.end(),
          outerCntAdapted.begin(), transformfun);
      std::transform(inner.cnt.begin(), inner.cnt.end(),
          innerCntAdapted.begin(), transformfun);

      Contours cnts;
      cnts.emplace_back(outerCntAdapted);
      cnts.emplace_back(innerCntAdapted);

      drawContours(dbg->dbgImage, cnts, -1, Scalar(0, 0, 255), 5);

    } else {
      initColorPlotSurfacte(dbg->dbgImage, dbg);

      Contours cnts;
      cnts.push_back(outer.cnt);
      cnts.push_back(inner.cnt);

      drawContours(dbg->dbgImage, cnts, -1, color, 2, 2);
    }

    imshow(dbg->windowName, dbg->dbgImage);
    waitKey(1);
  }

  return found;
}

void MarkerDetector_impl::fitEllipse(const Contour& cnt, Ellipse& out,
    DebugPlotConfig *dbg) {

  out = cv::fitEllipse(cnt);

}

void MarkerDetector_impl::normalizeSignal(std::vector<float>& sig_in) {

// depolarize

  Mat sig(1, sig_in.size(), CV_32FC1, sig_in.data());

  float u = mean(sig)[0];

  sig -= u;

  float M = -std::numeric_limits<float>::infinity();
  float m = std::numeric_limits<float>::infinity();

  for (unsigned int k = 0; k < sig_in.size(); k++) {
    if (sig_in[k] > M) {
      M = sig_in[k];
    } else if (sig_in[k] < m) {
      m = sig_in[k];
    }
  }

  if (m != M) {
    for (unsigned int k = 0; k < sig_in.size(); k++) {
      sig_in[k] = -1.0 + (sig_in[k] - m) * 2.0 / (M - m);
    }
  }
}

void MarkerDetector_impl::computeNormalizedxCorr(
    const std::vector<float>& sig_in, Mat &out) {

// prepare signals

  float raw[2 * sig_in.size()];

  for (unsigned int n = 0; n < 2; n++) {

    unsigned int segment = 0;
    float val = _cfg.markerSignalStartsWith;

    for (unsigned int k = 0; k < sig_in.size(); k++) {
      if (segment < _cfg.markerSignalModel.size()) {
        // test if I have to advance
        if (k > _cfg.markerSignalModel[segment] * sig_in.size()) {
          segment++;
          val = -val;
        }
      }

      raw[k + n * sig_in.size()] = val;
    }
  }

  Mat ref(1, 2 * sig_in.size(), CV_32FC1, raw);

  Mat sig(1, sig_in.size(), CV_32FC1, const_cast<float *>(sig_in.data()));

  /*
   cerr << "-------------" << endl;
   printVector(sig_in);
   for (int k = 0; k < 2 * sig_in.size(); k++) {
   cerr << raw[k] << ", ";
   }
   cerr << endl;
   //*/

// compute cross correlation
  matchTemplate(ref, sig, out, CV_TM_CCORR_NORMED);
}

void MarkerDetector_impl::getSignalInsideCircle(const Mat& image,
    const Circle& circle, float radiusPercentage, std::vector<float>& signal) {

  int N = ceil(2 * M_PI / (1.0 / circle.r));
//  int N = 360;

  signal.clear();
  signal.reserve(N);
//  signal.reserve(180);

  Mat px;
  float theta = -M_PI;

//  while (theta < M_PI) {
  for (int i = 0; i < N; ++i) {
    unsigned int x = round(
        circle.center.x + cos(theta) * circle.r * radiusPercentage);
    unsigned int y = round(
        circle.center.y + sin(theta) * circle.r * radiusPercentage);

    Scalar intensity = image.at<uchar>(y, x);
    signal.push_back(intensity[0]);

//    theta += 1.0 / circle.r;
//    theta += M_PI/90.0;

    theta += (2.0 * M_PI) / N;
  }
}

void MarkerDetector_impl::getSignalInsideEllipse(const Mat& image,
    const Ellipse& ellipse, float radiusPercentage, std::vector<float>& signal,
    float &increment, float thetasmall, float thetabig, vector<Point2f> *pts) {

  increment = 1.0
      / ((ellipse.size.width / 2.0 + ellipse.size.height / 2.0) / 2.0);

  int N = ceil((thetabig - thetasmall) / increment);

  signal.clear();
  signal.reserve(N);

  if (pts != NULL) {
    pts->clear();
    pts->reserve(N);
  }

  Mat px;

  float theta = thetasmall;

  for (int i = 0; i < N; ++i, theta += increment) {
    Point2f px = evalEllipse(theta, ellipse.center,
        ellipse.size.width * radiusPercentage / 2.0,
        ellipse.size.height * radiusPercentage / 2.0,
        ellipse.angle * M_PI / 180.0);

    Scalar intensity = image.at<uchar>(px.y, px.x);
    signal.push_back(intensity[0]);

    if (pts != NULL) {
      pts->push_back(px);
    }
  }
}

Point2f MarkerDetector_impl::evalEllipse(float at, const Point2f& c, float a,
    float b, float phi) {

  float offset = 0.0;
  if (a < b) {
    offset = M_PI / 2.0;
    std::swap(a, b);
  }

  Point2f ret;

  ret.x = c.x + a * cos(at - phi + offset - M_PI) * cos(phi + offset)
      - b * sin(at - phi + offset - M_PI) * sin(phi + offset);
  ret.y = c.y + a * cos(at - phi + offset - M_PI) * sin(phi + offset)
      + b * sin(at - phi + offset - M_PI) * cos(phi + offset);

  return ret;
}

void MarkerDetector_impl::initColorPlotSurfacte(Mat& in,
    const DebugPlotConfig *dbg) const {
  cvtColor(dbg->rawImage, in, CV_GRAY2BGR);
}

void MarkerDetector_impl::drawCube(Mat& image, const Mat& R,
    const Mat& t, const float l, const Scalar &color, const Point2i &basept,
    float ratio) {

  std::vector<Point2f> pointsInImage;
  std::vector<Point3f> points;

  float s = l / 2.0;

  points.push_back(Point3f(-s, -s, 0.0));
  points.push_back(Point3f(-s, s, 0.0));
  points.push_back(Point3f(s, s, 0.0));
  points.push_back(Point3f(s, -s, 0.0));

  points.push_back(Point3f(-s, -s, -2.0 * s));
  points.push_back(Point3f(-s, s, -2.0 * s));
  points.push_back(Point3f(s, s, -2.0 * s));
  points.push_back(Point3f(s, -s, -2.0 * s));

  projectPoints(points, R, t, _cfg.K, _cfg.distortion, pointsInImage);

  for (int i = 0; i < pointsInImage.size(); ++i) {
    pointsInImage[i] = transformPoint(pointsInImage[i], basept, ratio);
  }

  // TODO: the face on the marker for now it is red
  line(image, pointsInImage[0], pointsInImage[1], Scalar(0, 0, 255),
      ratio * 1.25);
  line(image, pointsInImage[1], pointsInImage[2], Scalar(0, 0, 255),
      ratio * 1.25);
  line(image, pointsInImage[2], pointsInImage[3], Scalar(0, 0, 255),
      ratio * 1.25);
  line(image, pointsInImage[3], pointsInImage[0], Scalar(0, 0, 255),
      ratio * 1.25);
  line(image, pointsInImage[0], pointsInImage[4], color, ratio * 1.25);
  line(image, pointsInImage[1], pointsInImage[5], color, ratio * 1.25);
  line(image, pointsInImage[2], pointsInImage[6], color, ratio * 1.25);
  line(image, pointsInImage[3], pointsInImage[7], color, ratio * 1.25);
  line(image, pointsInImage[4], pointsInImage[5], color, ratio * 1.25);
  line(image, pointsInImage[5], pointsInImage[6], color, ratio * 1.25);
  line(image, pointsInImage[6], pointsInImage[7], color, ratio * 1.25);
  line(image, pointsInImage[7], pointsInImage[4], color, ratio * 1.25);
}

void MarkerDetector_impl::getDistanceGivenCenter(const EllipsePoly& elps,
    const cv::Point2d& c, double r, double &mu, double &std, int N) {

  Point2d p1, p2;
  Eigen::Vector3d vr1, vr2, vrc;

  double sum = 0;
  double sumsq = 0;

  double f = _cfg.K.at<double>(0, 0); // focal lenght
  double cx = _cfg.K.at<double>(0, 2); // principal point x in pixels
  double cy = _cfg.K.at<double>(1, 2); // .. and y

  for (int i = 0; i < N; ++i) {
    double theta = (double) i / N * M_PI;

    getEllipseLineIntersections(elps, c.x, c.y, theta, p1, p2);

    // viewing rays
    vr1 << (p1.x - cx) / f, (p1.y - cy) / f, 1;
    vr2 << (p2.x - cx) / f, (p2.y - cy) / f, 1;
    vrc << (c.x - cx) / f, (c.y - cy) / f, 1;

    // angles

    double th1, th2;

    th1 = acos((vr1 / vr1.norm()).dot(vrc / vrc.norm()));
    th2 = acos((vr2 / vr2.norm()).dot(vrc / vrc.norm()));

    double curd = (sqrt(2) * r * sin(th1 + th2))
        / sqrt(3 - 2 * cos(2 * th1) - 2 * cos(2 * th2) + cos(2 * (th1 + th2)));

    sum += curd;
    sumsq += pow(curd, 2);
  }

  mu = sum / N;
  std = sqrt(sumsq / N - pow(mu, 2));

  if (std::isnan(std)) {
    std = 0.0;
  }
}

void MarkerDetector_impl::getPoseGivenCenter(const EllipsePoly& elps,
    const cv::Point2d& c, double r, double &d, double &phi, double &kappa,
    DebugPlotConfig *dbg) {
  Point2d p1, p2;
  Eigen::Vector3d vr1, vr2, vrc;

  double f = _cfg.K.at<double>(0, 0); // focal lenght
  double cx = _cfg.K.at<double>(0, 2); // principal point x in pixels
  double cy = _cfg.K.at<double>(1, 2); // .. and y

  double theta;

  double th1, th2;

  Eigen::Array4d G1, B2, res; // candidate solutions
  int i;
  double ad, gamma1, beta2; // final solutions

  Eigen::Vector3d CB_phi, CB_kappa;
  double OB;

  // AROUND X AXIS
  theta = M_PI / 2;

  getEllipseLineIntersections(elps, c.x, c.y, theta, p1, p2);

  if (p2.y < p1.y) {
    swap(p1, p2);
  }

  // viewing rays
  vr1 << (p1.x - cx) / f, (p1.y - cy) / f, 1;
  vr2 << (p2.x - cx) / f, (p2.y - cy) / f, 1;
  vrc << (c.x - cx) / f, (c.y - cy) / f, 1;

  th1 = acos((vr1 / vr1.norm()).dot(vrc / vrc.norm()));
  th2 = acos((vr2 / vr2.norm()).dot(vrc / vrc.norm()));

  d = (sqrt(2) * r * sin(th1 + th2))
      / sqrt(3 - 2 * cos(2 * th1) - 2 * cos(2 * th2) + cos(2 * (th1 + th2)));

  G1(0) = G1(1) = asin(d / r * sin(th1));
  G1(2) = G1(3) = M_PI - G1(0);

  B2(0) = B2(2) = asin(d / r * sin(th2));
  B2(1) = B2(3) = M_PI - B2(0);

  res = G1 + B2 + th1 + th2 - M_PI;
  ad = res.abs().minCoeff(&i);

//  cerr << "anglediff: " << res.transpose() << endl;

  if (ad > 1e-6) {
    cerr << "WARNING, geometrical inconsistence: " << ad << endl;
  }

  gamma1 = G1(i);
  beta2 = B2(i);

  OB = r * sin(M_PI - beta2 - th2) / sin(th2);
  CB_phi = OB * vr2 - d * vrc;
  phi = atan2(CB_phi(1), CB_phi(2));

//  cerr << CB_phi.transpose() << " ("
//      << sqrt(pow(CB_phi(1), 2) + pow(CB_phi(2), 2)) << ") " << d << endl;

  // AROUND Y AXIS
  theta = 0.0;

  getEllipseLineIntersections(elps, c.x, c.y, theta, p1, p2);

  if (p2.x < p1.x) {
    swap(p1, p2);
  }

  // viewing rays
  vr1 << (p1.x - cx) / f, (p1.y - cy) / f, 1;
  vr2 << (p2.x - cx) / f, (p2.y - cy) / f, 1;
  vrc << (c.x - cx) / f, (c.y - cy) / f, 1;

  th1 = acos((vr1 / vr1.norm()).dot(vrc / vrc.norm()));
  th2 = acos((vr2 / vr2.norm()).dot(vrc / vrc.norm()));

  d = (sqrt(2) * r * sin(th1 + th2))
      / sqrt(3 - 2 * cos(2 * th1) - 2 * cos(2 * th2) + cos(2 * (th1 + th2)));

  G1(0) = G1(1) = asin(d / r * sin(th1));
  G1(2) = G1(3) = M_PI - G1(0);

  B2(0) = B2(2) = asin(d / r * sin(th2));
  B2(1) = B2(3) = M_PI - B2(0);

  res = G1 + B2 + th1 + th2 - M_PI;
  ad = res.abs().minCoeff(&i);

//  cerr << "anglediff: " << res.transpose() << endl;

  if (ad > 1e-6) {
    cerr << "WARNING, geometrical inconsistence: " << ad << endl;
  }

  gamma1 = G1(i);
  beta2 = B2(i);

  OB = r * sin(M_PI - beta2 - th2) / sin(th2);
  CB_kappa = OB * vr2 - d * vrc;
  kappa = atan2(CB_kappa(0), CB_kappa(2));

//  cerr << CB_kappa.transpose() << " ("
//      << sqrt(pow(CB_kappa(0), 2) + pow(CB_kappa(2), 2)) << ") " << d << endl;

}

double MarkerDetector_impl::evalDistanceF(const EllipsePoly &outer,
    const EllipsePoly &inner, const cv::Point2d &x, const cv::Point2d &x0) {

  double mu, std;

  double ret;

  getDistanceGivenCenter(outer, x, _cfg.markerDiameter / 2.0, mu, std, 4);
  ret = std;

//  getDistanceGivenCenter(inner, x, _cfg.markerInnerDiameter / 2.0, mu, std);
//  ret += std;

  // regularization

//  float d0 = sqrt(pow(x.x - x0.x, 2) + pow(x.y - x0.y, 2));
//  ret += 10 * exp(-pow(d0 - 10, 2) / pow(0.75, 2));

  ret += 0.1 * (pow(x.x - x0.x, 2) + pow(x.y - x0.y, 2));

  return ret;
}

void MarkerDetector_impl::getDistanceWithGradientDescent(
    const EllipsePoly& outer, const EllipsePoly& inner, const cv::Point2d x0,
    double step, double lambda, cv::Point2d& x, double tolX, double tolFun) {

  Eigen::VectorXd D;
  Point2d newx;

  x = x0;
  int it = 1;
  double f, newf;

  while (true) {

    // estimate gradient
    double fplus, fminus;
    Point2d g;

    // .. along x
    fplus = evalDistanceF(outer, inner, x + Point2d(step, 0.0), x0);
    fminus = evalDistanceF(outer, inner, x + Point2d(-step, 0.0), x0);
    g.x = (fplus - fminus) / 2.0 / step;

    // and along y
    fplus = evalDistanceF(outer, inner, x + Point2d(0.0, step), x0);
    fminus = evalDistanceF(outer, inner, x + Point2d(0.0, -step), x0);
    g.y = (fplus - fminus) / 2.0 / step;

    bool hadToReduce = false, stepDone = false;

    f = evalDistanceF(outer, inner, x, x0);

    while (lambda * norm(g) > tolX && !stepDone) {

      newx = x - lambda * g;
      newf = evalDistanceF(outer, inner, newx, x0);

      if (newf < f) {
        x = newx;
        f = newf;

        if (!hadToReduce) {
          lambda = lambda * 2.0;
        }

        stepDone = true;

//        cerr << it << " " << x << " " << newf << " " << lambda << endl;
      } else {
        lambda = lambda / 2.0;
        hadToReduce = true;
      }
    }

    if (!stepDone) {
      break;
    }

    if (fabs(f - newf) < tolFun) {
      break;
    }

    ++it;
  }
}

void MarkerDetector_impl::buildCircleHolesImage(const EllipsePoly& elps,
    const cv::Point2f& c, double targetRaduis, const string &fname, int N,
    float dx) {

  cv::Mat img(2 * N + 1, 2 * N + 1, CV_32F), normalized;

  double mu, std;

  for (int x = -N; x <= N; ++x) {
    for (int y = -N; y <= N; ++y) {
      getDistanceGivenCenter(elps, c + Point2f((float) x * dx, (float) y * dx),
          targetRaduis, mu, std, 4);

      img.at<float>(x + N, y + N) = std;
    }
  }

  ofstream f(fname);
  f << cv::format(img, "csv") << std::endl;
  f.close();

}

void MarkerDetector_impl::initZoomedSubregionSurface(Point2i center, float r,
    Mat &in, Mat &out, int size, float &ratio, Point2i &basept) {

  int roundedr = roundf(r);

  // first select the interesting region
  int tx = center.x - roundedr;
  int ty = center.y - roundedr;

  int basex = std::max(0, tx);
  int basey = std::max(0, ty);
  int width = std::min(center.x, roundedr)
      + std::min(in.cols - center.x, roundedr);
  int height = std::min(center.y, roundedr)
      + std::min(in.rows - center.y, roundedr);

  // blit it on a 16:9 black image

  Rect r1 = Rect(basex, basey, width, height);
  Rect r2 = Rect(-std::min(0, tx), -std::min(0, ty), width, height);

  Mat roi(2 * roundedr, 2 * roundedr, CV_8UC1, Scalar(0, 0, 0));
  in(r1).copyTo(roi(r2));

  // convert it to color
  Mat roibgr;
  cvtColor(roi, roibgr, CV_GRAY2BGR);

  // resize the image
  resize(roibgr, out, Size(size, size), 0, 0, INTER_NEAREST);

  basept = Point2i(basex - r2.x, basey - r2.y);
  ratio = (float) size / 2.0 / roundedr;
}

void MarkerDetector_impl::subpixelEdgeWithLeastSquares(const cv::Mat &image,
    const Ellipse &elps, const EllipsePoly &poly, float theta, float a, float b,
    cv::Point2f &subpixedge, int N, DebugPlotConfig *dbg) {

  // evaluate preliminary edge position
  Point2f e = evalEllipse(theta, elps.center, elps.size.width / 2.0,
      elps.size.height / 2.0, elps.angle * M_PI / 180.0);

  // evaluate gradient at edge position
  float g = atan2(2.0 * poly(2) * e.y + poly(1) * e.x + poly(4),
      2.0 * poly(0) * e.x + poly(1) * e.y + poly(3));

  // compute a shift with respect to the preliminary edge position based on max gradient pixel
  float old = getSubpix(image,
      Point2f(e.x + (float) (-N - 1) * cos(g),
          e.y + (float) (-N - 1) * sin(g))), cur;
  float maxDelta = 0, delta;
  int maxDeltaI;
  for (int i = -N; i <= N + 1; ++i) {
    cur = getSubpix(image,
        Point2f(e.x + (float) i * cos(g), e.y + (float) i * sin(g)));
    delta = fabs(cur - old);
    if (delta > maxDelta) {
      maxDelta = delta;
      maxDeltaI = i;
    }

    old = cur;
  }

  // get pixel values along the orthogonal line
  vector<Point2f> orth(2 * N + 1);
  vector<float> subpixvalues(2 * N + 1);

  for (int i = -N; i <= N; ++i) {
    orth[i + N].x = e.x + (float) (i + maxDeltaI) * cos(g);
    orth[i + N].y = e.y + (float) (i + maxDeltaI) * sin(g);

    subpixvalues[i + N] = getSubpix(image, orth[i + N]);

//    cerr << subpixvalues[i + N] << endl;
  }

  // swap them in case order is inverted
  float swapmult = 1.0;
  if (subpixvalues[0] > subpixvalues[2 * N]) {
    swapmult = -1.0;
    for (int i = 0; i <= N; ++i) {
      swap(subpixvalues[i], subpixvalues[2 * N - i]);
    }
  }

  // toggle local limits estimation if not provided
  if (a < 0) {
    a = getSubpix(image,
        Point2f(e.x + (float) (-(N + 1) + maxDeltaI) * cos(g),
            e.y + (float) (-(N + 1) + maxDeltaI) * sin(g)));

    b = getSubpix(image,
        Point2f(e.x + (float) (N + 1 + maxDeltaI) * cos(g),
            e.y + (float) (N + 1 + maxDeltaI) * sin(g)));

    if (swapmult == -1.0) {
      swap(a, b);
    }

    b = b - a;
  }
  //

  // do the least square thing
  // TODO: sometimes the estimation fails (dx = nan)
  float mu = 0.25, logsigma = 0;

  auto Sign = [](float x) {return x>=0 ? 1.0: -1.0;};

  Eigen::MatrixXd H(2 * N + 1, 2);
  Eigen::VectorXd err(2 * N + 1);

  Eigen::VectorXd dx(2);

  int it = 0;

  do {
    it++;

    float sigma = exp(logsigma);

    for (int i = -N; i <= N; ++i) {
      float y = i;

      // error
      err(i + N) = (2.0 * a + b
          - b
              * sqrt(
                  1.0 - exp((-2.0 * pow(mu - y, 2)) / (M_PI * pow(sigma, 2))))
              * Sign(mu - y)) / 2.0 - subpixvalues[i + N];

      // TODO: they could be optimized
      // derivative with respect to mu
      H(i + N, 0) = (b * (-mu + y) * Sign(mu - y))
          / (exp((2.0 * pow(mu - y, 2)) / (M_PI * pow(sigma, 2)))
              * sqrt(
                  1.0 - exp((-2.0 * pow(mu - y, 2)) / (M_PI * pow(sigma, 2))))
              * M_PI * pow(sigma, 2));
      // derivative with respect to sigma
      H(i + N, 1) = sigma * (b * pow(mu - y, 2) * Sign(mu - y))
          / (exp((2.0 * pow(mu - y, 2)) / (M_PI * pow(sigma, 2)))
              * sqrt(
                  1.0 - exp((-2.0 * pow(mu - y, 2)) / (M_PI * pow(sigma, 2))))
              * M_PI * pow(sigma, 3));
    }

    dx = (H.transpose() * H).inverse() * H.transpose() * err;

    mu = mu - dx(0);
    logsigma = logsigma - dx(1);

    /* write some debug output
     cerr << "it: " << it << " mu:   " << mu << " logsigma: " << logsigma
     << " ||err||: "
     << err.norm() << endl;

     cerr << "predicted: ";
     for (int i = -N; i <= N; i++) {
     float y = i;
     cerr
     << (2.0 * a + b
     - b
     * sqrt(
     1.0
     - exp(
     (-2.0 * pow(mu - y, 2)) / (M_PI * pow(sigma, 2))))
     * Sign(mu - y)) / 2.0 << " ";
     }

     cerr << "\nmeasured: ";
     printVector(subpixvalues);

     cerr << "residuals: " << err.transpose() << endl;

     cerr << "dmu: " << dx(0) << " dsigma: " << dx(1) << endl;
     //*/

  } while (fabs(dx(0)) > 1e-2 && it <= 5);

  subpixedge = e
      + Point2f((swapmult * mu + maxDeltaI) * cos(g),
          (swapmult * mu + maxDeltaI) * sin(g));

  if (dbg != NULL) {

    float ratio;
    Point2i basept;

    initZoomedSubregionSurface(elps.center, 0.75 * elps.size.width,
        dbg->rawImage, dbg->dbgImage, 1280, ratio, basept);

    for (int i = -N; i <= N; ++i) {
      circle(dbg->dbgImage, transformPoint(orth[i + N], basept, ratio),
          0.5 * ratio, Scalar(0, 0, 255), -1);
    }

    circle(dbg->dbgImage, transformPoint(subpixedge, basept, ratio),
        0.5 * ratio, Scalar(0, 255, 0), -1);

    imshow(dbg->windowName, dbg->dbgImage);
    waitKey(1);
  }
}

void MarkerDetector_impl::refineEllipseCntWithSubpixelEdges(
    const cv::Mat &image, const Target &tg, const Ellipse &elps,
    bool ignoreSignalAreas, int N, std::vector<cv::Point2f> &cnt,
    std::vector<double> &angles, DebugPlotConfig *dbg) {

  EllipsePoly elpsPoly;
  getEllipsePolynomialCoeff(elps, elpsPoly);

  // decide the angular increment, avoid to reuse information
  float elpsSize = (elps.size.width + elps.size.height) / 4.0;
  float inc = 1.0 / (elpsSize - N);

  cnt.clear();
  cnt.reserve(ceil(2.0 * M_PI / inc));

  angles.clear();
  angles.reserve(ceil(2.0 * M_PI / inc));

  float theta = 0.0;
  float val = _cfg.markerSignalStartsWith;
  int segment = 0;

//  float maxPSFwidth = 1.0;
//  float safetyAngle = maxPSFwidth / elpsSize
//      / _cfg.markerSignalRadiusPercentage; // to avoid getting the white dots
  float safetyAngle = 1 / elpsSize;

  while (theta < 2 * M_PI) {
    bool considerThisAngle = !ignoreSignalAreas;

    if (segment < _cfg.markerSignalModel.size()
        && theta > 2 * M_PI * _cfg.markerSignalModel[segment]) {
      ++segment;
      val = -val;
    }
    if (val == -1.0) {
      float minAngle =
          2 * M_PI
              * _cfg.markerSignalModel[(segment - 1)
                  % _cfg.markerSignalModel.size()]
              + safetyAngle;
      float maxAngle = 2 * M_PI
          * _cfg.markerSignalModel[(segment) % _cfg.markerSignalModel.size()]
          - safetyAngle
          + (segment >= _cfg.markerSignalModel.size() ? 2 * M_PI : 0);

      if (theta >= minAngle && theta <= maxAngle) {
        considerThisAngle = true;
      }
    }

    if (considerThisAngle) {
      Point2f edge;
      //  subpixelEdgeWithLeastSquares(image, elps, elpsPoly, heading + theta,
      //      black, white - black, edge, N, dbg);

      // disable exposure hint
      subpixelEdgeWithLeastSquares(image, elps, elpsPoly, tg.heading + theta,
          -1, -1, edge, N);

      if (!isnan(edge.x) && !isnan(edge.y)) {
        cnt.push_back(edge);
        angles.push_back(tg.heading + theta);
      }
    }

    theta += inc;
  }

  if (dbg != NULL) {

    float ratio;
    Point2i basept;

    initZoomedSubregionSurface(elps.center, 1.5 * elpsSize, dbg->rawImage,
        dbg->dbgImage, 1280, ratio, basept);

    for (int i = 0; i < cnt.size(); ++i) {
      circle(dbg->dbgImage, transformPoint(cnt[i], basept, ratio),
          0.5 * ratio, Scalar(0, 0, 255), -1);
    }

    imshow(dbg->windowName, dbg->dbgImage);
    waitKey(1);
  }
}

bool MarkerDetector_impl::measureRough(const cv::Mat &image,
    std::shared_ptr<Target> tg, DebugPlotConfig *dbg) {

  if (!tg->detected) {
    return false;
  }

  if (dbg != NULL) {
    tic();
  }

  bool pointsFound = true;
  bool success = true;

  // check and eventually allocate mask
  if (_floodfillMask.rows != image.rows + 2
      || _floodfillMask.cols != image.cols + 2) {
    _floodfillMask = Mat::zeros(image.rows + 2, image.cols + 2, CV_8UC1);
  }

  // color the mask so floodfill cannot surpass the outer circle
  for (auto it = tg->outer.cnt.begin(); it != tg->outer.cnt.end(); ++it) {
    _floodfillMask.at<unsigned char>(it->y + 1, it->x + 1) = 255;
  }

  // generate or validate seedpoints
  const unsigned int NPTS = _worldPoints.size();

  if (tg->seedPoints.empty()) {
    // get the seed points for the floodfill and the true world points
    Ellipse outerElps;
    fitEllipse(tg->outer.cnt, outerElps);

    for (int cnt = 0; cnt < _cfg.markerSignalModel.size() / 2; ++cnt) {
      int i = (_cfg.markerSignalStartsWith == 1.0 ? 0 : 1) + 2 * cnt;

      float maxAngle, minAngle, angle;

      if (i == 0) {
        minAngle = 2 * M_PI
            * (_cfg.markerSignalModel[_cfg.markerSignalModel.size() - 1] - 1);
      } else {
        minAngle = 2 * M_PI * _cfg.markerSignalModel[i - 1];
      }

      maxAngle = 2 * M_PI * _cfg.markerSignalModel[i];
      angle = 0.5 * (maxAngle + minAngle);

      tg->seedPoints.push_back(
          evalEllipse(angle + tg->heading, outerElps.center,
              outerElps.size.width / 2.0 * _cfg.markerSignalRadiusPercentage,
              outerElps.size.height / 2.0 * _cfg.markerSignalRadiusPercentage,
              outerElps.angle * M_PI / 180.0));
    }
  } else {
    if (tg->seedPoints.size() != NPTS) {
      cerr << "ERROR: not enough or too much seedpoints provided" << endl;
      assert(false);
    }
  }

  // floodfill and mask exploration to compute centroids
  Rect bounds[NPTS];
  int times = floor(253 / NPTS);

  for (int i = 0; i < NPTS; ++i) {
    floodFill(image, _floodfillMask, tg->seedPoints[i], 255, &(bounds[i]),
        (tg->white - tg->black) * 0.4,
        255,
        4 | ((2 + i * times) << 8) | CV_FLOODFILL_FIXED_RANGE
            | CV_FLOODFILL_MASK_ONLY);
  }

  unsigned int cnt[NPTS];
  for (unsigned int i = 0; i < NPTS; ++i) {
    cnt[i] = 0;
    tg->codePoints.push_back(Point2f(0.0, 0.0));
  }

  // compute centroids employing the bounds regions

  for (unsigned int i = 0; i < NPTS; ++i) {
    if (bounds[i].width == 0 || bounds[i].height == 0) {
      pointsFound = false; // I need to keep on to finish uncolouring the mask
    } else {
      for (unsigned int x = bounds[i].x; x <= bounds[i].x + bounds[i].width - 1;
          ++x) {
        for (unsigned int y = bounds[i].y;
            y <= bounds[i].y + bounds[i].height - 1; ++y) {
          unsigned char maskval = _floodfillMask.at<unsigned char>(y + 1,
              x + 1);

          if ((maskval - 2) / times == i) {
            // clear this mask point
            _floodfillMask.at<unsigned char>(y + 1, x + 1) = 0;

            // weighted average
            cnt[i] += image.at<unsigned char>(y, x);

            tg->codePoints[i].x += image.at<unsigned char>(y, x) * x;
            tg->codePoints[i].y += image.at<unsigned char>(y, x) * y;
            //*/

            /* simple average
             cnt[i]++;

             tg->codePoints[i].x += x;
             tg->codePoints[i].y += y;
             //*/

            // const_cast<Mat &>(image).at<unsigned char>(y, x) = 0;
          }
        }
      }
    }
  }

  // un-color the mask so floodfill cannot surpass the outer circle
  for (auto it = tg->outer.cnt.begin(); it != tg->outer.cnt.end(); ++it) {
    _floodfillMask.at<unsigned char>(it->y + 1, it->x + 1) = 0;
  }
  //*/

  vector<Point2f> prj_points;

  if (pointsFound) {
    for (unsigned int i = 0; i < NPTS; ++i) {
      tg->codePoints[i].x /= cnt[i];
      tg->codePoints[i].y /= cnt[i];
    }

    // solve the PnP problem
    Mat rod;
    solvePnP(_worldPoints, tg->codePoints, _cfg.K, _cfg.distortion, rod,
        tg->rought);
    Rodrigues(rod, tg->roughR);

    // reproject points and compute error

    projectPoints(_worldPoints, rod, tg->rought, _cfg.K, _cfg.distortion,
        prj_points);

    tg->meanReprojectionError = 0;

    for (unsigned int i = 0; i < NPTS; i++) {
      float err = sqrt(
          pow(prj_points[i].x - tg->codePoints[i].x, 2)
              + pow(prj_points[i].y - tg->codePoints[i].y, 2));

      if (err > 5) {
        cerr << "WARNING, high reprojection error " << i << "-th code point: "
            << err << endl;

        success = false;
      }
      tg->meanReprojectionError += err;
    }

    tg->meanReprojectionError /= NPTS;
  } else {
    success = false;
  }

  if (success) {
    tg->roughlyMeasured = true;
  }

  if (dbg != NULL) {

    dbg->tstats.tMeasureRough = toc();

    if (dbg->enabled && dbg->enableRoughMeasure) {

      float ratio = 1.0;
      Point2i basept(0, 0);
      if (dbg->blitSubRegion) {
        initZoomedSubregionSurface(tg->outer.center,
            tg->outer.r * dbg->blitRegionWidthMultiplier, dbg->rawImage,
            dbg->dbgImage, 1280, ratio, basept);
      } else {
        initColorPlotSurfacte(dbg->dbgImage, dbg);
      }

      if (success & pointsFound) {

        for (int i = 0; i < NPTS; ++i) {
          circle(dbg->dbgImage,
              transformPoint(tg->codePoints[i], basept, ratio),
              ratio * 2, Scalar(0, i * 255.0 / NPTS, 255), -1);
        }

        drawCube(dbg->dbgImage, tg->roughR, tg->rought,
            _cfg.markerDiameter * 1.5,
            Scalar(255, 0, 0),
            basept, ratio);
      }

      for (int i = 0; i < NPTS; ++i) {
        circle(dbg->dbgImage,
            transformPoint(tg->seedPoints[i], basept, ratio),
            ratio * 0.5, Scalar(i * 255.0 / NPTS, 0, 0), -1);
      }

      imshow(dbg->windowName, dbg->dbgImage);
      waitKey(1);
    }
  }

  /* check that the mask is completely black
   for (int x = 1; x < floodfillMask.cols - 1; ++x) {
   for (int y = 1; y < floodfillMask.rows - 1; ++y) {
   if (floodfillMask.at<unsigned char>(y, x) != 0) {
   namedWindow("mask");
   imshow("mask", floodfillMask);
   waitKey(-1);
   assert(false);
   }
   }
   }
   */

  return success;
}

Point2f MarkerDetector_impl::distort(const Point2f& p) {
  // To relative coordinates <- this is the step you are missing.

  double cx = _cfg.K.at<double>(0, 2);
  double cy = _cfg.K.at<double>(1, 2);
  double fx = _cfg.K.at<double>(0, 0);
  double fy = _cfg.K.at<double>(1, 1);
  double k1 = _cfg.distortion.at<double>(0);
  double k2 = _cfg.distortion.at<double>(1);
  double p1 = _cfg.distortion.at<double>(2);
  double p2 = _cfg.distortion.at<double>(3);
  double k3 = _cfg.distortion.at<double>(4);

  double x = (p.x - cx) / fx;
  double y = (p.y - cy) / fy;

  double r2 = x * x + y * y;

  // Radial distorsion
  double xDistort = x * (1 + k1 * r2 + k2 * r2 * r2 + k3 * r2 * r2 * r2);
  double yDistort = y * (1 + k1 * r2 + k2 * r2 * r2 + k3 * r2 * r2 * r2);

  // Tangential distorsion
  xDistort = xDistort + (2 * p1 * x * y + p2 * (r2 + 2 * x * x));
  yDistort = yDistort + (p1 * (r2 + 2 * y * y) + 2 * p2 * x * y);

  // Back to absolute coordinates.
  xDistort = xDistort * fx + cx;
  yDistort = yDistort * fy + cy;

  return Point2f(xDistort, yDistort);
}

void MarkerDetector_impl::drawPyramid(cv::Mat& image, const cv::Mat& R,
    const cv::Mat& t, const float l, const cv::Scalar& color) {

  std::vector<Point2f> pointsInImage;
  std::vector<Point3f> points;

  float s = l / 2.0;

  points.push_back(Point3f(-s, -s, 0.0));
  points.push_back(Point3f(-s, s, 0.0));
  points.push_back(Point3f(s, s, 0.0));
  points.push_back(Point3f(s, -s, 0.0));

  points.push_back(Point3f(0.0, 0.0, -4.0 * s));

  projectPoints(points, R, t, _cfg.K, _cfg.distortion, pointsInImage);

  line(image, pointsInImage[0], pointsInImage[1], color, 2);
  line(image, pointsInImage[1], pointsInImage[2], color, 2);
  line(image, pointsInImage[2], pointsInImage[3], color, 2);
  line(image, pointsInImage[3], pointsInImage[0], color, 2);
  line(image, pointsInImage[0], pointsInImage[4], color, 2);

  line(image, pointsInImage[4], pointsInImage[0], color, 2);
  line(image, pointsInImage[4], pointsInImage[1], color, 2);
  line(image, pointsInImage[4], pointsInImage[2], color, 2);
  line(image, pointsInImage[4], pointsInImage[3], color, 2);

}

void MarkerDetector_impl::getEllipseMatrix(const Ellipse &elps,
    Eigen::Matrix3d &Q, DebugPlotConfig *dbg) {

  double a, b, cx, cy, theta;

  a = elps.size.width / 2;
  b = elps.size.height / 2;
  theta = elps.angle * M_PI / 180.0;

  cx = elps.center.x;
  cy = elps.center.y;

  if (a < b) {
    theta += M_PI / 2.0;
    std::swap(a, b);
  }

// Compute elements of conic matrix

  double A, B, C, D, E, F;

  A = pow(a, 2) * pow(sin(theta), 2) + pow(b, 2) * pow(cos(theta), 2);
  B = 2 * (pow(b, 2) - pow(a, 2)) * sin(theta) * cos(theta);
  C = pow(a, 2) * pow(cos(theta), 2) + pow(b, 2) * pow(sin(theta), 2);
  D = -2 * A * cx - B * cy;
  E = -B * cx - 2 * C * cy;
  F = A * pow(cx, 2) + B * cx * cy + C * pow(cy, 2) - pow(a, 2) * pow(b, 2);

// Store into Q
  Q << A, B / 2.0, D / 2.0, B / 2.0, C, E / 2.0, D / 2.0, E / 2.0, F;

}

void MarkerDetector_impl::getEllipsePolynomialCoeff(const Ellipse &elps,
    EllipsePoly &poly) {
  double a, b, cx, cy, theta;

  a = elps.size.width / 2;
  b = elps.size.height / 2;
  theta = elps.angle * M_PI / 180.0;

  cx = elps.center.x;
  cy = elps.center.y;

  if (a < b) {
    theta += M_PI / 2.0;
    std::swap(a, b);
  }

// Compute elements of conic matrix

  double A, B, C, D, E, F;

  A = pow(a, 2) * pow(sin(theta), 2) + pow(b, 2) * pow(cos(theta), 2);
  B = 2 * (pow(b, 2) - pow(a, 2)) * sin(theta) * cos(theta);
  C = pow(a, 2) * pow(cos(theta), 2) + pow(b, 2) * pow(sin(theta), 2);
  D = -2 * A * cx - B * cy;
  E = -B * cx - 2 * C * cy;
  F = A * pow(cx, 2) + B * cx * cy + C * pow(cy, 2) - pow(a, 2) * pow(b, 2);

  double k =
      1.0
          / sqrt(
              (pow(A, 2) + pow(B, 2) + pow(C, 2) + pow(D, 2) + pow(E, 2)
                  + pow(F, 2)));

  poly << A, B, C, D, E, F;
  poly *= k;

}

void MarkerDetector_impl::getEllipseLineIntersections(const EllipsePoly& em,
    double x0, double y0, double theta, cv::Point2d& p1, cv::Point2d& p2) {
  Eigen::Vector3d lm;

  if (fabs(theta - M_PI / 2.0) < 1e-6) {

    p1.x = x0;
    p1.y = (-(x0 * em(1)) - em(4)
        - sqrt(
            pow(x0 * em(1) + em(4), 2)
                - 4 * em(2) * (pow(x0, 2) * em(0) + x0 * em(3) + em(5))))
        / (2 * em(2));
    p2.x = x0;
    p2.y = (-(x0 * em(1)) - em(4)
        + sqrt(
            pow(x0 * em(1) + em(4), 2)
                - 4 * em(2) * (pow(x0, 2) * em(0) + x0 * em(3) + em(5))))
        / (2 * em(2));

  } else {
    lm << -tan(theta), 1.0, tan(theta) * x0 - y0;

#   include "generated/EllipseLineIntersection.cpp"
  }
}

void MarkerDetector_impl::getProjectedCircleCenter(const Circles &in,
    cv::Point2f &projCenter, DebugPlotConfig *dbg) {

// if only 2 concentric circles
  if (in.size() == 2) {

    // fit ellipses
    const Circle &c1 = in.at(0);
    const Circle &c2 = in.at(1);

    Ellipse e1;
    Ellipse e2;
    fitEllipse(c1.cnt, e1);
    fitEllipse(c2.cnt, e2);

    // compute ellipse matrices
    Eigen::Matrix3d Q1;
    Eigen::Matrix3d Q2;
    getEllipseMatrix(e1, Q1);
    getEllipseMatrix(e2, Q2);

    // compute lambda3 matrix
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> sas(Q2 * Q1.inverse());
    Eigen::Vector3d lam;
    lam(0) = sas.eigenvalues()(0);
    lam(1) = sas.eigenvalues()(1);
    lam(2) = sas.eigenvalues()(2);

    float lam3 = 0;
    for (int i = 0; i < 2; i++) {
      for (int j = i + 1; j < 3; j++) {

        // compute metric
        float metric = 2 * abs(abs(lam(i)) - abs(lam(j)))
            / (abs(lam(i)) + abs(lam(j)));

        if (metric < 0.1) {
          lam3 = (abs(lam(i)) + abs(lam(j))) / 2;
        }
      }
    }

    // compute M = (Q1^-1) - lambda3*(Q2^-1) and normalize wrt
    // 3rd column to get homogeneous coordinates
    Eigen::Matrix3d M = Q1.inverse() - lam3 * Q2.inverse();

    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        M(i, j) = M(i, j) / M(i, 2);
      }
    }

    // Take average of (non-zero) row vectors to compute homogeneous center
    projCenter.x = (M(0, 0) + M(1, 0) + M(2, 0)) / 3;
    projCenter.y = (M(0, 1) + M(1, 1) + M(2, 1)) / 3;

  }

  if (dbg != NULL) {
    initColorPlotSurfacte(dbg->dbgImage, dbg);

    circle(dbg->dbgImage, projCenter, 3, Scalar(255, 0, 0), -1);

    imshow(dbg->windowName, dbg->dbgImage);
    waitKey(1);
  }
}

}
