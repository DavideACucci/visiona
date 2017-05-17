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
 * batchdetector.cpp
 *
 *  Created on: May 13, 2015
 *      Author: Davide A. Cucci (davide.cucci@epfl.ch)
 */

#include <fstream>
#include <iostream>
#include <vector>

#include <dirent.h>
#include <unistd.h>
#include <getopt.h>

#include "opencv2/highgui/highgui.hpp"

#include <libconfig.h++>

#include "Visiona.h"
#include "Timer.h"

using namespace std;
using namespace cv;
using namespace libconfig;
using namespace visiona;

#define OPTPATHSET 1
#define OPTEXTSET 2
#define OPTONLYONEFRAME 4
#define OPTCFGFILESET 8
#define OPTDEBUG 16
#define OPTWAITKEYPRESS 32
#define OPTSTOPAT 64
#define OPTPREFIXSET 128
#define OPTSKIPDETECTION 256
#define OPTUSESEEDPOINTS 512

struct ImageDesc {
    string fileName;
    long int frameNumber;

    ImageDesc(const string &name, long int frame) :
        fileName(name), frameNumber(frame) {
    }

    bool operator <(const ImageDesc & b) const {
      return frameNumber < b.frameNumber;
    }
};

void loadCircle(const Setting &points, Circle &out);

template<typename PointType>
void load2DPointVector(const Setting &points, vector<PointType> &v) {
  int N = points.getLength();
  for (int i = 0; i < N; ++i) {
    PointType pt(points[i][0], points[i][1]);
    v.push_back(pt);
  }
}

int main(int argc, char *argv[]) {

  // --------------------- parse command line arguments ------------------------

  static struct option long_options[] = {
      { "path", required_argument, 0, 'p' },
      { "ext", required_argument, 0, 'e' },
      { "start-frame", required_argument, 0, 's' },
      { "stop-frame", required_argument, 0, 't' },
      { "only-frame", required_argument, 0, 'o' },
      { "wait", no_argument, 0, 'w' },
      { "config", required_argument, 0, 'c' },
      { "debug", no_argument, 0, 'd' },
      { "prefix", required_argument, 0, 'f' },
      { "skip-detection", required_argument, 0, 'k' },
      { 0, 0, 0, 0 }
  };

  unsigned int optionflag = 0;
  char *imagePath, *imgext, *configpath, *prefix, *detectioncfgpath,
      *seedpointspath;
  long int startFrom = 0, stopAt = -1, onlyFrame;

  opterr = 0;
  int c;
  while ((c = getopt_long_only(argc, argv, "", long_options, NULL)) != -1) {
    switch (c) {
    case 'p':

      if (optarg[strlen(optarg) - 1] == '/') {
        optarg[strlen(optarg) - 1] = 0;
      }

      imagePath = optarg;
      optionflag |= OPTPATHSET;
      break;
    case 'e':
      imgext = optarg;
      optionflag |= OPTEXTSET;
      break;
    case 's':
      startFrom = atoi(optarg);
      break;
    case 't':
      stopAt = atoi(optarg);
      optionflag |= OPTSTOPAT;
      break;
    case 'o':
      onlyFrame = atoi(optarg);
      optionflag |= OPTONLYONEFRAME;
      break;
    case 'w':
      optionflag |= OPTWAITKEYPRESS;
      break;
    case 'c':
      configpath = optarg;
      optionflag |= OPTCFGFILESET;
      break;
    case 'd':
      optionflag |= OPTDEBUG;
      break;
    case 'f':
      prefix = optarg;
      optionflag |= OPTPREFIXSET;
      break;
    case 'k':
      detectioncfgpath = optarg;
      optionflag |= OPTSKIPDETECTION;
      break;
    case '?':
      cerr << " * ERROR: unknown option or missing argument" << endl;
      return 1;

    default:
      abort();
    }
  }

  if ((optionflag & OPTPATHSET) == 0) {
    cerr << " * ERROR: image path not specified (-p)" << endl;
    return 1;
  }
  if ((optionflag & OPTEXTSET) == 0) {
    cerr << " * ERROR: image extension not specified (-e)" << endl;
    return 1;
  }

  DIR *dir;
  if ((dir = opendir(imagePath)) == NULL) {
    cerr << " * ERROR: could not open image folder" << endl;
    return 1;
  }

  // --------------------- configuration ---------------------------------------

  MarkerDetectorConfig cfg;

  if (optionflag & OPTCFGFILESET) {
    if (!cfg.loadConfig(configpath)) {
      return 1;
    }
  }

  MarkerDetector *d = MarkerDetectorFactory::makeMarkerDetector(cfg);

  shared_ptr<Target> tgfromcfg(new Target);

  if (optionflag & OPTSKIPDETECTION) {
    Config cfgfile;
    try {
      cfgfile.readFile(detectioncfgpath);
    } catch (const FileIOException &fioex) {
      cerr << " * ERROR: I/O error while reading " << detectioncfgpath << endl;
      return false;
    } catch (const ParseException &pex) {
      cerr << " * ERROR: malformed cfg file at " << pex.getFile() << ":"
          << pex.getLine() << " - " << pex.getError() << endl;
      return false;
    }

    const Setting &root = cfgfile.getRoot();

    tgfromcfg->detected = true;

    loadCircle(root["OuterPoints"], tgfromcfg->outer);
    loadCircle(root["InnerPoints"], tgfromcfg->inner);
    load2DPointVector(root["SeedPoints"], tgfromcfg->seedPoints);
  }

  // --------------------- generating image list -------------------------------

  vector<ImageDesc> images;

  struct dirent *file;

  while ((file = readdir(dir)) != NULL) {
    if (strcmp(file->d_name + strlen(file->d_name) - 3, imgext) == 0) {

      string fn_prefix = "_";
      if (optionflag & OPTPREFIXSET) {
        fn_prefix = string(prefix);
      }

      string fname(file->d_name);

      int sep = fname.find(fn_prefix);
      if (sep != string::npos) {
        stringstream s(fname.substr(sep + fn_prefix.length(), 16));

        unsigned int frame;
        s >> frame;

        if ((optionflag & OPTONLYONEFRAME) && frame != onlyFrame) {
          continue;
        }

        if (frame < startFrom) {
          continue;
        }

        if ((optionflag & OPTSTOPAT) && frame > stopAt) {
          continue;
        }

        images.push_back(ImageDesc(fname, frame));
      }

    }
  }

  closedir(dir);

  sort(images.begin(), images.end());

  // --------------------- prepare output files --------------------------------

  vector<ofstream *> imf;

  for (int i = 0; i < cfg.markerSignalModel.size() / 2; ++i) {
    stringstream s;
    s << "CP" << fixed << setfill('0') << setw(2) << i << ".txt";

    imf.push_back(new ofstream(s.str()));
  }

  // --------------------- process every image ---------------------------------

  for (auto it = images.begin(); it != images.end(); ++it) {

    cerr << it->frameNumber << " - " << it->fileName << " ..." << endl;

    string imgName = imagePath + string("/") + it->fileName;

    Mat raw = imread(imgName, CV_LOAD_IMAGE_GRAYSCALE);

    DebugPlotConfig *dbg = NULL;

    if (optionflag & OPTDEBUG) {
      dbg = new DebugPlotConfig;
      dbg->rawImage = raw;

      namedWindow("Debug", CV_GUI_EXPANDED);
      dbg->windowName = "Debug";

      dbg->blitSubRegion = true;
      dbg->blitRegionWidthMultiplier = 3.0;

      dbg->enabled = true;
      dbg->enableCirclesClusters = false;
      dbg->enableSelectedTargets = false;
      dbg->enableEsposure = false;
      dbg->enableRoughMeasure = false;
      dbg->enableSubPixelEllipses = true;

      dbg->frameNumber = it->frameNumber;
      dbg->writeContours = true;
      dbg->writeSubpixelContours = true;
      dbg->debugFilesPath = "contours";
    }

    if (optionflag & OPTDEBUG) {
      tic();
    }

    // --- real detection proces starts here

    shared_ptr<Target> tg;

    if (optionflag & OPTSKIPDETECTION) {
      tg = tgfromcfg;
    } else {
      vector<shared_ptr<Target>> ret = d->detect(raw, dbg);
      tg = ret[0];
    }

    if (tg->detected) {
      d->evaluateExposure(raw, tg, dbg);

      d->measureRough(raw, tg, dbg);

      d->measure(raw, tg, dbg);
    }

    // --- and ends here

    if (optionflag & OPTDEBUG) {
      dbg->tstats.tTotal = toc();

      stringstream fname;
      fname << "debug/dbg_" << setfill('0') << setw(6)
          << dbg->frameNumber << ".jpg";
      imwrite(fname.str(), dbg->dbgImage);

      waitKey((optionflag & OPTWAITKEYPRESS) ? 0 : 1);
    }

    // --- where output is produced

    // TODO: introduce a flag to enable/disable this

    if (tg->roughlyMeasured) {

      for (int i = 0; i < cfg.markerSignalModel.size() / 2; ++i) {
        double x = tg->codePoints[i].x, y = tg->codePoints[i].y;

        // convert to photogrammetry convention
        // TODO: put an option
        if (true) {
          swap(x, y);
          x = -(x - raw.rows / 2.0) * 4.7e-3;
          y = -(y - raw.cols / 2.0) * 4.7e-3;
        }

        (*imf[i]) << it->fileName.substr(0, it->fileName.length() - 4) << " ";
        (*imf[i]) << fixed << setprecision(6) << x << " ";
        (*imf[i]) << fixed << setprecision(6) << y << endl;

        imf[i]->flush();
      }
    }

    delete dbg;
  }

  return 0;
}

void loadCircle(const Setting &points, Circle &out) {
  int N = points.getLength();

  for (int i = 0; i < N; ++i) {
    Point2f pt(points[i][0], points[i][1]);
    out.cnt.push_back(pt);
    out.center += pt;
  }
  out.center.x /= N;
  out.center.y /= N;

  for (auto it = out.cnt.begin(); it != out.cnt.end();
      ++it) {
    out.r += sqrt(pow(it->x - out.center.x, 2) + pow(it->y - out.center.y, 2));
  }

  out.r /= N;
}

