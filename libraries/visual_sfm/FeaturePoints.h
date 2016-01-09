////////////////////////////////////////////////////////////////////////////
//	File:		FeaturePoints.h
//	Author:		Changchang Wu (ccwu1130@gmail.com)
//	Description :
//
//  Copyright (c) 2011  Changchang Wu (ccwu1130@gmail.com)
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation; either
//  Version 3 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  General Public License for more details.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef FEATURE_POINTS_H_INCLUDED
#define FEATURE_POINTS_H_INCLUDED
#pragma once

#include <fstream>
#include <vector>
#include <string>
using std::vector;
using std::ifstream;
using std::ofstream;
using std::string;

#include "points.h"

namespace visual_sfm {

///////////////////////////////////////////////////////////////////
// DESCRIPTOR TYPE
typedef unsigned char DTYPE;
// FEATURE LOCATION TYPE
typedef float LTYPE;

//#define DEBUG_DESCRIPTOR_BUG
////////////////////////////////////////////////////////////////////
////feature set with position info
//#define SIFT_LOCATION_NEW

class FeatureData {
  typedef struct sift_fileheader_v2 {
    int szFeature;
    int szVersion;
    int npoint;
    int nLocDim;
    int nDesDim;
  } sift_fileheader_v2;
  typedef struct sift_fileheader_v1 {
    int npoint;
    int nLocDim;
    int nDesDim;
  } sift_fileheader_v1;
  enum {
    //		READ_BUFFER_SIZE = 0x100000,
    SIFT_NAME = ('S' + ('I' << 8) + ('F' << 16) + ('T' << 24)),
    MSER_NAME = ('M' + ('S' << 8) + ('E' << 16) + ('R' << 24)),
    RECT_NAME = ('R' + ('E' << 8) + ('C' << 16) + ('T' << 24)),
    // SIFT_VERSION_2=('V'+('2'<<8)+('.'<<16)+('0'<<24)),
    // SIFT_VERSION_3=('V'+('3'<<8)+('.'<<16)+('0'<<24)),
    SIFT_VERSION_4 = ('V' + ('4' << 8) + ('.' << 16) + ('0' << 24)),
    SIFT_VERSION_5 = ('V' + ('5' << 8) + ('.' << 16) + ('0' << 24)),
    SIFT_EOF = (0xff + ('E' << 8) + ('O' << 16) + ('F' << 24)),
  };
  //	static char readBuf[READ_BUFFER_SIZE];
  //	static char sift_version[8];
  static inline int IsValidFeatureName(int value) {
    return value == SIFT_NAME || value == MSER_NAME;
  }
  static inline int IsValidVersionName(int value) {
    return value == SIFT_VERSION_4 || value == SIFT_VERSION_5;
  }

 public:
  class LocationData : public Points<LTYPE> {
   public:
    int _file_version;

   public:
    // each feature has x,y, z, size, orientation
    // for rect feature, there is x, y, width, height
    // for eclips feature, there is u,v,a,b,c +z
   public:
    void glPaint2D(int style);
    LocationData(int d, int n) : Points<LTYPE>(d, n), _file_version(0) {};
    LocationData(LocationData& sup, int index[], int n)
        : Points<LTYPE>(sup, index, n), _file_version(0) {};
    static void glPaintTexSIFT(const float* loc, float td[3]);
    static void glPaintTexFrontalVIP(const float* loc, float td[3]);
    static void glPaintSIFT(const LTYPE* loc);
    static void glPaintSIFTSQ(LTYPE* loc);
    static void glPaintELIPS(LTYPE* loc);
    static void glPaintRECT(LTYPE* loc);
    static void SetPaintColor(LTYPE sz);
  };

  class DescriptorData : public Points<DTYPE> {
    // generally d is 128
   public:
    DescriptorData(DescriptorData& sup, int index[], int n)
        : Points<DTYPE>(sup, index, n) {};
    DescriptorData(int d, int n) : Points<DTYPE>(d, n) {};
  };
  static float gSiftDisplayScale;
  static int gSiftVisualStyle;

 protected:
  void ReadFeatureDataA(ifstream& is, int npt, int locDim, int desDim,
                        int szFeatureName);
  // the set of feature descriptors
  DescriptorData* _desData;
  // the set of feature locations
  LocationData* _locData;
  int _npoint;
  int _updated;

 public:
  void SetUpdated() { _updated = 1; }
  int GetUpdated() { return _updated; }
  void CopyToFeatureData(FeatureData& fd);
  int appendSIFTB(const char* szFile, int pos);
  int validate() { return _locData && _desData; }
  void ResizeFeatureData(int npoint, int locDim = 5, int desDim = 128) {
    if (npoint == 0) {
      if (_locData) delete _locData;
      if (_desData) delete _desData;
      _locData = NULL;
      _desData = NULL;
    } else {
      if (_locData)
        _locData->resize(locDim, npoint);
      else
        _locData = new LocationData(locDim, npoint);
      if (_desData)
        _desData->resize(desDim, npoint);
      else
        _desData = new DescriptorData(desDim, npoint);
      _locData->_file_version = SIFT_VERSION_4;
    }
    _npoint = npoint;
  }
  void operator=(FeatureData& ref) { ref.CopyToFeatureData(*this); }
  void ResizeLocationData(int npoint, int locDim) {
    if (npoint == 0) {
      if (_locData) delete _locData;
      if (_desData) delete _desData;
      _locData = NULL;
      _desData = NULL;
    } else {
      if (_locData)
        _locData->resize(locDim, npoint);
      else
        _locData = new LocationData(locDim, npoint);
      if (_desData) {
        delete _desData;
        _desData = NULL;
      }
      _locData->_file_version = SIFT_VERSION_4;
    }
  }

  void ShrinkLocationData(int ndim = 2, int npoint = -1);
  void ReleaseDescriptorData() {
    if (_desData) {
      delete _desData;
      _desData = NULL;
    }
  }

 public:
  void SortSIFT();
  void SaveSIFTBClip(const char* szFileName, int x1, int x2, int y1, int y2);
  int ValidateIndex(int index[], int n);
  void ConvertA2B(const char* szFile);
  void SaveDescriptorFile(const char* szFeatureFile);
  void SaveLocationFile(const char* szFeatureFile);
  static int ReadRandomFeature(const char* path, double* pt, int nDim);
  static int ExportSIFT(const char* featurefile, Points<double>& points);
  static int OpenSeekSIFT(const char* featurefile, int& npoint, int& ndim);
  static int ExportSIFT(const char* featurefile, int fdx);
  int ReadSIFTB(const char* szFile);
  int ReadSIFTB_LOC(const char* szFile);
  int ReadSIFTB_LOCT(const char* szFile, int fmax);
  int ReadSIFTB_DES(const char* szFile, int fmax);
  static int ReadSIFTB_DES(const char* szFile, unsigned char* buf, int nmax);
  static int ReadSIFTA_DES(const char* szFile, unsigned char* buf, int nmax);
  static int ReadSIFTB_LOC(const char* szFile, float* buf, int nmax);
  static int ReadSIFTB(const char* szFile, float* locbuf, unsigned char* desbuf,
                       int nmax);
  void offsetFeatures(int ox, int oy);
  void saveSIFTB2(const char* szFile);
  void saveSIFTA(const char* szFile, const int nDesDim = 128);
  void Overwrite_SIFTLOC(const char* szFile);
  void CreateFeatureClip(FeatureData& src, int x1, int x2, int y1, int y2,
                         int smax = 100000);
  void SetFeatureClip(FeatureData& src, int nsf, int index[], int subset = 1);
  void ReleaseFeatureData();
  void FlipLocation(float imheight);
  int ReadSIFTA(const char* szFile);
  FeatureData();
  virtual ~FeatureData();
  DescriptorData& getDescriptorData() { return *_desData; }
  LocationData& getLocationData() const { return *_locData; }
  int IsValidFeatureData() { return getFeatureNum() > 0; }
  int getFeatureNum() { return _npoint; }
  int getLoadedFeatureNum() { return _locData ? _locData->npoint() : 0; }
};

typedef FeatureData::LocationData FeatureLocationData;
typedef FeatureData::DescriptorData FeatureDescriptorData;

}  // namespace visual_sfm

#endif  // !defined(FEATURE_POINTS_H_INCLUDED)
