////////////////////////////////////////////////////////////////////////////
//	File:		FeaturePoints.cpp
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

#include "stdio.h"
#include <algorithm>
#include <cassert>
#include <fcntl.h>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
using namespace std;

#include "FeaturePoints.h"

#ifndef O_BINARY
#define O_BINARY  0
#define O_TEXT    0
#endif

namespace visual_sfm {

float FeatureData::gSiftDisplayScale = 6.0f;
int FeatureData::gSiftVisualStyle = 0;

//////////////////////////////////////////////////////////////////////
// Construction/Destruction

FeatureData::FeatureData() {
  _desData = NULL;
  _locData = NULL;
  _updated = 0;
  _npoint = 0;
}

FeatureData::~FeatureData() {

  if (_desData) delete _desData;
  if (_locData) delete _locData;
}

int FeatureData::ReadSIFTA(const char* szFile) {
  int npt = 0, nd = 0;
  ifstream is(szFile);
  if (!is.is_open()) return 0;
  is >> npt >> nd;
  if (npt <= 0 || nd <= 0) {
    _npoint = 0;
    is.close();
    return 0;
  } else {
    ReadFeatureDataA(is, npt, 4, nd, SIFT_NAME);
    is.close();
    SetUpdated();
    return 1;
  }
}

void FeatureData::ReadFeatureDataA(ifstream& is, int npt, int locDim,
                                   int desDim, int szFeatureName) {
  int i, j;

  ///////////////////////////////////////
  _npoint = npt;
  _locData = new LocationData(locDim + 1, npt);  // one more column for z
  _desData = new DescriptorData(desDim, npt);
  LTYPE* pLoc = _locData->data();
  DTYPE* pDes = _desData->data();

  ////////////////////////////////////////////////
  for (i = 0; i < npt; i++) {

    if (locDim == 4) {
      is >> *(pLoc + 1);  // x
      is >> *pLoc;  //
    } else {
      is >> *pLoc;  // x
      is >> *(pLoc + 1);  // y
    }
    // SIFT_LOCATION_NEW
    pLoc[0] += 0.5f;
    pLoc[1] += 0.5f;
    pLoc[2] = 0;
    pLoc += 3;
    if (locDim == 4) {
      is >> pLoc[0];
      is >> pLoc[1];
      pLoc[1] = -pLoc[1];  // flip the angle
      pLoc += 2;
    } else {
      is >> *pLoc++;
      is >> *pLoc++;
      is >> *pLoc++;
    }

    for (j = 0; j < desDim; j++) {
      unsigned int value;
      is >> value;
      *pDes++ = value;
    }
  }
}

void FeatureData::CopyToFeatureData(FeatureData& fd) {
  if (getFeatureNum() > 0) {
    fd.ResizeFeatureData(getFeatureNum(), getLocationData().ndim(), 128);
    fd.getLocationData() = getLocationData();
    fd.getDescriptorData() = getDescriptorData();
    // memcpy(fd.getLocationData().data(), getLocationData().data(),
    // getLocationData().bsize());
    // memcpy(fd.getDescriptorData().data(), getDescriptorData().data(),
    // getDescriptorData().bsize());
  }
}

void FeatureData::FlipLocation(float imheight) {

  if (_locData == NULL) return;
  LTYPE* feature;
  if (_locData->ndim() == 5) {
    for (int i = 0; i < _locData->npoint(); i++) {
      feature = _locData->getpt(i);
      feature[1] = imheight - feature[1];
    }
  } else if (_locData->ndim() == 6) {
    for (int i = 0; i < _locData->npoint(); i++) {
      feature = _locData->getpt(i);
      //	SIFT_LOCATION_NEW
      feature[1] = imheight - feature[1];
      feature[4] = -feature[4];
    }
  }
}

void FeatureData::ReleaseFeatureData() {
  if (_locData) delete _locData;
  if (_desData) delete _desData;
  _locData = NULL;
  _desData = NULL;
  _npoint = 0;
}

void FeatureData::SetFeatureClip(FeatureData& src, int nsf, int index[],
                                 int subset) {
  ReleaseFeatureData();
  if (subset) {
    if (src._locData)
      _locData = new LocationData(src.getLocationData(), index, nsf);
    if (src._desData)
      _desData = new DescriptorData(src.getDescriptorData(), index, nsf);
  } else {
    FeatureData fd2;
    fd2.SetFeatureClip(src, nsf, index, 1);
    ResizeFeatureData(nsf, src._locData->width(), src._desData->width());
    getLocationData() = fd2.getLocationData();
    getDescriptorData() = fd2.getDescriptorData();
  }
  SetUpdated();
}

void FeatureData::offsetFeatures(int ox, int oy) {
  if (_locData == NULL) return;
  LTYPE dox = (LTYPE)ox, doy = (LTYPE)oy;
  LTYPE* feature;
  for (int i = 0; i < _locData->npoint(); i++) {
    feature = _locData->getpt(i);
    feature[0] = feature[0] + dox;
    feature[1] = feature[1] + doy;
  }
}

void FeatureData::saveSIFTA(const char* szFile, const int nDesDim) {
  if (_locData == NULL) return;
  int n = _locData->npoint();
  if (n <= 0) return;
  ofstream out(szFile);
  out << n << " " << nDesDim << "\n";
  //
  for (int i = 0; i < n; i++) {
    LTYPE* loc = _locData->getpt(i);
    unsigned char* des = _desData->getpt(i);
    out << loc[1] << " " << loc[0] << " " << loc[3] << " " << loc[4] << endl;
    for (int k = 0; k < nDesDim; k++) {
      out << ((unsigned int)des[k]) << " ";
      if ((k + 1) % 20 == 0) out << endl;
    }
    out << endl;
  }
}

void FeatureData::saveSIFTB2(const char* szFile) {

  int i, j, sift_eof = SIFT_EOF;
  sift_fileheader_v2 sfh;
  int fd = open(szFile, O_BINARY | O_CREAT | O_WRONLY | O_TRUNC,
                 S_IREAD | S_IWRITE);
  if (fd < 0) return;

  ///
  sfh.szFeature = SIFT_NAME;
  sfh.szVersion = SIFT_VERSION_4;
  sfh.npoint = _locData->npoint();
  sfh.nLocDim = _locData->ndim();
  sfh.nDesDim = _desData->ndim();
  write(fd, &sfh, sizeof(sfh));
  ////
  LTYPE* lp;
  DTYPE* dp;
  unsigned char* fph;
  float* fp;
  unsigned char* ucp;
  int Max, MemNum;
  lp = _locData->data();
  MemNum = sfh.nDesDim * sfh.npoint;
  Max = sfh.npoint * sfh.nLocDim;
  fph = new unsigned char[MemNum];  // MemCollect::Malloc(MemNum);
  fp = (float*)fph;

  for (i = 0; i < sfh.npoint; i++) {
    lp = (*_locData)[i];
    for (j = 0; j < sfh.nLocDim; j++) {
      *fp++ = (float)*lp++;
    }
  }
  write(fd, fph, sizeof(float) * Max);

  dp = _desData->data();
  Max = sfh.npoint * sfh.nDesDim;
  ucp = (unsigned char*)fph;
  for (i = 0; i < sfh.npoint; i++) {
    dp = (*_desData)[i];
    for (j = 0; j < sfh.nDesDim; j++, dp++) {
      *ucp++ = *dp;
    }
  }
  write(fd, fph, sizeof(unsigned char) * Max);
  write(fd, &sift_eof, sizeof(int));
  close(fd);
}

int FeatureData::ReadSIFTB(const char* szFile) {
  int name, version, npoint, nLocDim, nDesDim, sift_eof;
  int fd = open(szFile, O_BINARY | O_RDONLY, S_IREAD);
  if (fd < 0) return 0;
  ///
  read(fd, &name, sizeof(int));
  read(fd, &version, sizeof(int));
  if (IsValidFeatureName(name) && IsValidVersionName(version)) {
    // version 2 file
    read(fd, &npoint, sizeof(int));
    read(fd, &nLocDim, sizeof(int));
    read(fd, &nDesDim, sizeof(int));
    if (npoint > 0 && nLocDim > 0 && nDesDim == 128) {
      ResizeFeatureData(npoint, nLocDim, nDesDim);
      read(fd, _locData->data(), nLocDim * npoint * sizeof(float));
      read(fd, _desData->data(), nDesDim * npoint * sizeof(unsigned char));
      read(fd, &sift_eof, sizeof(int));
      // assert(sift_eof == SIFT_EOF);
      close(fd);
      _locData->_file_version = version;
      SetUpdated();
#ifdef DEBUG_DESCRIPTOR_BUG
      unsigned char* pb = desData->data();
      int bad_descriptor_count = 0;
      for (int i = 0; i < npoint; ++i, pb += 128) {
        int good_descriptor = 0;
        for (int j = 0; j < 128; ++j) {
          if (pb[j] != 45) {
            good_descriptor = 1;
            break;
          }
        }
        if (good_descriptor == 0) bad_descriptor_count++;
      }
      if (bad_descriptor_count)
        printf("xxxxxx %d bad descriptor xxxxxx\r\n", bad_descriptor_count);
#endif
    } else {
      ResizeFeatureData(0, 0, 0);
      close(fd);
      return 0;
    }
    return 1;
  } else {
    close(fd);
    return 0;
  }
}

int FeatureData::ReadSIFTB_DES(const char* szFile, int fmax) {
  sift_fileheader_v2 sfh;
  int fd = open(szFile, O_BINARY | O_RDONLY, S_IREAD);
  if (fd < 0) return 0;
  ///
  read(fd, &sfh, sizeof(sfh));

  int npoint = fmax > 0 ? min(sfh.npoint, fmax) : sfh.npoint;
  if (_desData)
    _desData->resize(sfh.nDesDim, npoint);
  else
    _desData = new DescriptorData(sfh.nDesDim, npoint);
  lseek(fd, sfh.nLocDim * sfh.npoint * sizeof(float), SEEK_CUR);
  read(fd, _desData->data(), sfh.nDesDim * npoint * sizeof(unsigned char));
  close(fd);
  return 1;
}

int FeatureData::ReadSIFTB_LOC(const char* szFile, float* buf, int nmax) {
  sift_fileheader_v2 sfh;
  int fd = open(szFile, O_BINARY | O_RDONLY, S_IREAD);
  if (fd < 0) return 0;
  ///
  read(fd, &sfh, sizeof(sfh));

  nmax = min(nmax, sfh.npoint);

  read(fd, buf, sfh.nLocDim * sfh.npoint * sizeof(float));
  close(fd);
  return nmax;
}

int FeatureData::ReadSIFTB_DES(const char* szFile, unsigned char* buf,
                               int nmax) {
  sift_fileheader_v2 sfh;
  int fd = open(szFile, O_BINARY | O_RDONLY, S_IREAD);
  if (fd < 0) return 0;
  ///
  read(fd, &sfh, sizeof(sfh));

  nmax = min(nmax, sfh.npoint);

  lseek(fd, sfh.nLocDim * sfh.npoint * sizeof(float), SEEK_CUR);
  read(fd, buf, sfh.nDesDim * nmax * sizeof(unsigned char));
  close(fd);
  return nmax;
}

int FeatureData::ReadSIFTA_DES(const char* szFile, unsigned char* buf,
                               int nmax) {
  FeatureData fd;
  fd.ReadSIFTA(szFile);
  if (fd.getLocationData().getpt(0)[3] <
      fd.getLocationData().getpt(fd.getFeatureNum() - 1)[3]) {
    // sort
    fd.SortSIFT();
  }
  nmax = min(nmax, fd.getFeatureNum());
  if (nmax > 0) {
    memcpy(buf, fd.getDescriptorData().data(), nmax * 128);
  }
  return nmax;
}

int FeatureData::ReadSIFTB(const char* szFile, float* locbuf,
                           unsigned char* desbuf, int nmax) {
  sift_fileheader_v2 sfh;
  int fd = open(szFile, O_BINARY | O_RDONLY, S_IREAD);
  if (fd < 0) return 0;
  ///
  read(fd, &sfh, sizeof(sfh));

  nmax = min(nmax, sfh.npoint);

  read(fd, locbuf, sfh.nLocDim * sfh.npoint * sizeof(float));
  read(fd, desbuf, sfh.nDesDim * nmax * sizeof(unsigned char));
  close(fd);
  return nmax;
}

int FeatureData::ReadSIFTB_LOC(const char* szFile) {
  sift_fileheader_v2 sfh;
  int fd = open(szFile, O_BINARY | O_RDONLY, S_IREAD);
  if (fd < 0) return 0;
  ///
  read(fd, &sfh, sizeof(sfh));
  if (IsValidFeatureName(sfh.szFeature) && IsValidVersionName(sfh.szVersion)) {
    // version 2 file
    if (sfh.npoint > 0 && sfh.nLocDim > 0 && sfh.nDesDim == 128) {
      _npoint = sfh.npoint;
      ResizeLocationData(sfh.npoint, sfh.nLocDim);
      if (_locData->data()) {
        read(fd, _locData->data(), sfh.nLocDim * sfh.npoint * sizeof(float));
        _locData->_file_version = sfh.szVersion;
      } else {
        printf("ERROR: ReadSIFTB_LOC allocatoin faled\r\n");
      }
      close(fd);
      SetUpdated();
    } else {
      ResizeFeatureData(0, 0, 0);
      close(fd);
      return 0;
    }
    return 1;
  } else {
    _npoint = 0;
    close(fd);
    return 0;
  }
}

int FeatureData::appendSIFTB(const char* szFile, int pos) {
  int npoint, nLocDim, nDesDim;
  int fd = open(szFile, O_BINARY | O_RDONLY, S_IREAD);
  if (fd < 0) return pos;
  ///
  ///
  read(fd, &npoint, sizeof(int));
  read(fd, &nLocDim, sizeof(int));

  if (npoint == SIFT_NAME && IsValidVersionName(nLocDim)) {
    read(fd, &npoint, sizeof(int));
    read(fd, &nLocDim, sizeof(int));
    // assert(nLocDim == 5);
    read(fd, &nDesDim, sizeof(int));  // assert(nDesDim == 128);
    read(fd, (float*)_locData->data() + pos * nLocDim,
          nLocDim * npoint * sizeof(float));
    read(fd, (unsigned char*)_desData->data() + pos * nDesDim,
          nDesDim * npoint * sizeof(unsigned char));
  }
  close(fd);

  return pos + npoint;
}

int FeatureData::OpenSeekSIFT(const char* featurefile, int& npoint, int& ndim) {
  sift_fileheader_v2 sfh = {SIFT_NAME, SIFT_VERSION_4, 0, 0, 0};

  int fd = open(featurefile, O_BINARY | O_RDONLY);
  if (fd < 0) return 0;
  ///
  read(fd, &sfh, sizeof(sfh));
  if (sfh.szFeature != SIFT_NAME || !IsValidVersionName(sfh.szVersion) ||
      sfh.npoint <= 0 || sfh.nLocDim <= 0 || sfh.nDesDim != 128) {
    // incorrect version
    close(fd);
    return 0;
  } else {
    lseek(fd, sfh.nLocDim * sfh.npoint * sizeof(float), SEEK_CUR);
    npoint = sfh.npoint;
    ndim = sfh.nDesDim;
    return fd;
  }
}

void FeatureData::SaveLocationFile(const char* szFeatureFile) {
  int nf = _locData->npoint();
  if (nf == 0) return;
  //	struct _stat buf;
  //	if(_stat(szFeatureFile, &buf)!=-1) return; //already saved
  vector<float> loc(nf * 2);
  float* p = &loc[0];
  int fid = open(szFeatureFile, O_CREAT | O_BINARY | O_WRONLY | O_TRUNC,
                  S_IREAD | S_IWRITE);
  if (fid == -1) return;

  for (int i = 0; i < nf; i++) {
    *p++ = (float)_locData->getpt(i)[0];
    *p++ = (float)_locData->getpt(i)[1];
  }

  write(fid, &loc[0], loc.size() * sizeof(float));
  close(fid);
}

void FeatureData::SaveDescriptorFile(const char* szFeatureFile) {
  int fid = open(szFeatureFile, O_CREAT | O_BINARY | O_WRONLY | O_TRUNC,
                  S_IREAD | S_IWRITE);

  int nf = _desData->npoint();
  DTYPE* dp = _desData->data();
  write(fid, dp, nf * 128 * sizeof(unsigned char));
  close(fid);
}

void FeatureData::ConvertA2B(const char* szFile) {

  // LET'S AUTOMATICALLY GENERATE A BINARY VRSION

  // make sure it is the same format of sift

  double normal = 0;
  DTYPE* d = _desData->data();

  for (int i = 0; i < 128; i++) normal += d[i] * d[i];
  normal = sqrt(normal);

  if (normal < 400) {
    //__asm int 03h;
    for (int i = 0; i < _desData->npoint(); i++) {
      d = _desData->getpt(i);
      normal = 0;
      for (int j = 0; j < 128; j++) {
        normal += d[j] * d[j];
      }

      normal = 512.0 / sqrt(normal);
      for (int k = 0; k < 128; k++) {
        d[k] = (DTYPE)floor(normal * d[k] + 0.5);
      }
    }

    for (int j = 0; j < _locData->npoint(); j++) {
      LTYPE* l = _locData->getpt(j);
      LTYPE temp = l[0];
      l[0] = l[1];
      l[1] = temp;
    }
  }

  int len = strlen(szFile);
  if (len > 10 && strcmp(szFile + len - 10, ".sift.sift") == 0) {
    char szNewFile[260];
    strcpy(szNewFile, szFile);
    szNewFile[len - 5] = 0;
    saveSIFTB2(szNewFile);
    remove(szFile);
  } else if (len > 5 && strcmp(szFile + len - 5, ".sift") == 0) {
    saveSIFTB2(szFile);
  }
}

int FeatureData::ValidateIndex(int index[], int n) {
  int num = getFeatureNum();
  for (int i = 0; i < n; i++) {
    if (index[i] >= num || index[i] < 0) {

      return 0;
    }
  }
  return 1;
}

void FeatureData::SortSIFT() {
  int i, j, nf = getFeatureNum();
  if (nf == 0 || _locData == NULL || _desData == NULL) return;
  LTYPE* loc = _locData->data();
  vector<int> orderv(nf);
  vector<LTYPE> sizes(nf);
  int* order = &orderv[0];
  LTYPE* size = &sizes[0];

  // bubble sort
  for (i = 0; i < nf; i++) order[i] = nf - i - 1;
  for (i = 0; i < nf; i++) size[i] = loc[5 * (nf - i - 1) + 3];
  for (i = 0; i < nf; i++) {
    int swaped = 0;
    for (j = nf - 1; j > i; j--) {
      if (size[j] > size[j - 1]) {
        // size
        LTYPE dtemp = size[j];
        size[j] = size[j - 1];
        size[j - 1] = dtemp;

        // order
        int itemp = order[j];
        order[j] = order[j - 1];
        order[j - 1] = itemp;
        swaped = 1;
      }
    }
    if (!swaped) break;
  }
  _locData->reorder(order);
  _desData->reorder(order);
  // printf("\r\n");
}

void FeatureData::ShrinkLocationData(int ndim, int npoint) {
  if (_locData)
    _locData->shrink(ndim, npoint < 0 ? _locData->npoint()
                                      : min(npoint, _locData->npoint()));
}

}  // namespace visual_sfm
