/*
  AstDetector - the interface class for the AGAST corner detector

  Copyright (c) 2010, Elmar Mair
  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are met:
  * Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.
  * Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer in the
  documentation and/or other materials provided with the distribution.
  * Neither the name of the owner nor the names of its contributors may be
  used to endorse or promote products derived from this software without
  specific prior written permission.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS BE LIABLE FOR ANY
  DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/


#ifndef ASTDETECTOR_H
#define ASTDETECTOR_H

#include <vector>
#include <iostream>
#include "cvWrapper.h"

namespace agast{

class AstDetector
{
 public:
  AstDetector():xsize(0),ysize(0),b(-1) {}
  AstDetector(int width, int height, int thr):xsize(width),ysize(height),b(thr) {}
  virtual ~AstDetector(){;}
  virtual void detect(const unsigned char* im, std::vector<OpenCVPoint>& corners_all)=0;
  virtual int get_borderWidth()=0;
  void nms(const unsigned char* im,
           const std::vector<OpenCVPoint>& corners_all, std::vector<OpenCVPoint>& corners_nms);
  void processImage(const unsigned char* im,
                    std::vector<OpenCVPoint>& keypoints_nms) {
    std::vector<OpenCVPoint> keypoints;
    detect(im,keypoints);
    nms(im,keypoints,keypoints_nms);}
  void set_threshold(int b_){b=b_;}
  void set_imageSize(int xsize_, int ysize_){xsize=xsize_; ysize=ysize_; init_pattern();}
  virtual int cornerScore(const unsigned char* p)=0;

 protected:
  virtual void init_pattern()=0;
  void score(const unsigned char* i, const std::vector<OpenCVPoint>& corners_all);
  void nonMaximumSuppression(const std::vector<OpenCVPoint>& corners_all,
                             std::vector<OpenCVPoint>& corners_nms);
  std::vector<int> scores;
  std::vector<int> nmsFlags;
  int xsize, ysize;
  int b;
};

}

#endif /* AGASTDETECTOR_H */
