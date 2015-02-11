/*
    oast9 - OAST, an optimal corner detector based on the
            accelerated segment test for a 16 pixel mask

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


#ifndef OAST9_16_H
#define OAST9_16_H

#include <stdint.h>
#include "AstDetector.h"

struct OpenCVPoint;

namespace agast{

	class OastDetector9_16 : public AstDetector
	{
		public:
			OastDetector9_16():AstDetector(){;}
			OastDetector9_16(int width, int height, int thr):AstDetector(width, height, thr){init_pattern();};
			~OastDetector9_16(){}
			void detect(const unsigned char* im,
					std::vector<OpenCVPoint>& keypoints);
			void nms(const unsigned char* im,
					const std::vector<OpenCVPoint>& keypoints, std::vector<OpenCVPoint>& keypoints_nms);
			int get_borderWidth(){return borderWidth;}
			int cornerScore(const unsigned char* p);

		private:
			static const int borderWidth=3;
			int_fast16_t s_offset0;
			int_fast16_t s_offset1;
			int_fast16_t s_offset2;
			int_fast16_t s_offset3;
			int_fast16_t s_offset4;
			int_fast16_t s_offset5;
			int_fast16_t s_offset6;
			int_fast16_t s_offset7;
			int_fast16_t s_offset8;
			int_fast16_t s_offset9;
			int_fast16_t s_offset10;
			int_fast16_t s_offset11;
			int_fast16_t s_offset12;
			int_fast16_t s_offset13;
			int_fast16_t s_offset14;
			int_fast16_t s_offset15;

			void init_pattern()
			{
				s_offset0=(-3)+(0)*xsize;
				s_offset1=(-3)+(-1)*xsize;
				s_offset2=(-2)+(-2)*xsize;
				s_offset3=(-1)+(-3)*xsize;
				s_offset4=(0)+(-3)*xsize;
				s_offset5=(1)+(-3)*xsize;
				s_offset6=(2)+(-2)*xsize;
				s_offset7=(3)+(-1)*xsize;
				s_offset8=(3)+(0)*xsize;
				s_offset9=(3)+(1)*xsize;
				s_offset10=(2)+(2)*xsize;
				s_offset11=(1)+(3)*xsize;
				s_offset12=(0)+(3)*xsize;
				s_offset13=(-1)+(3)*xsize;
				s_offset14=(-2)+(2)*xsize;
				s_offset15=(-3)+(1)*xsize;
			}
	};

}

#endif /* OAST9_16_H */
