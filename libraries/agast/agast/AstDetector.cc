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

#include "AstDetector.h"
#include "cvWrapper.h"

using namespace std;
using namespace agast;

void AstDetector::score(const unsigned char* i, const std::vector<OpenCVPoint>& corners_all)
{
	unsigned int n=0;
    unsigned int num_corners=corners_all.size();

	if(num_corners > scores.capacity())
	{
		if(scores.capacity()==0)
		{
			scores.reserve(512 > num_corners ? 512 : num_corners);
		}
		else
		{
			unsigned int nScores = scores.capacity()*2;
			if(num_corners > nScores)
				nScores = num_corners;
			scores.reserve(nScores);
		}
	}

    scores.resize(num_corners);

    for(; n < num_corners; n++)
        scores[n] = cornerScore(i + corners_all[n].y*xsize + corners_all[n].x);
}

void AstDetector::nms(const unsigned char* im, const std::vector<OpenCVPoint>& corners_all,
		std::vector<OpenCVPoint>& corners_nms)
{
	score(im,corners_all);
	nonMaximumSuppression(corners_all, corners_nms);
}
