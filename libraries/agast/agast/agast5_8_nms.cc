/*
    agast5 - AGAST, an adaptive and generic corner detector based on the
             accelerated segment test for a 8 pixel mask

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


#include <stdint.h>
#include <stdlib.h>
#include "cvWrapper.h"
#include "agast5_8.h"

using namespace std;
using namespace agast;

//using also bisection as propsed by Edward Rosten in FAST,
//but it is based on the OAST
int AgastDetector5_8::cornerScore(const unsigned char* p)
{
    int bmin = b;
    int bmax = 255;
    int b_test = (bmax + bmin)/2;

    register int_fast16_t offset0=s_offset0;
	register int_fast16_t offset1=s_offset1;
	register int_fast16_t offset2=s_offset2;
	register int_fast16_t offset3=s_offset3;
	register int_fast16_t offset4=s_offset4;
	register int_fast16_t offset5=s_offset5;
	register int_fast16_t offset6=s_offset6;
	register int_fast16_t offset7=s_offset7;

	while(1)
	{
		register const int cb = *p + b_test;
		register const int c_b = *p - b_test;
		if(p[offset0] > cb)
		  if(p[offset2] > cb)
		    if(p[offset3] > cb)
		      if(p[offset5] > cb)
		        if(p[offset1] > cb)
		          if(p[offset4] > cb)
		            goto is_a_corner;
		          else
		            if(p[offset7] > cb)
		              goto is_a_corner;
		            else
		              goto is_not_a_corner;
		        else
		          if(p[offset4] > cb)
		            if(p[offset6] > cb)
		              goto is_a_corner;
		            else
		              goto is_not_a_corner;
		          else
		            goto is_not_a_corner;
		      else
		        if(p[offset1] > cb)
		          if(p[offset4] > cb)
		            goto is_a_corner;
		          else
		            if(p[offset7] > cb)
		              goto is_a_corner;
		            else
		              goto is_not_a_corner;
		        else
		          goto is_not_a_corner;
		    else
		      if(p[offset7] > cb)
		        if(p[offset6] > cb)
		          if(p[offset5] > cb)
		            if(p[offset1] > cb)
		              goto is_a_corner;
		            else
		              if(p[offset4] > cb)
		                goto is_a_corner;
		              else
		                goto is_not_a_corner;
		          else
		            if(p[offset1] > cb)
		              goto is_a_corner;
		            else
		              goto is_not_a_corner;
		        else
		          goto is_not_a_corner;
		      else
		        if(p[offset5] < c_b)
		          if(p[offset3] < c_b)
		            if(p[offset7] < c_b)
		              if(p[offset4] < c_b)
		                if(p[offset6] < c_b)
		                  goto is_a_corner;
		                else
		                  goto is_not_a_corner;
		              else
		                goto is_not_a_corner;
		            else
		              goto is_not_a_corner;
		          else
		            goto is_not_a_corner;
		        else
		          goto is_not_a_corner;
		  else
		    if(p[offset5] > cb)
		      if(p[offset7] > cb)
		        if(p[offset6] > cb)
		          if(p[offset1] > cb)
		            goto is_a_corner;
		          else
		            if(p[offset4] > cb)
		              goto is_a_corner;
		            else
		              goto is_not_a_corner;
		        else
		          goto is_not_a_corner;
		      else
		        goto is_not_a_corner;
		    else
		      if(p[offset5] < c_b)
		        if(p[offset3] < c_b)
		          if(p[offset2] < c_b)
		            if(p[offset1] < c_b)
		              if(p[offset4] < c_b)
		                goto is_a_corner;
		              else
		                goto is_not_a_corner;
		            else
		              if(p[offset4] < c_b)
		                if(p[offset6] < c_b)
		                  goto is_a_corner;
		                else
		                  goto is_not_a_corner;
		              else
		                goto is_not_a_corner;
		          else
		            if(p[offset7] < c_b)
		              if(p[offset4] < c_b)
		                if(p[offset6] < c_b)
		                  goto is_a_corner;
		                else
		                  goto is_not_a_corner;
		              else
		                goto is_not_a_corner;
		            else
		              goto is_not_a_corner;
		        else
		          goto is_not_a_corner;
		      else
		        goto is_not_a_corner;
		else if(p[offset0] < c_b)
		  if(p[offset2] < c_b)
		    if(p[offset7] > cb)
		      if(p[offset3] < c_b)
		        if(p[offset5] < c_b)
		          if(p[offset1] < c_b)
		            if(p[offset4] < c_b)
		              goto is_a_corner;
		            else
		              goto is_not_a_corner;
		          else
		            if(p[offset4] < c_b)
		              if(p[offset6] < c_b)
		                goto is_a_corner;
		              else
		                goto is_not_a_corner;
		            else
		              goto is_not_a_corner;
		        else
		          if(p[offset1] < c_b)
		            if(p[offset4] < c_b)
		              goto is_a_corner;
		            else
		              goto is_not_a_corner;
		          else
		            goto is_not_a_corner;
		      else
		        if(p[offset5] > cb)
		          if(p[offset3] > cb)
		            if(p[offset4] > cb)
		              if(p[offset6] > cb)
		                goto is_a_corner;
		              else
		                goto is_not_a_corner;
		            else
		              goto is_not_a_corner;
		          else
		            goto is_not_a_corner;
		        else
		          goto is_not_a_corner;
		    else
		      if(p[offset7] < c_b)
		        if(p[offset3] < c_b)
		          if(p[offset5] < c_b)
		            if(p[offset1] < c_b)
		              goto is_a_corner;
		            else
		              if(p[offset4] < c_b)
		                if(p[offset6] < c_b)
		                  goto is_a_corner;
		                else
		                  goto is_not_a_corner;
		              else
		                goto is_not_a_corner;
		          else
		            if(p[offset1] < c_b)
		              goto is_a_corner;
		            else
		              goto is_not_a_corner;
		        else
		          if(p[offset6] < c_b)
		            if(p[offset5] < c_b)
		              if(p[offset1] < c_b)
		                goto is_a_corner;
		              else
		                if(p[offset4] < c_b)
		                  goto is_a_corner;
		                else
		                  goto is_not_a_corner;
		            else
		              if(p[offset1] < c_b)
		                goto is_a_corner;
		              else
		                goto is_not_a_corner;
		          else
		            goto is_not_a_corner;
		      else
		        if(p[offset3] < c_b)
		          if(p[offset5] < c_b)
		            if(p[offset1] < c_b)
		              if(p[offset4] < c_b)
		                goto is_a_corner;
		              else
		                goto is_not_a_corner;
		            else
		              if(p[offset4] < c_b)
		                if(p[offset6] < c_b)
		                  goto is_a_corner;
		                else
		                  goto is_not_a_corner;
		              else
		                goto is_not_a_corner;
		          else
		            if(p[offset1] < c_b)
		              if(p[offset4] < c_b)
		                goto is_a_corner;
		              else
		                goto is_not_a_corner;
		            else
		              goto is_not_a_corner;
		        else
		          goto is_not_a_corner;
		  else
		    if(p[offset5] > cb)
		      if(p[offset3] > cb)
		        if(p[offset2] > cb)
		          if(p[offset1] > cb)
		            if(p[offset4] > cb)
		              goto is_a_corner;
		            else
		              goto is_not_a_corner;
		          else
		            if(p[offset4] > cb)
		              if(p[offset6] > cb)
		                goto is_a_corner;
		              else
		                goto is_not_a_corner;
		            else
		              goto is_not_a_corner;
		        else
		          if(p[offset7] > cb)
		            if(p[offset4] > cb)
		              if(p[offset6] > cb)
		                goto is_a_corner;
		              else
		                goto is_not_a_corner;
		            else
		              goto is_not_a_corner;
		          else
		            goto is_not_a_corner;
		      else
		        goto is_not_a_corner;
		    else
		      if(p[offset5] < c_b)
		        if(p[offset7] < c_b)
		          if(p[offset6] < c_b)
		            if(p[offset1] < c_b)
		              goto is_a_corner;
		            else
		              if(p[offset4] < c_b)
		                goto is_a_corner;
		              else
		                goto is_not_a_corner;
		          else
		            goto is_not_a_corner;
		        else
		          goto is_not_a_corner;
		      else
		        goto is_not_a_corner;
		else
		  if(p[offset3] > cb)
		    if(p[offset5] > cb)
		      if(p[offset2] > cb)
		        if(p[offset1] > cb)
		          if(p[offset4] > cb)
		            goto is_a_corner;
		          else
		            goto is_not_a_corner;
		        else
		          if(p[offset4] > cb)
		            if(p[offset6] > cb)
		              goto is_a_corner;
		            else
		              goto is_not_a_corner;
		          else
		            goto is_not_a_corner;
		      else
		        if(p[offset7] > cb)
		          if(p[offset4] > cb)
		            if(p[offset6] > cb)
		              goto is_a_corner;
		            else
		              goto is_not_a_corner;
		          else
		            goto is_not_a_corner;
		        else
		          goto is_not_a_corner;
		    else
		      goto is_not_a_corner;
		  else
		    if(p[offset3] < c_b)
		      if(p[offset5] < c_b)
		        if(p[offset2] < c_b)
		          if(p[offset1] < c_b)
		            if(p[offset4] < c_b)
		              goto is_a_corner;
		            else
		              goto is_not_a_corner;
		          else
		            if(p[offset4] < c_b)
		              if(p[offset6] < c_b)
		                goto is_a_corner;
		              else
		                goto is_not_a_corner;
		            else
		              goto is_not_a_corner;
		        else
		          if(p[offset7] < c_b)
		            if(p[offset4] < c_b)
		              if(p[offset6] < c_b)
		                goto is_a_corner;
		              else
		                goto is_not_a_corner;
		            else
		              goto is_not_a_corner;
		          else
		            goto is_not_a_corner;
		      else
		        goto is_not_a_corner;
		    else
		      goto is_not_a_corner;

		is_a_corner:
			bmin=b_test;
			goto end;

		is_not_a_corner:
			bmax=b_test;
			goto end;

		end:

		if(bmin == bmax - 1 || bmin == bmax)
			return bmin;
		b_test = (bmin + bmax) / 2;
	}
}


