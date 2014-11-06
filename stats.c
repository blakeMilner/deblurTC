/* $Id$
 *
 *  (C) 2003-2007, Friedrich Woeger <woeger@kis.uni-freiburg.de>,
 *  Kiepenheuer-Institut für Sonnenphysik, Freiburg (Breisgau), Germany
 *  All rights reserved.
 *  
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are
 *  met:
 *  
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above copyright
 *     notice, this list of conditions and the following disclaimer in
 *     the documentation and/or other materials provided with the
 *     distribution.
 *   * Neither the name of the Kiepenheuer-Institut nor the
 *     names of its contributors may be used to endorse or promote
 *     products derived from this software without specific prior
 *     written permission.
 *  
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 *  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 *  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 *  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 *  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 *  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 *  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 *  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 */


/* =====================================================================
**
**	Supplementary speckle math routines
**
** =====================================================================
**
**	    uses fftw3 library (http://www.fftw.org)
**	    Library for fast fourier transform
**
** =====================================================================
**
** Author:  F. Woeger
**   Kiepenheuer-Institut fuer Sonnenphysik
**   Freiburg, Germany
**
** Written 17. 06. 2003
**
** =====================================================================
*/


#include "speckle_math.h"


void mean(float *in, int nx, int ny, int nfr, float *out)
{
    // Declare helper variables
    long i, l;
    long N = nx * ny;
    // set resultant matrix to zero!
    memset(out, 0.0, N * sizeof(float));

    // calculate the mean image
    for (i = 0; i < N; i++) {
		for (l = 0; l < nfr; l++) {
			out[i] += in[i + l * N];
		}
	out[i] /= (float) nfr;
    }
}


void stats(float *in, long N, float thres, double *out1, double *out2)
{
    /* Declare helper variable */
    long i, t, c = 0;

    if (out2 != NULL) {
	*out1 = 0;
	*out2 = 0;
	for (i = 0; i < N; i++) {
	    *out1 += in[i];
	    *out2 += in[i] * in[i];
	}
	/* calculate the mean value of image */
	*out1 /= N;
	/* calculate the variance of image */
	*out2 = *out2 / N - (*out1) * (*out1);
    } else {
	*out1 = 0;
	for (i = 0; i < N; i++) {
	    *out1 += in[i];
	}
	/* calculate the mean value of image */
	*out1 /= N;
    }
    /*
     *  if thres greater zero then compute mean only from pixels that are
     *  greater than thres percent of overall mean
     */
    if (thres > 0) {
	t = thres * (*out1);
	*out1 = 0;
	for (i = 0; i < N; i++) {
	    if (in[i] > t) {
		*out1 += in[i];
		c++;
	    }
	}
	*out1 /= c;
    }
}
