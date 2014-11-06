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
** core speckle routines
**
** =====================================================================
**
** uses:  fftw3 library:
**   http://www.fftw.org
**
**   Library for fast fourier transform
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

#include "speckle_core.h"


/****************************************************************************
***************TRIPLE CORRELATION PART FOR PHASE RECONSTRUCTION**************
****************************************************************************/
void bs_red(long *index, long bs_cnt, slaveinfo l, float *bsc, float *wc,
	    float **bs, float **w)
{
    long u1, u2, v1, v2;	// bs vector indices
    long cnt;			// loop variables
    long i1, i2;		// index variable
    long l2s;			// squared limit

    l2s = l.l2 * l.l2;

    (*bs) =
	(float *) calloc(2 * (l.l1 + 1) * (2 * (l.l2 + 1)) *
			 (2 * (l.max_rad + 1)) * (2 * (l.max_rad + 1)),
			 sizeof(float));
    (*w) =
	(float *) calloc((l.l1 + 1) * (2 * (l.l2 + 1)) *
			 (2 * (l.max_rad + 1)) * (2 * (l.max_rad + 1)),
			 sizeof(float));

    // create redundant matrices (eq. readin.f)
    for (cnt = 0; cnt < bs_cnt; cnt++) {
		// get vectors
		u1 = index[0 + cnt * 4];
		u2 = index[1 + cnt * 4];
		v1 = index[2 + cnt * 4];
		v2 = index[3 + cnt * 4];
		// use symmetries to fill bispectrum and weights matrices
		if ((u1 * u1 + u2 * u2) <= l2s) {
			mbs_idx(u1, u2, v1, v2, l.l1, l.l2, l.max_rad, i1);
			i2 = i1 / 2;
			(*bs)[i1] = bsc[2 * cnt];
			(*bs)[i1 + 1] = bsc[2 * cnt + 1];
			(*w)[i2] = wc[cnt];
		}
		if ((v1 * v1 + v2 * v2) <= l2s) {
			mbs_idx(v1, v2, u1, u2, l.l1, l.l2, l.max_rad, i1);
			i2 = i1 / 2;
			(*bs)[i1] = bsc[2 * cnt];
			(*bs)[i1 + 1] = bsc[2 * cnt + 1];
			(*w)[i2] = wc[cnt];
		}
		if (((u1 + v1) * (u1 + v1) + (u2 + v2) * (u2 + v2)) <= l2s) {
			mbs_idx(u1 + v1, u2 + v2, -u1, -u2, l.l1, l.l2, l.max_rad, i1);
			i2 = i1 / 2;
			(*bs)[i1] = bsc[2 * cnt];
			(*bs)[i1 + 1] = -bsc[2 * cnt + 1];
			(*w)[i2] = wc[cnt];
			mbs_idx(u1 + v1, u2 + v2, -v1, -v2, l.l1, l.l2, l.max_rad, i1);
			i2 = i1 / 2;
			(*bs)[i1] = bsc[2 * cnt];
			(*bs)[i1 + 1] = -bsc[2 * cnt + 1];
			(*w)[i2] = wc[cnt];
		}
    }
}
