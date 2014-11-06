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
void iwlspr(float *p1, float *p2, float *pc, float *bsc, float *wc,
	    long *index, long bs_cnt, slaveinfo l, int maxk)
{
    static long INIT = 0;	// initialisation flag
    static long *d = NULL;	// circle with radius l.bs_mr
    static float *pct = NULL;	// temporary phase consistency

    long i, j, cnt;		// loop variable
    long i1, i2, i3;		// index variable
    long nx, ny, nxh, nyh;	// (half of) field size
    long u1, u2, v1, v2;	// bispectrum vectors
    float ct[2], ct1[2], ct2[2];	// complex temporary variable
    float w1;			// real temporary variable
    float t1;		// real temporary variable
    long c1, c2, c3, c4;	// condition flags

    nx = l.sfsizex;
    ny = l.sfsizey;
    nxh = nx / 2;
    nyh = ny / 2;

    // set (quarter-)circle
    if (INIT != l.max_rad) {
		if (d != NULL)
			free(d);
		if (pct != NULL)
			free(pct);
		d = (long *) malloc((l.max_rad + 1) * sizeof(long));
		for (i = 0; i <= l.max_rad; i++) {
			d[i] = (long) sqrt((float) (l.max_rad * l.max_rad - i * i));
		}
		// setup temporary phase consistency
		pct = (float *) malloc(nx * ny * sizeof(float));
		INIT = l.max_rad;
    }
    // reset matrix
    memset(pct, 0.0, nx * ny * sizeof(float));

    for (cnt = 0; cnt < bs_cnt; cnt++) {
		// get vectors
		u1 = index[0 + cnt * 4];
		u2 = index[1 + cnt * 4];
		v1 = index[2 + cnt * 4];
		v2 = index[3 + cnt * 4];

		// don't go along the negative y axis with u
		if ((u1 == 0) && (u2 <= 0))
			continue;

		if (wc[cnt] >= l.snr) {
			ct[0] = bsc[2 * cnt];
			ct[1] = bsc[2 * cnt + 1];
		} else {
			ct[0] = 0;
			ct[1] = 0;
		}
		w1 = pow(wc[cnt], l.eps);

		// set condition flags
		//   the first three are there to not alter the initial guesses around the origin
		if ((u1 * u1 + u2 * u2) <= 1)
			c1 = 1;
		else
			c1 = 0;
		if ((v1 * v1 + v2 * v2) <= 1)
			c2 = 1;
		else
			c2 = 0;
		if (((u1 + v1) * (u1 + v1) + (u2 + v2) + (u2 + v2)) <= 1)
			c3 = 1;
		else
			c3 = 0;
		//  the last one is for fitting theory
		if ((u1 == v1) && (u2 == v2))
			c4 = 1;
		else
			c4 = 0;
		// compute indices for phases only once
		phs_idx(nxh - v1, nyh + v2, nx, i1);	// v position
		phs_idx(nxh - u1, nyh + u2, nx, i2);	// u position
		phs_idx(nxh - u1 - v1, nyh + u2 + v2, nx, i3);	// u+v position
		// finally calculate the values
		if (!c4) {
			if (!c1) {
				// compute weight
				t1 = pc[i1 / 2] * pc[i3 / 2] * w1;
				// phase at u
				conj(&p1[i3], ct1);
				cmul(&p1[i1], ct1, ct2);
				cdiv(ct, ct2, ct1);
				ct2[0] = t1 * ct1[0];
				ct2[1] = t1 * ct1[1];
				cadd(&p2[i2], ct2, &p2[i2]);
				pct[i2 / 2] += t1;
			}
			if (!c2) {
				// compute weight
				t1 = pc[i2 / 2] * pc[i3 / 2] * w1;
				// phase at v
				conj(&p1[i3], ct1);
				cmul(&p1[i2], ct1, ct2);
				cdiv(ct, ct2, ct1);
				ct2[0] = t1 * ct1[0];
				ct2[1] = t1 * ct1[1];
				cadd(&p2[i1], ct2, &p2[i1]);
				pct[i1 / 2] += t1;
			}
		} else {
			// compute weight
			t1 = 4.0 * pc[i3 / 2] * w1;
			// phase at u=v
			conj(&p1[i3], ct1);
			cdiv(ct, ct1, ct1);
			csqr(ct1, ct2);
			ct2[0] *= t1;
			ct2[1] *= t1;
			cadd(&p2[i2], ct2, &p2[i2]);
			pct[i2 / 2] += t1;
		}
		if (!c3) {
			// compute weight
			t1 = pc[i1 / 2] * pc[i2 / 2] * w1;
			// phase at u+v
			cmul(&p1[i1], &p1[i2], ct1);
			cdiv(ct, ct1, ct1);
			conj(ct1, ct2);
			ct2[0] *= t1;
			ct2[1] *= t1;
			cadd(&p2[i3], ct2, &p2[i3]);
			pct[i3 / 2] += t1;
		}
    }
// assign values to output phase p1
//for (i=0; i<=l.max_rad; i++) {
    for (i = maxk; i <= l.max_rad; i++) {
		i3 = (i == 0) ? 1 : -d[i];
		for (j = i3; j < d[i]; j++) {
			phs_idx(nxh - i, nyh + j, nx, i1);
			phs_idx(nxh + i, nyh - j, nx, i2);
	/*
				if ((p2[i1] == 0) && (p2[i1+1] == 0)) {
				  p1[i1] = 1.0; p1[i1+1] = 0.0;
				  p1[i2] = 1.0; p1[i2+1] = 0.0;
				} else {
	 */
			if ((p2[i1] != 0) && (p2[i1 + 1] != 0)) {
				t1 = cmod(&p2[i1]);
				// normalize
				p1[i1] = p2[i1] / t1;
				p1[i1 + 1] = p2[i1 + 1] / t1;
				pc[i1 / 2] = (pct[i1 / 2] != 0) ? t1 / pct[i1 / 2] : 0.0;
				// mirror
				conj(&p1[i1], &p1[i2]);
				pc[i2 / 2] = pc[i1 / 2];
			}
		}
    }
    // reset temp matrix
    memset(p2, 0, 2 * nx * ny * sizeof(float));
}
