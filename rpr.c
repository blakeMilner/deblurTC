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
void rpr(float *p1, float *p2, float *pc, long *index, long bs_cnt,
	 float *bsc, float *wc, slaveinfo l)
{
    static long INIT = 0;	// Initialisation flag
    static long *d = NULL;	// quarter circles

    float *bs, *w;		// redundant bispectrum/weights
    long i, j, u, v, r, rs;	// loop variables
    long ir;			// loop variables
    long i1, i2, id1, id2;	// index variable
    long l1, l2, l3, l4;	// limits
    long l2s;			// squared limits
    long nx, ny, nxh, nyh;	// (half of) field size
    float t, t1, t2;		// real temporary variables
    float ct[2], ct1[2], ct2[2], ct3[2];	// complex temporary variables

    nx = l.sfsizex;
    ny = l.sfsizey;
    nxh = nx / 2;
    nyh = ny / 2;
    l2s = l.l2 * l.l2;

    // create bi-spectrum/weights in redundant matrix
    bs_red(index, bs_cnt, l, bsc, wc, &bs, &w);

    // could use improvements later regarding matrix allocation
    if (INIT != l.max_rad) {
		if (d != NULL)
			free(d);
		d = (long *) calloc((l.max_rad + 1) * (l.max_rad + 1), sizeof(long));

		for (r = 0; r <= l.max_rad; r++) {
			rs = r * r;
			// set (quarter-)circle
			for (i = 0; i <= r; i++) {
				d[r * l.max_rad + i] = (long) sqrt((float) (rs - i * i));
			}
		}

		INIT = l.max_rad;
    }
    // begin reconstruction
    //   - move outwards with frequency
    for (ir = 1; ir <= l.max_rad; ir++) {
	//   - move in half-circles
		for (i = 0; i <= ir; i++) {
			idx(i, ir - 1, l.max_rad, id1);
			l1 = (i == ir) ? 0 : d[id1] + 1;
			idx(i, ir, l.max_rad, id2);

			for (j = l1; j <= d[id2]; j++) {
				// if this condition is true then skip
				if (((i == 0) && (j <= 1)) || ((j == 0) && (i <= 1)))
					continue;
				//   - find all frequency values within the above half-circle
				l2 = ir - i - 1;
				for (u = -l2; u <= -i - 1; u++) {
					idx(abs(u), ir - 1, l.max_rad, id1);
					idx(abs(u - i), ir - 1, l.max_rad, id2);
					l3 = (-d[id1] > j - d[id2]) ? -d[id1] : j - d[id2];
					l4 = (d[id1] < j + d[id2]) ? d[id1] : j + d[id2];
					for (v = l3; v <= l4; v++) {
						if ((u * u + v * v) <= l2s) {
							mbs_idx(-u, -v, i, j, l.l1, l.l2, l.max_rad, i1);
							i2 = i1 / 2;
							t1 = w[i2];
							if (t1 < l.snr) {
								ct1[0] = 0.0;
								ct1[1] = 0.0;
							} else {
								ct1[0] = bs[i1];
								ct1[1] = -bs[i1 + 1];
							}
							mbs_idx(-u, v, i, -j, l.l1, l.l2, l.max_rad, i1);
							i2 = i1 / 2;
							t2 = w[i2];
							if (t2 < l.snr) {
								ct2[0] = 0.0;
								ct2[1] = 0.0;
							} else {
								ct2[0] = bs[i1];
								ct2[1] = -bs[i1 + 1];
							}
						} else if (((i - u) * (i - u) + (j - v) * (j - v)) <= l2s) {
							mbs_idx(i - u, j - v, u, v, l.l1, l.l2,	l.max_rad, i1);
							i2 = i1 / 2;
							t1 = w[i2];
							if (t1 < l.snr) {
								ct1[0] = 0.0;
								ct1[1] = 0.0;
							} else {
								ct1[0] = bs[i1];
								ct1[1] = bs[i1 + 1];
							}
							mbs_idx(i - u, -j + v, u, -v, l.l1, l.l2, l.max_rad, i1);
							i2 = i1 / 2;
							t2 = w[i2];
							if (t2 < l.snr) {
								ct2[0] = 0.0;
								ct2[1] = 0.0;
							} else {
								ct2[0] = bs[i1];
								ct2[1] = bs[i1 + 1];
							}
						} else
							continue;

						phs_idx(nxh - u, nyh + v, nx, i1);
						phs_idx(nxh - i + u, nyh + j - v, nx, i2);
						cmul(&p1[i1], &p1[i2], ct3);
						cdiv(ct1, ct3, ct);
						conj(ct, ct);
						phs_idx(nxh - i, nyh + j, nx, i1);
						t = pow(t1, l.eps);
						p1[i1] += t * ct[0];
						p1[i1 + 1] += t * ct[1];
						pc[i1 / 2] += t;
						p2[i1] += t * ct[0];
						p2[i1 + 1] += t * ct[1];

						phs_idx(nxh - u, nyh - v, nx, i1);
						phs_idx(nxh - i + u, nyh - j + v, nx, i2);
						cmul(&p1[i1], &p1[i2], ct3);
						cdiv(ct2, ct3, ct);
						conj(ct, ct);
						phs_idx(nxh - i, nyh - j, nx, i1);
						t = pow(t2, l.eps);
						p1[i1] += t * ct[0];
						p1[i1 + 1] += t * ct[1];
						pc[i1 / 2] += t;
						p2[i1] += t * ct[0];
						p2[i1 + 1] += t * ct[1];
					}
				}
				l2 = (i < ir - i - 1) ? i : ir - i - 1;
				for (u = -l2; u <= -1; u++) {
					idx(abs(u), ir - 1, l.max_rad, id1);
					idx(abs(u - i), ir - 1, l.max_rad, id2);
					l3 = (-d[id1] > j - d[id2]) ? -d[id1] : j - d[id2];
					l4 = (d[id1] < j + d[id2]) ? d[id1] : j + d[id2];
					if (-u == i)
						l4 = (l4 < j - 1) ? l4 : j - 1;
						for (v = l3; v <= l4; v++) {
							if ((u * u + v * v) <= l2s) {
								mbs_idx(-u, -v, i, j, l.l1, l.l2, l.max_rad,
									i1);
								i2 = i1 / 2;
								t1 = w[i2];
								if (t1 < l.snr) {
								ct1[0] = 0.0;
								ct1[1] = 0.0;
								} else {
								ct1[0] = bs[i1];
								ct1[1] = -bs[i1 + 1];
								}
								mbs_idx(-u, v, i, -j, l.l1, l.l2, l.max_rad,
									i1);
								i2 = i1 / 2;
								t2 = w[i2];
								if (t2 < l.snr) {
								ct2[0] = 0.0;
								ct2[1] = 0.0;
								} else {
								ct2[0] = bs[i1];
								ct2[1] = -bs[i1 + 1];
								}
							}
							else if (((i - u) * (i - u) + (j - v) * (j - v)) <= l2s) {
								mbs_idx(i - u, j - v, u, v, l.l1, l.l2,
									l.max_rad, i1);
								i2 = i1 / 2;
								t1 = w[i2];
								if (t1 < l.snr) {
								ct1[0] = 0.0;
								ct1[1] = 0.0;
								} else {
								ct1[0] = bs[i1];
								ct1[1] = bs[i1 + 1];
								}
								mbs_idx(i - u, -j + v, u, -v, l.l1, l.l2,
									l.max_rad, i1);
								i2 = i1 / 2;
								t2 = w[i2];
								if (t2 < l.snr) {
								ct2[0] = 0.0;
								ct2[1] = 0.0;
								} else {
								ct2[0] = bs[i1];
								ct2[1] = bs[i1 + 1];
							}
						}
						else
							continue;

						phs_idx(nxh - u, nyh + v, nx, i1);
						phs_idx(nxh - i + u, nyh + j - v, nx, i2);
						cmul(&p1[i1], &p1[i2], ct3);
						cdiv(ct1, ct3, ct);
						conj(ct, ct);
						phs_idx(nxh - i, nyh + j, nx, i1);
						t = pow(t1, l.eps);
						p1[i1] += t * ct[0];
						p1[i1 + 1] += t * ct[1];
						pc[i1 / 2] += t;
						p2[i1] += t * ct[0];
						p2[i1 + 1] += t * ct[1];

						phs_idx(nxh - u, nyh - v, nx, i1);
						phs_idx(nxh - i + u, nyh - j + v, nx, i2);
						cmul(&p1[i1], &p1[i2], ct3);
						cdiv(ct2, ct3, ct);
						conj(ct, ct);
						phs_idx(nxh - i, nyh - j, nx, i1);
						t = pow(t2, l.eps);
						p1[i1] += t * ct[0];
						p1[i1 + 1] += t * ct[1];
						pc[i1 / 2] += t;
						p2[i1] += t * ct[0];
						p2[i1 + 1] += t * ct[1];
					}
				}

				l2 = (l.l1 <
					  l.l2) ? ((l.l1 < (long) i / 2) ? l.l1 : (long) i / 2)
					: ((l.l2 < (long) i / 2) ? l.l2 : (long) i / 2);
				for (u = 0; u <= l2; u++) {
					idx(abs(u), ir, l.max_rad, id1);
					idx(abs(u - i), ir, l.max_rad, id2);

					l3 = (-l.l1 > -d[id1]) ?
					((-l.l1 > j - d[id2]) ? -l.l1 : j - d[id2])
					: ((-d[id1] > j - d[id2]) ? -d[id1] : j - d[id2]);

					l4 = (l.l1 < d[id1]) ?
					((l.l1 < j + d[id2]) ? l.l1 : j + d[id2])
					: ((d[id1] < j + d[id2]) ? d[id1] : j + d[id2]);

					if (u == 0) {
						l3 = 1;
						l4 = (l4 < 2 * j - 1) ? l4 : 2 * j - 1;
					}

					if (u == (i - u)) {
						l4 = (l4 < (long) j / 2) ? l4 : (long) j / 2;
					}

					for (v = l3; v <= l4; v++) {
						if ((u * u + v * v) > l2s)
							continue;
						mbs_idx(u, v, i - u, j - v, l.l1, l.l2, l.max_rad,
							i1);
						i2 = i1 / 2;
						t1 = w[i2];
						if (t1 < l.snr) {
							ct1[0] = 0.0;
							ct1[1] = 0.0;
						} else {
							ct1[0] = bs[i1];
							ct1[1] = bs[i1 + 1];
						}
						mbs_idx(u, -v, i - u, -j + v, l.l1, l.l2,
							l.max_rad, i1);
						i2 = i1 / 2;
						t2 = w[i2];
						if (t2 < l.snr) {
							ct2[0] = 0.0;
							ct2[1] = 0.0;
						} else {
							ct2[0] = bs[i1];
							ct2[1] = bs[i1 + 1];
						}
						phs_idx(nxh - u, nyh + v, nx, i1);
						phs_idx(nxh - i + u, nyh + j - v, nx, i2);
						cmul(&p1[i1], &p1[i2], ct3);
						cdiv(ct1, ct3, ct);
						conj(ct, ct);
						phs_idx(nxh - i, nyh + j, nx, i1);
						t = pow(t1, l.eps);
						p1[i1] += t * ct[0];
						p1[i1 + 1] += t * ct[1];
						pc[i1 / 2] += t;
						p2[i1] += t * ct[0];
						p2[i1 + 1] += t * ct[1];

						phs_idx(nxh - u, nyh - v, nx, i1);
						phs_idx(nxh - i + u, nyh - j + v, nx, i2);
						cmul(&p1[i1], &p1[i2], ct3);
						cdiv(ct2, ct3, ct);
						conj(ct, ct);
						phs_idx(nxh - i, nyh - j, nx, i1);
						t = pow(t2, l.eps);
						p1[i1] += t * ct[0];
						p1[i1 + 1] += t * ct[1];
						pc[i1 / 2] += t;
						p2[i1] += t * ct[0];
						p2[i1 + 1] += t * ct[1];
					}
				}

				phs_idx(nxh - i, nyh - j, nx, i1);
				t1 = cmod(&p1[i1]);

				if (t1 != 0) {
					p1[i1] /= t1;
					p1[i1 + 1] /= t1;
				} else {
					p1[i1] = 1.0;
					p1[i1 + 1] = 0.0;
				}

				phs_idx(nxh + i, nyh + j, nx, i2);
				p1[i2] = p1[i1];
				p1[i2 + 1] = -p1[i1 + 1];
				pc[i2 / 2] = pc[i1 / 2];
				p2[i2] = p2[i1];
				p2[i2 + 1] = -p2[i1 + 1];

				phs_idx(nxh - i, nyh + j, nx, i1);
				t1 = cmod(&p1[i1]);

				if (t1 != 0) {
					p1[i1] /= t1;
					p1[i1 + 1] /= t1;
				} else {
					p1[i1] = 1.0;
					p1[i1 + 1] = 0.0;
				}

				phs_idx(nxh + i, nyh - j, nx, i2);
				p1[i2] = p1[i1];
				p1[i2 + 1] = -p1[i1 + 1];
				pc[i2 / 2] = pc[i1 / 2];
				p2[i2] = p2[i1];
				p2[i2 + 1] = -p2[i1 + 1];
			}
		}
    }

    for (i = 0; i < nx * ny; i++) {
		t1 = cmod(&p2[2 * i]);
		pc[i] = ((t1 != 0) && (pc[i] != 0.0)) ? t1 / pc[i] : 0.0;
    }

    // reset temp matrix
    memset(p2, 0, 2 * nx * ny * sizeof(float));
    // free mem
    free(bs);
    free(w);
}
