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
/*
 *  bs_init: initialise the bispectrum matrices
 *  1. determine the indices for the bispectrum elements to be used
 *    note: this depends on
 *    a) the diffraction limit in px (or a percentage thereof)
 *    b) the limit entered for the size of u in x direction
 *    c) the limit entered for the length of u, v and u+v
 *  2. determine the overall number of elements in bispectrum
 *  3. return the index matrix
 *
 *    bs_cnt: # of elements in bispectrum (return value)
 *    bs_mr: see 1. a)
 *    bs_l1: see 1. b)
 *    bs_l2: see 1. c)
 */

long *bs_init(long *bs_cnt, slaveinfo l)
{
    long i;			// loop variable
    long limit;			// helper variable
    long mrs, l2s;		// squared limits
    long u1, u2, v1, v2, uv1, uv2;	// u, v and u+v vector components
    long u1s, u2s, v1s, v2s, uv1s, uv2s;	// squared vector components
    long *d;			// circle with radius l.bs_mr
    long *index = NULL;		// returned matrix

    d = (long *) malloc((l.max_rad + 1) * sizeof(long));

    mrs = l.max_rad * l.max_rad;
    l2s = l.l2 * l.l2;

    // set (quarter-)circle
    for (i = 0; i <= l.max_rad; i++) {
    	d[i] = (long) sqrt((float) (mrs - i * i));
    }

    /*
     *  Define the parts of the bispectrum to be evaluated
     *  Store the indices in a one dimensional array for later use by the
     *    bispectrum algorithm
     */

    i = 0;
    for (u1 = 0; u1 <= l.l1; u1++) {
		u1s = u1 * u1;

		for (u2 = -d[u1]; u2 <= d[u1]; u2++) {
			u2s = u2 * u2;
			for (v1 = u1; v1 <= l.max_rad; v1++) {
				v1s = v1 * v1;
				uv1 = u1 + v1;
				uv1s = uv1 * uv1;
				// set negative limit for v2
				if ((u1 == v1) && (u2 != 0))
					limit = -abs(u2);
				else
					limit = -1;
				for (v2 = -l.max_rad; v2 <= limit; v2++) {
					v2s = v2 * v2;
					uv2 = u2 + v2;
					uv2s = uv2 * uv2;
					if (((uv1s + uv2s) <= mrs) &&
					((v1s + v2s) <= mrs) &&
					(((u1s + u2s) <= l2s) ||
					 ((v1s + v2s) <= l2s) || ((uv1s + uv2s) <= l2s))) {
						// reallocation of index in chunks of MEMCHUNK
						if (i % MEMCHUNK == 0)
							index =	(long *) realloc(index, 4 * (i + MEMCHUNK) * sizeof(long));
							index[0 + i * 4] = u1;
							index[1 + i * 4] = u2;
							index[2 + i * 4] = v1;
							index[3 + i * 4] = v2;
							i++;
						}
				}
				// set positive limit for v2
				if (u1 == v1)
					limit = abs(u2);
				else
					limit = 0;
				for (v2 = limit; v2 <= l.max_rad; v2++) {
					v2s = v2 * v2;
					uv2 = u2 + v2;
					uv2s = uv2 * uv2;
					if (((uv1s + uv2s) <= mrs) &&
						((v1s + v2s) <= mrs) &&
						(((u1s + u2s) <= l2s) ||
						 ((v1s + v2s) <= l2s) || ((uv1s + uv2s) <= l2s))) {
					// reallocation of index in chunks of MEMCHUNK
					if (i % MEMCHUNK == 0)
						index = (long *) realloc(index, 4 * (i + MEMCHUNK) * sizeof(long));
						index[0 + i * 4] = u1;
						index[1 + i * 4] = u2;
						index[2 + i * 4] = v1;
						index[3 + i * 4] = v2;
						i++;
					}
				}
			}
		}
    }
    // optimum memory usage
    *bs_cnt = i;
    index = (long *) realloc(index, 4 * (*bs_cnt) * sizeof(long));

    free(d);

    return (index);
}
