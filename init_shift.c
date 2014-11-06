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
*************EXTENDED KNOX-THOMPSON PART FOR PHASE RECONSTRUCTION************
****************************************************************************/
/*
** Compute shifts that will be utilized to calculate the Knox-Thompson
** Cross Spectra.
** Input:
**  nxh: half of field size in x direction
**  nyh: half of field size in y direction
**  len: maximum length of shift vector
**   NOTE: for KT this may not be higher than the
**    seeing limit!
**  maxk: number of shifts computed
**  shifts: pointer to resultant array
**
** written by Friedrich Woeger
**
** last change: 14. 10. 03
*/

void init_shift(int nxh, int nyh, int len, int *maxk, int **shifts)
{
    // helper variables
    int limit;
    int count;
    int rad;
    int i, j;
    int x, y;

    if (len == 1) {
	// set count value and allocate memory for original Knox-Thompson shifts
		count = 6;
		(*shifts) = (int *) realloc((*shifts), count * sizeof(int));
		(*shifts)[0] = 0;
		(*shifts)[1] = 0;
		(*shifts)[2] = 1;
		(*shifts)[3] = 0;
		(*shifts)[4] = 0;
		(*shifts)[5] = 1;
	} else {
		// set count value and allocate memory for zero shifts
		count = 2;
		(*shifts) = (int *) realloc((*shifts), count * sizeof(int));
		(*shifts)[0] = 0;
		(*shifts)[1] = 0;

		// limit = minimum of half of window size, since we start in the middle
		limit = (nxh <= nyh) ? nxh : nyh;

		// Begin computation
		for (rad = 1; rad <= len; rad++) {

			if (rad >= limit)
			break;

			for (i = 0; i <= rad; i++) {
				for (j = 0; j < 4 * i; j++) {
					if (j <= i) {
						x = i;
						y = j;
					}
					if ((j > i) && (j <= 3 * i)) {
						x = 2 * i - j;
						y = i;
					}
					if ((j > 3 * i) && (j <= 4 * i)) {
						x = -i;
						y = 4 * i - j;
					}

					if (rad < len) {
						if (((x * x + y * y) <= (rad * rad))
							&& ((x * x + y * y) > (rad - 1) * (rad - 1))) {
							(*shifts) =
							(int *) realloc((*shifts),
									(count + 2) * sizeof(int));
							(*shifts)[count] = x;
							(*shifts)[count + 1] = y;
							count += 2;
						}
					}
					else {
						if (((x * x + y * y) < (rad * rad))
							&& ((x * x + y * y) > (rad - 1) * (rad - 1))) {
							(*shifts) =
							(int *) realloc((*shifts),
									(count + 2) * sizeof(int));
							(*shifts)[count] = x;
							(*shifts)[count + 1] = y;
							count += 2;
						}
					}
				}
			}
		}
	}

    *maxk = (int) count / 2;
}
