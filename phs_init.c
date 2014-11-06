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
 *  phs_init: initialise the phase matrices
 *
 *    p1, p2: phases to be initialised
 *    nx, ny: size of images
 */

void phs_init(float *p1, float *p2, float *pc, slaveinfo l,
	      int maxk, int *shifts, float *origin)
{
    static int INIT = 0;
    static float *mask = NULL;
    static long N;
    long nx, ny;
    long i;
    long i1;

    nx = l.sfsizex;
    ny = l.sfsizey;

    // compute elliptic mask
    if ((INIT != l.max_rad) || ny != (N / nx)) {
		N = nx * ny;
		if (mask != NULL)
			free(mask);
		mask = ellmask(nx, ny, NULL, l.max_rad, l.max_rad);
		INIT = l.max_rad;
    }
    // set phases
    for (i = 0; i < maxk; i++) {
		phs_idx(nx / 2 + shifts[2 * i], ny / 2 + shifts[2 * i + 1], nx, i1);
		p1[i1] = origin[2 * i];
		p1[i1 + 1] = origin[2 * i + 1];
		p2[i1] = origin[2 * i];
		p2[i1 + 1] = origin[2 * i + 1];
		pc[i1 / 2] = 1.0;
		phs_idx(nx / 2 - shifts[2 * i], ny / 2 - shifts[2 * i + 1], nx,	i1);
		p1[i1] = origin[2 * i];
		p1[i1 + 1] = -origin[2 * i + 1];
		p2[i1] = origin[2 * i];
		p2[i1 + 1] = -origin[2 * i + 1];
		pc[i1 / 2] = 1.0;
    }

    // reset phase values to decent values
    for (i = 0; i < N; i++) {
		p1[2 * i] *= mask[i];
		p1[2 * i + 1] *= mask[i];
    }
}
