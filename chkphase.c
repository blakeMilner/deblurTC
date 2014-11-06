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


float chkphase(slaveinfo l, float *phs)
{
    static long INIT = 0;
    static float *old = NULL, *new = NULL;

    float *tmp = NULL;
    long i;
    float t, var = 0.0;

    if (INIT / l.sfsizex != l.sfsizey) {
		if (old != NULL)
			free(old);
		if (new != NULL)
			free(new);
		INIT = l.sfsizex * l.sfsizex;
		old = (float *) calloc(INIT, sizeof(float));
		new = (float *) calloc(INIT, sizeof(float));
    }

    for (i = 0; i < INIT; i++) {
		new[i] = atan2(phs[2 * i + 1], phs[2 * i]);
		t = (old[i] - new[i]);
		var += t * t;
    }

    // some fun stuff: pointer bending!
    tmp = old;
    old = new;
    new = tmp;
    tmp = NULL;

    return (var);
}
