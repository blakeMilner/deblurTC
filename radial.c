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



float *rad2im(vect * in, int nx, int ny, int *ori, float aspect)
{
    /* static helper variables */
    static int INIT = 0;
    static float *matrix = NULL;
    /* helper variables */
    vect *w_vect;
    int szh[2];
    long int i, j;

    /* output variable */
    float *out;

    /* memory allocation of static helper matrix */
    if (INIT / nx != ny) {
		if (matrix != NULL) {
			free(matrix);
			matrix = NULL;
		}
		matrix = (float *) malloc(nx * ny * sizeof(float));
		INIT = nx * ny;
    }
    /* memory allocation for output variable */
    out = (float *) calloc(nx * ny, sizeof(float));

    if (ori == NULL) {
		szh[0] = nx / 2;
		szh[1] = ny / 2;
    } else {
		szh[0] = ori[0];
		szh[1] = ori[1];
    }

    /* initialise helper matrix */
    for (j = -szh[1]; j < szh[1]; j++) {
		for (i = -szh[0]; i < szh[0]; i++) {
			matrix[(j + szh[1]) * nx + (i + szh[0])]
			=
			(float) ((int)
				 (sqrt(i * i + (j * aspect) * (j * aspect)) +
				  0.5));
		}
    }

    for (i = 0; i < (in->size); i++) {
		w_vect = w_func(matrix, "==", (float) i, nx * ny);
		for (j = 0; j < (w_vect->size); j++) {
			out[(int) (w_vect->res_vec[j])] += in->res_vec[i];
		}
		if (w_vect->size != -1) {
			free_vect(&w_vect);
		}
    }

    return (out);
}
