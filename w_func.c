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


vect *w_func(float *in, char *expr, float comp, int size)
{
    /* resultant vector */
    vect *out;

    /* helper variables */
    int i, count = 0;

    /*
     *  Allocate memory for resultant vector
     */
    out = (vect *) malloc(sizeof(vect));

    /*
     *  Allocate memory for array vector which has indices
     *  of numbers for which the condition is fulfilled
     */
    out->res_vec = (float *) malloc(size * sizeof(float));

    /* Parse through possible conditions */
    if (!strcmp(expr, "<")) {
		for (i = 0; i < size; i++) {
			if (in[i] < comp)
			out->res_vec[count++] = i;
		}
    }

    if (!strcmp(expr, "<=") || !strcmp(expr, "=<")) {
		for (i = 0; i < size; i++) {
			if (in[i] <= comp)
			out->res_vec[count++] = i;
		}
    }

    if (!strcmp(expr, "==")) {
		for (i = 0; i < size; i++) {
			if (in[i] == comp)
			out->res_vec[count++] = i;
		}
	}

    if (!strcmp(expr, ">=") || !strcmp(expr, "=>")) {
		for (i = 0; i < size; i++) {
			if (in[i] >= comp)
			out->res_vec[count++] = i;
		}
    }

    if (!strcmp(expr, ">")) {
		for (i = 0; i < size; i++) {
			if (in[i] > comp)
			out->res_vec[count++] = i;
		}
    }

    if (!strcmp(expr, "!=")) {
		for (i = 0; i < size; i++) {
			if (in[i] != comp)
			out->res_vec[count++] = i;
		}
    }

    if (!strcmp(expr, "|<|")) {
		for (i = 0; i < size; i++) {
			if (fabs(in[i]) < comp)
			out->res_vec[count++] = i;
		}
    }

    if (!strcmp(expr, "|<=|") || !strcmp(expr, "|=<|")) {
		for (i = 0; i < size; i++) {
			if (fabs(in[i]) <= comp)
			out->res_vec[count++] = i;
		}
    }

    if (!strcmp(expr, "|==|")) {
		for (i = 0; i < size; i++) {
			if (fabs(in[i]) == comp)
			out->res_vec[count++] = i;
		}
    }

    if (!strcmp(expr, "|>=|") || !strcmp(expr, "|=>|")) {
		for (i = 0; i < size; i++) {
			if (fabs(in[i]) >= comp)
			out->res_vec[count++] = i;
		}
    }

    if (!strcmp(expr, "|>|")) {
		for (i = 0; i < size; i++) {
			if (fabs(in[i]) > comp)
			out->res_vec[count++] = i;
		}
    }

    if (!strcmp(expr, "|!=|")) {
		for (i = 0; i < size; i++) {
			if (fabs(in[i]) != comp)
			out->res_vec[count++] = i;
		}
    }

    if (count == 0) {
		out->size = -1;
		free(out->res_vec);
		out->res_vec = NULL;
    }
    else {
		out->size = count;
		out->res_vec =
	    (float *) realloc(out->res_vec, count * sizeof(float));
    }

    return (out);
}
