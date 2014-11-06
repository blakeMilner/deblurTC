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


void surfit(float *in, int nx, int ny, int deg, float *out)
{
    // static variables
    static long n = 0;
    static int nrcoeffs = 0;
    static gsl_vector *a, *b, *c;

    // helper
    long i;
    // calculation matrices
    gsl_matrix *ut;
    gsl_matrix *kk;

    if (n == 0) {
	// set variables
	n = nx * ny;
	nrcoeffs = (deg + 1) * (deg + 2) / 2;
	// coefficient vector
	a = gsl_vector_alloc(n);
	b = gsl_vector_alloc(n);
	c = gsl_vector_alloc(nrcoeffs);
    }

    if ((n != 0)
	&& ((n / nx != ny) || (nrcoeffs != (deg + 1) * (deg + 2) / 2))) {
	// clear memory
	gsl_vector_free(a);
	gsl_vector_free(b);
	gsl_vector_free(c);
	// set variables
	n = nx * ny;
	nrcoeffs = (deg + 1) * (deg + 2) / 2;
	// coefficient vector
	a = gsl_vector_alloc(n);
	b = gsl_vector_alloc(n);
	c = gsl_vector_alloc(nrcoeffs);
    }
    // assign input
    for (i = 0; i < n; i++)
	gsl_vector_set(a, i, (double) in[i]);

    // Initialise matrix
    init_matrix(nx, ny, deg, &ut, &kk);

    // Get the fitted coefficients
    gsl_blas_dgemv(CblasNoTrans, 1.0, kk, a, 0.0, c);
    // Create fit of surface
    gsl_blas_dgemv(CblasTrans, 1.0, ut, c, 0.0, b);

    // assign output
    for (i = 0; i < n; i++)
	out[i] = (float) gsl_vector_get(b, i);
}
