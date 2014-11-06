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


void init_matrix(int nx, int ny, int deg, gsl_matrix ** ut, gsl_matrix ** kk)
{
    /* important static variables */
    static int init_done = 0;
    static gsl_matrix *out_ut;
    static gsl_matrix *out_kk;
    static int n1, n2, n3;

    /* helper variables */
    unsigned long nrcoeffs = (deg + 1) * (deg + 2) / 2;
    unsigned long n = nx * ny;
    long i, k, l;
    unsigned long j;
    int sig;
    gsl_vector *vx;
    gsl_vector *vy;
    gsl_matrix *tmp1;
    gsl_matrix *tmp2;
    gsl_permutation *p;

    /* Initialise fit matrix if needed */
    if ((init_done == 1) &&
	((n1 != nx) || (n2 != ny) || (n3 != deg) ||
	 ((*out_ut).size1 != nrcoeffs) || ((*out_ut).size2 != n) ||
	 ((*out_kk).size1 != nrcoeffs) || ((*out_kk).size2 != n))) {
		/* free memory */
		gsl_matrix_free(out_ut);
		gsl_matrix_free(out_kk);
		/* set init flag */
		init_done = 0;
    }
    if (init_done == 0) {
		// set control values
		n1 = nx;
		n2 = ny;
		n3 = deg;
		/* init helper variables */
		vx = gsl_vector_calloc(nx);
		vy = gsl_vector_calloc(ny);
		tmp1 = gsl_matrix_calloc(nrcoeffs, nrcoeffs);
		tmp2 = gsl_matrix_calloc(nrcoeffs, nrcoeffs);
		p = gsl_permutation_calloc(nrcoeffs);
		/* allocate memory for fit matrix */
		out_ut = gsl_matrix_calloc(nrcoeffs, n);
		out_kk = gsl_matrix_calloc(nrcoeffs, n);
		/* Initialise the base matrix which has the grid inside: *
		 * it has as many rows as paramters to be fitted and the *
		 * grid for each of the parameters                       */
		for (i = 0; i < nx; i++) {
			gsl_vector_set(vx, i, i);
		}
		for (i = 0; i < ny; i++) {
			gsl_vector_set(vy, i, i);
		}
		i = 0;
		for (k = 0; k <= deg; k++) {
			for (l = 0; l <= deg; l++) {
				if ((k + l) <= deg) {
					for (j = 0; j < n; j++) {
					gsl_matrix_set(out_ut, i, j,
							   pow(gsl_vector_get
							   (vx, (int) j % nx),
							   k) * pow(gsl_vector_get(vy,
										   (int) j
										   / nx),
									l));
					}
					i++;
				}
			}
		}
		/* Invert the base matrix squared and multiply. This     *
		 * will result in 1/a^2 * a = 1/a and gives us the       *
		 * base for the least squares fit                        */
		gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, out_ut, out_ut, 0.0,
				   tmp1);
		gsl_linalg_LU_decomp(tmp1, p, &sig);
		gsl_linalg_LU_invert(tmp1, p, tmp2);
		gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, tmp2, out_ut, 0.0,
				   out_kk);
		/* Free some memory */
		gsl_vector_free(vx);
		gsl_vector_free(vy);
		gsl_matrix_free(tmp1);
		gsl_matrix_free(tmp2);
		gsl_permutation_free(p);
		/* set init flag */
		init_done = 1;
    }
    /* Set the pointers right */
    *ut = out_ut;
    *kk = out_kk;
}
