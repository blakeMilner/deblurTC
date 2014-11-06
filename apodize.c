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

/* creates a either a hamming or hanning window based on the value of alpha fed in */
void hanming(float *res, int nx, int ny, float alpha)
{
    /* For a Hamming use alpha = 0.53836 */
    /* For a Hanning use alpha = 0.5     */
    /* Initialise helpers */
    int i, j;
    float fact;
    float *w1, *w2;

    /* Initialise constants */
    double pi = 3.141592654;
    double two_pi = 2.0 * pi;

    nx = (nx > 0) ? nx : 1;
    ny = (ny > 0) ? ny : 1;

    /* Initialize window arrays */
    w1 = (float *) (malloc(nx * sizeof(float)));
    for (i = 0; i < nx; i++) {
		w1[i] = alpha +
			(alpha - 1.) *
			((float)
			 (cos((double) (two_pi * ((float) i) / ((float) nx)))));
    }
    if (ny > 1) {
		w2 = (float *) (malloc(ny * sizeof(float)));
		for (i = 0; i < ny; i++) {
			w2[i] = alpha +
			(alpha - 1.) *
			((float)
			 (cos((double) (two_pi * ((float) i) / ((float) ny)))));
		}
    }

    /* Set window */
    for (j = 0; j < ny; j++) {
		fact = (ny > 1) ? w2[j] : 1.0;
		for (i = 0; i < nx; i++) {
			res[j * nx + i] = fact * w1[i];
		}
    }

    free(w1);
    if (ny > 1)
    	free(w2);
}


void frachamming(float *res, int nx, int ny, float frac, int maxk, int *sh)
{
    int i, j;
    int nxf, nyf, nxfh, nyfh;
    long N, l;
    float *wx, *wy;
    float *Hx, *Hy;
    float m = 0, maximum;

    fftw_complex *in, *out;
    fftw_plan fft, fftb;

    N = nx * ny;
    nxf = 2 * (int) ((frac * nx + 1) / 2);
    nyf = 2 * (int) ((frac * ny + 1) / 2);

    if (nxf == 0) {
		for (i = 0; i < nx * ny; i++)
			res[i] = 1.0;
	}
    else {
		wx = (float *) malloc(nx * sizeof(float));
		wy = (float *) malloc(ny * sizeof(float));
		Hx = (float *) malloc(nxf * sizeof(float));
		Hy = (float *) malloc(nyf * sizeof(float));
		nxfh = nxf / 2;
		nyfh = nyf / 2;
		/* get hamming values */
		hanming(Hx, nxf, 1, 0.5);
		hanming(Hy, nyf, 1, 0.5);
		/* init the vectors to 1.0 (for values in the middle!) */
		for (i = 0; i < nx; i++) {
			wx[i] = 1.0;
		}
		for (i = 0; i < ny; i++) {
			wy[i] = 1.0;
		}
		/* set hamming values */
		for (i = 0; i < nxfh; i++) {
			wx[i] = Hx[i];
			wx[(nx - nxfh) + i] = Hx[nxfh + i];
		}
		for (i = 0; i < nyfh; i++) {
			wy[i] = Hy[i];
			wy[(ny - nyfh) + i] = Hy[nyfh + i];
		}
		/* compute the tensor product */
		for (j = 0; j < ny; j++) {
			for (i = 0; i < nx; i++) {
			res[j * nx + i] = wy[j] * wx[i];
			}
		}
		free(wx);
		free(wy);
		free(Hx);
		free(Hy);
	}

    /* this isn't ever used in KISIP since maxk is always fed in as 0 */
	if ((maxk != 0) && (sh != NULL)) {
		// allocate memory for fourier transform
		in = fftw_malloc(N * sizeof(fftw_complex));
		out = fftw_malloc(N * sizeof(fftw_complex));
		fft = fftw_plan_dft_2d(nx, ny, in, out, -1, FFTW_ESTIMATE);
		fftb = fftw_plan_dft_2d(nx, ny, out, in, +1, FFTW_ESTIMATE);

		// set input matrix
		for (l = 0; l < N; l++) {
			in[l][0] = (double) res[l];
			in[l][1] = 0.0;
		}
		// execute FFT
		fftw_execute(fft);
		// Now set phases at shifts to zero. This follows:
		// Christoph Keller: Opt. Apod. for Speckle Imaging of ext. sourc.
		// ASP Conference Series #183 (1998)
		// Shifts are relative to centered spectrum!
		for (i = 1; i < maxk; i++) {
			shift_idx(nx / 2 + sh[2 * i], ny / 2 + sh[2 * i + 1], nx / 2,
				  ny / 2, nx, ny, l);
			out[l][0] = 0.0;
			out[l][1] = 0.0;
			shift_idx(nx / 2 - sh[2 * i], ny / 2 - sh[2 * i + 1], nx / 2,
				  ny / 2, nx, ny, l);
			out[l][0] = 0.0;
			out[l][1] = 0.0;
		}
		// execute reFFT
		fftw_execute(fftb);
		// set result matrix
		for (l = 0; l < N; l++) {
			res[l] = (float) in[l][0];
		}
		// free memory
		fftw_cleanup();
		//fftw_destroy_plan(fft);
		//fftw_destroy_plan(fftb);
		fftw_free(in);
		fftw_free(out);
	}

    // find max and normalize
    maximum = res[0];
    for (l = 0; l < N; l++) {
		if (res[l] > maximum)
			maximum = res[l];
		}
    for (l = 0; l < N; l++) {
		res[l] /= maximum;
		m += res[l] * res[l];
    }
    m /= N;
    m = sqrt(m);
    for (l = 0; l < N; l++) {
    	res[l] /= m;
    }
}


float *ellmask(int nx, int ny, int *ori, int len_x, int len_y)
{
    int i;
    float aspect;
    vect *v;
    float *result;

    v = (vect *) malloc(sizeof(vect));
    v->size = len_x;
    v->res_vec = (float *) malloc(len_x * sizeof(float));

    for (i = 0; i < len_x; i++) {
    	(v)->res_vec[i] = 1.0;
    }

    aspect = (float) len_y / len_x;
    result = rad2im(v, nx, ny, ori, aspect);

    free_vect(&v);
    return (result);
}


void assmask(float *res, int nx, int ny, float frac)
{
    // Initialise helpers
    int i, j;
    long offx, offy;
    float *temp;

    // set offsets for window
    offx = 0;
    offy = 0;
    //offx = (int) rint(frac / 2 * nx);
    //offy = (int) rint(frac / 2 * ny);

    temp = (float *) malloc((nx - 2 * offx) * (ny - 2 * offy) * sizeof(float));
    hanming(temp, (nx - 2 * offx), (ny - 2 * offy), 0.5);

    memset(res, 0, nx * ny * sizeof(float));

    for (j = offy; j < ny - offy; j++) {
		for (i = offx; i < nx - offx; i++) {
			res[j * nx + i] =
			temp[(j - offy) * (nx - 2 * offx) + (i - offx)];
		}
    }
    free(temp);
}
