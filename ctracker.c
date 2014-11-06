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


void ctracker(float *ref, float *in, int nx, int ny, int nfr, float *out)
{
/*
 * Declaration of local variables
 */
    static int nxh = 0, nyh = 0;	/* coordinates of origin */
    static float *win = NULL;	/* apodisation window */
    static float *fit = NULL;	/* fitted surface */
    static double *temp1 = NULL;	/* temporary data storage space */
    static double *temp2 = NULL;	/* temporary data storage space */
    static double *ccr = NULL;	/* cross correlation */
    static fftw_complex *fref = NULL;	/* FFT of reference img */
    static fftw_complex *fimage = NULL;	/* FFT of image */
    static fftw_complex *fccr = NULL;	/* FFT of cross cor */

    long k, l;			/* helper variables */
    long i;			/* index counter variable */
    long N;			/* total number of pixels in subfield */

    double maxval;		/* value of max of cross correlation */
    long maxpos = 0;		/* vector position of maximum */
    int mx, my;			/* x and y position of max */
    double mval;		/* mean intensity value */
    fftw_plan fftref;		/* planning to FFT mean */
    fftw_plan fftimage;		/* planning to FFT image */
    fftw_plan crosscor;		/* planning to FFT cross cor */

/*
 *	Initialize local variables
 */
    N = nx * ny;

    if ((nxh == 0) || (nyh == 0)) {
//      Allocate memory
		nxh = nx / 2;
		nyh = ny / 2;
		win = (float *) malloc(N * sizeof(float));
		fit = (float *) malloc(N * sizeof(float));
		temp1 = (double *) malloc(N * sizeof(double));
		temp2 = (double *) malloc(N * sizeof(double));
		ccr = (double *) malloc(N * sizeof(double));
		fref = fftw_malloc(nx * (nyh + 1) * sizeof(fftw_complex));
		fimage = fftw_malloc(nx * (nyh + 1) * sizeof(fftw_complex));
		fccr = fftw_malloc(nx * (nyh + 1) * sizeof(fftw_complex));
		hanming(win, nx, ny, 0.5);
    }
    else if ((nxh != nx / 2) || (nyh != ny / 2)) {
//      Clear all used memory
		free(win);
		free(fit);
		free(temp1);
		free(temp2);
		free(ccr);
		fftw_free(fref);
		fftw_free(fimage);
		fftw_free(fccr);
	//      Allocate memory
		nxh = nx / 2;
		nyh = ny / 2;
		win = (float *) malloc(N * sizeof(float));
		hanming(win, nx, ny, 0.5);
		fit = (float *) malloc(N * sizeof(float));
		temp1 = (double *) malloc(N * sizeof(double));
		temp2 = (double *) malloc(N * sizeof(double));
		ccr = (double *) malloc(N * sizeof(double));
		fref = fftw_malloc(nx * (nyh + 1) * sizeof(fftw_complex));
		fimage = fftw_malloc(nx * (nyh + 1) * sizeof(fftw_complex));
		fccr = fftw_malloc(nx * (nyh + 1) * sizeof(fftw_complex));
    }

    memset(temp1, 0.0, N * sizeof(double));
    memset(temp2, 0.0, N * sizeof(double));
    memset(ccr, 0.0, N * sizeof(double));
    memset(fref, 0.0, nx * (nyh + 1) * sizeof(fftw_complex));
    memset(fimage, 0.0, nx * (nyh + 1) * sizeof(fftw_complex));
    memset(fccr, 0.0, nx * (nyh + 1) * sizeof(fftw_complex));
    fftref = fftw_plan_dft_r2c_2d(nx, ny, temp1, fref, FFTW_ESTIMATE);
    fftimage = fftw_plan_dft_r2c_2d(nx, ny, temp2, fimage, FFTW_ESTIMATE);
    crosscor = fftw_plan_dft_c2r_2d(nx, ny, fccr, ccr, FFTW_ESTIMATE);

/*
 *	Set reference frame if present. If not present, then set first image
 *	as reference. Reference is stored in temp1. Then calculate FFT(Reference, -1)
 */
    if (ref != NULL) {
		surfit(ref, nx, ny, 2, fit);
		stats(fit, N, 0.0, &mval, NULL);
		for (i = 0; i < N; i++) {
			temp1[i] = (double) (win[i] * (ref[i] - fit[i])) + mval;
		}
    }
    else {
		surfit(in, nx, ny, 2, fit);
		stats(fit, N, 0.0, &mval, NULL);
		for (i = 0; i < N; i++) {
			temp1[i] = (double) (win[i] * (in[i] - fit[i])) + mval;
		}
    }
    fftw_execute(fftref);

/*
 *	Loop through burst
 */
	for (l = 0; l < nfr; l++) {
		/* index variable for image frame */
		k = l * N;
	/*
	 *	Initialize Images. Surface fit image. Store resultant image in dtemp2.
	 *
	 *	Surface fit
	 */
		surfit(&in[k], nx, ny, 2, fit);
		stats(fit, N, 0.0, &mval, NULL);
	/*	Counter variable for image frame */
		for (i = 0; i < N; i++) {
			temp2[i] = (double) (win[i] * (in[i + k] - fit[i])) + mval;
		}
	/*
	 *	Do FFT(Image,-1)
	 */
		fftw_execute(fftimage);
	/*
	 *	Calculation of FFT(Image,-1)*Conj(FFT(Reference,-1))
	 */
		for (i = 0; i < nx * (nyh + 1); i++) {
			fccr[i][0] = fimage[i][0] * fref[i][0]
			+ fimage[i][1] * fref[i][1];
			fccr[i][1] = fimage[i][1] * fref[i][0]
			- fimage[i][0] * fref[i][1];
		}
	/*
	 *        Calculation of FFT(FFT(Image,-1)*Conj(FFT(Reference,-1)),+1)
	 */
		fftw_execute(crosscor);
	/*
	 *        Find maximum
	 */
		maxval = ccr[0];
		for (i = 0; i < N; i++) {
			if (ccr[i] > maxval) {
			maxpos = i;
			maxval = ccr[i];
			}
		}
	/*
	 *        Convert information about max into relevant values
	 *        Assign displacements to output variable
	 */

		mx = (int) (maxpos % nx);
		if (mx > nxh){
			mx -= nx;
		}

		/* setting max displacement to zero */
		out[2 * l] = 0;
		//out[2 * l] = mx;
		my = (int) (maxpos / nx);
		if (my > nyh){
			my -= ny;
		}

		/* setting max displacement to zero */
		//out[2 * l + 1] = my;
		out[2 * l + 1] = 0;
	}

    fftw_cleanup();
}
