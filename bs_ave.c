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
 *  bs_ave: compute the average bispectrum of an image series
 *
 *    images: pointer to series of images
 *    win: apodization window for the fourier transform of the images
 *    nx, ny: size of images
 *    nfr: # of frames
 *    index: pointer to matrix with bispectrum elements to be used
 *    bs: BiSpectrum output
 *    w: Weight output
 */

float *bs_ave(float *images, float *win, slaveinfo l, long *index,
	      long bs_cnt, float *bsc, float *wc, float *amp, int maxk,
	      int *shifts)
{
    long i, j, k, m;		// loop & helper variables
    long i1, i2, i3;		// index variables
    double t1, t2, t3, t4;	// temporary real variables
    float ct1[2], ct2[2];	// temporary complex variables
    long nx, ny, nxh, nyh;	// half of subfield size
    long nfr, N, cnt;		// # of pixels in one image, loop variable
    double s = 0;		// mean intensity in image
    long u1, u2, v1, v2;	// bispectrum vectors
    fftw_complex *im1;		// temporary matrices
    fftw_complex *im2;
    float *bsr, *bsi;		// used for calculation of SNR
    float *origin;		// phase values around origin (needed for init)

    fftw_plan fftimage;		// plan for fourier transform

    nfr = l.nrofframes;
    nx = l.sfsizex;
    ny = l.sfsizey;
    nxh = nx / 2;
    nyh = ny / 2;
    N = nx * ny;       // subfield size in pixels

    origin = (float *) calloc(2 * maxk, sizeof(float));

    im1 = fftw_malloc(N * sizeof(fftw_complex));
    im2 = fftw_malloc(N * sizeof(fftw_complex));

    // setting up resultant matrices
    bsr = (float *) calloc(bs_cnt, sizeof(float));
    bsi = (float *) calloc(bs_cnt, sizeof(float));

    fftimage = fftw_plan_dft_2d(nx, ny, im1, im2, -1, FFTW_ESTIMATE);

	for (k = 0; k < nfr; k++) {
		// set a helper index
		m = k * N;
		// Calculate the mean of the image
		stats(&images[m], N, 0, &s, NULL);
		// store the windowed image in image1
		for (j = 0; j < ny; j++) {
			for (i = 0; i < nx; i++) {
				idx(i, j, nx, i1);
				im1[i1][0] = (double) ((images[i1 + m] - s) * win[i1] + s);  // problem here, with image
				im1[i1][1] = 0.0;

			}
		}
		// a) do FFT(imagei,-1)
		//    write the result into im2 (see fftw_plan fftimage)
		// b) shift the image, so low frequencies are in the origin
		//    write the results back into im1, scaled with
		//    1/(nx*ny) (not done by FFTW but by IDL does, thanks
		//    Alexandra Tritschler for the hint!)
		fftw_execute(fftimage);
		for (j = 0; j < ny; j++) {		// iterates through subfield pixels
			for (i = 0; i < nx; i++) {
				idx(i, j, nx, i1);		// jump to pixel position
				shift_idx(i, j, nxh, nyh, nx, ny, i2);
				im1[i1][0] = im2[i2][0] / N;		// complex part
				im1[i1][1] = im2[i2][1] / N;		// imaginary part
				amp[i1] +=
					(float) (sqrt
						 (im1[i1][0] * im1[i1][0] +			// why are these so large??
						  im1[i1][1] * im1[i1][1]) / nfr);
			}
		}
		// phase around origin
		for (cnt = 0; cnt < maxk; cnt++) {
			idx(nxh + shifts[2 * cnt], nyh + shifts[2 * cnt + 1], nx, i1);
			origin[2 * cnt] += im1[i1][0];
			origin[2 * cnt + 1] += im1[i1][1];
		}
		// calculate the mean raw bispectrum fast in non redundant matrix
		for (cnt = 0; cnt < bs_cnt; cnt++) {
			// get vectors
			u1 = index[0 + cnt * 4];
			u2 = index[1 + cnt * 4];
			v1 = index[2 + cnt * 4];
			v2 = index[3 + cnt * 4];
			// get matrix indices of vectors
			idx(nxh - u1, nyh + u2, nx, i1);
			idx(nxh - v1, nyh + v2, nx, i2);
			idx(nxh - u1 - v1, nyh + u2 + v2, nx, i3);
			// calculate bispectrum i(u)*i(v)*conj(i(u+v))
			//   (a+ib)(c+id)=((ac-bd)+i(bc+ad))
			ct1[0] = im1[i1][0] * im1[i2][0] - im1[i1][1] * im1[i2][1];
			ct1[1] = im1[i1][1] * im1[i2][0] + im1[i1][0] * im1[i2][1];
			//   (a+ib)(c-id)=((ac+bd)+i(bc-ad))
			ct2[0] = ct1[0] * im1[i3][0] + ct1[1] * im1[i3][1];
			ct2[1] = ct1[1] * im1[i3][0] - ct1[0] * im1[i3][1];
			// assign values
			bsc[2 * cnt] += ct2[0];
			//      bsc[2*cnt] += ct2[0]   // noise terms attached
			//         - ((im1[i1][0]*im1[i1][0] + im1[i1][1]*im1[i1][1])
			//           +(im1[i2][0]*im1[i2][0] + im1[i2][1]*im1[i2][1])
			//           +(im1[i3][0]*im1[i3][0] + im1[i3][1]*im1[i3][1]));
			bsc[2 * cnt + 1] += ct2[1];
			bsr[cnt] += bsc[2 * cnt] * bsc[2 * cnt];
			bsi[cnt] += bsc[2 * cnt + 1] * bsc[2 * cnt + 1];
		}
	}

    // phase around origin is averaged and normalized
    for (cnt = 0; cnt < maxk; cnt++) {
		t1 = cmod(&origin[2 * cnt]);
		origin[2 * cnt] /= t1;
		origin[2 * cnt + 1] /= t1;
    }

    // create bispectrum SNR (weights) and mean in non-redundant matrices
	for (cnt = 0; cnt < bs_cnt; cnt++) {
		// calculate SNR
		// sum of 2 normally - distributed random variables (real & imaginary part)
		t1 = bsc[2 * cnt] / nfr;
		t1 *= t1;		// (squared) mean real part
		t2 = bsc[2 * cnt + 1] / nfr;
		t2 *= t2;		// (squared) mean imaginary part
		t3 = fabs(bsr[cnt] - nfr * t1) / (nfr - 1);	// variance of real part
		t4 = fabs(bsi[cnt] - nfr * t2) / (nfr - 1);	// variance of imaginary part
		wc[cnt] = (float) sqrt(t1 + t2) * sqrt(nfr) / sqrt(t3 + t4);	// => SNR of MEAN absolute

			// make mean and normalize the phases if nonzero
		if ((bsc[2 * cnt] != 0.0) || (bsc[2 * cnt + 1] != 0.0)) {
			t1 = cmod(&bsc[2 * cnt]);
			bsc[2 * cnt] /= t1;
			bsc[2 * cnt + 1] /= t1;
		}
    }

    fftw_free(im1);
    fftw_free(im2);
    free(bsr);
    free(bsi);

    fftw_cleanup();

    return (origin);
}
