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


/* =====================================================================
**
** Reassemble the Image
** This is basically putting together the phases and amplitudes.
** Then a retransformation FFT is performed.
**
** =====================================================================
*/

void assemble(int nx, int ny, float *amp, float *phs, float *rec, maininfo info)
{
    float *mask = ellmask(nx, ny, NULL, info.tc.max_rad, info.tc.max_rad);

    static int INIT = 0;
    static long N;
    static fftw_complex *in;
    static fftw_complex *out;
    fftw_plan fft;

    long i, j;
    long r1, r2, r3;

    if (INIT == 0) {
		N = nx * ny;
		in = fftw_malloc(N * sizeof(fftw_complex));
		out = fftw_malloc(N * sizeof(fftw_complex));
		INIT = 1;
    } else if (N / nx != ny) {
		fftw_cleanup();
		fftw_free(in);
		fftw_free(out);
		N = nx * ny;
		in = fftw_malloc(N * sizeof(fftw_complex));
		out = fftw_malloc(N * sizeof(fftw_complex));
		INIT = 1;
    }
    fft = fftw_plan_dft_2d(nx, ny, in, out, +1, FFTW_ESTIMATE);

    for (j = 0; j < ny; j++) {
		for (i = 0; i < nx; i++) {
			idx(i, j, nx, r1);
			shift_idx(i, j, nx / 2, ny / 2, nx, ny, r2);
			phs_idx(i, j, nx, r3);
			in[r2][0] = amp[r1] * phs[r3] * mask[r1];
			in[r2][1] = amp[r1] * phs[r3 + 1] * mask[r1];
		}
    }
    fftw_execute(fft);
    for (j = 0; j < ny; j++) {
		for (i = 0; i < nx; i++) {
			idx(i, j, nx, r1);
			rec[r1] = (float) out[r1][0];
		}
    }

    fftw_cleanup();
    free(mask);
}
