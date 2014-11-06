/* $Id$
 *
 *  (C) 2003-2007, Friedrich Woeger <woeger@kis.uni-freiburg.de>,
 *  Kiepenheuer-Institut fuer Sonnenphysik, Freiburg (Breisgau), Germany
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


/*	=====================================================================
**
**	supplementary speckle math routines
**
**	=====================================================================
**
**	uses:		fftw3 library:
**			http://www.fftw.org
**
**			Library for fast fourier transform
**
**	=====================================================================
**
**	Author:		F. Woeger
**			Kiepenheuer-Institut fuer Sonnenphysik
**			Freiburg, Germany
**
**	Written 17. 06. 2003
**
**	=====================================================================
*/

#ifndef KISIP_SPECKLE_MATH_H
#define KISIP_SPECKLE_MATH_H

#include <fftw3.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multifit_nlin.h>
#include "inlines.h"

typedef struct vect {
    long size;
    float *res_vec;
} vect;

void free_vect(vect ** v);
void ctracker(float *ref, float *img, int nsx, int nsy, int nfr, float *out);
void surfit(float *in, int nx, int ny, int deg, float *out);
void init_matrix(int nx, int ny, int deg, gsl_matrix ** ut, gsl_matrix ** kk);
vect *w_func(float *in, char *expr, float comp, int size);
void mean(float *burst, int nx, int ny, int nfr, float *out);
void stats(float *in, long N, float thres, double *out1, double *out2);
float *rad2im(vect * in, int nx, int ny, int *ori, float aspect);
void hanming(float *res, int d1, int d2, float alpha);
void frachamming(float *res, int nx, int ny, float frac, int maxk, int *sh);
float *ellmask(int nx, int ny, int *ori, int len_x, int len_y);
void assmask(float *res, int nx, int ny, float frac);

#endif				/*      KISIP_SPECKLE_MATH_H    */
