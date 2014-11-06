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

#ifndef	 KISIP_INLINES_H
#define KISIP_INLINES_H


/*
 *      Definitions for inline complex functions
 *          cadd - complex addition             r = a + b
 *          csub - complex subtraction          r = a - b
 *          cmul - complex multiplication       r = a * b
 *          cdiv - complex division             r = a / b
 *          conj - complex conjugate            r = a*
 *          cmod - modulus                      |a|
 *          ccpy - complex copy                 r = a
 */

#define cadd(a,b,r) \
{ \
  *(r) = *(a) + *(b); \
  *((r)+1) = *((a)+1) + *((b)+1); \
}
#define csub(a,b,r) \
{ \
  *(r) = *(a) - *(b) ; \
  *((r)+1) = *((a)+1) - *((b)+1); \
}
#define cmul(a,b,r) \
{ \
  float t1, t2; \
  t1 = *(a) ; t2 = *((a) + 1); \
  *(r) = *(b)*t1 - *((b)+1)*t2; \
  *((r)+1) = *((b)+1)*t1 + *(b)*t2; \
}
#define cdiv(a,b,r) \
{ \
  float t1, t2, denom; \
  t1 = *(b) ; t2 = *((b) + 1); \
  denom = t1*t1 + t2*t2; \
  if (denom == 0.0) { \
    *(r) = 0.0; \
    *((r)+1) = 0.0; \
  } else { \
    *(r) = (*(a)*t1 + *((a)+1)*t2)/denom; \
    *((r)+1) = (*((a)+1)*t1 - *(a)*t2)/denom; \
  } \
}
#define csqr(a,r) \
{ \
  double t1, t2, t3, t4, t5, t6, t7; \
  t1 = *(a); t2 = *((a)+1); \
  if ((t1 == 0.0) && (t2 == 0.0)) { \
    *(r) = 0.0; \
    *((r)+1) = 0.0; \
  } else { \
    if ((t3 = fabs(t1)) >= (t4 = fabs(t2))) \
      {t5 = t4/t3; t6 = sqrt(t3)*sqrt(0.5*(1.0+sqrt(1.0+t5*t5)));} \
    else \
      {t5 = t3/t4; t6 = sqrt(t4)*sqrt(0.5*(t5+sqrt(1.0+t5*t5)));} \
    if (t1 >= 0.0) \
      {*(r) = (float)t6; *(r+1) = (float)(t2/(2.0*t6));} \
    else \
      {t7 = (t2 >= 0.0) ? t6 : -t6; *(r) = (float)(t2/(2.0*t7)); *(r+1) = (float)t7;} \
  } \
}
#define conj(a,r) \
{ \
  *(r) = *(a); \
  *((r)+1) = -(*((a)+1)); \
}
#define ccpy(a,r) \
{ \
  *(r) = *(a); \
  *((r)+1) = *((a)+1); \
}
#define cmod(a) ((float)(sqrt((double)(*(a)*(*(a)) + *((a)+1)*(*((a)+1))))))

/*
 *      Definitions for pointer arithmetic
 *
 *      sf_idx:         jump to correct vector-position in image to
 *                      split a matrix into subfields
 *                      (i,j)           -       requested pixel position
 *                      (m,n)           -       subfield index
 *                      fr		-	index of frame in burst
 *                      nxh,nyh		-	half of subfield height/width
 *                      ox, oy		-	offset from border
 *                      Nx,Ny		-	height/width of source image
 *      shift_idx:      jump to shifted pixel-position in image
 *                      (i,j)           -       requested pixel position
 *                      (s_x,s_y)       -       shift-vector
 *			nx,ny		-	image height/width
 *      mbs_idx:        vector-position for MeanBiSpectrum manipulation
 *                      (i,j,m,n)       -       requested (u,v) position
 *			(l1,l2,l3)	-	limits whithin BS is calculated
 *      w_idx:          vector-position for weight matrix manipulation
 *                      (i,j,m,n)       -       requested (u,v) position
 *			(l1,l2,l3)	-	limits whithin BS is calculated
 *      phs_idx:        vector-position for Phase manipulation
 *                      (i,j)           -       requested pixel position
 *			nx		-	image width
 *      idx:            jump to pixel-position in subfield image
 *                      (i,j)           -       requested pixel position
 *			nx		-	image width
 *      idx3d:          jump to pixel-position in datacube
 *          (i,j,k) -   requested pixel position
 *			nx		-	image width
 */

  // single backslash merges lines into one long line
#define sf_idx(i, j, m, n, fr, nxh, nyh, ox, oy, Nx, Ny, r) \
{ \
  long tx, ty; \
  tx = (ox) + (m)*(nxh) + (i) + (Nx); \
  ty = (oy) + (n)*(nyh) + (j) + (Ny); \
  tx %= (Nx); \
  ty %= (Ny); \
  r = ty*(Nx) + tx + (fr)*(Nx)*(Ny); \
}
#define shift_idx(i, j, s_x, s_y, nx, ny, r) \
{ \
  long tx, ty; \
  tx = (i) + (nx) + (s_x); \
  ty = (j) + (ny) + (s_y); \
  tx = (tx)%(nx); \
  ty = (ty)%(ny); \
  r = ty*(nx) + tx; \
}
#define mbs_idx(i, j, m, n, l1, l2, l3, r) \
{ \
  r = (2*(((((n)+(l3))*(2*((l3)+1)) + (m)+(l3))*(2*((l2)+1)) + (j)+(l2))*(l1+1) + (i)));\
}
#define w_idx(i, j, m, n, l1, l2, l3, r) \
{ \
  r = (((((n)+(l3))*(2*((l3)+1)) + (m)+(l3))*(2*((l2)+1)) + (j)+(l2))*(l1+1) + (i));\
}
#define phs_idx(i, j, nx, r) \
{ \
  r = (2*((j)*(nx) + (i)));\
}
#define idx(i, j, nx, r) {r = (j)*(nx) + (i);}
#define idx3d(i, j, k, nx, ny, r) {r = ((k)*(ny)+(j))*(nx) + (i);}

#endif

