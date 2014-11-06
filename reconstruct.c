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
** speckle administration routines
**
** =====================================================================
**
** uses:  fileops routines
**
** =====================================================================
**
** Author:  F. Woeger
**   Kiepenheuer-Institut fuer Sonnenphysik
**   Freiburg, Germany
**
** Written 13. 10. 2003
**
** =====================================================================
*/

#include "speckle_admin.h"


void reconstruct(float *srec, maininfo * in, float *out)
{
    long i, j;
    long nsf, l;
    long sfN, N;
    long p, q;
    float *weights;
    float *win;

    N = (*in).xsize * (*in).ysize;
    sfN = (*in).subfields.ssizex * (*in).subfields.ssizey;
    nsf = (*in).subfields.nfrx * (*in).subfields.nfry;
    win = malloc(sfN * sizeof(float));
    weights = calloc(N, sizeof(float));

    assmask(win, (*in).subfields.ssizex, (*in).subfields.ssizey, (*in).tc.limApod);

    for (l = 0; l < nsf; l++) {         // total number of subfields
                for (j = 0; j < (*in).subfields.ssizey; j++) {     // number of pixels in subfield
                        for (i = 0; i < (*in).subfields.ssizex; i++) {
                                sf_idx(i, j,
                                           (int) (l % (*in).subfields.nfrx),
                                           (int) (l / (*in).subfields.nfrx), 0,
                                           (int) (*in).subfields.ssizex / 2,
                                           (int) (*in).subfields.ssizey / 2,
                                           (*in).subfields.offx, (*in).subfields.offy,
                                           (*in).xsize, (*in).ysize, q);
                                idx(i, j, (*in).subfields.ssizex, p);     // jump to pixel-position in subfield image
                                out[q] += srec[p + l * sfN] * win[p];     // add result
                                weights[q] += win[p];
                        }
                }
    }


    // divide through by weight and set border pixels to zero
    for (i = 0; i < N; i++) {
        if( weights[i] != 0 ){
                out[i] /= weights[i];
        }
	else{
		out[i] = 0.0;
	}
    }

    free(weights);
    free(win);

}
