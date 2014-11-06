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


void subfielding(float *in, int count, maininfo * info, float *out)
{
    int nsx, nsy;		/* Subfield size in x/y                 */
    int nfrx;		/* No. of subfields in x/y              */
    int ox, oy;			/* offsets in x/y                       */
    int Nx, Ny;			/* Source image size in x/y             */
    int frames;			/* Nr. of frames in burst               */
    int i, j;			/* pixel index                          */
    int f;			/* frame counter index                  */
    int p;			/* helper index                         */

    nsx = (*info).subfields.ssizex;
    nsy = (*info).subfields.ssizey;
    nfrx = (*info).subfields.nfrx;
    ox = (*info).subfields.offx;
    oy = (*info).subfields.offy;
    Nx = (*info).xsize;
    Ny = (*info).ysize;
    frames = (*info).nrofframes;

    /* raw cutting */
	for (f = 0; f < frames; f++) {
		for (j = 0; j < nsy; j++) {
			for (i = 0; i < nsx; i++) {
				sf_idx(i, j, (int) count % nfrx, (int) count / nfrx, f,
				   (int) nsx / 2, (int) nsy / 2, ox, oy, Nx, Ny, p);
				out[j * nsx + i + f * nsx * nsy] = in[p];   // p is new pixel position
			}
		}
	}
}
