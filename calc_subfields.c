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


void calc_subfields(int szx, int szy, float sfs, maininfo * result)
{
    /* helper variables */
    int pot;			/* power for subfield size     */

    /* start subfielding if wanted */
    if (sfs != 0) {
		/* raw subfield size in x/y */
		(*result).subfields.ssizex = (int) sfs;
		(*result).subfields.ssizey = (int) sfs;

		/* rescale raw size to a power of two */
		pot = 0;
		while ((*result).subfields.ssizex > 0) {
			(*result).subfields.ssizex =
			(int) (((float) (*result).subfields.ssizex) / 2);
			pot += 1;
		}
		(*result).subfields.ssizex = pow(2, pot);

		pot = 0;
		while ((*result).subfields.ssizey > 0) {
			(*result).subfields.ssizey =
			(int) (((float) (*result).subfields.ssizey) / 2);
			pot += 1;
		}
		(*result).subfields.ssizey = pow(2, pot);

		/* find number of subfields, overall size reconstr. and border offsets */

		(*result).subfields.nfrx = (int) (2 * szx / (*result).subfields.ssizex - 1);
		(*result).subfields.nsizex = (int) ((float) ((*result).subfields.ssizex) *
			   (float) ((*result).subfields.nfrx + 1) / 2);
		(*result).subfields.offx = (int) (szx - (*result).subfields.nsizex) / 2;

		(*result).subfields.nfry = (int) (2 * szy / (*result).subfields.ssizey - 1);
		(*result).subfields.nsizey = (int) ((float) ((*result).subfields.ssizey) *
				(float) ((*result).subfields.nfry + 1) / 2);
		(*result).subfields.offy = (int) (szy - (*result).subfields.nsizey) / 2;	// offset makes the border!
		/* if not wanted set subfield size to frame size */
    } else {
		(*result).subfields.ssizex = szx;
		(*result).subfields.ssizey = szy;
		(*result).subfields.nsizex = szx;
		(*result).subfields.nsizey = szy;
		(*result).subfields.nfrx = 1;
		(*result).subfields.nfry = 1;
		(*result).subfields.offx = 0;
		(*result).subfields.offy = 0;
    }
}
