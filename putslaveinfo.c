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


void putslaveinf(slaveinfo * out, maininfo * in, int i)
{
    /* set position of subfield in whole image */
    (*out).position = i;
    /* set subfield size */
    (*out).sfsizex = (*in).subfields.ssizex;
    (*out).sfsizey = (*in).subfields.ssizey;
    /* set # of frames */
    (*out).nrofframes = (*in).nrofframes;

    /* set reconstr. limit (don't go too far) */
    (*out).rad_x = (*in).subfields.ssizex / 2 - 2;
    if ((*out).rad_x > (*in).subfields.ssizex / 2 - 2)
    	(*out).rad_x = (int) (*in).subfields.ssizex / 2 - 2;
    (*out).rad_y = (*in).subfields.ssizey / 2 - 2;
    if ((*out).rad_y > (*in).subfields.ssizey / 2 - 2)
    	(*out).rad_y = (int) (*in).subfields.ssizey / 2 - 2;

	/* set maximum reconstruction settings */
    (*out).max_rad = (*in).tc.max_rad;
	// *** check if this is too far
	if ((*out).max_rad > (*out).rad_x)
	    (*out).max_rad = (*out).rad_x;
	if ((*out).max_rad > (*out).rad_y)
	    (*out).max_rad = (*out).rad_y;

	// ****************************
	(*out).l1 = (*in).tc.l1; // hardcode
	// *** check if this is too far
	if ((*out).l1 > (*out).rad_x)
	    (*out).l1 = (*out).rad_x;
	if ((*out).l1 > (*out).rad_y)
	    (*out).l1 = (*out).rad_y;

	// ****************************
	(*out).l2 = (*in).tc.l2;
	// *** check if this is too far
	if ((*out).l2 > (*out).rad_x)
	    (*out).l2 = (*out).rad_x;
	if ((*out).l2 > (*out).rad_y)
	    (*out).l2 = (*out).rad_y;

	// ****************************
	(*out).max_it = (*in).tc.max_iterations;
	(*out).snr = (*in).tc.snr;
	(*out).eps = (*in).tc.eps;
	/* set phase reconstruction apodisation */
	(*out).limApod = (*in).tc.limApod;
}
