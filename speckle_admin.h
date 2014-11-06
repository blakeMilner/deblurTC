/* $Id$
 *
 *  (C) 2003-2007, Friedrich Woeger <woeger@kis.uni-freiburg.de>,
 *  Kiepenheuer-Institut f??r Sonnenphysik, Freiburg (Breisgau), Germany
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
**	speckle administration routines
**
**	=====================================================================
**
**	uses:		fileops routines
**
**	=====================================================================
**
**	Author:		F. Woeger
**			Kiepenheuer-Institut fuer Sonnenphysik
**			Freiburg, Germany
**
**	Written 13. 10. 2003
**
**	=====================================================================
*/

#ifndef KISIP_SPECKLE_ADMIN_H
#define KISIP_SPECKLE_ADMIN_H

#include <string.h>
#include <math.h>
#include <getopt.h>

#include "inlines.h"
#include "speckle_math.h"

/* max string length for reading data from file */
#define  ZEILENLAENGE   128
#define  MAXLINES   20
#define  LL 1024
#define  BS 100

/* Flags for calculations */
#define  GO   1
#define  NOGO 0

typedef struct maininfo {
    /* names of image files */
    char **imagefiles;
    /* name of result file */
    char savefiles[ZEILENLAENGE];

    int headeroffset;

    /* data properties */
    int xsize, ysize;		/* x/y size in px of data     */
    int nrofframes;			/* number of frames in burst  */

    /* properties for the reconstruction */
    struct {
		int ssizex, ssizey;	/* subfield size in px in x/y */
		int nsizex, nsizey;	/* overall size reconstr. x/y */
		int nfrx, nfry;		/* nr of subframes in x/y     */
		int offx, offy;		/* offset border in px in x/y */
    } subfields;

    /* parameters for triple correlation */
	struct {
		float sfs;			/* sfs in arcsec            */
		float lim1;			/* max length of u in x-dir in %age of diff. lim. */
		float lim2;			/* max length of u, v, u+v in %age of diff. lim. */
		int max_iterations;	/* max number of iterations */
		float snr;			/* SNR threshold for rec    */
		float eps;			/* weighting factor for rec */
		float limApod;		/* Percantage of reconstr.  apodisation */
		int max_rad;
		int l1, l2;
	} tc;
} maininfo;

/* IF YOU ALTER THIS, YOU MUST ADJUST MPI_SETSLAVETYPE() ACCORDINGLY OR ELSE  *
 * YOU WILL EXPERIENCE HORRIBLY CONFOUNDING ERRORS!                           */
typedef struct slaveinfo {
    int position;			/* subfield position          */
    int sfsizex, sfsizey;	/* subfield size in x/y       */
    int nrofframes;			/* number of frames in burst  */
    int max_rad;			/* phase int/it radius in px  */
    int max_it;				/* max number of iterations   */
    int l1, l2;				/* Limits for bispectrum      */
    float rad_x, rad_y;		/* DiffLim Radius in x/y      */
    float snr, eps;			/* params for triple corr.    */
    float limApod;			/* %age of apod.win for phrec */
} slaveinfo;

int getinfo(maininfo * result, int argc, char** argv);
void calc_subfields(int szx, int szy, float sfs, maininfo * result);
void subfielding(float *in, int count, maininfo * info, float *out);
void putslaveinf(slaveinfo * out, maininfo * in, int i);
void reconstruct(float *srec, maininfo * in, float *out);

#endif				/*      KISIP_SPECKLE_ADMIN_H   */
