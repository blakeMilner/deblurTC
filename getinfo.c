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
#include "opts.h"

int getinfo(maininfo* out, int argc, char** argv)
{
	char ch;
	int i;
	int xPresent = 0;
	int yPresent = 0;

	/* Create default values for when no argument is given */
	strcpy(out->savefiles, 		DEFAULT_SAVEFILE);
	out->xsize = 				DEFAULT_XSIZE;
	out->ysize = 				DEFAULT_YSIZE;
	out->headeroffset = 		DEFAULT_HEADERSIZE;
	out->tc.sfs = 				DEFAULT_SFS;
	out->tc.l1 = 				DEFAULT_BSAXISLENGTH;
	out->tc.l2 = 				DEFAULT_BSLENGTH;
	out->tc.max_rad = 			DEFAULT_MAXRAD;
	out->tc.max_iterations =	DEFAULT_MAXITER;
	out->tc.snr = 				DEFAULT_SNR;
	out->tc.eps = 				DEFAULT_EXP;
	out->tc.limApod = 			DEFAULT_APODLIM;

	while ((ch = getopt_long(argc, argv, "oxyfsvwrpnea:", options, NULL)) != -1){
		 switch (ch) { // make most of these not required...
			 case 'o':
				 strcpy(out->savefiles, optarg);
				 break;
			 case 'x':
				 out->xsize = atoi(optarg);
				 xPresent = 1;
				 break;
			 case 'y':
				 out->ysize = atoi(optarg);
				 yPresent = 1;
				 break;
			 case 'h':
				 out->headeroffset = atoi(optarg);
				 break;
			 case 's':
				 out->tc.sfs = (float) atof(optarg);
				 break;
			 case 'v':
				 out->tc.l1 = (int) atoi(optarg);
				 break;
			 case 'w':
				 out->tc.l2 = (int) atoi(optarg);
				 break;
			 case 'r':
				 out->tc.max_rad = (int) atoi(optarg);
				 break;
			 case 'p':
				 out->tc.max_iterations = atoi(optarg);
				 break;
			 case 'n':
				 out->tc.snr = (float) atof(optarg) / 100.0;
				 break;
			 case 'e':
				 out->tc.eps = (float) atof(optarg);
				 break;
			 case 'a':
				 out->tc.limApod = (float) atof(optarg) / 100.0;
				 break;
		 }
	 }
	 argc -= optind;
	 argv += optind;

	 if(argc == 0){
		 fprintf(stderr, "ERROR: no image files provided.\n");
		 return(NOGO);
	 }
	 else{
		 out->imagefiles = malloc(argc * sizeof(char*));
	 }

	 // interpret arguments without options as input image files
	 out->nrofframes = argc;
	 for(i = 0; i < out->nrofframes; i++){
		 out->imagefiles[i] = malloc(ZEILENLAENGE * sizeof(char));
		 strcpy(out->imagefiles[i], argv[i]);
	 }

#ifndef EMAN2
	 if( !xPresent || !yPresent ){
		 fprintf(stderr, "ERROR: pixels in x- and y-direction required.\n");
		 return(NOGO);
	 }
#else
	 if( xPresent || yPresent ){
		 fprintf(stderr, "Pixels in x- and y-direction not needed, ignored.\n");
	 }
#endif

	 return(GO);
}
