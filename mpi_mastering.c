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


/*
 *  KI SIP Graphical User Interface
 *  -------------------------------
 *  The Graphical User Interface for the KI SIP v6
 *  Speckle Library in C
 *
 *  created by Friedrich Woeger
 *  Kiepenheuer-Institut fuer Sonnenphysik, Freiburg
 *  last changes: 21.07.2003
 *
 *  ----------------------------------------------
 *  Main functions for the MPI control of jobs for
 *  the KI SIP v6 Speckle Library in C
 *  ----------------------------------------------
 *
 *  This is a parallized version, implemented using MPI.
 *  An MPI Library can be downloaded at various places,
 *  since MPI has standardised functions, e.g.
 *  Internet:  http://www-unix.mcs.anl.gov/mpi/mpich
 *
 */

#include "mpifuncts.h"
#include "opts.h"

//---------------------------------------------------------------------------
// mpi_mastering:
//---------------------------------------------------------------------------
void mpi_mastering(int argc, char **argv)
{
    // +++ Declare variables +++
    // MPI variables
    int jobs;				/* # of jobs                            */
    int used_size;			/* # of processors used                 */

    // data variables
    float *alldata = NULL;		/* data variables                       */
    float *famp = NULL;			/* result of calculations (amplitudes)  */
    float *fphs = NULL;			/* result of calculations (phases)      */
    float *result = NULL;		/* result of calculations (transf.)     */
    float *rec = NULL;			/* assembled reconstruction             */
    float *flipped_recon = NULL;
    long sfN;				/* number of pixels in a subfield       */
    int i, j;				/* counters for loops                   */

#ifdef DEBUG
    time_t start, end;			/* timer values				*/
#endif

    maininfo info;

//---------------------------------------------------------------------------
//	Get the vital data for processing from the GUI entries.
//	Set total number of jobs equal to the number of subfields
//	and set the number of processors used.
//
//	Attempt to read in image data. If successful, start the
//	reconstruction process.
//---------------------------------------------------------------------------


    if( getinfo(&info, argc, argv) ){
    	alldata = malloc(info.xsize * info.ysize * info.nrofframes * sizeof(float));
    }

    if ( (alldata != NULL) && readims(info.imagefiles, &info, alldata) ) {

    	calc_subfields(info.xsize, info.ysize, info.tc.sfs, &info);

	// assume maximum reconstruction radius if no command line option given
        if( info.tc.max_rad == DEFAULT_MAXRAD ){
                info.tc.max_rad = (info.subfields.ssizex / 2 - 2) ;
        }

	jobs = info.subfields.nfrx * info.subfields.nfry;
	sfN = info.subfields.ssizex * info.subfields.ssizey;
	if (jobs < proc_nr - 1){
		used_size = (jobs + 1);
	}
	else{
		used_size = proc_nr;
	}

//---------------------------------------------------------------------------
//	Allocate memory shared with slaves:
//	alldata - the entire input time series, original data format
//	famp - Fourier amplitudes. Real valued, subfield times jobs
//	fphs - Fourier phases. Real valued, subfield times jobs
//	result - reconstructed image. subfield times jobs
//	rec - assembled reconstruction. Original image size.
//---------------------------------------------------------------------------
	famp = (float *) malloc(jobs * sfN * sizeof(float));
	fphs = (float *) malloc(jobs * 2 * sfN * sizeof(float));
	result = (float *) malloc(jobs * sfN * sizeof(float));
	rec = (float *) malloc(info.xsize * info.ysize * sizeof(float));
	flipped_recon = malloc(info.xsize * info.ysize * sizeof(float));
//---------------------------------------------------------------------------
//	File exists and contains enough data. Start timer. Awake the slaves.
//	Broadcast # of processors used and let slave decide, whether he is
//	needed to join the calculation.
//---------------------------------------------------------------------------
#ifdef	DEBUG
	    time(&start);
	    fprintf(stderr, "\nWaking up slaves...");
	    fflush(stdout);
#endif
	    MPI_Bcast(&used_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
#ifdef	DEBUG
	    fprintf(stderr, " done.\r\n");
	    fprintf(stderr, "Starting calculation.\r\n");
#endif

//---------------------------------------------------------------------------
//	This function distributes the tasks and manages the subfield
//	reconstruction process. Slave work is done when MPI_MASTER_REC
//	returns.
//---------------------------------------------------------------------------

	    mpi_master_rec(alldata, &info, used_size, famp, fphs);

//---------------------------------------------------------------------------
//	Reassemble subimages with Fourier amplification (noise filter in af)
//---------------------------------------------------------------------------
#ifdef	DEBUG
	    fprintf(stderr, "\n\nAssembling... ");
	    fflush(stdout);
#endif

	    for (i = 0; i < jobs; i++) {

			/*  Combine Fourier phases and amplitudes of subimages
			 *  Calculate inverse Fourier transform to obtain reconstructed subimages
			 */
			assemble(info.subfields.ssizex, info.subfields.ssizey, &famp[i * sfN],
				 &fphs[(2 * i) * sfN], &result[i * sfN], info);
	    }

//---------------------------------------------------------------------------
//	Compute full field reconstruction.
//	End the timer.
//---------------------------------------------------------------------------
	if( jobs != 1 ){
		reconstruct(result, &info, rec);
	}
	else{
		rec = result;
	}

#ifdef DEBUG
	    time(&end);
	    fprintf(stderr, "done.\n\n");
	    fprintf(stderr, "Master done! Used time: %.2lf seconds\r\n", difftime(end, start));
#endif

//---------------------------------------------------------------------------
//      Flip and save the results
//---------------------------------------------------------------------------
            for(i = 0; i < info.xsize; i++){
                for(j = 0; j < info.ysize; j++){
                        flipped_recon[(j * info.xsize) + i] = rec[(info.ysize - j - 1) * info.xsize + i];
                }
            }

            savefloat(info.savefiles, (long) info.xsize * info.ysize, info, flipped_recon);

//---------------------------------------------------------------------------
//	Search for input data not successful. Loop to next file.
//---------------------------------------------------------------------------
	} else {
	    fprintf(stderr, "Read failed!\r\n");
	}

//---------------------------------------------------------------------------
//	Free Memory.
//---------------------------------------------------------------------------

    if (alldata != NULL)
    	free(alldata);
    if (famp != NULL)
    	free(famp);
    if (fphs != NULL)
    	free(fphs);
    if (rec != NULL)
    	free(rec);
    if (info.imagefiles != NULL){
    	for(i = 0; i < info.nrofframes; i++){
    		free(info.imagefiles[i]);
    	}
		free(info.imagefiles);
    }
}
