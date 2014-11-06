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

//---------------------------------------------------------------------------
// mpi_slaving:
//---------------------------------------------------------------------------
int mpi_slaving()
{
    // Declare local variables
    MPI_Datatype MPI_SLAVEINF;	/* Dataype for MPI communications       */
    slaveinfo sl_info;		/* Information structure for slaves     */
    MPI_Status *status;		/* Recv status handler                  */
    int used_size;		/* # of processors used                 */

    long sfN;			/* # of pixels in subfield              */
    float *data = NULL;		/* subfield data burst                  */
    float *winH = NULL;		/* apodisation window for SR calc.      */
    float *winF = NULL;		/* apodisation window for phase rec.    */
    float *mask = NULL;		/* DiffLim Mask                         */
    float *pc = NULL;		/* Phase Consistency                    */

    /* for reallocation: needs to be NULL   */
    int *shifts = NULL;		/* shifts for KT Cross Spectra          */
    /* for reallocation: needs to be NULL   */
    int maxk;			/* number of shifts                     */

    // Triple correlation part
    long *index = NULL;		/* index list for bispectrum vectors    */
    long bs_cnt;		/* # of vectors used                    */
    float *bsc = NULL;		/* complex non-red/red bispectrum       */
    float *wc = NULL;		/* complex non-red/red bs weights       */
    float *p1 = NULL;		/* phase matrix for iterativ reconstr.  */
    float *aphs = NULL;		/* average of low frequencies           */

    // General part
    float *phs = NULL;		/* Phase of reconstructed image         */
    float *amp = NULL;		/* Amplitude of reconstructed image     */
    int c = GO;			/* helpers                              */
    long i;



    // set up status variable for sends and receives
    status = (MPI_Status *) malloc(sizeof(MPI_Status));

    // Set Slaveinfo datatype for MPI sends and receives
    mpi_setslavetype(&MPI_SLAVEINF);

    // Receive number of jobs ...
    MPI_Bcast(&used_size, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // ...and if used_size = NOGO, we are quitting the slave!
    if (used_size == NOGO) {
		// free memory
		MPI_Type_free(&MPI_SLAVEINF);
		free(status);
		// returning NOGO quits the slaves in entry.c
		return NOGO;
    }
    // ...and decide if I am needed
    if (used_size > proc_id) {
	// work until TAG says NOGO
		while (GO) {
			// Receive the data information
			MPI_Recv(&sl_info, 1, MPI_SLAVEINF, 0, MPI_ANY_TAG, MPI_COMM_WORLD, status);

			// Break if shutdown signal was sent
			if (status->MPI_TAG == NOGO)
				break;

			if (status->MPI_TAG == TC) {
			/* First time allocation of memory */
				if (c) {
					// find number of pixels in subfield
					sfN = sl_info.sfsizex * sl_info.sfsizey;
					// allocate memory for data
					data = (float *) malloc(sfN * sl_info.nrofframes * sizeof(float));
					winH = (float *) malloc(sfN * sizeof(float));
					winF = (float *) malloc(sfN * sizeof(float));
					pc = (float *) malloc(sfN * sizeof(float));

					/* Initialise hamming window, fractional hamming & mask */
					hanming(winH, sl_info.sfsizex, sl_info.sfsizey, 0.5);
					frachamming(winF, sl_info.sfsizex, sl_info.sfsizey, sl_info.limApod, 0, NULL);
					// mask is no longer used...
					mask = ellmask(sl_info.sfsizex, sl_info.sfsizey, NULL, sl_info.rad_x, sl_info.rad_y);

					// allocate appropriate memory
					index = bs_init(&bs_cnt, sl_info);

					bsc = (float *) malloc(2 * bs_cnt * sizeof(float));
					wc = (float *) malloc(bs_cnt * sizeof(float));

					// Allocate memory for reconstructed amplitudes & phases
					amp = (float *) malloc(sfN * sizeof(float));
					phs = (float *) malloc(2 * sfN * sizeof(float));
					p1 = (float *) malloc(2 * sfN * sizeof(float));

					// set flag so as to not allocate new memory later
					c = NOGO;
				}
			// initialise amps, phases & phase consistency to zero
				memset(bsc, 0.0, 2 * bs_cnt * sizeof(float));
				memset(wc, 0.0, bs_cnt * sizeof(float));
				memset(amp, 0.0, sfN * sizeof(float));
				memset(phs, 0.0, 2 * sfN * sizeof(float));
				memset(p1, 0.0, 2 * sfN * sizeof(float));
				memset(pc, 0.0, sfN * sizeof(float));

				// finally receive the data
				MPI_Recv(data, sfN * sl_info.nrofframes, MPI_FLOAT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, status);

				// compute mean of burst
				//mean(data, sl_info.sfsizex, sl_info.sfsizey, sl_info.nrofframes, temp);

				/* compute the position of the 'good' average phase parts
				 	 3 is an argument here because we deleted:
				   ((int) (alpha[0] * i_rad) > 3 ?
					(int) (alpha[0] * i_rad) : 3),  	*/
				init_shift(sl_info.sfsizex / 2, sl_info.sfsizey / 2, 3, &maxk, &shifts);

				// create bi-spectrum/weights in non-redundant matrix
				aphs = bs_ave(data, winF, sl_info, index, bs_cnt, bsc, wc, amp, maxk, shifts);

				// set snr-threshold
				bs_snrt(wc, bs_cnt, &sl_info);

				// phase reconstruction
				//   init phase matrices
				phs_init(phs, p1, pc, sl_info, maxk, shifts, aphs);

				//   recursive approach
				rpr(phs, p1, pc, index, bs_cnt, bsc, wc, sl_info);

				//   iterative approach
				for (i = 0; i < sl_info.max_it; i++) {
					iwlspr(phs, p1, pc, bsc, wc, index, bs_cnt, sl_info, maxk);
					if (chkphase(sl_info, phs) < 1.0e-5)
						break;
				}

				// Send back info, so master knows, which subfield burst this was
				// MPI_TAG will indicate whether we are at the end of the process
				MPI_Send(&sl_info, 1, MPI_SLAVEINF, 0, NOGO, MPI_COMM_WORLD);

				// Send back processed data
				MPI_Send(amp, sfN, MPI_FLOAT, 0, NOGO, MPI_COMM_WORLD);
				// send back phase
				MPI_Send(phs, 2 * sfN, MPI_FLOAT, 0, NOGO, MPI_COMM_WORLD);

				// free some memory
				free(shifts);
				shifts = NULL;
				free(aphs);
				aphs = NULL;
			}
		}

		// free memory
		if (data != NULL)
			free(data);
		if (winH != NULL)
			free(winH);
		if (winF != NULL)
			free(winF);
		if (mask != NULL)
			free(mask);
		if (pc != NULL)
			free(pc);
		if (shifts != NULL)
			free(shifts);
		if (amp != NULL)
			free(amp);
		if (phs != NULL)
			free(phs);
		if (index != NULL)
			free(index);
		if (p1 != NULL)
			free(p1);
		if (bsc != NULL)
			free(bsc);
		if (wc != NULL)
			free(wc);
		
		free(status);
		MPI_Type_free(&MPI_SLAVEINF);
		// idle until new data or shutdown Broadcast is sent
		return GO;
    }
    else {
		// free memory
		free(status);
		MPI_Type_free(&MPI_SLAVEINF);
		// idle until new data or shutdown Broadcast is sent
		return GO;
    }
}
