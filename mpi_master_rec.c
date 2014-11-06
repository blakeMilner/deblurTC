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
// mpi_master_rec:
//---------------------------------------------------------------------------
void mpi_master_rec(float *alldata, maininfo * info, int used_size, float *amp, float *phs)
{
    //  Declare variables
    MPI_Datatype MPI_SLAVEINF;	/* Dataype for MPI communications       */
    MPI_Status *status;		/* Recv status handler                  */
    slaveinfo sl_info;		/* Information structure for slaves     */

    float *data;		/* subfield data burst                  */
    int jobs, jobs_done;	/* # of jobs, # of jobs done            */
    int sfN;			/* number of pixels in a subfield       */
    int i;			/* counter for loops                    */

    putslaveinf(&sl_info, info, 0);

    // Set Slaveinfo datatype for MPI sends and receives
    mpi_setslavetype(&MPI_SLAVEINF);

    // Initialise status variable for MPI data exchange
    status = (MPI_Status *) malloc(sizeof(MPI_Status));

    // Initialise variables
    jobs = (*info).subfields.nfrx * (*info).subfields.nfry;
    sfN = (*info).subfields.ssizex * (*info).subfields.ssizey;
    data = (float *) malloc(sfN * (*info).nrofframes * sizeof(float));

	/*************************************
     **** INITIALIZE SLAVES WITH WORK ****
     *************************************/

    for (i = 1; i < used_size; i++) {
    	// divide data into subfields
		subfielding(alldata, i - 1, info, data);
		// now set the data for the slaves
		putslaveinf(&sl_info, info, i - 1);

	    // Send information to the slaves
	    MPI_Send(&sl_info, 1, MPI_SLAVEINF, i, TC, MPI_COMM_WORLD);

	    MPI_Send(data, (*info).subfields.ssizex * (*info).subfields.ssizey *
		     (*info).nrofframes, MPI_FLOAT, i, TC, MPI_COMM_WORLD);
    }

	/********************************************************
     **** CONTINUE GIVING TASKS TO SLAVES FOR COMPLETION ****
     ********************************************************/

    // set counter for jobs done
    jobs_done = used_size - 1;
    //  if we are not done, send next job to slave that responds first
    while (jobs_done < jobs) {
	// receive tracked data info
		MPI_Recv(&sl_info, 1, MPI_SLAVEINF, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, status);

		/**** MPI_TAG indicates whether we are at the end of the process, ****
		 **** or just getting the mean image and tracking data ***************/
		if (status->MPI_TAG != GO) {
			// receive processed data
			MPI_Recv(&amp[sl_info.position * sfN], sfN, MPI_FLOAT,
				 status->MPI_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD,
				 status);
			MPI_Recv(&phs[2 * sl_info.position * sfN], 2 * sfN, MPI_FLOAT,
				 status->MPI_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD,
				 status);


			// Cut next subfield to be analyzed
			subfielding(alldata, jobs_done, info, data);
			// set the data for the slaves
			putslaveinf(&sl_info, info, jobs_done);

			// Send information to the slaves
			MPI_Send(&sl_info, 1, MPI_SLAVEINF, status->MPI_SOURCE, TC, MPI_COMM_WORLD);

			// send data to slaves
			MPI_Send(data, sfN * (*info).nrofframes, MPI_FLOAT,
				 status->MPI_SOURCE, TC, MPI_COMM_WORLD);

			// increase jobs_done counter
			jobs_done += 1;
		}
    }

	/*********************************************************************************
     **** WAIT FOR THE REST OF THE SLAVES TO COMPLETE AFTER THERE IS NO MORE WORK ****
     *********************************************************************************/

    // collect the remaining data
    i = 1;
    while (i < used_size) {
		// receive tracked data info
		MPI_Recv(&sl_info, 1, MPI_SLAVEINF, MPI_ANY_SOURCE, MPI_ANY_TAG,
			 MPI_COMM_WORLD, status);

		/**** MPI_TAG indicates whether we are at the end of the process, ****
		 **** or just getting the mean image and tracking data ***************/
		if (status->MPI_TAG != GO) {
			// receive processed data
			MPI_Recv(&amp[sl_info.position * sfN], sfN, MPI_FLOAT,
				 status->MPI_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, status);
			MPI_Recv(&phs[2 * sl_info.position * sfN], 2 * sfN, MPI_FLOAT,
				 status->MPI_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, status);

			// send idle signal to slaves
			MPI_Send(&sl_info, 1, MPI_INT, status->MPI_SOURCE, NOGO, MPI_COMM_WORLD);

			// set counter for slaves finished off
			i += 1;
		}
    }

    free(data);
    free(status);
    MPI_Type_free(&MPI_SLAVEINF);
}
