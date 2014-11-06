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
 * Look at Header file for more details
 */

#include "entry.h"

int main(int argc, char *argv[])
{

/* Shutdown Flag */
    int go_flag = GO;

/* Declare important MPI variables here					*/
    int rank;			/* Rank of processor                    */
    int number;			/* Number of available processors       */

/* Initialize MPI */
    MPI_Init(&argc, &argv);
/* Get rank */
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
/* Get the total number of processors */
    MPI_Comm_size(MPI_COMM_WORLD, &number);

/* Tell mpifuncts rank and number */
#ifdef	DEBUG
	fprintf(stderr, "Initializing MPI [processors = %3d, rank = %3d].\n",number,rank);
#endif

    mpi_setmoddata(number, rank);

/*
 * 	MASTER PROCESS
 */
    if (rank == 0) {
/* start MPI */
		mpi_mastering(argc, argv);
		mpi_shutdown();
/*
 *	SLAVES PROCESSES
 */
    } else {
		while (go_flag != NOGO) {
			go_flag = mpi_slaving();
		}
    }

/* Terminate MPI */
    MPI_Finalize();

    return 0;
}
