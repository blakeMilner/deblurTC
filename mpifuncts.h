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
 *  -------------------------------------------
 *  Header for the KI SIP v6 MPI Functions in C
 *  -------------------------------------------
 *
 *  uses mpi.h for MPI functions
 *  uses speckle.h for the actual calculation functions
 *  uses visual.h for visualizing the data
 *
 */

#ifndef KISIP_MPIFUNCTS_H
#define KISIP_MPIFUNCTS_H

#include <mpi.h>
#include <time.h>
#include "speckle_admin.h"
#include "speckle_core.h"
#include "speckle_math.h"
#include "fileops.h"

#ifndef GO
#define  GO   1
#endif
#ifndef NOGO
#define  NOGO 0
#endif
#ifndef KT
#define  KT 10
#endif
#ifndef TC
#define  TC 11
#endif
#ifndef NOISE
#define  NOISE 12
#endif

void mpi_setmoddata(int number, int rank);
void mpi_setslavetype(MPI_Datatype * typ);
void mpi_mastering(int argc, char **argv);
void mpi_master_rec(float *alldata, maininfo * info, int used_size, float *amp, float *phs);
void mpi_shutdown();
int mpi_slaving();

// These varaibles are declared here and defined in "mpi_setmoddata"
extern int proc_id;		/* Rank of processor               */
extern int proc_nr;		/* Number of available processors  */
extern char ptmodel[];		/* Path to model                   */

#endif				/* KISIP_MPIFUNCTS_H    */
