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

#define SLAVEINFO_NUM 14

//---------------------------------------------------------------------------
//      setslavetype:
//---------------------------------------------------------------------------
void mpi_setslavetype(MPI_Datatype * typ)
{
    // Declaration of the slaveinfo MPI Type here for communication purposes
    int bl[SLAVEINFO_NUM] =
	{ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
    MPI_Aint ind[SLAVEINFO_NUM];
    MPI_Datatype otyp[SLAVEINFO_NUM] =
    { 	MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT,
		MPI_INT, MPI_INT, MPI_INT,
		MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT,
		MPI_FLOAT,
		MPI_UB
    };

    ind[0] = 0;
    ind[1] = ind[0] + (MPI_Aint) sizeof(int);
    ind[2] = ind[1] + (MPI_Aint) sizeof(int);
    ind[3] = ind[2] + (MPI_Aint) sizeof(int);
    ind[4] = ind[3] + (MPI_Aint) sizeof(int);
    ind[5] = ind[4] + (MPI_Aint) sizeof(int);
    ind[6] = ind[5] + (MPI_Aint) sizeof(int);
    ind[7] = ind[6] + (MPI_Aint) sizeof(int);
    ind[8] = ind[7] + (MPI_Aint) sizeof(float);
    ind[9] = ind[8] + (MPI_Aint) sizeof(float);
    ind[10] = ind[9] + (MPI_Aint) sizeof(float);
    ind[11] = ind[10] + (MPI_Aint) sizeof(float);
    ind[12] = ind[11] + (MPI_Aint) sizeof(float);
    ind[SLAVEINFO_NUM - 1] = sizeof(slaveinfo);

    MPI_Type_struct(SLAVEINFO_NUM, bl, ind, otyp, typ);
    MPI_Type_commit(typ);
}
