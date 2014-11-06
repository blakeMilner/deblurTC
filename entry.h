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
 *  last changes: 13.10.2003
 *
 *  ---------------------------------------------
 *  Header for the KI SIP v6 Speckle Library in C
 *  ---------------------------------------------
 *
 *  The GUI uses GraphApp.
 *  Internet:  http://enchantia.com/software/graphapp
 *
 *  This is a parallized version using MPI.
 *  An MPI Library can be downloaded at various places,
 *  since MPI is just defines a standard for functions.
 *
 *  Download for example at
 *  Internet:  http://www-unix.mcs.anl.gov/mpi/mpich
 *
 */

#ifndef KISIP_ENTRY_H
#define	 KISIP_ENTRY_H

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include "mpifuncts.h"

#ifndef GO
#define	 GO	1
#endif

#ifndef NOGO
#define	 NOGO 0
#endif

#endif				/*      KISIP_ENTRY_H   */
