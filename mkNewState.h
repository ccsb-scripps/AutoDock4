/*

 $Id: mkNewState.h,v 1.6 2009/05/08 23:02:14 rhuey Exp $

 AutoDock 

Copyright (C) 2009 The Scripps Research Institute. All rights reserved.

 AutoDock is a Trade Mark of The Scripps Research Institute.

 This program is free software; you can redistribute it and/or
 modify it under the terms of the GNU General Public License
 as published by the Free Software Foundation; either version 2
 of the License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

 */

#ifndef MKNEWSTATE
#define MKNEWSTATE

#include "constants.h"
#include "qmultiply.h"
#include "cnv_state_to_coords.h"

void  mkNewState( State *now,
		  State *last,
		  State *change,
		  /*
		  ** Real qtn[NQTN],
                  ** Real tor[MAX_TORS],
                  ** Real qtnLast[NQTN],
                  ** Real torLast[MAX_TORS],
                  ** Real qtnChange[NQTN],
                  ** Real torChange[MAX_TORS],
		  */
                  Real vt[MAX_TORS][NTRN],
                  int   tlist[MAX_TORS][MAX_ATOMS],
                  int   ntor,
                  Real crd[MAX_ATOMS][NTRN],
                  Real crdpdb[MAX_ATOMS][NTRN],
                  int   natom,
                  Real trnStep,
                  Real qtwStep,
                  Real torStep,
	          Real F_TorConRange[MAX_TORS][MAX_TOR_CON][2],
		  int   N_con[MAX_TORS]);
#endif
