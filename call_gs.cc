/*

 $Id: call_gs.cc,v 1.6 2007/04/27 06:01:48 garrett Exp $

 AutoDock 

 Copyright (C) 1989-2007,  Garrett M. Morris, David S. Goodsell, Ruth Huey, Arthur J. Olson, 
 All Rights Reserved.

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


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/********************************************************************
     Call_gs:  Invokes a Global Searcher object on a randomly
               generated population of solution to the docking 
               problem.

				rsh 3/12/96
********************************************************************/
// possibly unnecessary // #include <iostream.h>
#include "gs.h"
#include "support.h"
#include "eval.h"
#include "hybrids.h"

   #include "constants.h"
   #include "structs.h"

extern Eval evaluate;

State call_gs(Global_Search *global_method, State now, unsigned int num_evals, unsigned int pop_size,
              Molecule *mol,
              int extOutputEveryNgens,
              GridMapSetInfo *info)
{
   register unsigned int i;

   evaluate.reset();
   global_method->reset(extOutputEveryNgens);

   Population thisPop(pop_size);

   for (i=0; i<pop_size; i++) {
      thisPop[i] = random_ind(now.ntor, info);
      thisPop[i].mol = mol;
   }

   do {
      global_method->search(thisPop);
   } while ((evaluate.evals() < num_evals) && (!global_method->terminate()));

   thisPop.msort(3);
   return( thisPop[0].state(now.ntor) );
}
