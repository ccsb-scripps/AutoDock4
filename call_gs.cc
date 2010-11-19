/*

 $Id: call_gs.cc,v 1.6.2.1 2010/11/19 20:09:30 rhuey Exp $

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
#include "pso.h"

#include "constants.h"
#include "structs.h"

extern Eval evaluate;
unsigned int maxEvalNum;
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
   } while (( (unsigned)evaluate.evals() < num_evals) && (!global_method->terminate()));

   // TSRI 20101101 changed 3 to 1 in next line, added printing of final energy M Pique
   thisPop.msort(1);
    (void)fprintf(logFile,"Final-Value: %.3f\n", thisPop[0].value(Normal_Eval));
   return( thisPop[0].state(now.ntor) );
}


//////////////////////////////////////////////////////////
// Wrapper to call PSO  
// -Huameng 07/09/2008
//////////////////////////////////////////////////////////
State call_pso(
	Global_Search *global_method, 
	State now, 
	unsigned int num_evals, 
	unsigned int pop_size,
    Molecule *mol,
    int extOutputEveryNgens,
    GridMapSetInfo *info
    )
{
	unsigned int i;
   	int allEnergiesEqual = 1, numTries = 0;
 
   	double firstEnergy = 0.0;
   	double indvEnergy = 0.0;
   	EvalMode localEvalMode = Normal_Eval;
   	evaluate.reset();
   	global_method->reset(extOutputEveryNgens);

	maxEvalNum = num_evals;
   	Population thisPop(pop_size);
   
   	fprintf(logFile, "start initializing %d particles\n", pop_size);
   	do { 	  	
   	   ++numTries;
   	   
	   // initialize particles in population
	   for (i=0; i < pop_size; i++) {
	      thisPop[i] = random_ind(now.ntor, info);	      
	      thisPop[i].mol = mol;
	      //fprintf(logFile, "Done initializing particle %d\n\n", i+1);
	      //fflush(logFile);
	   }
	   // Now ensure that there is some variation in the energies...
       firstEnergy = thisPop[0].value(localEvalMode);
           
#ifdef DEBUG
   	(void)fprintf(logFile,"\n\ncall_pso() ensuring there is variation in the energies, firstEnergy=%lf\n\n", firstEnergy);	
#endif
        for (i=1; i<pop_size; i++) {       	
			 indvEnergy = thisPop[i].value(localEvalMode);			           
             allEnergiesEqual = allEnergiesEqual && (indvEnergy == firstEnergy);                          
        }
        if (allEnergiesEqual) {
            (void)fprintf(logFile,"NOTE: All energies are equal in population; re-initializing. (Try Number %d)\n", numTries);
        }
   	} while (allEnergiesEqual);
   
   	fprintf(logFile, "Beginning PSO search... \n");
   	fflush(logFile);

   	do {
      	global_method->search(thisPop);
   	} while (( (unsigned)evaluate.evals() < num_evals) && (!global_method->terminate()));

   	//thisPop.msort(3);
   	//return( thisPop[0].state(now.ntor) ); 
   // TSRI 20101101 changed 3 to 1 in next line, added printing of final energy M Pique
   thisPop.msort(1);
    (void)fprintf(logFile,"Final-Value: %.3f\n", thisPop[0].value(Normal_Eval));
   	return (((ParticleSwarmGS *)global_method)->getBest().state(now.ntor));
}


