/*

 $Id: call_gs.cc,v 1.3 2003/02/26 00:27:20 garrett Exp $

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
              FloatOrDouble xlo, FloatOrDouble xhi, FloatOrDouble ylo, FloatOrDouble yhi, FloatOrDouble zlo, FloatOrDouble zhi, Molecule *mol,
              int extOutputEveryNgens)
{
   int i;

   evaluate.reset();
   global_method->reset(extOutputEveryNgens);

   Population thisPop(pop_size);

   for (i=0; i<pop_size; i++)
   {
      thisPop[i] = random_ind(now.ntor, double(xlo), double(xhi), double(ylo), double(yhi), double(zlo), double(zhi));
      thisPop[i].mol = mol;
   }

   do
   {
      global_method->search(thisPop);
   } while ((evaluate.evals() < num_evals) && (!global_method->terminate()));

   thisPop.msort(3);
   return( thisPop[0].state(now.ntor) );
}
