/********************************************************************
     Call_ls:  Invokes a local searcher on a docking to try and 
               find the locally optimal solution.  So, the docking
               must be specified BEFORE calling this routine.
               Assumes a population size of 1.

				rsh 2/5/96
********************************************************************/
// possibly unnecessary // #include <iostream.h>
#include "ls.h"
#include "support.h"
#include "eval.h"

   #include "constants.h"
   #include "structs.h"
//   extern FILE *logFile;
   #include "qmultiply.h"

extern Eval evaluate;

Representation **cnv_state_to_rep(const State &state)
{
   register int i;
   Representation **retval;

   retval = new Representation *[5];
   retval[0] = new RealVector(1);
   retval[0]->write(state.T.x, 0);
   retval[1] = new RealVector(1);
   retval[1]->write(state.T.y, 0);
   retval[2] = new RealVector(1);
   retval[2]->write(state.T.z, 0);
   retval[3] = new RealVector(3);
   retval[3]->write(state.Q.nx, 0);
   retval[3]->write(state.Q.ny, 1);
   retval[3]->write(state.Q.nz, 2);
   retval[4] = new RealVector(1+state.ntor);
   retval[4]->write(state.Q.ang, 0);
   for(i=1; i<=state.ntor; i++)
   {
      retval[4]->write(state.tor[i-1], i);
   }

   return(retval);
}

Individual cnv_state_to_ind(const State &original)
{
   // BEGIN DELETION
   // return(Individual(Genotype(5, cnv_state_to_rep(original)), Phenotype(5, cnv_state_to_rep(original))));
   // END DELETION

   // BEGIN ADDITION
   // Added by gmm, 27-MAR-97, to solve these compiler warnings:
   //
   // call_ls.cc:59: warning: In this statement, the initialization of a non-const reference requires a temporary for "Genotype(5,cnv_state_to_rep(original))". (reftemporary)
   // call_ls.cc:59: warning: In this statement, the initialization of a non-const reference requires a temporary for "Phenotype(5,cnv_state_to_rep(original))". (reftemporary)
   // call_ls.cc:59: warning: In this statement, the initialization of a non-const reference requires a temporary for "(Individual(Genotype(5,cnv_state_to_rep(original)),Phenotype(5,cnv_state_to_rep(original))))". (reftemporary)

   Genotype temp_Gtype;
   Phenotype temp_Ptype;

   temp_Gtype = Genotype(5, cnv_state_to_rep(original));
   temp_Ptype = Phenotype(5, cnv_state_to_rep(original));

   Individual temp(temp_Gtype, temp_Ptype);

   return(temp);
   // END ADDITION

}

State call_ls(Local_Search *local_method, State now, unsigned int pop_size, Molecule *mol) 
{
   int i;

   evaluate.reset();
   local_method->reset();

   Population thisPop(pop_size);
   for(i=0; i<pop_size; i++)
   {
      thisPop[i] = cnv_state_to_ind(now); 
      thisPop[i].mol = mol;
   }

   for(i=0; i<pop_size; i++)
   {
      local_method->search( thisPop[i] );
   }

   thisPop.msort(1);
   return(thisPop[0].state(now.ntor));
}
