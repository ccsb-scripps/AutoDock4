/********************************************************************
     Call_glss:  Invokes a GA-LS hybrid to try and solve the
                 docking problem.

                                rsh 9/95
********************************************************************/
/*
** $Log: call_glss.cc,v $
** Revision 1.1  2001/08/13 22:05:52  gillet
** *** empty log message ***
**
*/

#include <iostream.h>
#include "gs.h"
#include "ls.h"
#include "support.h"
#include "eval.h"
#include "hybrids.h"

   #include "constants.h"
   #include "structs.h"
   extern FILE *logFile;

int global_ntor;

Eval evaluate;

Representation **generate_R(int num_torsions, double xlo, double xhi, double ylo, double yhi, double zlo, double zhi)
{
   Representation **retval;

   retval = new Representation *[5];
   retval[0] = new RealVector(1, xlo, xhi);
   retval[1] = new RealVector(1, ylo, yhi);
   retval[2] = new RealVector(1, zlo, zhi);
   retval[3] = new RealVector(3);
   retval[4] = new RealVector(num_torsions+1);

   return(retval);
}

Genotype generate_Gtype(int num_torsions, double xlo, double xhi, double ylo, double yhi, double zlo, double zhi)
{
   Genotype temp(5, generate_R(num_torsions, xlo, xhi, ylo, yhi, zlo, zhi));

   return(temp);
}

Phenotype generate_Ptype(int num_torsions, double xlo, double xhi, double ylo, double yhi, double zlo, double zhi)
{
   Phenotype temp(5, generate_R(num_torsions, xlo, xhi, ylo, yhi, zlo, zhi));

   return(temp);
}

Individual random_ind(int num_torsions, double xlo, double xhi, double ylo, double yhi, double zlo, double zhi)
{
   Genotype temp_Gtype;
   Phenotype temp_Ptype;

   temp_Gtype = generate_Gtype(num_torsions, xlo, xhi, ylo, yhi, zlo, zhi);
   temp_Ptype = generate_Ptype(num_torsions, xlo, xhi, ylo, yhi, zlo, zhi);

   Individual temp(temp_Gtype, temp_Ptype);

   return(temp);
}

State call_glss(Global_Search *global_method, Local_Search *local_method, 
                State now, unsigned int num_evals, unsigned int pop_size, 
                float xlo, float xhi, float ylo, 
                float yhi, float zlo, float zhi,
                int outlev, unsigned int extOutputEveryNgens,
                Molecule *mol, Boole B_template)
{
    register int i;
    int num_iterations = 0, num_loops = 0, allEnergiesEqual = 1, numTries = 0;
    double firstEnergy = 0.0;
    EvalMode localEvalMode = Normal_Eval;

    global_method->reset(extOutputEveryNgens);
    local_method->reset();
    evaluate.reset();

    (void)fprintf( logFile, "\nCreating a population of %u individuals.\n", pop_size);
    Population thisPop(pop_size);

    (void)fprintf( logFile, "\nAssigning a random translation, a random orientation and %d random torsions to each of the %u individuals.\n\n", now.ntor, pop_size);
    global_ntor = now.ntor;//debug
    do {
        ++numTries;
        // Create a population of pop_size random individuals...
        for (i=0; i<pop_size; i++) {
            thisPop[i] = random_ind( now.ntor, double(xlo), double(xhi), double(ylo), double(yhi), double(zlo), double(zhi));
            thisPop[i].mol = mol;
            thisPop[i].age = 0L;
        }
        // Now ensure that there is some variation in the energies...
        firstEnergy = thisPop[0].value(localEvalMode);
        for (i=1; i<pop_size; i++) {
             allEnergiesEqual = allEnergiesEqual && (thisPop[i].value(localEvalMode) == firstEnergy);
        }
        if (allEnergiesEqual) {
            (void)fprintf(logFile,"NOTE: All energies are equal in population; re-initializing. (Try Number %d)\n", numTries);
        }
    } while (allEnergiesEqual);

    if (outlev > 2) { thisPop.printPopulationAsStates(logFile, pop_size, now.ntor); }

    (void)fprintf( logFile, "Beginning Lamarckian Genetic Algorithm (LGA), with a maximum of %u\nenergy evaluations.\n\n", num_evals);

    do {
        if (outlev > 1) { (void)fprintf( logFile, "Global-local search iteration %d,  Starting global search.\n", ++num_loops); }
        global_method->search(thisPop);
        if (outlev > 1) { (void)fprintf( logFile, "\tEnding global search.\n"); }
        if (outlev > 2) { thisPop.printPopulationAsStates(logFile, pop_size, now.ntor); }
        if (outlev > 3) { minmeanmax( logFile, thisPop, ++num_iterations ); }

        if (outlev > 1) { (void)fprintf( logFile, "\tStarting local search.\n"); }
        for (i=0; i<pop_size; i++) {
            local_method->search(thisPop[i]);
        }
        if (outlev > 1) { (void)fprintf( logFile, "\tEnding local search.\n"); }
        if (outlev > 2) { thisPop.printPopulationAsStates(logFile, pop_size, now.ntor); }
    } while ((evaluate.evals() < num_evals) && (!global_method->terminate()));

    thisPop.msort(3);
    return( thisPop[0].state(now.ntor) );
}
