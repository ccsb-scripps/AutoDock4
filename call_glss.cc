/********************************************************************
     Call_glss:  Invokes a GA-LS hybrid to try and solve the
                 docking problem.

                                rsh 9/95
********************************************************************/
/*
** $Log: call_glss.cc,v $
** Revision 1.3  2002/04/17 06:04:01  garrett
** Introduced the int "indiv" which can later be used to set arbitrary
** individuals in a population, not just the first one.
**
** Also added the mapping() to ensure the phenotype matches the genotype.
** This was causing an erroneous report of the initial state variables
** in the very first listing of the top %d individuals in the population
** when outlev was set to 3.                                       (gmm)
**
** Revision 1.2  2002/04/17 00:40:26  lindy
** added Booleans to call_glss for tran0, quat0, dihe0
** dpf parameters. Initilized the thisPop[0] to the
** values placed in the State now.
**
** Revision 1.1.1.1  2001/08/13 22:05:52  gillet
**  import initial of autodock sources
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

#ifdef FALSE
Individual set_ind(int num_torsions, double xlo, double xhi, double ylo, double yhi, double zlo, double zhi, State state)
{
   Genotype temp_Gtype;
   Phenotype temp_Ptype;
   int i;

   temp_Gtype = generate_Gtype(num_torsions, xlo, xhi, ylo, yhi, zlo, zhi);
   temp_Ptype = generate_Ptype(num_torsions, xlo, xhi, ylo, yhi, zlo, zhi);

   // use the state to generate a Genotype
   temp_Gtype.write(state.T.x, 0);
   temp_Gtype.write(state.T.y, 1);
   temp_Gtype.write(state.T.z, 2);
   temp_Gtype.write(state.Q.nx, 3);
   temp_Gtype.write(state.Q.ny, 4);
   temp_Gtype.write(state.Q.nz, 5);
   temp_Gtype.write(state.Q.ang, 6);
   for (i=0;i<state.ntor; i++) {
       temp_Gtype.write(state.tor[i], 7+i);
   };

   Individual temp(temp_Gtype, temp_Ptype);   

   // use mapping to generate a Phenotype
   temp.phenotyp =  temp.mapping();

   return(temp);
}
#endif

State call_glss(Global_Search *global_method, Local_Search *local_method, 
                State now, unsigned int num_evals, unsigned int pop_size, 
                float xlo, float xhi, float ylo, 
                float yhi, float zlo, float zhi,
                int outlev, unsigned int extOutputEveryNgens,
                Molecule *mol, Boole B_template,
		Boole B_RandomTran0, Boole B_RandomQuat0, Boole B_RandomDihe0)
{
    register int i;
    int num_iterations = 0, num_loops = 0, allEnergiesEqual = 1, numTries = 0;
    int indiv = 0; // Number of Individual in Population to set initial state variables for.
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

	// If initial values were supplied, put them in thisPop[0]
	if (!B_RandomTran0) {
          if (outlev > 1) { (void)fprintf(logFile, "Setting the initial translation (tran0) for individual number %d to %lf %lf %lf\n\n", indiv+1, now.T.x, now.T.y, now.T.z); }
	  thisPop[indiv].genotyp.write( now.T.x, 0);
	  thisPop[indiv].genotyp.write( now.T.y, 1);
	  thisPop[indiv].genotyp.write( now.T.z, 2);
	  // Remember to keep the phenotype up-to-date
	  thisPop[indiv].phenotyp = thisPop[indiv].mapping();
	};
	if (!B_RandomQuat0) {
          if (outlev > 1) { (void)fprintf(logFile, "Setting the initial quaternion (quat0) for individual number %d to %lf %lf %lf  %lf deg\n\n", indiv+1, now.Q.nx, now.Q.ny, now.Q.nz, Deg(now.Q.ang)); }
	  thisPop[indiv].genotyp.write( now.Q.nx, 3);
	  thisPop[indiv].genotyp.write( now.Q.ny, 4);
	  thisPop[indiv].genotyp.write( now.Q.nz, 5);
	  thisPop[indiv].genotyp.write( now.Q.ang, 6);
	  // Remember to keep the phenotype up-to-date
	  thisPop[indiv].phenotyp = thisPop[indiv].mapping();
	};
	if (!B_RandomDihe0) {
          if (outlev > 1) { (void)fprintf(logFile, "Setting the initial torsions (dihe0) for individual number %d to ", indiv+1); }
	  for (i=0;i<now.ntor; i++) {
	    thisPop[indiv].genotyp.write( now.tor[i], 7+i);
            if (outlev > 1) { (void)fprintf(logFile, "%lf ", Deg(now.tor[i])); }
	  };
          if (outlev > 1) { (void)fprintf(logFile, " deg\n\n"); }
	  // Remember to keep the phenotype up-to-date
	  thisPop[indiv].phenotyp = thisPop[indiv].mapping();
	};

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
