/*

 $Id: call_glss.cc,v 1.14 2006/01/30 04:25:10 billhart Exp $

*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/********************************************************************
     Call_glss:  Invokes a GA-LS hybrid to try and solve the
                 docking problem.

                                rsh 9/95
********************************************************************/

#include <string.h>
#include "gs.h"
#include "ls.h"
#include "support.h"
#include "eval.h"
#include "hybrids.h"
#include "constants.h"
#include "structs.h"
#include "openfile.h"

extern FILE *logFile;
extern char *programname;

int global_ntor;

Eval evaluate;

Representation **generate_R(int num_torsions, GridMapSetInfo *info)
{
   Representation **retval;

#ifdef DEBUG
    // gmm 20-FEB-2003
    (void)fprintf(logFile,"call_glss.cc/Representation **generate_R()  about to create a new Representation with 5 elements, retval...\n");
#endif
   retval = new Representation *[5];
#ifdef DEBUG
    // gmm 20-FEB-2003
    (void)fprintf(logFile,"call_glss.cc/Representation **generate_R()  done creating   a new Representation with 5 elements, retval...\n");
#endif
   retval[0] = new RealVector(1, info->lo[X], info->hi[X]);
   retval[1] = new RealVector(1, info->lo[Y], info->hi[Y]);
   retval[2] = new RealVector(1, info->lo[Z], info->hi[Z]);
   retval[3] = new RealVector(3);
   retval[4] = new RealVector(num_torsions+1);
#ifdef DEBUG
    // gmm 20-FEB-2003
    (void)fprintf(logFile,"call_glss.cc/Representation **generate_R()  done assigning each of the retval[0-5] elements...\n");
#endif

   return(retval);
}

Genotype generate_Gtype(int num_torsions, GridMapSetInfo *info)
{
#ifdef DEBUG
    // gmm 20-FEB-2003
    (void)fprintf(logFile,"call_glss.cc/Genotype generate_Gtype() about to call Genotype temp(5, generate_R())...\n");
#endif
   Genotype temp((unsigned int)5, generate_R(num_torsions, info));
#ifdef DEBUG
    // gmm 20-FEB-2003
    (void)fprintf(logFile,"call_glss.cc/Genotype generate_Gtype() done calling  Genotype temp(5, generate_r())...\n");
#endif

   return(temp);
}

Phenotype generate_Ptype(int num_torsions, GridMapSetInfo *info) 
{
#ifdef DEBUG
    // gmm 20-FEB-2003
    (void)fprintf(logFile,"call_glss.cc/Genotype generate_Ptype() about to call Phenotype temp(5, generate_R())...\n");
#endif
   Phenotype temp((unsigned int)5, generate_R(num_torsions, info));
#ifdef DEBUG
    // gmm 20-FEB-2003
    (void)fprintf(logFile,"call_glss.cc/Genotype generate_Ptype() done calling  Phenotype temp(5, generate_R())...\n");
#endif

   return(temp);
}

Individual random_ind(int num_torsions,  GridMapSetInfo *info) 
{
   Genotype temp_Gtype;
   Phenotype temp_Ptype;

#ifdef DEBUG
    // gmm 20-FEB-2003
    (void)fprintf(logFile,"call_glss.cc/Individual random_ind() about to generate_Gtype()...\n");
#endif
   temp_Gtype = generate_Gtype(num_torsions, info);
#ifdef DEBUG
    // gmm 20-FEB-2003
    (void)fprintf(logFile,"call_glss.cc/Individual random_ind() about to generate_Ptype()...\n");
#endif
   temp_Ptype = generate_Ptype(num_torsions, info); // differs on Linux gcc 2.96 and Mac OS X gcc 3.1

#ifdef DEBUG
    // gmm 20-FEB-2003
    (void)fprintf(logFile,"call_glss.cc/Individual random_ind() about to Individual temp(temp_Gtype, temp_Ptype)...\n");
#endif
   Individual temp(temp_Gtype, temp_Ptype);
#ifdef DEBUG
    // gmm 20-FEB-2003
    (void)fprintf(logFile,"call_glss.cc/Individual random_ind() done     Individual temp(temp_Gtype, temp_Ptype)...\n");
#endif

   return(temp);
}

#ifdef FALSE
Individual set_ind(int num_torsions,  GridMapSetInfo *info, State state)
{
   Genotype temp_Gtype;
   Phenotype temp_Ptype;
   int i;

   temp_Gtype = generate_Gtype(num_torsions, info);
   temp_Ptype = generate_Ptype(num_torsions, info);

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
                State now, 
                unsigned int num_evals, unsigned int pop_size, 
                int outlev, 
                unsigned int extOutputEveryNgens, Molecule *mol, 
                Boole B_template,
                Boole B_RandomTran0, Boole B_RandomQuat0, Boole B_RandomDihe0,
                GridMapSetInfo *info, char FN_pop_file[MAX_CHARS] )
{
    register unsigned int i;
    register int j;
    int num_iterations = 0, num_loops = 0, allEnergiesEqual = 1, numTries = 0;
    int indiv = 0; // Number of Individual in Population to set initial state variables for.
    double firstEnergy = 0.0;
    EvalMode localEvalMode = Normal_Eval;
    FILE *pop_fileptr;

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
#ifdef DEBUG
    // gmm 20-FEB-2003
    (void)fprintf(logFile,"call_glss.cc/State call_glss(): Creating individual i= %d in thisPop[i];  about to call random_ind()...\n", i);
#endif
            thisPop[i] = random_ind( now.ntor, info);
#ifdef DEBUG
    // gmm 20-FEB-2003
    (void)fprintf(logFile,"call_glss.cc/State call_glss(): Created  individual i= %d in thisPop[i]\n", i);
#endif
            thisPop[i].mol = mol;
            thisPop[i].age = 0L;
        }

    // If initial values were supplied, put them in thisPop[0]
    if (!B_RandomTran0) {
          if (outlev > 1) { (void)fprintf(logFile, "Setting the initial translation (tran0) for individual number %d to %.2lf %.2lf %.2lf\n\n", indiv+1, now.T.x, now.T.y, now.T.z); }
      thisPop[indiv].genotyp.write( now.T.x, 0);
      thisPop[indiv].genotyp.write( now.T.y, 1);
      thisPop[indiv].genotyp.write( now.T.z, 2);
      // Remember to keep the phenotype up-to-date
      thisPop[indiv].phenotyp = thisPop[indiv].mapping();
    };
    if (!B_RandomQuat0) {
          if (outlev > 1) { (void)fprintf(logFile, "Setting the initial quaternion (quat0) for individual number %d to %.2lf %.2lf %.2lf  %.2lf deg\n\n", indiv+1, now.Q.nx, now.Q.ny, now.Q.nz, Deg(now.Q.ang)); }
      thisPop[indiv].genotyp.write( now.Q.nx, 3);
      thisPop[indiv].genotyp.write( now.Q.ny, 4);
      thisPop[indiv].genotyp.write( now.Q.nz, 5);
      thisPop[indiv].genotyp.write( now.Q.ang, 6);
      // Remember to keep the phenotype up-to-date
      thisPop[indiv].phenotyp = thisPop[indiv].mapping();
    };
    if (!B_RandomDihe0) {
          if (outlev > 1) { (void)fprintf(logFile, "Setting the initial torsions (dihe0) for individual number %d to ", indiv+1); }
            for (j=0; j<now.ntor; j++) {
                thisPop[indiv].genotyp.write( now.tor[j], 7+j);
                if (outlev > 1) { (void)fprintf(logFile, "%.2lf ", Deg(now.tor[j])); }
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
        if (outlev > 1) { (void)fprintf( logFile, "Global-local search iteration %d\n", ++num_loops); }
        
        global_method->search(thisPop);

        if (outlev > 2) { thisPop.printPopulationAsStates( logFile, pop_size, now.ntor); }
        if (outlev > 3) { minmeanmax( logFile, thisPop, ++num_iterations ); }

        for (i=0; i<pop_size; i++) {
            local_method->search(thisPop[i]);
        }
        if (outlev > 2) { thisPop.printPopulationAsStates( logFile, pop_size, now.ntor); }

        if (strcmp (FN_pop_file, "") != 0) { // YES, do print!
            if ((pop_fileptr = ad_fopen( FN_pop_file, "w")) == NULL) {
                pr(logFile, "\n%s: ERROR:  I'm sorry, I cannot create\"%s\".\n\n", programname, FN_pop_file);
            } else {
                thisPop.printPopulationAsCoordsEnergies( pop_fileptr, pop_size, now.ntor); 
                fclose( pop_fileptr );
            }
        }

        (void)fflush(logFile);
    } while ((evaluate.evals() < num_evals) && (!global_method->terminate()));

    thisPop.msort(3);
    (void)fprintf(logFile,"Final-Value: %.3f\n", thisPop[0].value(Normal_Eval));
    return( thisPop[0].state(now.ntor) );
}
