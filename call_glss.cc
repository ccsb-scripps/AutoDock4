/*

 $Id: call_glss.cc,v 1.46 2010/01/08 20:13:47 mp Exp $

 AutoDock  

Copyright (C) 2009 The Scripps Research Institute. All rights reserved.
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
#include "qmultiply.h"

extern FILE *logFile;
extern char *programname;

int global_ntor;

Eval evaluate;

Representation **generate_R(int num_torsions, GridMapSetInfo *info)
{
   Representation **retval;
   Quat q;

#ifdef DEBUG
    (void)fprintf(logFile,"call_glss.cc/Representation **generate_R()  about to create a new Representation with 5 elements, retval...\n");
#endif
   retval = new Representation *[5];
#ifdef DEBUG
    (void)fprintf(logFile,"call_glss.cc/Representation **generate_R()  done creating   a new Representation with 5 elements, retval...\n");
#endif
   // Set the x-translation
   retval[0] = new RealVector( 1, info->lo[X], info->hi[X] );
   // Set the y-translation
   retval[1] = new RealVector( 1, info->lo[Y], info->hi[Y] );
   // Set the z-translation
   retval[2] = new RealVector( 1, info->lo[Z], info->hi[Z] );

   // Generate a uniformly-distributed random quaternion for a random rotation (UDQ)
   q = uniformQuat();
   q = convertQuatToRot( q );
#ifdef DEBUG
   printQuat( logFile, q );
#endif

   // Set the unit vector components (the "axis"), for the rotation about axis
   retval[3] = new RealVector( 3, -1., 1., q.nx, q.ny, q.nz ); // uniformly-distributed quaternion (UDQ)

   // Set the angle (the "rotation") for the rotation about axis, 
   // and any torsion angles
   retval[4] = new RealVector( num_torsions+1, -PI, PI, q.ang ); // uniformly-distributed quaternion (UDQ)
   // retval[4] = new RealVector( num_torsions+1, -PI, PI );  // rotation-about-axis angle is uniformly distributed, -PI to PI, not UDQ

#ifdef DEBUG
    (void)fprintf(logFile,"call_glss.cc/Representation **generate_R()  done assigning each of the retval[0-5] elements...\n");
#endif

   return(retval);
}

Representation **generate_R_quaternion(int num_torsions, GridMapSetInfo *info)
{
   Representation **retval;
   Quat q;

#ifdef DEBUG
    (void)fprintf(logFile,"\ncall_glss.cc/Representation **generate_R_quaternion()  about to create a new Representation with 5 elements, retval...\n");
#endif
   retval = new Representation *[5];
#ifdef DEBUG
    (void)fprintf(logFile,"call_glss.cc/Representation **generate_R_quaternion()  done creating   a new Representation with 5 elements, retval...\n");
#endif
   // Set the x-translation
   retval[0] = new RealVector( 1, info->lo[X], info->hi[X] );
   // Set the y-translation
   retval[1] = new RealVector( 1, info->lo[Y], info->hi[Y] );
   // Set the z-translation
   retval[2] = new RealVector( 1, info->lo[Z], info->hi[Z] );

   // Generate a uniformly-distributed random quaternion for a random rotation (UDQ)
   q = uniformQuat();
   q = convertQuatToRot( q );
#ifdef DEBUG
   printQuat( logFile, q );
#endif

#ifdef DEBUG_QUAT
#ifdef DEBUG_QUAT_PRINT
    pr( logFile, "DEBUG_QUAT: generate_R_quaternion()\n" );
    (void) fflush(logFile);
#endif
    //  Make sure the quaternion is suitable for 3D rotation
    assertQuatOK( q );
#endif

   // Set the quaternion (x,y,z,w) genes
   retval[3] = new RealVector( 4, -1., 1., q.x, q.y, q.z, q.w ); // uniformly-distributed quaternion (UDQ)
   // TODO retval[3] = new ConstrainedRealVector( 4, -1., 1., q.x, q.y, q.z, q.w ); // uniformly-distributed quaternion (UDQ)

   // Set the torsion angles
   retval[4] = new RealVector( num_torsions, -PI, PI );

#ifdef DEBUG
    (void)fprintf(logFile,"call_glss.cc/Representation **generate_R_quaternion()  done assigning each of the retval[0-5] elements...\n\n");
#endif

   return(retval);
}

Genotype generate_Gtype(int num_torsions, GridMapSetInfo *info)
{
#ifdef DEBUG
    // (void)fprintf(logFile,"\ncall_glss.cc/Genotype generate_Gtype() about to call Genotype temp(5, generate_R())...\n");
    (void)fprintf(logFile,"\ncall_glss.cc/Genotype generate_Gtype() about to call Genotype temp(5, generate_R_quaternion())...\n");
#endif
   // Genotype temp((unsigned int)5, generate_R(num_torsions, info));
   Genotype temp((unsigned int)5, generate_R_quaternion(num_torsions, info));
#ifdef DEBUG
   // (void)fprintf(logFile,"call_glss.cc/Genotype generate_Gtype() done calling  Genotype temp(5, generate_R())...\n\n");
   (void)fprintf(logFile,"call_glss.cc/Genotype generate_Gtype() done calling  Genotype temp(5, generate_R_quaternion())...\n\n");
#endif

   return(temp);
}

Phenotype generate_Ptype(int num_torsions, GridMapSetInfo *info) 
{
#ifdef DEBUG
    // (void)fprintf(logFile,"\ncall_glss.cc/Genotype generate_Ptype() about to call Phenotype temp(5, generate_R())...\n");
    (void)fprintf(logFile,"\ncall_glss.cc/Genotype generate_Ptype() about to call Phenotype temp(5, generate_R_quaternion())...\n");
#endif
   // Phenotype temp((unsigned int)5, generate_R(num_torsions, info));
   Phenotype temp((unsigned int)5, generate_R_quaternion(num_torsions, info));
#ifdef DEBUG
   // (void)fprintf(logFile,"call_glss.cc/Genotype generate_Ptype() done calling  Phenotype temp(5, generate_R())...\n\n");
   (void)fprintf(logFile,"call_glss.cc/Genotype generate_Ptype() done calling  Phenotype temp(5, generate_R_quaternion())...\n\n");
#endif

   return(temp);
}

Individual random_ind(int num_torsions,  GridMapSetInfo *info) 
{

#ifdef DEBUG
    (void)fprintf(logFile,"\ncall_glss.cc/Individual random_ind()  About to generate_Gtype()...\n");
#endif
   Genotype temp_Gtype = generate_Gtype(num_torsions, info);
#ifdef DEBUG
   (void)fprintf(logFile,"call_glss.cc/Individual random_ind()  About to generate_Ptype()...\n");
#endif
   Phenotype temp_Ptype = generate_Ptype(num_torsions, info); 

#ifdef DEBUG
   (void)fprintf(logFile,"call_glss.cc/Individual random_ind()  About to Individual temp(temp_Gtype, temp_Ptype)...\n");
#endif
   //shotgun wedding: does not map genotype to phenotype
   Individual temp(temp_Gtype, temp_Ptype);
   temp.mapping();

#ifdef DEBUG
    (void)fprintf(logFile,"call_glss.cc/Individual random_ind()  Done     Individual temp(temp_Gtype, temp_Ptype)...\n\n");
#endif

   return(temp);
}

#ifdef FALSE
Individual set_ind(int num_torsions,  GridMapSetInfo *info, State state)
{
   Genotype temp_Gtype;
   Phenotype temp_Ptype;
   Quat q;
   int i;

   temp_Gtype = generate_Gtype(num_torsions, info);
   temp_Ptype = generate_Ptype(num_torsions, info);

   // use the state to generate a Genotype
   temp_Gtype.write(state.T.x, 0);
   temp_Gtype.write(state.T.y, 1);
   temp_Gtype.write(state.T.z, 2);

   q = convertRotToQuat( state.Q );

#ifdef DEBUG_QUAT
#ifdef DEBUG_QUAT_PRINT
    pr( logFile, "DEBUG_QUAT: set_ind()\n" );
    (void) fflush(logFile);
#endif
    //  Make sure the quaternion is suitable for 3D rotation
    assertQuatOK( q );
#endif

   temp_Gtype.write( q.x, 3);
   temp_Gtype.write( q.y, 4);
   temp_Gtype.write( q.z, 5);
   temp_Gtype.write( q.w, 6);

   for (i=0;i<state.ntor; i++) {
       temp_Gtype.write(state.tor[i], 7+i);
   };

   Individual temp(temp_Gtype, temp_Ptype);   

   // use mapping to generate a Phenotype
   //temp.phenotyp =  temp.mapping();
   temp.mapping();

   return(temp);
}
#endif

State call_glss(Global_Search *global_method, Local_Search *local_method, 
                State sInit, 
                unsigned int num_evals, unsigned int pop_size, 
                int outlev, 
                unsigned int extOutputEveryNgens, Molecule *mol, 
                Boole B_RandomTran0, Boole B_RandomQuat0, Boole B_RandomDihe0,
                GridMapSetInfo *info, char *FN_pop_file,
                int end_of_branch[MAX_TORS])
{
    register unsigned int i;
    register int j;
    int num_generations = 0, allEnergiesEqual = 1, numTries = 0;
    int indiv = 0; // Number of Individual in Population to set initial state variables for.
    int max_numTries = 1000;
    double firstEnergy = 0.0;
    EvalMode localEvalMode = Normal_Eval;
    FILE *pop_fileptr;

    global_method->reset(extOutputEveryNgens);
    local_method->reset();
    evaluate.reset();

    (void)fprintf( logFile, "\nCreating an initial population of %u individuals.\n", pop_size);
    Population thisPop(pop_size);
    //  Pass in the end_of_branch tree for Branch Crossover Mode.
    thisPop.set_eob( end_of_branch );

    if (sInit.ntor > 0) {
        (void)fprintf( logFile, "\nAssigning a random translation, a random orientation and %d random torsions to each of the %u individuals.\n\n", sInit.ntor, pop_size);
    } else {
        (void)fprintf( logFile, "\nAssigning a random translation and a random orientation to each of the %u individuals.\n\n", pop_size);
    }
    global_ntor = sInit.ntor; //debug
#ifdef DEBUG
    (void)fprintf(logFile,"\ncall_glss.cc/State call_glss():  {\n");
#endif
    do {
        ++numTries;
        // Create a population of pop_size random individuals...
        for (i=0; i<pop_size; i++) {
#ifdef DEBUG
    (void)fprintf(logFile,"\ncall_glss.cc/State call_glss():  Creating individual thisPop[i=%d] using random_ind(%d,info)...\n", i, sInit.ntor);
#endif
            thisPop[i] = random_ind( sInit.ntor, info );
#ifdef DEBUG
    (void)fprintf(logFile,"call_glss.cc/State call_glss(): Created  individual i= %d in thisPop[i]\n\n", i);
#endif
            thisPop[i].mol = mol;
            thisPop[i].age = 0L;
        }

        // If initial values were supplied, put them in thisPop[0] and remap
        if (!B_RandomTran0) {
            if (outlev > 1) { (void)fprintf(logFile, "Setting the initial translation (tran0) for individual number %d to %.2lf %.2lf %.2lf\n\n", indiv+1, sInit.T.x, sInit.T.y, sInit.T.z); }
            thisPop[indiv].genotyp.write( sInit.T.x, 0 );
            thisPop[indiv].genotyp.write( sInit.T.y, 1 );
            thisPop[indiv].genotyp.write( sInit.T.z, 2 );
            // Remember to keep the phenotype up-to-date
            thisPop[indiv].mapping();
        }
        if (!B_RandomQuat0) {
            if (outlev > 1) { 
                (void)fprintf(logFile, "Setting the initial orientation using axis-angle values for individual number %d to %.2lf %.2lf %.2lf  %.2lf deg\n\n", indiv+1, sInit.Q.nx, sInit.Q.ny, sInit.Q.nz, RadiansToDegrees(sInit.Q.ang)); 
                (void)fprintf(logFile, "which corresponds to the quaternion (x,y,z,w) values:  %.2lf %.2lf %.2lf %.2lf\n\n", sInit.Q.x, sInit.Q.y, sInit.Q.z, sInit.Q.w); 
            }
            thisPop[indiv].genotyp.write( sInit.Q.x, 3 );
            thisPop[indiv].genotyp.write( sInit.Q.y, 4 );
            thisPop[indiv].genotyp.write( sInit.Q.z, 5 );
            thisPop[indiv].genotyp.write( sInit.Q.w, 6 );
            // Remember to keep the phenotype up-to-date
            thisPop[indiv].mapping();
        }
        if (sInit.ntor > 0) {
            if (!B_RandomDihe0) {
                if (outlev > 1) { (void)fprintf(logFile, "Setting the initial torsions (dihe0) for individual number %d to ", indiv+1); }
                for (j=0; j<sInit.ntor; j++) {
                    thisPop[indiv].genotyp.write( sInit.tor[j], 7+j );
                    if (outlev > 1) { (void)fprintf(logFile, "%.2lf ", RadiansToDegrees(sInit.tor[j])); }
                };
                if (outlev > 1) { (void)fprintf(logFile, " deg\n\n"); }
                // Remember to keep the phenotype up-to-date
                thisPop[indiv].mapping();
            }
        }

#ifdef DEBUG
    (void)fprintf(logFile,"\n\ncall_glss.cc  // ensuring there is variation in the energies\n\n");
#endif
        // Now ensure that there is some variation in the energies...
        firstEnergy = thisPop[0].value(localEvalMode);
#ifdef DEBUG
    (void)fprintf(logFile,"\n\ncall_glss.cc  // ensuring there is variation in the energies, firstEnergy=%lf\n\n", firstEnergy);
#endif
        for (i=1; i<pop_size; i++) {
#ifdef DEBUG
    (void)fprintf(logFile,"\n\ncall_glss.cc  // ensuring there is variation in the energies, i=%d, thisPop[i].value=%lf\n\n", i, thisPop[i].value(localEvalMode));
#endif
             allEnergiesEqual = allEnergiesEqual && (thisPop[i].value(localEvalMode) == firstEnergy);
        }
        if ( pop_size>1 && allEnergiesEqual) {
            (void)fprintf(logFile,"NOTE: All energies are equal in population; re-initializing. (Try Number %d)\n", numTries);
        }
        if (numTries > max_numTries) {
            (void)fprintf(logFile,"WARNING: the number of tries has exceeded the maximum number of tries permitted.\nWARNING: AutoDock will attempt continue with the currently-generated random population.\n");
            break;
        }
    } while (pop_size>1 && allEnergiesEqual);


#ifdef DEBUG
    (void)fprintf(logFile,"\ncall_glss.cc/State call_glss():  }\n");
    if (outlev > 2) { 
    thisPop.printPopulationAsCoordsEnergies( logFile, pop_size, sInit.ntor );
    }
#endif

    if (outlev > 2) { 
        (void)fprintf( logFile, "The initial population consists of the following %d individuals:\n\n", pop_size);
        (void)fprintf( logFile, "<generation t=\"%d\" after_performing=\"initialisation of population\">\n", num_generations);
        (void)fprintf( logFile, "</generation>\n\n\n");
    }

    if (pop_size > 1 && outlev > 3 ) { minmeanmax( logFile, thisPop, num_generations, info ); }

//We now have a mapped and evaluated population suitable for global search

    (void)fprintf( logFile, "Beginning Lamarckian Genetic Algorithm (LGA), with a maximum of %u\nenergy evaluations.\n\n", num_evals);

 // M Pique 28 Oct 2009 - adding a (harmless...) thisPop.msort(1) here
 // broke the unbound-extended test. Investigating why
 // Conclusion: the test is sensitive to population order since it does
 // only one generation of glss.  So I'm removing the msort(1) and doing
 // the search for lowest-energy individual by hand.
#ifdef DEBUGMSORT
     double *beforesort = new double[pop_size];
     double *aftersort = new double[pop_size];
     int *beforeorder = new int[pop_size];
     int *afterorder = new int[pop_size];

    (void)fprintf(logFile,"@@ #evals before msort %-6u  DEBUGMSORT call_glss.cc\n", evaluate.evals());
     for (i=0; i<pop_size; i++) beforesort[i] = thisPop[i].value(Normal_Eval);
     //thisPop.msort(1);
     for (i=0; i<pop_size; i++) aftersort[i] = thisPop[i].value(Normal_Eval);

     // set beforeorder[i] to location j of each aftersort value
     for (i=0; i<pop_size; i++) beforeorder[i] = -1 ; // flag unset values
     for (i=0; i<pop_size; i++) \
     for (j=0;j<pop_size && aftersort[i]!=beforesort[beforeorder[i]=j];j++) ;
     // set afterorder[i] to location j of each beforesort value
     for (i=0; i<pop_size; i++) afterorder[i] = -1 ; // flag unset values
     for (i=0; i<pop_size; i++) \
     for (j=0;j<pop_size && beforesort[i]!=aftersort[afterorder[i]=j];j++) ;
     // report 
    (void)fprintf(logFile,"@@ #evals after msort %-6u  DEBUGMSORT call_glss.cc\n", evaluate.evals());
    (void)fprintf(logFile,"@@ pop before     after sort   locn   DEBUGMSORT call_glss.cc\n");
     for (i=0; i<pop_size; i++) 
       (void)fprintf(logFile,"%3d %9.2f %3d   %9.2f %3d\n", 
         i, beforesort[i], afterorder[i], aftersort[i], beforeorder[i]);

      delete[] beforesort;
      delete[] aftersort;
      delete[] beforeorder;
      delete[] afterorder;
#endif
    if(pop_size>0) {
       double bestenergy = thisPop[0].value(Normal_Eval);
       for (i=1; i<pop_size; i++) if(bestenergy>thisPop[i].value(Normal_Eval)) \
         bestenergy=thisPop[i].value(Normal_Eval);
       (void)fprintf(logFile,"Initial-Value: %.3f\n", bestenergy);

       // print "Population at Generation:" line with low/high/mean/median/stddev...
       // and search counts (expected to be zero)
       if(outlev>0) {
         (void) thisPop.printPopulationStatisticsVerbose(logFile, 
         num_generations, evaluate.evals(), "");
	 fprintf(logFile, " cg_count: %u", global_method->cg_count);
	 fprintf(logFile, " ci_count: %u", global_method->ci_count);
	 fprintf(logFile, " mg_count: %u", global_method->mg_count);
	 fprintf(logFile, " mi_count: %u", global_method->mi_count);
	 fprintf(logFile, " ls_count: %u", local_method->count);
	 fprintf(logFile, "\n");
	 }
     }

    do {
        ++num_generations;

        if (outlev > 1) { (void)fprintf( logFile, "Global-Local Search Iteration: %d\n", num_generations); }
        
        if (outlev > 1) { (void)fprintf( logFile, "Performing Global Search.\n"); }

        global_method->search(thisPop);

        if (outlev > 2) {
            (void)fprintf( logFile, "<generation t=\"%d\" after_performing=\"global search\">\n", num_generations);
            thisPop.printPopulationAsStates( logFile, pop_size, sInit.ntor );
            (void)fprintf( logFile, "</generation>\n\n\n");
        }

        if (pop_size > 1 && outlev > 3) { minmeanmax( logFile, thisPop, num_generations, info ); }

        if (outlev > 1) { (void)fprintf( logFile, "Performing Local Search.\n"); }

        for (i=0; i<pop_size; i++) {
            if (outlev > 1) {
                (void)fprintf( logFile, "LS: %d",num_generations); 
                (void)fprintf( logFile, " %d",i+1); 
                (void)fprintf( logFile, " %f",thisPop[i].value(localEvalMode)); 
            }
            local_method->search(thisPop[i]);
            if (outlev > 1) {
                (void)fprintf( logFile, " %f",thisPop[i].value(localEvalMode)); 
                (void)fprintf( logFile, " \n"); 
            }
        }

        if (outlev > 2) {
            (void)fprintf( logFile, "<generation t=\"%d\" after_performing=\"local search\">\n", num_generations);
            thisPop.printPopulationAsStates( logFile, pop_size, sInit.ntor );
            (void)fprintf( logFile, "</generation>\n\n\n");
        }

       // Print extended generational statistics 
       //  every generation 1 to 20
       //  every tenth generation 20 to extOutputEveryNgens (100)
       //  every "extOutputEveryNgens" (100) if greater than that
       // This is purely for studying population convergence -  M Pique 2010
       if (outlev > 1 && extOutputEveryNgens != 0 && (
	 (num_generations <= 20) ||
	 (num_generations <= extOutputEveryNgens && num_generations%10 == 0) ||
         (num_generations%extOutputEveryNgens == 0 )
	 )) {
	       // print "Population at Generation:" line with low/high/mean/median/stddev...
	       // followed by global search stats (crossover count, mutation count)
	       // and local search stats (invocation count)
	       (void) thisPop.printPopulationStatisticsVerbose(logFile, 
		 num_generations, evaluate.evals(), "");
		 fprintf(logFile, " cg_count: %u", global_method->cg_count);
		 fprintf(logFile, " ci_count: %u", global_method->ci_count);
		 fprintf(logFile, " mg_count: %u", global_method->mg_count);
		 fprintf(logFile, " mi_count: %u", global_method->mi_count);
		 fprintf(logFile, " ls_count: %u", local_method->count);
		 fprintf(logFile, "\n");
		 }
        if (pop_size > 1 && outlev > 3) { minmeanmax( logFile, thisPop, num_generations, info ); }

        if (strcmp (FN_pop_file, "") != 0) { // YES, do print!
            if ((pop_fileptr = ad_fopen( FN_pop_file, "w")) == NULL) {
                pr(logFile, "\n%s: ERROR:  I'm sorry, I cannot create\"%s\".\n\n", programname, FN_pop_file);
            } else {
                thisPop.printPopulationAsCoordsEnergies( pop_fileptr, pop_size, sInit.ntor); 
                fclose( pop_fileptr );
            }
        }

        (void)fflush(logFile);
    } while ((evaluate.evals() < num_evals) && (!global_method->terminate()));

    thisPop.msort(1);
    (void)fprintf(logFile,"Final-Value: %.3f\n", thisPop[0].value(Normal_Eval));

    // print "Population at Generation:" line with low/high/mean/median/stddev...
    if(outlev>0) {
       (void) thisPop.printPopulationStatisticsVerbose(logFile, 
        num_generations, evaluate.evals(), "");
	 fprintf(logFile, " cg_count: %u", global_method->cg_count);
	 fprintf(logFile, " ci_count: %u", global_method->ci_count);
	 fprintf(logFile, " mg_count: %u", global_method->mg_count);
	 fprintf(logFile, " mi_count: %u", global_method->mi_count);
	 fprintf(logFile, " ls_count: %u", local_method->count);
	 fprintf(logFile, "\n");
	 }

    return( thisPop[0].state(sInit.ntor) );
}
