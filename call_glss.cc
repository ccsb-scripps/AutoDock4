/*

 $Id: call_glss.cc,v 1.29.2.1 2010/11/19 20:09:29 rhuey Exp $

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
extern int nlig;   // number of ligands value assigned in main.cc
extern int ntor_lig[MAX_LIGANDS];  //no of torsions in each ligand assigned in readPDBQT.cc
extern int gene_index_lig[MAX_LIGANDS][2];  //gene num start_point & end_point of a ligand.
extern int global_ntor; // defined and assigned in main.cc

Eval evaluate;
/*
Representation **generate_R(int num_torsions, GridMapSetInfo *info)
{
   Representation **retval;
   Quat q;
   int num_rep = 4 * nlig + 1;
   int i =0;
     
#ifdef DEBUG
    (void)fprintf(logFile,"call_glss.cc/Representation **generate_R()  about to create a new Representation with 5 elements, retval...\n");
#endif
   // Modify to handle multi-ligand docking -Huameng 11/08/2007
   //retval = new Representation *[5];  
   retval = new Representation *[num_rep];
#ifdef DEBUG
    (void)fprintf(logFile,"call_glss.cc/Representation **generate_R()  done creating   a new Representation with 5 elements, retval...\n");
#endif
   // Modify to handle multi-ligand docking -Huameng 11/08/2007
   for(i = 0; i < nlig; i++) {
   	  
	   // Set the x-translation
	   retval[4*i + 0] = new RealVector( 1, info->lo[X], info->hi[X] );
	   // Set the y-translation
	   retval[4*i + 1] = new RealVector( 1, info->lo[Y], info->hi[Y] );
	   // Set the z-translation
	   retval[4*i + 2] = new RealVector( 1, info->lo[Z], info->hi[Z] );
	
	   // Generate a uniformly-distributed random quaternion for a random rotation (UDQ)
	   q = uniformQuat();
	   q = convertQuatToRot( q );
	#ifdef DEBUG
	   printQuat( logFile, q );
	#endif
	
	   // Set the unit vector components (the "axis"), for the rotation about axis
	   retval[4*i + 3] = new RealVector( 3, -1., 1., q.nx, q.ny, q.nz ); // uniformly-distributed quaternion (UDQ)
   } //nlig
   
   // Set the angle (the "rotation") for the rotation about axis, 
   // and any torsion angles
   retval[4*nlig] = new RealVector( num_torsions+1, -PI, PI, q.ang ); // uniformly-distributed quaternion (UDQ)
   // retval[4] = new RealVector( num_torsions+1, -PI, PI );  // rotation-about-axis angle is uniformly distributed, -PI to PI, not UDQ

#ifdef DEBUG
    (void)fprintf(logFile,"call_glss.cc/Representation **generate_R()  done assigning each of the retval[0-5] elements...\n");
#endif

   return(retval);
}
*/

Representation **generate_R_quaternion(int num_torsions, GridMapSetInfo *info)
{
   Representation **retval;
   Quat q;
   register int num_reps = 4*nlig + 1; // each ligand uses 4 vector representations: x,y,z and quaternion
   									   // last one is for torsions
   register int i = 0;
#ifdef DEBUG
    (void)fprintf(logFile,"\ncall_glss.cc/Representation **generate_R_quaternion()  about to create a new Representation with 5 elements, retval...\n");
#endif
   // Modify to handle multi-ligand docking -Huameng 11/08/2007
   //retval = new Representation *[5];  
   retval = new Representation *[num_reps];
#ifdef DEBUG
    (void)fprintf(logFile,"call_glss.cc/Representation **generate_R_quaternion()  done creating   a new Representation with 5 elements, retval...\n");
#endif
  // Modify to handle multi-ligand docking -Huameng 11/08/2007
   for( i = 0; i < nlig; i++) {
	   // Set the x-translation
	   retval[4*i + 0] = new RealVector( 1, info->lo[X], info->hi[X] );
	   // Set the y-translation
	   retval[4*i + 1] = new RealVector( 1, info->lo[Y], info->hi[Y] );
	   // Set the z-translation
	   retval[4*i + 2] = new RealVector( 1, info->lo[Z], info->hi[Z] );
	
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
	   retval[4*i + 3] = new RealVector( 4, -1., 1., q.x, q.y, q.z, q.w ); // uniformly-distributed quaternion (UDQ)
	   // TODO retval[3] = new ConstrainedRealVector( 4, -1., 1., q.x, q.y, q.z, q.w ); // uniformly-distributed quaternion (UDQ)
   } //end nlig
   
   // Set the torsion angles
   // Note last rep is for torsion, index = 4*nlig
   //retval[4] = new RealVector( num_torsions, -PI, PI );
   retval[4*nlig] = new RealVector( num_torsions, -PI, PI );

#ifdef DEBUG
    (void)fprintf(logFile,"call_glss.cc/Representation **generate_R_quaternion()  done assigning each of the retval[0-%d] elements...\n\n", num_reps);
#endif

   return(retval);
}

Genotype generate_Gtype(int num_torsions, GridMapSetInfo *info)
{
	//Modify to handle multi-ligand docking -Huameng 11/08/2007
	int num_reps = 4*nlig + 1; // each ligand uses 4 vector representations: x,y,z and quaternion
#ifdef DEBUG
    // (void)fprintf(logFile,"\ncall_glss.cc/Genotype generate_Gtype() about to call Genotype temp(5, generate_R())...\n");
    (void)fprintf(logFile,"\ncall_glss.cc/Genotype generate_Gtype() about to call Genotype temp(5, generate_R_quaternion())...\n");
#endif
   // Genotype temp((unsigned int)5, generate_R(num_torsions, info));
   //Genotype temp((unsigned int)5, generate_R_quaternion(num_torsions, info));
   Genotype temp((unsigned int)num_reps, generate_R_quaternion(num_torsions, info));
#ifdef DEBUG
   // (void)fprintf(logFile,"call_glss.cc/Genotype generate_Gtype() done calling  Genotype temp(5, generate_R())...\n\n");
   (void)fprintf(logFile,"call_glss.cc/Genotype generate_Gtype() done calling  Genotype temp(%d, generate_R_quaternion())...\n\n", num_reps);
#endif

   return(temp);
}

Phenotype generate_Ptype(int num_torsions, GridMapSetInfo *info) 
{
	int num_reps = 4*nlig + 1; // each ligand uses 4 vector representations: x,y,z and quaternion
#ifdef DEBUG
    // (void)fprintf(logFile,"\ncall_glss.cc/Genotype generate_Ptype() about to call Phenotype temp(5, generate_R())...\n");
    (void)fprintf(logFile,"\ncall_glss.cc/Genotype generate_Ptype() about to call Phenotype temp(5, generate_R_quaternion())...\n");
#endif
   // Phenotype temp((unsigned int)5, generate_R(num_torsions, info));
   //Phenotype temp((unsigned int)5, generate_R_quaternion(num_torsions, info));
   Phenotype temp((unsigned int)num_reps, generate_R_quaternion(num_torsions, info));
#ifdef DEBUG
   // (void)fprintf(logFile,"call_glss.cc/Genotype generate_Ptype() done calling  Phenotype temp(5, generate_R())...\n\n");
   (void)fprintf(logFile,"call_glss.cc/Genotype generate_Ptype() done calling  Phenotype temp(%d, generate_R_quaternion())...\n\n", num_reps);
   fflush(logFile);
#endif

   return(temp);
}

// set random individual
Individual random_ind(int num_torsions,  GridMapSetInfo *info) 
{
   Genotype temp_Gtype;
   Phenotype temp_Ptype;

#ifdef DEBUG
    (void)fprintf(logFile,"\ncall_glss.cc/Individual random_ind()  About to generate_Gtype()...\n");
#endif
   temp_Gtype = generate_Gtype(num_torsions, info);
#ifdef DEBUG
    (void)fprintf(logFile,"call_glss.cc/Individual random_ind()  About to generate_Ptype()...\n");
#endif
   temp_Ptype = generate_Ptype(num_torsions, info); // differs on Linux gcc 2.96 and Mac OS X gcc 3.1

#ifdef DEBUG
    (void)fprintf(logFile,"call_glss.cc/Individual random_ind()  About to Individual temp(temp_Gtype, temp_Ptype)...\n");
#endif
   Individual temp(temp_Gtype, temp_Ptype);
#ifdef DEBUG
    (void)fprintf(logFile,"call_glss.cc/Individual random_ind()  Done     Individual temp(temp_Gtype, temp_Ptype)...\n\n");
	fflush(logFile);
#endif

   return(temp);
}


#ifdef FALSE
// set individual
Individual set_ind(int num_torsions,  GridMapSetInfo *info, State state)
{
   Genotype temp_Gtype;
   Phenotype temp_Ptype;
   Quat q;
   int i;

   temp_Gtype = generate_Gtype(num_torsions, info);
   temp_Ptype = generate_Ptype(num_torsions, info);
   
	// handle multi-ligands  -Huameng Li 11/08/2007
	for(int n = 0; n < nlig; n++) {
		// use the state to generate a Genotype
	   temp_Gtype.write(state.T[n].x, 7*n + 0);
	   temp_Gtype.write(state.T[n].y, 7*n + 1);
	   temp_Gtype.write(state.T[n].z, 7*n + 2);
	
	   q = convertRotToQuat( state.Q[n] );
	
	#ifdef DEBUG_QUAT
	#ifdef DEBUG_QUAT_PRINT
	    pr( logFile, "DEBUG_QUAT: set_ind()\n" );
	    (void) fflush(logFile);
	#endif
	    //  Make sure the quaternion is suitable for 3D rotation
	    assertQuatOK( q );
	#endif
	
	   temp_Gtype.write( q.x, 7*n + 3);
	   temp_Gtype.write( q.y, 7*n + 4);
	   temp_Gtype.write( q.z, 7*n + 5);
	   temp_Gtype.write( q.w, 7*n + 6);
	} // nlig
   
   for (i=0;i<state.ntor; i++) {
       temp_Gtype.write(state.tor[i], 7*nlig + i);
   };

   Individual temp(temp_Gtype, temp_Ptype);   

   // use mapping to generate a Phenotype
   temp.phenotyp =  temp.mapping();

   return(temp);
}
#endif

State call_glss(Global_Search *global_method, Local_Search *local_method, 
                State sInit, 
                unsigned int num_evals, unsigned int pop_size, 
                int outlev, 
                unsigned int extOutputEveryNgens, Molecule *mol, 
                Boole B_RandomTran0, Boole B_RandomQuat0, Boole B_RandomDihe0,
                GridMapSetInfo *info, char FN_pop_file[MAX_CHARS] )
{
    register unsigned int i;
    register int j;
    int num_generations = 0, allEnergiesEqual = 1, numTries = 0;
    int indiv = 0; // Number of Individual in Population to set initial state variables for.
    double firstEnergy = 0.0;
    double indvEnergy = 0.0;
    EvalMode localEvalMode = Normal_Eval;
    FILE *pop_fileptr;

    global_method->reset(extOutputEveryNgens);
    local_method->reset();
    evaluate.reset();
      
    (void)fprintf( logFile, "\nCreating an initial population of %u individuals.\n", pop_size);
    Population thisPop(pop_size);

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
			///////////////////////////////////////////////			
            // Here call generate_Gtype(), generate_Ptype()
            // --> generate_R(), generate_R_quaternion() 
            // --> create STATE for each individual
            ///////////////////////////////////////////////
            thisPop[i] = random_ind( sInit.ntor, info );   
                   
#ifdef DEBUG
    		(void)fprintf(logFile,"call_glss.cc/State call_glss(): Created  individual i= %d in thisPop[i]\n\n", i);    
#endif
            // assign ligands to this individual
            thisPop[i].mol = mol;
            thisPop[i].age = 0L;
            // Make sure the phenotype corresponds to the genotype.
            /// gmm 2006-10-18 thisPop[i].phenotyp = thisPop[i].mapping();
        } // for pop_size
        
	    fprintf(logFile, "call_glss.cc/State call_glss(): After random_ind(), created %d individuals.\n", i);
	   	  
        // If initial values were supplied, put them in thisPop[0]
        // Note here indiv = 0, replace [indiv] with [0] so it is not confusing
        // Handle multi-ligand docking  -Huameng Li 11/08/2007
        for( i = 0; i < (unsigned int) nlig; i++) {
        	if (!B_RandomTran0) {
	            if (outlev > 1) { (void)fprintf(logFile, "Setting the initial translation (tran0) for individual number %d to %.2lf %.2lf %.2lf\n\n", indiv+1, sInit.T[i].x, sInit.T[i].y, sInit.T[i].z); }
	            thisPop[0].genotyp.write( sInit.T[i].x, 7*i + 0 );
	            thisPop[0].genotyp.write( sInit.T[i].y, 7*i + 1 );
	            thisPop[0].genotyp.write( sInit.T[i].z, 7*i + 2 );
	            // Remember to keep the phenotype up-to-date
	            thisPop[0].phenotyp = thisPop[0].mapping();
	        }
	        if (!B_RandomQuat0) {	        	
	            if (outlev > 1) { 
	                (void)fprintf(logFile, "Setting the initial orientation using axis-angle values for individual number %d to %.2lf %.2lf %.2lf  %.2lf deg\n\n", indiv+1, sInit.Q[i].nx, sInit.Q[i].ny, sInit.Q[i].nz, RadiansToDegrees(sInit.Q[i].ang)); 
	                (void)fprintf(logFile, "which corresponds to the quaternion (x,y,z,w) values:  %.2lf %.2lf %.2lf %.2lf\n\n", sInit.Q[i].x, sInit.Q[i].y, sInit.Q[i].z, sInit.Q[i].w); 
	            }           
	            thisPop[0].genotyp.write( sInit.Q[i].x, 7*i + 3 );
	            thisPop[0].genotyp.write( sInit.Q[i].y, 7*i + 4 );
	            thisPop[0].genotyp.write( sInit.Q[i].z, 7*i + 5 );
	            thisPop[0].genotyp.write( sInit.Q[i].w, 7*i + 6 );
	            // Remember to keep the phenotype up-to-date
	            thisPop[0].phenotyp = thisPop[0].mapping();
	        }	        
        } // nlig
        
        if (sInit.ntor > 0) {
            if (!B_RandomDihe0) {
                if (outlev > 1) { (void)fprintf(logFile, "Setting the initial torsions (dihe0) for individual number %d to ", indiv+1); }
                for (j=0; j<sInit.ntor; j++) {
                    thisPop[0].genotyp.write( sInit.tor[j], 7*nlig + j );
                    if (outlev > 1) { (void)fprintf(logFile, "%.2lf ", RadiansToDegrees(sInit.tor[j])); }
                };
                if (outlev > 1) { (void)fprintf(logFile, " deg\n\n"); }
                // Remember to keep the phenotype up-to-date
                thisPop[0].phenotyp = thisPop[0].mapping();
            }
        }
       
#ifdef DEBUG  
		pr(logFile, "Torsion in sInit.ntor for population is %d\n", sInit.ntor);     
        pr(logFile, "Before firstEnergy, thisPop[0].genotyp.number_of_genes=%d\tphenotyp.num_dimensions=%d\tphenotyp.number_of_points= %d\n", 
        	thisPop[0].genotyp.num_genes(), 
        	thisPop[0].phenotyp.num_dimensions(),       	
        	thisPop[0].phenotyp.num_pts());
#endif       	
        // Now ensure that there is some variation in the energies...
        firstEnergy = thisPop[0].value(localEvalMode);
           
#ifdef DEBUG
    (void)fprintf(logFile,"\n\ncall_glss.cc  // ensuring there is variation in the energies, firstEnergy=%lf\n\n", firstEnergy);	
#endif
        for (i=1; i<pop_size; i++) {
        	
			 indvEnergy = thisPop[i].value(localEvalMode);
			 
             //allEnergiesEqual = allEnergiesEqual && (thisPop[i].value(localEvalMode) == firstEnergy);
             allEnergiesEqual = allEnergiesEqual && (indvEnergy == firstEnergy);
#ifdef DEBUG
    //(void)fprintf(logFile,"\n\ncall_glss.cc Initial individual i = %d, ENERGY (thisPop[i].value)=%lf\n\n", i, indvEnergy );	
#endif                          
        }
        if (allEnergiesEqual) {
            (void)fprintf(logFile,"NOTE: All energies are equal in population; re-initializing. (Try Number %d)\n", numTries);
        }
    } while (allEnergiesEqual); 
    
    //////////////////////////////////////////////////////////////////// 
    // End creating individuals in this population and engergy checking
	////////////////////////////////////////////////////////////////////
	(void)fprintf(logFile,"\ncall_glss.cc/State call_glss():  Done initial population.\nfirstEnergy = %.3f\n", firstEnergy);
	
    if (outlev > 2) { 
        (void)fprintf( logFile, "The initial population consists of the following %d individuals:\n\n", pop_size);
        (void)fprintf( logFile, "<generation t=\"%d\" after_performing=\"initialisation of population\">\n", num_generations);
        thisPop.printPopulationAsStates( logFile, pop_size, sInit.ntor );
        (void)fprintf( logFile, "</generation>\n\n\n");
    }

    if (outlev > 3) { minmeanmax( logFile, thisPop, num_generations, info ); }
   
    
	///////////////////////////////////////////////
	// Beginning Lamarckian Genetic Algorithm (LGA)
	///////////////////////////////////////////////
    (void)fprintf( logFile, "Beginning Lamarckian Genetic Algorithm (LGA), with a maximum of %u energy evaluations.\n\n", num_evals);
	fflush(logFile);
	
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

        if (outlev > 3) { minmeanmax( logFile, thisPop, num_generations, info ); }

        if (outlev > 1) { (void)fprintf( logFile, "Performing Local Search.\n"); }

        for (i=0; i<pop_size; i++) {
            local_method->search(thisPop[i]);
        }
        
		
        if (outlev > 2) {
            (void)fprintf( logFile, "<generation t=\"%d\" after_performing=\"local search\">\n", num_generations);
            thisPop.printPopulationAsStates( logFile, pop_size, sInit.ntor );
            (void)fprintf( logFile, "</generation>\n\n\n");
        }
		
        if (outlev > 3) { minmeanmax( logFile, thisPop, num_generations, info ); }

        if (strcmp (FN_pop_file, "") != 0) { // YES, do print!
            if ((pop_fileptr = ad_fopen( FN_pop_file, "w")) == NULL) {
                pr(logFile, "\n%s: ERROR:  I'm sorry, I cannot create\"%s\".\n\n", programname, FN_pop_file);
            } else {
                thisPop.printPopulationAsCoordsEnergies( pop_fileptr, pop_size, sInit.ntor); 
                fclose( pop_fileptr );
            }
        }		        
    } while (( (unsigned) evaluate.evals() < num_evals) && (!global_method->terminate()));

    thisPop.msort(3);
    (void)fprintf(logFile,"Final-Value: %.3f\n", thisPop[0].value(Normal_Eval));
    
    //State finalState = thisPop[0].state(sInit.ntor);   
    //pr (logFile, "in call_gals:, finalState E = %.3lf\n", evaluate.evalpso(&finalState));
    //return( finalState );
    return( thisPop[0].state(sInit.ntor) );
}
