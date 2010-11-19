/*

 $Id: eval.cc,v 1.22.2.1 2010/11/19 20:09:29 rhuey Exp $

 AutoDock 

 Copyright (C) 1989-2007,  Scott Halliday, Rik Belew, Garrett M. Morris, David S. Goodsell, Ruth Huey, Arthur J. Olson, 
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
     These are the functions associated with the evaluation object.

                                rsh 9/95
Modification: HL 11/06/2007 add state for multi-ligands 		
********************************************************************/


#include <math.h>
#include "eval.h"
#include "stateLibrary.h"
#include "assert.h"

extern FILE *logFile;
extern int nlig;    //global variable number of ligand set in main.cc

#include <stdio.h>
#include <string.h>

#ifdef sgi
    #include <ieeefp.h>
#endif

#ifdef sun
    #include <ieeefp.h>
#endif

/*  The chromosome is assumed to have a layout like this -

       | x | y | z | qx | qy | qz | qw | tor1 | ... | tor N |

 	For multiple ligands:  -Huameng
 	The chromosome is assumed to have a layout like this -

     	| x1 | y1 | z1 | qx1 | qy1 | qz1 | qw1 
     	| x2 | y2 | z2 | qx2 | qy2 | qz2 | qw2 | ... 
     	| tor1 | ... | tor N |
    where:
       x is the x translation
       y is the y translation
       z is the z translation
       qx, qy, qz, qw are the components of a 4D-normalized quaternion
       tor 1, ..., tor N are the ntor torsion angles
*/

void make_state_from_rep(Representation **rep, State *stateNow)
/*
    This routine modifies the various components of stateNow to correspond
    to the chromosome.  
*/
{
   register int i = 0;
   register int start_idx = 0;
      
#ifdef DEBUG
   (void)fprintf(logFile, "eval.cc/make_state_from_rep(Representation **rep, State *stateNow)\n");
   (void)fprintf(logFile, "global nlig=%d Total torsions=%d\n", nlig, stateNow->ntor);
#endif /* DEBUG */
  
   // modify to handle multi-ligand docking -Huameng 11/08/2007 
   stateNow->nlig = nlig;     
   for(i = 0; i < nlig; i++) {  	
   	   //each ligand uses 4 representation
	   start_idx = 4*i;  
	    	
	   //  Do the translations
	   assert(    !ISNAN( rep[start_idx + 0]->gene(0).real ) );
	   stateNow->T[i].x = rep[start_idx + 0]->gene(0).real;
	   assert(    !ISNAN( rep[start_idx + 1]->gene(0).real ) );
	   stateNow->T[i].y = rep[start_idx + 1]->gene(0).real;
	   assert(    !ISNAN( rep[start_idx + 2]->gene(0).real ) );
	   stateNow->T[i].z = rep[start_idx + 2]->gene(0).real;
	
	   //  Set up the quaternion
	   assert(    !ISNAN( rep[start_idx + 3]->gene(0).real ) );
	   stateNow->Q[i].x = rep[start_idx + 3]->gene(0).real;
	   assert(    !ISNAN( rep[start_idx + 3]->gene(1).real ) );
	   stateNow->Q[i].y = rep[start_idx + 3]->gene(1).real;
	   assert(    !ISNAN( rep[start_idx + 3]->gene(2).real ) );
	   stateNow->Q[i].z = rep[start_idx + 3]->gene(2).real;
	   assert(    !ISNAN( rep[start_idx + 3]->gene(3).real ) );
	   stateNow->Q[i].w = rep[start_idx + 3]->gene(3).real;
	
	   // Generate the corresponding axis-angle ("rotation")
	   Quat q_axis_angle;
	   q_axis_angle = convertQuatToRot( stateNow->Q[i] );
	
	   // Update the axis-angle values in stateNow
	   stateNow->Q[i].nx = q_axis_angle.nx;
	   stateNow->Q[i].ny = q_axis_angle.ny;
	   stateNow->Q[i].nz = q_axis_angle.nz;
	   stateNow->Q[i].ang = q_axis_angle.ang;
   
   } //nlig
   
   
   // Copy the angles,  4*nlig 
   for (i=0; i < stateNow->ntor; i++) {
      assert( !ISNAN( rep[4*nlig]->gene(i).real ) );
      stateNow->tor[i] = rep[4*nlig]->gene(i).real;
   }

   // mkUnitQuat(&(stateNow->Q));
}

double Eval::operator()(Representation **rep)
{
   make_state_from_rep(rep, &stateNow);
   return eval();
}

double Eval::operator()(Representation **rep, int term)
{
   make_state_from_rep(rep, &stateNow);
   return eval(term);
}


double Eval::eval()
{
   register int i;
   int   B_outside = 0;
   int   I_tor = 0;
   int   indx = 0;
   double energy = 0.0L;
   // Real emap[MAX_ATOMS] = { 0.0L };
   // Real elec[MAX_ATOMS] = { 0.0L };

#ifdef DEBUG
    (void)fprintf(logFile,"eval.cc/double Eval::eval()\n");
#endif /* DEBUG */

#ifdef DEBUG
	// multip Ligands Huameng Li 10/20/2007
	for(i = 0; i < stateNow.nlig; i++ ) {
	    if (is_out_grid_info(stateNow.T[i].x, stateNow.T[i].y, stateNow.T[i].z)) {
	       (void)fprintf(logFile,"eval.cc/stateNow.T is outside grid!\n");
	    }
	}
#endif /* DEBUG */

#ifdef DEBUG
    (void)fprintf(logFile,"eval.cc/eval()  Converting state to coordinates...\n");
    printState( logFile, stateNow, 2 );
#endif /* DEBUG */

   // Ligand could be inside or could still be outside, check all the atoms...
   // cnv_state_to_coords(stateNow, vt, tlist, stateNow.ntor, crdreo, crd, natom);
   cnv_state_to_coords(stateNow, vt, tlist, stateNow.ntor, crdpdb, crd, natom);

#ifdef DEBUG
(void)fprintf(logFile,"eval.cc/Checking to see if all coordinates are inside grid...\n");
#endif /* DEBUG */

   //  Check to see if crd is valid
   for (i=0; (i<natom)&&(!B_outside); i++) {
      B_outside = is_out_grid_info(crd[i][0], crd[i][1], crd[i][2]);
   } // i


    if (B_compute_intermol_energy) {
        energy = trilinterp( 0, natom, crd, charge, abs_charge, type, map, 
                             info, B_outside?SOME_ATOMS_OUTSIDE_GRID:ALL_ATOMS_INSIDE_GRID, 
                             ignore_inter, NULL_ELEC, NULL_EVDW, NULL_ELEC_TOTAL, NULL_EVDW_TOTAL, NULL);
    }
#ifdef DEBUG
(void)fprintf(logFile,"eval.cc/double Eval::eval() after trilinterp, energy= %.5lf\n",energy);
#endif /* DEBUG */
    energy += eintcal( nonbondlist, ptr_ad_energy_tables, crd, Nnb, B_calcIntElec, B_include_1_4_interactions, scale_1_4, qsp_abs_charge, parameterArray, B_use_non_bond_cutoff, B_have_flexible_residues ) - unbound_internal_FE;
#ifdef DEBUG
(void)fprintf(logFile,"eval.cc/double Eval::eval() after eintcal, energy= %.5lf\n",energy);
#endif /* DEBUG */
 
    if (B_isGaussTorCon) {
        for (I_tor = 0; I_tor <= stateNow.ntor; I_tor++) {
            if (B_isTorConstrained[I_tor] == 1) {
                indx = RadiansToDivs( WrpModRad(stateNow.tor[I_tor]) );
                if (B_ShowTorE) {
                    energy += (double)(US_TorE[I_tor] = US_torProfile[I_tor][indx]);
                } else {
                    energy += (double)US_torProfile[I_tor][indx];
                }
            }
        } // I_tor
    }/*if*/

   num_evals++;

   if (!finite(energy)) {
      (void)fprintf( logFile, "eval.cc:  ERROR!  energy is infinite!\n\n");
      for (i=0; i<natom; i++) {
           // (void)fprintf( logFile, "ATOM  %5d  C   INF     1    %8.3f%8.3f%8.3f %+8.2f %+6.2f  %+6.3f\n", i+1, crd[i][X], crd[i][Y], crd[i][Z], eval_emap[i], eval_elec[i], charge[i]); 
          (void)fprintf(logFile, FORMAT_PDBQ_ATOM_RESSTR, "", i+1, "C   INF     1", crd[i][X], crd[i][Y], crd[i][Z], 0.0, 0.0, charge[i]); 
          (void)fprintf(logFile, "\n");
      } // i
   }
   if (ISNAN(energy)) {
      (void)fprintf( logFile, "eval.cc:  ERROR!  energy is not a number!\n\n");
      for (i=0; i<natom; i++) {
          // (void)fprintf( logFile, "ATOM  %5d  C   NaN     1    %8.3f%8.3f%8.3f %+8.2f %+6.2f  %+6.3f\n", i+1, crd[i][X], crd[i][Y], crd[i][Z], eval_emap[i], eval_elec[i], charge[i]); 
          (void)fprintf(logFile, FORMAT_PDBQ_ATOM_RESSTR, "", i+1, "C   NaN     1", crd[i][X], crd[i][Y], crd[i][Z], 0.0, 0.0, charge[i]); 
          (void)fprintf(logFile, "\n");
      } // i
   }
#ifdef DEBUG
    (void)fprintf(logFile,"eval.cc/double Eval::eval() returns energy= %.5lf\n", energy);
#endif /*DEBUG*/

   return(energy);
   
}



/*
// Use this method, eval(int term), to compute just one particular term of the total energy
//
// we define term=0 as total energy
//           term=1 as total non-bonded energy, i.e. vdW+Hb+desolv
//           term=2 as total electrostatic energy
*/
double Eval::eval(int term)
{
   register int i;
   int   B_outside = 0;
   int   I_tor = 0;
   int   indx = 0;
   double energy = 0.0L;
   double retval = 0.0L;

	Real emap_total = 0.0L;
	Real elec_total = 0.0L;
	double desol_total = 0.0L;
	Real emap[MAX_ATOMS] = { 0.0L };
	Real elec[MAX_ATOMS] = { 0.0L };


#ifdef DEBUG
    (void)fprintf(logFile,"eval.cc/double Eval::eval(int term=%d)\n", term);
#endif /* DEBUG */

#ifdef DEBUG
	// multi-Ligand Huameng Li 10/20/2007
	for(i = 0; i < stateNow.nlig; i++) {
	    if (is_out_grid_info(stateNow.T[i].x, stateNow.T[i].y, stateNow.T[i].z)) {
	       (void)fprintf(logFile,"eval.cc/stateNow.T is outside grid!\n");
	    }
	}
#endif /* DEBUG */

#ifdef DEBUG
    (void)fprintf(logFile,"eval.cc/eval(int term)  Converting state to coordinates...\n");
    printState( logFile, stateNow, 2 );
#endif /* DEBUG */
 
   // Ligand could be inside or could still be outside, check all the atoms...
   // cnv_state_to_coords(stateNow, vt, tlist, stateNow.ntor, crdreo, crd, natom);
   cnv_state_to_coords(stateNow, vt, tlist, stateNow.ntor, crdpdb, crd, natom);

#ifdef DEBUG
(void)fprintf(logFile,"eval.cc/Checking to see if all coordinates are inside grid...\n");
#endif /* DEBUG */

   //  Check to see if crd is valid
   for (i=0; (i<natom)&&(!B_outside); i++) {
      B_outside = is_out_grid_info(crd[i][0], crd[i][1], crd[i][2]);
   }

   // Use standard energy function

#ifdef DEBUG
(void)fprintf(logFile,"eval.cc/All coordinates are inside grid...\n");
#endif /* DEBUG */

    if (B_compute_intermol_energy) {
        energy = trilinterp( 0, natom, crd, charge, abs_charge, type, map, 
                             info, B_outside?SOME_ATOMS_OUTSIDE_GRID:ALL_ATOMS_INSIDE_GRID, 
                             ignore_inter, elec, emap, &elec_total, &emap_total, &desol_total);
    }
    
#ifdef DEBUG
(void)fprintf(logFile,"eval.cc/double Eval::eval(int term=%d) after trilinterp, energy= %.5lf\n", term, energy);
#endif /* DEBUG */
    energy += eintcal( nonbondlist, ptr_ad_energy_tables, crd, Nnb, B_calcIntElec, B_include_1_4_interactions, scale_1_4, qsp_abs_charge, parameterArray, B_use_non_bond_cutoff, B_have_flexible_residues ) - unbound_internal_FE;
#ifdef DEBUG
(void)fprintf(logFile,"eval.cc/double Eval::eval(int term=%d) after eintcal, energy= %.5lf\n", term, energy);
#endif /* DEBUG */
 
    if (B_isGaussTorCon) {
        for (I_tor = 0; I_tor <= stateNow.ntor; I_tor++) {
            if (B_isTorConstrained[I_tor] == 1) {
                indx = RadiansToDivs( WrpModRad(stateNow.tor[I_tor]) );
                if (B_ShowTorE) {
                    energy += (double)(US_TorE[I_tor] = US_torProfile[I_tor][indx]);
                } else {
                    energy += (double)US_torProfile[I_tor][indx];
                }
            }
        } // I_tor
    }/*if*/

   // num_evals++;

   if (!finite(energy)) {
      (void)fprintf( logFile, "eval.cc:  ERROR!  energy is infinite!\n\n");
      for (i=0; i<natom; i++) {
           // (void)fprintf( logFile, "ATOM  %5d  C   INF     1    %8.3f%8.3f%8.3f %+8.2f %+6.2f  %+6.3f\n", i+1, crd[i][X], crd[i][Y], crd[i][Z], eval_emap[i], eval_elec[i], charge[i]); 
          (void)fprintf(logFile, FORMAT_PDBQ_ATOM_RESSTR, "", i+1, "C   INF     1", crd[i][X], crd[i][Y], crd[i][Z], 0.0, 0.0, charge[i]); 
          (void)fprintf(logFile, "\n");
      } // i
   }
   if (ISNAN(energy)) {
      (void)fprintf( logFile, "eval.cc:  ERROR!  energy is not a number!\n\n");
      for (i=0; i<natom; i++) {
          // (void)fprintf( logFile, "ATOM  %5d  C   NaN     1    %8.3f%8.3f%8.3f %+8.2f %+6.2f  %+6.3f\n", i+1, crd[i][X], crd[i][Y], crd[i][Z], eval_emap[i], eval_elec[i], charge[i]); 
          (void)fprintf(logFile, FORMAT_PDBQ_ATOM_RESSTR, "", i+1, "C   NaN     1", crd[i][X], crd[i][Y], crd[i][Z], 0.0, 0.0, charge[i]); 
          (void)fprintf(logFile, "\n");
      } // i
   }
    switch (term) {
    default:
    case 0:
        // Return the total energy.
        retval = energy;
        break;
    case 1:
        // Return the non-bonded energy, vdW+Hb+desolv.
        retval = (double)emap_total;
        break;
    case 2:
        // Return the electrostatics energy.
        retval = (double)elec_total;
        break;
    }

#ifdef DEBUG
    (void)fprintf(logFile,"eval.cc/double Eval::eval(int term=%d) returns retval= %.5lf\n", term, retval);
#endif /*DEBUG*/
   return(retval);
}

int Eval::write(FILE *out_file, Representation **rep)
{
    int i=0, retval=0;
    //char rec14[14];

#ifdef DEBUG
    (void)fprintf(logFile,"eval.cc/int Eval::write(FILE *out_file, Representation **rep)\n");
#endif /*DEBUG*/

    make_state_from_rep(rep, &stateNow);
    // cnv_state_to_coords(stateNow, vt, tlist, stateNow.ntor, crdreo, crd, natom);
       
    cnv_state_to_coords(stateNow, vt, tlist, stateNow.ntor, crdpdb, crd, natom);
    
    for (i=0; i<natom; i++) {
        // strncpy( rec14, &atomstuff[i][13], (size_t)13);
        // rec14[13]='\0';
        //strncpy(rec14, "C   RES     1\0", (size_t)14);
        //retval = fprintf( out_file, "ATOM  %5d  %13s    %8.3f%8.3f%8.3f %+8.2f %+6.2f  %+6.3f\n", i+1, rec14, crd[i][X], crd[i][Y], crd[i][Z], 0., 0., charge[i]); 
        retval = fprintf( out_file, FORMAT_PDBQ_ATOM_RESSTR, "", i+1, "C   RES     1", crd[i][X], crd[i][Y], crd[i][Z], 0., 0., charge[i]); 
        (void)fprintf(out_file, "\n");
    } // i
    return retval;
}


/********************************************************
 *  added wrapper function to compute free Energy for PSO
 *  -Huameng 02/23/08
 *
 ********************************************************/ 

double Eval::evalpso(State *stateObj)
{
   double energy = 0.0L;
     
   // need to conversion between x,y,z,w and nx,ny,nz,ang before calculate E.  
   for(int i = 0; i < nlig; i++) {   		 
   		//stateObj->Q[i] = convertQuatToRot(stateObj->Q[i]);	 
   		stateObj->Q[i] = convertRotToQuat(stateObj->Q[i]); 
   		stateObj->Q[i] = normQuat(stateObj->Q[i]);  		   	     	  	 
   }
 
   // need to set stateNow for Eval object
   // use memcpy may be faster than copyState(*dest, source)   
   //memcpy(&stateNow, stateObj, sizeof(State));
   copyState(&stateNow, *stateObj);
   
   //do Energy evalulation
   energy = eval();
   
   return(energy);
}


// we define term=0 as total energy
//           term=1 as total non-bonded energy, i.e. vdW+Hb+desolv
//           term=2 as total electrostatic energy
double Eval::evalpso(State *stateObj, int term)
{
   double energy = 0.0L;
   // need to conversion between x,y,z,w and nx,ny,nz,ang before calculate E.  
   for(int i = 0; i < nlig; i++) {   		
   		//stateObj->Q[i] = convertQuatToRot(stateObj->Q[i]);	 
   		stateObj->Q[i] = convertRotToQuat(stateObj->Q[i]);
   		stateObj->Q[i] = normQuat(stateObj->Q[i]);    		   	     	  	 
   }
   //memcpy(&stateNow, stateObj, sizeof(State));
   copyState(&stateNow, *stateObj);
   energy = eval(term);
   
   //fprintf(logFile,"eval.cc/double Eval::evalpso() returns energy= %.5lf\n", energy);
   return(energy);
}


/*
 *  code for ANT-COLINY
 * 
 */
#if defined(USING_COLINY) // {
double Eval::operator()(double* vec, int len)
{
   make_state_from_rep(vec, len, &stateNow);
   return eval();
}


void make_state_from_rep(double *rep, int len, State *now)
{
	#ifdef DEBUG
    	(void)fprintf(logFile, "USING_COLINY. eval.cc/make_state_from_rep(double *rep, int len, State *now)\n");
	#endif // DEBUG
	//  modify for multi-ligands  -Huameng 11/07/2007	
	int i;
	for(i = 0; i < nlig; i++) {	
			
	    //  Do the translations
	    now->T[i].x = rep[7*i + 0];
	    now->T[i].y = rep[7*i + 1];
	    now->T[i].z = rep[7*i + 2];
	
	    //  Set up the quaternion
	    now->Q[i].x = rep[7*i + 3];
	    now->Q[i].y = rep[7*i + 4];
	    now->Q[i].z = rep[7*i + 5];
	    now->Q[i].w = rep[7*i + 6];
	
	    now->Q[i] = convertQuatToRot( now->Q[i] );
	}
    //  Copy the angles
    now->ntor = len - 7*nlig;
    for (i=0, j=7*nlig; j < len; i++, j++) {
      now->tor[i] = rep[j];
    }

    //mkUnitQuat(&(now->Q));
}

extern Eval evaluate;

double ADEvalFn(double* x, int n)
{
	//
	// Normalize the data
	//
	//
	// Quaternion vector
	/*
	double sum=0.0;
	if (x[3] < 0.0) x[3] = 1e-16;
	if (x[4] < 0.0) x[4] = 1e-16;
	if (x[5] < 0.0) x[5] = 1e-16;
	*/
	double sum = sqrt(x[3]*x[3]+x[4]*x[4]+x[5]*x[5]);
	if (sum < 1e-8)
	   x[3]=x[4]=x[5]=1.0L/sqrt(3.0L);
	   else {
	      x[3] /= sum;
	      x[4] /= sum;
	      x[5] /= sum;
	      }
	
	// torsion angles
	for (int i=6; i<n; i++)
	  x[i] = WrpModRad(x[i]);
	
	return ::evaluate(x,n);
}
//
#endif // USING_COLINY // }
