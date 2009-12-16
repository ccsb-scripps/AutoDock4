/*

 $Id: eval.cc,v 1.28 2009/12/16 23:35:54 rhuey Exp $

 AutoDock  

Copyright (C) 2009 The Scripps Research Institute. All rights reserved.

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
********************************************************************/


#include <math.h>
#include "eval.h"
#include "stateLibrary.h"
#include "assert.h"

extern FILE *logFile;

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
   register int i;

#ifdef DEBUG
   (void)fprintf(logFile, "eval.cc/make_state_from_rep(Representation **rep, State *stateNow)\n");
#endif /* DEBUG */

   //  Do the translations
   assert( !ISNAN( rep[0]->gene(0).real ) );
   stateNow->T.x = rep[0]->gene(0).real;
   assert( !ISNAN( rep[1]->gene(0).real ) );
   stateNow->T.y = rep[1]->gene(0).real;
   assert( !ISNAN( rep[2]->gene(0).real ) );
   stateNow->T.z = rep[2]->gene(0).real;

   //  Set up the quaternion
   assert( !ISNAN( rep[3]->gene(0).real ) );
   stateNow->Q.x = rep[3]->gene(0).real;
   assert( !ISNAN( rep[3]->gene(1).real ) );
   stateNow->Q.y = rep[3]->gene(1).real;
   assert( !ISNAN( rep[3]->gene(2).real ) );
   stateNow->Q.z = rep[3]->gene(2).real;
   assert( !ISNAN( rep[3]->gene(3).real ) );
   stateNow->Q.w = rep[3]->gene(3).real;

   // Generate the corresponding axis-angle ("rotation")
   Quat q_axis_angle;
   q_axis_angle = convertQuatToRot( stateNow->Q );

   // Update the axis-angle values in stateNow
   stateNow->Q.nx = q_axis_angle.nx;
   stateNow->Q.ny = q_axis_angle.ny;
   stateNow->Q.nz = q_axis_angle.nz;
   stateNow->Q.ang = q_axis_angle.ang;
   
   //  Copy the angles
   for (i=0; i<stateNow->ntor; i++) {
      assert( !ISNAN( rep[4]->gene(i).real ) );
      stateNow->tor[i] = rep[4]->gene(i).real;
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
#ifdef DEBUG
   (void) fprintf(logFile,"eval.cc eval() calling eval(3)\n");
#endif /* DEBUG */
   return eval(3); // default is total energy
}


double Eval::eval(int term)

// Use this method, eval(int term), to compute just one particular term of the total energy
//
// we define term=0 as total energy
//           term=1 as total non-bonded energy, i.e. vdW+Hb+desolv
//           term=2 as total electrostatic energy
//           term=3 as total energy if invoked by eval()

{
   register int i;
   int   B_outside = 0;
   int   I_tor = 0;
   int   indx = 0;
   double energy = 0.0L;
   double retval = 0.0L;

	Real emap_total = 0.0L;
	Real elec_total = 0.0L;
	Real emap[MAX_ATOMS] = { 0.0L };
	Real elec[MAX_ATOMS] = { 0.0L };


#ifdef DEBUG
    (void)fprintf(logFile,"eval.cc/double Eval::eval(int term=%d)\n", term);
#endif /* DEBUG */

#ifdef DEBUG
    if (is_out_grid_info(stateNow.T.x, stateNow.T.y, stateNow.T.z)) {
       (void)fprintf(logFile,"eval.cc/stateNow.T is outside grid!\n");
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
    if(B_outside) (void)fprintf(logFile,"eval.cc/Some coordinates are outside grid...\n");
    else (void)fprintf(logFile,"eval.cc/All coordinates are inside grid...\n");
#endif /* DEBUG */

    if (B_compute_intermol_energy) {
        if(term==3) // do not need energy breakdown in this eval() case
        energy = trilinterp( 0, natom, crd, charge, abs_charge, type, map, 
                             info, B_outside?SOME_ATOMS_OUTSIDE_GRID:ALL_ATOMS_INSIDE_GRID, 
                             ignore_inter, NULL_ELEC, NULL_EVDW, NULL_ELEC_TOTAL, NULL_EVDW_TOTAL);
        else
        energy = trilinterp( 0, natom, crd, charge, abs_charge, type, map, 
                             info, B_outside?SOME_ATOMS_OUTSIDE_GRID:ALL_ATOMS_INSIDE_GRID, 
                             ignore_inter, elec, emap, &elec_total, &emap_total);
    }
    
#ifdef DEBUG
(void)fprintf(logFile,"eval.cc/double Eval::eval(int term=%d) after trilinterp, energy= %.5lf\n", term, energy);
#endif /* DEBUG */
    energy += eintcal( nonbondlist, ptr_ad_energy_tables, crd, Nnb,
                 B_calcIntElec, B_include_1_4_interactions, 
                 scale_1_4, qsp_abs_charge,
                 B_use_non_bond_cutoff, B_have_flexible_residues);
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

   // increment evaluation counter only for "total energy" case
   if(term==3) num_evals++;

   if ((!finite(energy)) || ISNAN(energy)) {
      (void)fprintf( logFile, "eval.cc:  ERROR!  energy is %s!\n\n",
       (!finite(energy))?"infinite":"not a number");
      for (i=0; i<natom; i++) {
            print_PDBQ_atom_resstr( logFile, "", i,   " C   INF     1 ", crd, 
             0.0, 0.0, charge[i],"\n"); 
      } // i
   }
    switch (term) {
    default:
    case 0:
    case 3:
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
            print_PDBQ_atom_resstr( logFile, "", i,   " C   RES     1 ", crd, 
             0.0, 0.0, charge[i],"\n"); 
    } // i
    return retval;
}

#if defined(USING_COLINY) // {
double Eval::operator()(double* vec, int len)
{
   make_state_from_rep(vec, len, &stateNow);
   return eval();
}


void make_state_from_rep(double *rep, int n, State *now)
{
#   ifdef DEBUG
    (void)fprintf(logFile, "eval.cc/make_state_from_rep(double *rep, int n, State *now)\n");
#   endif /* DEBUG */

    //  Do the translations
    now->T.x = rep[0];
    now->T.y = rep[1];
    now->T.z = rep[2];

    //  Set up the quaternion
    now->Q.x = rep[3];
    now->Q.y = rep[4];
    now->Q.z = rep[5];
    now->Q.w = rep[6];

    now->Q = convertQuatToRot( now->Q );

    //  Copy the angles
    now->ntor = n - 7;
    for (int i=0, j=7; j<n; i++, j++) {
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

double Eval::evalpso(State *state)
{
    register int i;
    int   B_outside = 0;
    int   I_tor = 0;
    int   indx = 0;
    double energy = 0.0;

	//double einter = 0.0; 
	//double eintra = 0.0; 
    Real emap_total = 0.0L;
    Real elec_total = 0.0L;
    Real emap[MAX_ATOMS] = { 0.0L };
    Real elec[MAX_ATOMS] = { 0.0L };

	
	mkUnitQuat( &(state->Q) );
	copyState(&stateNow, *state);

#ifdef DEBUG
    if (is_out_grid_info(stateNow.T.x, stateNow.T.y, stateNow.T.z)) {
       (void)fprintf(logFile,"eval.cc/evalpso: stateNow.T is outside grid!\n");
    }
    (void)fprintf(logFile,"eval.cc/evalpso: Converting state to coordinates...\n");
#endif /* DEBUG */
 
   // Ligand could be inside or could still be outside, check all the atoms...
   cnv_state_to_coords(stateNow, vt, tlist, stateNow.ntor, crdpdb, crd, natom);

#ifdef DEBUG
(void)fprintf(logFile,"eval.cc/Checking to see if all coordinates are inside grid...\n");
#endif /* DEBUG */

   //  Check to see if crd is valid
   for (i=0; (i<natom)&&(!B_outside); i++) {
      B_outside = is_out_grid_info(crd[i][0], crd[i][1], crd[i][2]);
   } // i
   
   //if (!B_template){
   // Use standard energy function
   if (!B_outside) {

#ifdef DEBUG
(void)fprintf(logFile,"eval.cc/All coordinates are inside grid...\n");
#endif /* DEBUG */

        // formerly quicktrilinterp->last 4 arguments are NULL:
        // use NULL_ELEC, NULL_EVDW, NULL_ELEC_TOTAL, NULL_EVDW_TOTAL);
        energy = trilinterp( 0, natom, crd, charge, abs_charge, type, map, 
                             info, B_outside?SOME_ATOMS_OUTSIDE_GRID:ALL_ATOMS_INSIDE_GRID, 
                             ignore_inter, elec, emap, &elec_total, &emap_total);
                             //ignore_inter, NULL_ELEC, NULL_EVDW, NULL_ELEC_TOTAL, NULL_EVDW_TOTAL);

        energy += eintcal( nonbondlist, ptr_ad_energy_tables, crd, Nnb,
                       B_calcIntElec, B_include_1_4_interactions,
                       scale_1_4, qsp_abs_charge,
                       B_use_non_bond_cutoff, B_have_flexible_residues);
         
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
        }/*if isGaussTorCon*/
    } else {  //not B_compute_intermol_energy
        //gmm: outsidetrilinterp-> ??set SOME_ATOMS_OUTSIDE_GRID to TRUE
        //energy = outsidetrilinterp( crd, charge, type, natom, map,
        //                            info->inv_spacing, // eval_elec, eval_emap, 
        //                            info->lo[0], info->lo[1], info->lo[2],
        //                            info->hi[0], info->hi[1], info->hi[2],  
        //                            info->cen[0], info->cen[1], info->cen[2] )
        energy = trilinterp( 0, natom, crd, charge, abs_charge, type, map, 
                         info, B_outside?SOME_ATOMS_OUTSIDE_GRID:ALL_ATOMS_INSIDE_GRID, 
                         ignore_inter, elec, emap, &elec_total, &emap_total);
        
        energy += eintcal( nonbondlist, ptr_ad_energy_tables, crd, Nnb,
                   B_calcIntElec, B_include_1_4_interactions,
                   scale_1_4, qsp_abs_charge,
                   B_use_non_bond_cutoff, B_have_flexible_residues);


                    // + eintcal( nonbondlist, e_internal, crd, type, Nnb,
                    //            B_calcIntElec, q1q2);
            if (B_isGaussTorCon) {
                for (I_tor = 0; I_tor <= stateNow.ntor; I_tor++) {
                    if (B_isTorConstrained[I_tor] == 1) {
                        indx = RadiansToDivs( WrpModRad(stateNow.tor[I_tor]) );
                        if (B_ShowTorE) {
                            energy += (double)(US_TorE[I_tor] = US_torProfile[I_tor][indx
    ]);
                        } else {
                            energy += (double)US_torProfile[I_tor][indx];
                        }
                    }
                } // I_tor
            }; // if
       }; //else
    //} else {
    //    // Use template scoring function
    //    if (!B_outside) {
    //        energy = template_trilinterp( crd, charge, type, natom, map, inv_spacing,
    //                              xlo, ylo, zlo, template_energy, template_stddev);
    //    } else {
    //        energy = outside_templ_trilinterp( crd, charge, type, natom, map,
    //                                           inv_spacing, 
    //                                           xlo, ylo, zlo,
    //                                           xhi, yhi, zhi,  xcen, ycen, zcen,
    //                                           template_energy, template_stddev);
    //    }
    //}

   num_evals++;

   if (!finite(energy)) {
      (void)fprintf( logFile, "eval.cc:  ERROR!  energy is infinite!\n\n");
      for (i=0; i<natom; i++) {
            print_PDBQ_atom_resstr( logFile, "", i,   " C   INF     1 ", crd, 
             0.0, 0.0, charge[i],"\n"); 
      } // i
   }
   if (ISNAN(energy)) {
      (void)fprintf( logFile, "eval.cc:  ERROR!  energy is not a number!\n\n");
      for (i=0; i<natom; i++) {
            print_PDBQ_atom_resstr( logFile, "", i,   " C   NaN     1 ", crd, 
             0.0, 0.0, charge[i],"\n"); 
      } // i
   }

   return(energy);
}


