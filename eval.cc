/********************************************************************
     These are the functions associated with the evaluation object.

                                rsh 9/95
********************************************************************/
// #include <iostream.h>
// #include <fstream.h>
#include <math.h>
#include "eval.h"

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

       | x | y | z | nx | ny | nz | ang | tor1 | ... | tor N |

    where
       x is the x translation
       y is the y translation
       z is the z translation
       nx is the x component of the quaternion
       ny is the y component of the quaternion
       nz is the z component of the quaternion
       ang is the angle portion of the quaternion
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
   stateNow->T.x = rep[0]->gene(0).real;
   stateNow->T.y = rep[1]->gene(0).real;
   stateNow->T.z = rep[2]->gene(0).real;

   //  Set up the quaternion
   stateNow->Q.nx = rep[3]->gene(0).real;
   stateNow->Q.ny = rep[3]->gene(1).real;
   stateNow->Q.nz = rep[3]->gene(2).real;
   stateNow->Q.ang = rep[4]->gene(0).real;
   
   //  Copy the angles
   for (i=1; i<=stateNow->ntor; i++) {
      stateNow->tor[i-1] = rep[4]->gene(i).real;
   }

   mkUnitQuat(&(stateNow->Q));
}

double Eval::operator()(Representation **rep)
{
   register int i;
   int   B_outside = 0;
   int   I_tor = 0;
   int   indx = 0;
   double energy = 0.0;

#ifdef DEBUG
    (void)fprintf(logFile,"eval.cc/double Eval::operator()(Representation **rep)\n");
    if (B_template) {
        (void)fprintf(logFile,"eval.cc/double Eval::operator() -- B_template is true.\n");
    }
#endif /* DEBUG */

   make_state_from_rep(rep, &stateNow);

#ifdef DEBUG
    if (is_out_grid(stateNow.T.x, stateNow.T.y, stateNow.T.z)) {
       (void)fprintf(logFile,"eval.cc/stateNow.T is outside grid!\n");
    }
#endif /* DEBUG */

#ifdef DEBUG
    (void)fprintf(logFile,"eval.cc/Converting state to coordinates...\n");
#endif /* DEBUG */
 
   // Ligand could be inside or could still be outside, check all the atoms...
   cnv_state_to_coords(stateNow, vt, tlist, stateNow.ntor, crdpdb, crd, natom);

#ifdef DEBUG
(void)fprintf(logFile,"eval.cc/Checking to see if all coordinates are inside grid...\n");
#endif /* DEBUG */

   //  Check to see if crd is valid
   for (i=0; (i<natom)&&(!B_outside); i++) {
      B_outside = is_out_grid(crd[i][0], crd[i][1], crd[i][2]);
   } // i

   if (!B_template) {
       // Use standard energy function
       if (!B_outside) {

#ifdef DEBUG
(void)fprintf(logFile,"eval.cc/All coordinates are inside grid...\n");
#endif /* DEBUG */

    /*---------------------------------------------------------------------------
    * I have removed this call to save a stack-push, and inlined the code 
    * instead...
    *                 -- Garrett
    *
    *    energy = evaluate_energy(crd, charge, type, natom, map, inv_spacing, 
    *                              xlo, ylo, zlo, nonbondlist, 
    *                              e_internal, Nnb, B_calcIntElec, q1q2, 
    *                              B_isGaussTorCon, B_isTorConstrained, stateNow, 
    *                              B_ShowTorE, US_TorE, US_torProfile);
    * --------------------------------------------------------------------------*/

            energy = quicktrilinterp( crd, charge, type, natom, map, inv_spacing,
                                      xlo, ylo, zlo)
                     + eintcal( nonbondlist, e_internal, crd, type, Nnb, 
                                B_calcIntElec, q1q2);
            /*
            energy = trilinterp(      crd, charge, type, natom, map, inv_spacing, 
                                      eval_elec, eval_emap, 
                                      xlo, ylo, zlo)
                     + eintcal( nonbondlist, e_internal, crd, type, Nnb, 
                                B_calcIntElec, q1q2);
            */
         
            if (B_isGaussTorCon) {
                for (I_tor = 0; I_tor <= stateNow.ntor; I_tor++) {
                    if (B_isTorConstrained[I_tor] == 1) {
                        indx = Rad2Div( WrpModRad(stateNow.tor[I_tor]) );
                        if (B_ShowTorE) {
                            energy += (double)(US_TorE[I_tor] = US_torProfile[I_tor][indx]);
                        } else {
                            energy += (double)US_torProfile[I_tor][indx];
                        }
                    }
                } // I_tor
            }/*if*/
           } else {
            /*
             * This confuses the GA and GA-LS, because there is no gradient
             * information when all outside conformations are given the same
             * energy.
             *
             * energy = BIG_ENERGY;  / / A really big number defined in autocomm.h
             */
            /*
             * Instead...
             *
             * Penalise atoms outside grid based on the square of the 
             * distance from centre of grid map, otherwise use the normal 
             * trilinear interpolation.
             */
            energy = outsidetrilinterp( crd, charge, type, natom, map,
                                        inv_spacing, // eval_elec, eval_emap, 
                                        xlo, ylo, zlo,
                                        xhi, yhi, zhi,  xcen, ycen, zcen )
                     + eintcal( nonbondlist, e_internal, crd, type, Nnb,
                                B_calcIntElec, q1q2);
            if (B_isGaussTorCon) {
                for (I_tor = 0; I_tor <= stateNow.ntor; I_tor++) {
                    if (B_isTorConstrained[I_tor] == 1) {
                        indx = Rad2Div( WrpModRad(stateNow.tor[I_tor]) );
                        if (B_ShowTorE) {
                            energy += (double)(US_TorE[I_tor] = US_torProfile[I_tor][indx
    ]);
                        } else {
                            energy += (double)US_torProfile[I_tor][indx];
                        }
                    }
                } // I_tor
            } // if
        }
    } else {
        // Use template scoring function
        if (!B_outside) {
            energy = template_trilinterp( crd, charge, type, natom, map, inv_spacing,
                                  xlo, ylo, zlo, template_energy, template_stddev);
        } else {
            energy = outside_templ_trilinterp( crd, charge, type, natom, map,
                                               inv_spacing, 
                                               xlo, ylo, zlo,
                                               xhi, yhi, zhi,  xcen, ycen, zcen,
                                               template_energy, template_stddev);
        }
    }

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
   return(energy);
}

int Eval::write(FILE *out_file, Representation **rep)
{
    int i, retval;
    //char rec14[14];

#ifdef DEBUG
    (void)fprintf(logFile,"eval.cc/int Eval::write(FILE *out_file, Representation **rep)\n");
#endif /*DEBUG*/

    make_state_from_rep(rep, &stateNow);
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
