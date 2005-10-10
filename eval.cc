/*

 $Id: eval.cc,v 1.11.6.1 2005/10/10 16:45:52 alther Exp $

*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/********************************************************************
     These are the functions associated with the evaluation object.

                                rsh 9/95
********************************************************************/
#include <stdio.h>

#include "eval.h"
#include "cnv_state_to_coords.h"
#include "eintcal.h"
#include "qmultiply.h"
#include "trilinterp.h"

extern FILE *logFile;

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
   make_state_from_rep(rep, &stateNow);
   return eval();
}

double Eval::eval()
{
   register int i;
   int   B_outside = 0;
   int   I_tor = 0;
   int   indx = 0;
   double energy = 0.0L;

#ifdef DEBUG
    (void)fprintf(logFile,"eval.cc/double Eval::eval()\n");
    if (B_template) {
        (void)fprintf(logFile,"eval.cc/double Eval::eval() -- B_template is true.\n");
    }
#endif /* DEBUG */

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
      B_outside = is_out_grid_info(crd[i][0], crd[i][1], crd[i][2]);
   } // i

   if (!B_template) {
       // Use standard energy function
       if (!B_outside) {

#ifdef DEBUG
(void)fprintf(logFile,"eval.cc/All coordinates are inside grid...\n");
#endif /* DEBUG */

            energy = quicktrilinterp4( crd, charge, abs_charge, type, natom, map,
                             		   ignore_inter, info);
#ifdef DEBUG
    (void)fprintf(logFile,"eval.cc/double Eval::eval() after quicktrilinterp, energy= %.5lf\n",energy);
#endif /* DEBUG */
            energy += eintcal( nonbondlist, ptr_ad_energy_tables, crd, Nnb, B_calcIntElec, q1q2,
                               B_include_1_4_interactions, scale_1_4,
                               qsp_abs_charge, parameterArray,
                               unbound_internal_FE);
#ifdef DEBUG
    (void)fprintf(logFile,"eval.cc/double Eval::eval() after eintcal, energy= %.5lf\n",energy);
#endif /* DEBUG */

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
            energy = outsidetrilinterp4( crd, charge, abs_charge, type, natom, map, ignore_inter, info );
#ifdef DEBUG
    (void)fprintf(logFile,"eval.cc/double Eval::eval() after outsidetrilinterp, energy= %.5lf\n",energy);
#endif /* DEBUG */
            energy += eintcal( nonbondlist, ptr_ad_energy_tables, crd, Nnb, B_calcIntElec, q1q2, B_include_1_4_interactions, scale_1_4, qsp_abs_charge, parameterArray, unbound_internal_FE);
#ifdef DEBUG
    (void)fprintf(logFile,"eval.cc/double Eval::eval() after eintcal, energy= %.5lf\n",energy);
#endif /* DEBUG */
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
            energy = template_trilinterp( crd, charge, abs_charge, type, natom, map,
                                  template_energy, template_stddev, info);
        } else {
            energy = outside_templ_trilinterp( crd, charge, abs_charge, type, natom, map,
                                               template_energy, template_stddev, info);
        }
    }

   num_evals++;

   if (!FINITE(energy)) {
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
    (void)fprintf(logFile,"eval.cc/double Eval::eval() returns energy= %.5lf\n",energy);
#endif /*DEBUG*/
   return(energy);
}

int Eval::write(FILE *out_file, Representation **rep)
{
    int i=0, retval=0;
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

#if defined(USING_COLINY)
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
now->Q.nx = rep[3];
now->Q.ny = rep[4];
now->Q.nz = rep[5];
now->Q.ang = rep[6];

//  Copy the angles
now->ntor = n - 7;
for (int i=0, j=7; j<n; i++, j++)
  now->tor[i] = rep[j];

mkUnitQuat(&(now->Q));
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
#endif // USING_COLINY
