/*

 $Id: eintcal.cc,v 1.21 2009/07/14 23:42:33 rhuey Exp $

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

#include <stdio.h>
#include <math.h>
#include "eintcal.h"
#include "constants.h"
#include "distdepdiel.h"


extern Linear_FE_Model AD4;
extern int Nnb_array[3];
extern Real nb_group_energy[3];

#ifndef EINTCALPRINT

// Calculate internal energy
Real eintcal( NonbondParam * const nonbondlist,
              const EnergyTables  *ptr_ad_energy_tables,
              const Real tcoord[MAX_ATOMS][SPACE],
              const int           Nnb,
              const Boole         B_calcIntElec,
              const Boole         B_include_1_4_interactions,
              const Real scale_1_4,
              const Real qsp_abs_charge[MAX_ATOMS],
              const Boole B_use_non_bond_cutoff,
              const Boole B_have_flexible_residues  // if the receptor has flexibile residues, this will be set to TRUE
             )

#else 

// eintcalPrint [

extern FILE *logFile;

// Calculate internal energy and print out a detailed report
Real eintcalPrint( NonbondParam * const nonbondlist,
                   const EnergyTables  *ptr_ad_energy_tables,
                   const Real tcoord[MAX_ATOMS][SPACE],
                   const int           Nnb,
                   const Boole         B_calcIntElec,
                   const Boole         B_include_1_4_interactions,
                   const Real scale_1_4,
                   const Real qsp_abs_charge[MAX_ATOMS],
                   const Boole B_use_non_bond_cutoff,
                   const Boole B_have_flexible_residues  // if the receptor has flexibile residues, this will be set to TRUE
                  )
// eintcalPrint ]

#endif 

/* *****************************************************************************/
/*       Name: eintcal                                                         */
/*   Function: Calculate the Internal Energy of the Small Molecule.            */
/*             Accelerated non-square-rooting, dx,dy,dz version.               */
/*Copyright (C) 2009 The Scripps Research Institute. All rights reserved. */
/* ____________________________________________________________________________*/
/*    Authors: Garrett M. Morris, TSRI                                         */
/*             David Goodsell, UCLA                                            */
/*       Date: 16/03/94                                                        */
/* ____________________________________________________________________________*/
/*     Inputs: nonbondlist, ptr_ad_energy_tables, tcoord, type, Nnb            */
/*    Returns: total_e_internal                                                */
/*    Globals: NEINT, MAX_ATOMS, SPACE                                         */
/* ____________________________________________________________________________*/
/*  Modification Record                                                        */
/*  Date     Inits   Comments                                                  */
/*  07/05/92 DSG     Original FORTRAN                                          */
/*  15/05/92 GMM     Translated into C                                         */
/*  15/05/92 GMM     hypotenuse macro                                          */
/*  19/11/93 GMM     Accelerated non-square-rooting version.                   */
/*  16/03/94 GMM     Accelerated dx,dy,dz version.                             */
/*  10/02/04 GMM     Reduced NBC from 64.0 to 8.0                              */
/*  04/03/05 GMM     Added the new internal desolvation term                   */
/* *****************************************************************************/

{

  // if r is less than the non-bond-cutoff, 
  //  -OR-
  // If we are computing the unbound conformation then we ignore the non bond cutoff, NBC
#ifndef  EINTCALPRINT //  eintcal [
#   ifndef  NOSQRT
        register double r=0.0L; //  SQRT
        // if we have defined USE_8A_CUTOFF, then NBC = 8
        const double nbc  = B_use_non_bond_cutoff ? NBC  : 999;
#   else
        // if we have defined USE_8A_CUTOFF, then NBC = 8 // Xcode-gmm
        const double nbc2 = B_use_non_bond_cutoff ? NBC2 : 999 * 999;
#   endif
    //  eintcal ]
#else // eintcalPrint [
#   ifndef  NOSQRT
        register double d=0.0L; //  SQRT 
        //const double nbc  = NBC;
        // if we have defined USE_8A_CUTOFF, then NBC = 8
        const double nbc  = B_use_non_bond_cutoff ? NBC  : 999;
#   else
        //const double nbc2 = NBC2;
        // if we have defined USE_8A_CUTOFF, then NBC = 8 // Xcode-gmm
        const double nbc2 = B_use_non_bond_cutoff ? NBC2 : 999 * 999;
#   endif
#endif // eintcalPrint ]

    register double dx=0.0L, dy=0.0L, dz=0.0L;
    register double r2=0.0L;

    register double total_e_internal=0.0L; // total_e_internal = eint

    register double e_elec=0.0L;

#ifdef EINTCALPRINT
    double total_e_elec=0.0L;
    double total_e_vdW_Hb=0.0L;
    double e_vdW_Hb=0.0L;
    double total_e_desolv=0.0L;
#endif

    register int inb=0;
    register int a1=0, a2=0;
    register int t1=0, t2=0; // Xcode-gmm
    register int nonbond_type=0; // if = 4, it is a 1_4;  otherwise it is another kind of nonbond
    register int index_lt_NEINT=0;
    register int index_lt_NDIEL=0;
    register int nb_group=0;
    int inb_from=0;
    int inb_to=0;
    int nb_group_max = 1;  // By default, we have one nonbond group, (1) intramolecular in the ligand

    if (B_have_flexible_residues) {
        // If we have flexible residues, we need to consider three groups of nonbonds:
        // (1) intramolecular in the ligand, (2) intermolecular and (3) intramolecular in the receptor
        nb_group_max = 3;
    }

    // Loop over the nonbonding groups --
    // Either (intramolecular ligand nonbonds)
    // or (intramolecular ligand nonbonds, intermolecular nonbonds, and intramolecular receptor nonbonds)
    for (nb_group = 0;  nb_group < nb_group_max;  nb_group++) {

#ifdef EINTCALPRINT
        if (nb_group == 0) {
            pr(logFile, "\n\n\t\tLigand Intramolecular Energy Analysis\n");
            pr(logFile,     "\t\t=====================================\n\n");
        }
        if (nb_group == 1) {
            pr(logFile, "\n\n\t\tLigand-Receptor Moving-Atom Intermolecular Energy Analysis\n");
            pr(logFile,     "\t\t==========================================================\n\n");
        }
        if (nb_group == 2) {
            pr(logFile, "\n\n\t\tReceptor Moving-Atom Intramolecular Energy Analysis\n");
            pr(logFile,     "\t\t===================================================\n\n");
        }
        if (B_calcIntElec) {
            pr( logFile, "Non-bond  Atom1-Atom2  Distance   Total     Elec      vdW+Hb    Desolv     Sol_fn   Type Dielectric\n"); // eintcalPrint 
            pr( logFile, "________  ___________  ________   ______  ________  ________  ________   ________   ____ __________\n"); // eintcalPrint 
        } else {
            pr( logFile, "Non-bond  Atom1-Atom2  Distance   Total     vdW+Hb    Desolv     Sol_fn   Type Dielectric\n"); // eintcalPrint 
            pr( logFile, "________  ___________  ________   ______  ________  ________   ________   ____ __________\n"); // eintcalPrint 
        }
#endif

        if (nb_group == 0) {
            inb_from = 0;
        } else {
            inb_from = Nnb_array[nb_group-1];
        }
        inb_to   = Nnb_array[nb_group];

        // Loop over the non-bonds in this nonbond "group", "inb",
        for (inb = inb_from;  inb < inb_to;  inb++) {

            double e_internal=0;  // e_internal = epair
            double e_desolv=0;    // e_desolv = dpair

            a1 = nonbondlist[inb].a1;
            a2 = nonbondlist[inb].a2;
            t1 = nonbondlist[inb].t1; // Xcode-gmm  // t1 is a map_index
            t2 = nonbondlist[inb].t2; // Xcode-gmm  // t2 is a map_index
            nonbond_type = nonbondlist[inb].nonbond_type;
            double nb_desolv = nonbondlist[inb].desolv;
            double q1q2 = nonbondlist[inb].q1q2;

            dx = tcoord[a1][X] - tcoord[a2][X];
            dy = tcoord[a1][Y] - tcoord[a2][Y];
            dz = tcoord[a1][Z] - tcoord[a2][Z];

            // Calculate the van der Waals and/or H-bonding energy & the desolvation energy.
            //|
            //| desolvation energy = sol_fn[dist] * ( rec.vol * (lig.solpar + qsolpar * |lig.charge|)
            //|                                     + lig.vol * (rec.solpar + qsolpar * |rec.charge|) );
            //|
#ifndef NOSQRT 
            // Use square-root, slower...

            // r = the separation between the atoms a1 and a2 in this non-bond, inb, 
            r = clamp(hypotenuse(dx,dy,dz), RMIN_ELEC); // clamp prevents electrostatic potential becoming too high when shorter than RMIN_ELEC
            r2 = r*r;
            register const int index = Ang_to_index(r); // convert real-valued distance r to an index for energy lookup tables

#else   
            //  Non-square-rooting version, faster...
            r2 = sqhypotenuse(dx,dy,dz); // r2, the square of the separation between the atoms a1 and a2 in this non-bond, inb, 
            r2 = clamp(r2, RMIN_ELEC2);
            register const int index = SqAng_to_index(r2);

#endif  // NOSQRT ]

            index_lt_NEINT = BoundedNeint(index);  // guarantees that index_lt_NEINT is never greater than (NEINT - 1)
            index_lt_NDIEL = BoundedNdiel(index);  // guarantees that index_lt_NDIEL is never greater than (NDIEL - 1)

            if (B_calcIntElec) {
                //  Calculate  Electrostatic  Energy
                double r_dielectric = ptr_ad_energy_tables->r_epsilon_fn[index_lt_NDIEL];
                e_elec = q1q2 * r_dielectric;
                e_internal = e_elec;
            }
            if  ( r2 < nbc2 ) {   
                e_desolv = ptr_ad_energy_tables->sol_fn[index_lt_NEINT] * nb_desolv;
                if (B_include_1_4_interactions != 0 && nonbond_type == 4) {
                    //| Compute a scaled 1-4 interaction,
                    e_internal += scale_1_4 * (ptr_ad_energy_tables->e_vdW_Hb[index_lt_NEINT][t2][t1] + e_desolv);
                } else {
                    e_internal += ptr_ad_energy_tables->e_vdW_Hb[index_lt_NEINT][t2][t1] + e_desolv;
                }
            }


          total_e_internal += e_internal;
#ifdef EINTCALPRINT // eintcalPrint [

          total_e_desolv   += e_desolv;
          total_e_elec     += e_elec;
          double dielectric = ptr_ad_energy_tables->epsilon_fn[index_lt_NDIEL];

          if (B_calcIntElec) {

              e_vdW_Hb = e_internal - e_desolv - e_elec,
              pr( logFile, " %6d   %5d-%-5d  %7.2lf  %+8.3lf  %+8.3lf  %+8.3lf  %+8.3lf   %+8.3lf   %d  %8.3lf\n", 
                    (int)(inb+1), (int)(a1+1), (int)(a2+1), (double)sqrt(r2), 
                    (double)e_internal, (double)e_elec, (double)e_vdW_Hb, (double)e_desolv, 
                    (double)ptr_ad_energy_tables->sol_fn[index_lt_NEINT], (int)nonbond_type, (double)dielectric 
                 );
          } else {

              e_vdW_Hb = e_internal - e_desolv,
              pr( logFile, " %6d   %5d-%-5d  %7.2lf  %+8.3lf  %+8.3lf  %+8.3lf   %+8.3lf   %d  %8.3lf\n", 
                    (int)(inb+1), (int)(a1+1), (int)(a2+1), (double)sqrt(r2), 
                    (double)e_internal, (double)e_vdW_Hb, (double)e_desolv, 
                    (double)ptr_ad_energy_tables->sol_fn[index_lt_NEINT], (int)nonbond_type, (double)dielectric 
                 );
          }

          total_e_vdW_Hb += e_vdW_Hb;

#endif // eintcalPrint ]

        } //  inb -- next non-bond interaction

        if (nb_group == INTRA_LIGAND) { // [0]
            // Intramolecular energy of ligand
            nb_group_energy[INTRA_LIGAND] = total_e_internal;
        } else if (nb_group == INTER) { // [1]
            // intermolecular energy
            nb_group_energy[INTER] = total_e_internal - nb_group_energy[INTRA_LIGAND];
        } else if (nb_group == INTRA_RECEPTOR) { // [2]
            // intramolecular energy of receptor
            nb_group_energy[INTRA_RECEPTOR] = total_e_internal - nb_group_energy[INTRA_LIGAND] - nb_group_energy[INTER];
        }

    } // nb_group -- intra lig, inter, intra rec


#ifdef EINTCALPRINT
    if (B_calcIntElec) {
        pr( logFile, "                                ________  ________  ________  ________\n");
        pr( logFile, "Total                           %+8.3lf  %+8.3lf  %+8.3lf  %+8.3lf\n", total_e_internal, total_e_elec, total_e_vdW_Hb, total_e_desolv);
        pr( logFile, "                                ________  ________  ________  ________\n");
        pr( logFile, "                                   Total      Elec    vdW+Hb    Desolv\n");
    } else {
        pr( logFile, "                                ________  ________  ________\n");
        pr( logFile, "Total                           %+8.3lf  %+8.3lf  %+8.3lf\n", total_e_internal, total_e_vdW_Hb, total_e_desolv);
        pr( logFile, "                                ________  ________  ________\n");
        pr( logFile, "                                   Total    vdW+Hb    Desolv\n");
    }
#endif

#ifdef EINTCALPRINT
    pr( logFile, "\nTotal Intramolecular Interaction Energy   = %+.3lf kcal/mol\n", (double)total_e_internal); // eintcalPrint
#endif

    return (Real) total_e_internal;
}
/*  EOF */
