/*

 $Id: eintcal.cc,v 1.10.6.1 2005/10/10 16:35:51 alther Exp $

*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "eintcal.h"

extern Linear_FE_Model AD4;

/*
 * Calculate internal energy
 */

FloatOrDouble eintcal( const int           nonbondlist[MAX_NONBONDS][MAX_NBDATA],
                       const EnergyTables  *ptr_ad_energy_tables,
                       const FloatOrDouble tcoord[MAX_ATOMS][SPACE],
                       const int           Nnb,
                       const Boole         B_calcIntElec,
                       const FloatOrDouble q1q2[MAX_NONBONDS],
                       const Boole         B_include_1_4_interactions,
                       const FloatOrDouble scale_1_4,
                       const FloatOrDouble qsp_abs_charge[MAX_ATOMS],
                       const ParameterEntry parameterArray[MAX_MAPS],
                       const FloatOrDouble unbound_internal_FE
                       )
/* *****************************************************************************/
/*       Name: eintcal                                                         */
/*   Function: Calculate the Internal Energy of the Small Molecule.            */
/*             Accelerated non-square-rooting, dx,dy,dz version.               */
/*  Copyright: (C) 1994-2004, TSRI                                             */
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

#   ifndef  NOSQRT
        register double r=0.0L; //  SQRT
#   endif
    //  eintcal ]

    register double dx=0.0L, dy=0.0L, dz=0.0L;
    register double r2=0.0L;

    register double total_e_internal=0.0L; // total_e_internal = eint

    register double e_desolv=0.0L; // e_desolv = dpair
    register double e_elec=0.0L;

    register double dielectric = 1.0L;
    register double r_dielectric = 1.0L;

    register int inb=0;
    register int a1=0, a2=0;
    register int t1=0, t2=0; // Xcode-gmm
    register int nonbond_type=0; // if = 4, it is a 1_4;  otherwise it is another kind of nonbond
    register int index=0;

    register int index_lt_NEINT=0;
    register int index_lt_NDIEL=0;

    // Loop over all the non-bonds, "inb",
    for (inb = 0;  inb < Nnb;  inb++) {
        a1 = nonbondlist[inb][ATM1];
        a2 = nonbondlist[inb][ATM2];
        t1 = nonbondlist[inb][TYPE1]; // Xcode-gmm  // t1 is a map_index
        t2 = nonbondlist[inb][TYPE2]; // Xcode-gmm  // t2 is a map_index
        nonbond_type = nonbondlist[inb][NBTYPE];

        dx = tcoord[a1][X] - tcoord[a2][X];
        dy = tcoord[a1][Y] - tcoord[a2][Y];
        dz = tcoord[a1][Z] - tcoord[a2][Z];

#ifndef NOSQRT
// SQRT  [

        // Use square-root, slower...

        // r = the separation between the atoms a1 and a2 in this non-bond, inb,
        r = clamp(hypotenuse(dx,dy,dz), RMIN_ELEC); // clamp prevents electrostatic potential becoming too high when shorter than RMIN_ELEC

        index = Ang_to_index(r); // convert real-valued distance r to an index for energy lookup tables
        index_lt_NEINT = BoundedNeint(index);  // guarantees that index_lt_NEINT is never greater than (NEINT - 1)
        index_lt_NDIEL = BoundedNdiel(index);  // guarantees that index_lt_NDIEL is never greater than (NDIEL - 1)

        dielectric = ptr_ad_energy_tables->epsilon_fn[index_lt_NDIEL];
        // r_dielectric = ptr_ad_energy_tables->r_epsilon_fn[index_lt_NDIEL];
        r_dielectric = r * ptr_ad_energy_tables->epsilon_fn[index_lt_NDIEL];

        if (B_calcIntElec) {
            //  Calculate  Electrostatic Energy
            e_elec = (double) q1q2[inb] / (r_dielectric);
            total_e_internal += e_elec;     // eintcal
        }

        // if r is less than the non-bond-cutoff,
        if (r < NBC) {   // if we have defined USE_8A_CUTOFF, then NBC = 8

            // Calculate the van der Waals and/or H-bonding energy & the desolvation energy.

            //| Calculate the desolvation energy
            //|
            //| desolvation energy = sol_fn[dist] * ( rec.vol * (lig.solpar + qsolpar * |lig.charge|)
            //|                                     + lig.vol * (rec.solpar + qsolpar * |rec.charge|) );
            //|
            //| lig.solpar = parameterArray[t1].solpar;
            //| lig.vol    = parameterArray[t1].vol;
            //| lig.charge = qsp_abs_charge[a1]/qsolpar;
            //| rec.solpar = parameterArray[t2].solpar;
            //| rec.vol    = parameterArray[t2].vol;
            //| rec.charge = qsp_abs_charge[a2]/qsolpar;

            e_desolv = ptr_ad_energy_tables->sol_fn[index_lt_NEINT] *
                       ( parameterArray[t2].vol * (parameterArray[t1].solpar + qsp_abs_charge[a1])
                       + parameterArray[t1].vol * (parameterArray[t2].solpar + qsp_abs_charge[a2]) );

            if (B_include_1_4_interactions != 0 && nonbond_type == 4) {
                // Compute a scaled 1-4 interaction, multiply by scale_1_4
                total_e_internal += scale_1_4 * (ptr_ad_energy_tables->e_vdW_Hb[index_lt_NEINT][t2][t1] + e_desolv);
            } else {
                total_e_internal += ptr_ad_energy_tables->e_vdW_Hb[index_lt_NEINT][t2][t1] + e_desolv;
            }
        } else {
            // otherwise, the atoms are too far apart, so don't add anything to the total energy, & save some time
            e_desolv = 0.0L;
        }

// SQRT  ]
#else   // NOSQRT [
        //  Non-square-rooting version, faster...

        r2 = sqhypotenuse(dx,dy,dz); // r2, the square of the separation between the atoms a1 and a2 in this non-bond, inb,
        r2 = clamp(r2, RMIN_ELEC2);

        index = SqAng_to_index(r2);
        index_lt_NEINT = BoundedNeint(index);  // guarantees that index_lt_NEINT is never greater than (NEINT - 1)
        index_lt_NDIEL = BoundedNdiel(index);  // guarantees that index_lt_NDIEL is never greater than (NDIEL - 1)

        dielectric = ptr_ad_energy_tables->epsilon_fn[index_lt_NDIEL];
        // r_dielectric = ptr_ad_energy_tables->r_epsilon_fn[index_lt_NDIEL];
        r_dielectric = sqrt(r2) * ptr_ad_energy_tables->epsilon_fn[index_lt_NDIEL];  // Using sqrt() in NOSQRT is a no-no!  FIXME!

        if (B_calcIntElec) {
            //  Calculate  Electrostatic  Energy
            e_elec = (double) q1q2[inb] / r_dielectric;
            total_e_internal += e_elec; // eintcal
        }

        //| if r-squared is less than non-bond-cutoff-squared,
        if (r2 < NBC2) {  // Xcode-gmm

            // Calculate the van der Waals and/or H-bonding energy & the desolvation energy.

            //| Calculate the desolvation energy
            //|
            //| desolvation energy = sol_fn[dist] * ( rec.vol * (lig.solpar + qsolpar * |lig.charge|)
            //|                                     + lig.vol * (rec.solpar + qsolpar * |rec.charge|) );
            //|
            //| lig.solpar = parameterArray[t1].solpar;
            //| lig.vol    = parameterArray[t1].vol;
            //| lig.charge = qsp_abs_charge[a1]/qsolpar;
            //| rec.solpar = parameterArray[t2].solpar;
            //| rec.vol    = parameterArray[t2].vol;
            //| rec.charge = qsp_abs_charge[a2]/qsolpar;

            e_desolv = ptr_ad_energy_tables->sol_fn[index_lt_NEINT] *
                       ( parameterArray[t2].vol * (parameterArray[t1].solpar + qsp_abs_charge[a1])
                       + parameterArray[t1].vol * (parameterArray[t2].solpar + qsp_abs_charge[a2]) );

            if (B_include_1_4_interactions != 0 && nonbond_type == 4) {
                //| Compute a scaled 1-4 interaction,
                total_e_internal += scale_1_4 * (ptr_ad_energy_tables->e_vdW_Hb[index_lt_NEINT][t2][t1] + e_desolv);
            } else {
                total_e_internal += ptr_ad_energy_tables->e_vdW_Hb[index_lt_NEINT][t2][t1] + e_desolv;
            }
        } else {
            // otherwise, the atoms are too far apart, so don't add anything to the total energy, & save some time
            e_desolv = 0.0L;
        }
#endif  // NOSQRT ]

    } //  inb -- next non-bond interaction

    // Subtract the internal free energy of the unbound state
    total_e_internal = total_e_internal - unbound_internal_FE;

    return (FloatOrDouble) total_e_internal;
}
/*  EOF */
