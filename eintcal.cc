/*

 $Id: eintcal.cc,v 1.6 2005/03/11 02:11:29 garrett Exp $

*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>

#ifdef EINTCALPRINT
#include <stdio.h>
#endif

#include "eintcal.h"
#include "constants.h"

extern Linear_FE_Model AD4;

#ifndef EINTCALPRINT
/*
 * Calculate internal energy
 */
FloatOrDouble eintcal( const int           nonbondlist[MAX_NONBONDS][MAX_NBDATA],
                       const FloatOrDouble e_internal[NEINT][ATOM_MAPS][ATOM_MAPS],
                       const FloatOrDouble tcoord[MAX_ATOMS][SPACE],
                       const int           Nnb,
                       const Boole         B_calcIntElec,
                       const FloatOrDouble q1q2[MAX_NONBONDS],
                       const Boole         B_include_1_4_interactions,
                       const FloatOrDouble scale_1_4,
                       const FloatOrDouble qsp_abs_charge[MAX_ATOMS],
                       const FloatOrDouble sol_fn[NEINT],
                       const ParameterEntry parameterArray[MAX_MAPS]
                       )

#else
// eintcalPrint [

extern FILE *logFile;

FloatOrDouble eintcalPrint( const int           nonbondlist[MAX_NONBONDS][MAX_NBDATA],
                            const FloatOrDouble e_internal[NEINT][ATOM_MAPS][ATOM_MAPS],
                            const FloatOrDouble tcoord[MAX_ATOMS][SPACE],
                            const int           Nnb,
                            const Boole         B_calcIntElec,
                            const FloatOrDouble q1q2[MAX_NONBONDS],
                            const Boole         B_include_1_4_interactions,
                            const FloatOrDouble scale_1_4,
                            const FloatOrDouble qsp_abs_charge[MAX_ATOMS],
                            const FloatOrDouble sol_fn[NEINT],
                            const ParameterEntry parameterArray[MAX_MAPS]
                            )

// eintcalPrint ]
#endif

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
/*     Inputs: nonbondlist, e_internal, tcoord, type, Nnb                      */
/*    Returns: eint                                                            */
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

#ifndef  EINTCALPRINT
//  eintcal [
#ifndef  NOSQRT
    register double r; //  SQRT
#endif
//  eintcal ]
#else
// eintcalPrint [
    register double epair=0.0L;
#ifndef  NOSQRT
    register double d; //  SQRT 
#endif
// eintcalPrint ]
#endif

    register int inb;
    register double eint=0.0L, dx, dy, dz;
    register double dint=0.0L;
    register double r2 = 0.0L;
    register int a1, a2;
    register int t1, t2; // Xcode-gmm
    register int nonbond_type; // if = 4, it is a 1_4;  otherwise it is another kind of nonbond
    register int index;


#ifdef EINTCALPRINT
    pr( logFile, "Non-bond  Atom1-Atom2  Distance  Energy\n"); // eintcalPrint 
#endif

    // Loop over all the non-bonds, "inb",
    for (inb = 0;  inb < Nnb;  inb++) {

#ifdef EINTCALPRINT
        epair = 0.0L;
#endif
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
        // ______________________________________________________________________
        // Square-rooting version, slower...

        r = clamp(hypotenuse(dx,dy,dz), RMIN_ELEC); // clamp prevents electrostatic potential becoming too high

        if (B_calcIntElec) {
            //  Calculate  Electrostatic Energy
#   ifndef EINTCALPRINT
            eint += q1q2[inb] / (r*r);     // eintcal
#   else
            epair = q1q2[inb] / (r*r);     // eintcalPrint
#   endif
        }

        // if we are defining USE_8A_CUTOFF, then NBC = 8
        // so, if r is less than the non-bond-cutoff,
        if (r < NBC) {

            // r, the separation between the atoms a1 and a2 in this
            // non-bond, inb, are within the non-bond cutoff, NBC;
            // so calculate the van der Waals and/or H-bonding energy;
            // also, calculate the desolvation energy.

            index = Ang_to_index(r);

            //| Look up the desolvation parameters
            //|
            //| const double qsolpar = 0.01097L;
            //|
            //| desolvation energy = sol_fn[dist] * { rec.vol * (lig.solpar + qsolpar * |lig.charge|)
            //|                                       +
            //|                                       lig.vol * (rec.solpar + qsolpar * |rec.charge|) };
            //|
            //| lig.solpar = parameterArray[t1].solpar;
            //| lig.vol    = parameterArray[t1].vol;
            //| lig.charge = qsp_abs_charge[a1]/qsolpar;
            //| rec.solpar = parameterArray[t2].solpar;
            //| rec.vol    = parameterArray[t2].vol;
            //| rec.charge = qsp_abs_charge[a2]/qsolpar;
            
            dint = sol_fn[index] * ( parameterArray[t2].vol * (parameterArray[t1].solpar + qsp_abs_charge[a1]) 
                                    + parameterArray[t1].vol * (parameterArray[t2].solpar + qsp_abs_charge[a2]) );

#   ifndef EINTCALPRINT
            //|  eintcal [
            //| Xcode-gmm -- only do the double-to-int conversion if within nonbond cutoff
            if (B_include_1_4_interactions != 0 && nonbond_type == 4) {
                //| Compute a scaled 1-4 interaction, multiply by scale_1_4
                eint += scale_1_4 * (e_internal[SqAng_to_index_Int(r2)][t2][t1] + dint);
            } else {
                eint += e_internal[SqAng_to_index_Int(r2)][t2][t1] + dint;
            }
            //  eintcal ]
#   else
            // eintcalPrint [
            if (B_include_1_4_interactions != 0 && nonbond_type == 4) {
                // Compute a scaled 1-4 interaction, multiply by scale_1_4
                epair = scale_1_4 * (e_internal[SqAng_to_index_Int(r2)][t2][t1] + dint);
            } else {
                epair = e_internal[SqAng_to_index_Int(r2)][t2][t1] + dint;
            }
            // eintcalPrint ]
#   endif

#   ifndef EINTCALPRINT
            if (B_include_1_4_interactions != 0 && nonbond_type == 4) {
                // Compute a scaled 1-4 interaction, multiply by scale_1_4
                eint += (scale_1_4 * e_internal[BoundedNeint(index)][t2][t1]); //  eintcal 
            } else {
                eint += (e_internal[BoundedNeint(index)][t2][t1]); //  eintcal 
            }
#   else
            if (B_include_1_4_interactions != 0 && nonbond_type == 4) {
                // Compute a scaled 1-4 interaction, multiply by scale_1_4
                epair += (scale_1_4 * e_internal[BoundedNeint(index)][t2][t1]); // eintcalPrint 
            } else {
                epair += (e_internal[BoundedNeint(index)][t2][t1]); // eintcalPrint 
            }
#   endif
        } // if (r < NBC)
        // ______________________________________________________________________
// SQRT  ]

#else

// NOSQRT [
        // ______________________________________________________________________
        //  Non-square-rooting version, faster...

        r2 = sqhypotenuse(dx,dy,dz);
        r2 = clamp(r2, RMIN_ELEC2);

        if (B_calcIntElec) {
            //  Calculate  Electrostatic  Energy
#   ifndef EINTCALPRINT
            eint += q1q2[inb] / (r2);     // eintcal
#   else
            epair = q1q2[inb] / (r2);     // eintcalPrint
#   endif
        }

        //| If we are defining USE_8A_CUTOFF, then NBC = 8 and NBC2 = 64
        //| so, if r-squared is less than non-bond-cutoff-squared,
        //| i.e.  r2 < NBC2,  so  r < NBC
        if (r2 < NBC2) {  // Xcode-gmm

            //| r, the separation between the atoms a1 and a2 in this
            //| non-bond, inb, are within the non-bond cutoff, NBC;
            //| so calculate the van der Waals and/or H-bonding energy;
            //| also, calculate the desolvation energy.

            index = SqAng_to_index(r2);

            //| Look up the desolvation parameters
            //|
            //| desolvation energy = sol_fn[dist] * { rec.vol * (lig.solpar + qsolpar * |lig.charge|)
            //|                                       +
            //|                                       lig.vol * (rec.solpar + qsolpar * |rec.charge|) };
            //|
            //| lig.solpar = parameterArray[t1].solpar;
            //| lig.vol    = parameterArray[t1].vol;
            //| lig.charge = qsp_abs_charge[a1]/qsolpar;
            //| rec.solpar = parameterArray[t2].solpar;
            //| rec.vol    = parameterArray[t2].vol;
            //| rec.charge = qsp_abs_charge[a2]/qsolpar;
            
            dint = sol_fn[index] * ( parameterArray[t2].vol * (parameterArray[t1].solpar + qsp_abs_charge[a1]) 
                                    + parameterArray[t1].vol * (parameterArray[t2].solpar + qsp_abs_charge[a2]) );

#   ifndef EINTCALPRINT
            //|  eintcal [
            //| Xcode-gmm -- only do the double-to-int conversion if within nonbond cutoff
            if (B_include_1_4_interactions != 0 && nonbond_type == 4) {
                //| Compute a scaled 1-4 interaction, multiply by scale_1_4
                eint += scale_1_4 * (e_internal[SqAng_to_index_Int(r2)][t2][t1] + dint);
            } else {
                eint += e_internal[SqAng_to_index_Int(r2)][t2][t1] + dint;
            }
            //  eintcal ]
#   else
            // eintcalPrint [
            if (B_include_1_4_interactions != 0 && nonbond_type == 4) {
                // Compute a scaled 1-4 interaction, multiply by scale_1_4
                epair = scale_1_4 * (e_internal[SqAng_to_index_Int(r2)][t2][t1] + dint);
            } else {
                epair = e_internal[SqAng_to_index_Int(r2)][t2][t1] + dint;
            }
            // eintcalPrint ]
#   endif

        }   // otherwise, the atoms are too far apart, so don't add anything to the total energy, & save some time
        // ______________________________________________________________________

// NOSQRT ]
#endif


#ifdef EINTCALPRINT
// eintcalPrint [
        eint += epair;
        pr( logFile, " %6d   %5d-%-5d  %7.2lf  %+8.3lf\n", 
                (int)(inb+1), (int)(a1+1), (int)(a2+1), (double)sqrt(r2), (double)epair);
// eintcalPrint ]
#endif

    } //  inb -- next non-bond interaction

#ifdef EINTCALPRINT
    pr( logFile, "\n\nIntramolecular Interaction Energy = %+8.3lf\n", (double)eint); // eintcalPrint
#endif

    return (FloatOrDouble)eint;
}
/*  EOF */
