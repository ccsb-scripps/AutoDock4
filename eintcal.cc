/*

 $Id: eintcal.cc,v 1.3 2004/02/12 04:32:15 garrett Exp $

*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/* eintcal.cc */

#include <math.h>

#ifdef EINTCALPRINT        /*EINTCALPRINT[*/
    #include <stdio.h>
#endif                     /*EINTCALPRINT]*/
    #include "eintcal.h"
    #include "constants.h"

#ifndef EINTCALPRINT        /*!EINTCALPRINT[*/

FloatOrDouble eintcal( int           nonbondlist[MAX_NONBONDS][2],
                       FloatOrDouble e_internal[NEINT][ATOM_MAPS][ATOM_MAPS],
               FloatOrDouble tcoord[MAX_ATOMS][SPACE],
               int type[MAX_ATOMS],
               int Nnb,
               Boole B_calcIntElec,
               FloatOrDouble q1q2[MAX_NONBONDS])

/******************************************************************************/
/*      Name: eintcal                                                         */
/*  Function: Calculate the Internal Energy of the Small Molecule.            */
/*            Accelerated non-square-rooting, dx,dy,dz version.               */
/*  Copyright: (C) 1994-2004, TSRI                                            */
/*____________________________________________________________________________*/
/*   Authors: Garrett M. Morris, TSRI                                         */
/*            David Goodsell, UCLA                                            */
/*      Date: 16/03/94                                                        */
/* ____________________________________________________________________________*/
/*     Inputs: nonbondlist, e_internal, tcoord, type, Nnb                      */
/*   Returns: eint                                                            */
/*   Globals: NEINT, MAX_ATOMS, SPACE                                         */
/*____________________________________________________________________________*/
/* Modification Record                                                        */
/* Date     Inits   Comments                                                  */
/* 07/05/92 DSG     Original FORTRAN                                          */
/* 15/05/92 GMM     Translated into C                                         */
/* 15/05/92 GMM     hypotenuse macro                                          */
/* 19/11/93 GMM     Accelerated non-square-rooting version.                   */
/* 16/03/94 GMM     Accelerated dx,dy,dz version.                             */
/*  10/02/04 GMM     Reduced NBC from 64.0 to 8.0 Å							   */
/* *****************************************************************************/

                        /*!EINTCALPRINT]*/
#else                        /*EINTCALPRINT[*/

extern FILE *logFile;

FloatOrDouble eintcalPrint( int nonbondlist[MAX_NONBONDS][2],
                            FloatOrDouble e_internal[NEINT][ATOM_MAPS][ATOM_MAPS],
                    FloatOrDouble tcoord[MAX_ATOMS][SPACE],
                    int type[MAX_ATOMS],
                    int Nnb,
                    Boole B_calcIntElec,
                    FloatOrDouble q1q2[MAX_NONBONDS])

#endif                       /*EINTCALPRINT]*/


{
#ifndef  EINTCALPRINT        
/*  not EINTCALPRINT [ */
#ifndef  NOSQRT              
    double r; /*  SQRT */                     
#endif                       
/*  not EINTCALPRINT ] */
#else
/* EINTCALPRINT [ */
    double epair=0.0L;
#ifndef  NOSQRT
    double d; /*  SQRT  */
#endif
/* EINTCALPRINT ] */
#endif

    register int inb;
    double eint=0.0L, dx, dy, dz;
    double r2 = 0.0L;
    int a1, a2;
    int t1, t2; // Xcode-gmm
    int index;
        
#ifdef EINTCALPRINT          
    pr( logFile, "Non-bond  Atom1-Atom2  Distance  Energy\n"); /* EINTCALPRINT  */
#endif                        

    for (inb = 0;  inb < Nnb;  inb++) {

#pragma function_align 32
		
        a1 = nonbondlist[inb][ATM1];
        a2 = nonbondlist[inb][ATM2];
        t1 = type[a1]; // Xcode-gmm
        t2 = type[a2]; // Xcode-gmm

        dx = tcoord[a1][X] - tcoord[a2][X];
        dy = tcoord[a1][Y] - tcoord[a2][Y];
        dz = tcoord[a1][Z] - tcoord[a2][Z];

#ifndef NOSQRT              
/* SQRT  [ */
        /*  NOSQRT is _not_ defined, i.e. SQRTing version, which is slower... */
        if (B_calcIntElec) {
            /*** Calculate internal electrostatic energy too ***/
            /* r = hypotenuse(dx,dy,dz); */

            r = clamp(hypotenuse(dx,dy,dz), RMIN_ELEC);

            index = Ang_to_index(r);
#ifndef EINTCALPRINT        
            eint += (e_internal[ BoundedNeint(index) ][t2][t1] + q1q2[inb]/(r*r)); /*  not EINTCALPRINT  */
#else                       
            epair = (e_internal[ BoundedNeint(index) ][t2][t1] + q1q2[inb]/(r*r)); /* EINTCALPRINT  */
#endif                      

        } else {
            /*** Calculate van der Waals/H-bond energy only ***/

#ifdef BOUNDED              /*BOUNDED[*/
            index = Ang_to_index(hypotenuse(dx,dy,dz));
#ifndef EINTCALPRINT        /*!EINTCALPRINT[*/
            eint += eint_table[ BoundedNeint(index) ][type[a2]][type[a1]];
                            /*!EINTCALPRINT]*/
#else                       /*EINTCALPRINT[*/
            epair = eint_table[ BoundedNeint(index) ][type[a2]][type[a1]];
#endif                      /*EINTCALPRINT]*/
                            /*BOUNDED]*/
#else                       /*!BOUNDED[*/
            assert( Ang_to_index(hypotenuse(dx,dy,dz)) >= 0 );
            assert( Ang_to_index(hypotenuse(dx,dy,dz)) < NEINT );
#ifndef EINTCALPRINT        /*!EINTCALPRINT[*/
            eint += eint_table[ Ang_to_index(hypotenuse(dx,dy,dz)) ][type[a2]][type[a1]];
                            /*!EINTCALPRINT]*/
#else                       /*EINTCALPRINT[*/
            epair = eint_table[ Ang_to_index(hypotenuse(dx,dy,dz)) ][type[a2]][type[a1]];
#endif                      /*EINTCALPRINT]*/
#endif                      /*!BOUNDED]*/

        }
                            /*!NOSQRT]*/
#else                       /*NOSQRT[*/
        /* NOSQRTing version, faster... */
        if (B_calcIntElec) {
            r2 = sqhypotenuse(dx,dy,dz);
            r2 = clamp(r2, RMIN_ELEC2);
            index = SqAng_to_index(r2);
#ifndef EINTCALPRINT
            eint += (e_internal[ BoundedNeint(index) ][t2][t1] + q1q2[inb]/r2); /*  not EINTCALPRINT  */
#else
            epair = (e_internal[ BoundedNeint(index) ][t2][t1] + q1q2[inb]/r2); /* EINTCALPRINT  */
#endif
        } else {

            r2 = sqhypotenuse(dx,dy,dz);

#ifndef EINTCALPRINT
/*  not EINTCALPRINT [ */
            if (r2 < NBC2) {  // Xcode-gmm
                // atom pair is close enough to interact
                // Xcode-gmm -- only do the double-to-int conversion if within nonbond cutoff
                eint += e_internal[SqAng_to_index_Int(r2)][t2][t1];
            }   // otherwise, the atoms are too far apart, so don't add anything to the total energy, & save some time
/*  not EINTCALPRINT ] */
#else
/* EINTCALPRINT [ */
            if (r2 < NBC2) {  // Xcode-gmm
                epair = e_internal[SqAng_to_index_Int(r2)][t2][t1];
            } else {
                epair = 0.0L;
            }
/* EINTCALPRINT ] */
#endif
        }
#endif                      /*NOSQRT]*/

#ifdef EINTCALPRINT         /*EINTCALPRINT[*/
        eint += epair;
        pr( logFile, " %6d   %5d-%-5d  %7.2lf  %+8.3lf\n", (int)(inb+1), (int)(a1+1),
        (int)(a2+1), (double)sqrt(r2), (double)epair);
#endif                      /*EINTCALPRINT]*/

    } /* next non-bond interaction */

#ifdef EINTCALPRINT         /*EINTCALPRINT[*/
    pr( logFile, "\n\nIntramolecular Interaction Energy = %+8.3lf\n", (double)eint);
#endif                      /*EINTCALPRINT]*/

    return (FloatOrDouble)eint;
}
/* EOF */
