/* eintcal.cc */

#include <math.h>

#ifdef EINTCALPRINT        /*EINTCALPRINT[*/
    #include <stdio.h>
#endif                     /*EINTCALPRINT]*/
    #include "eintcal.h"
    #include "constants.h"

#ifndef EINTCALPRINT        /*!EINTCALPRINT[*/

float eintcal( int nonbondlist[MAX_NONBONDS][2],
               float eint_table[NEINT][ATOM_MAPS][ATOM_MAPS],
               float tcoord[MAX_ATOMS][SPACE],
               int type[MAX_ATOMS],
               int Nnb,
               Boole B_calcIntElec,
               float q1q2[MAX_NONBONDS])

/******************************************************************************/
/*      Name: eintcal                                                         */
/*  Function: Calculate the Internal Energy of the Small Molecule.            */
/*            Accelerated non-square-rooting, dx,dy,dz version.               */
/* Copyright: (C) 1994, TSRI                                                  */
/*____________________________________________________________________________*/
/*   Authors: Garrett M. Morris, TSRI                                         */
/*            David Goodsell, UCLA                                            */
/*      Date: 16/03/94                                                        */
/*____________________________________________________________________________*/
/*    Inputs: nonbondlist, eint_table, tcoord, type, Nnb                      */
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
/******************************************************************************/

                        /*!EINTCALPRINT]*/
#else                        /*EINTCALPRINT[*/

extern FILE *logFile;

float eintcalPrint( int nonbondlist[MAX_NONBONDS][2],
                    float eint_table[NEINT][ATOM_MAPS][ATOM_MAPS],
                    float tcoord[MAX_ATOMS][SPACE],
                    int type[MAX_ATOMS],
                    int Nnb,
                    Boole B_calcIntElec,
                    float q1q2[MAX_NONBONDS])

#endif                       /*EINTCALPRINT]*/


{
#ifndef  EINTCALPRINT        /*!EINTCALPRINT[*/

#ifndef  NOSQRT              /*!NOSQRT[*/
    double r;
                             /*!NOSQRT]*/
#else                        /*NOSQRT[*/
    double r2;
#endif                       /*NOSQRT]*/
                             /*!EINTCALPRINT]*/
#else                        /*EINTCALPRINT[*/

    double epair=0.;
    double r2;

#ifndef  NOSQRT              /*!NOSQRT[*/
    double d;
#endif                       /*!NOSQRT]*/
#endif                       /*EINTCALPRINT]*/

    double eint=0., dx,dy,dz;
    int a1, a2;
    register int inb;

#ifdef BOUNDED               /*BOUNDED[*/
    int index;
#endif                       /*BOUNDED]*/

#ifdef EINTCALPRINT          /*EINTCALPRINT[*/
    pr( logFile, "Non-bond  Atom1-Atom2  Distance  Energy\n");
#endif                        /*EINTCALPRINT]*/

    for (inb = 0;  inb < Nnb;  inb++) {

        dx = tcoord[(a1 = nonbondlist[inb][ATM1])][X] - tcoord[(a2 = nonbondlist[inb][ATM2])][X];
        dy = tcoord[a1][Y] - tcoord[a2][Y];
        dz = tcoord[a1][Z] - tcoord[a2][Z];

#ifndef NOSQRT              /*!NOSQRT[*/
        /* NOSQRT is not defined, i.e. SQRTing version, which is slower... */
        if (B_calcIntElec) {
            /*** Calculate internal electrostatic energy too ***/
            /* r = hypotenuse(dx,dy,dz); */

            r = clamp(hypotenuse(dx,dy,dz), RMIN_ELEC);

#ifdef BOUNDED              /*BOUNDED[*/
            index = Ang_to_index(r);
#ifndef EINTCALPRINT        /*!EINTCALPRINT[*/
            eint += (eint_table[ BoundedNeint(index) ][type[a2]][type[a1]] + q1q2[inb]/(r*r));
                            /*!EINTCALPRINT]*/
#else                       /*EINTCALPRINT[*/
            epair = (eint_table[ BoundedNeint(index) ][type[a2]][type[a1]] + q1q2[inb]/(r*r));
#endif                      /*EINTCALPRINT]*/
                            /*BOUNDED]*/
#else                       /*!BOUNDED[*/
            assert( Ang_to_index(r) >= 0 );
            assert( Ang_to_index(r) < NEINT );
#ifndef EINTCALPRINT        /*!EINTCALPRINT[*/
            eint += (eint_table[ Ang_to_index(r) ][type[a2]][type[a1]] + q1q2[inb]/(r*r));
                            /*!EINTCALPRINT]*/
#else                       /*EINTCALPRINT[*/
            epair = (eint_table[ Ang_to_index(r) ][type[a2]][type[a1]] + q1q2[inb]/(r*r));
#endif                      /*EINTCALPRINT]*/
#endif                      /*!BOUNDED]*/

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

#ifdef BOUNDED              /*BOUNDED[*/
            index = SqAng_to_index(r2);
#ifndef EINTCALPRINT        /*!EINTCALPRINT[*/
            eint += (eint_table[ BoundedNeint(index) ][type[a2]][type[a1]] + q1q2[inb]/r2);
                            /*!EINTCALPRINT]*/
#else                       /*EINTCALPRINT[*/
            epair = (eint_table[ BoundedNeint(index) ][type[a2]][type[a1]] + q1q2[inb]/r2);
#endif                      /*EINTCALPRINT]*/
                            /*BOUNDED]*/
#else                       /*!BOUNDED[*/
            assert( SqAng_to_index(r2) >= 0 );
            assert( SqAng_to_index(r2) < NEINT );
#ifndef EINTCALPRINT        /*!EINTCALPRINT[*/
            eint += (eint_table[ SqAng_to_index(r2) ][type[a2]][type[a1]] + q1q2[inb]/r2);
                            /*!EINTCALPRINT]*/
#else                       /*EINTCALPRINT[*/
            epair = (eint_table[ SqAng_to_index(r2) ][type[a2]][type[a1]] + q1q2[inb]/r2);
#endif                      /*EINTCALPRINT]*/
#endif                      /*!BOUNDED]*/

        } else {

#ifdef BOUNDED              /*BOUNDED[*/
            index = SqAng_to_index(sqhypotenuse(dx,dy,dz));
#ifndef EINTCALPRINT        /*!EINTCALPRINT[*/
            eint += eint_table[ BoundedNeint(index) ][type[a2]][type[a1]];
                            /*!EINTCALPRINT]*/
#else                       /*EINTCALPRINT[*/
            epair = eint_table[ BoundedNeint(index) ][type[a2]][type[a1]];
#endif                      /*EINTCALPRINT]*/
                            /*BOUNDED]*/
#else                       /*!BOUNDED[*/
            assert( SqAng_to_index(sqhypotenuse(dx,dy,dz)) >= 0 );
            assert( SqAng_to_index(sqhypotenuse(dx,dy,dz)) < NEINT );
#ifndef EINTCALPRINT        /*!EINTCALPRINT[*/
            eint += eint_table[ SqAng_to_index(sqhypotenuse(dx,dy,dz)) ][type[a2]][type[a1]];
                            /*!EINTCALPRINT]*/
#else                       /*EINTCALPRINT[*/
            epair = eint_table[ SqAng_to_index(sqhypotenuse(dx,dy,dz)) ][type[a2]][type[a1]];
#endif                      /*EINTCALPRINT]*/
#endif                      /*!BOUNDED]*/

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

    return (float)eint;
}
/* EOF */
