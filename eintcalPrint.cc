/* eintcalPrint.cc */

#include <math.h>

#include <stdio.h>
#include "eintcalPrint.h"


extern FILE *logFile;

Real eintcalPrint( NonbondParam * nonbondlist,
                            Real eint_table[NEINT][ATOM_MAPS][ATOM_MAPS],
                            Real tcoord[MAX_ATOMS][SPACE],
                            int type[MAX_ATOMS],
                            int Nnb,
                            Boole B_calcIntElec,
                            Real abs_charge[MAX_ATOMS]
                            )

/******************************************************************************/
/*      Name: eintcal                                                         */
/*  Function: Calculate the Internal Energy of the Small Molecule.            */
/*            Accelerated non-square-rooting, dx,dy,dz version.               */
/* Copyright: (C) 1994, TSRI                                                  */
/*____________________________________________________________________________*/
/*    Author: Garrett M. Morris, TSRI                                         */
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
{
    double eint = 0., epair=0., dx,dy,dz,r2;
#ifndef NOSQRT
    double d;
#endif /* NOSQRT */
    int a1, a2;
    register int inb;

    pr( logFile, "Non-bond  Atom1-Atom2  Distance  Energy\n");

    for (inb = 0;  inb < Nnb;  inb++) {

	dx = tcoord[(a1 = nonbondlist[inb].a1)][X] - tcoord[(a2 = nonbondlist[inb].a2)][X];
	dy = tcoord[a1][Y] - tcoord[a2][Y];
	dz = tcoord[a1][Z] - tcoord[a2][Z];

	r2 = sqhypotenuse(dx,dy,dz);

#ifndef NOSQRT 
	/*
	** NOSQRT not defined, i.e. sqrt-ing is OK,
	*/
	d = sqrt( r2 );
        if (B_calcIntElec) {
            epair = eint_table[ Ang_to_index(d) ][type[a2]][type[a1]] + nonbondlist[inb].q1q2/r2;
        } else {
	    epair = eint_table[ Ang_to_index(d) ][type[a2]][type[a1]];
        }
#else
	/*
	**  NOSQRT *was* defined,
	*/
        if (B_calcIntElec) {
            epair = eint_table[ SqAng_to_index(r2) ][type[a2]][type[a1]] + nonbondlist[inb].q1q2/r2;
        } else {
	    epair = eint_table[ SqAng_to_index(r2) ][type[a2]][type[a1]];
        }
#endif

        eint += epair;
        pr( logFile, " %6d   %5d-%-5d  %7.2lf  %+8.3lf\n", (int)(inb+1), (int)(a1+1), 
	(int)(a2+1), (double)sqrt(r2), (double)epair);

    }

    pr( logFile, "\n\n Total internal non-bonded energy = %+8.3lf\n", (double)eint);

    return (Real)eint;
}
/* EOF */
