/*

 $Id: intnbtable.cc,v 1.4 2005/03/11 02:11:30 garrett Exp $

*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/* intnbtable.cc */

#include <math.h>
#include <stdio.h>
#include <time.h>
#include <sys/types.h>
#include <sys/times.h>
#include "intnbtable.h"

#ifdef NOSQRT

/*  ACCELERATED NON SQUARE-ROOTING VERSION;
 *  Look-up internal non-bond energy based on square-of-the-distance,
 *  in square Angstroms. This saves a square-root operation for each
 *  non-bonded interaction.
 */

#define        LookUpProc(i)        sqrt( index_to_SqAng( i ) )

#else

/*  SQUARE-ROOTING VERSION;
 *  Look-up internal non-bond energy based on distance,
 *  in Angstroms.
 */

#define        LookUpProc(i)        index_to_Ang( i )

// #define index_to_Ang(i)         (( (FloatOrDouble) (i) ) * INV_A_DIV)
// constants.h:#define INV_A_DIV   0.01      // every 1/100-th of an Angstrom

#endif

extern FILE *logFile;
extern int debug;

void intnbtable( Boole *P_B_havenbp,
                 int *P_a1,
                 int *P_a2, 
                 int num_atm_maps,
                 char ligand_atom_types[MAX_MAPS][3],
                 FloatOrDouble cA, 
                 FloatOrDouble cB, 
                 int xA, 
                 int xB,
                 FloatOrDouble e_internal[NEINT][ATOM_MAPS][ATOM_MAPS],
                 FloatOrDouble sol_fn[NEINT],
                 double coeff_desolv,
                 double sigma)
{
    /* Local variables: */

    Clock  nbeEnd;
    Clock  nbeStart;

    double dxA;
    double dxB;
    double rA;
    double rB;
    double r;
    double minus_inv_two_sigma_sqd = -0.5L / (sigma * sigma);

    int Ais2B = FALSE;
    int Bis2A = FALSE;
    register int i;

    struct tms tms_nbeEnd;
    struct tms tms_nbeStart;

    *P_B_havenbp = TRUE;

    if (*P_a1 != *P_a2) {
        pr( logFile, "\nNon-bonded parameters for %s-%s and %s-%s interactions, used in internal energy calculations:\n", ligand_atom_types[*P_a1], ligand_atom_types[*P_a2], ligand_atom_types[*P_a2], ligand_atom_types[*P_a1] );
    } else {
        pr( logFile, "\nNon-bonded parameters for %s-%s interactions, used in internal energy calculations:\n", ligand_atom_types[*P_a1], ligand_atom_types[*P_a2] );
    }
    pr( logFile, "\n               %9.1lf       %9.1lf \n", cA, cB );
    pr( logFile, "    E      =  -----------  -  -----------\n");
    pr( logFile, "     %2s,%-2s         %2d              %2d\n", ligand_atom_types[*P_a1], ligand_atom_types[*P_a2], xA, xB );
    pr( logFile, "                  r               r \n\n");
    pr( logFile, "Calculating %s-%-s interaction energy versus atomic separation (%d data points).\n", ligand_atom_types[*P_a1], ligand_atom_types[*P_a2], NEINT );
    flushLog;
    nbeStart = times( &tms_nbeStart );
    dxA = (double) xA;
    dxB = (double) xB;
    if (xB == (2*xA)) {
        Bis2A = TRUE;
        Ais2B = FALSE;
    } else if (xA == (2*xB)) {
        Ais2B = TRUE;
        Bis2A = FALSE;
    } else {
        Ais2B = FALSE;
        Bis2A = FALSE;
    }
    // loop up to a maximum distance of  (NEINT * INV_A_DIV), 
    //                          usually    2048 * 0.01,       or 20.48 Angstroms

    for ( i = 1;  i < NEINT;  i++ ) {
        // i is essentially the distance (in integer form)

        r = LookUpProc(i); // r is the distance

        // Compute the distance-dependent gaussian component of the desolvation energy, sol_fn
        // Weight this by the coefficient for desolvation
        sol_fn[i] = coeff_desolv * exp( minus_inv_two_sigma_sqd * sq(r) );

        if (Bis2A) {
            rA = pow( r, dxA );
            rB = rA * rA;
        } else if (Ais2B) {
            rB = pow( r, dxB );
            rA = rB * rB;
        } else {
            rA = pow( r, dxA );
            rB = pow( r, dxB );
        }

        e_internal[i][*P_a1][*P_a2] = e_internal[i][*P_a2][*P_a1] = min( EINTCLAMP, (cA/rA - cB/rB) );

        if (debug > 1) {
            pr( logFile, "i=%6d  e_internal = %.3f,   r=%.4lf\n",i, e_internal[i][*P_a1][*P_a2], r ); // Xcode-gmm
        }

    } // 1 <= i < NEINT
    
    nbeEnd = times( &tms_nbeEnd );
    pr( logFile, "Time taken: ");
    timesys( nbeEnd - nbeStart, &tms_nbeStart, &tms_nbeEnd );

    if (++(*P_a2) >= num_atm_maps) {
        *P_a2 = ++(*P_a1);
    }
}
/* EOF of intnbtable.c */
