/*

 $Id: intnbtable.cc,v 1.3 2004/02/12 04:32:15 garrett Exp $

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

#endif

extern FILE *logFile;


void intnbtable(Boole *P_B_havenbp,
                int *P_a1,
                int *P_a2, 
                int num_atm_maps,
                char atm_tyP_str[ATOM_MAPS],
                FloatOrDouble cA, 
                FloatOrDouble cB, 
                int xA, 
                int xB,
                FloatOrDouble e_internal[NEINT][ATOM_MAPS][ATOM_MAPS] ) 
{ 
    /* Local variables: */

    Clock  nbeEnd;
    Clock  nbeStart;

    double dxA;
    double dxB;
    double rA;
    double rB;
    double r;

    int Ais2B = FALSE;
    int Bis2A = FALSE;
    register int i;

    struct tms tms_nbeEnd;
    struct tms tms_nbeStart;

    *P_B_havenbp = TRUE;

    pr( logFile, "Non-bonded parameters for %c-%c interactions, used in internal energy calculations:\n", atm_tyP_str[*P_a1], atm_tyP_str[*P_a2] );
    pr( logFile, "\n             %9.1lf       %9.1lf \n", cA, cB );
    pr( logFile, "    E    =  -----------  -  -----------\n");
    pr( logFile, "     %c,%c         %2d              %2d\n", atm_tyP_str[*P_a1], atm_tyP_str[*P_a2], xA, xB );
    pr( logFile, "                r               r \n\n");
    pr( logFile, "Calculating %c-%c interaction energy versus atomic separation (%d data points).\n", atm_tyP_str[*P_a1], atm_tyP_str[*P_a2], NEINT );
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
    for ( i = 1;  i < NEINT;  i++ ) {
        r = LookUpProc( i );
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
        e_internal[i][*P_a1][*P_a2] = e_internal[i][*P_a2][*P_a1] 
                        = min( EINTCLAMP, (cA/rA - cB/rB) );
        // pr( logFile, "i=%6d  e_internal = %.3f,   r=%.4lf\n",i, e_internal[i][*P_a1][*P_a2], r ); // Xcode-gmm
    }
    nbeEnd = times( &tms_nbeEnd );
    pr( logFile, "Time taken: ");
    timesys( nbeEnd - nbeStart, &tms_nbeStart, &tms_nbeEnd );
    if (++(*P_a2) >= num_atm_maps) {
        *P_a2 = ++(*P_a1);
    }
}
/* EOF of intnbtable.c */
