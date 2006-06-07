/*

 $Id: nbe.cc,v 1.5 2006/06/07 00:06:26 garrett Exp $

*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/* nbe.cc */

#include <math.h>

    #include <stdio.h>
    #include "nbe.h"


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

void nbe( GridMapSetInfo *info,
          EnergyTables *ptr_ad_energy_tables,
          int num_atm_maps )

{
 
    static int NUMPTS = 640;
    register int i = 0;
    register int j = 0;
    register int k = 0;
    Real r = 0.;

    pr( logFile,"SUMMARY OF PAIRWISE-ATOMIC NON-BONDED INTERNAL ENERGIES\n" );
    pr( logFile,"________________________________________________________\n\n");
    pr( logFile,"Clamp pairwise-atomic interaction energies at:\t%.2f\n", EINTCLAMP );
 
    pr( logFile, "    \t\n r  \tLook-up\t" );
    for ( i = 0; i < num_atm_maps; i++) {
        for ( j = i; j < num_atm_maps; j++) {
            pr( logFile, "  E    " );
        }
    }
    pr( logFile, "\n /Ang\tIndex\t" );
    for ( i = 0; i < num_atm_maps; i++) {
        for ( j = i; j < num_atm_maps; j++) {
            pr( logFile, "   %2s,%-2s ", info->atom_type_name[i], info->atom_type_name[j] );
        }
    }
    pr( logFile, "\n______\t_____\t" );
    for ( i = 0; i < num_atm_maps; i++) {
        for ( j = i; j < num_atm_maps; j++) {
            pr( logFile, " ______" );
        }
    }
    pr( logFile, "\n" );
    for ( k = 10;  k <= NUMPTS;  k += 10 ) {
        r = LookUpProc( k );
        pr( logFile, "%6.3f\t%5d\t", r, k );
        for ( i = 0;  i < num_atm_maps; i++) {
            for ( j = i;  j < num_atm_maps; j++) {
                pr( logFile, "%7.2f", ptr_ad_energy_tables->e_vdW_Hb[k][j][i] );
            } /*  j  */
        } /*  i  */
        pr( logFile, "\n" );
    } /*  k  */
    flushLog;
}
/* EOF */
