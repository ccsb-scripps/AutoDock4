/*

 $Id: nbe.cc,v 1.3.6.1 2005/10/10 16:48:02 alther Exp $

*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/* nbe.cc */

#ifdef __INTEL_COMPILER
   #include <mathimf.h>
#else
   #include <math.h>
#endif

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

void nbe( char atm_typ_str[ATOM_MAPS],
          EnergyTables *ptr_ad_energy_tables,
          int num_atm_maps )

{

    static int NUMPTS = 640;
    register int i = 0;
    register int j = 0;
    register int k = 0;
    FloatOrDouble r = 0.;

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
	    pr( logFile, "   %c,%c ", atm_typ_str[i], atm_typ_str[j] );
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
