/*

 $Id: intnbtable.cc,v 1.6 2006/04/17 05:41:04 garrett Exp $

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
#include <string.h>
#include "intnbtable.h"
#include "structs.h"
#include "distdepdiel.h"

#ifdef NOSQRT
    /*  ACCELERATED NON-SQUARE-ROOTING VERSION  *  Look-up internal non-bond energy based on square-of-the-distance, in square Angstroms. */
#   define IndexToDistance(i) sqrt( index_to_SqAng( i ) )

#else
    /*  SQUARE-ROOTING VERSION  *  Look-up internal non-bond energy based on distance, in Angstroms.  */
#   define IndexToDistance(i) index_to_Ang( i )

#endif

extern FILE *logFile;
extern int debug;

void intnbtable( Boole *P_B_havenbp,
                 int a1,
                 int a2, 
                 GridMapSetInfo *info,
                 FloatOrDouble cA, 
                 FloatOrDouble cB, 
                 int xA, 
                 int xB,
                 double coeff_desolv,
                 double sigma,
                 EnergyTables *ad_tables)
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

    char calc_type[128];

    if ((xA == 1) && (xB == 2)) {
        strcpy(calc_type, "unbound");
    } else {
        strcpy(calc_type, "internal");
    }

    *P_B_havenbp = TRUE;

    if (a1 != a2) {
        pr( logFile, "\nNon-bonded parameters for %s-%s and %s-%s interactions, used in %s energy calculations:\n", info->atom_type_name[a1], info->atom_type_name[a2], info->atom_type_name[a2], info->atom_type_name[a1], calc_type );
    } else {
        pr( logFile, "\nNon-bonded parameters for %s-%s interactions, used in %s energy calculations:\n", info->atom_type_name[a1], info->atom_type_name[a2], calc_type );
    }
    pr( logFile, "\n               %9.1lf       %9.1lf \n", cA, cB );
    pr( logFile, "    E      =  -----------  -  -----------\n");
    pr( logFile, "     %2s,%-2s         %2d              %2d\n", info->atom_type_name[a1], info->atom_type_name[a2], xA, xB );
    pr( logFile, "                  r               r \n\n");
    pr( logFile, "Calculating %s-%-s interaction energy versus atomic separation (%d data points).\n", info->atom_type_name[a1], info->atom_type_name[a2], NEINT );
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
        // i is the lookup-table index that corresponds to the distance

        r = IndexToDistance(i); // r is the distance that corresponds to index

        // Compute the distance-dependent gaussian component of the desolvation energy, sol_fn[i];
        // Weight this by the coefficient for desolvation, coeff_desolv.

        ad_tables->sol_fn[i] = coeff_desolv * exp( minus_inv_two_sigma_sqd * sq(r) );

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

        ad_tables->e_vdW_Hb[i][a1][a2]  =  ad_tables->e_vdW_Hb[i][a2][a1]  =  min( EINTCLAMP, (cA/rA - cB/rB) );

        if (debug > 1) {
            pr( logFile, "i=%6d  ad_tables->e_vdW_Hb = %.3f,   r=%.4lf\n",i, ad_tables->e_vdW_Hb[i][a1][a2], r ); // Xcode-gmm
        }

    } // 1 <= i < NEINT
    
    nbeEnd = times( &tms_nbeEnd );
    pr( logFile, "Time taken: ");
    timesys( nbeEnd - nbeStart, &tms_nbeStart, &tms_nbeEnd );
}
/* end of intnbtable */


void setup_distdepdiel( int outlev, 
                        EnergyTables *ptr_ad_energy_tables  // Holds vdw+Hb, desolvation & dielectric lookup tables
                      )
{
    register int i=0;
    register double distance=0.0L;

    if (outlev > 0) {
        pr(logFile, "Calculating distance-dependent dielectric function using the method of Mehler & Solmajer\n\n\n");
    }

    ptr_ad_energy_tables->epsilon_fn[0] = 1.0L;
    if (outlev > 1) {
        pr(logFile, "i, ptr_ad_energy_tables->epsilon_fn[i] = %d, %8.4lf\n", i, ptr_ad_energy_tables->epsilon_fn[i]);
    }
    for (i = 1;  i < NDIEL;  i++) {
        distance = IndexToDistance(i);
        ptr_ad_energy_tables->epsilon_fn[i] = calc_ddd_Mehler_Solmajer( distance, APPROX_ZERO );
        ptr_ad_energy_tables->r_epsilon_fn[i] = distance * calc_ddd_Mehler_Solmajer( distance, APPROX_ZERO );
        if (outlev > 1) {
            if (i%10 == 0) {
                pr(logFile, "i, ptr_ad_energy_tables->epsilon_fn[i] = %d, %8.4lf\n", i, ptr_ad_energy_tables->epsilon_fn[i]);
            }
        }
    } // next i
}
/* end of setup_distdepdiel */

/* EOF */
