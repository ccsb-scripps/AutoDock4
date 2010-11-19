/*

 $Id: intnbtable.cc,v 1.10 2007/04/27 06:01:49 garrett Exp $

 AutoDock 

 Copyright (C) 1989-2007,  Garrett M. Morris, David S. Goodsell, Ruth Huey, Arthur J. Olson, 
 All Rights Reserved.

 AutoDock is a Trade Mark of The Scripps Research Institute.

 This program is free software; you can redistribute it and/or
 modify it under the terms of the GNU General Public License
 as published by the Free Software Foundation; either version 2
 of the License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

 */


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>
#include <stdio.h>
#include <time.h>
#include <sys/types.h>
#include <sys/times.h>
#include <string.h>
#include "intnbtable.h"
#include "structs.h"
#include "distdepdiel.h"
#include "autocomm.h"

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
                 Real cA, 
                 Real cB, 
                 int xA, 
                 int xB,
                 double coeff_desolv,
                 double sigma,
                 EnergyTables *ad_tables,
                 Boole B_is_unbound_calculation )
{
    /* Local variables: */

    Clock  nbeEnd;
    Clock  nbeStart;

    double rA;
    double rB;
    double r;
    double minus_inv_two_sigma_sqd = -0.5L / (sigma * sigma);

    register int i;

    struct tms tms_nbeEnd;
    struct tms tms_nbeStart;

    char calc_type[128];

    if (B_is_unbound_calculation == TRUE) {
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
    // Output the form of the potential energy equation:
    if ( B_is_unbound_calculation ) {
        pr( logFile, "\n               %9.1lf\n", cA );
        pr( logFile, "    E      =  -----------  -  r\n");
        pr( logFile, "     %2s,%-2s         %2d\n", info->atom_type_name[a1], info->atom_type_name[a2], xA );
        pr( logFile, "                  r\n\n");
    } else {
        pr( logFile, "\n               %9.1lf       %9.1lf \n", cA, cB );
        pr( logFile, "    E      =  -----------  -  -----------\n");
        pr( logFile, "     %2s,%-2s         %2d              %2d\n", info->atom_type_name[a1], info->atom_type_name[a2], xA, xB );
        pr( logFile, "                  r               r \n\n");
    }
    pr( logFile, "Calculating %s-%-s interaction energy versus atomic separation (%d data points).\n", info->atom_type_name[a1], info->atom_type_name[a2], NEINT );
    flushLog;

    nbeStart = times( &tms_nbeStart );

    // loop up to a maximum distance of  (NEINT * INV_A_DIV), 
    //                          usually    2048 * 0.01,       or 20.48 Angstroms

    for ( i = 1;  i < NEINT;  i++ ) {
        // i is the lookup-table index that corresponds to the distance

        // r is the distance that corresponds to the lookup index
        r = IndexToDistance(i); 

        // Compute the distance-dependent gaussian component of the desolvation energy, sol_fn[i];
        // Weight this by the coefficient for desolvation, coeff_desolv.
        ad_tables->sol_fn[i] = coeff_desolv * exp( minus_inv_two_sigma_sqd * sq(r) );

        // Compute r^xA and r^xB:
        rA = pow( r, (double)xA );
        rB = pow( r, (double)xB );

        if ( B_is_unbound_calculation ) {
            // Calculate the unbound potential for computing the 
            // unbound extended conformation of the ligand:
            // E = -|r|
            // ad_tables->e_vdW_Hb[i][a1][a2]  =  ad_tables->e_vdW_Hb[i][a2][a1]  = -1. * fabs( r );
            // Calculate the interaction energy at this distance, r, using an equation 
            // of the form E  =  cA / r^xA  i.e. just the repulsive term
            // minus r, to make the potential long range
            ad_tables->e_vdW_Hb[i][a1][a2]  =  ad_tables->e_vdW_Hb[i][a2][a1]  =  min( EINTCLAMP, (cA/rA) ) - r;
        } else {
            // Calculate the bound potential for docking:

            // Calculate the interaction energy at this distance, r, using an equation 
            // of the form E  =  cA / r^xA  -  cB / r^xB
            ad_tables->e_vdW_Hb[i][a1][a2]  =  ad_tables->e_vdW_Hb[i][a2][a1]  =  min( EINTCLAMP, (cA/rA - cB/rB) );

            if (debug > 1) {
                pr( logFile, "i=%6d  ad_tables->e_vdW_Hb = %.3f,   r=%.4lf\n",i, ad_tables->e_vdW_Hb[i][a1][a2], r ); // Xcode-gmm
            }
        }

    } // next i // for ( i = 1;  i < NEINT;  i++ )
    
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
            if (i%1000 == 0) {
                pr(logFile, "i = %5d,  distance = %7.2lf,  epsilon_fn[i] = %8.4lf,  r_epsilon_fn[i] = %8.4lf\n", 
                        i, distance, ptr_ad_energy_tables->epsilon_fn[i], ptr_ad_energy_tables->r_epsilon_fn[i]);
            }
        }
        // pre-compute reciprocal to avoid having to do it later in eintcal.
        ptr_ad_energy_tables->r_epsilon_fn[i] = 1.0 /  ptr_ad_energy_tables->r_epsilon_fn[i];
    } // next i
}
/* end of setup_distdepdiel */

/* EOF */
