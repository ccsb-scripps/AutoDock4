/*

 $Id: intnbtable.cc,v 1.34 2020/05/07 21:23:14 mp Exp $

 AutoDock 

Copyright (C) 2009 The Scripps Research Institute. All rights reserved.

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
#include <string.h>
#include <stdlib.h> // for exit codes
#include "intnbtable.h"
#include "structs.h"
#include "distdepdiel.h"
#include "autocomm.h"
#include "stop.h"

#ifdef NOSQRT
    /*  ACCELERATED NON-SQUARE-ROOTING VERSION  *  Look-up internal non-bond energy based on square-of-the-distance, in square Angstroms. */
#   define IndexToDistance(i) sqrt( index_to_SqAng( i ) )

#else
    /*  SQUARE-ROOTING VERSION  *  Look-up internal non-bond energy based on distance, in Angstroms.  */
#   define IndexToDistance(i) index_to_Ang( i )

#endif

static void printminvalue( char *label, ETableType *e_vdW_Hb/*[NEINT]*/[MAX_ATOM_TYPES][MAX_ATOM_TYPES],
 const int neint, const int a1, const int a2, FILE *logFile);  // see end of this file


void intnbtable( Boole *const P_B_havenbp,
                 const int a1,
                 const int a2, 
                 const GridMapSetInfo *const info,
                 ConstReal cA, 
                 ConstReal cB, 
                 const int xA, 
                 const int xB,
		 Boole is_hbond,
                 ConstReal nb_cutoff, /* soft cutoff honored by vdw & hb only */
                 ConstReal r_smooth,
                 const Linear_FE_Model AD4,
                 ConstDouble sigma,
                 /* not const */ EnergyTables *const ad_tables,
                 const Boole B_is_unbound_calculation,
		 FILE *logFile,
		 const int outlev)
{
    /* Local variables: */


    double rA;
    double rB;
    double r;
    double minus_inv_two_sigma_sqd = -0.5L / (sigma * sigma);

    register int i;

    Clock  nbeStart; // can be made static to keep lint quiet
    struct tms tms_nbeStart; // can be made static to keep lint quiet
    Clock  nbeEnd;
    struct tms tms_nbeEnd;

    char *calc_type;

    if (B_is_unbound_calculation) calc_type="unbound";
    else calc_type="internal";

    *P_B_havenbp = TRUE;
    ad_tables->is_hbond[a1][a2]  =  ad_tables->is_hbond[a2][a1]  = is_hbond;
    ad_tables->nbc[a1][a2]  =  ad_tables->nbc[a2][a1]  = nb_cutoff;

    if( outlev >= LOGETABLES ) {
    if (a1 != a2) {
        pr( logFile, "\nNon-bonded parameters for %s-%s and %s-%s interactions, used in %s energy calculations:\n", info->atom_type_name[a1], info->atom_type_name[a2], info->atom_type_name[a2], info->atom_type_name[a1], calc_type );
    } else {
        pr( logFile, "\nNon-bonded parameters for %s-%s interactions, used in %s energy calculations:\n", info->atom_type_name[a1], info->atom_type_name[a2], calc_type );
    }
    // Output the form of the potential energy equation:
    if ( B_is_unbound_calculation ) {
        pr( logFile, "\n            %12.5lg\n", cA );
        pr( logFile, "    E      =  -----------  -  r\n");
        pr( logFile, "     %2s,%-2s         %2d\n", info->atom_type_name[a1], info->atom_type_name[a2], xA );
        pr( logFile, "                  r\n\n");
    } else {
	Real weight = is_hbond ? AD4.coeff_hbond: AD4.coeff_vdW;
	const char * wname = is_hbond ? "coeff_hbond" : "coeff_vdW";
	Real cA_unw = cA/weight;
	Real cB_unw = cB/weight;
        pr( logFile, "\n                            %12.5lg    %12.5lg    %12.5lg   %12.5lg \n", cA_unw, cB_unw, cA, cB );
	 
	pr( logFile, "    E      =  %11s * (-----------  -  -----------) =  -----------  -  -----------\n",
	wname);
        pr( logFile, "     %2s,%-2s        %6.4f           %-2d               %-2d             %-2d            %-2d\n", info->atom_type_name[a1], info->atom_type_name[a2], weight, xA, xB , xA,xB);
        pr( logFile, "                                  r                r              r             r \n\n");
//      pr( logFile, "\n            %12.5lg   %12.5lg \n", cA, cB );
//      pr( logFile, "    E      =  -----------  -  -----------\n");
//      pr( logFile, "     %2s,%-2s         %2d              %2d\n", info->atom_type_name[a1], info->atom_type_name[a2], xA, xB );
//      pr( logFile, "                  r               r \n\n");
    }
    pr( logFile, "\n Distance cutoff for %s evaluation = %.2f\n", 
	is_hbond ? "hbond" : "vdW", nb_cutoff);
    pr( logFile, "\n Smoothing range  = %.2f\n", r_smooth);
    if(r_smooth==0) pr( logFile, "\n No smoothing will be applied.\n");
    pr( logFile, "Calculating %s-%-s interaction energy versus atomic separation (%.3f to %.3f Ang, %d intervals: widths %.3f to %.3f).\n",
    info->atom_type_name[a1], info->atom_type_name[a2], 
    IndexToDistance(1), IndexToDistance(NEINT-1),
    NEINT, 
    IndexToDistance(2)-IndexToDistance(1),
    IndexToDistance(NEINT-1)-IndexToDistance(NEINT-2) );
    flushLog;

    } // end if outlev
    nbeStart = times( &tms_nbeStart );

  if ( outlev >= LOGETABLES ) {
 // MPique 2011 produce table of distances and values for judging cutoffs
    fprintf(logFile, "DISTANCES SQA_DIV= %.5f\n", SQA_DIV);
    fprintf(logFile, "DISTANCES INV_SQA_DIV= %.5f\n", INV_SQA_DIV);
    fprintf(logFile, "DISTANCES %s %.2f %.2f %.2f  ... %.2f %.2f %.2f [%d]\n",
       "NEINT 1,2,3,...,NEINT-1",
       IndexToDistance(1),
       IndexToDistance(2),
       IndexToDistance(3),
       IndexToDistance(NEINT-3),
       IndexToDistance(NEINT-2),
       IndexToDistance(NEINT-1),
       NEINT-1);
    }
    // Allocate the e_vdW_Hb [a1][a2] 1-D array if needed, 
    //  put pointer to it in [a2][a1] as well;
    //  TODO check if not empty upon first call here
    if ( ad_tables->e_vdW_Hb[a1][a2] == NULL ) {
	    ad_tables->e_vdW_Hb[a1][a2] = (ETableType *) calloc( NEINT, sizeof (ETableType));
	    ad_tables->e_vdW_Hb[a2][a1] = ad_tables->e_vdW_Hb[a1][a2];
	}
    if (ad_tables->e_vdW_Hb[a1][a2]==NULL) stop ( "unable to allocate e_vdW_Hb table");
	    

    // loop from index 1 up to a maximum distance of  (NEINT-1) * INV_A_DIV, 
    //                          usually    2048 * 0.01,       or 20.48 Angstroms
    // Note the zero-th entry is set to EINTCLAMP

    ad_tables->e_vdW_Hb[a1][a2][0]  =  ad_tables->e_vdW_Hb[a2][a1][0]  =   EINTCLAMP;

    for ( i = 0;  i < NEINT;  i++ ) {
        // i is the lookup-table index that corresponds to the distance

        // r is the distance that corresponds to the lookup index
        r = IndexToDistance(i); 

        // Compute r^xA and r^xB:
        rA = pow( r, (double)xA );
        rB = pow( r, (double)xB );

        if ( B_is_unbound_calculation ) {
            // Calculate the unbound potential for computing the 
            // unbound extended conformation of the ligand:
            // E = -|r|
            // ad_tables->e_vdW_Hb[a1][a2][i]  =  ad_tables->e_vdW_Hb[a2][a1][i]  = -1. * fabs( r );
            // Calculate the interaction energy at this distance, r, using an equation 
            // of the form E  =  cA / r^xA  i.e. just the repulsive term
            // minus r, to make the potential long range
            if(i>0) ad_tables->e_vdW_Hb[a1][a2][i]  =  ad_tables->e_vdW_Hb[a2][a1][i]  =  min( EINTCLAMP, (cA/rA) ) - r;
        } else {
            // Calculate the bound potential for docking:

            // Calculate the interaction energy at this distance, r, using an equation 
            // of the form E  =  cA / r^xA  -  cB / r^xB
            if(i>0) ad_tables->e_vdW_Hb[a1][a2][i]  =  ad_tables->e_vdW_Hb[a2][a1][i]  =  min( EINTCLAMP, (cA/rA - cB/rB) );

        }
    // optionally dump non-zero-ish values to log file for scoring function development
    if( outlev >= LOGETABLES && (!B_is_unbound_calculation) && (i==0 || fabs(ad_tables->e_vdW_Hb[a1][a2][i]) >= 0.001 )) {
                pr( logFile, "i=%6d  ad_tables->e_vdW_Hb = %12.8g,   r=%7.4lf\n",i, ad_tables->e_vdW_Hb[a1][a2][i], r );
            }

    } // next i // for ( i = 1;  i < NEINT;  i++ )
  
    //ad_tables->e_vdW_Hb[a1][a2][NEINT-1]  =  ad_tables->e_vdW_Hb[a2][a1][NEINT-1] = 0;


    // report range of minimum values before smoothing
    if( outlev >= LOGETABLES ) {
       char label[40];
       sprintf(label, "before smoothing %s %s-%-s ", calc_type, info->atom_type_name[a1], info->atom_type_name[a2]);
       printminvalue( label, ad_tables->e_vdW_Hb, NEINT, a1, a2, logFile);
    }

    /* smooth with min function; 
      r_smooth is Angstrom range of "smoothing" */
    if (r_smooth > 0) {
        Real energy_smooth[NEINT];
        for (i = 0;  i <= NEINT-1;  i++) {
	    double r = IndexToDistance(i);
	    double rlow  = r - r_smooth/2;
	    double rhigh = r + r_smooth/2;
            energy_smooth[i] = 100000.;
#ifdef NOSQRT
            for (int j = max(0, BoundedSqAng_to_index(rlow*rlow));  \
	      j <= min(NEINT-1, BoundedSqAng_to_index(rhigh*rhigh));  j++)
#else
            for (int j = max(0, BoundedAng_to_index(rlow)); \
	      j <= min(NEINT-1, BoundedAng_to_index(rhigh));  j++) 
#endif
              energy_smooth[i] = min(energy_smooth[i], ad_tables->e_vdW_Hb[a1][a2][j]);
        }
        for (i = 1;  i < NEINT-1;  i++) {
            ad_tables->e_vdW_Hb[a1][a2][i]  =  ad_tables->e_vdW_Hb[a2][a1][i] = energy_smooth[i];
        }
       // report range of minimum values after smoothing
       if( outlev >= LOGETABLES ) {
          char label[40];
          sprintf(label, "after  smoothing %s %s-%-s ", calc_type, info->atom_type_name[a1], info->atom_type_name[a2]);
          printminvalue( label, ad_tables->e_vdW_Hb, NEINT, a1, a2, logFile);
       }
    } /* endif smoothing */

    if ( outlev >= LOGETABLES ) {
    // MPique fall 2011
    fprintf(logFile, "DISTANCES %s %.3f %.3f %.3f %.3f  ... %.3f %.3f %.3f [%d]\n",
       "NDIEL 0,1,2,3,...NDIEL-1",
       IndexToDistance(0),
       IndexToDistance(1),
       IndexToDistance(2),
       IndexToDistance(3),
       IndexToDistance(NDIEL-3),
       IndexToDistance(NDIEL-2),
       IndexToDistance(NDIEL-1),
       NDIEL-1);
    }

    // loop from 0 up to a maximum distance of  (NDIEL-1) * INV_A_DIV, 
    //                          usually    16384 * 0.01,       or 163.84 Angstroms
    for ( i = 0;  i < NDIEL;  i++ ) {
        // i is the lookup-table index that corresponds to the distance

        // r is the distance that corresponds to the lookup index
        r = IndexToDistance(i); 

        // Compute the distance-dependent gaussian component of the desolvation energy, sol_fn[i];
        // Weight this by the coefficient for desolvation, AD4.coeff_desolv.
        ad_tables->sol_fn[i] = AD4.coeff_desolv * exp( minus_inv_two_sigma_sqd * sq(r) );
    } // next i // for ( i = 1;  i < NDIEL;  i++ )
    
    if( outlev >= LOGETABLES ) {
	    nbeEnd = times( &tms_nbeEnd );
	    pr( logFile, "Time taken: ");
	    timesys( nbeEnd - nbeStart, &tms_nbeStart, &tms_nbeEnd, logFile);
    }
}
/* end of intnbtable */


void setup_distdepdiel( FILE *logFile,
			const int outlev, 
                        EnergyTables *const ptr_ad_energy_tables  // Holds vdw+Hb, desolvation & dielectric lookup tables
                      )
{
    register int i=0;
    register double distance=0.0L;

    if (outlev >= LOGETABLES ) {
        pr(logFile, "Calculating distance-dependent dielectric function using the method of Mehler & Solmajer\n\n\n");
    }

    ptr_ad_energy_tables->epsilon_fn[0] = 1.0L;
    if (outlev >= LOGETABLES) {
        pr(logFile, "i, ptr_ad_energy_tables->epsilon_fn[i] = %d, %8.4lf\n", i, ptr_ad_energy_tables->epsilon_fn[i]);
    }
    for (i = 1;  i < NDIEL;  i++) {
        distance = IndexToDistance(i);
        ptr_ad_energy_tables->epsilon_fn[i] = calc_ddd_Mehler_Solmajer( distance, APPROX_ZERO );
        ptr_ad_energy_tables->r_epsilon_fn[i] = distance * calc_ddd_Mehler_Solmajer( distance, APPROX_ZERO );
        if (outlev >= LOGETABLES &&
             ((i<1000&&(i%100)==0) || (i%1000) == 0)) {
                pr(logFile, "i = %5d,  distance = %7.3lf,  epsilon_fn[i] = %9.4lf,  r_epsilon_fn[i] = %9.4lf\n", 
                        i, distance, ptr_ad_energy_tables->epsilon_fn[i], ptr_ad_energy_tables->r_epsilon_fn[i]);
            }
        // pre-compute reciprocal to avoid having to do it later in eintcal.
        ptr_ad_energy_tables->r_epsilon_fn[i] = 1.0 /  ptr_ad_energy_tables->r_epsilon_fn[i];
    } // next i
}
/* end of setup_distdepdiel */

static void printminvalue( char *label, ETableType *e_vdW_Hb/*[NEINT]*/[MAX_ATOM_TYPES][MAX_ATOM_TYPES],
 const int neint, const int a1, const int a2, FILE *logFile) 
{
   /* report distance range over which minimum value appears - a 'smoothing' debug tool mostly */
int i;
double minval;
int ifirstmin=-1, ilastmin=-1; // indices

	// find minimum value in list
	minval=EINTCLAMP; // big
	for(i=0;i<neint;i++) minval=min(minval, e_vdW_Hb[a1][a2][i]);

	// find first appearance of minval
	for(i=0;i<neint;i++) if(e_vdW_Hb[a1][a2][i]==minval) break;
	ifirstmin=i;

	// find last appearance of minval
	for(i=0;i<neint;i++) if(e_vdW_Hb[a1][a2][i]==minval) ilastmin=i;

	if(ifirstmin<0||ifirstmin>=neint||ilastmin<0||ilastmin>=neint) {
		// bug check
		fprintf(logFile, " bug check in intnbtable:printminvalue\n");
		exit(EXIT_FAILURE);
		}
	fprintf(logFile, "%s minval=%10.6f over range d= %8.5f to %8.5f\n",
	   label, minval, IndexToDistance(ifirstmin), IndexToDistance(ilastmin) );
}


/* EOF */
