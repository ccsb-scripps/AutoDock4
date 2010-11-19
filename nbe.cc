/*

 $Id: nbe.cc,v 1.6 2007/04/27 06:01:50 garrett Exp $

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
