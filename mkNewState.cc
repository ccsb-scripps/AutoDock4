/*

 $Id: mkNewState.cc,v 1.10 2007/04/27 06:01:49 garrett Exp $

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
#include <stdlib.h>
#include "mkNewState.h"


void mkNewState( State *now,
                 State *last,        /* ...must be a normalized quaternion! */
                 State *change,

                Real vt[MAX_TORS][NTRN],
                int   tlist[MAX_TORS][MAX_ATOMS],
                int   ntor,
                Real crd[MAX_ATOMS][NTRN],
                Real crdpdb[MAX_ATOMS][NTRN],
                int   natom,
                Real trnStep,
                Real qtwStep,
                Real torStep,
                Real F_TorConRange[MAX_TORS][MAX_TOR_CON][2],
                int N_con[MAX_TORS])
    // Create a new state, based on the current state
{
    register int i;
    double t;
    int I_ranCon;
    Real a, b;

    /*
    ** Translation
    */
    now->T.x = last->T.x + (change->T.x = random_pm( trnStep ));
    now->T.y = last->T.y + (change->T.y = random_pm( trnStep ));
    now->T.z = last->T.z + (change->T.z = random_pm( trnStep ));

    /*
    ** Quaternion angular displacement
    */
    if (qtwStep > APPROX_ZERO) {
        /*
        // (This code probably does *not* produce a uniformly distributed quaternion)
        change->Q.nx  = Randpm1; 
        change->Q.ny  = Randpm1; 
        change->Q.nz  = Randpm1; 
        change->Q.ang = random_pm( qtwStep );
        mkUnitQuat( &(change->Q) );
        */

        /*
        **  This should produce a uniformly distributed quaternion, according to
        **  Shoemake, Graphics Gems III.6, pp.124-132, "Uniform Random Rotations",
        **  published by Academic Press, Inc., (1992)
        */
        change->Q = uniformQuat();

        /*
        **  Apply random change, to Last Quaternion
        */
        qmultiply( &(now->Q), &(last->Q), &(change->Q) );
    }

    for (i=0; i<ntor; i++) {
        if (N_con[i] > 0) {
            if (N_con[i] > 1) {
                /* If N_con was 2, I_ranCon could be 0 or 1, never 2
                 * Select a random constraint */
                I_ranCon = (int)((double)N_con[i] * local_random());  
            } else {
                /* Hobson's choice...  */
                I_ranCon = 0;
            }
            a = F_TorConRange[i][I_ranCon][LOWER];
            b = F_TorConRange[i][I_ranCon][UPPER];
            t = random_range(a,b);
            now->tor[i]    = WrpModRad(t);
            change->tor[i] = now->tor[i] - last->tor[i];
        } else {
            t = last->tor[i] + (change->tor[i] = random_pm( torStep ));
            now->tor[i]    = WrpModRad(t);
        }
    }/*i*/

    cnv_state_to_coords( *now,  vt, tlist, ntor,  crdpdb, crd, natom );
} 
/* EOF */
