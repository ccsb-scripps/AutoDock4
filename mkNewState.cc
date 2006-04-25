/*

 $Id: mkNewState.cc,v 1.3 2006/04/25 22:32:32 garrett Exp $

*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/* mkNewState.cc */

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
                /*Real qtwStep,*/
                Real torStep,
                Real F_TorConRange[MAX_TORS][MAX_TOR_CON][2],
                int N_con[MAX_TORS])

{
    register int i;
    double t;
    int I_ranCon;
    double x0, r1, r2, t1, t2;
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
    /* 
    **  This code probably does *not* produce a uniformly distributed quaternion
    **
    ** change->Q.nx  = Randpm1; 
    ** change->Q.ny  = Randpm1; 
    ** change->Q.nz  = Randpm1; 
    ** change->Q.ang = random_pm( qtwStep );
    ** mkUnitQuat( &(change->Q) );
    */
    /*
    **  This should produce a uniformly distributed quaternion, according to
    **  Shoemake, Graphics Gems III.6, pp.124-132, "Uniform Random Rotations",
    **  published by Academic Press, Inc., (1992)
    */
             t1 = TWOPI * local_random();
    change->Q.x = sin( t1 ) * (  r1 = random_sign * sqrt( 1 - (x0 = local_random()) )  );
    change->Q.y = cos( t1 ) * r1;
             t2 = TWOPI * local_random();
    change->Q.z = sin( t2 ) * (  r2 = random_sign * sqrt( x0 )  );
    change->Q.w = cos( t2 ) * r2;

    /*
    **  Apply random change, to Last Quaternion
    */
    qmultiply( &(now->Q), &(last->Q), &(change->Q) );

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
