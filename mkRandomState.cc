/*

 $Id: mkRandomState.cc,v 1.4 2006/04/25 22:32:34 garrett Exp $

*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/* mkRandomState.cc */

#include <math.h>

    #include <stdlib.h>
    #include "structs.h"
    #include "mkRandomState.h"


State mkRandomState( int   ntor,
		     Real F_TorConRange[MAX_TORS][MAX_TOR_CON][2],
		     int   N_con[MAX_TORS],
                     GridMapSetInfo *info)

{
    State now;
    register int i;
    double t;
    int I_ranCon;
    double x0, r1, r2, t1, t2;
    Real a, b;

    now.ntor = ntor;

    /*
    ** Translation
    */
    now.T.x = random_range( info->lo[X], info->hi[X]);
    now.T.y = random_range( info->lo[Y], info->hi[Y]);
    now.T.z = random_range( info->lo[Z], info->hi[Z]);

    /*
    ** Quaternion angular displacement
    */
    /*
    **  This should produce a uniformly distributed quaternion, according to
    **  Shoemake, Graphics Gems III.6, pp.124-132, "Uniform Random Rotations",
    **  published by Academic Press, Inc., (1992)
    */
         x0 = local_random();
         r1 = random_sign * sqrt( 1 - x0 );
         t1 = TWOPI * local_random();
    now.Q.x = sin( t1 ) * r1;
    now.Q.y = cos( t1 ) * r1;
         r2 = random_sign * sqrt(     x0 );
         t2 = TWOPI * local_random();
    now.Q.z = sin( t2 ) * r2;
    now.Q.w = cos( t2 ) * r2;

    for (i=0; i<ntor; i++) {
	if (N_con[i] > 0) {
	    if (N_con[i] > 1) {
		/* If N_con was 2, I_ranCon could be 0 or 1, never 2 */
		/* Select a random constraint */
		I_ranCon = (int)((double)N_con[i] * local_random());  
	    } else {
		/* Hobson's choice...
		*/
		I_ranCon = 0;
	    }
	    a = F_TorConRange[i][I_ranCon][LOWER];
	    b = F_TorConRange[i][I_ranCon][UPPER];
	    t = random_range(a,b);
	    now.tor[i] = WrpModRad(t);
	} else {
	    now.tor[i] = random_range( -PI, PI );
	}
    }/*i*/

    return( now );
} 
/* EOF */
