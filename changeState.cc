/* changeState.cc */

#include <math.h>

    #include <stdlib.h>
    #include "constants.h"
    #include "qmultiply.h"
    #include "changeState.h"


State  changeState( State last,      /* ...must be a normalized quaternion! */
		    float trnStep,
		    float torStep,
		    int   ntor,
		    float F_TorConRange[MAX_TORS][MAX_TOR_CON][2],
		    int   N_con[MAX_TORS])

{
    register int i;
    double t;
    int I_ranCon;
    double x0, r1, r2, t1, t2;
    Quat changeQuat;
    State now;

    /*
    ** Translation
    */
    now.T.x = last.T.x + random_pm( trnStep );
    now.T.y = last.T.y + random_pm( trnStep );
    now.T.z = last.T.z + random_pm( trnStep );

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
    changeQuat.x = sin( t1 ) * r1;
    changeQuat.y = cos( t1 ) * r1;
	      r2 = random_sign * sqrt(     x0 );
              t2 = TWOPI * local_random();
    changeQuat.z = sin( t2 ) * r2;
    changeQuat.w = cos( t2 ) * r2;

    /*
    **  Apply random change, to Last Quaternion
    */
    qmultiply( &(now.Q), &(last.Q), &(changeQuat) );

    now.ntor = ntor;
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
	    t = random_range(F_TorConRange[i][I_ranCon][LOWER], 
			     F_TorConRange[i][I_ranCon][UPPER]);
	    now.tor[i] = WrpModRad(t);
	} else {
	    t = last.tor[i] + random_pm( torStep );
	    now.tor[i] = WrpModRad(t);
	}
    }/*i*/

    return(now);
} 
/* EOF */
