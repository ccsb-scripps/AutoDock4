/* mkRandomState.h */


#ifndef MKRANDOMSTATE
#define MKRANDOMSTATE

#include "constants.h"

State mkRandomState( float xlo,
		     float xhi,
		     float ylo,
		     float yhi,
		     float zlo,
		     float zhi,
		     int   ntor,
		     float F_TorConRange[MAX_TORS][MAX_TOR_CON][2],
		     int   N_con[MAX_TORS]);
#endif
