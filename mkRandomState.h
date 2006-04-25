/* mkRandomState.h */


#ifndef MKRANDOMSTATE
#define MKRANDOMSTATE

#include "constants.h"

State mkRandomState( int   ntor,
		     Real F_TorConRange[MAX_TORS][MAX_TOR_CON][2],
		     int   N_con[MAX_TORS],
                     GridMapSetInfo *info);
#endif
