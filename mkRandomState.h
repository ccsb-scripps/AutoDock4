/* mkRandomState.h */


#ifndef MKRANDOMSTATE
#define MKRANDOMSTATE

#include "constants.h"

State mkRandomState( FloatOrDouble xlo,
		     FloatOrDouble xhi,
		     FloatOrDouble ylo,
		     FloatOrDouble yhi,
		     FloatOrDouble zlo,
		     FloatOrDouble zhi,
		     int   ntor,
		     FloatOrDouble F_TorConRange[MAX_TORS][MAX_TOR_CON][2],
		     int   N_con[MAX_TORS]);
#endif
