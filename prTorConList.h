#ifndef PRTORCONLIST
#define PRTORCONLIST
#include "constants.h"

void  prTorConList( int   ntor,
		    Boole B_isTorConstrained[MAX_TORS],
		    unsigned short US_torProfile[MAX_TORS][NTORDIVS],
		    float F_TorConRange[MAX_TORS][MAX_TOR_CON][2],
		    int   N_con[MAX_TORS]);
#endif
