#ifndef MKNEWSTATE
#define MKNEWSTATE

#include "constants.h"
#include "qmultiply.h"
#include "cnv_state_to_coords.h"

void  mkNewState( State *now,
		  State *last,
		  State *change,
		  /*
		  ** Real qtn[NQTN],
                  ** Real tor[MAX_TORS],
                  ** Real qtnLast[NQTN],
                  ** Real torLast[MAX_TORS],
                  ** Real qtnChange[NQTN],
                  ** Real torChange[MAX_TORS],
		  */
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
		  int   N_con[MAX_TORS]);
#endif
