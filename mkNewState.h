#ifndef MKNEWSTATE
#define MKNEWSTATE

#include "constants.h"
#include "qmultiply.h"
#include "cnv_state_to_coords.h"

void  mkNewState( State *now,
		  State *last,
		  State *change,
		  /*
		  ** FloatOrDouble qtn[NQTN],
                  ** FloatOrDouble tor[MAX_TORS],
                  ** FloatOrDouble qtnLast[NQTN],
                  ** FloatOrDouble torLast[MAX_TORS],
                  ** FloatOrDouble qtnChange[NQTN],
                  ** FloatOrDouble torChange[MAX_TORS],
		  */
                  FloatOrDouble vt[MAX_TORS][NTRN],
                  int   tlist[MAX_TORS][MAX_ATOMS],
                  int   ntor,
                  FloatOrDouble crd[MAX_ATOMS][NTRN],
                  FloatOrDouble crdpdb[MAX_ATOMS][NTRN],
                  int   natom,
                  FloatOrDouble trnStep,
                  /*FloatOrDouble qtwStep,*/
                  FloatOrDouble torStep,
	          FloatOrDouble F_TorConRange[MAX_TORS][MAX_TOR_CON][2],
		  int   N_con[MAX_TORS]);
#endif
