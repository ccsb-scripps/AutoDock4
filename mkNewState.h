#ifndef MKNEWSTATE
#define MKNEWSTATE

#include "constants.h"
#include "qmultiply.h"
#include "cnv_state_to_coords.h"

void  mkNewState( State *now,
		  State *last,
		  State *change,
		  /*
		  ** float qtn[NQTN],
                  ** float tor[MAX_TORS],
                  ** float qtnLast[NQTN],
                  ** float torLast[MAX_TORS],
                  ** float qtnChange[NQTN],
                  ** float torChange[MAX_TORS],
		  */
                  float vt[MAX_TORS][NTRN],
                  int   tlist[MAX_TORS][MAX_ATOMS],
                  int   ntor,
                  float crd[MAX_ATOMS][NTRN],
                  float crdpdb[MAX_ATOMS][NTRN],
                  int   natom,
                  float trnStep,
                  /*float qtwStep,*/
                  float torStep,
	          float F_TorConRange[MAX_TORS][MAX_TOR_CON][2],
		  int   N_con[MAX_TORS]);
#endif
