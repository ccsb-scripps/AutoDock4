#ifndef INPUT_STATE
#define INPUT_STATE
#include "constants.h"
#include "qmultiply.h"

int input_state( State *S,
		 FILE  *fp, 
                 char  line[LINE_LEN], 
                 int   ntor, 
		 int   *P_istep, 
                 Real *P_energy, 
		 Real *P_eint, 
                 char  *P_lastmove );
#endif
