#ifndef INPUT_STATE
#define INPUT_STATE
#include "constants.h"
#include "qmultiply.h"

int input_state( State *S,
		 FILE  *fp, 
                 char  line[LINE_LEN], 
                 int   ntor, 
		 int   *P_istep, 
                 float *P_energy, 
		 float *P_eint, 
                 char  *P_lastmove );
#endif
