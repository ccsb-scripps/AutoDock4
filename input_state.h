#ifndef INPUT_STATE
#define INPUT_STATE
#include "constants.h"
#include "qmultiply.h"

int input_state( State *S,
		 FILE  *fp, 
                 char  line[LINE_LEN], 
                 int   ntor, 
		 int   *P_istep, 
                 FloatOrDouble *P_energy, 
		 FloatOrDouble *P_eint, 
                 char  *P_lastmove );
#endif
