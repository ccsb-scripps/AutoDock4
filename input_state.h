/*
 $Id: input_state.h,v 1.2.8.1 2005/10/11 00:05:46 alther Exp $
*/

#ifndef INPUT_STATE
#define INPUT_STATE

#include <stdio.h>
#include "autocomm.h"
#include "structs.h"
#include "typedefs.h"

int input_state(State*         S,
                FILE*          fp,
                char           line[LINE_LEN],
                int            ntor,
		          int*           P_istep,
                FloatOrDouble* P_energy,
		          FloatOrDouble* P_eint,
                char*          P_lastmove );
#endif
