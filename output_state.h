
#ifndef OUTPUT_STATE
#define OUTPUT_STATE
#include "constants.h"
void  output_state( FILE  *fp,
		    State S,
                    int   ntor,
                    int   istep,
                    FloatOrDouble energy,
                    FloatOrDouble eint,
                    char  lastmove,
                    Boole B_watch,
                    char  FN_watch[MAX_CHARS],
                    char  atomstuff[MAX_ATOMS][MAX_CHARS],
                    int   natom,
                    FloatOrDouble crd[MAX_ATOMS][SPACE]);
#endif
