
#ifndef OUTPUT_STATE
#define OUTPUT_STATE
#include "constants.h"
void  output_state( FILE  *fp,
		    State S,
                    int   ntor,
                    int   istep,
                    float energy,
                    float eint,
                    char  lastmove,
                    Boole B_watch,
                    char  FN_watch[MAX_CHARS],
                    char  atomstuff[MAX_ATOMS][MAX_CHARS],
                    int   natom,
                    float crd[MAX_ATOMS][SPACE]);
#endif
