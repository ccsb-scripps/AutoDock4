
#ifndef OUTPUT_STATE
#define OUTPUT_STATE
#include "constants.h"
void  output_state( FILE  *fp,
		    State S,
                    int   ntor,
                    int   istep,
                    Real energy,
                    Real eint,
                    char  lastmove,
                    Boole B_watch,
                    char  FN_watch[MAX_CHARS],
                    char  atomstuff[MAX_ATOMS][MAX_CHARS],
                    int   natom,
                    Real crd[MAX_ATOMS][SPACE]);
#endif
