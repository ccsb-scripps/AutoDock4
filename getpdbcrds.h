#ifndef GETPDBCRDS
#define GETPDBCRDS

#include "constants.h"
#include "openfile.h"

int   getpdbcrds( char  rms_ref_crds_FN[MAX_CHARS],
		Real ref_crds[MAX_ATOMS][SPACE] );
#endif
