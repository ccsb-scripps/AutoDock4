#ifndef NBE
#define NBE

#include "constants.h"

void  nbe( char  atm_tyP_str[ATOM_MAPS],
	   FloatOrDouble e_internal[NEINT][ATOM_MAPS][ATOM_MAPS],
           int   num_atm_maps );
#endif
