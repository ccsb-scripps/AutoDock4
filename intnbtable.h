#ifndef INTNBTABLE
#define INTNBTABLE
#include "constants.h"
#include "timesys.h"

void intnbtable(Boole *P_B_havenbp,
                int   *P_a1,
                int   *P_a2,
		int   num_atm_maps,
                char  atm_tyP_str[ATOM_MAPS],
                float cA,
                float cB,
                int   xA,
                int   xB,
                float e_internal[NEINT][ATOM_MAPS][ATOM_MAPS] );
#endif
