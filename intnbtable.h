#ifndef INTNBTABLE
#define INTNBTABLE
#include "constants.h"
#include "timesys.h"

void intnbtable(Boole *P_B_havenbp,
                int   *P_a1,
                int   *P_a2,
		        int   num_atm_maps,
                char ligand_atom_types[MAX_MAPS][3],
                FloatOrDouble cA,
                FloatOrDouble cB,
                int   xA,
                int   xB,
                FloatOrDouble e_internal[NEINT][ATOM_MAPS][ATOM_MAPS],
                FloatOrDouble sol_fn[NEINT],
                double coeff_desolv,
                double sigma);
#endif
