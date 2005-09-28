
#ifndef EINTCALPRINT
#define EINTCALPRINT
#include "constants.h"
FloatOrDouble  eintcalPrint( int   nonbondlist[MAX_NONBONDS][2], 
                             EnergyTables *ptr_ad_energy_tables,
                             FloatOrDouble tcoord[MAX_ATOMS][SPACE], 
                             int   atmtyp[MAX_ATOMS], 
                             int   Nnb,
                             Boole B_calcIntElec,
                             FloatOrDouble q1q2[MAX_NONBONDS],
                             FloatOrDouble abs_charge[MAX_ATOMS]);
#endif
