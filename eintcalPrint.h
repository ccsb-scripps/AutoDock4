
#ifndef EINTCALPRINT
#define EINTCALPRINT
#include "constants.h"
Real  eintcalPrint( int   **nonbondlist, 
                             EnergyTables *ptr_ad_energy_tables,
                             Real tcoord[MAX_ATOMS][SPACE], 
                             int   atmtyp[MAX_ATOMS], 
                             int   Nnb,
                             Boole B_calcIntElec,
                             Real abs_charge[MAX_ATOMS]);
#endif
