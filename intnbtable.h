#ifndef INTNBTABLE
#define INTNBTABLE
#include "constants.h"
#include "timesys.h"
#include "structs.h"

void intnbtable(Boole *P_B_havenbp,
                int   a1,
                int   a2,
                GridMapSetInfo *info,
                Real cA,
                Real cB,
                int   xA,
                int   xB,
                double coeff_desolv,
                double sigma,
                EnergyTables *ad_tables,
                Boole B_is_unbound_calculation );
#endif
