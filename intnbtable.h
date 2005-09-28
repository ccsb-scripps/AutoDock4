#ifndef INTNBTABLE
#define INTNBTABLE
#include "constants.h"
#include "timesys.h"
#include "structs.h"

void intnbtable(Boole *P_B_havenbp,
                int   *P_a1,
                int   *P_a2,
                GridMapSetInfo *info,
                FloatOrDouble cA,
                FloatOrDouble cB,
                int   xA,
                int   xB,
                double coeff_desolv,
                double sigma,
                EnergyTables *ad_tables);
#endif
