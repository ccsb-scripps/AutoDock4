
#ifndef EINTCAL
#define EINTCAL

#include "constants.h"
#include "structs.h"

#ifndef EINTCALPRINT

FloatOrDouble  eintcal( 
                        int   ** const nonbondlist, 
                        const EnergyTables  *ad_energy_tables,
                        const FloatOrDouble tcoord[MAX_ATOMS][SPACE], 
                        const int           Nnb,
                        const Boole         B_calcIntElec,
                        const FloatOrDouble q1q2[MAX_NONBONDS],
                        const Boole         B_include_1_4_interactions,
                        const FloatOrDouble scale_1_4,
                        const FloatOrDouble qsp_abs_charge[MAX_ATOMS],
                        const ParameterEntry parameterArray[MAX_MAPS],
                        const FloatOrDouble unbound_internal_FE
                      );

#else        /*EINTCALPRINT*/

FloatOrDouble  eintcalPrint( 
                             int  ** const nonbondlist,
                             const EnergyTables  *ad_energy_tables,
                             const FloatOrDouble tcoord[MAX_ATOMS][SPACE],
                             const int   Nnb,
                             const Boole B_calcIntElec,
                             const FloatOrDouble q1q2[MAX_NONBONDS],
                             const Boole B_include_1_4_interactions,
                             const FloatOrDouble scale_1_4,
                             const FloatOrDouble qsp_abs_charge[MAX_ATOMS],
                             const ParameterEntry parameterArray[MAX_MAPS],
                             const FloatOrDouble unbound_internal_FE
                           );

#endif        /*EINTCALPRINT*/

#endif        /*!EINTCAL*/
