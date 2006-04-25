
#ifndef EINTCAL
#define EINTCAL

#include "constants.h"
#include "structs.h"

#ifndef EINTCALPRINT

Real  eintcal( 
                        int   ** const nonbondlist, 
                        const EnergyTables  *ad_energy_tables,
                        const Real tcoord[MAX_ATOMS][SPACE], 
                        const int           Nnb,
                        const Boole         B_calcIntElec,
                        const Real q1q2[MAX_NONBONDS],
                        const Boole         B_include_1_4_interactions,
                        const Real scale_1_4,
                        const Real qsp_abs_charge[MAX_ATOMS],
                        const ParameterEntry parameterArray[MAX_MAPS]
                      );

#else        /*EINTCALPRINT*/

Real  eintcalPrint( 
                             int  ** const nonbondlist,
                             const EnergyTables  *ad_energy_tables,
                             const Real tcoord[MAX_ATOMS][SPACE],
                             const int   Nnb,
                             const Boole B_calcIntElec,
                             const Real q1q2[MAX_NONBONDS],
                             const Boole B_include_1_4_interactions,
                             const Real scale_1_4,
                             const Real qsp_abs_charge[MAX_ATOMS],
                             const ParameterEntry parameterArray[MAX_MAPS]
                           );

#endif        /*EINTCALPRINT*/

#endif        /*!EINTCAL*/
