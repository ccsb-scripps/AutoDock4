
#ifndef EINTCAL
#define EINTCAL

#include "constants.h"
#include "structs.h"

#ifndef EINTCALPRINT

FloatOrDouble  eintcal( 
                        const int           nonbondlist[MAX_NONBONDS][MAX_NBDATA], 
                        const FloatOrDouble e_internal[NEINT][ATOM_MAPS][ATOM_MAPS], 
                        const FloatOrDouble tcoord[MAX_ATOMS][SPACE], 
                        const int           Nnb,
                        const Boole         B_calcIntElec,
                        const FloatOrDouble q1q2[MAX_NONBONDS],
                        const Boole         B_include_1_4_interactions,
                        const FloatOrDouble scale_1_4,
                        const FloatOrDouble qsp_abs_charge[MAX_ATOMS],
                        const FloatOrDouble sol_fn[NEINT],
                        const ParameterEntry parameterArray[MAX_MAPS],
                        const FloatOrDouble unbound_internal_FE
                      );

#else        /*EINTCALPRINT*/

FloatOrDouble  eintcalPrint( 
                             const int   nonbondlist[MAX_NONBONDS][MAX_NBDATA],
                             const FloatOrDouble e_internal[NEINT][ATOM_MAPS][ATOM_MAPS],
                             const FloatOrDouble tcoord[MAX_ATOMS][SPACE],
                             const int   Nnb,
                             const Boole B_calcIntElec,
                             const FloatOrDouble q1q2[MAX_NONBONDS],
                             const Boole B_include_1_4_interactions,
                             const FloatOrDouble scale_1_4,
                             const FloatOrDouble qsp_abs_charge[MAX_ATOMS],
                             const FloatOrDouble sol_fn[NEINT],
                             const ParameterEntry parameterArray[MAX_MAPS],
                             const FloatOrDouble unbound_internal_FE
                           );

#endif        /*EINTCALPRINT*/

#endif        /*!EINTCAL*/
