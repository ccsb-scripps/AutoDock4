
#ifndef EINTCAL
#define EINTCAL

#include "constants.h"
#include "structs.h"

#ifndef EINTCALPRINT

Real  eintcal( NonbondParam * const nonbondlist, 
               const EnergyTables  *ad_energy_tables,
               const Real tcoord[MAX_ATOMS][SPACE], 
               const int  Nnb,
               const Boole B_calcIntElec,
               const Boole B_include_1_4_interactions,
               const Real scale_1_4,
               const Real qsp_abs_charge[MAX_ATOMS],
               const ParameterEntry parameterArray[MAX_MAPS],
               const Boole B_use_non_bond_cutoff,
               const Boole B_have_flexible_residues  );

#else        /*EINTCALPRINT*/

Real  eintcalPrint( NonbondParam * const nonbondlist,
                     const EnergyTables  *ad_energy_tables,
                     const Real tcoord[MAX_ATOMS][SPACE],
                     const int   Nnb,
                     const Boole B_calcIntElec,
                     const Boole B_include_1_4_interactions,
                     const Real scale_1_4,
                     const Real qsp_abs_charge[MAX_ATOMS],
                     const ParameterEntry parameterArray[MAX_MAPS],
                     const Boole B_use_non_bond_cutoff,
                     const Boole B_have_flexible_residues);  // if the receptor has flexibile residues, this will be set to TRUE

#endif        /*EINTCALPRINT*/

#endif        /*!EINTCAL*/
