#ifndef WRITEPDBQ
#define WRITEPDBQ
#include "constants.h"
#include "printEnergies.h"
void  writePDBQ( int   irun,
                 char  smFileName[MAX_CHARS],
                 char  dpfFN[MAX_CHARS],
                 Real sml_center[SPACE],
                 State stat,
                 int   ntor,
                 Real eintra,
                 Real einter,
                 int   natom,
                 char  atomstuff[MAX_ATOMS][MAX_CHARS],
                 Real crd[MAX_ATOMS][SPACE],
                 Real emap[MAX_ATOMS],
                 Real elec[MAX_ATOMS],
                 Real charge[MAX_ATOMS],
                 Real abs_charge[MAX_ATOMS],
                 Real qsp_abs_charge[MAX_ATOMS],
                 int ligand_is_inhibitor,
                 Real torsFreeEnergy,
                 int outlev,
                 int   ignore_inter[MAX_ATOMS],
                 const Boole         B_include_1_4_interactions,
                 const Real scale_1_4,

                 const ParameterEntry parameterArray[MAX_MAPS],
                 const Real unbound_internal_FE
                 );
#endif
