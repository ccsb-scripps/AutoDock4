#ifndef WRITEPDBQ
#define WRITEPDBQ
#include "constants.h"
#include "printEnergies.h"
void  writePDBQ( int   irun,
                 char  smFileName[MAX_CHARS],
                 char  dpfFN[MAX_CHARS],
                 FloatOrDouble sml_center[SPACE],
                 State stat,
                 int   ntor,
                 FloatOrDouble eintra,
                 FloatOrDouble einter,
                 int   natom,
                 char  atomstuff[MAX_ATOMS][MAX_CHARS],
                 FloatOrDouble crd[MAX_ATOMS][SPACE],
                 FloatOrDouble emap[MAX_ATOMS],
                 FloatOrDouble elec[MAX_ATOMS],
                 FloatOrDouble charge[MAX_ATOMS],
                 FloatOrDouble abs_charge[MAX_ATOMS],
                 FloatOrDouble qsp_abs_charge[MAX_ATOMS],
                 int ligand_is_inhibitor,
                 FloatOrDouble torsFreeEnergy,
                 int outlev,
                 int   ignore_inter[MAX_ATOMS],
                 const Boole         B_include_1_4_interactions,
                 const FloatOrDouble scale_1_4,
                 const FloatOrDouble sol_fn[NEINT],
                 const ParameterEntry parameterArray[MAX_MAPS],
                 const FloatOrDouble unbound_internal_FE
                 );
#endif
