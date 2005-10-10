/*
 $Id: cmdmode.h,v 1.7.6.1 2005/10/10 16:44:18 alther Exp $
*/

#ifndef CMDMODE
#define CMDMODE

#include "autocomm.h"
#include "grid.h"
#include "typedefs.h"
#include "structs.h"
//#include "constants.h"
//#include "set_cmd_io_std.h"
//#include "print_2x.h"
//#include "parse_com_line.h"
//#include "strindex.h"
//#include "print_avsfld.h"
//#include "printEnergies.h"
//#include "success.h"
//#include "readPDBQT.h"
//#include "get_atom_type.h"
//#include "timesys.h"
//#include "eintcal.h"
//#include "trilinterp.h"
//#include "qmultiply.h"
//#include "cnv_state_to_coords.h"
//#include "parse_trj_line.h"
//#include "input_state.h"
//#include "openfile.h"
int   cmdmode(int           natom,
              Clock         jobStart,
              struct tms    tms_jobStart,
              FloatOrDouble map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS],

              EnergyTables *ptr_ad_energy_tables,

              FloatOrDouble WallEnergy,
              FloatOrDouble vt[MAX_TORS][SPACE],
              int   tlist[MAX_TORS][MAX_ATOMS],
              int   ntor,
              int   Nnb,
              int   nonbondlist[MAX_NONBONDS][MAX_NBDATA],
              char  atomstuff[MAX_ATOMS][MAX_CHARS],
              FloatOrDouble crdpdb[MAX_ATOMS][SPACE],
              char  hostnm[MAX_CHARS],
              int   type[MAX_ATOMS],
              FloatOrDouble charge[MAX_ATOMS],
              FloatOrDouble abs_charge[MAX_ATOMS],
              FloatOrDouble qsp_abs_charge[MAX_ATOMS],
              Boole B_calcIntElec,
              FloatOrDouble q1q2[MAX_NONBONDS],
              char  atm_typ_str[ATOM_MAPS],
              FloatOrDouble torsFreeEnergy,
              int ligand_is_inhibitor,
              int ignore_inter[MAX_ATOMS],
              const Boole         B_include_1_4_interactions,
              const FloatOrDouble scale_1_4,

              const ParameterEntry parameterArray[MAX_MAPS],
              const FloatOrDouble unbound_internal_FE,

              GridMapSetInfo *info
             );
#endif
