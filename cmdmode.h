#ifndef CMDMODE
#define CMDMODE

#include "constants.h"
#include "set_cmd_io_std.h"
#include "print_2x.h"
#include "parse_com_line.h"
#include "strindex.h"
#include "print_avsfld.h"
#include "printEnergies.h"
#include "success.h"
#include "readPDBQT.h"
#include "get_atom_type.h"
#include "timesys.h"
#include "eintcal.h"
#include "trilinterp.h"
#include "qmultiply.h"
#include "cnv_state_to_coords.h"
#include "parse_trj_line.h"
#include "input_state.h"
#include "openfile.h"

int   cmdmode( int natom,
             Clock jobStart,
             struct tms tms_jobStart,
             Real map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS],

             EnergyTables *ptr_ad_energy_tables,

             Real WallEnergy,
             Real vt[MAX_TORS][SPACE],
             int   tlist[MAX_TORS][MAX_ATOMS],
             int   ntor,
             int   Nnb,
             int   **nonbondlist,
             char  atomstuff[MAX_ATOMS][MAX_CHARS],
             Real crdpdb[MAX_ATOMS][SPACE],
             char  hostnm[MAX_CHARS],
             int   type[MAX_ATOMS],
             Real charge[MAX_ATOMS],
             Real abs_charge[MAX_ATOMS],
             Real qsp_abs_charge[MAX_ATOMS],
             Boole B_calcIntElec,
             Real q1q2[MAX_NONBONDS],
             char  atm_typ_str[ATOM_MAPS],
             Real torsFreeEnergy,
             int ligand_is_inhibitor,
             int ignore_inter[MAX_ATOMS],
             const Boole         B_include_1_4_interactions,
             const Real scale_1_4,

             const ParameterEntry parameterArray[MAX_MAPS],
             const Real unbound_internal_FE,

             GridMapSetInfo *info,
             Boole B_have_flexible_residues
             );
#endif
