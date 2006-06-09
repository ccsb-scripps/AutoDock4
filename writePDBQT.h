#ifndef _WRITEPDBQT
#define _WRITEPDBQT

#include "structs.h"
#include "constants.h"
#include "printEnergies.h"
#include "trilinterp.h"
#include "eintcal.h"
#include "cnv_state_to_coords.h"
#include "stateLibrary.h"

void writePDBQT(int irun,FourByteLong seed[2],
                    char  smFileName[MAX_CHARS],
                    char  dpfFN[MAX_CHARS],
                    Real sml_center[SPACE],
                    State state,
                    int   ntor,
                    Real (*Ptr_eintra),
                    Real (*Ptr_einter),
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
                    Real vt[MAX_TORS][SPACE],
                    int   tlist[MAX_TORS][MAX_ATOMS],
                    Real crdpdb[MAX_ATOMS][SPACE],
                    int   **nonbondlist,
                    EnergyTables *ptr_ad_energy_tables,
                    int   type[MAX_ATOMS],
                    int   Nnb,
                    Boole B_calcIntElec,
                    Real q1q2[MAX_NONBONDS],
                    Real map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS],
                    int outlev,
                    int   ignore_inter[MAX_ATOMS],
                    const Boole         B_include_1_4_interactions,
                    const Real scale_1_4,
                    const ParameterEntry parameterArray[MAX_MAPS],
                    const Real unbound_internal_FE,
                    GridMapSetInfo *info,
                    int state_type,  // 0 means unbound, 1 means docked
                    char PDBQT_record[MAX_RECORDS][LINE_LEN],
                    Boole B_use_non_bond_cutoff,
                    Boole B_have_flexible_residues
                    );

#endif
