#ifndef PRINTENERGIES
#define PRINTENERGIES

#include "autocomm.h"
#include "constants.h"
#include "structs.h"

void printEnergies( EnergyBreakdown *eb,
                    char  *prefixString, 
                    int ligand_is_inhibitor,
                    Real emap_total,
                    Real elec_total,
                    Boole B_have_flexible_residues);

void printStateEnergies( EnergyBreakdown *eb,
			 char  *prefixString, 
			 int ligand_is_inhibitor);
#endif
