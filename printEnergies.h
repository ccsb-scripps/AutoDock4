#ifndef PRINTENERGIES
#define PRINTENERGIES
#include <stdio.h>
#include "autocomm.h"
void printEnergies( Real einter,
		    Real eintra,
		    Real torsFreeEnergy,
		    char  *prefixString, 
            int ligand_is_inhibitor,
			 Real elec_total,
			 Real emap_total,
			 Real unbound_internal_FE
		    );

void printStateEnergies( Real einter,
			 Real eintra,
			 Real torsFreeEnergy,
			 char  *prefixString, 
			 int ligand_is_inhibitor,
			 Real unbound_internal_FE
			 );
#endif
