#ifndef PRINTENERGIES
#define PRINTENERGIES
#include <stdio.h>
#include "autocomm.h"
void printEnergies( FloatOrDouble einter,
		    FloatOrDouble eintra,
		    FloatOrDouble torsFreeEnergy,
		    char  *prefixString, 
            int ligand_is_inhibitor,
			 FloatOrDouble elec_total,
			 FloatOrDouble emap_total,
			 FloatOrDouble unbound_internal_FE
		    );

void printStateEnergies( FloatOrDouble einter,
			 FloatOrDouble eintra,
			 FloatOrDouble torsFreeEnergy,
			 char  *prefixString, 
			 int ligand_is_inhibitor,
			 FloatOrDouble unbound_internal_FE
			 );
#endif
