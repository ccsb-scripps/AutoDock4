#ifndef PRINTENERGIES
#define PRINTENERGIES
#include <stdio.h>
#include "autocomm.h"
void printEnergies( float einter,
		    float eintra,
		    float torsFreeEnergy,
		    char  *prefixString, 
        int ligand_is_inhibitor
		    );
#endif
