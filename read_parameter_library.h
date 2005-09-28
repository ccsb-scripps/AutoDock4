#ifndef _READ_PARAMETER_LIBRARY
#define _READ_PARAMETER_LIBRARY

#include "autocomm.h"

void read_parameter_library(
        char FN_parameter_library[MAX_CHARS],
        int outlev
        );

void setup_parameter_library(
        int outlev
        );

void setup_distdepdiel( int outlev, 
                        EnergyTables *ptr_ad_energy_tables  // Holds vdw+Hb, desolvation & dielectric lookup tables
                      );


#endif
