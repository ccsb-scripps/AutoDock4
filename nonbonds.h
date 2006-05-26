
#ifndef NONBONDS
#define NONBONDS
#include "constants.h"
void  nonbonds( const Real crdpdb[MAX_ATOMS][SPACE],  
			    int         nbmatrix[MAX_ATOMS][MAX_ATOMS],
			    const int   natom, 
			    const int   bond_index[MAX_ATOMS],
                int         B_include_1_4_interactions,
                int         bonded[MAX_ATOMS][6]);
#endif

#ifndef GETBONDS
#define GETBONDS
#include "constants.h"
void getbonds(const Real crdpdb[MAX_ATOMS][SPACE], 
			  const int natom, 
			  const int bond_index[MAX_ATOMS],
              int bonded[MAX_ATOMS][6]);
#endif

#ifndef PRINTBONDS
#define PRINTBONDS
#include "constants.h"
void printbonds(const int natom, const int bonded[MAX_ATOMS][6], const char *message, const int B_print_all_bonds);
#endif

#ifndef PRINT14
#define PRINT14
#include "constants.h"
#include <stdio.h>
void print_1_4_message(FILE *file, Boole B_include_1_4_interactions,  Real scale_1_4);
#endif
