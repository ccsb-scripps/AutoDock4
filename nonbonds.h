
#ifndef NONBONDS
#define NONBONDS
#include "constants.h"
void  nonbonds( float crdpdb[MAX_ATOMS][SPACE],
                int   nbmatrix_binary[MAX_ATOMS][MAX_ATOMS],
                int   natom,
                int   atomnumber[MAX_RECORDS], 
                int   nrecord,
                char  record[MAX_RECORDS][LINE_LEN],
                int   piece[MAX_ATOMS],
                int   Htype,
                int   type[MAX_ATOMS] );
#endif
