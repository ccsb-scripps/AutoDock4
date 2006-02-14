#ifndef WEEDBONDS
#define WEEDBONDS
#include "constants.h"
#include "stop.h"

void  weedbonds( int   natom,
                 char  pdbaname[MAX_ATOMS][5],
                 int   piece[MAX_ATOMS],
                 int   ntor,
                 int   tlist[MAX_TORS][MAX_ATOMS],
                 int   nbmatrix_binary[MAX_ATOMS][MAX_ATOMS],
                 int   *P_Nnb,
                 int   **nonbondlist,
                 int   outlev,
				 int   type[MAX_ATOMS]);
#endif


#ifndef PRINT_NONBONDS
#define PRINT_NONBONDS

#include "constants.h"
#include "stop.h"
void print_nonbonds(
                int natom,
                char pdbaname[MAX_ATOMS][5],
                int piece[MAX_ATOMS],
                int ntor,
                int tlist[MAX_TORS][MAX_ATOMS],
                int nbmatrix[MAX_ATOMS][MAX_ATOMS],
                int Nnb,
                int **nonbondlist,
                int outlev,
                int type[MAX_ATOMS]);
#endif

