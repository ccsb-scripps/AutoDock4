#ifndef WEEDBONDS
#define WEEDBONDS
#include "constants.h"
#include "stop.h"

void  weedbonds( int   natom,
                 char  pdbaname[MAX_ATOMS][5],
                 int   piece[MAX_ATOMS],
                 int   ntor,
                 int   tlist[MAX_TORS][MAX_ATOMS],
                 int   *P_Nnb,
                 int   Nnbonds[MAX_ATOMS],
                 int   nbmatrix_binary[MAX_ATOMS][MAX_ATOMS],
                 int   nonbondlist[MAX_NONBONDS][4],
                 int   outlev,
				 int   type[MAX_ATOMS]);
#endif
