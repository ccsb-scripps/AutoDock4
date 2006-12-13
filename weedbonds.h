#ifndef WEEDBONDS
#define WEEDBONDS

#include "constants.h"
#include "stop.h"
#include "structs.h"

void  weedbonds( int   natom,
                 char  pdbaname[MAX_ATOMS][5],
                 int   piece[MAX_ATOMS],
                 int   ntor,
                 int   tlist[MAX_TORS][MAX_ATOMS],
                 int   nbmatrix_binary[MAX_ATOMS][MAX_ATOMS],
                 int   *P_Nnb,
                 NonbondParam *nonbondlist,
                 int   outlev,
				 int   type[MAX_ATOMS]);
#endif


#ifndef PRINT_NONBONDS
#define PRINT_NONBONDS

#include "constants.h"
#include "stop.h"
#include "structs.h"

void print_nonbonds(
                int natom,
                char pdbaname[MAX_ATOMS][5],
                int piece[MAX_ATOMS],
                int ntor,
                int tlist[MAX_TORS][MAX_ATOMS],
                int nbmatrix[MAX_ATOMS][MAX_ATOMS],
                int Nnb,
                NonbondParam *nonbondlist,
                int outlev,
                int type[MAX_ATOMS]);
#endif

