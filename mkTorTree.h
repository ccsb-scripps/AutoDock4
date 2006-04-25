#ifndef MKTORTREE
#define MKTORTREE

#include "constants.h"
#include "parse_PDBQT_line.h"
#include "stop.h"

void  mkTorTree(int   atomnumber[MAX_RECORDS],
                char  record[MAX_RECORDS][LINE_LEN],
                int   nrecord,
                int   tlist[MAX_TORS][MAX_ATOMS],
                int   *P_ntor,
                char  smFileName[MAX_CHARS],
                char  pdbaname[MAX_ATOMS][5],
                Boole *P_B_constrain,
                int   *P_atomC1,
                int   *P_atomC2,
                Real *P_sqlower,
                Real *P_squpper,
                int   *P_ntorsdof,
                int   ignore_inter[MAX_ATOMS]);
#endif
