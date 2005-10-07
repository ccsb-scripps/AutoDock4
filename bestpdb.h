#ifndef BESTPDB
#define BESTPDB

#include "constants.h"
#include "print_rem.h"
#include "strindex.h"
#include "print_avsfld.h"

void  bestpdb( int   ncluster, 
               int   num_in_clu[MAX_RUNS], 
               int   cluster[MAX_RUNS][MAX_RUNS], 
               FloatOrDouble econf[MAX_RUNS], 
               FloatOrDouble crd[MAX_RUNS][MAX_ATOMS][SPACE], 
               char  atomstuff[MAX_ATOMS][MAX_CHARS], 
               int   natom, 
               Boole B_write_all_clusmem, 
               FloatOrDouble ref_rms[MAX_RUNS]);
#endif
