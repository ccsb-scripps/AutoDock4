
#ifndef CLUSTER_ANALYSIS
#define CLUSTER_ANALYSIS
#include "constants.h"
#include "getrms.h"
int  cluster_analysis( float clus_rms_tol, 
                       int   cluster[MAX_RUNS][MAX_RUNS], 
                       int   num_in_clus[MAX_RUNS], 
                       int   isort[MAX_RUNS], 
                       int   nconf, 
                       int   natom, 
                       int   type[MAX_ATOMS], 
                       float crd[MAX_RUNS][MAX_ATOMS][SPACE], 
                       float crdpdb[MAX_ATOMS][SPACE], 
                       float sml_center[SPACE], 
                       float clu_rms[MAX_RUNS][MAX_RUNS], 
                       Boole B_symmetry_flag,
                       float ref_crds[MAX_ATOMS][SPACE],
                       int   ref_natoms,
                       float ref_rms[MAX_RUNS]);
#endif
