
#ifndef CLUSTER_ANALYSIS
#define CLUSTER_ANALYSIS
#include "constants.h"
#include "getrms.h"
int  cluster_analysis( Real clus_rms_tol, 
                       int   cluster[MAX_RUNS][MAX_RUNS], 
                       int   num_in_clus[MAX_RUNS], 
                       int   isort[MAX_RUNS], 
                       int   nconf, 
                       int   natom, 
                       int   type[MAX_ATOMS], 
                       Real crd[MAX_RUNS][MAX_ATOMS][SPACE], 
                       Real crdpdb[MAX_ATOMS][SPACE], 
                       Real sml_center[SPACE], 
                       Real clu_rms[MAX_RUNS][MAX_RUNS], 
                       Boole B_symmetry_flag,
                       Real ref_crds[MAX_ATOMS][SPACE],
                       int   ref_natoms,
                       Real ref_rms[MAX_RUNS]);
#endif
