
#ifndef CLUSTER_ANALYSIS
#define CLUSTER_ANALYSIS
#include "constants.h"
#include "getrms.h"
int  cluster_analysis( FloatOrDouble clus_rms_tol, 
                       int   cluster[MAX_RUNS][MAX_RUNS], 
                       int   num_in_clus[MAX_RUNS], 
                       int   isort[MAX_RUNS], 
                       int   nconf, 
                       int   natom, 
                       int   type[MAX_ATOMS], 
                       FloatOrDouble crd[MAX_RUNS][MAX_ATOMS][SPACE], 
                       FloatOrDouble crdpdb[MAX_ATOMS][SPACE], 
                       FloatOrDouble sml_center[SPACE], 
                       FloatOrDouble clu_rms[MAX_RUNS][MAX_RUNS], 
                       Boole B_symmetry_flag,
                       FloatOrDouble ref_crds[MAX_ATOMS][SPACE],
                       int   ref_natoms,
                       FloatOrDouble ref_rms[MAX_RUNS]);
#endif
