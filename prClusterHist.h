
#ifndef PRCLUSTERHIST
#define PRCLUSTERHIST
#include "constants.h"
void  prClusterHist(int   ncluster,
                    int   irunmax,
                    Real clus_rms_tol,
                    int   num_in_clu[MAX_RUNS],
                    int   cluster[MAX_RUNS][MAX_RUNS],
                    Real econf[MAX_RUNS],
                    Real clu_rms[MAX_RUNS][MAX_RUNS],
		    Real ref_rms[MAX_RUNS]);
#endif
