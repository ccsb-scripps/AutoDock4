
#ifndef PRCLUSTERHIST
#define PRCLUSTERHIST
#include "constants.h"
void  prClusterHist(int   ncluster,
                    int   irunmax,
                    FloatOrDouble clus_rms_tol,
                    int   num_in_clu[MAX_RUNS],
                    int   cluster[MAX_RUNS][MAX_RUNS],
                    FloatOrDouble econf[MAX_RUNS],
                    FloatOrDouble clu_rms[MAX_RUNS][MAX_RUNS],
		    FloatOrDouble ref_rms[MAX_RUNS]);
#endif
