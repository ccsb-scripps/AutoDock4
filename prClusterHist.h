
#ifndef PRCLUSTERHIST
#define PRCLUSTERHIST
#include "constants.h"
void  prClusterHist(int   ncluster,
                    int   irunmax,
                    float clus_rms_tol,
                    int   num_in_clu[MAX_RUNS],
                    int   cluster[MAX_RUNS][MAX_RUNS],
                    float econf[MAX_RUNS],
                    float clu_rms[MAX_RUNS][MAX_RUNS],
		    float ref_rms[MAX_RUNS]);
#endif
