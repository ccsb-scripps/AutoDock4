#ifndef _CALL_CPSO
#define _CALL_CPSO

#include "stateLibrary.h"
#include "cnv_state_to_coords.h"
#include "structs.h"
#include "constants.h"
#include "alea.h"
#include "dimLibrary.h"



//State call_cpso(int n_exec, State sInit, int S, int D, double *xmin, double *xmax, int eval_max, int K, double c1, double c2, int outlev);
State call_cpso(int n_exec, State sInit, int S, int D, double *xmin, double *xmax, int eval_max, int K, double c1, double c2, int outlev, 
	int ntor, unsigned int max_its,
        unsigned int max_succ, unsigned int max_fail,
	float expansion, float contraction,
	float search_freq, float *rho,
	float *lb_rho);

State SWLocalSearch(State *state, int ntor, unsigned int max_its, unsigned int max_succ, unsigned int max_fail,
	float expansion, float contraction, float *rho, float *lb_rho_ptr);

void addDeviatesBias(State* srcstate, State* deststate, float* deviates, float* bias);

#endif
