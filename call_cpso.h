#ifndef _CALL_CPSO
#define _CALL_CPSO

#include "stateLibrary.h"
#include "cnv_state_to_coords.h"
#include "structs.h"
#include "constants.h"

State call_cpso(
	int n_exec, 
	State sInit, 
	int S, int D, 
	double *xmin, 
	double *xmax, 
	long eval_max, 
	int K, 
	double c1, 
	double c2, 
	int outlev, 
	//int ntor, 
	unsigned int max_its,
    unsigned int max_succ, 
    unsigned int max_fail,
	float expansion, 
	float contraction,
	float search_freq, 
	float *rho_ptr,
	float *lb_rho_ptr,
	float rho,
	float lb_rho);

void make_state_from_position(State *S, Position P);
void copy_state_to_position(Position *R , State S);
void initialiseDimension(GridMapSetInfo *info, double *xmin, double *xmax, int D);

void initialiseParticle(
			int i, 
			int D, 
			Position *Xi, 
			Velocity *Vi, 
			double *xmin,
			double *xmax, 
			double *Vmin, 
			double *Vmax);
	
void initialize_random_particle ( 
			int i,  			// particle i
			int D,    			// degree of freedom
 			Position *Xi, 		// pointer to Positions
 			Velocity *Vi, 		// pointer to Velocity of particles
 			double *xmin,		// Min. x,y,z of gridInfoMap
 			double *xmax,       // Max. x,y,z of gridInfoMap
 			double *Vmin, 
 			double *Vmax );
		 
void swarmActivity(int S, int D, Position *Xi, int nb_eval, int outlev);


State PSWLocalSearch(
	State &state, 
	int ntor,  
	int max_its,  
	int max_succ,  
	int max_fail,
	int ilig, 
	int from_idx, 
	int to_idx,	 
	float *rho, 
	float *lb_rho_ptr);

State SWLocalSearch(
	State &state,
	int ilig,
	int from_idx,
	int to_idx, 
	int ntor, 
	int max_its, 
	int max_succ, 
	int max_fail,		
	float expansion, 
	float contraction, 
	float  rho,
	float  lb_rho  
	);
			
State addDeviatesBiasAD4(
	State &srcstate,
	float* deviates, 
	float* bias, 
	int ilig, 
	double sign);

#endif
