/* call_cpso.cc */

#include <stdio.h>
#include <math.h>


#include <time.h>      // time_t time(time_t *tloc);
#include <sys/times.h>
#include <stdlib.h>
#include <string.h>
#include <sys/param.h>
#include <ctype.h> // tolower

#include "call_cpso.h"
#include "printdate.h"
#include "timesyshms.h"
#include "assert.h"
#include "eval.h"
#include "qmultiply.h"
#include "ranlib.h"

extern FILE *logFile;
extern Eval evaluate;
extern int nlig;    // number of ligands defined in autoglobal.h
extern int ntor_lig[MAX_LIGANDS];
extern int gene_index_lig[MAX_LIGANDS][2];  //gene num start_point & end_point of a ligand.
extern int global_ntor; // set to current s.Init.ntor in main.cc

 
/*******************************************************************
 * Particle Swarms Optimization (PSO) implementation
 * 
 * Input Parameters: sInit the initial State of ligands
 * 
 * 
 * 
 * 
 * Reference:   cPSO, Vigneshwaran Namasivayam & Robert Guenther
 * 				University of Leipzig
 * 
 * Author: 		Huameng Li		04/26/2008
 * 		   		The Ohio State University 
 *******************************************************************/ 
State call_cpso(
	     int n_runs, 
	     State sInit, 
	     int swarm_size, 
	     int D, 
	     double *xmin, 
	     double *xmax, 
	     long eval_max, 
	     int K, 				//max links for each particle
	     double c1, 
	     double c2, 
	     int outlev,
		 unsigned int max_its, 
		 unsigned int max_succ, 
		 unsigned int max_fail, 
		 float expansion, 
		 float contraction, 
		 float search_freq, 
		 float *rho_ptr,
	     float *lb_rho_ptr,
	     float rho,
	     float lb_rho)
{

	int n_swarm_move = 0;

	Position Xi[S_max];    // array of particle positions
    Velocity Vi[S_max];    // array of particle Velocity
    Position P[S_max];	   // array of peronal Best of positions
    State sNew[S_max];
    State sTemp;
    State sbNew;  			// State variable for the best particle in current swarm Xi[best]
    int sb = 0;

	int i; 					// particle i index
	int j, k;
	int d; 
	//int m;  
	int best = 0;
	int gbest = 0;	
    double pso_energy;
    double energy_prev;
    double energyLocalSearch;      // Energy after local search
    
    //double min = 0.0;
	int init_links;
	int LINKS[S_max][S_max];
						 
	double prev_value[S_max];

	long int evaluations = 0;
	/*
	float wmax = 0.9;
	float wmin = 0.4;
	float w = 0.0;
	float ratio;
    */
	double Vmin[D_max], Vmax[D_max];
	
	//Constants for the Constriction PSO 
	double chi, phi;
	
	int ntor = sInit.ntor;	
	evaluate.reset();


	// Calculatin chi for CPSO
	phi = c1 + c2;
	
	//Use formula from Clerc, M. Kennedy, J., Evolutionary
	//Computation, IEEE Transactions on , vol. 6, no. 1, 58-73 (2002).	
	chi = 2.0 / (2.0 - phi - sqrt( (phi * phi) - (4 * phi) ));
	//chi = 2.0 / (phi + sqrt( (phi * phi) - (4 * phi) ));
	chi = fabs(chi);
	fprintf(logFile, "\nConstants used for the Constriction PSO \nc1: %lf c2: %lf phi: %lf chi: %lf\n", c1, c2, phi, chi);
	fflush(logFile);	
	
    ///////////////////////////////////////////////////
	// Initialization of the Particles
	////////////////////////////////////////////////////
	for( i = 0; i < swarm_size ; i++) {
				
		Xi[i].size = D; 
		Vi[i].size = D;
			
		///
		//Initialise the Particles 
		///
		//initialize_swarm_particle(i, D, Xi, Vi, xmin, xmax, Vmin, Vmax);	
		initialiseParticle(i, D, Xi, Vi, xmin, xmax, Vmin, Vmax);
		       
		sNew[i]=sInit;
		if(i == 0) {
			fprintf(logFile, "sNew[i] nlig=%d\n", sNew[i].nlig);
	    	fflush(logFile);
		}		
		
		//Transform position vector to STATE for E eval
		make_state_from_position(&sNew[i], Xi[i]);
		    	
		//Evaluate fitness function (binding Energy)
		Xi[i].f = evaluate.evalpso(&sNew[i]); 
		
		if (outlev >= 2) {
			pr(logFile, "\n(Initialization) Particle %d has Binding_E= %8.2lf\n", i, Xi[i].f);
			printState(logFile, sNew[i], 2);
			pr(logFile, "\n");
		}
		
		// set personal best positions so far   -Huameng				
		P[i] = Xi[i];// Best Postion = current One
					
		//store the current energy value of particle 
		prev_value[i] = BIG;
		//prev_value[i] = Xi[i].f;
		
		///
		// Find the best Energy of all particles afer intialization
		///
		if( P[i].f < P[best].f )	
			best = i;   		
							
	}
	
	
	// best fit Energy after intialization
	pso_energy = P[best].f;	
		
	energy_prev = pso_energy;


	//debug
	fprintf(logFile, "Done PSO Initialization of %d particles for run %d. lowestEnergy=%8.2lf\n\n", swarm_size, n_runs, P[best].f);		
	fflush(logFile);
	
	
	
	/////////////////////////////////////////////////////
	//
	// Start CPSO search	
	//
	/////////////////////////////////////////////////////
	init_links = 1;
	n_swarm_move=0;	
	
	do {
		/////////////////////////////////////
		//initialize links among particles
		/////////////////////////////////////
		if (init_links == 1) {
			//Who informs who, at random
			for (i = 0; i < swarm_size; i++) {
				
				// initialize 0 links to all other particles ( 0 neighbours)
				for( j = 0; j < swarm_size; j++) {
					if(i == j)
						LINKS[i][j] = 1; // Each particle informs itself
					else
						LINKS[i][j] = 0; // Init to "no link"
				}
				
				// set random K links (neigbours) for each particle
				j = 0;
				for(k = 0; k < K; k++) {
					j = (int) random_range(0, swarm_size -1);
					LINKS[i][j] = 1;
				}				
			}										
		} // end initialization of links
		
		
		/////////////////////////////////
		// Particle Swarm Moves
		/////////////////////////////////
		sb = 0;
		for (i = 0; i < swarm_size; i++) {		
			//Updating evaluations
			evaluations++;
										
			// find group best, gBest in neighborhood   -Huameng li				
			gbest = i;												
			for (j = 0; j < swarm_size; j++) {
				if ( LINKS[i][j] == 1 && P[j].f < P[gbest].f ) {
					gbest = j;
				}
			}
			///////////////////////////////////////////											
			//Update velocity and position
			///////////////////////////////////////////
			
			//Implement time-varying acceleration coefficients c1, c2
			float ratio = (float)evaluate.evals()/eval_max;
			//float c1_w = c1 + 8*ratio;    	// c1 increases with time
			//float c2_w = c2 * ratio;   // c2 deincreases with time
			//float w = 0.9 - (0.9 - 0.4) * ratio;
			for ( d = 0; d < D; d++ ) {											
				// -----Constriction -PSO ------Begin-								
			    if (Xi[i].f >= prev_value[i]) {
			    	// Relaxation Velocity Update
			    	if(	ratio < 1.0/4) {
			    		Vi[i].v[d] =  chi * (Vi[i].v[d] + c1 * random_range(0,1) * (P[i].x[d] - Xi[i].x[d]) 
								  + 0.2 * c2 * random_range(0,1) * (P[gbest].x[d] - Xi[i].x[d]));
			    	} else {	    			    	
						Vi[i].v[d] =  chi * (Vi[i].v[d] + c1 * random_range(0,1) * (P[i].x[d] - Xi[i].x[d]) 
								  + c2 * random_range(0,1) * (P[gbest].x[d] - Xi[i].x[d]));	
			    	}														
				}
														
				if(Vi[i].v[d] > Vmax[d]) {	
					Vi[i].v[d] = Vmax[d]; 
				} else if(Vi[i].v[d] < Vmin[d]) {
					Vi[i].v[d] = Vmin[d]; 
				}						
							
				Xi[i].x[d] = Xi[i].x[d] + Vi[i].v[d];
				
				//Interval confinement (keep in the gridInfo box  x, y, z )
				if (xmax[d] == PI) {
					// qw and torsions
					if (Xi[i].x[d] > PI) {
						//printf("\n Value Initial = %lf", Xi[i].x[d]);
						Xi[i].x[d] = Xi[i].x[d] - ((int)(Xi[i].x[d] / PI) * PI);
						//Vi[i].v[d]= 0.0 ;
						Vi[i].v[d]= random_range(Vmin[d], Vmax[d]) ; 						
					}
					
					if (Xi[i].x[d] < -PI) {
						//Xi[i].x[d] %= -PI;
						Xi[i].x[d] = Xi[i].x[d] - ((int)(Xi[i].x[d] / -PI) * -PI);
						//Vi[i].v[d]= 0.0 ;	
						Vi[i].v[d]= random_range(Vmin[d], Vmax[d]) ; 										  
					}					
				} else {				
					// reset x, y, z, qx,qy,qz for x[d] and v[d] 
					if( (Xi[i].x[d] < xmin[d]) || (Xi[i].x[d] > xmax[d]) ) {
						Xi[i].x[d]= random_range(xmin[d], xmax[d]); 
						Vi[i].v[d]= random_range(Vmin[d], Vmax[d]) ; 
					}
				}		
									
								
			} // -----Constriction -PSO ------End---

					
			//copy the dimension of Position to State Variable 
			make_state_from_position(&sNew[i], Xi[i]);
			
			//...evaluate the new position
			prev_value[i] = Xi[i].f;
			Xi[i].f = evaluate.evalpso(&sNew[i]); // compute bindng energy 
			
			if (outlev > 2) {
				pr(logFile,"\nSwarmMove: (%d) \tParticle:  %d \tEnergy= %8.2lf\n", n_swarm_move + 1, i + 1, Xi[i].f);
				printState(logFile, sNew[i], 0);
				fflush(logFile);
			}
			
			// get the Best of all current Xi[i]
			if( Xi[i].f <= Xi[sb].f ) 
				sb = i;
						
		} // end i swarm_size particles
		
		//fprintf(logFile, "Done swarm movement %d.\n, n_swarm_move");		
		//fflush(logFile);
		
		//Have to initialize sbNew  -Huameng 03/08/08
		sbNew = sNew[sb];    
		make_state_from_position(&sbNew, Xi[sb]);
					
		/////////////////////////////
		// do Solis Wets local search for the best of current swarm
		/////////////////////////////				
		if (ranf() < search_freq) {
			//fprintf(logFile, "Current swarm move, BEFORE LS gBest E = %8.2lf, ", Xi[sb].f);
			//fflush(logFile);
				
			//SW local search all ligands together			
			int	from_idx; 
			int to_idx;	   									   	
			for(int n = 0; n < nlig ; n++) {
		   		//get the segment if gene points for this ligand 	   	   
		   	  	from_idx = gene_index_lig[n][0];	   
		   	  	to_idx = gene_index_lig[n][1];
		   	  	//SWLocalSearch(sbNew, n, from_idx, to_idx, ntor, max_its, max_succ, max_fail, 2.0, 0.5, rho, lb_rho);			  	     	   	 		   	    		
		      	PSWLocalSearch(sbNew, ntor, max_its, max_succ, max_fail, i, from_idx, to_idx, rho_ptr, lb_rho_ptr);  				   		
		   	}	 	   	   	
			// all torsion
			from_idx = nlig*7;
			to_idx = nlig*7 + global_ntor;
			//SWLocalSearch(sbNew, nlig, from_idx, to_idx, ntor, max_its, max_succ, max_fail, 2.0, 0.5, rho, lb_rho);		  
			PSWLocalSearch(sbNew, ntor, max_its, max_succ, max_fail, nlig, from_idx, to_idx, rho_ptr, lb_rho_ptr);	   		   	   				
			
			// check if we get lower binding energy
			energyLocalSearch = evaluate.evalpso(&sbNew);
			if(energyLocalSearch < Xi[sb].f ) {			
				copy_state_to_position(&Xi[sb], sbNew);			
				Xi[sb].f = energyLocalSearch; //Energy after local search			
			}					
			//fprintf(logFile, " AFTER LS local Best E = %8.2lf\n", energyLocalSearch);				
		}							 			
		//////////////////////////////////////////
		// update pBest and gBest (Pi and Gi)
		/////////////////////////////////////////
		for (i = 0; i < swarm_size; i++) {
			//... Update the personal best so far i, e.g., pBest
			if (Xi[i].f < P[i].f) {
				P[i] = Xi[i];
				//...update the best of the bests, final global Best
				if(P[i].f < P[best].f) 
					best = i;
			}		
		}	// end for i
		
		n_swarm_move++;
		
#ifdef DEBUG
		//the code seems just for debuggging and does not improve the E    -Huameng		
		//Swarm - Begin         Activity moved as a function 		
		//swarmActivity(S, D, Xi, n_swarm_move, outlev);
		//Swarm  - End
#endif			
		
		pso_energy = P[best].f;
		
		// If no improvement, information links will be reinitialised
		if(pso_energy >= energy_prev) 
			init_links = 1;
		else 
			init_links = 0;
			
		energy_prev = pso_energy;
		if(n_swarm_move%100 == 0) {
			pr(logFile, "Swarm move %d, PSO gbest E= %8.2lf** (num_evals=%ld)\n", n_swarm_move,  pso_energy,  evaluate.evals());					
		}
		
		if (outlev > 1 ){
			pr(logFile, "PSO: Updating the swarm at move %d (= %ld and %ld evaluations)\n", n_swarm_move+1, evaluations, evaluate.evals());
		}				
		
	} while(evaluate.evals() < eval_max);  // max evaluations
   
 
	//Result is stored in the History state for writePDBQ and for clustering ....
	sTemp = sInit;
	make_state_from_position( &sTemp, P[best]);
		
	pso_energy = evaluate.evalpso(&sTemp);
	pr(logFile, "\n\nFinal State:  for run %d  E = %8.2lf \n", n_runs + 1 , pso_energy );			
	//printState(logFile, sTemp, 2);	
		
	return (sTemp);
	
} 


/**
 *  The Position of N demension for a particle is the representaion of a ligand's state
 *  P.x[1..Nd] has following sectors:
 *  	x([0..2])n is the translation
 *  	x([3..6])n is rotation
 *  	x([7..ntor])n is the torsion section
 */
void make_state_from_position( State *S, Position P)
{
    register int i, j;
    int idx;
    
    //handle multi-ligand  -Huameng
    for(i = 0; i < nlig; i++) {
    	idx = 7*i;
    	//translation
	    S->T[i].x = P.x[idx + 0];
	    S->T[i].y = P.x[idx + 1];
	    S->T[i].z = P.x[idx + 2];
	    // quarterion
	    /*    	    
	    S->Q[i].x = P.x[idx + 3];
	    S->Q[i].y = P.x[idx + 4];
	    S->Q[i].z = P.x[idx + 5];
	    S->Q[i].w = P.x[idx + 6];
	    */  	
	   // rotation	 
	   S->Q[i].nx = P.x[idx + 3];
	   S->Q[i].ny = P.x[idx + 4];
	   S->Q[i].nz = P.x[idx + 5];
	   S->Q[i].ang = P.x[idx + 6];		   
	      	         		   
    }
    
    // all torsion part
    for(j = 7*nlig, i = 0; i < S->ntor; i++, j++) {
		S->tor[i] = P.x[j];
	}	
}


void copy_state_to_position(Position *R , State S)
{
    register int i, j, idx;
    //handle multi-ligand  -Huameng
    for(i = 0; i < nlig; i++) {
    	idx = 7*i;  	
    		 	
		R->x[idx + 0] = S.T[i].x;
		R->x[idx + 1] = S.T[i].y;
	    R->x[idx + 2] = S.T[i].z;
	        
	    // quarterion
	    /*	   
	    R->x[idx + 3] = S.Q[i].x;
	    R->x[idx + 4] = S.Q[i].y;
	    R->x[idx + 5] = S.Q[i].z;
	    R->x[idx + 6] = S.Q[i].w;
	   	*/
	   	// rotation
	    R->x[idx + 3] = S.Q[i].nx;
	    R->x[idx + 4] = S.Q[i].ny;
	    R->x[idx + 5] = S.Q[i].nz;
	    R->x[idx + 6] = S.Q[i].ang;
	       	 	        
    }
    
    for(j= 7*nlig, i=0; i < S.ntor; i++, j++) {
		R->x[j] = S.tor[i];
    }
}


/**
 * Initialize the min and max of each dimensions of position/state variable 
 */
void initialiseDimension(GridMapSetInfo *info, double *xmin, double *xmax, int D)
{
	int i, d, offset;
			
	//handle multi-ligand -Huameng
	for(i = 0; i < nlig; i++)
	{								
		// 0, 1, 2 (x,y,z) for translation
		offset = 7*i;
		xmin[offset + 0] = info->lo[X];
    	xmax[offset + 0] = info->hi[X];
    						
		xmin[offset + 1] = info->lo[Y];
    	xmax[offset + 1] = info->hi[Y];	
    			
		xmin[offset + 2] = info->lo[Z];
    	xmax[offset + 2] = info->hi[Z];		
		fprintf(logFile, "In initialise X: d=%d  xmin=%.3f	xmax=%.3f\n", offset + 0, xmin[offset + 0], xmax[offset + 0]);
		fprintf(logFile, "In initialise Y: d=%d  xmin=%.3f	xmax=%.3f\n", offset + 1, xmin[offset + 1], xmax[offset + 1]);		
		fprintf(logFile, "In initialise Z: d=%d  xmin=%.3f	xmax=%.3f\n", offset + 2, xmin[offset + 2], xmax[offset + 2]);
		fflush(logFile);
					
		// 7*i + 3,4,5,6 are for orientation
		xmin[offset + 3] = 0;
    	xmax[offset + 3] = 1.0;					
		xmin[offset + 4] = 0;
    	xmax[offset + 4] = 1.0;			
		xmin[offset + 5] = 0;
    	xmax[offset + 5] = 1.0;		    	
    	xmin[offset + 6] = -PI;
    	xmax[offset + 6] = PI;
#ifdef DEBUG
		fprintf(logFile, "In initialise Q: d=%d  xmin=%.3f	xmax=%.3f\n", offset + 3, xmin[offset + 3], xmax[offset + 3]);
		fprintf(logFile, "In initialise Q: d=%d  xmin=%.3f	xmax=%.3f\n", offset + 4, xmin[offset + 4], xmax[offset + 4]);		
		fprintf(logFile, "In initialise Q: d=%d  xmin=%.3f	xmax=%.3f\n", offset + 5, xmin[offset + 5], xmax[offset + 5]);			    	
	    fprintf(logFile, "In initialise Q: d=%d  xmin=%.3f	xmax=%.3f\n", offset + 6, xmin[offset + 6], xmax[offset + 6]);
		fflush(logFile);
#endif					    												 		    
	}
	
	// the rest are for torsion
	for(d = 7*nlig; d < D; d++)
	{
		xmin[d] = -PI;
	    xmax[d] = PI;	    
	}
	
	fprintf(logFile, "End initialiseDimension() of PSO \n");
	fflush(logFile);
	/*
    for (d=0; d < D ; d++)
    {
    	switch(d)
        {
	        case 0:
				xmin[d] = xlo;
	            xmax[d] = xhi;
				break;
			case 1:
				xmin[d] = ylo;
				xmax[d] = yhi;
				break;
			case 2:
				xmin[d] = zlo;
				xmax[d] = zhi;
				break;
			case 3: // For Quaternion angles
				xmin[d] = 0;
				xmax[d] = 1;
				break;
			case 4:
				xmin[d] = 0;
				xmax[d] = 1;
				break;
			case 5:
				xmin[d] = 0;
				xmax[d] = 1;
				break;				
			case 6:
				xmin[d] = -PI;
				xmax[d] = PI;
				break;
			default: // For Torsional Angles
				xmin[d] = -PI;
				xmax[d] = PI;
				break;		
		}
	}
	*/
}


void initialiseParticle(int s, int D, Position *Xi, Velocity *Vi, double *xmin, 
						double *xmax, double *Vmin, double *Vmax)
{
    int d = 0, i, j;
	double temp;
	
	// handle multi-ligand -Huameng	02/24/08
	for(i = 0; i < nlig; i++) {
		for(j = 0; j < 7; j++) {			
			if( j <= 2 ) {
				temp = xmin[d] - xmax[d];
				Vmax[d] = fabs(temp) / 2;
				Vmin[d] = -Vmax[d];
			} else if( j > 2 && j < 6) {
				temp = xmin[d] - xmax[d];
				Vmax[d] = fabs(temp);
				Vmin[d] = -Vmax[d];				
			} else {				
				Vmin[d] = -PI;
				Vmax[d] = PI;				
			}
			Xi[s].x[d] = random_range(xmin[d], xmax[d]);
			Xi[s].prev_x[d] = Xi[s].x[d];
			Vi[s].v[d] = random_range(Vmin[d], Vmax[d] ); 
			d++;
#ifdef DEBUG			
			pr(logFile, "Particle first 7 Dim: Xi[%d].x[%d]=%.3f\n", s, d, Xi[s].x[d]);	
			fflush(logFile);			
#endif																
		}						
	}
		
	// torsion part
	for( d = 7*nlig; d < D; d++ )
	{
		Vmin[d] = -PI;
		Vmax[d] = PI;
				
		Xi[s].x[d] = random_range(xmin[d], xmax[d]);
		Xi[s].prev_x[d] = Xi[s].x[d];
		Vi[s].v[d] = random_range(Vmin[d], Vmax[d] ); 							
	}
	//pr(logFile, "initialise torsions  Xi[%d].x[%d]=%.3f\n", s, d-1, Xi[s].x[d-1]);				
}

/*********************************************************************
 * Generate a swarm of particles with random position Xi
 * Xi has translation, quarterinon and torsion as D dimensions
 * in format: (x, y, z, qx, qy, qz, qw)*nligand + tor1,tor2,...torN
 * 
 * Huameng Li OSU   02/26/2008
 * 
 ********************************************************************/
void  initialize_random_particle ( 
 			int i,    			// particle i
 			int D,    			// degree of freedom
 			Position *Xi, 		// pointer to Positions
 			Velocity *Vi, 		// pointer to Velocity of particles
 			double *xmin,		// Min. x,y,z of gridInfoMap
 			double *xmax,     // Max. x,y,z of gridInfoMap
 			double *Vmin,  
 			double *Vmax )
{	
 	int d = 0, n, j;
	double temp;
	Quat q;
	
	// handle multi-ligand -Huameng	02/24/08
	for(n = 0; n < nlig; n++) {
		for(j = 0; j < 7; j++) {			
			d = 7*n + j;
			//pr(logFile, "Particle %d Dimension %d: \n", s, d);
			//fflush(logFile);			
			if( j < 3) {
				//Xi[i].x[d] = random_range(xmin[d], xmax[d]);
				Xi[i].x[d] = double(genunf(xmin[d], xmax[d]));
				temp = xmin[d] - xmax[d];
				Vmax[d] = fabs(temp) / 2;
				Vmin[d] = -Vmax[d];
			} else { 			
				// Generate a uniformly-distributed random quaternion for a random rotation (UDQ)
		   		q = uniformQuat();
		   		q = convertQuatToRot( q );      //return rotation
		   				   		
		   		if(j == 3)
					Xi[i].x[d] = q.nx;
				if(j == 4) 
					Xi[i].x[d] = q.ny;
				if(j == 5)
					Xi[i].x[d] = q.nz;
				if(j == 6)
					Xi[i].x[d] = q.ang;
				
				if( j < 6) {
					Vmax[d] = fabs(xmin[d] - xmax[d]);
					Vmin[d] = -Vmax[d];
				} else {
					Vmax[d] = PI;
					Vmin[d] = -PI;
				}
				Xi[i].prev_x[d] = Xi[i].x[d];					
				Vi[i].v[d] = random_range( Vmin[d], Vmax[d] );		
							
				//if(j == 3)
				//	Xi[i].x[d] = q.x;
				//if(j == 4) 
				//	Xi[i].x[d] = q.y;
				//if (j == 5)
				//	Xi[i].x[d] = q.z;
				//if (j == 6)
				//	Xi[i].x[d] = q.w;
				//Vi[i].v[d] = double(genunf(xmin[d], xmax[d])); 					
				//Vmax[d] = 1.0;
				//Vmin[d] = -1.0;																																
			}															
#ifdef DEBUG			
			pr(logFile, "Particle %d Dimension %d: Xi[i].x[d]=%.3f  Vi[i].v[d] = %.3f\n", i, d, Xi[i].x[d], Vi[i].v[d]);	
#endif																
		} // end for 					
	} // end for ligands
	
	// torsion part
	for( d = 7*nlig; d < D; d++ ) {
		Vmax[d] = PI;
		Vmin[d] = -PI;
		
		//Xi[i].x[d] = random_range(xmin[d], xmax[d]);
		Xi[i].x[d] = double(genunf(-PI, PI));		
		Xi[i].prev_x[d] = Xi[i].x[d];
		//Vi[i].v[d] = random_range(Vmin[d], Vmax[d] );
		Vi[i].v[d] =  double(genunf(-PI, PI));  							
	}
	//pr(logFile, "initialise Particle: last dimension Xi[%d].x[%d]=%.3f\n", s, d-1, Xi[i].x[d-1]);			
 }

			
void swarmActivity(int S, int D, Position *Xi, int nb_eval, int outlev)
{
	int s, d;
	double swarm_activity;
	double position_dist[S_max];    // Euclidian distance of position and prev position
	double diff;                    // just a helper variable for calculating the distance
	// RG Swarm activity begin
    swarm_activity = 0;
    for(s=0; s < S; s++)
	{
    	for(d = 0; d < D; d++)
        {
			diff = Xi[s].x[d] - Xi[s].prev_x[d];
			//pr(logFile, "swarm_move %d particle= %d dim= %d diff= %f\n", nb_eval, s, d, diff);
			// Differences of angles
			if (d > 5)
			{
				if (diff < -PI) 
					diff = diff + PI;
				if (diff > PI)  
					diff = diff - PI;
				
			}
			position_dist[s] += (diff * diff);
		}
		position_dist[s] = sqrt(position_dist[s]);
		swarm_activity += position_dist[s];
	}
	
	if (outlev >1)
	{
		swarm_activity = (swarm_activity) / (S * D);
		pr(logFile, "swarm_move %d SA= %f\n", nb_eval+1, swarm_activity);
	}
	// RG Swarm activity end
}


/**
 * 
 * SW local search for multi-ligand PSO
 *  
 * Solis & Wets algorithm adds random deviates to every real number in the STATE variable.
 * SWLocalSearch based on Solis & Wets local search 
 */
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
	) 
{
	float bias[7*nlig + ntor];
	float deviates[7*nlig + ntor];
	
	State newstate;
	float temp_rho = rho;
	
	
	int i, j, num_successes = 0, num_failures = 0;
		
	//int size= (unsigned)ntor+7;
	//int size = (unsigned) 7*nlig + ntor ;
	
	//Have to copy best state to tmpstate, otherwise, it causes segment fault -Huameng		
	
	//  Reset bias	   
    //for (i=0; i < size; i++) {	  
    //  bias[i] = 0.0;
    //}
    //  Reset bias
	for (i= from_idx; i < to_idx; i++) {
	    bias[i] = 0.0;
	}
	for (i=0; i < max_its; i++) {	  	
	    // Generate deviates
	    for (j=from_idx; j < to_idx; j++) {
	      deviates[j] = gennor(0.0, temp_rho);
	    }
	    
	    // zeta = x + bias + deviates
	 	newstate = addDeviatesBiasAD4(state, deviates, bias, ilig, +1.0);
	 		 
	    if (evaluate.evalpso(&newstate) < evaluate.evalpso(&state)) {
	      	num_successes++;
	      	num_failures = 0;
	      	//memcpy(&state, &newstate, sizeof(State));  
	        copyState(&state, newstate);		 	
	      	for (j=from_idx; j < to_idx; j++) {
				//bias[j] = 0.20*bias[j] + 0.40*deviates[j];
				bias[j] = 0.60*bias[j] + 0.40*deviates[j];   // strict SW
	      	}
	    } else {
	    	
		    newstate = addDeviatesBiasAD4(state,  deviates, bias, ilig, -1.0);
	 		    
		    if (evaluate.evalpso(&newstate) < evaluate.evalpso(&state)) {
				num_successes++;
				num_failures = 0;
				//memcpy(&state, &newstate, sizeof(State));
				copyState(&state, newstate);
				for (j=from_idx; j < to_idx; j++) {
			  		//bias[j] -= 0.40*deviates[j];
			  		bias[j] = 0.60*bias[j] - 0.40*deviates[j];   // strict SW
				}
	      	} else {
				num_failures++;
				num_successes = 0;
				for (j=from_idx; j < to_idx; j++) {
				  	bias[j] *= 0.50;
				}
	      	}
	    }
	        	 	
	    // Check to see if we need to expand or contract
      	if (num_successes >= max_succ) {
         	temp_rho *= expansion;
         	num_successes = num_failures = 0;
      	} else if (num_failures >= max_fail) {
         	temp_rho *= contraction;
         	num_successes = num_failures = 0;
      	}
	    
	   if (temp_rho < lb_rho);
	    	break;
	    
	} //  i-loop
	  //pr(logFile, "lowest E after SW local search = %.3f\n", evaluate.evalpso(state));
	  //fflush(logFile);
	 
	  return (state);
} // void Solis_Wets::SW(Phenotype &vector)



/**
 * 
 * pseudo-SW local search for multi-ligand PSO
 *  
 * pseudo-Solis & Wets algorithm adds random deviates to every real number in the STATE variable.
 * SWLocalSearch based on Solis & Wets local search 
 */
State PSWLocalSearch(
	State &state, 
	int ntor, 
	 int max_its, 
	 int max_succ, 
	 int max_fail,
	 int ilig,		
     int from_idx,
     int to_idx,	
	float *rho_ptr, 
	float *lb_rho_ptr) 
{	
	float bias[7*nlig + ntor];
	float deviates[7*nlig + ntor];
	
	State newstate;
	float temp_rho[7*nlig + ntor];
	
	float expansion = 2.0;
	float contraction = 0.5;
	int i, j, num_successes = 0, num_failures = 0,  all_rho_stepsizes_too_small = 1;
		
	//int size = (unsigned) 7*nlig + ntor ; // toal dimensions for all ligands		
	  
	//  Initialize the temp_rho's
	//for (i = 0; i < size; i++) {
	for (i = from_idx; i < to_idx; i++) {
	    temp_rho[i] = rho_ptr[i];
	}
	   
	//  Reset bias
	for (i= from_idx; i < to_idx; i++) {
	    bias[i] = 0.0;
	}
		
	for (i=0; i < max_its; i++) {	  	
	    // Generate deviates
	    for (j= from_idx; j < to_idx; j++) {
	      deviates[j] = gennor(0.0, temp_rho[j]);
	    }
	    	    
	 	newstate = addDeviatesBiasAD4(state, deviates, bias, ilig, 1.0);
	 		 
	    if (evaluate.evalpso(&newstate) < evaluate.evalpso(&state)) {
	      	num_successes++;
	      	num_failures = 0;
	      	//memcpy(state, &newstate, sizeof(State));  
	        copyState(&state, newstate);
	        		 	
	      	for (j= from_idx; j < to_idx; j++) {
				//bias[j] = 0.20*bias[j] + 0.40*deviates[j];
				bias[j] = 0.60*bias[j] + 0.40*deviates[j];   // strict SW
	      	}
	    } else {
	      	//for (j=0; j < size; j++) {
			//	deviates[j] *= -1.0;
	      	//}	      			    
		    newstate = addDeviatesBiasAD4(state, deviates, bias, ilig, -1.0);	    	 	  
		    if (evaluate.evalpso(&newstate) < evaluate.evalpso(&state)) {
				num_successes++;
				num_failures = 0;
				//memcpy(state, &newstate, sizeof(State));
				copyState(&state, newstate);
				
				for (j = from_idx; j < to_idx; j++) {
			  		//bias[j] -= 0.40*deviates[j];
			  		bias[j] = 0.60*bias[j] - 0.40*deviates[j];   // strict SW
				}
	      	} else {
				num_failures++;
				num_successes = 0;
				for (j = from_idx; j < to_idx; j++) {
				  	bias[j] *= 0.50;
				}
	      	}
	    }
	        	 	
	    // Check to see if we need to expand or contract
	    if (num_successes >= max_succ) {
	      for(j = from_idx; j < to_idx; j++) {
			temp_rho[j] *= expansion;
	      }
	      num_successes = num_failures = 0;
	    } else if (num_failures >= max_fail) {
	      for(j = from_idx; j < to_idx; j++) {
			temp_rho[j] *= contraction;
	      }
	      num_successes = num_failures = 0;
	    }
	    
	    //  WEH - Scott's code doesn't do anything!!! no stopping based upon step scale!!!
	    //  GMM - corrected Scott's code; this does now stop correctly, based upon step scale.
	    //  GMM - This version only exits if all the step sizes are too small...	    
	    for(j = from_idx; j < to_idx; j++) {   
	      all_rho_stepsizes_too_small = all_rho_stepsizes_too_small & (temp_rho[j] < lb_rho_ptr[j]);
	    } //  j-loop
	    if (all_rho_stepsizes_too_small) {
	      break; //  GMM - THIS breaks out of i loop, which IS what we want...
	    }
	} //  i-loop
	//pr(logFile, "lowest E after SW local search = %.2f\n", evaluate.evalpso(state));
	
	return (state);
} // void Pseudo_Solis_Wets::SW(Phenotype &vector)


/**
 * -Huameng Li 03/08/08
 */
State addDeviatesBiasAD4(State &srcstate,  float* deviates, float* bias, int i, double sign) 
{	
	State deststate;
	
	deststate.nlig = srcstate.nlig;
	deststate.ntor = srcstate.ntor;
	
	// handle multi-ligand -huameng 02/23/08
	int idx = 7*i;   //idx of the start dimension for this ligand
	if(idx < 7*nlig) { 	
		// tanslation for ligand i			  		
		deststate.T[i].x = srcstate.T[i].x + sign * (deviates[idx + 0]+ bias[idx + 0]);
		deststate.T[i].y = srcstate.T[i].y + sign * (deviates[idx + 1]+ bias[idx + 1]);
		deststate.T[i].z = srcstate.T[i].z + sign * (deviates[idx + 2]+ bias[idx + 2]);
		
		// rotation angle			
		deststate.Q[i].nx = srcstate.Q[i].nx + sign * (deviates[idx + 3]+ bias[idx + 3]);
		deststate.Q[i].ny = srcstate.Q[i].ny + sign * (deviates[idx + 4] + bias[idx + 4]);
		deststate.Q[i].nz = srcstate.Q[i].nz + sign * (deviates[idx + 5]+ bias[idx + 5]);
		deststate.Q[i].ang = srcstate.Q[i].ang + sign * (deviates[idx + 6]+ bias[idx + 6]);
			
		//use quaternion in stead of rotation
		/*
		deststate.Q[i].x = srcstate.Q[i].x +  (deviates[idx + 3]+ bias[idx + 3]);
		deststate.Q[i].y = srcstate.Q[i].y +  (deviates[idx + 4] + bias[idx + 4]);
		deststate.Q[i].z = srcstate.Q[i].z +  (deviates[idx + 5]+ bias[idx + 5]);
		deststate.Q[i].w = srcstate.Q[i].w + (deviates[idx + 6]+ bias[idx + 6]);
		*/
		//pr(logFile, "\nIn rotation: deviates=%.3f\t bias=%.3f, sign=%.2f\n", deviates[idx+6], bias[idx+6], sign);			  
		//pr(logFile, "index=%d, destDateQ->ang = %.3f   srcstate->Q[i].ang = %.3f\n", idx+6, deststate->Q[i].ang, srcstate->Q[i].ang);
		//fflush(logFile);
	} else {
		//torsion part			  
		for (int j=0; j < srcstate.ntor; j++) {
		    deststate.tor[j]= srcstate.tor[j] + sign *(deviates[7*nlig + j] + bias[7*nlig + j]);		    
		    //pr(logFile, "\nIn torsion:  deviates[7*nlig + j] =%.3f,  bias[7*nlig + j]=%.3f\n", deviates[7*nlig + j], bias[7*nlig + j]);
		    //pr(logFile, "torsion j=%d, deststate->tor[j] = %.3f   srcstate->tor[j] = %.3f\n", j, deststate->tor[j], srcstate->tor[j]);
		    //fflush(logFile);
		}
	}
	
	return deststate;
	
}

