#include "pso.h"
#include "ranlib.h"
#include <math.h>

extern FILE *logFile;
extern Eval evaluate;


inline float Norm(float *x, int n)
{
	float s = 0;
	for(int i=0;i<n;i++)
		s += (x[i]*x[i]);
	return sqrt(s);
}


/***********************************************************************
 * Particle Swarm Optimization (PSO) with time varying inertia weight w.
 * Global PSO and local S&W search.
 * 
 * Based on PSO in SODOCK: J. Comput Chem, 28, 2, 612-623, 2007.
 * 
 * Huameng Li, 07/2008
 * 
 *
 * Key modifications spring-summer 2011 by R Huey & M Pique at TRSI:
 *   No longer does local search here, the caller will do that
 *     after each generation using our localsearch method, below.
 * Caution: do not reorder the Pop array outside of these functions. MP/RH
 ***********************************************************************/
int ParticleSwarmGS::search(Population &Pop)
{
	int i, j, g;
	double r1, r2;
	double curVal;
	double newVal;
	double piCurE = 0.0;
	//float phi, c;
	float ratio = 0.0;
	
	// on first call per run, allocate Pi array, initialize velocity vectors
	if(_Pi == NULL) {
		
		pr(logFile, "PSO max_generations = %d\n", max_generations);
		pr(logFile, "PSO size = %d\n", size); // MP debug
				
		fprintf(logFile, "PSO Allocate initial velocity of particles...\n");
		_Pi = new Population(Pop); // copy of Pop (? MP)
		
		pop_size = Pop.num_individuals();
		pr(logFile, "PSO pop_size = %d\n", pop_size); // MP debug
						
		best = 0;
		for(i = 1; i < pop_size; i++)
			if( (*_Pi)[i].value(Normal_Eval) < (*_Pi)[best].value(Normal_Eval) )
				best = i;
		_Pg = new Individual( (*_Pi)[best] );
		
		// allocate velocity 
		v = new float * [pop_size];
		for(i=0; i < pop_size; i++) {
			v[i] = new float [size];
		}
			
		// initial velocity of particles
		for(i = 0; i < pop_size; i++) {			
			for(j = 0; j < size; j++) {
				v[i][j] = random_range(vmin[j], vmax[j]);				
			}			
		}
		/* analysis
		pr(logFile, "Initial velocity V:\n");
		for(i=0;i<pop_size;i++) {
			pr(logFile, " V[%d] ", i);
			for(j=0;j<size;j++)
				pr(logFile, "%+6.2f ", v[i][j]);
			pr(logFile, "\n");
		}				
		pr(logFile, "Pi:\n");
		for(i=0;i<pop_size;i++) {
			pr(logFile, "Pi[%d] ", i);
			for(j=0;j<size;j++)
				pr(logFile, "%+6.2f ", (*_Pi)[i].phenotyp.gread(j).real);
			pr(logFile, "\n");
		}
		*/
								
		// Display field title
		//pr(logFile, "Generation NumEvals     Pg       PopAvg     PiAvg       w      |V| Avg\n");
		//pr(logFile, "---------- -------- ---------- ---------- ---------- -------- --------\n");					
		//pr(logFile, "Generation\tNumEvals\tpi X_best\tgBest Pg\tImproved\n");	
		pr(logFile, "Generation NumEvals  gBest E   \n");
		pr(logFile, "---------- -------- ---------- \n");					
									
	}

	// Update weights as the run progresses
	// Huameng had:  ratio = (float)evaluate.evals() / num_evals;
	ratio = (float)generations / max_generations;  // MP TODO since we dont have num_evals
	pso_w = wmin +(wmax - wmin) * (1.0 - ratio);
	
	// Set shorthand names for best for each indiv, best individual in pop
	Population &Pi = (Population &)(*_Pi);
	Individual &Pg = (Individual &)(*_Pg);
	
	//double piBestE[pop_size];	    // E array for personal Best	
	//double prevE[pop_size];       // E array for particles in previous run 
	//double curE[pop_size];        // E array for particles in current run

	// update Pi from Pop, set Global Best Pg, which points to _Pg
	for(i = 0; i < pop_size; i++ ) {		
		if( Pop[i].value(Normal_Eval) < Pi[i].value(Normal_Eval) )		
		        Pi[i] = Pop[i];
	 // assert energy of Pi[i] <= energy of Pop[i]
	}
	best = 0;	
	for(i = 1; i < pop_size; i++ ) {		
		if( Pi[i].value(Normal_Eval) < Pi[best].value(Normal_Eval) )		
			best = i;			
	}
	Pg = Pi[best];
	// assert energy of Pg <= energy of Pi[*] <= energy of Pop[*]
		   		   	 	     	   	 	   	   	   	  
	//update position variable related to the translation and rotation of each particle			
	for(i = 0; i < pop_size; i++) {
				
		//r2 = 1.0 - r1;
		//get the best in neighborhood, e.g., Neighborhood Best (nbBest)		
		g = i;
		for(j = i + 1; j < i + pso_K; j++) {
			if(Pi[j % pop_size].value(Normal_Eval) < Pi[g].value(Normal_Eval))
				g = (j % pop_size);
		}
				
		// size is the dimension of docking search (e.g, n*nlig + ntor)			
		for(j = 0; j < size; j++) {
			r1 = ranf();
			r2 = ranf();
			// note that phenotype "gread" is x,y,z,qx,qy,qz,qw,t1...
			curVal = Pop[i].phenotyp.gread(j).real;
									
			
			// -Huameng
			// update velocity
			if(ratio <= 0.8) {				
				v[i][j] = pso_w * v[i][j] + c1 * r1 * (Pi[i].phenotyp.gread(j).real - curVal)
			                      + c2 * r2 * (Pi[g].phenotyp.gread(j).real - curVal);
			} else {
				//Incorporate stage 2 constriction PSO 
				v[i][j] = 0.1 * (v[i][j] + 6.05 * r1 * (Pi[i].phenotyp.gread(j).real - curVal)
			                    + 6.05 * r2 * (Pi[g].phenotyp.gread(j).real - curVal));
			}
			                      																			
			// verify and restrict the new velocity
			if(v[i][j] > vmax[j])
				v[i][j] = vmax[j];
				//v[i][j] = random_range(vmin[j], vmax[j]) ;
			else if(v[i][j] < vmin[j])
				v[i][j] = vmin[j];
				//v[i][j] = random_range(vmin[j], vmax[j]) ;
				
			newVal = curVal + v[i][j];
										
			// update x,y,z, quaternion, torsion of particle i
			Pop[i].phenotyp.write(newVal, j);												
		}
		// must normalize Quaternion after modifying its components
		Quat q = Pop[i].phenotyp.readQuat();
		Pop[i].phenotyp.writeQuat( normQuat( q ));

		Pop[i].inverse_mapping(); // MP@@ copy phenotype to genotype (may not be necc)
		
	}	// end current swarm	
	
	// Update personal best Pi after this swarm search generation.
	// Find the best in Pop after this swarm search generation.
	// Note each could be better or worse than former best, stored in Pg
	// This 'best' value will be used by a subsequent call to LocalSearch,
	//   see below. M Pique June 2011
	best = 0;
	double X_best_value = 99999999;
	for(i = 0; i < pop_size; i++) {		
		piCurE = Pop[i].value(Normal_Eval);
		if(piCurE < X_best_value) {
			best = i;
			X_best_value = piCurE;
		}
		
		// Update Pi, personal Best in history				
		if(piCurE < Pi[i].value(Normal_Eval)) {
			Pi[i] = Pop[i];			
		}										
	}

#ifdef SWAPBEST
// MP June 2011 - unworkable try to unify GA/LGA and PSO
	// exchange the best Pop individual into location Pop[0]
	// and exchange the corresponding Pi history entry
	pr(logFile, "PSO pre-swap  best=%d Pop[%d]=%.2f X_best_value=%.2f, Pop[0,1,..] = ", 
	  best, best, Pop[best].value(Normal_Eval), X_best_value);
	for(i=0;i<8; i++) pr(logFile, "%8.2f ", Pop[i].value(Normal_Eval));
	pr(logFile, "\n");



	if(best!=0) {
		Individual temp;
		temp = Pop[0];
		Pop[0] = Pop[best];
		Pop[best] = temp;

		temp = Pi[0];
		Pi[0] = Pi[best];
		Pi[best] = temp;
	}

	// MP debug - look again for best after swap, should be at [0]
	// Find the best in Pop	
	best = 0;
	X_best_value = 99999999; // MP debug Pop[0].value(Normal_Eval);	
	for(i = 0; i < pop_size; i++) {
		piCurE = Pop[i].value(Normal_Eval);
		if(piCurE < X_best_value) {
			best = i;
			X_best_value = piCurE;
		}
		
		// Update Pi, personal Best in history				
		if(piCurE < Pi[i].value(Normal_Eval)) {
			Pi[i] = Pop[i];			
		}										
	}
	pr(logFile, "PSO post-swap best=%d Pop[%d]=%.2f X_best_value=%.2f, Pop[0,1,..] = ", 
	  best, best, Pop[best].value(Normal_Eval), X_best_value);
	for(i=0;i<8; i++) pr(logFile, "%8.2f ", Pop[i].value(Normal_Eval));
	pr(logFile, "\n");
#endif

	////////////////////////////////////////////////////////
	// Local Search
	////////////////////////////////////////////////////////

	generations++;
	
	// comments out the swarm analysis   -Huameng
	/********************************************************
	float Pop_avg = 0;
	float Pi_avg = 0;
	float v_avg = 0;
	for(i = 0; i < pop_size;i++) {
		Pop_avg += Pop[i].value(Normal_Eval);
		Pi_avg += Pi[i].value(Normal_Eval);
		v_avg += Norm(v[i], size);
	}
	Pop_avg /= pop_size;
	Pi_avg /= pop_size;
	v_avg /= pop_size;
		
	float diff[50];
	float dist = 0;
	int k;	
	int count = 0;
	for(i = 0; i < pop_size-1; i++) {
		for(j = i+1;j < pop_size; j++)
			for(k = 0; k < size; k++)
				diff[k] = Pop[i].phenotyp.gread(k).real - Pop[j].phenotyp.gread(k).real;
		dist += Norm(diff, size);
		count++;
	}
	dist /= count;
	****************************************************************/
	
	// Output PSO statistics
	if(outputEveryNgens > 0 && 
	  (generations % outputEveryNgens == 0||generations==1)) {
		//pr(logFile, "%d %8d %10.2f %10.2f %10.2f %6.2f %8.2f\n", generations, evaluate.evals(), Pg.value(Normal_Eval), Pop_avg, Pi_avg, w, v_avg);
		//ORIG pr(logFile, "%8d %10ld %10.2f  \n", generations, evaluate.evals(), Pg.value(Normal_Eval));		
		fflush(logFile);
	}
	//pr(logFile, "%8d\t%6d\t%6.3f\t%6.3f\t%6.3f\n", generations, evaluate.evals(), X_best_value, Pg.value(Normal_Eval), dist);	
	//pr(logFile, "%6d\t%6d\t%6.3f\t%6.3f\n\n", generations, evaluate.evals(), X_best_value, Pg.value(Normal_Eval));		
	//fflush(logFile);
	
	return (0);
}

int ParticleSwarmGS::localsearch(Population &Pop, Local_Search *local_method)
{
// Set shorthand names for best for each indiv, best individual in pop
Population &Pi = (Population &)(*_Pi);
Individual &Pg = (Individual &)(*_Pg);
	if(local_method == NULL) return(-1);
	if(best<0||best>=Pop.num_individuals()) {
		stop("PSO LocalSearch bad best value");
	}
#ifdef PSODEBUG
 // disabled by MP Jun 2011 after code move TODO
	// MP DEBUG 2011
	(void)fprintf( logFile, "PSO LS only Pop[best=%d] was %f\n", 
	 best, Pop[0].value(localEvalMode));
#endif
	// If local seach method is defined, apply local search			
	local_method->search(Pop[best]);

	if(Pop[best].value(Normal_Eval) < Pi[best].value(Normal_Eval)) {
		Pi[best] = Pop[best];			
		Pg = Pi[best]; // mp maybe not needed
	}					
#ifdef PSODEBUG
 // disabled by MP Jun 2011 after code move TODO
	// MP DEBUG 2011:
	(void)fprintf( logFile, "PSO LS only Pop[best=%d] now %f\n", 
	 best, Pop[0].value(localEvalMode));
   if (outlev > 1) {
	(void)fprintf( logFile, " %f",Pop[best].value(localEvalMode)); 
	(void)fprintf( logFile, " \n"); 
	    }
#endif
	return(0);
}
