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
 ***********************************************************************/
int ParticleSwarmGS::search(Population &Pop)
{
	int i, j, g;
	int best;
	double r1, r2;
	double curVal;
	double newVal;
	double piCurE = 0.0;
	//float phi, c;
	float ratio = 0.0;
	
	// initialize velocity		
	if(_Pi == NULL) {
		
		pr(logFile, "PSO max_generations = %d\n", max_generations);
		pr(logFile, "PSO size = %d\n", size); // MP debug
				
		fprintf(logFile, "PSO Allocate initial velocity of particles...\n");
		_Pi = new Population(Pop);
		
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
	
	Population &Pi = (Population &)(*_Pi);
	Individual &Pg = (Individual &)(*_Pg);
	
	//double piBestE[pop_size];	    // E array for personal Best	
	//double prevE[pop_size];       // E array for particles in previous run 
	//double curE[pop_size];        // E array for particles in current run
		   		   	 	     	   	 	   	   	   	  
	//update position variable related to the translation and rotation of each particle			
	for(i = 0; i < pop_size; i++) {
				
		//r2 = 1.0 - r1;
		//get the best in neighborhood, e.g., Neighborhood Best (nbBest)		
		g = i;
		for(j = i + 1; j < i + K; j++) {
			if(Pi[j % pop_size].value(Normal_Eval) < Pi[g].value(Normal_Eval))
				g = (j % pop_size);
		}
				
		// size is the dimension of docking search (e.g, n*nlig + ntor)			
		for(j = 0; j < size; j++) {		
			r1 = ranf();
			r2 = ranf();
			// update velocity
			curVal = Pop[i].phenotyp.gread(j).real;
									
			
			// -Huameng
			if(ratio <= 0.8) {				
				v[i][j] = w * v[i][j] + c1 * r1 * (Pi[i].phenotyp.gread(j).real - curVal)
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
										
			// update x,y,z, quarternion, torsion of particle i
			Pop[i].phenotyp.write(newVal, j);												
		}
		// need to normalize Quarternion
		Quat q;
		q = Pop[i].phenotyp.readQuat();
		Pop[i].phenotyp.writeQuat( normQuat( q ));
		
		//evaluate new solution
		//this will be done in step: 'Find the best in Pop'						
		//Pop[i].value(Normal_Eval);
		
	}	// end current swarm	
	
	// Find the best in Pop	
	best = 0;
	double X_best_value = 99999999; // MP debug Pop[0].value(Normal_Eval);	
	for(i = 0; i < pop_size; i++) {		
		Pop[i].inverse_mapping(); // MP@@ copy phenotype to genotype
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
	// MP debug - verifying msort is bringing best individual to Pop[0]
	pr(logFile, "PSO pre-msort  best=%d Pop[%d]=%.2f X_best_value=%.2f, Pop[0,1,..] = ", 
	  best, best, Pop[best].value(Normal_Eval), X_best_value);
	for(i=0;i<8; i++) pr(logFile, "%8.2f ", Pop[i].value(Normal_Eval));
	pr(logFile, "\n");
	Pop.msort(1); // MP@@ - bring best individual to position 0 (for local search)
	// MP debug - look again for best after msort, should be at [0]
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
	pr(logFile, "PSO post-msort best=%d Pop[%d]=%.2f X_best_value=%.2f, Pop[0,1,..] = ", 
	  best, best, Pop[best].value(Normal_Eval), X_best_value);
	for(i=0;i<8; i++) pr(logFile, "%8.2f ", Pop[i].value(Normal_Eval));
	pr(logFile, "\n");

	////////////////////////////////////////////////////////
	// Local Search
	////////////////////////////////////////////////////////
	// If local seach method is defined, apply local search			
//	if(LocalSearchMethod) {							
//		LocalSearchMethod->search(Pop[best]);
					
//		if( Pop[best].value(Normal_Eval) < Pi[best].value(Normal_Eval) ) {
//			Pi[best] = Pop[best];			
//		}					
		//pr(logFile, "NumEvals=%d\tX_best before LS=%.2f\tX_best after LS=%.2f\n",	
		//			evaluate.evals(), X_best_value, Pop[best].value(Normal_Eval));		
//	}		
	
	
	// Update Pi, personal Best in history
	//for(i = 0; i < pop_size; i++) {					
	//	if(Pop[i].value(Normal_Eval) < Pi[i].value(Normal_Eval))	
	//		Pi[i] = Pop[i];			
	//}
				
	// update Global Best Pg, which points to _Pg
	best = 0;	
	for(i = 1; i < pop_size; i++ ) {		
		if( Pi[i].value(Normal_Eval) < Pi[best].value(Normal_Eval) )		
			best = i;			
	}
	
	Pg = Pi[best];
	
	// Update weight
	// Huameng:  ratio = (float)evaluate.evals() / num_evals;

	ratio = (float)generations / max_generations;  // MP since we dont have num_evals
	w = wmin +(wmax - wmin) * (1.0 - ratio);
			 		
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
	
	// Output Information
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
