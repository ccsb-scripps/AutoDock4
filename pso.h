#ifndef _pso_h
#define _pso_h

#include "gs.h"
#include "ls.h"

extern FILE *logFile;

class ParticleSwarmGS : public Global_Search 
{
	private:
		Population *_Pi;	// best solution for each individual in its own searching history
		Individual	*_Pg;	// current best solution
		int pop_size;	// population size
		int size;	// number of dimensions (7*nlig + num_torsion)
		float **v;	// velocity
		float *vmax;	//max velocity 
		float *vmin;	// min velocity
		float w;	   // inertia weight
		float wmax;
		float wmin;
		float c1;	// cognitive
		float c2;	// social
	    int K;      // number of neighbor particles
	    float c;    // constriction factor for cPSO
	    
		int generations;
		int num_evals;
		int outputEveryNgens;
	
		Local_Search *LocalSearchMethod;
	
	public:	
		~ParticleSwarmGS();
		ParticleSwarmGS(
			int init_size, 
			float *init_vmax, 
			float *init_vmin, 
			float init_wmax, 
			float init_wmin,
			Local_Search *init_LS, 
			int init_num_evals,
			float pso_c1,
			float pso_c2,
			int pso_k, 
			int init_ouput_gens); 			
		
		Individual& getBest();	
		
		// The following part are derived virtual functions
		void reset(void);
        void reset(unsigned int);
        int terminate(void);
        int search(Population &);  			
};

			
inline ParticleSwarmGS::~ParticleSwarmGS()
{
	if(_Pi)
		delete _Pi;
	if(_Pg)
		delete _Pg;
	if(v) {
		for(int i=0;i < pop_size; i++)	
			delete [] v[i];
		delete [] v;
	}
}

inline ParticleSwarmGS::ParticleSwarmGS(
			int init_size, 
			float *init_vmax, 
			float *init_vmin, 
			float init_wmax, 
			float init_wmin,
			Local_Search *init_LS, 
			int init_num_evals, 
			float pso_c1,
			float pso_c2,
			int pso_k,
			int init_ouput_gens)			
{
	size = init_size;
	vmax = init_vmax; 
	vmin = init_vmin;			  
	wmax = init_wmax; 
	wmin = init_wmin;
	LocalSearchMethod = init_LS;
	num_evals = init_num_evals; 
	outputEveryNgens = init_ouput_gens;
	_Pi =NULL; 
	_Pg = NULL ;
	 v = NULL; 
	c1 = pso_c1;	      
	c2 = pso_c1; 
	K = pso_k;
    c = 0.01;
}

inline Individual& ParticleSwarmGS:: getBest()
{
	return *_Pg;
}

// The following part are derived virtual functions
//int search(Population &);
inline int ParticleSwarmGS::terminate(void)
{
	return 0;
}
	
inline void ParticleSwarmGS::reset(void)
{
	generations = 0;
	w = wmax;
	if(_Pi)
		delete _Pi;
	if(_Pg)
		delete _Pg;
	if(v) {
		for(int i=0;i < pop_size; i++)	
			delete [] v[i];
		delete [] v;
	}
	
	_Pi = NULL;
	_Pg = NULL;
	v = NULL;
}
	
inline void ParticleSwarmGS::reset(unsigned int extOutputEveryNgens)
{
	reset();
}
	


#endif
