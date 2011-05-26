#ifndef _pso_h
#define _pso_h

#include "gs.h"
#include "ls.h"
#include "structs.h"

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
		int max_generations;
		int outputEveryNgens;
        Output_pop_stats output_pop_stats;
	
		Local_Search *LocalSearchMethod;
	
	public:	
		~ParticleSwarmGS();
		ParticleSwarmGS(
			float *init_vmax, 
			float *init_vmin, 
			float init_wmax, 
			float init_wmin,
			Local_Search *init_LS, 
			float pso_c1,
			float pso_c2,
			int pso_k, 
			int max_generations,
			Output_pop_stats output_pop_stats); 			
		
		Individual& getBest();	

      void initialize(const unsigned int, const unsigned int);
      unsigned int num_generations(void) const;
		
		// The following part are derived virtual functions
        char *shortname(void);
        char *longname(void);
		void reset(void);
        void reset(const Output_pop_stats&);
        int terminate(void);
        int search(Population &);  			
};

inline char * ParticleSwarmGS::shortname(void)
{
        return "PSO";
}

inline char * ParticleSwarmGS::longname(void)
{
        return "PARTICLE SWARM OPTIMIZATION";
}

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
			float *init_vmax, 
			float *init_vmin, 
			float init_wmax, 
			float init_wmin,
			Local_Search *init_LS, 
			float pso_c1,
			float pso_c2,
			int pso_k,
			const int init_max_generations, 
			Output_pop_stats init_output_pop_stats)
{
	vmax = init_vmax; 
	vmin = init_vmin;			  
	wmax = init_wmax; 
	wmin = init_wmin;
	LocalSearchMethod = init_LS;
    generations = 0;
    max_generations = init_max_generations;
	output_pop_stats = init_output_pop_stats;
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
inline void ParticleSwarmGS::initialize(const unsigned int init_pop_size, const unsigned int ndims)
{

    pop_size = init_pop_size;
    size = ndims;
}


// The following part are derived virtual functions
//int search(Population &);

inline int ParticleSwarmGS::terminate(void)
{
   if (max_generations>0) {
      return(generations>=max_generations); 
   } else {
      return(0);  //  Don't terminate
   }
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
	
inline void ParticleSwarmGS::reset(const Output_pop_stats &init_output_pop_stats)
{
    output_pop_stats = init_output_pop_stats; 
	reset();
}
	


#endif
