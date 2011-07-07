#ifndef _pso_h
#define _pso_h

#include "gs.h"
#include "ls.h" 
#include "structs.h"

extern FILE *logFile;

// define class for PSO algorithmic options
// to avoid having to change dozens of method signatures whenever
// the option set is changed. MP TSRI 2011
struct PSO_Options {
	float pso_w;	   // inertia weight
	float wmax;	// pso_w at beginning of run
	float wmin;	// pso_w at conclusion of run, see pso.cc
	float c1;	// cognitive
	float c2;	// social
	int pso_K;      // number of neighbor particles
	float c;    // constriction factor for cPSO
	Boole pso_neighbors_dynamic; // MP
	Boole pso_random_by_dimension; // MP
	Boole pso_interpolate_as_scalars; // MP
  public:
    PSO_Options () :
        pso_w(1.), // w
        wmax(1.), // wmax
        wmin(1.), // wmin
        c1(6.), // c1
        c2(6.), // c2
        pso_K(4),  // pso_K
        c(0.01),  // c
        pso_neighbors_dynamic(false), // 
        pso_random_by_dimension(true),  // 
        pso_interpolate_as_scalars(true)  //
        { }
	};

class ParticleSwarmGS : public Global_Search 
{
	private:
		Population *_Pi;	// best solution for each individual in its own searching history
		Individual	*_Pg;	// current best solution
		int best; // index in Pi of current global best solution
		int pop_size;	// population size
		int size;	// number of dimensions (7*nlig + num_torsion)
		float **v;	// velocity
		float *vmax;	//max velocity 
		float *vmin;	// min velocity
		PSO_Options pso_options;
	    
		int generations;
		int outputEveryNgens;
        Output_pop_stats output_pop_stats;
	
		Local_Search *LocalSearchMethod;
	
	public:	
		~ParticleSwarmGS();
		ParticleSwarmGS(
			float *init_vmax, 
			float *init_vmin, 
			PSO_Options init_pso_options, 
			Local_Search *init_LS, 
			unsigned int init_max_evals, 
			unsigned int init_max_generations, 
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
	int localsearch(Population &, Local_Search *);
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
			PSO_Options init_pso_options, 
			Local_Search *init_LS, 
			const unsigned int init_max_evals, 
			const unsigned int init_max_generations, 
			Output_pop_stats init_output_pop_stats) :
    Global_Search(  init_max_evals, init_max_generations)
{
	vmax = init_vmax; 
	vmin = init_vmin;			  
	pso_options = init_pso_options;
	LocalSearchMethod = init_LS;
    generations = 0;
	output_pop_stats = init_output_pop_stats;
	_Pi =NULL; 
	_Pg = NULL ;
	 v = NULL; 
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
   if (max_generations>(unsigned) 0) {
      return((unsigned)generations>=max_generations); 
   } else {
      return(0);  //  Don't terminate
   }
}

	
inline void ParticleSwarmGS::reset(void)
{
	generations = 0;
	//MP pso_w = wmax;
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
