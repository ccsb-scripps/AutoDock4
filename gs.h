//  These are the classes in the global search hierarchy.  Notice that
//  Global_Search is an abstract base class and as such it can never
//  be instantiated.  For now, the only global search operator is the
//  Genetic Algorithm.
//  rsh 07/08/95

#ifndef _GLOBAL_SEARCH_H
#define _GLOBAL_SEARCH_H

#include "support.h"

enum M_mode {ERR = -1, BitFlip, CauchyDev, IUniformSub};
enum Selection_Mode {Proportional=0, Tournament=1, Boltzmann=2};
enum Xover_Mode {TwoPt=0, OnePt=1, Uniform=2};
enum Worst_Mode {AverageOfN, OfN, Ever};

class Global_Search
{
   public:
      Global_Search(void);
      virtual ~Global_Search(void);
      virtual int search(Population &) = 0;
      virtual int terminate(void) = 0;
      virtual void reset(void) = 0;
      virtual void reset(unsigned int) = 0;
};

// The class Genetic_Algorithm is a Global_Search method,
// 
class Genetic_Algorithm : public Global_Search
{
//   friend void debug(Genetic_Algorithm &, Population &);
   private:
      unsigned int generations, max_generations;
      unsigned int converged; // gmm 7-jan-98
      unsigned int *ordering;
      unsigned int outputEveryNgens; // gmm 2000.11.1
      float *alloc, *mutation_table;
      EvalMode e_mode;
      Selection_Mode s_mode;
      Worst_Mode w_mode;
      int window_size, elitism;
      Xover_Mode c_mode;
      unsigned int m_table_size;
      float m_rate, c_rate, tournament_prob;
      float alpha, beta;
      float tranStep, quatStep, torsStep;
      int low, high;
      double avg, worst;
      double *worst_window;

      double worst_this_generation(Population &);
      void set_worst(Population &);
      void make_table(int, float);
      int check_table(float);
      M_mode m_type(RepType);
      void mutate(Genotype &, int);
      void mutation(Population &);
      void crossover(Population &);
      void crossover_2pt(Genotype &, Genotype &, unsigned int, unsigned int);
      void selection_proportional(Population &, Individual *);
      void selection_tournament(Population &, Individual *);
      Individual *selection(Population &);

   public:
      Genetic_Algorithm(void);
      // Genetic_Algorithm(EvalMode, Selection_Mode, Xover_Mode, Worst_Mode, int, float, float, int, unsigned int); // before 2000.11.1
      Genetic_Algorithm(EvalMode, Selection_Mode, Xover_Mode, Worst_Mode, int, float, float, int, unsigned int, unsigned int); // after 2000.11.1
      ~Genetic_Algorithm(void);
      void initialize(unsigned int, unsigned int);
      void mutation_values(int, int, float, float);
      unsigned int num_generations(void);
      void reset(void);
      void reset(unsigned int);
      int terminate(void);
      int search(Population &);
};

//  Inline Functions
inline Global_Search::Global_Search(void)
{
}

inline Global_Search::~Global_Search(void)
{
}

// Default values set in this constructor.
inline Genetic_Algorithm::Genetic_Algorithm(void)
: m_table_size(0), mutation_table(NULL), ordering(NULL), alloc(NULL), worst_window(NULL)
{
   generations = 0;
   elitism = window_size = low = high = 0;
   m_rate = 0.02;
   c_rate = 0.80;
   alpha = beta = 0.0;
   tranStep = 2.0;
   quatStep = torsStep = Rad( 50.0 );
   worst = avg = 0.0L;
   converged = 0; // gmm 7-jan-98
   outputEveryNgens = 100; // gmm 2000-nov-1
}

inline Genetic_Algorithm::~Genetic_Algorithm(void)
{
   if (worst_window!=NULL) {
      delete [] worst_window;
   }

   if (alloc!=NULL) {
      delete [] alloc;
   }

   if (ordering!=NULL) {
      delete [] ordering;
   }

   if (mutation_table!=NULL) {
      delete [] mutation_table;
   }
}

inline void Genetic_Algorithm::mutation_values(int init_low, int init_high, float init_alpha, float init_beta)
{
   low = init_low;
   high = init_high;
   alpha = init_alpha;
   beta = init_beta;
}

inline unsigned int Genetic_Algorithm::num_generations(void)
{
   return(generations);
}

inline int Genetic_Algorithm::terminate(void)
{
   if (max_generations>0) {
      // before 7-jan-98, was: return(generations>=max_generations);
      return((generations>=max_generations)||(converged==1)); // gmm 7-jan-98
   } else {
      return(0);  //  Don't terminate
   }
}

inline void Genetic_Algorithm::reset(void)
{
   generations = 0;
   converged = 0; // gmm 7-jan-98
}

inline void Genetic_Algorithm::reset(unsigned int extOutputEveryNgens) // gmm 2000.11.1
{
   outputEveryNgens = extOutputEveryNgens; // gmm 2000.11.1
   generations = 0;
   converged = 0; // gmm 7-jan-98
}

#endif
