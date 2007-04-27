/*

 $Id: gs.h,v 1.11 2007/04/27 06:01:49 garrett Exp $

 AutoDock 

 Copyright (C) 1989-2007,  Garrett M. Morris, David S. Goodsell, Ruth Huey, Arthur J. Olson, 
 All Rights Reserved.

 AutoDock is a Trade Mark of The Scripps Research Institute.

 This program is free software; you can redistribute it and/or
 modify it under the terms of the GNU General Public License
 as published by the Free Software Foundation; either version 2
 of the License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

 */

//  These are the classes in the global search hierarchy.  Notice that
//  Global_Search is an abstract base class and as such it can never
//  be instantiated.  For now, the only global search operator is the
//  Genetic Algorithm.
//  rsh 07/08/95

#ifndef _GLOBAL_SEARCH_H
#define _GLOBAL_SEARCH_H

#include "support.h"

enum M_mode { ERR = -1, BitFlip, CauchyDev, IUniformSub };
enum Selection_Mode { Proportional=0, Tournament=1, Boltzmann=2 };
enum Xover_Mode { TwoPt=0, OnePt=1, Uniform=2, Arithmetic=3 };
enum Worst_Mode { AverageOfN, OfN, Ever };

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
      EvalMode e_mode;
      Selection_Mode s_mode;
      Xover_Mode c_mode;
      Worst_Mode w_mode;
      unsigned int elitism;
	  Real c_rate;
	  Real m_rate;
      unsigned int window_size;
      Real alpha;
	  Real beta;
      Real tranStep, quatStep, torsStep;
      int low, high; // should these be int or Real?
      unsigned int generations; 
	  unsigned int max_generations;
      unsigned int outputEveryNgens; // gmm 2000.11.1,2003.08.18
      unsigned int converged; // gmm 7-jan-98
 	  Real *alloc;
      Real *mutation_table;
      unsigned int *ordering;	  
      unsigned int m_table_size;
      double worst, avg;
      double *worst_window;
	  Real tournament_prob;

      double worst_this_generation(Population &);
      void set_worst(Population &);
      void make_table(int, Real);
      int check_table(Real);
      M_mode m_type(RepType);
      void mutate(Genotype &, int);
      void mutation(Population &);
      void crossover(Population &);
      void crossover_2pt(Genotype &, Genotype &, unsigned int, unsigned int);
      void crossover_uniform(Genotype &, Genotype &, unsigned int);
      void crossover_arithmetic(Genotype &, Genotype &, Real);
      void selection_proportional(Population &, Individual *);
      void selection_tournament(Population &, Individual *);
      Individual *selection(Population &);

   public:
      Genetic_Algorithm(void);
      // Genetic_Algorithm(EvalMode, Selection_Mode, Xover_Mode, Worst_Mode, int, Real, Real, int, unsigned int); // before 2000.11.1
      Genetic_Algorithm(EvalMode, Selection_Mode, Xover_Mode, Worst_Mode, int, Real, Real, int, unsigned int, unsigned int); // after 2000.11.1
      ~Genetic_Algorithm(void);
      void initialize(unsigned int, unsigned int);
      void mutation_values(int, int, Real, Real,  Real, Real, Real );
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
: alloc(NULL), mutation_table(NULL), ordering(NULL), m_table_size(0), worst_window(NULL)
{
   generations = 0;
   elitism = window_size = low = high = 0;
   m_rate = 0.02;
   c_rate = 0.80;
   alpha = beta = 0.0;
   tranStep = 2.0;
   quatStep = torsStep = DegreesToRadians( 30.0 );
   worst = avg = 0.0L;
   converged = 0; // gmm 7-jan-98
   outputEveryNgens = OUTLEV1_GENS; // gmm 2000-nov-1
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

inline void Genetic_Algorithm::mutation_values(int init_low, int init_high, 
        Real init_alpha, Real init_beta, 
        Real init_tranStep, Real init_quatStep, Real init_torStep )
{
   low = init_low;
   high = init_high;
   alpha = init_alpha;
   beta = init_beta;
   tranStep = init_tranStep;
   quatStep = init_quatStep;
   torsStep = init_torStep;
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
