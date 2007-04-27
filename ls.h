/*

 $Id: ls.h,v 1.5 2007/04/27 06:01:49 garrett Exp $

 AutoDock 

 Copyright (C) 1989-2007,  Scott Halliday, Rik Belew, Garrett M. Morris, David S. Goodsell, Ruth Huey, Arthur J. Olson, 
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

//  These are the classes associated with the local search operator hierarchy.
//  Notice that Local_Search is an abstract base class and as such cannot be
//  instantiated.  For now, local search is just embodied by Solis & Wets 
//  Algorithms 1 & 2.
//  rsh 07/08/95
//  At the suggestion of Bruce Duncan, the Pseudo_Solis_Wets classes were
//  added.  The major difference between these classes and the regular 
//  Solis_Wets classes is that the variance varies across the dimensions
//  rsh 02/16/96

#ifndef _LOCAL_SEARCH_H
#define _LOCAL_SEARCH_H

#include "support.h"
#include "ranlib.h"
#include <stdlib.h>

class Local_Search
{
   public:
      Local_Search(void);
      virtual ~Local_Search(void);
      virtual void reset(void) = 0;
      virtual int terminate(void) = 0;
      virtual int search(Individual &) = 0;
};

class Pattern_Search : public Local_Search
{
   protected:
      unsigned int size; 
			unsigned int max_success;
      Real step_size, current_step_size;
      Real step_threshold, expansion, contraction;
			Real search_frequency;
			Real *pattern;
			unsigned int *index;
			unsigned int successes;
      //int HJ_bias; // not implemented
      //Real * step_scales; // not implemented
			Phenotype exploratory_move(const Phenotype&);
			Phenotype pattern_explore(const Phenotype&);
			Phenotype pattern_move(const Phenotype&);
			void reset_pattern(void);
			void reset_indexes(void);
			void shuffle_indexes(void);
   public:
      Pattern_Search(void);
      Pattern_Search(unsigned int, unsigned int, Real, Real, Real, Real, Real);
      ~Pattern_Search(void);
      void reset(void);
      int terminate(void);
      int search(Individual &);
};

class Solis_Wets_Base : public Local_Search
{
   protected:
      unsigned int size, max_its, max_successes, max_failures;
      Real expansion, contraction;
      Real search_frequency;
      Real *deviates, *bias;

   public:
      Solis_Wets_Base(void);
      Solis_Wets_Base(unsigned int, unsigned int, unsigned int, unsigned int, Real, Real, Real);
      virtual ~Solis_Wets_Base(void);
      virtual double gen_deviates(Real) = 0;
      virtual void SW(Phenotype &) = 0;
      virtual void reset(void);
      virtual int terminate(void);
      int search(Individual &);
};

class Solis_Wets : public Solis_Wets_Base
{
   protected:
      Real rho, lower_bound_on_rho;

   public:
      Solis_Wets(void);
      Solis_Wets(unsigned int, unsigned int, unsigned int, unsigned int, Real, Real, Real, Real, Real);
      virtual ~Solis_Wets(void);
      virtual double gen_deviates(Real) = 0;
      void SW(Phenotype &);
};

class Pseudo_Solis_Wets : public Solis_Wets_Base
{
   protected:
      Real *rho, *lower_bound_on_rho;
      Real *temp_rho;

   public:
      Pseudo_Solis_Wets(void);
      Pseudo_Solis_Wets(unsigned int, unsigned int, unsigned int, unsigned int, Real, Real, Real);
      Pseudo_Solis_Wets(unsigned int, unsigned int, unsigned int, unsigned int, Real, Real, Real, Real *, Real *);
      virtual ~Pseudo_Solis_Wets(void);
      virtual double gen_deviates(Real) = 0;
      void SW(Phenotype &);
};

class Solis_Wets1 : public Solis_Wets
{
   public:
      Solis_Wets1(void);
      Solis_Wets1(unsigned int, unsigned int, unsigned int, unsigned int, Real, Real, Real, Real, Real);
      ~Solis_Wets1(void);
      double gen_deviates(Real);
};

class Solis_Wets2 : public Solis_Wets
{
   public:
      Solis_Wets2(void);
      Solis_Wets2(unsigned int, unsigned int, unsigned int, unsigned int, Real, Real, Real, Real, Real);
      ~Solis_Wets2(void);
      double gen_deviates(Real);
};

class Pseudo_Solis_Wets1 : public Pseudo_Solis_Wets
{
   public:
      Pseudo_Solis_Wets1(void);
      Pseudo_Solis_Wets1(unsigned int, unsigned int, unsigned int, unsigned int, Real, Real, Real);
      Pseudo_Solis_Wets1(unsigned int, unsigned int, unsigned int, unsigned int, Real, Real, Real, Real *, Real *);
      ~Pseudo_Solis_Wets1(void);
      double gen_deviates(Real);
};

class Pseudo_Solis_Wets2 : public Pseudo_Solis_Wets
{
   public:
      Pseudo_Solis_Wets2(void);
      Pseudo_Solis_Wets2(unsigned int, unsigned int, unsigned int, unsigned int, Real, Real, Real);
      Pseudo_Solis_Wets2(unsigned int, unsigned int, unsigned int, unsigned int, Real, Real, Real, Real *, Real *);
      ~Pseudo_Solis_Wets2(void);
      double gen_deviates(Real);
};

//  Inline Functions
inline Local_Search::Local_Search(void)
{
}

inline Local_Search::~Local_Search(void)
{
}

inline Solis_Wets_Base::Solis_Wets_Base(void)
:  size(0), max_its(10), max_successes(5), max_failures(5), expansion(2.0), contraction(0.5),
   search_frequency(1.0), deviates(NULL), bias(NULL)
{
}

inline Solis_Wets_Base::Solis_Wets_Base(unsigned int init_size, unsigned int init_max_its, 
                                        unsigned int init_max_succ, unsigned int init_max_fail, 
                                        Real init_expansion, Real init_contraction, Real init_search_freq)
:  size(init_size), max_its(init_max_its), max_successes(init_max_succ), max_failures(init_max_fail),
   expansion(init_expansion), contraction(init_contraction), search_frequency(init_search_freq)
{
   bias = new Real[size];
   deviates = new Real[size];
}

inline Solis_Wets_Base::~Solis_Wets_Base(void)
{
   if(deviates!=NULL)
   {
      delete [] deviates;
   }

   if(bias!=NULL)
   {
      delete [] bias;
   }
}

inline void Solis_Wets_Base::reset(void)
{
   //  Do nothing
}

inline int Solis_Wets_Base::terminate(void)
{
   return(0);  //  Don't terminate
}

inline Solis_Wets::Solis_Wets(void)
:  Solis_Wets_Base(), rho(1.0), lower_bound_on_rho(0.0)
{
}

inline Solis_Wets::Solis_Wets(unsigned int init_size, unsigned int init_max_its, unsigned int init_max_succ, 
                              unsigned int init_max_fail, Real init_rho, Real init_lb_on_rho, 
                              Real init_expansion, Real init_contraction, Real init_search_freq)
:  Solis_Wets_Base(init_size, init_max_its, init_max_succ, init_max_fail, init_expansion, init_contraction, 
                   init_search_freq), rho(init_rho), lower_bound_on_rho(init_lb_on_rho)
{
}

inline Solis_Wets::~Solis_Wets(void)
{
}

inline Pseudo_Solis_Wets::Pseudo_Solis_Wets(void)
:  Solis_Wets_Base(), rho(NULL), lower_bound_on_rho(NULL), temp_rho(NULL)
{
}

inline Pseudo_Solis_Wets::Pseudo_Solis_Wets(unsigned int init_size, unsigned init_max_its, 
                                            unsigned int init_max_succ, unsigned int init_max_fail, 
                                            Real init_expansion, Real init_contraction, Real init_search_freq)
:  Solis_Wets_Base(init_size, init_max_its, init_max_succ, init_max_fail, init_expansion, init_contraction, 
                   init_search_freq), rho(NULL), lower_bound_on_rho(NULL), temp_rho(NULL)
{
}

inline Pseudo_Solis_Wets::Pseudo_Solis_Wets(unsigned int init_size, unsigned init_max_its, 
                                            unsigned int init_max_succ, unsigned int init_max_fail, 
                                            Real init_expansion, Real init_contraction, Real init_search_freq, 
                                            Real *init_rho, Real *init_lb_on_rho)
:  Solis_Wets_Base(init_size, init_max_its, init_max_succ, init_max_fail, init_expansion, init_contraction, 
                   init_search_freq), rho(init_rho), lower_bound_on_rho(init_lb_on_rho)
{
   temp_rho = new Real[init_size];
}

inline Pseudo_Solis_Wets::~Pseudo_Solis_Wets(void)
{
   if (rho!=NULL)
   {
      delete [] rho;
   }

   if (lower_bound_on_rho!=NULL)
   {
      delete [] lower_bound_on_rho;
   }

   if (temp_rho!=NULL)
   {
      delete [] temp_rho;
   }
}

inline Solis_Wets1::Solis_Wets1(void)
: Solis_Wets()
{
}

inline Solis_Wets1::Solis_Wets1(unsigned int init_size, unsigned int init_max_its, unsigned int init_max_succ, 
                                unsigned int init_max_fail, Real init_rho, Real init_lb_on_rho, 
                                Real init_expansion, Real init_contraction, Real init_search_freq)
:  Solis_Wets(init_size, init_max_its, init_max_succ, init_max_fail, init_rho, init_lb_on_rho, init_expansion, 
              init_contraction, init_search_freq)
{
}

inline Solis_Wets1::~Solis_Wets1(void)
{
}

inline double Solis_Wets1::gen_deviates(Real rho)
{
   return(gennor(0.0, rho));
}

inline Solis_Wets2::Solis_Wets2(void)
:  Solis_Wets()
{
}

inline Solis_Wets2::Solis_Wets2(unsigned int init_size, unsigned int init_max_its, unsigned int init_max_succ, 
                              unsigned int init_max_fail, Real init_rho, Real init_lb_on_rho, 
                              Real init_expansion, Real init_contraction, Real init_search_freq)
:  Solis_Wets(init_size, init_max_its, init_max_succ, init_max_fail, init_rho, init_lb_on_rho, init_expansion,
              init_contraction, init_search_freq)
{
}

inline Solis_Wets2::~Solis_Wets2(void)
{
}

inline double Solis_Wets2::gen_deviates(Real rho)
{
   return(genunf(-rho/2.0, rho/2.0));
}

inline Pseudo_Solis_Wets1::Pseudo_Solis_Wets1(void)
:  Pseudo_Solis_Wets()
{
}

inline Pseudo_Solis_Wets1::Pseudo_Solis_Wets1(unsigned int init_size, unsigned int init_max_its, 
                                              unsigned int init_max_succ, unsigned int init_max_fail,  
                                              Real init_expansion, Real init_contraction, 
                                              Real init_search_freq)
:  Pseudo_Solis_Wets(init_size, init_max_its, init_max_succ, init_max_fail, init_expansion,
                     init_contraction, init_search_freq)
{
}

inline Pseudo_Solis_Wets1::Pseudo_Solis_Wets1(unsigned int init_size, unsigned int init_max_its, 
                                              unsigned int init_max_succ, unsigned int init_max_fail,  
                                              Real init_expansion, Real init_contraction, 
                                              Real init_search_freq, Real *init_rho,
                                              Real *init_lb_on_rho)
:  Pseudo_Solis_Wets(init_size, init_max_its, init_max_succ, init_max_fail, init_expansion,
                     init_contraction, init_search_freq, init_rho, init_lb_on_rho)
{
}

inline Pseudo_Solis_Wets1::~Pseudo_Solis_Wets1(void)
{
}

inline double Pseudo_Solis_Wets1::gen_deviates(Real rho)
{
   return(gennor(0.0, rho));
}

inline Pseudo_Solis_Wets2::Pseudo_Solis_Wets2(void)
:  Pseudo_Solis_Wets()
{
}

inline Pseudo_Solis_Wets2::Pseudo_Solis_Wets2(unsigned int init_size, unsigned int init_max_its, 
                                              unsigned int init_max_succ, unsigned int init_max_fail,  
                                              Real init_expansion, Real init_contraction, 
                                              Real init_search_freq)
:  Pseudo_Solis_Wets(init_size, init_max_its, init_max_succ, init_max_fail, init_expansion,
                     init_contraction, init_search_freq)
{
}

inline Pseudo_Solis_Wets2::Pseudo_Solis_Wets2(unsigned int init_size, unsigned int init_max_its, 
                                              unsigned int init_max_succ, unsigned int init_max_fail,  
                                              Real init_expansion, Real init_contraction, 
                                              Real init_search_freq, Real *init_rho,
                                              Real *init_lb_on_rho)
:  Pseudo_Solis_Wets(init_size, init_max_its, init_max_succ, init_max_fail, init_expansion,
                     init_contraction, init_search_freq, init_rho, init_lb_on_rho)
{
}

inline Pseudo_Solis_Wets2::~Pseudo_Solis_Wets2(void)
{
}

inline double Pseudo_Solis_Wets2::gen_deviates(Real rho)
{
   return(genunf(-rho/2.0, rho/2.0));
}

#endif
