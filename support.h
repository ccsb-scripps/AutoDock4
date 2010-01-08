/*

 $Id: support.h,v 1.16 2010/01/08 20:13:47 mp Exp $

 AutoDock 

Copyright (C) 2009 The Scripps Research Institute. All rights reserved.
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

//  These are the class used to support the Representation classes.
//  Genotypes are the representations that the Global_Search class 
//  and its derivations acts on.  Local_Search and its children act
//  on the Phenotype classes.  Phenotypes are basically what results
//  from mapping the Genotype to the solution domain.  It has the 
//  fundamental characteristic of fitness.  We need to make sure that 
//  the const pointers are never used to indirectly change values!
//  We can factor Genotype and Phenotype into 
//  a common base class Chromosome
//  rsh 07/08/95

/*
** $Log: support.h,v $
** Revision 1.16  2010/01/08 20:13:47  mp
** Extended population generations statistics, turned on if outlev>1.
**  Modified Files: RELEASENOTES call_glss.cc support.cc support.h
**
** Revision 1.15  2009/12/12 18:44:20  mp
** Added "Population at Generation:" statistics.
**  Modified Files: call_glss.cc configure.ac gs.cc support.cc support.h
**
** Revision 1.14  2009/05/08 23:02:18  rhuey
** Updated copyright notice in 188 source files
**
** Revision 1.13  2009/05/08 21:46:11  rhuey
** removed debugging comments and print-out
**
** Revision 1.12  2009/04/28 21:12:19  rhuey
** Changed so now Individual does mapping of its genotype into its phenotype and inverse_mapping of its phenotype into its genotype; in both cases returns a reference to itself; added a check for self-assignment
**
** Revision 1.11  2008/06/09 22:27:51  garrett
** Added "end_of_branch[MAX_TORS]" to the Population class, the logic being that every Individual in the Population is the same, so rather than add the overhead to all the Individuals, we added it to the Population.  Also introducted two new methods, set_eob() and get_eob(), to set the end_of_branch[] array, and get values given a key torsion number.  These changes are to support the new "Branch Crossover mode".
**
** Revision 1.10  2007/04/27 06:01:51  garrett
** Added the files necessary for GNU Autotools and the "dot-slash-configure dance"...
**
** Revision 1.9  2007/03/21 06:30:56  garrett
** Created a branch of AutoDock 4 with internal representation of orientations changed from axis-angle nx,ny,nz,ang to quaternion-components qx,qy,qz,qw.  This is intended to avoid rotation singularities of axis-angles near ((1,0,0),0 radians), and to avoid orientational bias in dockings.
**
** Revision 1.8  2006/11/03 02:10:48  garrett
** Significant change.  The initial population is now generated in a different way; previously, the axis that defined the rotation was created by generating uniformly-distributed random numbers in the range REALV_LOW to REALV_HIGH.  The same for the rotation angle (and torsion angles).  Now, we use a range of +/- 1 for the initial unit vector, and +/- PI for the rotation angle (and torsion angles).\
**
** Revision 1.7  2005/10/14 03:10:01  garrett
** Completed the "printPopulationAsCoordsEnergies" member function of the "Population" class, so that it now prints the nonbonded energy and the electrostatics energy, in addition to the translation and total energy for each member of the population.  These numbers are written to the population file at the end of each generation.  The DPF keyword "output_pop_file" expects the name of this population file; if this keyword is not given before a given "ga_run" command, then no population file will be written.
**
** Revision 1.6  2005/09/29 03:34:42  garrett
** Added a new method to the Population class, called "printPopulationAsCoordsEnergies", which is used to print out the translation of the centre of each individual and its total interaction energy.
**
** Revision 1.5  2004/12/07 02:07:53  gillet
** -- fix problem of compilation with g++ 3.3.2 on Linux:
** 	added Genotype(Genotype const &); in support.h
** 	it s definition in support.cc
**
** 	added Individual(Individual const &)
**
** Use the following message to resolve problem:
** http://gcc.gnu.org/ml/gcc-help/2003-10/msg00121.html
**
** You should put a copy constructor in your Anton class.
**
** "Anton(Anton& a)" isn't a copy constructor.  You need an "Anton(Anton const& a)".
**
** If Anton objects cannot be used in a copy constructor, then there are certain operations which Anton objects cannot perform.
**
** It appears that you hit upon one of them.
**
** If the "some code which _needs_ to modify a" is doing so in such a way that the LOGICAL state of the object is not affected, then those data members which are being modified should be marked as "mutable" in the class itself.  For example, certain reference counting schemes.  Another example is the std::string's c_str() method, which is a const method.
**
** If the "some code which _needs_ to modify a" does modify the LOGICAL state of the Anton object being copied from, then that's not kosher for use in a copy constructor.  C'est la vie.
**
** The Standard C++ Library auto_ptr<> is an example of a template class which modifies the state of the copied-from object.  That's one of the reasons that auto_ptr<>'s and STL don't mix (by-and-large).
**
** Or to say it another way, auto_ptr<> doesn't satisfy the contract requirements of STL containers.  Generally speaking.  If someone is REALLY careful, they may be able to use auto_ptr<> in a STL container... but I tend to recommend against it.   The BOOST <www.boost.org> folks have some smart pointer classes that are STL friendly.
**
** Revision 1.4  2004/11/16 23:42:56  garrett
** This is the result of merging the existing CVS respository with the AutoDock 4.0 code.  We have tested the code with a variety of problems: rigid ligand-rigid protein, rigid ligand-flexible protein, flexible ligand-rigid protein and flexible ligand-flexible protein: all four tests passed.  There was a bug fix regarding the flexible ligand-rigid protein case, to do with the absence of a BEGIN_RES record in the PDBQ file. -- GMM & RH
**
** Revision 1.3  2004/02/12 04:32:16  garrett
**
** After a first round of compilation using Apple's IDE, Xcode, various
** warnings have been eliminated (mainly unsigned ints and ints being
** interchanged); and
**
** After using Apple's Shark tool for profiling source code, the
** internal energy calculation has been optimized.
**
** The non-bonded cutoff used in the internal energy calculation has been
** reduced from 64 Angstroms to 8 Angstroms.  Most contributions beyond
** 8 Angstroms are very small, less than -0.001 kcal/mol, even for the
** largest atoms. Also, the conversion from double to int used to
** be done before the if to decide if we were within the cutoff; now
** the square of the distance is used in the comparison, and only if
** we are within the cutoff, do we convert from the double to int.
**
** The version checked in here still uses the type array to lookup
** the energy of interaction for a nonbond; this level of indirection
** can be pre-computed, and this should appear in my next round of checkins
**
** -- Garrett
**
** Revision 1.2  2002/10/30 01:49:15  garrett
** Commented out the #include <iostream.h> lines, since these appeared
** to conflict with <stdio.h>.  Also, added -lsupc++ to the linker
** options for Mac OS X 10.2, which now uses GCC 3.1; this may be
** necessary on GNU/Linux systems that use GCC 3.1.
**
** -- Lindy and Garrett
**
** Revision 1.1.1.1  2001/08/13 22:05:53  gillet
**  import initial of autodock sources
**
*/

#ifndef _SUPPORT_H
#define _SUPPORT_H

#include <stdio.h>
#include "rep.h"
#include "eval.h"
#include "structs.h"

enum EvalMode { Reset, Always_Eval, Normal_Eval, Always_Eval_Nonbond, Always_Eval_Elec };

typedef struct
{
   unsigned int vector;
   unsigned int index;
} Lookup;

//  For class Genotype, right now assume the user implements the
//  default constructor.
class Genotype
{
   //friend void debug(Genotype &);
   protected:
      //  Could some of these be made static?
      unsigned int number_of_genes;
      unsigned int number_of_vectors; // #vectors in rep_vector
      Lookup *lookup;		      // a table that helps in looking up a gene
      Representation **rep_vector; /* the actual representation of the genotype
				      like arrays of reals, bits, ints */
      unsigned modified : 1; /* used in caching for genotype operators, 
				e.g. crossover */

   public:
      Genotype(void);
      Genotype(Genotype &); /* copy constructor */
      Genotype(Genotype const &);
      Genotype(unsigned int, Representation **); /* creates a genotype from the
					     representation & total # vectors */
      ~Genotype(void); /* destructor */
      Genotype &operator=(const Genotype &);
      unsigned int num_vectors(void); /* e.g. "real,bit,bit,int" would = 4 */
      unsigned int num_genes(void); /* returns number_of_genes (see above) */
      RepType gtype(int); /* returns the type (real,bit,int) for 
							    a particular gene */
      const Element gread(int);
      const Representation *vread(int);
      void write(Element, int);
      void write(unsigned char, int);
      void write(FourByteLong, int);
      void write(double, int);
      void write(const Representation &, int);
      Quat readQuat();
      void writeQuat( Quat q );
};

//  Should Phenotype automatically evaluate itself upon construction?
class Phenotype
{
   //friend void debug(Phenotype &);
   protected:
      unsigned int number_of_dimensions, number_of_points;
      Lookup *lookup;
      Representation **value_vector;
      double value;
      unsigned evalflag : 1;  //  =1 means that this is the current evaluation
      unsigned modified : 1;  //  =1 means that this has been modified

   public:
      Phenotype(void);
      Phenotype(const Phenotype &);
      //Phenotype(const Genotype &);//to do
      Phenotype(unsigned int, Representation **);
      ~Phenotype(void);
      Phenotype &operator=(const Phenotype &);
      RepType gtype(int);
      const Element gread(int);
      const Representation *vread(int);
      void write(Element, int);
      void write(unsigned char, int);
      void write(FourByteLong, int);
      void write(double, int);
      void write(const Representation &, int);
      double evaluate(EvalMode);  //  This should return evaluation if that's the right answer, and it should evaluate otherwise.
      State make_state(int);
      unsigned int num_dimensions(void);
      unsigned int num_pts(void);
      Quat readQuat();
      void writeQuat( Quat q );
};

//  This should be an encapsulated class within Population
class Individual
{
   //friend void debug(Individual &);
   public:
      Genotype genotyp;   /* Genotype  is operated upon by *global search* operators */
      Phenotype phenotyp; /* Phenotype  "     "      "   " *local search*  operators, eg SW */
      Molecule *mol;		/* molecule */
      unsigned long age;	/* age of this individual; gmm, 1998-07-10 */

      Individual(void);
      Individual(Individual &); /* copy constructor */
      Individual(Individual const &);
      Individual(Genotype &, Phenotype &);
      ~Individual(void); /* destructor */
      Individual &operator=(const Individual &); /* assignment function for
						    individuals */
      Individual &mapping(void);         //updates phenotype from current genotype values 
      Individual &inverse_mapping(void); //updates genotype from current phenotype values 
      //Phenotype mapping(void); /* takes the genotype and converts it into a phenotype.  */
      //Genotype inverse_mapping(void);  // Scott should do: Also copy Phenotype's value
      double value(EvalMode); /* evaluation of the individual gives its value */
      State state(int); /* state variables in AutoDock */
      void  getMol(Molecule *); /* converts phenotype to mol's state and returns this individual's mol data */
      void printIndividualsState(FILE *, int, int); /* print out the state of this individual */
      void incrementAge(); /* make individual grow 1 generation older */
      int serial; // serial number of this individual
};

class Population
{
   //friend void debug(Population &);
   protected:
      int lhb;  //  These keep track of the lower & upper heap bounds
      int size; /* the number of individuals in the population */
      Individual *heap; /* a heap of individuals -- special binary tree */
      void swap(Individual &, Individual &); /* for maintaining the heap order*/
      void SiftUp(void); /* for maintaining the heap order*/
      void SiftDown(void); /* for maintaining the heap order*/
      int end_of_branch[MAX_TORS]; // For Branch Crossover Mode

   public:
      Population(void);
      Population(int); /* create a pop. with this many individuals */
      Population(int, Individual *); /* takes an array of ind's and turns into pop. */
      Population(Population &); /* copy constructor */
      ~Population(void); /* destructor */
      Individual &operator[](int);  /* for accessing a particular indiv.in pop*/
      Population &operator=(const Population &);
      unsigned int num_individuals(void); /* returns the size of the pop. */
      void msort(int); /* sorts the first m individuals using heap properties */
      // void print(ostream &, int); /* prints top int energies */
      void print(FILE *, int); /* like above */
      int printPopulationStatistics(FILE *, int, Boole); /* prints best, worse, mean, etc energies */
      int printPopulationStatisticsVerbose(FILE *, unsigned int, long int, const char []); /* print with generations & #evals */
      void printPopulationAsStates(FILE *, int, int); /*prints energies,states of top energies */
      void printPopulationAsCoordsEnergies(FILE *, int, int); /*prints energies,states of top energies */
      void set_eob(int init_end_of_branch[MAX_TORS]); // For Branch Crossover Mode
      int get_eob(int init_tor); // For Branch Crossover Mode
};

/**************************************************************************
      Inline Functions
**************************************************************************/

//  The following should be the user's default constructor.  For now,
//  we'll deal with just RealVectors
inline Genotype::Genotype(void)
{
   number_of_genes = number_of_vectors = 0;
   modified = 0;
   rep_vector = (Representation **)NULL;
   lookup = (Lookup *)NULL;
}

inline unsigned int Genotype::num_genes(void)
{
   return(number_of_genes);
}

inline unsigned int Genotype::num_vectors(void)
{
   return(number_of_vectors);
}

inline RepType Genotype::gtype(int gene_number)
{
   return(rep_vector[lookup[gene_number].vector]->type());
}

inline const Element Genotype::gread(int gene_number)
{
   return(rep_vector[lookup[gene_number].vector]->gene(lookup[gene_number].index));
}

inline const Representation *Genotype::vread(int vector_number)
{
   return(rep_vector[vector_number]);
}

//  More user definable stuff
inline Phenotype::Phenotype(void)
{
   value_vector = (Representation **)NULL;
   lookup = (Lookup *)NULL;
   number_of_dimensions = 0;
   number_of_points = 0;
   value = 0;
   evalflag = 0;
}

inline RepType Phenotype::gtype(int gene_number)
{
   return(value_vector[lookup[gene_number].vector]->type());
}

inline const Element Phenotype::gread(int gene_number)
{
   return(value_vector[lookup[gene_number].vector]->gene(lookup[gene_number].index));
}

inline const Representation *Phenotype::vread(int vector_number)
{
   return(value_vector[vector_number]);
}

inline unsigned int Phenotype::num_pts(void)
{
   return(number_of_points);
}

//  Constructs an Individual using the default constructors
inline Individual::Individual(void)
{
   age = 0L;
}


inline Individual::Individual(Individual &original)
: genotyp(original.genotyp), phenotyp(original.phenotyp)
{
}

// caution, does not do mapping
inline Individual::Individual(Genotype &init_genotyp, Phenotype &init_phenotyp)
: genotyp(init_genotyp), phenotyp(init_phenotyp)
{
}

inline Individual::~Individual(void)
{
}

inline Individual &Individual::operator=(const Individual &original)
{
   if (this == &original) {//Prevent self assignment
      return *this ;
   }
   genotyp = original.genotyp;
   phenotyp = original.phenotyp;
   mol = original.mol;
   age = original.age;
   return(*this);
}

inline double Individual::value(EvalMode mode)
{ // TO DO: check if mapping from genotyp to phenotyp is up-to-date
  // note that phenotyp.evaluate only does evaluation if evalflag is false
   return(phenotyp.evaluate(mode));
}

inline Population::Population(void)
:lhb(-1), size(0), heap((Individual *)NULL)
{
    for (int i=0; i<MAX_TORS; i++) {
        end_of_branch[i] = -1;
    }
}

inline Population::Population(int num_inds)
: lhb(num_inds-1), size(num_inds)
{
   heap = new Individual[num_inds];
    for (int i=0; i<MAX_TORS; i++) {
        end_of_branch[i] = -1;
    }
}

inline Population::Population(int newpopsize, Individual *newpop)
: size(newpopsize), heap(newpop)
{
   //  Do initialization stuff
    for (int i=0; i<MAX_TORS; i++) {
        end_of_branch[i] = -1;
    }
}

inline Population::~Population(void)
{
   if(heap != (Individual *)NULL)
   {
      delete [] heap;
   }
}

inline unsigned int Population::num_individuals(void)
{
   return(size);
}

#endif
