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
** Revision 1.1  2001/08/13 22:05:53  gillet
** *** empty log message ***
**
*/

#ifndef _SUPPORT_H
#define _SUPPORT_H

#include <stdio.h>
#include "rep.h"
#include "eval.h"
#include "structs.h"

enum EvalMode {Reset, Always_Eval, Normal_Eval};

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
      Genotype(unsigned int, Representation **); /* creates a genotype from the
					     representation & total # vectors */
      ~Genotype(void);
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
};

//  Should Phenotype automatically evaluate itself upon construction?
class Phenotype
{
   //friend void debug(Phenotype &);
   protected:
      unsigned int number_of_dimensions, number_of_points;
      double value;
      Lookup *lookup;
      Representation **value_vector;
      unsigned evalflag : 1;  //  =1 means that this is the current evaluation
      unsigned modified : 1;  //  =1 means that this has been modified

   public:
      Phenotype(void);
      Phenotype(const Phenotype &);
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
      Individual(Genotype &, Phenotype &);
      ~Individual(void); /* destructor */
      Individual &operator=(const Individual &); /* assignment function for
						    individuals */
      Phenotype mapping(void); /* takes the genotype and converts it into a
				  phenotype.  */
      Genotype inverse_mapping(void);  // Scott should do: Also copy Phenotype's value
      double value(EvalMode); /* evaluation of the individual gives its value */
      State state(int); /* state variables in AutoDock */
      void  getMol(Molecule *); /* converts phenotype to mol's state and returns this individual's mol data */
      void printIndividualsState(FILE *, int); /* print out the state of this individual */
      void incrementAge(); /* make individual grow 1 generation older */
};

class Population
{
   //friend void debug(Population &);
   protected:
      Individual *heap; /* a heap of individuals -- special binary tree */
      int lhb;  //  These keep track of the lower & upper heap bounds
      int size; /* the number of individuals in the population */
      void swap(Individual &, Individual &); /* for maintaining the heap order*/
      void SiftUp(void); /* for maintaining the heap order*/
      void SiftDown(void); /* for maintaining the heap order*/

   public:
      Population(void);
      Population(int); /* create a pop. with this many individuals */
      Population(int, Individual *); /* takes an array of ind's and turns into pop. */
      Population(Population &); /* copy constructor */
      ~Population(void);
      Individual &operator[](int);  /* for accessing a particular indiv.in pop*/
      Population &operator=(const Population &);
      unsigned int num_individuals(void); /* returns the size of the pop. */
      void msort(int); /* sorts the first m individuals using heap properties */
      void print(ostream &, int); /* prints top int energies */
      void print(FILE *, int); /* like above */
      void printPopulationAsStates(FILE *, int, int); /*prints energies,states of top energies */
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

inline Individual::Individual(Genotype &init_genotyp, Phenotype &init_phenotyp)
: genotyp(init_genotyp), phenotyp(init_phenotyp)
{
}

inline Individual::~Individual(void)
{
}

inline Individual &Individual::operator=(const Individual &original)
{
   genotyp = original.genotyp;
   phenotyp = original.phenotyp;
   mol = original.mol;
   age = original.age;

   return(*this);
}

inline double Individual::value(EvalMode mode)
{
   return(phenotyp.evaluate(mode));
}

inline Population::Population(void)
:heap((Individual *)NULL), lhb(-1), size(0)
{
}

inline Population::Population(int num_inds)
: lhb(num_inds-1), size(num_inds)
{
   heap = new Individual[num_inds];
}

inline Population::Population(int newpopsize, Individual *newpop)
: size(newpopsize), heap(newpop)
{
   //  Do initialization stuff
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
