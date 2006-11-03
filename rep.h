//  These are the classes  associated with the Representation class
//  hierarchy.  The Representation class is meant to be a generic
//  place holder for any type of Representation that a user might
//  need to build a problem out of.  By the way, these derived 
//  (Representation) classes look like perfect candidates for templates 
//  We need to make sure that the const pointers
//  are never used to indirectly change values!
//  rsh 07/08/95
//  
//  The type Element was added to take care of all of the void *
//  rsh 02/23/96

#ifndef _REP_H
#define _REP_H

// appears not to be used //  #include <iostream.h>

#include <stdio.h>
#include "structs.h"

#define ACCURACY 0.001
#define REALV_LOW -999.999 //gmm, 2003-11-11
#define REALV_HIGH 999.999 //gmm, 2003-11-11
// #define REALV_LOW -100.0 //gmm, 2003-10-13
// #define REALV_HIGH 100.0 //gmm, 2003-10-13
// #define REALV_LOW -3.14159265358979323846 //gmm, 1998-07-08
// #define REALV_HIGH 3.14159265358979323846 //gmm, 1998-07-08

enum RepType { T_BASE, T_IntV, T_RealV, T_CRealV, T_BitV };

typedef union 
{
   double real;
   FourByteLong integer;
   unsigned char bit;
} Element;

class Representation
{
   protected:
      unsigned int number_of_pts;
      RepType mytype;
      unsigned char normalized; // =1 means the vector's normalized

   public:
      Representation(void);
      Representation(unsigned int);
      virtual ~Representation(void);
      virtual Representation &operator=(const Representation &) = 0;
      unsigned int number_of_points(void) const;
      int is_normalized(void) const;
      virtual RepType type(void) const; 
      virtual void write(unsigned char, int) = 0;
      virtual void write(FourByteLong, int) = 0;
      virtual void write(double, int) = 0;
//      virtual void write(const void *, int) = 0;
//      virtual const void *gene(unsigned int) const = 0;
      virtual void write(const Element, int) = 0;
      virtual const Element gene(unsigned int) const = 0;
      virtual Representation *clone(void) const = 0;
      virtual const void *internals(void) const = 0;
};

class IntVector : public Representation
{
//   friend void debug(IntVector &);
   protected:
      static FourByteLong low, high;
      FourByteLong *vector;

      const void *internals(void) const;
      Representation *clone(void) const;

   public:
      IntVector(void);
      IntVector(int);
      IntVector(int, FourByteLong *);
      IntVector(int, FourByteLong, FourByteLong);
      IntVector(const IntVector &);
      ~IntVector(void);
      void write(unsigned char, int);
      void write(FourByteLong, int);
      void write(double, int);
//      void write(const void *, int);
      void write(const Element, int);
      Representation &operator=(const Representation &);
//      const void *gene(unsigned int) const;
      const Element gene(unsigned int) const;
};

class RealVector : public Representation
{
//   friend void debug(RealVector &);
   protected:
      Real high, low;
      double *vector;

      const void *internals(void) const;
      Representation *clone(void) const;

   public:
      RealVector(void);
      RealVector(int);
      RealVector(int, double *);
      RealVector(int, double, double);
      RealVector(int, double, double, double); // use this to set the first value in the vector--useful for random quaternions
      RealVector(double, double, double, int); // use this to create a vector of length 3 with these values--useful for random quaternions
      RealVector(const RealVector &);
      ~RealVector(void);
      void write(unsigned char, int);
      void write(FourByteLong, int);
      void write(double, int);
//      void write(const void *, int);
      void write(const Element, int);
      Representation &operator=(const Representation &);
//      const void *gene(unsigned int) const;
      const Element gene(unsigned int) const;
};

//  Maybe this should be derived from RealVector
class ConstrainedRealVector : public Representation
{
//   friend debug(ConstrainedRealVector &);
   protected:
      static Real high, low;
      static double sum;
      double *vector;

      const void *internals(void) const;
      Representation *clone(void) const;
      void normalize(void) const;

   public:
      ConstrainedRealVector(void);
      ConstrainedRealVector(int);
      ConstrainedRealVector(int, double *);
      ConstrainedRealVector(int, double, double);
      ConstrainedRealVector(const ConstrainedRealVector &);
      ~ConstrainedRealVector(void);
      void write(unsigned char, int);
      void write(FourByteLong, int);
      void write(double, int);
//      void write(const void *, int);
      void write(const Element, int);
      Representation &operator=(const Representation &);
//      const void *gene(unsigned int) const;
      const Element gene(unsigned int) const;
};

class BitVector : public Representation
{
//   friend void debug(BitVector &);
   protected:
      static Real one_prob;
      unsigned char *vector;

      const void *internals(void) const;
      Representation *clone(void) const;

   public:
      BitVector(void);
      BitVector(int);
      BitVector(int, unsigned char *);
      BitVector(int, Real);
      BitVector(const BitVector &);
      ~BitVector(void);
      void write(unsigned char, int);
      void write(FourByteLong, int);
      void write(double, int);
//      void write(const void *, int);
      void write(const Element, int);
      Representation &operator=(const Representation &);
//      const void *gene(unsigned int) const;
      const Element gene(unsigned int) const;
};

/**************************************************************************
      Inline Functions
**************************************************************************/

inline Representation::Representation(void)
{
   number_of_pts = 0;
   mytype = T_BASE;
}

inline Representation::Representation(unsigned int pts)
{
   number_of_pts = pts;
}

inline Representation::~Representation(void)
{
}

inline unsigned int Representation::number_of_points(void) const
{
   return(number_of_pts);
}

inline int Representation::is_normalized(void) const
{
   return(normalized);
}

inline RepType Representation::type(void) const
{
   return(mytype);
}

inline Representation *Representation::clone(void) const
{
   return(NULL);
}

inline IntVector::IntVector(void)
: Representation(0)
{
   vector = (FourByteLong *)NULL;
   mytype = T_IntV;
}

//  This constructor does a shallow copy of the array
inline IntVector::IntVector(int num_els, FourByteLong *array)
: Representation(num_els)
{
   vector = array;
   mytype = T_IntV;
}

inline IntVector::~IntVector(void)
{
   if(vector!=(FourByteLong *)NULL)
   {
      delete [] vector;
   }
}

inline Representation *IntVector::clone(void) const
{
   return(new IntVector(*this));
}

inline RealVector::RealVector(void)
: Representation(0)
{
   low = REALV_LOW;
   high = REALV_HIGH;
   vector = (double *)NULL;
   mytype = T_RealV;
}

//  This performs a shallow copy of the array
inline RealVector::RealVector(int num_els, double *array)
: Representation(num_els)
{
   low = REALV_LOW;
   high = REALV_HIGH;
   vector = array;
   mytype = T_RealV;
}

inline RealVector::~RealVector(void)
{
   if(vector!=(double *)NULL)
   {
      delete [] vector;
   }
}

inline Representation *RealVector::clone(void) const
{
   return(new RealVector(*this));
}

inline ConstrainedRealVector::ConstrainedRealVector(void)
: Representation(0)
{
   normalized = 1;
   vector = (double *)NULL;
   mytype = T_CRealV;
}

inline ConstrainedRealVector::ConstrainedRealVector(int num_els, double *array)
:  Representation(num_els)
{
   normalized = 0;
   vector = array;
   mytype = T_CRealV;
}

inline ConstrainedRealVector::~ConstrainedRealVector(void)
{
   if(vector!=(double *)NULL)
   { 
      delete [] vector;
   }
}

inline Representation *ConstrainedRealVector::clone(void) const
{
   return(new ConstrainedRealVector(*this));
}

inline BitVector::BitVector(void)
: Representation(0)
{
   vector = (unsigned char *)NULL;
   mytype = T_BitV;
}

//  Do a shallow copy
inline BitVector::BitVector(int num_els, unsigned char *array)
: Representation(num_els)
{
   vector = array;
   mytype = T_BitV;
}

inline BitVector::~BitVector(void)
{
   if(vector!=(unsigned char *)NULL)
   {
      delete [] vector;
   }
}

inline Representation *BitVector::clone(void) const
{
   return(new BitVector(*this));
}

#endif
