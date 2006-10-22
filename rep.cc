/*

 $Id: rep.cc,v 1.9 2006/10/22 20:35:35 garrett Exp $

*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/* rep.cc */

/********************************************************************
     The methods associated with the Representation class hierarchy.

				rsh 9/95
********************************************************************/
 
// possibly unnecessary // #include <iostream.h>
#include <stdio.h>
#include <math.h>
#include <limits.h>
#include "rep.h"
#include "ranlib.h"
#include "structs.h"

extern FILE *logFile;
extern int debug;

//  Initializations
FourByteLong IntVector::low = -INT_MAX/4;
FourByteLong IntVector::high = INT_MAX/4;
/* A nonstatic data member cannot be defined outside its class:
 * Real RealVector::low = REALV_LOW;
 * Real RealVector::high = REALV_HIGH;
 */
//  For now assume that normalize handles this constraint
Real ConstrainedRealVector::low = REALV_LOW;
Real ConstrainedRealVector::high = REALV_HIGH;
double ConstrainedRealVector::sum = 1.0;
Real BitVector::one_prob = 0.5;

//  The member functions for the canonical base classes

//  This constructor is used to generate the initial (random) instances
//  of an integer vector.
IntVector::IntVector(int number_of_els)
: Representation(number_of_els)
{
   register int i;

#ifdef DEBUG
    (void)fprintf(logFile, "rep.cc/IntVector::IntVector(int number_of_els=%d) \n",number_of_els);
#endif /* DEBUG */

   mytype = T_IntV;
   vector = new FourByteLong[number_of_els];
   for (i=0; i<number_of_els; i++) {
      vector[i] = ignuin(low, high);
   }
}

IntVector::IntVector(int num_els, FourByteLong init_low, FourByteLong init_high)
: Representation(num_els)
{
   register int i;

#ifdef DEBUG
    (void)fprintf(logFile, "rep.cc/IntVector::IntVector(int num_els=%d, FourByteLong init_low=%ld, FourByteLong init_high=%ld) \n",num_els,init_low,init_high);
#endif /* DEBUG */


   mytype = T_IntV;
   vector = new FourByteLong[num_els];
   for (i=0; i<num_els; i++) {
      vector[i] = ignuin(init_low, init_high);
   }
}

//  This constructor does an actual copy of the vector.
//  We could make gains by doing reference counting, but
//  that's for the future.
IntVector::IntVector(const IntVector &original)
: Representation(original.number_of_pts)
{

#ifdef DEBUG
    (void)fprintf(logFile, "rep.cc/IntVector::IntVector(const IntVector &original) \n");
#endif /* DEBUG */

   mytype = T_IntV;
   if (original.vector!=NULL) {
      vector = new FourByteLong[number_of_pts];
   } else {
      vector = NULL;
   }

   for (register unsigned int i=0; i<number_of_pts; i++) {
      vector[i] = original.vector[i];
   }
}

void IntVector::write(unsigned char value, int gene)
{

#ifdef DEBUG
    (void)fprintf(logFile, "rep.cc/void IntVector::write(unsigned char value=%c, int gene=%d) \n",value,gene);
#endif /* DEBUG */

   (void)fprintf(logFile,"Writing a Bit to an Int!\n"); // used to be "stderr"
   (void)fprintf(logFile,"value= \"%c\", gene= %d\n", value, gene); // used to be "stderr"
}

void IntVector::write(FourByteLong value, int gene)
{

#ifdef DEBUG
    (void)fprintf(logFile, "rep.cc/void IntVector::write(FourByteLong value=%ld, int gene=%d)` \n",value,gene);
#endif /* DEBUG */

   if (value<low) {
      vector[gene] = low;
   } else if (value>high) {
      vector[gene] = high;
   } else {
      vector[gene] = value;
   }
}

void IntVector::write(double value, int gene)
{

#ifdef DEBUG
    (void)fprintf(logFile, "rep.cc/void IntVector::write(double value=%lf, int gene=%d) \n",value,gene);
#endif /* DEBUG */

   (void)fprintf(logFile,"Writing a Real to an Int!\n"); // used to be "stderr"
   (void)fprintf(logFile,"value= %lf, gene= %d\n", value, gene); // used to be "stderr"
}

/*
void IntVector::write(const void *value, int gene)
{

#ifdef DEBUG
    (void)fprintf(logFile, "rep.cc/void IntVector::write(const void *value, int gene=%d) \n",gene);
#endif // * DEBUG * /

   if (*((FourByteLong *)value)<low) {
      (void)fprintf(logFile,"Writing out-of-bounds Int!\n"); // used to be "stderr"
      vector[gene] = low;
   } else if (*((FourByteLong *)value)>high) {
      (void)fprintf(logFile,"Writing out-of-bounds Int!\n"); // used to be "stderr"
      vector[gene] = high;
   } else {
      vector[gene] = *((FourByteLong *)value);
   }
}
*/

void IntVector::write(const Element value, int gene)
{

#ifdef DEBUG
    (void)fprintf(logFile, "rep.cc/void IntVector::write(const Element value, int gene=%d) \n",gene);
#endif /* DEBUG */

   if (value.integer<low) {
      (void)fprintf(logFile,"Writing out-of-bounds Int!\n"); // used to be "stderr"
      vector[gene] = low;
   } else if (value.integer>high) {
      (void)fprintf(logFile,"Writing out-of-bounds Int!\n"); // used to be "stderr"
      vector[gene] = high;
   } else {
      vector[gene] = value.integer;
   }
}

/*
const void *IntVector::gene(unsigned int gene_number) const
{

#ifdef DEBUG
    (void)fprintf(logFile, "rep.cc/const void *IntVector::gene(unsigned int gene_number=%d) const \n",gene_number);
#endif // * DEBUG * /

   if (gene_number>=number_of_pts) {
      (void)fprintf(logFile,"Trying to access an out-of-bounds gene!\n"); // used to be "stderr"
      return(NULL);
   } else {
      return((void *)(&vector[gene_number]));
   }
}
*/

const Element IntVector::gene(unsigned int gene_number) const
{
   Element retval;

#ifdef DEBUG
    (void)fprintf(logFile, "rep.cc/const Element IntVector::gene(unsigned int gene_number=%d) const \n",gene_number);
#endif /* DEBUG */


   if (gene_number>=number_of_pts) {
      (void)fprintf(logFile,"Trying to access an out-of-bounds gene!\n"); // used to be "stderr"
      retval.integer = 0;
      return(retval);
   } else {
      retval.integer = vector[gene_number];
      return(retval);  // typecast int as Element
   }
}

const void *IntVector::internals(void) const
{

#ifdef DEBUG
    (void)fprintf(logFile, "rep.cc/const void *IntVector::internals(void) const \n");
#endif /* DEBUG */

   return((void *)(&vector[0]));
}

Representation &IntVector::operator=(const Representation &original)
{
   register unsigned int i;
   FourByteLong *array;

#ifdef DEBUG
    (void)fprintf(logFile, "rep.cc/Representation &IntVector::operator=(const Representation &original) \n");
#endif /* DEBUG */


   array = (FourByteLong *)original.internals();
   if (original.type()==T_IntV) {
      number_of_pts = original.number_of_points();
      if (vector!=NULL) {
         delete [] vector;
      }

      if (array!=NULL) {
         vector = new FourByteLong[number_of_pts];
      } else {
         vector = NULL;
      }

      for (i=0; i<number_of_pts; i++) {
         vector[i] = array[i];
      }
   } else {
      (void)fprintf(logFile,"Unable to invoke operator= because Representations don't match!\n"); // used to be "stderr"
   }

   return(*this);
}

//  This is the constructor used to initialize the (random)
//  starting population.
RealVector::RealVector(int num_els)
: Representation(num_els)
{

#ifdef DEBUG
    (void)fprintf(logFile, "rep.cc/RealVector::RealVector(int num_els=%d) \n",num_els);
#endif /* DEBUG */

   mytype = T_RealV;
   low = REALV_LOW;
   high = REALV_HIGH;
   vector = new double[num_els];
   for (; --num_els>=0;) {
      vector[num_els] = double(genunf(low, high));
#ifdef DEBUG
      (void)fprintf(logFile, "rep.cc/RealVector::RealVector(num_els)   vector[num_els] = %.3f\n", vector[num_els] );
#endif /* DEBUG */
   }
}

RealVector::RealVector(int num_els, double init_low, double init_high)
: Representation(num_els)
{

#ifdef DEBUG
    (void)fprintf(logFile, "rep.cc/RealVector::RealVector(int num_els=%d, double init_low=%lf, double init_high=%lf) \n",num_els,init_low,init_high);
#endif /* DEBUG */

   mytype = T_RealV;
   low = init_low;
   high = init_high;
   vector = new double[num_els];
   for (; --num_els>=0;) {
      vector[num_els] = double(genunf(init_low, init_high));
#ifdef DEBUG
      (void)fprintf(logFile, "rep.cc/RealVector::RealVector(num_els, init_low, init_high)   vector[num_els] = %.3f\n", vector[num_els] );
#endif /* DEBUG */
   }
}

//  Do a deep copy of the original
RealVector::RealVector(const RealVector &original)
: Representation(original.number_of_pts)
{

#ifdef DEBUG
    (void)fprintf(logFile, "rep.cc/RealVector::RealVector(const RealVector &original) \n");
#endif /* DEBUG */

   mytype = T_RealV;
   low =  original.low;
   high = original.high;
   if (original.vector!=NULL) {
      vector = new double[original.number_of_pts];
   } else {
      vector = NULL;
   }

   for (register unsigned int i=0; i<original.number_of_pts; i++) {
      vector[i] = original.vector[i];
#ifdef DEBUG
    (void)fprintf(logFile, "rep.cc/i=%d, original.number_of_pts=%d, vector[%d]= %.3f\n",i, original.number_of_pts, i, vector[i]);
#endif /* DEBUG */
   }
}

void RealVector::write(unsigned char value, int gene)
{

#ifdef DEBUG
    (void)fprintf(logFile, "rep.cc/void RealVector::write(unsigned char value=%c, int gene=%d) \n",value,gene);
#endif /* DEBUG */

   (void)fprintf(logFile,"Writing a Bit to a Real!\n"); // used to be "stderr"
   (void)fprintf(logFile,"value= \"%c\", gene= %d\n", value, gene); // used to be "stderr"
}

void RealVector::write(FourByteLong value, int gene)
{

#ifdef DEBUG
    (void)fprintf(logFile, "rep.cc/void RealVector::write(FourByteLong value=%ld, int gene=%d) \n",value,gene);
#endif /* DEBUG */

   (void)fprintf(logFile,"Writing an Int to a Real!\n"); // used to be "stderr"
   (void)fprintf(logFile,"value= %ld, gene= %d\n", value, gene); // used to be "stderr"
}

void RealVector::write(double value, int gene)
{

#ifdef DEBUG
    (void)fprintf(logFile, "rep.cc/void RealVector::write(double value=%lf, int gene=%d) \n",value,gene);
#endif /* DEBUG */

   if (value<low) {
      // if (debug > 0) {
          // (void)fprintf(logFile,"WARNING:  Writing out of bounds Real!  value (%lf) too low (%lf)\n",value,low); // used to be "stderr"
      // }
      vector[gene] = low;
   } else if (value>high) {
      // if (debug > 0) {
          // (void)fprintf(logFile,"WARNING:  Writing out of bounds Real!  value (%lf) too high (%lf)\n",value,high); // used to be "stderr"
      // }
      vector[gene] = high;
   } else {
      vector[gene] = value;
   }
}

/*
void RealVector::write(const void *value, int gene)
{

#ifdef DEBUG
    (void)fprintf(logFile, "rep.cc/void RealVector::write(const void *value, int gene=%d) \n",gene);
#endif // * DEBUG * /

   if (*((double *)value)<low) {
      vector[gene] = low;
   } else if (*((double *)value)>high) {
      vector[gene] = high;
   } else {
      vector[gene] = *((double *)value);
   }
}
*/

void RealVector::write(const Element value, int gene)
{

#ifdef DEBUG
    (void)fprintf(logFile, "rep.cc/void RealVector::write(const Element value, int gene=%d) \n",gene);
#endif /* DEBUG */

   if (value.real<low) {
      vector[gene] = low;
   } else if (value.real>high) {
      vector[gene] = high;
   } else {
      vector[gene] = value.real;
   }
}

/*
 *const void *RealVector::gene(unsigned int gene_number) const
 *{
 *
 *#ifdef DEBUG
     *(void)fprintf(logFile, "rep.cc/const void *RealVector::gene(unsigned int gene_number=%d) const \n",gene_number);
 *#endif // * DEBUG * /
 *
    *if (gene_number>=number_of_pts) {
       *(void)fprintf(logFile,"Trying to access out-of-bounds gene\n"); // used to be "stderr"
       *return(NULL);
    *} else {
       *return((void *)(&vector[gene_number]));
    *}
 *}
 */

const Element RealVector::gene(unsigned int gene_number) const
{

#ifdef DEBUG
    (void)fprintf(logFile, "rep.cc/const Element RealVector::gene(unsigned int gene_number=%d) const \n",gene_number);
#endif /* DEBUG */

   Element retval;

   if (gene_number>=number_of_pts) {
      (void)fprintf(logFile,"Trying to access out-of-bounds gene\n"); // used to be "stderr"
      retval.real = 0.0;
      return(retval);
   } else {
      retval.real = vector[gene_number];
      return(retval);
   }
}

const void *RealVector::internals(void) const
{

#ifdef DEBUG
    (void)fprintf(logFile, "rep.cc/const void *RealVector::internals(void) const \n");
#endif /* DEBUG */

   return((void *)(&vector[0]));
}

Representation &RealVector::operator=(const Representation &original)
{

#ifdef DEBUG
    (void)fprintf(logFile, "rep.cc/Representation &RealVector::operator=(const Representation &original) \n");
#endif /* DEBUG */

   register unsigned int i;
   double *array;

   if (original.type()==T_RealV) {
      low = ((const RealVector &)original).low;
      high = ((const RealVector &)original).high;
      array = (double *)original.internals();
      number_of_pts = original.number_of_points();
      if (vector!=NULL) {
         delete [] vector;
      }

      if (array!=NULL) {
         vector = new double[number_of_pts];
      } else {
         vector = NULL;
      }

      for (i=0; i<number_of_pts; i++) {
         vector[i] = array[i];
      }
   } else {
      (void)fprintf(logFile,"Unable to invoke operator= because Representations don't match!\n"); // used to be "stderr"
   }

   return(*this);
}

ConstrainedRealVector::ConstrainedRealVector(int num_els)
:  Representation(num_els)
{

#ifdef DEBUG
    (void)fprintf(logFile, "rep.cc/ConstrainedRealVector::ConstrainedRealVector(int num_els) \n");
#endif /* DEBUG */

   mytype = T_CRealV;
   normalized = 0;
   vector = new double[num_els];
   for (; --num_els>=0;) {
      vector[num_els] = double(genunf(low, high));
   }

//   normalize();
}

ConstrainedRealVector::ConstrainedRealVector(int num_els, double init_low, double init_high)
:  Representation(num_els)
{

#ifdef DEBUG
    (void)fprintf(logFile, "rep.cc/ConstrainedRealVector::ConstrainedRealVector(int num_els=%d, double init_low=%lf, double init_high=%lf) \n",num_els,init_low,init_high);
#endif /* DEBUG */

   mytype = T_CRealV;
   normalized = 0;
   vector = new double[num_els];
   for (; --num_els>=0;) {
      vector[num_els] = double(genunf(init_low, init_high));
   }

//   normalize();
}

ConstrainedRealVector::ConstrainedRealVector(const ConstrainedRealVector &original)
:  Representation(original.number_of_pts)
{

#ifdef DEBUG
    (void)fprintf(logFile, "rep.cc/ConstrainedRealVector::ConstrainedRealVector(const ConstrainedRealVector &original) \n");
#endif /* DEBUG */

   mytype = T_CRealV;
   normalized = original.normalized;
   if (original.vector != NULL) {
      vector = new double[original.number_of_pts];
   } else {
      vector = NULL;
   }

   for (register unsigned int i=0; i<original.number_of_pts; i++) {
      vector[i] = original.vector[i];
   }
}

void ConstrainedRealVector::write(unsigned char value, int gene)
{

#ifdef DEBUG
    (void)fprintf(logFile, "rep.cc/void ConstrainedRealVector::write(unsigned char value=%c, int gene=%d) \n",value,gene);
#endif /* DEBUG */

   (void)fprintf(logFile,"Writing a Bit to a Constrained Real\n"); // used to be "stderr"
   (void)fprintf(logFile,"value= \"%c\",  gene= %d\n", value, gene); // used to be "stderr"
}

void ConstrainedRealVector::write(FourByteLong value, int gene)
{

#ifdef DEBUG
    (void)fprintf(logFile, "rep.cc/void ConstrainedRealVector::write(FourByteLong value=%ld, int gene=%d) \n",value,gene);
#endif /* DEBUG */

   (void)fprintf(logFile,"Writing an Integer to a Constrained Real\n"); // used to be "stderr"
   (void)fprintf(logFile,"value= %ld, gene= %d\n",value,gene); // used to be "stderr"
}

void ConstrainedRealVector::write(double value, int gene)
{

#ifdef DEBUG
    (void)fprintf(logFile, "rep.cc/void ConstrainedRealVector::write(double value=%lf, int gene=%d) \n",value,gene);
#endif /* DEBUG */

   if (value<low) {
      (void)fprintf(logFile,"Writing out-of-bounds Constrained Real\n"); // used to be "stderr"
      vector[gene] = low;
   } else if (value>high) {
      (void)fprintf(logFile,"Writing out-of-bounds Constrained Real\n"); // used to be "stderr"
      vector[gene] = high;
   } else {
      vector[gene] = value;
   }

   normalized = 0;
}

/*
void ConstrainedRealVector::write(const void *value, int gene)
{

#ifdef DEBUG
    (void)fprintf(logFile, "rep.cc/void ConstrainedRealVector::write(const void *value, int gene=%d) \n",gene);
#endif // * DEBUG * /

   if (*((double *)value)<low) {
      vector[gene] = low;
   } else if (*((double *)value)>high) {
      vector[gene] = high;
   } else {
      vector[gene] = *((double *)value);
   }

   normalized = 0;
}
*/

void ConstrainedRealVector::write(const Element value, int gene)
{

#ifdef DEBUG
    (void)fprintf(logFile, "rep.cc/void ConstrainedRealVector::write(const Element value, int gene=%d) \n",gene);
#endif /* DEBUG */

   if (value.real<low) {
      vector[gene] = low;
   } else if (value.real>high) {
      vector[gene] = high;
   } else {
      vector[gene] = value.real;
   }

   normalized = 0;
}

void ConstrainedRealVector::normalize(void) const
{
   //unsigned char *kluge; commented out by gmm, 9-17-97

#ifdef DEBUG
    (void)fprintf(logFile, "rep.cc/void ConstrainedRealVector::normalize(void) const \n");
#endif /* DEBUG */


//   kluge = &normalized;
   if (!normalized) {
      register unsigned int i;
      register double tempsum = 0.0, hypotenuse;

      for (i=0; i<number_of_pts; i++) {
         tempsum+=vector[i]*vector[i];
      }

      if ((tempsum-sum>ACCURACY)||(sum-tempsum>ACCURACY)) {
         hypotenuse = sqrt(tempsum);
         for (i=0; i<number_of_pts; i++) {
            vector[i] /= hypotenuse;
         }
      }

      //normalized = 1;
      //*kluge = 1; commented out by gmm, 9-17-97
   }
}

/*
const void *ConstrainedRealVector::gene(unsigned int gene_number) const
{

#ifdef DEBUG
    (void)fprintf(logFile, "rep.cc/const void *ConstrainedRealVector::gene(unsigned int gene_number=%d) const \n",gene_number);
#endif // * DEBUG * /

   if (gene_number>=number_of_pts) {
      (void)fprintf(logFile,"Trying to access an out-of-bounds gene\n"); // used to be "stderr"
      return(NULL);
   } else {
//      normalize();
      return((void *)(&vector[gene_number]));
   }
}
*/

const Element ConstrainedRealVector::gene(unsigned int gene_number) const
{
   Element retval;

#ifdef DEBUG
    (void)fprintf(logFile, "rep.cc/const Element ConstrainedRealVector::gene(unsigned int gene_number=%d) const \n",gene_number);
#endif /* DEBUG */


   if (gene_number>=number_of_pts) {
      (void)fprintf(logFile,"Trying to access an out-of-bounds gene\n"); // used to be "stderr"
      retval.real = 0.0;
      return(retval);
   } else {
//      normalize();
      retval.real = 0.0;
      return(retval);
   }
}

const void *ConstrainedRealVector::internals(void) const
{

#ifdef DEBUG
    (void)fprintf(logFile, "rep.cc/const void *ConstrainedRealVector::internals(void) const \n");
#endif /* DEBUG */

   return((void *)(&vector[0]));
}

Representation &ConstrainedRealVector::operator=(const Representation &original)
{
   register unsigned int i;
   double *array;

#ifdef DEBUG
    (void)fprintf(logFile, "rep.cc/Representation &ConstrainedRealVector::operator=(const Representation &original) \n");
#endif /* DEBUG */


   array = (double *)original.internals();
   if (original.type()==T_CRealV) {
      number_of_pts = original.number_of_points();
      normalized = original.is_normalized();
      if (vector!=NULL) {
         delete [] vector;
      }
      
      if (array!=NULL) {
         vector = new double[number_of_pts];
      } else {
         vector = NULL;
      }
      
      for (i=0; i<number_of_pts; i++) {
         vector[i] = array[i];
      }
   } else {
      (void)fprintf(logFile,"Unable to invoke operator= because Representations don't match!\n"); // used to be "stderr"
   }

   return(*this);
}

//  This constructor is used to initialize the first
//  generation of any particular bitvector.  Right
//  now bits are assumed to be unsigned chars.
BitVector::BitVector(int num_els)
: Representation(num_els)
{

#ifdef DEBUG
    (void)fprintf(logFile, "rep.cc/BitVector::BitVector(int num_els=%d) \n",num_els);
#endif /* DEBUG */

   mytype = T_BitV;
   vector = new unsigned char[num_els];
   for (; --num_els>=0;) {
      vector[num_els] = ((ranf()<one_prob)? 1 : 0);
   }
}

BitVector::BitVector(int num_els, Real prob)
: Representation(num_els)
{

#ifdef DEBUG
    (void)fprintf(logFile, "rep.cc/BitVector::BitVector(int num_els=%d, Real prob=%f) \n",num_els,prob);
#endif /* DEBUG */

   mytype = T_BitV;
   vector = new unsigned char[num_els];
   for (; --num_els>=0;) {
      vector[num_els] = ((ranf()<prob)? 1 : 0);
   }
}

//  There are probably better ways of doing this, e.g.
//  using memcpy()
BitVector::BitVector(const BitVector &original)
: Representation(original.number_of_pts)
{

#ifdef DEBUG
    (void)fprintf(logFile, "rep.cc/BitVector::BitVector(const BitVector &original) \n");
#endif /* DEBUG */

   mytype = T_BitV;
   if (original.vector!=NULL) {
      vector = new unsigned char[number_of_pts];
   } else {
      vector = NULL;
   }

   for (register unsigned int i=0; i<number_of_pts; i++) {
      vector[i] = original.vector[i];
   }
}

void BitVector::write(unsigned char value, int gene)
{

#ifdef DEBUG
    (void)fprintf(logFile, "rep.cc/void BitVector::write(unsigned char value=%c, int gene=%d) \n",value,gene);
#endif /* DEBUG */

   vector[gene] = value;
}

void BitVector::write(FourByteLong value, int gene)
{

#ifdef DEBUG
    (void)fprintf(logFile, "rep.cc/void BitVector::write(FourByteLong value=%ld, int gene=%d) \n",value,gene);
#endif /* DEBUG */

   (void)fprintf(logFile,"Writing Int to Bit!\n"); // used to be "stderr"
   (void)fprintf(logFile,"value= %ld, gene= %d\n",value,gene); // used to be "stderr"
}

void BitVector::write(double value, int gene)
{

#ifdef DEBUG
    (void)fprintf(logFile, "rep.cc/void BitVector::write(double value=%lf, int gene=%d) \n",value,gene);
#endif /* DEBUG */

   (void)fprintf(logFile,"Writing Real to Bit!\n"); // used to be "stderr"
   (void)fprintf(logFile,"value= %lf, gene= %d\n",value,gene); // used to be "stderr"
}

/*
void BitVector::write(const void *value, int gene)
{

#ifdef DEBUG
    (void)fprintf(logFile, "rep.cc/void BitVector::write(const void *value, int gene=%d) \n",gene);
#endif // * DEBUG * /

   vector[gene] = *((unsigned char *)value);
}
*/

void BitVector::write(const Element value, int gene)
{

#ifdef DEBUG
    (void)fprintf(logFile, "rep.cc/void BitVector::write(const Element value, int gene=%d) \n",gene);
#endif /* DEBUG */

   vector[gene] = value.bit;
}

/*
const void *BitVector::gene(unsigned int gene_number) const
{

#ifdef DEBUG
    (void)fprintf(logFile, "rep.cc/const void *BitVector::gene(unsigned int gene_number=%d) const \n",gene_number);
#endif // * DEBUG * /

   if (gene_number>=number_of_pts) {
      (void)fprintf(logFile,"Trying to access an out-of-bounds gene\n"); // used to be "stderr"
      return(NULL);
   } else {
      return((void *)(&vector[gene_number]));
   }
}
*/

const Element BitVector::gene(unsigned int gene_number) const
{
   Element retval;

#ifdef DEBUG
    (void)fprintf(logFile, "rep.cc/const Element BitVector::gene(unsigned int gene_number=%d) const \n",gene_number);
#endif /* DEBUG */


   if (gene_number>=number_of_pts) {
      (void)fprintf(logFile,"Trying to access an out-of-bounds gene\n"); // used to be "stderr"
      retval.bit = 0;
      return(retval);
   } else {
      retval.bit = vector[gene_number];
      return(retval);
   }
}

const void *BitVector::internals(void) const
{

#ifdef DEBUG
    (void)fprintf(logFile, "rep.cc/const void *BitVector::internals(void) const \n");
#endif /* DEBUG */

   return((void *)(&vector[0]));
}

Representation &BitVector::operator=(const Representation &original)
{
   register unsigned int i;
   unsigned char *array;

#ifdef DEBUG
    (void)fprintf(logFile, "rep.cc/Representation &BitVector::operator=(const Representation &original) \n");
#endif /* DEBUG */


   array = (unsigned char *)original.internals();
   if (original.type()==T_BitV) {
      if (vector!=NULL) {
         delete [] vector;
      }

      number_of_pts = original.number_of_points();
      if (array!=NULL) {
         vector = new unsigned char[number_of_pts];
      } else {
         vector = NULL;
      }

      for (i=0; i<number_of_pts; i++) {
         vector[i] = array[i];
      }
   } else {
      (void)fprintf(logFile,"Unable to invoke operator= because Representations don't match!\n"); // used to be "stderr"
   }

   return(*this);
}

