/*

 $Id: support.cc,v 1.33 2010/01/08 20:13:47 mp Exp $

 AutoDock 

Copyright (C) 2009 The Scripps Research Institute. All rights reserved.

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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "eval.h"
#include "support.h"
#include "stateLibrary.h"
#include "structs.h"

//#define DEBUG 
extern FILE *logFile;

extern class Eval evaluate;


//  These are the member functions for the support classes.


Population::Population(Population &original)
: lhb(original.lhb), size(original.size)
{
   register int i;

#ifdef DEBUG
   (void)fprintf(logFile, "support.cc/Population::Population(Population &original)\n");
#endif /* DEBUG */

   heap = new Individual[size];
   for (i=0; i<size; i++) {
      heap[i] = original.heap[i];
      heap[i].age = 0L; // gmm, 1998-07-14
   }
}

/*  Heap Functions:  In this case, the heap condition means the maximal 
    element wrt fitness (i.e. the best) is at the top of the heap.  lhb 
    is the index of the last element to be inserted into the heap.  Some 
    the standard functions on the heap can be accomplished in the following 
    manner:
       Root = size - 1  (Note: that the root is fixed)
       LeftChild(I) = 2I - size
       RightChild(I) = 2I - size - 1
       Parent(I) = (I + size + 1)/2
    It is important to notice that the heap condition is maintained from
    lhb to size-1 *in reverse order*.
*/

void Population::swap(Individual &individual1, Individual &individual2)
{
   Individual temp;

#ifdef DEBUG
   (void)fprintf(logFile, "support.cc/void Population::swap(Individual &individual1, Individual &individual2)\n");
#endif /* DEBUG */


   temp = individual1;
   individual1 = individual2;
   individual2 = temp;
}

/*  This routine assumes that the heap condition is satisfied
    between lhb and size-1 and that the new individual is in
    position lhb-1.
*/
void Population::SiftUp(void)
{
   int i, parent;

#ifdef DEBUG
   (void)fprintf(logFile, "support.cc/void Population::SiftUp(void)\n");
#endif /* DEBUG */


   i = lhb-1;
   while (((parent=(i+size+1)/2)<size)&&(heap[parent].value(Normal_Eval)>heap[i].value(Normal_Eval))) {
      swap(heap[parent], heap[i]);
      i = parent;
   }
   lhb--;
}

/*  This routine assumes that the heap condition is satisfied
    between lhb & size-2 initially, and that the individual at size-1
    needs to be accomodated.*/
void Population::SiftDown(void)
{
   int i, child;

#ifdef DEBUG
   (void)fprintf(logFile, "support.cc/void Population::SiftDown(void)\n");
#endif /* DEBUG */


   i = size-1;
   while ((child=2*i-size)>=lhb) {
      if (child-1>=lhb) {
         if (heap[child-1].value(Normal_Eval)<heap[child].value(Normal_Eval)) {
            child--;
         }
      }
      /*  Now child holds the index of the best child of i  */
      if (heap[i].value(Normal_Eval)<heap[child].value(Normal_Eval)) {
         break;
      } else {
         swap(heap[child], heap[i]);
         i = child;
      }
   }
}

void Population::msort(int m)
{
   register int i;

#ifdef DEBUG
   (void)fprintf(logFile, "support.cc/void Population::msort(int m=%d)\n", m);
#endif /* DEBUG */

   

   //  First make a heap of the whole array, i.e lhb = 0 & uhb = size
   lhb = size-1;
   while (lhb>0) {
      SiftUp();
   }

   assert(lhb==0);
 
   //  Now place the m best members at the beginning of the array
   if (m==size) return; //we're done
   if (m>size) {
        (void)fprintf(logFile, "support.cc/Population::msort(int m=%d) -- ERROR!  m > size!\n\n", m);
        exit(-1);
   }

   for (i=0; i<m && i<size-1; i++) {
#ifdef DEBUG
   (void)fprintf(stderr, "support.cc/placing %d of %d best of %d, lhb=%d )\n", i, m, size, lhb);
#endif /* DEBUG */
      swap(heap[i], heap[size-1]);
      lhb++;
#ifdef DEBUG
   (void)fprintf(stderr, "support.cc/calling SiftDown lhb=%d )\n", lhb);
#endif /* DEBUG */
      SiftDown();
   }
   
   //  Assert: heap[0..m-1] sorted
}

//void Population::print(ostream &output, int num)
//{
   //register int i;

//#ifdef DEBUG
   //(void)fprintf(logFile, "support.cc/void Population::print(ostream &output, int num=%d)\n", num);
//#endif /* DEBUG */


   //(void)fprintf(logFile, "The top %d individuals in the population:\n", num);
   //for (i=0; i<num; i++) {
      //(void)fprintf(logFile, "%lf\n", heap[i].value(Normal_Eval));
   //}
//}

void Population::print(FILE *output, int num) {
   register int i;

#ifdef DEBUG
   (void)fprintf(logFile, "support.cc/void Population::print(FILE *output, int num=%d)\n", num);
#endif /* DEBUG */

   (void)fprintf( output, "The top %d individuals in the population:\n\n", num);
   for (i=0; i<num; i++) {
      (void)fprintf( output, "(%d):\t %8.2f\n", i+1, heap[i].value(Always_Eval));
   }
}

static int
doublecompare(const void *p1, const void *p2)
{
	double i = *((double *)p1);
	double j = *((double *)p2);
	if ( i > j ) return 1;
	if ( i < j ) return -1;
	return 0;
	}
static double
compute_k_quantile(const int k, const int q, 
 const double energy[], const int size)
{
 // return k-th (0-origin) "q"-tile of array (assumed sorted ascending)
 // example: median is k=1  q=2
 // example: first quartile is k=1 q=4
 //
 // If the computed index falls between two input sample values, the
 //   weighted average of the two is returned.
 //
 // M Pique December 2009 from wikipedia "Quantile" 
 //    "Estimating the quantiles of a population", "Weighted average"
 //    http://en.wikipedia.org/wiki/Quantile 
double p = (size-1) * k / q;
double j_double; // integer part of p
double g; // fractional part of p
int ilow, ihigh; // indices to be weighted
  g = modf(p, &j_double);
  ilow = (int) j_double;
  assert(ilow>=0 && ilow<size);
  if (g == 0.0) return energy[ilow];
  ihigh=ilow+1;
  assert(ihigh<size);
  return energy[ilow] + g*(energy[ihigh]-energy[ilow]);
  }


  
   

int Population::printPopulationStatistics(FILE *output, int level, Boole appendNewline) {
// write best energy, etc, depending on level
// return 0 if OK, non-zero if error
// Code adapted from gs.cc  - M Pique  December 2009
int returnCode=0;
double sum, best_e, worst_e;
unsigned long oldest;
int oldestIndividual, best_i;
   if(size<=0) return 2; // error, give up
   sum=best_e=worst_e=heap[best_i=0].value(Normal_Eval);
   oldest=heap[oldestIndividual=0].age;
   for (int i=1; i<size; i++) {
      double e = heap[i].value(Normal_Eval);
      sum+=e;
      if(e>worst_e) worst_e = e;
      if(e<best_e) {
          best_e = e;
	  best_i = i;
	  }
      if (heap[i].age >= oldest) {
          oldest = heap[i].age;
          oldestIndividual = i;
      }
   }

   // level == 1 - short output as gs.cc did previously: print oldest and best (in exactly same format)
   // level == 2 - add age info (supporting DEBUG3 option in gs.cc)
   // level == 3 - no age info, but add worst, mean, median, quartiles, 
   //    standard deviation to (1) in simpler format
switch (level) {
    case 1:
	    (void)fprintf(output, " Oldest's energy: %.3f    Lowest energy: %.3f", 
               heap[oldestIndividual].value(Normal_Eval), best_e);
	       break;
    case 2:
       (void)fprintf(output, " Oldest ind.: %u/%u, age: %lu, energy: %.3f    Lowest energy individual: %u/%u, age: %lu, energy: %.3f", 
               oldestIndividual+1, size, 
	       heap[oldestIndividual].age, 
               heap[oldestIndividual].value(Normal_Eval), 
	       best_i+1, size, 
               heap[best_i].age, best_e);
	       break;
    default:
	{
	double mean, sum_squares, median, stddev;
	//double q14, q34; // 1st and 3rd quartiles
	//double q15, q45; // 1st and 4th quintiles
	double energy[size]; // array for sorting
	mean = sum/size;
	sum_squares=0;
	for (int i=0; i<size; i++) energy[i] = heap[i].value(Normal_Eval);
	for (int i=0; i<size; i++) sum_squares += (energy[i]-mean)*(energy[i]-mean);
	if(size<2) stddev=0; else stddev=sqrt(sum_squares/(size-1));
	// compute median and quartiles 
	//  Avoid using msort because reordering the actual population
	//  causes runs to produce different results depending on statistics output level
	qsort( (void *) energy, size, sizeof(*energy), doublecompare); // ascend

	if( 1 == (size%2) ) median = energy[size/2]; // odd size
	else median = (energy[size/2] + energy[size/2-1])/2.0; // even size

	(void) fprintf(output, "Lowest: %.3f Highest: %.3f Mean: %.3f Median: %.3f Std.Dev: %.3f", 
	  best_e, worst_e, mean, median, stddev);

	// quartiles and quintiles
	//q15 = compute_k_quantile(1, 5, energy, size);
	//q14 = compute_k_quantile(1, 4, energy, size);
	//q34 = compute_k_quantile(3, 4, energy, size);
	//q45 = compute_k_quantile(4, 5, energy, size);
#define quantile(k, q) \
	(void) fprintf(output, " Q%d/%d: %.3f", \
	  k, q, compute_k_quantile(k, q, energy, size))

	quantile(1,10);
	quantile(1, 5);
	quantile(1, 4);
	quantile(3, 4);
	quantile(4, 5);
	quantile(9,10);

#undef quantile
	  // debug print every energy:
	  if(level>3) for(int i=0; i<size; i++) fprintf(output, " %.3f", energy[i]);
	  }
	break;
	}
    if(appendNewline) fprintf(output, "\n");
    return returnCode;
 }
int Population::printPopulationStatisticsVerbose(FILE * output, 
 unsigned int generations, long int nevals, const char suffix[]){ /* print with generations & #evals */
int returnCode=0;
   // print "Population at Generation:" line with low/high/mean/median/stddev...
   (void) fprintf(output, "Population at Generation: %3u ", generations);
   // highest level, with newline at end:
   returnCode= printPopulationStatistics(logFile, 3, FALSE);
   (void) fprintf(logFile, " Num.evals: %ld%s", nevals, suffix );
    return returnCode; 
}
 	

void Population::printPopulationAsStates(FILE *output, int num, int ntor) {
   register int i;
#ifdef DEBUG2
   register int j;
   char resstr[LINE_LEN];
#endif /* DEBUG2 */
   double thisValue;

#ifdef DEBUG
   (void)fprintf(logFile, "support.cc/void Population::printPopulationAsStates(FILE *output, int num=%d, int ntor=%d)\n", num, ntor);
#endif /* DEBUG */

   // Print an XML-like tag indicating this is a population, with attribute size
   // being the number of individuals in the population
   (void)fprintf( output, "<population size=\"%d\">\n", num);
   for (i=0; i<num; i++) {
      thisValue = heap[i].value(Always_Eval);
      (void)fprintf( output, "%4d\t%9.4lg\t", i+1, thisValue);
      (void)fprintf( output, "%4lu\t", heap[i].age );
      heap[i].printIndividualsState(output, ntor, 0);

#ifdef DEBUG2
      // to print only infinite or NaN structures // if (!finite(thisValue) || ISNAN(thisValue)) {//debug
      // Convert state to coords and print it out...//debug
      cnv_state_to_coords(heap[i].state(ntor), heap[i].mol->vt,  heap[i].mol->tlist,  ntor, heap[i].mol->crdpdb,  heap[i].mol->crd,  heap[i].mol->natom);//debug
      (void)fprintf(logFile, "MODEL     %4d\n", i+1);
      for (j=0; j<heap[i].mol->natom; j++) {//debug
            print_PDBQ_atom_resnum( logFile, "", j, " C   RES  ",   1,  // replace 1 with i+1 for incrementing residue numbers.
            heap[i].mol->crd, 0.0, 0.0, 0.0, "\n"); //debug
      }/*j*///debug
      (void)fprintf(logFile, "ENDMDL\n");
      // to print only infinite or NaN structures // }// thisValue is either infinite or not-a-number.//debug
#endif /* DEBUG2 */

   }// i
   (void)fprintf( output, "</population>\n");
}

void Population::printPopulationAsCoordsEnergies(FILE *output, int num, int ntor) {
   register int i;
   double thisValue;

#ifdef DEBUG
   (void)fprintf(logFile, "support.cc/void Population::printPopulationAsCoordsEnergies(FILE *output, int num=%d, int ntor=%d)\n", num, ntor);
#endif // DEBUG

   //(void)fprintf( output, "The top %d individuals in the population:\n\n", num);
   for (i=0; i<num; i++) {

      // Print the number of this individual in the population (counting from 1, not 0)
      (void)fprintf( output, "%d\t", i+1);

      // Print the translation
      heap[i].printIndividualsState(output, ntor, 3);  // 3 means print just the translation
      //heap[i].printIndividualsState(output, ntor, 0);  // 0 means print the whole state

      // Print the energy
      thisValue = heap[i].value(Normal_Eval); // was Always_Eval before 13-Jan-2006
      (void)fprintf( output, "\t%9.2lf", thisValue);

      // Print the non-bonded energy, i.e. vdW+Hb+desolv
      thisValue = heap[i].value(Always_Eval_Nonbond);
      (void)fprintf( output, "\t%9.2lf", thisValue);
   
      // Print the electrostatic energy, i.e. elec
      thisValue = heap[i].value(Always_Eval_Elec);
      (void)fprintf( output, "\t%9.2lf", thisValue);
     
      // Write a newline at the end
      (void)fprintf( output, "\n");

      // We need the coordinates of this individual to compute the electrostatic and nonbond energies
      //cnv_state_to_coords( heap[i].state(ntor), heap[i].mol->vt,  heap[i].mol->tlist,  ntor, heap[i].mol->crdpdb,  heap[i].mol->crd,  heap[i].mol->natom);

   }// i
   (void)fprintf( output, "\n");
}

void Population::set_eob(int init_end_of_branch[MAX_TORS])
// Set the end_of_branch[MAX_TORS] array
{
   for (register int i=0; i<MAX_TORS; i++) {
       end_of_branch[i] = init_end_of_branch[i];
   }
}

int Population::get_eob(int init_tor)
// Get the end_of_branch[] value for the supplied torsion number, init_tor
{
    if ((init_tor >= 0) && (init_tor < MAX_TORS)) {
        return end_of_branch[init_tor];
    } else {
        (void)fprintf(logFile, "support.cc/Population::get_eob(int init_tor=%d) -- ERROR!  Attempt to access out-of-bounds torsion!\n\n", init_tor);
        exit(-1);
    }
}

Genotype::Genotype(unsigned int init_number_of_vectors, Representation **
init_rep_vector)
: number_of_vectors(init_number_of_vectors), rep_vector(init_rep_vector), 
  modified(0)
{
   register unsigned int i, j, k;

#ifdef DEBUG
   (void)fprintf(logFile, "support.cc/Genotype::Genotype(unsigned int init_number_of_vectors=%d, Representation **init_rep_vector)\n", init_number_of_vectors);
#endif /* DEBUG */


   number_of_genes = 0;
   for (i=0; i<number_of_vectors; i++) {
      number_of_genes += rep_vector[i]->number_of_points();
#ifdef DEBUG
      (void)fprintf(logFile, "support.cc/Genotype::Genotype(init_number_of_vectors=%d, **init_rep_vector) number_of_genes=%d   rep_vector[%d]->number_of_points()=%d\n", init_number_of_vectors, number_of_genes, i, rep_vector[i]->number_of_points());
#endif /* DEBUG */
   }
 
   i=0;
   lookup = new Lookup[number_of_genes];
   for (j=0; j<number_of_vectors; j++) {
      for (k=0; k<rep_vector[j]->number_of_points(); k++) {
         lookup[i].vector = j;
         lookup[i].index = k;
         i++;
      }
   }
}

Genotype::Genotype(Genotype &original)
{
   register unsigned int i;

#ifdef DEBUG
   (void)fprintf(logFile, "support.cc/Genotype::Genotype(Genotype &original)\n");
#endif /* DEBUG */
   number_of_genes = original.number_of_genes;
   number_of_vectors = original.number_of_vectors;
   modified = original.modified;
   if (original.rep_vector!=NULL) {
      rep_vector = new Representation*[number_of_vectors];
      lookup = new Lookup[number_of_genes];
   } else {
      rep_vector = NULL;
      lookup = NULL;
   }

   for (i=0; i<number_of_vectors; i++) {
      rep_vector[i] = original.rep_vector[i]->clone();
   }

   for (i=0; i<number_of_genes; i++) {
      lookup[i] = original.lookup[i];
   }
}

Genotype::Genotype(Genotype const &original)
{
   register unsigned int i;

#ifdef DEBUG
   (void)fprintf(logFile, "support.cc/Genotype::Genotype(Genotype const &original)\n");
#endif /* DEBUG */

   number_of_genes = original.number_of_genes;
   number_of_vectors = original.number_of_vectors;
   modified = original.modified;
   if (original.rep_vector!=NULL) {
      rep_vector = new Representation*[number_of_vectors];
      lookup = new Lookup[number_of_genes];
   } else {
      rep_vector = NULL;
      lookup = NULL;
   }

   for (i=0; i<number_of_vectors; i++) {
      rep_vector[i] = original.rep_vector[i]->clone();
   }

   for (i=0; i<number_of_genes; i++) {
      lookup[i] = original.lookup[i];
   }
}

Genotype::~Genotype(void)
{
   register unsigned int i;

#ifdef DEBUG
   (void)fprintf(logFile, "support.cc/Genotype::~Genotype(void)\n");
#endif /* DEBUG */


   if (rep_vector!=NULL) {
      for (i=0; i<number_of_vectors; i++) {
         delete rep_vector[i];
      }
      delete [] rep_vector;
      delete [] lookup;
   }
}

Genotype &Genotype::operator=(const Genotype &original)
{
   register unsigned int i;

#ifdef DEBUG
   (void)fprintf(logFile, "\nsupport.cc/Genotype &Genotype::operator=(const Genotype &original): this==original is %d\n\n", this==&original);
#endif /* DEBUG */

   if (this==&original){ //Prevent self assignment
      return *this;
   }

   if (rep_vector!=NULL) {
      for (i=0; i<number_of_vectors; i++) {
         delete rep_vector[i];
      }
      delete [] rep_vector;
      delete [] lookup;
   }

   number_of_vectors = original.number_of_vectors;
   number_of_genes = original.number_of_genes;
//   modified = original.modified;
   modified = 1;

   if (original.rep_vector!=NULL) {
      rep_vector = new Representation *[number_of_vectors];
      lookup = new Lookup[number_of_genes];
   } else {
      rep_vector = NULL;
      lookup = NULL;
   }

   for (i=0; i<number_of_vectors; i++) {
      rep_vector[i] = original.rep_vector[i]->clone();
   }

   for (i=0; i<number_of_genes; i++) {
      lookup[i] = original.lookup[i];
   }
   return(*this);
}


void Genotype::write(Element value, int gene_number)
{

#ifdef DEBUG
   (void)fprintf(logFile, "support.cc/void Genotype::write(Element value, int gene_number=%d)\n", gene_number);
#endif /* DEBUG */

   modified = 1;
   rep_vector[lookup[gene_number].vector]->write(value, lookup[gene_number].index);
}

void Genotype::write(unsigned char value, int gene_number)
{

#ifdef DEBUG
   (void)fprintf(logFile, "support.cc/void Genotype::write(unsigned char value=%c, int gene_number=%d)\n", value, gene_number);
#endif /* DEBUG */

   modified = 1;
   rep_vector[lookup[gene_number].vector]->write(value, lookup[gene_number].index);
}

void Genotype::write(FourByteLong value, int gene_number)
{

#ifdef DEBUG
   (void)fprintf(logFile, "support.cc/void Genotype::write(FourByteLong value=%ld, int gene_number=%d)\n", value, gene_number);
#endif /* DEBUG */

   modified = 1;
   rep_vector[lookup[gene_number].vector]->write(value, lookup[gene_number].index);
}

void Genotype::write(double value, int gene_number)
{

#ifdef DEBUG
   (void)fprintf(logFile, "support.cc/void Genotype::write(double value=%lf, int gene_number=%d)\n", value, gene_number);
#endif /* DEBUG */

   modified = 1;
   rep_vector[lookup[gene_number].vector]->write(value, lookup[gene_number].index);
}

void Genotype::write(const Representation &value, int gene_number)
{

#ifdef DEBUG
   (void)fprintf(logFile, "support.cc/void Genotype::write(const Representation &value, int gene_number=%d)\n", gene_number);
#endif /* DEBUG */

   modified = 1;
//   *rep_vector[lookup[gene_number].vector] = value;
   *(rep_vector[gene_number]) = value;
}

Quat Genotype::readQuat()
{
    Quat q;
    q.x = gread(3).real;
    q.y = gread(4).real;
    q.z = gread(5).real;
    q.w = gread(6).real;
    // q = convertQuatToRot( q );
    return q;
}

void Genotype::writeQuat( Quat q )
{
    write( q.x, 3 );
    write( q.y, 4 );
    write( q.z, 5 );
    write( q.w, 6 );
}

Quat Phenotype::readQuat()
{
    Quat q;
    q.x = gread(3).real;
    q.y = gread(4).real;
    q.z = gread(5).real;
    q.w = gread(6).real;
    // q = convertQuatToRot( q );
    return q;
}

void Phenotype::writeQuat( Quat q )
{
    write( q.x, 3 );
    write( q.y, 4 );
    write( q.z, 5 );
    write( q.w, 6 );
}

//  Maybe we should evaluate the Phenotype?
Phenotype::Phenotype(unsigned int init_number_of_dimensions, Representation **init_value_vector)
: number_of_dimensions(init_number_of_dimensions), value_vector(init_value_vector), 
  value(0.0), evalflag(0)
{
   register unsigned int i, j, k;

#ifdef DEBUG
   (void)fprintf(logFile, "support.cc/Phenotype::Phenotype(unsigned int init_number_of_dimensions=%d, Representation **init_value_vector)\n", init_number_of_dimensions);
#endif /* DEBUG */


   number_of_points = 0;
   for (i=0; i<number_of_dimensions; i++) {
      number_of_points += value_vector[i]->number_of_points();
   }

   i = 0;
   lookup = new Lookup[number_of_points];
   for (j=0; j<number_of_dimensions; j++) {
      for (k=0; k<value_vector[j]->number_of_points(); k++) {
         assert ( i < number_of_points ); // mp!
         lookup[i].vector = j;
         lookup[i].index = k;
         i++;
      }
   }
}

Phenotype::Phenotype(const Phenotype &original)
{
   register unsigned int i;

#ifdef DEBUG
   (void)fprintf(logFile, "support.cc/Phenotype::Phenotype(const Phenotype &original)\n");
#endif /* DEBUG */


   number_of_dimensions = original.number_of_dimensions;
   number_of_points = original.number_of_points;
   evalflag = original.evalflag;
   value = original.value;

   if (original.value_vector!=NULL) {
      value_vector = new Representation *[number_of_dimensions];
      lookup = new Lookup[number_of_points];
   } else {
      value_vector = NULL;
      lookup = NULL;
   }

   for (i=0; i<number_of_dimensions; i++) {
      value_vector[i] = original.value_vector[i]->clone();
   }
   
   for (i=0; i<number_of_points; i++) {
      lookup[i] = original.lookup[i];
   }
}

Phenotype &Phenotype::operator=(const Phenotype &original)
{
   register unsigned int i;

#ifdef DEBUG
   (void)fprintf(logFile, "\nsupport.cc/Phenotype &Phenotype::operator=(const Phenotype &original): this==original is %d\n\n", this==&original);
#endif /* DEBUG */

   if (this==&original){ //Prevent self assignment
      return *this;
   }

   //  Do the destructors get called on each element of value_vector?
   if (value_vector!=NULL) {
      for (i=0; i<number_of_dimensions; i++) {
         delete value_vector[i];
      }
      delete [] value_vector;
      delete [] lookup;
   }

   number_of_dimensions = original.number_of_dimensions;
   number_of_points = original.number_of_points;
   value = original.value;
   evalflag = original.evalflag;

   if (original.value_vector!=NULL) {
      value_vector = new Representation *[number_of_dimensions];
      lookup = new Lookup[number_of_points];
   } else {
      value_vector = NULL;
      lookup = NULL;
   }

   for (i=0; i<number_of_dimensions; i++) {
      value_vector[i] = original.value_vector[i]->clone();
   }

   for (i=0; i<number_of_points; i++) {
      lookup[i] = original.lookup[i];
   }

   return(*this);
}

Phenotype::~Phenotype(void)
{
   register unsigned int i;

#ifdef DEBUG
   (void)fprintf(logFile, "support.cc/Phenotype::~Phenotype(void)\n");
#endif /* DEBUG */


   if (value_vector!=NULL) {
      for (i=0; i<number_of_dimensions; i++) {  
         delete value_vector[i];
      }
      delete [] value_vector;
      delete [] lookup;
   }
}

void Phenotype::write(Element value, int gene_number)
{

#ifdef DEBUG
   (void)fprintf(logFile, "support.cc/void Phenotype::write(Element value, int gene_number=%d)\n", gene_number);
#endif /* DEBUG */

   evalflag = 0;
   value_vector[lookup[gene_number].vector]->write(value, lookup[gene_number].index);
}

void Phenotype::write(unsigned char value, int gene_number)
{

#ifdef DEBUG
   (void)fprintf(logFile, "support.cc/void Phenotype::write(unsigned char value=%c, int gene_number=%d)\n", value, gene_number);
#endif /* DEBUG */

   evalflag = 0;
   value_vector[lookup[gene_number].vector]->write(value, lookup[gene_number].index);
}

void Phenotype::write(FourByteLong value, int gene_number)
{

#ifdef DEBUG
   (void)fprintf(logFile, "support.cc/void Phenotype::write(FourByteLong value=%ld, int gene_number=%d)\n", value, gene_number);
#endif /* DEBUG */

   evalflag = 0;
   value_vector[lookup[gene_number].vector]->write(value, lookup[gene_number].index);
}

void Phenotype::write(double value, int gene_number)
{

#ifdef DEBUG
   (void)fprintf(logFile, "support.cc/void Phenotype::write(double value=%lf, int gene_number=%d)\n", value, gene_number);
#endif /* DEBUG */

   evalflag = 0;
   value_vector[lookup[gene_number].vector]->write(value, lookup[gene_number].index);
}

void Phenotype::write(const Representation &value, int gene_number)
{

#ifdef DEBUG
   (void)fprintf(logFile, "support.cc/void Phenotype::write(const Representation &value, int gene_number=%d)\n", gene_number);
#endif /* DEBUG */

   evalflag = 0;
//   *(value_vector[lookup[gene_number].vector]) = value;
   *(value_vector[gene_number]) = value;
}

double Phenotype::evaluate(EvalMode mode)
{

#ifdef DEBUG
   (void)fprintf(logFile, "support.cc/double Phenotype::evaluate(EvalMode mode)\n");
#endif /* DEBUG */

   switch(mode)
   {
      case Always_Eval:
         value = ::evaluate(value_vector);
         evalflag = 1;
         break;
      case Always_Eval_Nonbond:
         value = ::evaluate(value_vector, 1); // term=1 as total non-bonded energy, i.e. vdW+Hb+desolv
         evalflag = 0;
         break;
      case Always_Eval_Elec:
         value = ::evaluate(value_vector, 2); // term=2 as total electrostatic energy
         evalflag = 0;
         break;
      case Normal_Eval:
         if (!evalflag) {
            value = ::evaluate(value_vector);
            evalflag = 1;
         }
         break;
      case Reset:
         evalflag = 0;
         break;
      default:
         (void)fprintf(logFile, "Unknown Evaluation Mode!\n");
         break;
   }

   return(value);
}

State Phenotype::make_state(int ntor)
{
   State retval;

#ifdef DEBUG
   (void)fprintf(logFile, "support.cc/State Phenotype::make_state(int ntor=%d)\n", ntor);
#endif /* DEBUG */


   retval.ntor = ntor;
   make_state_from_rep(value_vector, &retval);
   return(retval);
}

Individual &Population::operator[](int ind_num)
{

#ifdef DEBUG
   (void)fprintf(logFile, "support.cc/Individual &Population::operator[](int ind_num=%d)\n", ind_num);
#endif /* DEBUG */

   if ((ind_num<0)||(ind_num>=size)) {
      (void)fprintf(logFile, "ERROR: support.cc/Trying to access %d, an out of bounds individual! (0<i<%d)\n", ind_num, size);
      return(heap[0]);
   } else {
      return(heap[ind_num]);
   }
}

State Individual::state(int ntor)
{

#ifdef DEBUG
   (void)fprintf(logFile, "support.cc/State Individual::state(int ntor=%d)\n", ntor);
#endif /* DEBUG */

   return(phenotyp.make_state(ntor));
}

void Individual::getMol(Molecule *returnedMol)
{
// Converts phenotype to mol's state and returns this individual's mol data.

    State molState;
    Molecule molcopy;

    molState = phenotyp.make_state(mol->S.ntor);
    molcopy = copyStateToMolecule(&molState, mol);
    returnedMol = &molcopy;
}

void Individual::printIndividualsState(FILE *filePtr, int ntor, int detail) 
{
#ifdef DEBUG
   (void)fprintf(logFile, "support.cc/void Individual::printIndividualsState(FILE *filePtr, int ntor=%d, int detaiil=%d)\n", ntor, detail);
#endif /* DEBUG */

    printState( filePtr, state(ntor), detail ); 
    fprintf( filePtr, "\n" );
}

void Individual::incrementAge(void)
{
    ++age;
}

Population &Population::operator=(const Population &original)
{
   register int i;

#ifdef DEBUG
   (void)fprintf(logFile, "\nsupport.cc/Population &Population::operator=(const Population &original):this=original is %d\n\n", this==&original);
#endif /* DEBUG */

   if (this==&original){ //Prevent self assignement
      return *this;
   }

   if (heap!=NULL) {
      delete [] heap;
   }

   size = original.size;
   heap = new Individual[size];
   lhb = original.lhb;
   for (i=0; i<size; i++) {
      heap[i] = original.heap[i];
   }

   return(*this);
}
