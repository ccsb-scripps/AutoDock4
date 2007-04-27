/*

 $Id: gs.cc,v 1.19 2007/04/27 06:01:48 garrett Exp $

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


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/********************************************************************
     These are the methods for Global_Search and its derivations.

                                rsh 9/95
*********************************************************************/

#include <stdio.h>
#include <sys/types.h>      /*time_t time(time_t *tloc); */
#include <time.h>           /*time_t time(time_t *tloc); */
#include <sys/times.h>

#include <math.h>
#include "gs.h"
#include "ranlib.h"
#include "eval.h"
#include "rep.h"
#include "rep_constants.h"
#include "assert.h"

#ifdef sgi
    #include <ieeefp.h>
#endif
#ifdef sun
    #include <ieeefp.h>
#endif

#include "constants.h"
#include "autocomm.h"
#include "timesyshms.h"
#include "writePDBQT.h"


extern FILE *logFile;
extern class Eval evaluate;
extern int sel_prop_count;//debug
extern int global_ntor;//debug
extern int debug;//debug


double worst_in_window(double *window, int size)
{
   register int i;
   double worst;

   worst = window[0];

#ifdef DEBUG2
   (void)fprintf(logFile, "gs.cc/double worst_in_window(double *window, int size)_________________________\n");//debug
#endif

   for (i=1; i<size; i++) {

#ifdef DEBUG2
      (void)fprintf(logFile, "gs.cc/window[%d]= %.3f\tworst= %.3f\n", i, window[i], worst);//debug
#endif

      if (window[i]>worst) {
         worst = window[i];

#ifdef DEBUG2
         (void)fprintf(logFile, "gs.cc/i= %d\t(window[i]>worst)\tUpdating: worst= %.3f\n", i, worst);//debug
#endif

      }
   }// for i

#ifdef DEBUG2
   (void)fprintf(logFile, "gs.cc/Returning: worst= %.3f\n\n", worst);//debug
#endif

   return(worst);
}

double avg_in_window(double *window, int size)
{
   register int i;
   double mysum = 0.0, myavg = 0.0;

#ifdef DEBUG2
   (void)fprintf(logFile, "gs.cc/avg_in_window(double *window, int size)_________________________\n");//debug
#endif
   for (i=0; i<size; i++) {
      mysum += window[i];
#ifdef DEBUG2
      (void)fprintf(logFile, "gs.cc/mysum= %.3f\twindow[%d]= %.3f\n",mysum, i, window[i]);//debug
#endif
   }
   myavg = mysum / size;
#ifdef DEBUG2
   (void)fprintf(logFile, "gs.cc/Returning: myavg= %.3f\n\n",myavg);//debug
#endif

   return(myavg);
}

//  Also set avg
double Genetic_Algorithm::worst_this_generation(Population &pop)
{
   register unsigned int i;
   double worstval, avgval;

#ifdef DEBUG2
   (void)fprintf(logFile, "gs.cc/worst_this_generation(Population &pop)_________________________\n");
#endif

#ifdef DEBUG
   (void)fprintf(logFile, "gs.cc/double Genetic_Algorithm::worst_this_generation(Population &pop)\n");
#endif /* DEBUG */

   avgval = worstval = pop[0].value(Normal_Eval);
#ifdef DEBUG2
   (void)fprintf(logFile, "gs.cc/avgval= %.3f\tworstval= %.3f\n", avgval, worstval);
#endif
   for (i=1; i<pop.num_individuals(); i++) {
      avgval += pop[i].value(Normal_Eval);
#ifdef DEBUG2
      (void)fprintf(logFile, "gs.cc/avgval= %.3f\tpop[%d].value(Normal_Eval)= %.3f\n", avgval, i, pop[i].value(Normal_Eval));
#endif
      if (pop[i].value(Normal_Eval)>worstval) {
         worstval = pop[i].value(Normal_Eval);
#ifdef DEBUG2
         (void)fprintf(logFile, "gs.cc/(pop[i].value(Normal_Eval)>worstval): Updating: worstval= %.3f\n", worstval);
#endif
      }
   }

   avg = avgval/pop.num_individuals();
#ifdef DEBUG2
   (void)fprintf(logFile, "gs.cc/Returning: avg= %.3f, worstval= %.3f\n\n", avg, worstval);
#endif
   return(worstval);
}

//  This could be made inline

Genetic_Algorithm::Genetic_Algorithm( EvalMode init_e_mode, 
                                      Selection_Mode init_s_mode, 
                                      Xover_Mode init_c_mode,
                                      Worst_Mode init_w_mode, 
                                      int init_elitism, 
                                      Real init_c_rate, 
                                      Real init_m_rate, 
                                      int init_window_size, 
                                      unsigned int init_max_generations,
                                      unsigned int outputEveryNgens)

:  e_mode(init_e_mode),
s_mode(init_s_mode),
c_mode(init_c_mode),
w_mode(init_w_mode),
elitism(init_elitism),
c_rate(init_c_rate),
m_rate(init_m_rate),
window_size(init_window_size),
alpha(1.0),
beta(0.0),
tranStep(2.0),                         //  2 Angstroms
quatStep( DegreesToRadians( 30.0 ) ),  // 30 degrees
torsStep( DegreesToRadians( 30.0 ) ),  // 30 degrees
low(-100),
high(100),
generations(0),
max_generations(init_max_generations),
outputEveryNgens(100),
converged(0),
alloc(NULL),
mutation_table(NULL),
ordering(NULL),
m_table_size(0),
worst(0.0L),
avg(0.0L)

{
#ifdef DEBUG
   (void)fprintf(logFile, "gs.cc/Genetic_Algorithm::Genetic_Algorithm(EvalMode init_e_mode,...\n");
#endif /* DEBUG */

   worst_window = new double[window_size];
}

void Genetic_Algorithm::set_worst(Population &currentPop)
{
   double temp = 0.0;

#ifdef DEBUG
   (void)fprintf(logFile, "gs.cc/void Genetic_Algorithm::set_worst(Population &currentPop)\n");
#endif /* DEBUG */

   worst_window[generations%window_size] = worst_this_generation(currentPop);
   switch(w_mode)
   {
      //  Assume for this case that there's a window_size of 1
      case Ever:
         if (generations!=0) {
            if (temp>worst)
               worst = worst_window[0];
         } else {
            worst = worst_window[0];
         }
         break;
      case OfN:
         if (generations>=window_size) {
            worst = worst_in_window(worst_window, window_size);
         } else {
            worst = worst_in_window(worst_window, generations+1);
         }
         break;
      case AverageOfN:
         if (generations>=window_size) {
            worst = avg_in_window(worst_window, window_size);
         } else {
            worst = avg_in_window(worst_window, generations+1);
         }
         break;
      default:
         (void)fprintf(logFile,"gs.cc/Unable to set the individual with the worst fitness!\n");
   }
}

M_mode Genetic_Algorithm::m_type(RepType type)
{

#ifdef DEBUG
   (void)fprintf(logFile, "gs.cc/M_mode Genetic_Algorithm::m_type(RepType type)\n");
#endif /* DEBUG */
   switch(type)
   {
      case T_BitV:
         return(BitFlip);
      case T_RealV:
      case T_CRealV:
         return(CauchyDev);
      case T_IntV:
         return(IUniformSub);
      default:
         (void)fprintf(logFile,"gs.cc/Unrecognized Type (The allowed types are:  T_BitV, T_RealV, T_CRealV and T_IntV)!\n");
         return(ERR);
   }
}

void Genetic_Algorithm::make_table(int size, Real prob)
{
   register int i, j;
   double L = 0.0L;

#ifdef DEBUG
   (void)fprintf(logFile, "gs.cc/void Genetic_Algorithm::make_table(int size=%d, Real prob=%f)\n",size, prob);
#endif /* DEBUG */

#ifdef DEBUG_MUTATION
   (void)fprintf(logFile, "\ngs.cc/void Genetic_Algorithm::make_table(int size=%d, Real prob=%f)\n",size, prob);
#endif /* DEBUG */

   m_table_size = size;
   mutation_table = new Real[size+1];

   mutation_table[0] = pow(1-prob, size);
   mutation_table[size] = 1;

   i = 1;
   while (i<=(int)size*prob) {
      L = 0.0;
      for (j=1; j<=i; j++) {
         L += log(size+1-j) - log(j);
#ifdef DEBUG_MUTATION
         fprintf(logFile,"j= %d (size+1-j)= %d log(_)= %.4f log(j)=%.4f L= %.4f\n", j, (size+1-j), log(size+1-j), log(j), L);
#endif
      }
      L += i*log(prob) + (size-i)*log(1-prob);
#ifdef DEBUG_MUTATION
      fprintf(logFile,"j= %d (size+1-j)= %d log(_)= %.4f log(j)=%.4f L= %.4f\n", j, (size+1-j), log(size+1-j), log(j), L);
#endif

      mutation_table[i] = mutation_table[i-1]+exp(L);
#ifdef DEBUG_MUTATION
      fprintf(logFile,"mutation_table[%d] = %.3f\n", i, mutation_table[i]);
#endif
      i++;
   } // end while

#ifdef DEBUG_MUTATION
      fprintf(logFile,"L= %.3f  exp(L)= %.3f\n", L, exp(L) );
#endif
   L = exp(L);
   for (; i<size; i++) {
#ifdef DEBUG_MUTATION
      fprintf(logFile,"i= %d  L= %.3f  prob= %.3f size= %d  (size+1-i)= %d  i*(1-prob)= %.3f\n", i, L, prob, size, (size+1-i), i*(1-prob) );
#endif
      L = (L*prob*(size+1-i))/(i*(1-prob));
#ifdef DEBUG_MUTATION
      fprintf(logFile,"L= %.3f\n", L );
#endif
      mutation_table[i] = mutation_table[i-1]+L;
#ifdef DEBUG_MUTATION
      fprintf(logFile,"mutation_table[%d] = %.3f\n", i, mutation_table[i]);
#endif
   }
#ifdef DEBUG_MUTATION
      fflush(logFile);
#endif
}

int Genetic_Algorithm::check_table(Real prob)
{
   int low, high;

#ifdef DEBUG
   (void)fprintf(logFile, "gs.cc/int Genetic_Algorithm::check_table(Real prob=%f)\n",prob);
#endif /* DEBUG */

#ifdef DEBUG_MUTATION
   (void)fprintf(logFile, "\ngs.cc/int Genetic_Algorithm::check_table(Real prob=%f)\n",prob);
#endif /* DEBUG */

   low = 0; high = m_table_size;
#ifdef DEBUG_MUTATION
      fprintf(logFile,"prob= %.3f  low= %d  high= %d\n", prob, low, high );
#endif

   while (high-low>1) {
#ifdef DEBUG_MUTATION
      fprintf(logFile,"high-low= %d\n", high-low );
      fprintf(logFile,"(high+low)/2= %d  mutation_table[(high+low)/2]= %.3f  prob= %.3f\n", (high+low)/2,  mutation_table[(high+low)/2], prob );
#endif
      if (mutation_table[(high+low)/2]<prob) {
         low = (high+low)/2;
#ifdef DEBUG_MUTATION
         fprintf(logFile,"low= %d\n", low );
#endif
      } else if (mutation_table[(high+low)/2]>prob) {
         high = (high+low)/2;
#ifdef DEBUG_MUTATION
         fprintf(logFile,"high= %d\n", high );
#endif
      } else {
         high = low = (high+low)/2;
#ifdef DEBUG_MUTATION
         fprintf(logFile,"low= %d  high= %d\n", low, high );
#endif
      }
   }
   return(low);
}

void Genetic_Algorithm::initialize(unsigned int pop_size, unsigned int num_poss_mutations)
{
   register unsigned int i;

#ifdef DEBUG
   (void)fprintf(logFile, "gs.cc/void Genetic_Algorithm::initialize(unsigned int pop_size=%d, ",pop_size);
   (void)fprintf(logFile, "unsigned int num_poss_mutations=%d)\n",num_poss_mutations);
#endif /* DEBUG */

   if (alloc!=NULL) {
      delete [] alloc;
   }

   if (ordering!=NULL) {
      delete [] ordering;
   }

   if (mutation_table!=NULL) {
      delete [] mutation_table;
   }

   alloc = new Real[pop_size];

   ordering = new unsigned int[pop_size];
   for (i=0; i<pop_size; i++) {
      ordering[i] = i;
      assert(ordering[i] < pop_size);//debug
      alloc[i] = 1.0; // changed by gmm, 12-sep-1997.
   }

   make_table(pop_size*num_poss_mutations, m_rate);
}

void Genetic_Algorithm::mutate(Genotype &mutant, int gene_number)
{
   Element tempvar;

#ifdef DEBUG
   (void)fprintf(logFile, "gs.cc/void Genetic_Algorithm::mutate(Genotype &mutant, int gene_number=%d)\n",gene_number);
#endif /* DEBUG */

#ifdef DEBUG_MUTATION
   (void)fprintf(logFile, "gs.cc/void Genetic_Algorithm::mutate(Genotype &mutant, int gene_number=%d)\n",gene_number);
#endif /* DEBUG */

   switch(m_type(mutant.gtype(gene_number))) {

      case BitFlip:
         //((unsigned char *)gene)[point] = 1 - ((unsigned char *)gene)[point];
         //  Read the bit
         tempvar = mutant.gread(gene_number);
         //  Flip it
         tempvar.bit = 1 - tempvar.bit;
         //  write it
         mutant.write(tempvar, gene_number);
         break;

      case CauchyDev:
         // gene_numbers 3, 4, 5 and 6 correspond to the 
         // four components of the quaternion, (x,y,z,w)
         if ( is_rotation_index( gene_number ) ) {
            // Mutate all four comopnents of the quaternion, (x,y,z,w) simultaneously:
            // Generate a uniformly-distributed quaternion
            Quat q_change;
            q_change = uniformQuat();
#ifdef DEBUG_MUTATION
            fprintf( logFile, "q_change -- after uniformQuat\n" );
            printQuat_q( logFile, q_change );
#endif
            Quat q_current = mutant.readQuat();
#ifdef DEBUG_MUTATION
            fprintf( logFile, "q_current -- after reading mutant.gread\n" );
            printQuat_q( logFile, q_current );
#endif
            Quat q_new;
#ifdef DEBUG_MUTATION
            fprintf( logFile, "q_current -- about to call qmultiply\n" );
#endif
#ifdef DEBUG_QUAT
#ifdef DEBUG_QUAT_PRINT
            pr( logFile, "DEBUG_QUAT: mutate() -- q_current\n" );
            (void) fflush(logFile);
#endif // DEBUG_QUAT_PRINT
            //  Make sure the quaternion is suitable for 3D rotation
            assertQuatOK( q_current );
#ifdef DEBUG_QUAT_PRINT
            pr( logFile, "DEBUG_QUAT: mutate() -- q_change\n" );
            (void) fflush(logFile);
#endif // DEBUG_QUAT_PRINT
            //  Make sure the quaternion is suitable for 3D rotation
            assertQuatOK( q_change );
#endif // DEBUG_QUAT
            // Multiply the quaternions, applying the rotation to the current orientation
            qmultiply( &q_new, &q_current, &q_change );
#ifdef DEBUG_QUAT
#ifdef DEBUG_QUAT_PRINT
            pr( logFile, "DEBUG_QUAT: mutate() -- q_current\n" );
            (void) fflush(logFile);
#endif
            //  Make sure the quaternion is suitable for 3D rotation
            assertQuatOK( q_new );
#endif // DEBUG_QUAT
#ifdef DEBUG_MUTATION
            fprintf( logFile, "q_new - after qmultiply\n" );
            printQuat_q( logFile, q_new );
#endif
            mutant.writeQuat( q_new );
         } else {
             //  Read the real
             tempvar = mutant.gread(gene_number);
#ifdef DEBUG_MUTATION
             (void)fprintf(logFile, "   ---CauchyDev---\n" );
             (void)fprintf(logFile, "   Before mutating:        tempvar= %.3f\n", tempvar.real );
             (void)fprintf(logFile, "   gene_number= %d\n", gene_number );
             (void)fprintf(logFile, "   tempvar.real += scauchy2()\n" );
#endif
             //  Add deviate
             //  We never vary alpha and beta, so just use the faster "scauchy2()" function:
             tempvar.real += scauchy2();
#ifdef DEBUG_MUTATION
             (void)fprintf(logFile, "   Add Cauchy deviate:     tempvar= %.3f\n", tempvar.real );
             (void)fflush(logFile );
#endif
             //  Write it
             mutant.write(tempvar, gene_number);
         }
         break;

      case IUniformSub:
         //((int *)gene)[point] = ignuin(low, high);
         //  Generate the new integer
         tempvar.integer = ignuin(low, high);
         //  Write it
         //mutant.write((void *)(&tempvar.integer), gene_number);
         mutant.write(tempvar, gene_number);
         break;

      default:
         (void)fprintf(logFile,"gs.cc/Unrecognized mutation Mode!\n");
         break; 
   }
}

void Genetic_Algorithm::mutation(Population &pure)
{
   int num_mutations, individual, gene_number;

#ifdef DEBUG
   (void)fprintf(logFile, "gs.cc/void Genetic_Algorithm::mutation(Population &pure)\n");
#endif /* DEBUG */

#ifdef DEBUG_MUTATION
   (void)fprintf(logFile, "gs.cc/void Genetic_Algorithm::mutation(Population &pure)\n");
#endif /* DEBUG */

   num_mutations = check_table(ranf());

   //  Note we don't check to see if we mutate the same gene twice.
   //  So, effectively we are lowering the mutation rate, etc...
   //  Might want to check out Bentley's chapter on selection.
   for (; num_mutations>0; num_mutations--) {
      individual = ignlgi()%pure.num_individuals();
      gene_number = ignlgi()%pure[individual].genotyp.num_genes();
      mutate(pure[individual].genotyp, gene_number);
      pure[individual].age = 0L;
#ifdef DEBUG_MUTATION
      (void)fprintf(logFile, "num_mutations= %d\n", num_mutations );
      (void)fprintf(logFile, "   individual= %d\n", individual );
      (void)fprintf(logFile, "   gene_number= %d\n", gene_number );
#endif /* DEBUG */
   }
}

void Genetic_Algorithm::crossover(Population &original_population)
{
   register unsigned int i;
   int first_point, second_point, temp_index, temp_ordering;
   Real alpha = 0.5;
   
#ifdef DEBUG
   (void)fprintf(logFile, "gs.cc/void Genetic_Algorithm::crossover(Population &original_population)\n");
#endif /* DEBUG */

   //  Shuffle the population
   //  Permute the ordering of the population, "original_population"
   for (i=0; i<original_population.num_individuals(); i++) {
      temp_ordering = ordering[i];
      // ignlgi is GeNerate LarGe Integer and is in com.cc
      temp_index = ignlgi()%original_population.num_individuals();
      ordering[i] = ordering[temp_index];
      assert(ordering[i] < original_population.num_individuals());//debug
      ordering[temp_index] = temp_ordering;
      assert(ordering[temp_index] < original_population.num_individuals());//debug
   }
                                                                    
   //  How does Goldberg implement crossover?

   // Loop over individuals in population
   for (i=0; i<original_population.num_individuals()-1; i=i+2) {
      // The two individuals undergoing crossover are original_population[ordering[i]] and original_population[ordering[i+1]]

      if (ranf() < c_rate) {
         // Perform crossover with a probability of c_rate
         
         switch(c_mode) {
            case TwoPt:
                // First crossover point is a random integer from 0 to the number of genes minus 1
#ifdef DO_NOT_CROSSOVER_IN_QUAT
                // Make sure we do not crossover inside a quaternion...
                do {
                    first_point = ignuin(0, original_population[i].genotyp.num_genes()-1);
                } while ( is_rotation_index( first_point ) ) ;
                do {
                    second_point = first_point+ignuin(0, original_population[i].genotyp.num_genes()-first_point-1);
                } while ( is_rotation_index( second_point ) );
#else
                first_point = ignuin(0, original_population[i].genotyp.num_genes()-1);
                second_point = first_point+ignuin(0, original_population[i].genotyp.num_genes()-first_point-1);
#endif
                // Do two-point crossover, with the crossed-over offspring replacing the parents in place:
                crossover_2pt( original_population[ordering[i]].genotyp, 
                               original_population[ordering[i+1]].genotyp, 
                               first_point, 
                               second_point );
                original_population[ordering[i]].age = 0L;
                original_population[ordering[i+1]].age = 0L;
                break;
            case OnePt:
                // First crossover point is a random integer from 0 to the number of genes minus 1
                // Make sure we do not crossover inside a quaternion...
                do {
                    first_point = ignlgi()%original_population[i].genotyp.num_genes();
                } while ( is_rotation_index( first_point ) ) ;
                //  We can accomplish one point crossover by using the 2pt crossover operator
                crossover_2pt( original_population[ordering[i]].genotyp, 
                               original_population[ordering[i+1]].genotyp,
                               first_point, 
                               original_population[ordering[i]].genotyp.num_genes()-1 );
                original_population[ordering[i]].age = 0L;
                original_population[ordering[i+1]].age = 0L;
                break;
            case Uniform:
                crossover_uniform( original_population[ordering[i]].genotyp, 
                                   original_population[ordering[i+1]].genotyp,
                                   original_population[ordering[i]].genotyp.num_genes() - 1);

                break;
            case Arithmetic:
               // select the parents A and B
               // create new offspring, a and b, where
               // a = x*A + (1-x)*B, and b = (1-x)*A + x*B    -- note: x is alpha in the code
               alpha = (Real) ranf();
#ifdef DEBUG
               (void)fprintf(logFile, "gs.cc/  alpha = " FDFMT "\n", alpha);
               (void)fprintf(logFile, "gs.cc/ About to call crossover_arithmetic with original_population[%d] & [%d]\n", i, i+1);
#endif /* DEBUG */
               crossover_arithmetic( original_population[ i ].genotyp, 
                                     original_population[i+1].genotyp, 
                                     alpha );
               break;
            default:
                (void)fprintf(logFile,"gs.cc/ Unrecognized crossover mode!\n");
         }
      }
   }
}


void Genetic_Algorithm::crossover_2pt(Genotype &father, Genotype &mother, unsigned int pt1, unsigned int pt2)
{
    /*  Assumes that 0<=pt1<pt2<=number_of_pts  
     *  There are four cases to consider:-
     *  (1) the copied area is contained entirely within the gene
     *  (2) the gene is contained entirely within the copied area
     *  (3) the copied area is partially contained within the gene
     *  (4) there's no intersection between the copied area and the gene
     */
   register unsigned int i;
   Element temp;

#ifdef DEBUG
   (void)fprintf(logFile, "gs.cc/void Genetic_Algorithm::crossover_2pt(Genotype");
   (void)fprintf(logFile, "&father, Genotype &mother, unsigned int pt1, unsigned int pt2)\n");
   (void)fprintf(logFile,"gs.cc/Performing crossover from %d to %d \n", pt1,pt2);
#endif /* DEBUG */

   // loop over genes to be crossed over
   for (i=pt1; i<=pt2; i++) {
#ifdef DEBUG
      //(void)fprintf(logFile,"gs.cc/1::At pt %d   father: %.3lf   mother: %.3lf\n",
      //i, *((double *)father.gread(i)), *((double *)mother.gread(i)) );
#endif /* DEBUG */
      temp = father.gread(i);
      father.write(mother.gread(i), i);
      mother.write(temp, i);
#ifdef DEBUG
      //(void)fprintf(logFile,"gs.cc/1::At pt %d   father: %.3lf   mother: %.3lf\n",
      //i, *((double *)father.gread(i)), *((double *)mother.gread(i)) );
#endif /* DEBUG */
   }

#ifdef DEBUG_QUAT
   Quat q_father, q_mother;

#ifdef DEBUG_QUAT_PRINT
   pr( logFile, "DEBUG_QUAT: crossover_2pt()  q_father\n" );
   pr( logFile, "DEBUG_QUAT: crossover_2pt()  pt1=%d, pt2=%d\n", pt1, pt2 );
#endif // endif DEBUG_QUAT_PRINT

   //  Make sure the quaternion is suitable for 3D rotation
   q_father = father.readQuat();
#ifndef DO_NOT_CROSSOVER_IN_QUAT
   q_father = normQuat( q_father );
   father.writeQuat( q_father );
#endif
#ifdef DEBUG_QUAT_PRINT
   printQuat( logFile, q_father );
   (void) fflush(logFile);
#endif // endif DEBUG_QUAT_PRINT
   assertQuatOK( q_father );
#ifdef DEBUG_QUAT_PRINT
   pr( logFile, "DEBUG_QUAT: crossover_2pt()  q_mother\n" );
   (void) fflush(logFile);
#endif // endif DEBUG_QUAT_PRINT
   //  Make sure the quaternion is suitable for 3D rotation
   q_mother = mother.readQuat();
#ifndef DO_NOT_CROSSOVER_IN_QUAT
   q_mother = normQuat( q_mother );
   mother.writeQuat( q_mother );
#endif
#ifdef DEBUG_QUAT_PRINT
   printQuat( logFile, q_mother );
   (void) fflush(logFile);
#endif // endif DEBUG_QUAT_PRINT
   assertQuatOK( q_mother );
#endif // endif DEBUG_QUAT
}


void Genetic_Algorithm::crossover_uniform(Genotype &father, Genotype &mother, unsigned int num_genes)
{
    register unsigned int i;
    Element temp;

#ifdef DEBUG
    (void)fprintf(logFile, "gs.cc/void Genetic_Algorithm::crossover_uniform(Genotype");
    (void)fprintf(logFile, "&father, Genotype &mother, unsigned int num_genes)\n");
#endif /* DEBUG */

    for (i=0; i<num_genes; i++) {
        // Choose either father's or mother's gene, with a 50/50 probability
        if (ranf() > 0.5) {
            temp = father.gread(i);
            father.write(mother.gread(i), i);
            mother.write(temp, i);
        }
    }
}

void Genetic_Algorithm::crossover_arithmetic(Genotype &A, Genotype &B, Real alpha)
{
   register unsigned int i;
   Element temp_A, temp_B;
   Real one_minus_alpha;

   one_minus_alpha = 1.0 - alpha;

#ifdef DEBUG
   (void)fprintf(logFile, "gs.cc/void Genetic_Algorithm::crossover_arithmetic(Genotype");
   (void)fprintf(logFile, "&A, Genotype &B, Real alpha)\n");
   (void)fprintf(logFile, "gs.cc/void Genetic_Algorithm::crossover_arithmetic");
   (void)fprintf(logFile, "/Trying to perform arithmetic crossover using alpha = " FDFMT "\n", alpha);
   (void)fflush(logFile);
#endif /* DEBUG */

   // loop over genes to be crossed over
   for (i=0; i<A.num_genes(); i++) {
#ifdef DEBUG
       (void)fprintf(logFile, "gs.cc/void Genetic_Algorithm::crossover_arithmetic");
       (void)fprintf(logFile, "/looping over genes to be crossed over, i = %d\n", i);
       (void)fflush(logFile);
#endif /* DEBUG */
      temp_A = A.gread(i);
      temp_B = B.gread(i);
#ifdef DEBUG
       (void)fprintf(logFile, "gs.cc/void Genetic_Algorithm::crossover_arithmetic");
       (void)fprintf(logFile, "/temp_A = %.3f  &  temp_B = %.3f\n", temp_A.real, temp_B.real);
       (void)fflush(logFile);
#endif
      // a = alpha*A + (1-alpha)*B
      // b = (1-alpha)*A + alpha*B
      A.write( (alpha * temp_A.real  +  one_minus_alpha * temp_B.real), i);
      B.write( (one_minus_alpha * temp_A.real  +  alpha * temp_B.real), i);
#ifdef DEBUG
       (void)fprintf(logFile, "gs.cc/void Genetic_Algorithm::crossover_arithmetic");
       (void)fprintf(logFile, "/A = %.3f  &  B = %.3f\n", A.gread(i).real, B.gread(i).real);
       (void)fflush(logFile);
#endif
   }
}

/* * */

/*
 * Proportional Selection
 *
 *
 *  We want to perform minimization on a function.  To do so, we
 *      take fitness values in a given generation and normalize them s.t.
 *      they are non-negative, and such that smaller f's are given a higher
 *      weighting.
 *  If f in [a,b], then f' in [A,B] s.t. f(a) => f'(B) and f(b) => f'(A) is
 *
 *      f' = (B-A) * (1 - (f-a)/(b-a)) + A
 *
 *         = (B-A) * (b-f)/(b-a) + A
 *
 *  Note that it suffices to map f into [0,1], giving
 *
 *      f' = (b-f)/(b-a)
 *
 *  Now the expected number of samples generated from a given point is
 *
 *      N * f'_i / \sum f'_i  =  N * ((b-f_i)/(b-a)) / \sum ((b-f_i)/(b-a))
 *
 *                            =  N * (b-f_i) / (N*b - \sum f_i)
 *
 *  Within a given generation, let 'b' be the maximal (worst) individual and 
 *      let 'a' be the minimal individual. 
 *
 *  Note:
 *      (1) This calculation is invariant to the value of B, but _not_ 
 *              invariant to the value of A.
 *      (2) This selection strategy works fine for functions bounded above,
 *              but it won't necessarily work for functions which are unbounded
 *              above.
 *
 *  The 'b' parameter is represented as 'Worst' and is selected by the 
 *    scale_fitness method.  If a value is greater than Worst, it is given a 
 *    value of zero for it's expectation.
 */

void Genetic_Algorithm::selection_proportional(Population &original_population, Individual *new_pop)
{
   register unsigned int i=0, start_index = 0;
   int temp_ordering, temp_index;
#ifdef DEBUG2
   Real debug_ranf;
   int allzero = 1;//debug
   Molecule *individualMol;//debug
#endif

#ifdef CHECK_ISNAN
   int allEnergiesEqual = 1;
   double diffwa = 0.0, invdiffwa = 0.0, firstEnergy = 0.0;
#endif

#ifdef DEBUG
   (void)fprintf(logFile, "gs.cc/void Genetic_Algorithm::");
   (void)fprintf(logFile, "selection_proportional(Population &original_population, Individual *new_pop)\n");
#endif /* DEBUG */

#ifdef DEBUG2
   (void)fprintf(logFile, "gs.cc/At the start of sel_prop:  sel_prop_count= %d, start_index= %d\n\n",sel_prop_count, start_index); //debug

   original_population.printPopulationAsStates(logFile, original_population.num_individuals(), global_ntor);//debug
#endif

#ifdef CHECK_ISNAN

   /*
    * This is the new code to check for the NaN case that could
    * arise; this will preempt any fatal errors.
    */

   // Calculate expected number of children for each individual

   assert(finite(worst));
   assert(finite(avg));
   assert(!ISNAN(worst));
   assert(!ISNAN(avg));

   diffwa = worst - avg;

   assert(finite(diffwa));
   assert(!ISNAN(diffwa));

   if (diffwa != 0.0) { // added by gmm, 4-JUN-1997

      invdiffwa = 1.0 / diffwa;

      if (ISNAN(invdiffwa)) {
          (void)fprintf(logFile,"WARNING!  While doing proportional selection, not-a-number was detected (NaN).\n");
          (void)fprintf(logFile,"All members of the population will be arbitrarily allocated 1 child each.\n\n");
          for (i=0;  i < original_population.num_individuals();  i++) {
             alloc[i] = 1.0;  // arbitrary
          }
          // Assume run has converged:
          converged = 1;
          (void)fprintf(stderr,"WARNING!  The population appears to have converged, so this run will shortly terminate.\n");
      } else {
          assert(finite(invdiffwa));
          assert(finite(worst));
          assert(finite(original_population.num_individuals()));
          assert(!ISNAN(invdiffwa));
          assert(!ISNAN(worst));
          assert(!ISNAN(original_population.num_individuals()));

          for (i=0;  i < original_population.num_individuals();  i++) {
             alloc[i] = (worst - original_population[i].value(e_mode)) * invdiffwa;

             assert(finite(original_population[i].value(e_mode)));
             assert(finite(alloc[i]));
             assert(!ISNAN(original_population[i].value(e_mode)));
             assert(!ISNAN(alloc[i]));
          }// for i
      }// endif (ISNAN(invdiffwa))
   } else {
      // diffwa = 0.0,  so worst = avg
      // This implies the population may have converged.
      converged = 1;

      // Write warning message to both stderr and logFile...
      (void)fprintf(stderr,"WARNING!  While doing proportional selection, worst (%6.2le) = avg (%6.2le).\n", worst, avg);
      (void)fprintf(stderr,"          This would cause a division-by-zero error.\n");
      (void)fprintf(stderr,"          All members of the population will be arbitrarily allocated 1 child each.\n");
      (void)fprintf(stderr,"WARNING!  The population appears to have converged, so this run will shortly terminate.\n\n");

      (void)fprintf(logFile,"WARNING!  While doing proportional selection, worst (%6.2le) = avg (%6.2le).\n", worst, avg);
      (void)fprintf(logFile,"          This would cause a division-by-zero error.\n");
      (void)fprintf(logFile,"          All members of the population will be arbitrarily allocated 1 child each.\n");
      (void)fprintf(logFile,"WARNING!  The population appears to have converged, so this run will shortly terminate.\n\n");

      alloc[0] = 1.0; // Added by gmm, 2-APR-1997
      firstEnergy = original_population[0].value(e_mode);
      allEnergiesEqual = 1;
      for (i=1;  i < original_population.num_individuals();  i++) {
         alloc[i] = 1.0; // Added by gmm, 2-APR-1997
         allEnergiesEqual = allEnergiesEqual && (firstEnergy == original_population[i].value(e_mode));
      }
      if (allEnergiesEqual) {
          (void)fprintf(logFile,"          All individuals in the population have the same fitness (%6.2le)\n", firstEnergy);
      } else {
          (void)fprintf(logFile,"          Here are the fitness values of the population:\n\n");
          for (i=0;  i < original_population.num_individuals();  i++) {
              (void)fprintf(logFile,"%3d = %6.2le,  ", i+1, original_population[i].value(e_mode) );
          }
      }
      (void)fprintf(logFile,"\n");
   } // diffwa = 0.0,  so worst = avg

   /*
    * Commented out because this writes too many
    * errors out -- (worst-avg) is constant inside loop,
    * so can be brought outside for-loop.
    * for (i=0; i<original_population.num_individuals(); i++) {
    * if ((worst - avg) != 0.0) { // added by gmm, 4-JUN-1997
    * //  In function minimization, the max energy is the worst
    * alloc[i] = (worst - original_population[i].value(e_mode))/(worst - avg);
    * } else {
    * (void)fprintf(logFile,"gs.cc/WARNING!  While doing proportional selection,
    * worst (%6.2le) and avg (%6.2le) were found equal, which would cause a 
    * division-by-zero error; this was just prevented, for individual %d\n", 
    * worst, avg, i); // Added by gmm, 4-JUN-1997
    * alloc[i] = 1.0; // Added by gmm, 2-APR-1997
    * }
    * }
    */

#else 

   /*
    * This is how the code used to be, before the worst=avg problem
    * was observed.
    */

   // Calculate expected number of children for each individual
   for (i=0;  i < original_population.num_individuals();  i++)
   {
      //  In our case of function minimization, the max individual is the worst
      alloc[i] = (worst - original_population[i].value(e_mode))/(worst - avg);
   }

#endif /* not CHECK_ISNAN */

#ifdef DEBUG2
   allzero = 1; //debug
   unsigned int J;//debug
   (void)fprintf(logFile, "gs.cc: checking that all alloc[] variables are not all zero...\n"); //debug
   for (J=0;  J < original_population.num_individuals();  J++) {//debug
       allzero = allzero & (alloc[J] == (Real)0.0);//debug
   }//debug
   if (allzero) {//debug
       (void)fprintf(logFile, "gs.cc:  W A R N I N G !  all alloc variables are zero!\n"); //debug
   }//debug
#endif

   //  Permute the individuals

#ifdef DEBUG2
   (void)fprintf(logFile, "gs.cc:  Permuting beginning: original_population.num_individuals()= %d\n", original_population.num_individuals()); //debug
#endif

   for (i=0;  i < original_population.num_individuals();  i++) {
      temp_ordering = ordering[i];

      assert(ordering[i] < original_population.num_individuals());//debug

      temp_index = ignlgi()%(original_population.num_individuals());

      assert(ordering[temp_index] < original_population.num_individuals());//debug

      ordering[i] = ordering[temp_index];
      ordering[temp_index] = temp_ordering;

#ifdef DEBUG2
      //(void)fprintf(logFile, "gs.cc: permuting in sel_prop: ordering check: ordering[i=%d]= %d, ordering[temp_index=%d]= %d\n", ordering[i],i, ordering[temp_index], temp_index);//debug
#endif
   }// endfor i

   //  We might get some savings here if we sorted the individuals before calling
   //  this routine
   for (i=0; (i < original_population.num_individuals()) && (start_index < original_population.num_individuals()); i++) {
      for (; (alloc[i] >= 1.0) && (start_index < original_population.num_individuals());  alloc[i]-= 1.0) {
         new_pop[start_index] = original_population[i];
         //new_pop[start_index].incrementAge();
         ++start_index;
      }
   }

#ifdef DEBUG2
   (void)fprintf(stderr, "gs.cc/void Genetic_Algorithm::"); //debug
   (void)fprintf(stderr, "selection_proportional(Population &original_population, Individual *new_pop)\n"); //debug
#endif

   i = 0;

#ifdef DEBUG2
   int count = 0;//debug
   (void)fprintf(stderr, "gs.cc/beginning \"while(start_index < original_population.num_individuals()) {\" loop\n"); //debug
#endif

   // ??? start_index = 0; // gmm, 1998-07-13 ???

   while (start_index < original_population.num_individuals()) {
#ifdef DEBUG2
      (void)fprintf(stderr, "gs.cc:596/inside \"while(start_index(=%d) < original_population.num_individuals()(=%d)) \" loop:  count= %d\n", start_index, original_population.num_individuals(), ++count); //debug
#endif
      assert(ordering[i] < original_population.num_individuals());//debug
#ifdef DEBUG2
      debug_ranf = ranf();
      (void)fprintf(stderr, "gs.cc:599/inside debug_ranf= %.3f, alloc[ordering[i]]= %.3e, ordering[i]= %d,  i= %d\n", debug_ranf, alloc[ordering[i]], ordering[i], i); // debug
      if (debug_ranf < alloc[ordering[i]]) {
#else
      if (ranf() < alloc[ordering[i]]) {                        // non-debug
#endif //  not DEBUG2
#ifdef DEBUG2
         (void)fprintf(stderr, "gs.cc:603/inside (debug_ranf < alloc[ordering[i]]) is true!\n"); //debug
         (void)fprintf(stderr, "gs.cc:604/inside about to increment start_index in:  \"new_pop[start_index++] = original_population[ordering[i]];\"; right now, start_index= %d\n", start_index); //debug
#endif
         new_pop[start_index] = original_population[ordering[i]];
         //new_pop[start_index].incrementAge();
         start_index++;
#ifdef DEBUG2
         (void)fprintf(stderr, "gs.cc:605/inside just incremented start_index, now start_index= %d\n", start_index); //debug
#endif
      }// endif (ranf() < alloc[ordering[i]])

#ifdef DEBUG2
      (void)fprintf(stderr, "gs.cc:608/inside i= %d, original_population.num_individuals()= %d\n", i, original_population.num_individuals()); //debug
      (void)fprintf(stderr, "gs.cc:609/inside about to \"i = (i+1)%%original_population.num_individuals();\"\n"); //debug
#endif
      i = (i+1)%original_population.num_individuals();
#ifdef DEBUG2
      (void)fprintf(stderr, "gs.cc:611/inside just done \"i = (i+1)%%original_population.num_individuals();\"\n"); //debug
      (void)fprintf(stderr, "gs.cc:612/inside i= %d  _____________________________________________________\n\n", i); //debug

       allzero = 1;//debug
       for (J=0;  J < original_population.num_individuals();  J++) {//debug
           allzero = allzero & (alloc[J] == (Real)0.0);//debug
       }//debug
       if (allzero) {//debug
           (void)fprintf(logFile, "gs.cc:  W A R N I N G !  all alloc variables are zero!\n"); //debug
       }//debug
#endif

   }// endwhile (start_index < original_population.num_individuals())

#ifdef DEBUG2
  (void)fprintf(stderr, "gs.cc/finished \"while(start_index < original_population.num_individuals()) \" loop\n"); //debug
#endif

}

/* Probabilistic Binary Tournament Selection
 *
 * This type of tournament selection was described in Goldberg and Deb (1991),
 * and performs the same type of selection as is seen in linear rank
 * selection.  Two individuals are chosen at random, and the better individual
 * is selected with probability P, 0.5 <= P <= 1.  If we let 2P = C, then
 * we see that the expected number of samples generated by rank r_i is
 *
 *      N * [2P - 2(2P-1) r_i]            , 1 >= P >= 0.5
 *
 * We can again parameterize this ranking with K, the relative probability
 * between the best and worst individual.  Since 2P = C,
 * P = K/(1+K).
 */
void Genetic_Algorithm::selection_tournament(Population &original_population, Individual *new_pop)
{
   register unsigned int i = 0, start_index = 0;
   int temp_ordering, temp_index;

#ifdef DEBUG
   (void)fprintf(logFile, "gs.cc/void Genetic_Algorithm::");
   (void)fprintf(logFile, "selection_tournament(Population &original_population, Individual *new_pop)\n");
#endif /* DEBUG */

   original_population.msort(original_population.num_individuals());
   for (i=0; i<original_population.num_individuals(); i++) {
      alloc[i] = original_population.num_individuals()*(2*tournament_prob - i*(4*tournament_prob - 2));
   }

   for (i=0;  (i < original_population.num_individuals()) && (start_index < original_population.num_individuals());  i++) {
      for (; (alloc[i] >= 1.0) && (start_index < original_population.num_individuals());  alloc[i] -= 1.0) {
         new_pop[start_index++] = original_population[i];
      }
   }

   for (i=0; i < original_population.num_individuals(); i++) {
      temp_ordering = ordering[i];
      temp_index = ignlgi()%original_population.num_individuals();
      ordering[i] = ordering[temp_index];
      ordering[temp_index] = temp_ordering;
   }

   i = 0;
   while (start_index < original_population.num_individuals()) {
      if (ranf() < alloc[ordering[i]]) {
         new_pop[start_index++] = original_population[ordering[i]];
      }
      i = (i+1)%original_population.num_individuals();
   }
}

Individual *Genetic_Algorithm::selection(Population &solutions)
{
   Individual *next_generation;

#ifdef DEBUG
   (void)fprintf(logFile, "gs.cc/Individual *Genetic_Algorithm::selection(Population &solutions)\n");
#endif /* DEBUG */

   next_generation = new Individual[solutions.num_individuals()];
   
   set_worst(solutions);
   switch(s_mode)
   {
      case Proportional:
         selection_proportional(solutions, next_generation);
         break;
      case Tournament:
         selection_tournament(solutions, next_generation);
         break;
      case Boltzmann:
         (void)fprintf(logFile,"gs.cc/Unimplemented Selection Method - using proportional\n");
         selection_proportional(solutions, next_generation);
         break;
      default:
         (void)fprintf(logFile,"gs.cc/Unknown Selection Mode!\n");
   }

   return(next_generation);
}

//  For right now global search is taken to be a GA
//
//  This is where the action is... SEARCH!
//
int Genetic_Algorithm::search(Population &solutions)
{
   register unsigned int i;
   unsigned int oldest = 0, oldestIndividual = 0, fittestIndividual = 0;
   double fittest = BIG;

   struct tms tms_genStart;
   struct tms tms_genEnd;

   Clock genStart;
   Clock genEnd;

#ifdef DEBUG /* DEBUG { */
   (void)fprintf(logFile, "gs.cc/int Genetic_Algorithm::search(Population &solutions)\n");
#endif /* } DEBUG */

   genStart = times( &tms_genStart );

#ifdef DEBUG3 /* DEBUG3 { */
   (void)fprintf(logFile,"About to perform Mapping on the solutions.\n");
   for (i=0; i<solutions.num_individuals(); i++) {
       (void)fprintf(logFile,"%ld ", solutions[i].age);
   }
   (void)fprintf(logFile,"\n");
#endif /* } DEBUG3 */

   //
   // Map from genotype to phenotype
   //
   for (i=0; i<solutions.num_individuals(); i++) {
      solutions[i].mapping();
   }
   
#ifdef DEBUG3 /* DEBUG3 { */
   (void)fprintf(logFile,"About to perform Selection on the solutions.\n");
   for (i=0; i<solutions.num_individuals(); i++) {
       (void)fprintf(logFile,"%ld ", solutions[i].age);
   }
   (void)fprintf(logFile,"\n");
#endif /* } DEBUG3 */

   //
   // Perform selection
   //
   Population newPop(solutions.num_individuals(), selection(solutions));

#ifdef DEBUG3 /* DEBUG3 { */
   (void)fprintf(logFile,"About to perform Crossover on the population, newPop.\n");
   for (i=0; i<solutions.num_individuals(); i++) {
       (void)fprintf(logFile,"%ld ", newPop[i].age);
   }
   (void)fprintf(logFile,"\n");
#endif /* } DEBUG3 */

   //
   // Perform crossover
   // 
   crossover(newPop);

#ifdef DEBUG3 /* DEBUG3 { */
   (void)fprintf(logFile,"About to perform mutation on the population, newPop.\n");
   for (i=0; i<solutions.num_individuals(); i++) {
       (void)fprintf(logFile,"%ld ", newPop[i].age);
   }
   (void)fprintf(logFile,"\n");
#endif /* } DEBUG3 */

   //
   // Perform mutation
   // 
   mutation(newPop);

#ifdef DEBUG3 /* DEBUG3 { */
   (void)fprintf(logFile,"About to perform elitism, newPop.\n");
   for (i=0; i<solutions.num_individuals(); i++) {
       (void)fprintf(logFile,"%ld ", newPop[i].age);
   }
   (void)fprintf(logFile,"\n");
#endif /* } DEBUG3 */

   //
   // Copy the n best individuals to the next population, if the elitist flag is set, where n is the value of elitism.
   //
   if (elitism > 0) {
      solutions.msort(elitism);
      for (i=0; i<elitism; i++) {
         newPop[solutions.num_individuals()-1-i] = solutions[i];
      }
   }

#ifdef DEBUG3 /* DEBUG3 { */
   (void)fprintf(logFile,"About to Update the Current Generation, newPop.\n");
   for (i=0; i<solutions.num_individuals(); i++) {
       (void)fprintf(logFile,"%ld ", newPop[i].age);
   }
   (void)fprintf(logFile,"\n");
#endif /* } DEBUG3 */

   //
   // Update current generation 
   //  
   solutions = newPop;

   //
   // Increase the number of generations
   //
   generations++;

   //
   // Increment the age of surviving individuals...
   // 
   for (i=0; i<solutions.num_individuals(); i++) {
       solutions[i].incrementAge();
   }

   if (debug > 0) {
       (void)fprintf(logFile,"DEBUG:  Generation: %3u, outputEveryNgens = %3u, generations%%outputEveryNgens = %u\n",
                     generations, outputEveryNgens, generations%outputEveryNgens);
   }
   if (generations%outputEveryNgens == 0) {
       oldest  = 0L;
       fittest = BIG;
       for (i=0; i<solutions.num_individuals(); i++) {
          if (solutions[i].age >= oldest) {
              oldest = solutions[i].age;
              oldestIndividual = i;
          }
          if (solutions[i].value(Normal_Eval) <= fittest) {
              fittest = solutions[i].value(Normal_Eval);
              fittestIndividual = i;
          }
       }
       /* Only output if the output level is not 0. */
       if (outputEveryNgens != OUTLEV0_GENS) {
           // (void)fprintf(logFile, "___\noutputEveryNgens = %d, OUTLEV0_GENS=%d\n___\n", outputEveryNgens, OUTLEV0_GENS);
           if (outputEveryNgens > 1) {
    #ifndef DEBUG3
               (void)fprintf(logFile,"Generation: %3u   Oldest's energy: %.3f    Lowest energy: %.3f    Num.evals.: %ld   Timing: ", 
               generations, solutions[oldestIndividual].value(Normal_Eval), solutions[fittestIndividual].value(Normal_Eval), 
               evaluate.evals() );
    #else
               (void)fprintf(logFile,"Generation: %3u   Oldest ind.: %u/%u, age: %lu, energy: %.3f    Lowest energy individual: %u/%u, age: %lu, energy: %.3f    Num.evals.: %ld    Timing: ", 
               generations, oldestIndividual+1, solutions.num_individuals(), solutions[oldestIndividual].age, 
               solutions[oldestIndividual].value(Normal_Eval), fittestIndividual+1, solutions.num_individuals(), 
               solutions[fittestIndividual].age, solutions[fittestIndividual].value(Normal_Eval), 
               evaluate.evals() );
    #endif /* DEBUG3 */
           } else {
    #ifndef DEBUG3
               (void)fprintf(logFile,"Generation: %3u   Oldest's energy: %.3f    Lowest energy: %.3f    Num.evals.: %ld   Timing: ", 
               generations, solutions[oldestIndividual].value(Normal_Eval), solutions[fittestIndividual].value(Normal_Eval), evaluate.evals() );
    #else
               (void)fprintf(logFile,"Generation: %3u   Oldest: %u/%u, age: %lu, energy: %.3f    Lowest energy individual: %u/%u, age: %lu, energy: %.3f    Num.evals.: %ld   Timing: ", 
               generations, oldestIndividual+1, solutions.num_individuals(), solutions[oldestIndividual].age, 
               solutions[oldestIndividual].value(Normal_Eval), fittestIndividual+1, solutions.num_individuals(), 
               solutions[fittestIndividual].age, solutions[fittestIndividual].value(Normal_Eval), evaluate.evals() );
    #endif /* DEBUG3 */
           }
       }
       genEnd = times( &tms_genEnd );
       timesyshms( genEnd - genStart, &tms_genStart, &tms_genEnd );
       genStart = times( &tms_genStart );
   }

   return(0);
}
