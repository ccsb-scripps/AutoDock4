/*

 $Id: ls.cc,v 1.4 2004/11/16 23:42:53 garrett Exp $

*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/********************************************************************
      These are the methods of the local searchers

				rsh 9/95

      Modifications to the class heirarchy made 2/21/96 rsh
********************************************************************/
#include "ls.h"
extern class Eval evaluate;

extern FILE *logFile;

//  This function adds array1 + array2 to all the reals in the representation
Phenotype genPh(const Phenotype &original, FloatOrDouble *array1, FloatOrDouble *array2)
{
   RepType genetype;
   register unsigned int i, index = 0;
   Phenotype retval(original);

#ifdef DEBUG
   (void)fprintf(logFile, "ls.cc/Phenotype genPh(const Phenotype &original, FloatOrDouble *array1, FloatOrDouble *array2)\n");
#endif /* DEBUG */


   for (i=0; i < retval.num_pts(); i++) {
      genetype = retval.gtype(i);
      if ((genetype == T_RealV)||(genetype == T_CRealV)) {
         retval.write(retval.gread(i).real + array1[index] + array2[index], i);
         index++;
      }
   }

   return(retval);
}

//  What Solis & Wetts does is add random deviates to every
//  real number in the Phenotype.
//  
//  This has only one value of rho, for all genes.
//
void Solis_Wets::SW(Phenotype &vector)
{
   register unsigned int i, j, num_successes = 0, num_failures = 0;
   register FloatOrDouble temp_rho = rho;
   Phenotype newPh;

#ifdef DEBUG
   (void)fprintf(logFile, "ls.cc/void Solis_Wets::SW(Phenotype &vector)\n");
#endif /* DEBUG */

   //  Reset bias
   for (i=0; i < size; i++) {
      bias[i] = 0.0;
   }

   for (i=0; i < max_its; i++) {
      // Generate deviates
      for (j=0; j < size; j++) {
         deviates[j] = gen_deviates(temp_rho);
      }

      newPh = genPh(vector, deviates, bias);
      // Evaluate
      if (newPh.evaluate(Normal_Eval) < vector.evaluate(Normal_Eval)) {
         num_successes++;
         num_failures = 0;
         vector = newPh;
         for (j=0; j < size; j++) {
            bias[j] = 0.20*bias[j] + 0.40*deviates[j];
         }
      } else  {
         //  We need to check if the opposite deviates do any good
         for (j=0; j < vector.num_pts(); j++) {
            deviates[j] *= -1.0;
         }

         newPh = genPh(vector, deviates, bias);
         if (newPh.evaluate(Normal_Eval) < vector.evaluate(Normal_Eval)) {
            num_successes++;
            num_failures = 0;
            vector = newPh;
            for (j=0; j < size; j++) {
               bias[j] -= 0.40*deviates[j];
            }
         } else {
            num_failures++;
            num_successes = 0;
            for (j=0; j < size; j++) {
               bias[j] *= 0.50;
            }
         }
      }

      // Check to see if we need to expand or contract
      if (num_successes >= max_successes) {
         temp_rho *= expansion;
         num_successes = num_failures = 0;
      } else if (num_failures >= max_failures) {
         temp_rho *= contraction;
         num_successes = num_failures = 0;
      }
         
      if (temp_rho < lower_bound_on_rho)
         break;  // GMM - this breaks out of the i loop...
   } // i-loop
} // void Solis_Wets::SW(Phenotype &vector)


//  This is pseudo-Solis & Wets in that it adds random deviates to every dimension
//  of the current solution, but the variances vary across dimensions.
//
//  This has a different value of rho for each gene.
//
void Pseudo_Solis_Wets::SW(Phenotype &vector)
{
   register unsigned int i, j, num_successes = 0, num_failures = 0,  all_rho_stepsizes_too_small = 1;
    
   Phenotype newPh;

#ifdef DEBUG
   (void)fprintf(logFile, "ls.cc/void Pseudo_Solis_Wets::SW(Phenotype &vector)\n");
#endif /* DEBUG */

   //  Initialize the temp_rho's
   for (i=0; i < size; i++) {
      temp_rho[i] = rho[i];
   }
   //  Reset bias
   for (i=0; i < size; i++) {
      bias[i] = 0.0;
   }

   for (i=0; i < max_its; i++) {
      // Generate deviates
      for (j=0; j < size; j++) {
         deviates[j] = gen_deviates(temp_rho[j]);
      }

      newPh = genPh(vector, deviates, bias);
      // Evaluate
      if (newPh.evaluate(Normal_Eval) < vector.evaluate(Normal_Eval)) {
         num_successes++;
         num_failures = 0;
         vector = newPh;
         for (j=0; j < size; j++) {
            bias[j] = 0.20*bias[j] + 0.40*deviates[j];
         }
      } else  {
         //  We need to check if the opposite deviates do any good
         for (j=0; j < vector.num_pts(); j++) {
            deviates[j] *= -1.0;
         }

         newPh = genPh(vector, deviates, bias);
         if (newPh.evaluate(Normal_Eval) < vector.evaluate(Normal_Eval)) {
            num_successes++;
            num_failures = 0;
            vector = newPh;
            for (j=0; j < size; j++) {
               bias[j] -= 0.40*deviates[j];
            }
         } else {
            num_failures++;
            num_successes = 0;
            for (j=0; j < size; j++) {
               bias[j] *= 0.50;
            }
         }
      }

      // Check to see if we need to expand or contract
      if (num_successes >= max_successes) {
         for(j=0; j < size; j++) {
            temp_rho[j] *= expansion;
         }
         num_successes = num_failures = 0;
      } else if (num_failures >= max_failures) {
         for(j=0; j < size; j++) {
            temp_rho[j] *= contraction;
         }
         num_successes = num_failures = 0;
      }
      
      //  WEH - Scott's code doesn't do anything!!! no stopping based upon step scale!!!
      //  GMM - corrected Scott's code; this does now stop correctly, based upon step scale.
      //  GMM - This version only exits if all the step sizes are too small...
      for(j=0; j < size; j++) {   
         all_rho_stepsizes_too_small = all_rho_stepsizes_too_small & (temp_rho[j] < lower_bound_on_rho[j]);
      } //  j-loop
      if (all_rho_stepsizes_too_small) {
         break; //  GMM - THIS breaks out of i loop, which IS what we want...
      }
   } //  i-loop
} // void Pseudo_Solis_Wets::SW(Phenotype &vector)


int Solis_Wets_Base::search(Individual &solution)
{

#ifdef DEBUG
   (void)fprintf(logFile, "ls.cc/int Solis_Wets_Base::search(Individual &solution)\n");
#endif /* DEBUG */

   if (ranf() < search_frequency) {
      SW(solution.phenotyp);
      solution.inverse_mapping();
   }

   return(0);
}
