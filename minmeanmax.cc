/* minmeanmax.cc */

/* Calculates the min, mean and max of each gene for a given population */

/* N.B.: the macro is defined in "hybrids.h" */

#include <stdio.h>
#include "rep.h"
#include "support.h" /* Population defined here */

void minmeanmax( FILE *fp, Population &pop, int num_its )
{
   register int i=0, g=0;
   int     num_indvs=1, num_genes=1;
   double  state_var=0., *minimum, *sum, *maximum, *best;
   // double  best_energy=1e20;
   // double energy;
   Element temp;

   num_indvs = pop.num_individuals();
   num_genes = pop[0].genotyp.num_genes();

   minimum = new double[ num_genes ];
   sum     = new double[ num_genes ];
   maximum = new double[ num_genes ];
   best    = new double[ num_genes ];

   // energy  = best_energy;

   for (g=0; g < num_genes; g++)
   {
      temp       = pop[0].genotyp.gread(g);
      minimum[g] = temp.real;
      sum[g]     = temp.real;
      maximum[g] = temp.real;
      best[g]    = temp.real;
   }

   for (i=1; i < num_indvs; i++) // i-th individual in population
   {
      for (g=0; g < num_genes; g++) // g-th gene in individual
      {
         temp = pop[i].genotyp.gread(g);
         state_var = temp.real;
         if (minimum[g] > state_var) {
            minimum[g] = state_var;
         }
         if (maximum[g] < state_var) {
            maximum[g] = state_var;
         }
         sum[g] += state_var;
      }
   }
   for (g=0; g < num_genes; g++) // g-th gene in individual 
   {
       fprintf( fp, "%d  %2d  %8.3f  %8.3f  %8.3f\n", num_its, (g+1), minimum[g], sum[g]/num_indvs, maximum[g] ); fflush( fp );
   }

   delete [] minimum;
   delete [] sum;
   delete [] maximum;
}
