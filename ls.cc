/*

 $Id: ls.cc,v 1.10.2.1 2010/11/19 20:09:29 rhuey Exp $

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
extern int nlig;		// assign value in mian.cc
extern int ntor_lig[MAX_LIGANDS];
extern int gene_index_lig[MAX_LIGANDS][2];  //gene num start_point & end_point of a ligand.
extern int global_ntor; // set to current s.Init.ntor in main.cc

//  This function adds sign * (array1 + array2) to all the reals in the representation
Phenotype genPhByLigand(const Phenotype &original, 
		Real sign, Real *array1, Real *array2, 
		int from_gene_point, int to_gene_point, int ilig)
{
   RepType genetype;
   register int i, index = 0;
   Phenotype retval(original);

#ifdef DEBUG
   (void)fprintf(logFile, "ls.cc/Phenotype genPh(const Phenotype &original, Real *array1, Real *array2)\n");
   //(void)fprintf(logFile, "ls.cc/Phenotype retval.num_pts() = %d\n", retval.num_pts());
   (void)fprintf(logFile, "ls.cc/Phenotype from_gene_point = %d, to_gene_point = %d\n", from_gene_point, to_gene_point);
#endif // DEBUG 

   // handle multi-ligand fragment  -Huameng 12/06/2007
   //for (i=0; i < retval.num_pts(); i++) {
   for (i = from_gene_point; i < to_gene_point; i++) {
      //genetype = retval.gtype(i);
      //pr( logFile, "DEBUG_genetype = %d \n", genetype);
      //if ((genetype == T_RealV)||(genetype == T_CRealV)) {
         retval.write(retval.gread(i).real + sign * (array1[index] + array2[index]), i);
         index++;
      //}
   }
   
   // multiple ligand Quat  -Huameng 11/18/2007
   if(ilig < nlig) {
   	 Quat q;    	
   	 q = retval.readQuat(ilig);
#ifdef DEBUG_QUAT
   	 pr( logFile, "DEBUG_QUAT: genPh()  q\n" );
   	 printQuat( logFile, q );
   	 assertQuatOK( q );
#endif // endif DEBUG_QUAT
   	 retval.writeQuat( normQuat( q ), ilig);
   }
  
   return(retval);
}


//  This function adds sign * (array1 + array2) to all the reals in the representation
Phenotype genPh(const Phenotype &original, Real sign, Real *array1, Real *array2)
{
   RepType genetype;
   register unsigned int i, index = 0;
   Phenotype retval(original);

#ifdef DEBUG
   (void)fprintf(logFile, "ls.cc/Phenotype genPh(const Phenotype &original, Real *array1, Real *array2)\n");
   (void)fprintf(logFile, "ls.cc/Phenotype retval.num_pts() = %d\n", retval.num_pts());
#endif // DEBUG 

   for (i=0; i < retval.num_pts(); i++) { 
      genetype = retval.gtype(i);
      if ((genetype == T_RealV)||(genetype == T_CRealV)) {
         retval.write(retval.gread(i).real + sign * (array1[index] + array2[index]), i);
         index++;
      }
   }
   // multiple ligand Quat  -Huameng 11/18/2007
   Quat q;
   for (int n = 0; n < nlig; n++) {  	
   	q = retval.readQuat(n);
#ifdef DEBUG_QUAT
	   pr( logFile, "DEBUG_QUAT: genPh()  q\n" );
	   printQuat( logFile, q );
	   assertQuatOK( q );
#endif // endif DEBUG_QUAT
   	retval.writeQuat( normQuat( q ), n );
   	
   }
   return(retval);
}


/**
 *   
 * What Solis & Wets does is add random deviates to every
 *  real number in the Phenotype.
 *   
 *  This has only one value of rho, for all genes.
 */
void Solis_Wets::SWMultiLigand(Phenotype &vector, int from_gene_point, int to_gene_point, int ilig)
{
   register unsigned int i, j, num_successes = 0, num_failures = 0;
   register Real temp_rho = rho;
   Phenotype newPh;
   
#ifdef DEBUG
   (void)fprintf(logFile, "ls.cc/void Solis_Wets::SW(Phenotype &vector)\n");
#endif /* DEBUG */
  
   // handle multi-ligand to do separate local search -Huameng 12/01/2007    
   unsigned int gene_size = to_gene_point - from_gene_point;
       	   	   	    	   
   //  Reset bias	   
   for (i=0; i < gene_size; i++) {	  
      bias[i] = 0.0;
   }
   num_successes = 0;
   num_failures = 0;
   temp_rho = rho;
   for (i=0; i < max_its; i++) {
      // Generate deviates
      for (j=0; j < gene_size; j++) {
         deviates[j] = gen_deviates(temp_rho);
      }
      	  	  
      // zeta = x + bias + deviates
      newPh = genPhByLigand(vector, +1., deviates, bias, from_gene_point, to_gene_point, ilig); // zeta
      // Evaluate
      if (newPh.evaluate(Normal_Eval) < vector.evaluate(Normal_Eval)) {
         num_successes++;
         num_failures = 0;
         vector = newPh;
         for (j=0; j < gene_size; j++) {
            // bias[j] = 0.20*bias[j] + 0.40*deviates[j];  // original & questionable
            bias[j] = 0.60*bias[j] + 0.40*deviates[j]; // strict Solis+Wets
         }
      } else {
         // We need to check if the opposite move does any good (move = bias[j] + deviates[j])
         newPh = genPhByLigand(vector, -1., deviates, bias, from_gene_point, to_gene_point, ilig); // 2x - zeta = x - move
         if (newPh.evaluate(Normal_Eval) < vector.evaluate(Normal_Eval)) {
            num_successes++;
            num_failures = 0;
            vector = newPh;
            for (j=0; j < gene_size; j++) {
               // bias[j] -= 0.40*deviates[j]; // incorrect
               bias[j] = 0.60*bias[j] - 0.40*deviates[j]; // correct if deviates is not changed
            }
         } else {
            num_failures++;
            num_successes = 0;
            // vector is unchanged  // x
            for (j=0; j < gene_size; j++) {
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
   } // i-loop	 (i < max_its)  	  
  
         
} // void Solis_Wets::SWMultiLigand(Phenotype &vector,  )


//  What Solis & Wets does is add random deviates to every
//  real number in the Phenotype.
//  
//  This has only one value of rho, for all genes.
//
void Solis_Wets::SW(Phenotype &vector)
{
   register unsigned int i, j, num_successes = 0, num_failures = 0;
   register Real temp_rho = rho;
   Phenotype newPh;
   
#ifdef DEBUG
   (void)fprintf(logFile, "ls.cc/void Solis_Wets::SW(Phenotype &vector)\n");
#endif /* DEBUG */
  
   // move all points together 
   for (i=0; i < max_its; i++) {
      // Generate deviates
      for (j=0; j < size; j++) {
         deviates[j] = gen_deviates(temp_rho);
      }
      	  	  
      // zeta = x + bias + deviates
      newPh = genPh(vector, +1., deviates, bias); // zeta
      // Evaluate
      if (newPh.evaluate(Normal_Eval) < vector.evaluate(Normal_Eval)) {
         num_successes++;
         num_failures = 0;
         vector = newPh;
         for (j=0; j < size; j++) {
            // bias[j] = 0.20*bias[j] + 0.40*deviates[j];  // original & questionable
            bias[j] = 0.60*bias[j] + 0.40*deviates[j]; // strict Solis+Wets
         }
      } else {
         // We need to check if the opposite move does any good (move = bias[j] + deviates[j])
         newPh = genPh(vector, -1., deviates, bias); // 2x - zeta = x - move
         if (newPh.evaluate(Normal_Eval) < vector.evaluate(Normal_Eval)) {
            num_successes++;
            num_failures = 0;
            vector = newPh;
            for (j=0; j < size; j++) {
               // bias[j] -= 0.40*deviates[j]; // incorrect
               bias[j] = 0.60*bias[j] - 0.40*deviates[j]; // correct if deviates is not changed
            }
         } else {
            num_failures++;
            num_successes = 0;
            // vector is unchanged  // x
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
void Pseudo_Solis_Wets::SWMultiLigand(Phenotype &vector, int from_gene_point, int to_gene_point, int ilig)
{
   register unsigned int i, j, num_successes = 0, num_failures = 0,  all_rho_stepsizes_too_small = 1;
    
   Phenotype newPh;

#ifdef DEBUG
   (void)fprintf(logFile, "ls.cc/void Pseudo_Solis_Wets::SW(Phenotype &vector)\n");
#endif // DEBUG 

   //  Initialize the temp_rho's
   for (i=0; i < size; i++) {
      temp_rho[i] = rho[i];
   }
   //  Reset bias
   for (i=0; i < size; i++) {
      bias[i] = 0.0;
   }
   // handle multi-ligand to do separate local search -Huameng 12/01/2007
   unsigned int gene_size = to_gene_point - from_gene_point;
  
   num_successes = 0;
   num_failures = 0;  	   
   for (i=0; i < max_its; i++) {
      // Generate deviates
      for (j= from_gene_point; j < (unsigned)to_gene_point; j++) {
         deviates[j] = gen_deviates(temp_rho[j]);
      }
	  // handle multi-ligands -Huameng
      //newPh = genPh(vector, +1., deviates, bias);
      newPh = genPhByLigand(vector, +1., deviates, bias, from_gene_point, to_gene_point, ilig);
      // Evaluate
      if (newPh.evaluate(Normal_Eval) < vector.evaluate(Normal_Eval)) {
         num_successes++;
         num_failures = 0;
         vector = newPh;
         for (j=0; j < gene_size; j++) {
            // bias[j] = 0.20*bias[j] + 0.40*deviates[j];
            bias[j] = 0.60*bias[j] + 0.40*deviates[j]; // strict Solis+Wets
         }
      } else  {
         // We need to check if the opposite move does any good (move = bias[j] + deviates[j])
         // handle multi-ligands -Huameng
         //newPh = genPh(vector, -1., deviates, bias);
         newPh = genPhByLigand(vector, -1., deviates, bias, from_gene_point, to_gene_point, ilig);
         if (newPh.evaluate(Normal_Eval) < vector.evaluate(Normal_Eval)) {
            num_successes++;
            num_failures = 0;
            vector = newPh;
            for (j=0; j < gene_size; j++) {
               // bias[j] -= 0.40*deviates[j];
               bias[j] = 0.60*bias[j] - 0.40*deviates[j]; // correct if deviates is not changed
            }
         } else {
            num_failures++;
            num_successes = 0;
            // vector is unchanged  // x
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
      all_rho_stepsizes_too_small = 1;
      for(j=0; j < gene_size; j++) {   
         all_rho_stepsizes_too_small = all_rho_stepsizes_too_small & (temp_rho[j] < lower_bound_on_rho[j]);
      } //  j-loop
      if (all_rho_stepsizes_too_small) {
         break; //  GMM - THIS breaks out of i loop, which IS what we want...
      }
   } //  i-loop
	   
} // void Pseudo_Solis_Wets::SW(Phenotype &vector)



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
#endif // DEBUG 

   //  Initialize the temp_rho's
   for (i=0; i < size; i++) {
      temp_rho[i] = rho[i];
   }
   //  Reset bias
   for (i=0; i < size; i++) {
      bias[i] = 0.0;
   }
   
   num_successes = 0;
   num_failures = 0;  	   
   for (i=0; i < max_its; i++) {
      // Generate deviates
      for (j=0; j < size; j++) {
         deviates[j] = gen_deviates(temp_rho[j]);
      }	  
      newPh = genPh(vector, +1., deviates, bias);
      // Evaluate
      if (newPh.evaluate(Normal_Eval) < vector.evaluate(Normal_Eval)) {
         num_successes++;
         num_failures = 0;
         vector = newPh;
         for (j=0; j < size; j++) {
            // bias[j] = 0.20*bias[j] + 0.40*deviates[j];
            bias[j] = 0.60*bias[j] + 0.40*deviates[j]; // strict Solis+Wets
         }
      } else  {
         //  We need to check if the opposite move does any good (move = bias[j] + deviates[j])        
         newPh = genPh(vector, -1., deviates, bias);
         if (newPh.evaluate(Normal_Eval) < vector.evaluate(Normal_Eval)) {
            num_successes++;
            num_failures = 0;
            vector = newPh;
            for (j=0; j < size; j++) {
               // bias[j] -= 0.40*deviates[j];
               bias[j] = 0.60*bias[j] - 0.40*deviates[j]; // correct if deviates is not changed
            }
         } else {
            num_failures++;
            num_successes = 0;
            // vector is unchanged  // x
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
      all_rho_stepsizes_too_small = 1;
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
   int from_gene_point, to_gene_point;
#ifdef DEBUG
   (void)fprintf(logFile, "ls.cc/int Solis_Wets_Base::search(Individual &solution)\n");
#endif /* DEBUG */

   if (ranf() < search_frequency) {
   	
   	   // handle multi-ligand searching -Huameng
   	   // search rotation & quarternion for each ligand (0, nlig -1),
	   // at last seach rotation segment (n = nlig)
	   for(int n = 0; n < nlig ; n++) {
	   	  //get the segment if gene points for this ligand 	   	   
	   	  from_gene_point = gene_index_lig[n][0];	   
	   	  to_gene_point = gene_index_lig[n][1];  	     	   	 
	   	  // translation and rotation   	   
	      SWMultiLigand(solution.phenotyp, from_gene_point, to_gene_point, n);	      	      	           
	   }
	   
	   // all torsion
	   from_gene_point = nlig*7;
	   to_gene_point = from_gene_point + global_ntor;
	   SWMultiLigand(solution.phenotyp, from_gene_point, to_gene_point, nlig);  
	   /*
	   // Handle the torsion of flexible residue
	   // need to skip this in UNBOUND STATE
	   if(to_gene_point < (global_ntor - 7*nlig)) {	   
		   from_gene_point = gene_index_lig[2*nlig][0];
		   to_gene_point = gene_index_lig[2*nlig][1]; 
		   //pr(logFile, "in side chain local search from_gene_point=%d, to_gene_point=%d\n",
		   //   		from_gene_point, to_gene_point);
		   if((to_gene_point - from_gene_point) > 0 ) { 	        	   
		   	 SWMultiLigand(solution.phenotyp, from_gene_point, to_gene_point, 2*nlig);
		   }	 
	   }
	   */
	   
	   //inverse mapping
       solution.inverse_mapping();
   }

   return(0);
}


/******************************************************
 * 
 * Huameng -7/27/2008 
 * Local search using partical swarm optimization PSO
 * NO inverse_mapping compared to genetic algorithms
 * 
 ****************************************************/
int Solis_Wets_Base::searchByPSO(Individual &solution)
{
#ifdef DEBUG
   (void)fprintf(logFile, "ls.cc-> Start Solis_Wets_Base::searchByPSO(Individual &solution)\n");
#endif /* DEBUG */
   int from_point, to_point;
   if (ranf() < search_frequency) {
   	   	
   	   for (int i = 0; i < nlig; i++) {
   	   	from_point = gene_index_lig[i][0];
   	   	to_point = gene_index_lig[i][1];
   	   	// handle multi-ligand searching -Huameng  	   	   	       	   
	   	SWMultiLigand(solution.phenotyp, from_point, to_point, i);	      	      	           
	  	   	   	   
#ifdef DEBUG
	   	fprintf(logFile, "ls.cc->  End local searchByPSO(Individual &solution)\n"); 
	   	fflush(logFile); 
#endif /* DEBUG */
   	   } 
   	   
   	   // all torsion
	   from_point = nlig*7;
	   to_point = from_point + global_ntor;
	   SWMultiLigand(solution.phenotyp, from_point, to_point, nlig);  
   }

   return(0);
}


Pattern_Search::Pattern_Search(void)
{
}

Pattern_Search::Pattern_Search(unsigned int init_size, unsigned int init_max_success, Real init_step_size, Real init_step_threshold, Real init_expansion, Real init_contraction, Real init_search_frequency)
: size(init_size), max_success(init_max_success), step_size(init_step_size), step_threshold(init_step_threshold), expansion(init_expansion), contraction(init_contraction), search_frequency(init_search_frequency)
{
  current_step_size = step_size;
  pattern = new Real[size];
	index = new unsigned int[size];
  reset_pattern();
	reset_indexes();
  successes = 0;
}

Pattern_Search::~Pattern_Search(void)
{
	delete []pattern;
	delete []index;
}

void Pattern_Search::reset()
{
  current_step_size = step_size;
  reset_pattern();
  reset_indexes();
  successes = 0;
}

void Pattern_Search::reset_pattern() {
  for (unsigned int i=0; i < size; i++) {
    pattern[i] = 0.0;
  }
}

void Pattern_Search::reset_indexes() {
	for (unsigned int i=0; i < size; i++) {
		index[i] = i;
	}
}

void Pattern_Search::shuffle_indexes() {
	int select;
	unsigned int temp;
	for (unsigned int i=size; i > 1; i--) {
		select = rand() % i;
		temp = index[select];
		index[select] = index[i-1];
		index[i-1] = temp;
	}
}

int Pattern_Search::terminate(void)
{
   return (0);
}

int Pattern_Search::search(Individual &solution)
{
  // TODO: implement scaling?

  if (ranf() >= search_frequency) {
    return(0);
  }

	reset();
  Phenotype base = solution.phenotyp;
  Phenotype newPoint;
  // evaluate function at base point
  while (current_step_size > step_threshold) {
    // do exploratory moves
    //fprintf(stderr, "base point energy: %f\n", base.evaluate(Normal_Eval));
    newPoint = exploratory_move(base);
    //fprintf(stderr, "newPoint energy: %f\n", newPoint.evaluate(Normal_Eval));
    if (newPoint.evaluate(Normal_Eval) < base.evaluate(Normal_Eval)) {
      // new point is more favorable than base point
      // set new point as base point
      base = newPoint;

      while (true) {
        newPoint = pattern_explore(base);
        if (newPoint.evaluate(Normal_Eval) < base.evaluate(Normal_Eval)) {
			successes++;
          base = newPoint;
        } else {
			break;
			successes = 0;
		}
		if (successes > max_success) {
			//fprintf(stderr, "Expanding step size\n");
			successes = 0;
			current_step_size *= expansion;
		}
      } //while
    }
    else {
      current_step_size *= contraction;
			successes = 0;
      reset_pattern();
      //fprintf(stderr, "Contracted to %f after %ld evaluations.\n", current_step_size, evaluate.evals());
    }
  }
  
  solution.phenotyp = base;
  solution.inverse_mapping();
  return (0);
}

int Pattern_Search::searchByPSO(Individual &solution)
{
  // TODO: implement scaling?

  if (ranf() >= search_frequency) {
    return(0);
  }

	reset();
  Phenotype base = solution.phenotyp;
  Phenotype newPoint;
  // evaluate function at base point
  while (current_step_size > step_threshold) {
    // do exploratory moves
    //fprintf(stderr, "base point energy: %f\n", base.evaluate(Normal_Eval));
    newPoint = exploratory_move(base);
    //fprintf(stderr, "newPoint energy: %f\n", newPoint.evaluate(Normal_Eval));
    if (newPoint.evaluate(Normal_Eval) < base.evaluate(Normal_Eval)) {
      // new point is more favorable than base point
      // set new point as base point
      base = newPoint;

      while (true) {
        newPoint = pattern_explore(base);
        if (newPoint.evaluate(Normal_Eval) < base.evaluate(Normal_Eval)) {
			successes++;
          base = newPoint;
        } else {
			break;
			successes = 0;
		}
		if (successes > max_success) {
			//fprintf(stderr, "Expanding step size\n");
			successes = 0;
			current_step_size *= expansion;
		}
      } //while
    }
    else {
      current_step_size *= contraction;
			successes = 0;
      reset_pattern();
      //fprintf(stderr, "Contracted to %f after %ld evaluations.\n", current_step_size, evaluate.evals());
    }
  }
  
  solution.phenotyp = base;
  return (0);
}

Phenotype Pattern_Search::exploratory_move(const Phenotype& base) {
  Phenotype newBase(base);
	shuffle_indexes();
	unsigned int current_index;
	int direction;

  for (unsigned int i=0; i < size; i++) {
    Phenotype trialPoint(newBase);

		current_index = index[i];
		// pick a random direction
		if (rand()%2 == 0) {
			direction = 1;
		}
		else {
			direction = -1;
		}
    // try first coordinate direction
    trialPoint.write(trialPoint.gread(current_index).real+current_step_size*direction, current_index);
    // if successful, keep new point
    if (trialPoint.evaluate(Normal_Eval) < newBase.evaluate(Normal_Eval)) {
      newBase = trialPoint;
      pattern[current_index] += current_step_size*direction;
    }
    // otherwise, try opposite coordinate and test again
    else {
      trialPoint.write(trialPoint.gread(current_index).real-2.0*current_step_size*direction, current_index);
      if (trialPoint.evaluate(Normal_Eval) < newBase.evaluate(Normal_Eval)) {
        newBase = trialPoint;
        pattern[current_index] -= current_step_size*direction;
      }
    }
  }
  return newBase;
}

Phenotype Pattern_Search::pattern_explore(const Phenotype& base) {
  Phenotype newPoint = pattern_move(base);
  reset_pattern();
  Phenotype newBase = exploratory_move(newPoint);
  return newBase;
}

Phenotype Pattern_Search::pattern_move(const Phenotype& base) {
  Phenotype newPoint(base);
  for (unsigned int i=0; i < size; i++) {
    newPoint.write(newPoint.gread(i).real + pattern[i] , i);
  }
  return newPoint;
}
