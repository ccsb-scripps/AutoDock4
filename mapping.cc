/********************************************************************
     These are the user defined functions that perform the
     mapping between Genotype and Phenotype and its inverse

				rsh 9/95
********************************************************************/
#include "support.h"

extern FILE *logFile;

//  This should be made more efficient.  As it is now, we (de facto) AlwaysEval!!!!
Phenotype Individual::mapping(void)
{

#ifdef DEBUG
   (void)fprintf(logFile, "mapping.cc/Phenotype Individual::mapping(void)\n");
#endif /* DEBUG */

   phenotyp.write(*genotyp.vread(0), 0);
   phenotyp.write(*genotyp.vread(1), 1);
   phenotyp.write(*genotyp.vread(2), 2);
   phenotyp.write(*genotyp.vread(3), 3);
   phenotyp.write(*genotyp.vread(4), 4);
   value(Normal_Eval);

   return(phenotyp);
}

Genotype Individual::inverse_mapping(void)
{

#ifdef DEBUG
   (void)fprintf(logFile, "mapping.cc/Genotype Individual::inverse_mapping(void)\n");
#endif /* DEBUG */

   genotyp.write(*phenotyp.vread(0), 0);
   genotyp.write(*phenotyp.vread(1), 1);
   genotyp.write(*phenotyp.vread(2), 2);
   genotyp.write(*phenotyp.vread(3), 3);
   genotyp.write(*phenotyp.vread(4), 4);

   return(genotyp);
}
