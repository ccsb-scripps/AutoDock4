/* sort_enrg.cc */

    #include "sort_enrg.h"


extern FILE *logFile;

void sort_enrg( float econf[MAX_RUNS],
                int isort[MAX_RUNS],
		int nconf )

{
/*__________________________________________________________________________
 | Sort conformations on energy                                             |
 |__________________________________________________________________________|
 | Searches through all conformations;  puts in isort[0] the index of the   |
 | lowest energy, in isort[1] the next lowest energy's index, and so on.    |
 |__________________________________________________________________________|
 | WARNING: Fails if any 2 or more econf[] energies are equal.              |
 |__________________________________________________________________________|*/

    quicksort( econf, isort, 0, nconf-1 );
}
/* EOF */
