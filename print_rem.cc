/*

 $Id: print_rem.cc,v 1.2 2003/02/26 01:29:27 garrett Exp $

*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/* print_rem.cc */ 
    #include <stdio.h>
    #include "print_rem.h"


void print_rem( FILE *outFile,
		int Rank,
		int NumMem,
		int Run,
		FloatOrDouble ref_rms)
{
    fprintf( outFile, "MODEL     %4d\n", Run );
    fprintf( outFile, "USER    Run = %d\n", Run );
    fprintf( outFile, "USER    Cluster Rank = %d\n", Rank );
    fprintf( outFile, "USER    Number of conformations in this cluster = %d\n", NumMem );
    fprintf( outFile, "USER  \n");
    fprintf( outFile, "USER    RMSD from reference structure       = %.3f A\n", ref_rms );
    fprintf( outFile, "USER  \n");
}
/* EOF */
