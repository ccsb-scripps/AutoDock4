/* summarizegrids.cc */

    #include <stdio.h>
    #include "summarizegrids.h"


extern FILE *logFile;

void summarizegrids( char  atm_typ_str[ATOM_MAPS],
		     float mapmax[MAX_MAPS],
		     float mapmin[MAX_MAPS],
		     int   num_all_maps,
		     int   num_atm_maps )
{
    /*
     Summarize in the log-file what was found:
    */

    int i = 0;

    fprintf( logFile, "\n" );
    fprintf( logFile, UnderLine );
    fprintf( logFile, "SUMMARY OF GRID MAPS' ENERGIES:\n_______________________________\n\n\t\t\tvan der Waals or H-bond\n\tGrid\tAtom\t    Potential Energy\n\tMap\tType\tMinimum\t\tMaximum\n\t___\t____\t_________\t_________\n" );

    for ( i=0; i < num_atm_maps ; i++ ) {

	fprintf( logFile, "\t %d\t %c\t", (i+1), atm_typ_str[i]);

	/*fprintf( logFile, "<DEBUG 2>\n");  */

	fprintf( logFile, "%6.2f\t\t", mapmin[i]); 

	fflush( logFile );

	/* fprintf( logFile, "<DEBUG 3>\n");  */

	fflush( logFile );

	fprintf( logFile, "%6.2e\n", mapmax[i]);

	fflush( logFile );

	/* fprintf( logFile, "<DEBUG 4>\n");  */

	fflush( logFile );
    }

    fprintf( logFile, "\n\t\t\t    Electrostatic\n\tGrid\tGrid\t    Potential Energy:\n\tMap\tType\tMinimum\t\tMaximum\n\t___\t____\t_________\t_________\n" );

    fprintf( logFile, "\t %d\telec\t%+8.2f\t%+8.2f\n\n", num_all_maps, mapmin[num_atm_maps],mapmax[num_atm_maps]);

    /* fprintf( logFile, UnderLine ); */

    fflush( logFile );
}
/* EOF */
