/* print_atomic_energies.cc */

    #include <stdio.h>
    #include <string.h>
    #include "print_atomic_energies.h"


extern FILE *logFile;

/*----------------------------------------------------------------------------*/
void print_atomic_energies( int natom, 
			    char atomstuff[MAX_ATOMS][MAX_CHARS],
			    int type[MAX_ATOMS],
			    float emap[MAX_ATOMS],
			    float elec[MAX_ATOMS],
			    float charge[MAX_ATOMS] )

/*----------------------------------------------------------------------------*/
{
    char rec[16];
    register int i;

    fprintf( logFile, 
     "Small Molecule  Atom  Non-bonded   Electrostatic  Partial\n" );
    fprintf( logFile, 
     "Atom & Residue  Type    Energy        Energy      Charge\n" );
    fprintf( logFile, 
     "______________  ____  ___________  _____________  ______\n" );

    for (i = 0;  i < natom;  i++) {
        strncpy( rec, &atomstuff[i][6], (size_t)14 );
        fprintf( logFile, "%.14s   %1d  %+11.2f  %+11.2f      %+6.3f\n", rec, (type[i]+1), emap[i], elec[i], charge[i] );
    }
}
/* EOF */
