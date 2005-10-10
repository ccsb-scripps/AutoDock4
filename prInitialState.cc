/*

 $Id: prInitialState.cc,v 1.3.6.1 2005/10/10 23:47:05 alther Exp $

*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/* prInitialState.cc */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "prInitialState.h"


extern int keepresnum;
extern FILE *logFile;
extern char *programname;


void prInitialState(

    FloatOrDouble einter,
    FloatOrDouble eintra,
    FloatOrDouble torsFreeEnergy,
    int natom,
    FloatOrDouble crd[MAX_ATOMS][SPACE],
    char atomstuff[MAX_ATOMS][MAX_CHARS],
    int type[MAX_ATOMS],
    FloatOrDouble emap[MAX_ATOMS],
    FloatOrDouble elec[MAX_ATOMS],
    FloatOrDouble charge[MAX_ATOMS],
    int ligand_is_inhibitor)

{
    char rec8[10];
    char rec13[15];
    char descriptor[17];
    register int i = 0;
    int a = 0;
    FloatOrDouble emap_total = 0.0;
    FloatOrDouble elec_total = 0.0;

    strncpy(descriptor, "INITIAL STATE:  ", (size_t)16);

    pr( logFile, "\n\t\t%s\n\t\t______________\n\n\n", descriptor );

    pr( logFile, "%sUSER    Transformed Initial Coordinates\n", descriptor );
    for (i = 0;  i < natom;  i++) {
        pr( logFile, "%s", descriptor);
	if (keepresnum > 0) {
	    strncpy( rec13, &atomstuff[i][13], (size_t)13);
	    pr(logFile, FORMAT_PDBQ_ATOM_RESSTR, "", i+1, rec13,   crd[i][X], crd[i][Y], crd[i][Z], 1.0, 0.0, charge[i]);
	    pr(logFile, "\n");
	} else {
	    strncpy( rec8, &atomstuff[i][13], (size_t)8);
	    pr(logFile, FORMAT_PDBQ_ATOM_RESNUM, "", i+1, rec8, 0, crd[i][X], crd[i][Y], crd[i][Z], 1.0, 0.0, charge[i]);
	    pr(logFile, "\n");
	}
    } /* i */
    pr( logFile, "%sTER\n\n\n", descriptor );

    pr( logFile, "\t\tINITIAL ENERGY BREAKDOWN\n" );
    pr( logFile, "\t\t________________________\n" );
    pr( logFile, "\n\nEnergy of starting position of Small Molecule by atom: \n\n" );

    print_atomic_energies( natom, atomstuff, type, emap, elec, charge );

    emap_total = 0.0;
    elec_total = 0.0;
    for (a=0; a<natom; a++) {
        emap_total += emap[a];
        elec_total += elec[a];
    }
    printEnergies( einter, eintra, torsFreeEnergy, "Initial ", ligand_is_inhibitor, emap_total, elec_total);

    flushLog;
}
/* EOF */
