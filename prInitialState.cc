/*

 $Id: prInitialState.cc,v 1.8 2007/04/27 06:01:50 garrett Exp $

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

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "prInitialState.h"


extern int keepresnum;
extern FILE *logFile;
extern char *programname;


void prInitialState(
    EnergyBreakdown *p_eb,
    int natom,
    Real crd[MAX_ATOMS][SPACE],
    char atomstuff[MAX_ATOMS][MAX_CHARS],
    int type[MAX_ATOMS],
    Real emap[MAX_ATOMS],
    Real elec[MAX_ATOMS],
    Real charge[MAX_ATOMS],
    int ligand_is_inhibitor,
    Boole B_have_flexible_residues )

{
    char rec8[10];
    char rec13[15];
    char descriptor[17];
    register int i = 0;
    int a = 0;
    Real emap_total = 0.0;
    Real elec_total = 0.0;

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
    
	pr( logFile, "\n\n" );
    printEnergies( p_eb, "Initial ", ligand_is_inhibitor, emap_total, elec_total, B_have_flexible_residues );
    pr( logFile, "\n\n" );

    flushLog;
}
/* EOF */
