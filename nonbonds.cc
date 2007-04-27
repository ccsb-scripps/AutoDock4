/*

 $Id: nonbonds.cc,v 1.10 2007/04/27 06:01:50 garrett Exp $

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
#include <stdio.h>
#include "nonbonds.h"
#include "mdist.h"  /* mindist and maxdist are here */


extern int debug;
extern  FILE    *logFile;

#ifdef DEBUG
extern int debug;
extern  FILE    *logFile;
#endif /* DEBUG */

using namespace std;

void nonbonds(const Real  crdpdb[MAX_ATOMS][SPACE],
		      int         nbmatrix[MAX_ATOMS][MAX_ATOMS],
		      const int   natom, 
              const int   bond_index[MAX_ATOMS],
              int         B_include_1_4_interactions,
              int         bonded[MAX_ATOMS][6])
{
	int i,j,k,l;
    int nonbond_type;

    if (debug>0) {
        printbonds(natom, bonded, "\nDEBUG:  3. INSIDE nonbonds, bonded[][] array is:\n\n", 0);
    }

    //
    // in "nbmatrix", the values 1 (and 4) mean this pair of atoms will be included in the internal, non-bonded list
    //                           0                                         ignored
 
    // set all nonbonds in nbmatrix to 1, except "1-1 interactions" (self interaction)
    for (i = 0; i<MAX_ATOMS; i++) {
        for (j = 0; j<MAX_ATOMS; j++){
            nbmatrix[i][j] = 1;
        } // j
        nbmatrix[i][i] = 0; /* 2005-01-10 RH & GMM */
    } // i

	for (i = 0; i<natom; i++) {  // loop over atoms, "i", from 1 to natom
        for (j=0; j<bonded[i][5]; j++) {  // loop over this atom's "j" bonds
            // Ignore 1-2 Interactions
            nbmatrix[ i            ][ bonded[i][j] ] = 0;
            nbmatrix[ bonded[i][j] ][ i            ] = 0;
        } // j
    } // i

	for (i=0; i<natom; i++) {  // loop over all atoms, "i"
		for (j=0; j<bonded[i][5]; j++) {  // loop over each atom "j" bonded to the current atom, "i"
			for (k=0; k<bonded[bonded[i][j]][5]; k++) {  // loop over each atom "k" bonded to the current atom "j"

				// Ignore "1-3 Interactions"
				nbmatrix[bonded[bonded[i][j]][k]][i] = 0;
				nbmatrix[i][bonded[bonded[i][j]][k]] = 0;

                if (B_include_1_4_interactions == FALSE) {
                    // include is FALSE,  i.e.: ignore
                    // Ignore "1-4 Interactions"
                    nonbond_type = 0;
                } else {
                    // remember that this is a 1-4 interaction so we can scale it later...
                    // Remember "1-4 Interactions"
                    nonbond_type = 4;
                } //  end if we are ignoring "1-4 interactions"

                // Loop over each atom "l" bonded to the atom "k" bonded to the atom "j" bonded to the atom "i"...
                for (l=0; l<bonded[bonded[bonded[i][j]][k]][5]; l++) { 
                    nbmatrix[i][bonded[bonded[bonded[i][j]][k]][l]] = nonbond_type;
                    nbmatrix[bonded[bonded[bonded[i][j]][k]][l]][i] = nonbond_type;
                } //  l


			} //  k
		} //  j
	} //  i
	return;
} // end of nonbonds

/*----------------------------------------------------------------------------*/

void getbonds(const Real crdpdb[MAX_ATOMS][SPACE],
              const int from_atom,
              const int to_atom,
              const int bond_index[MAX_ATOMS],
              int   bonded[MAX_ATOMS][6])
{
	int i,j;
	double dist,dx,dy,dz;

    // set up all the minimum and maximum possible distances for bonds
	mdist();

    // determine the bonded atoms or "1-2 interactions"

    // loop over atoms, "i", from "from_atom" to "to_atom"
	for (i = from_atom;  i < to_atom;  i++) {

        // while on atom"i", loop over atoms "j", from "i+1" to "to_atom"
		for (j = i+1;  j < to_atom;  j++) {
	
			dx = crdpdb[i][X] - crdpdb[j][X];
			dy = crdpdb[i][Y] - crdpdb[j][Y];
			dz = crdpdb[i][Z] - crdpdb[j][Z];

			dist = sqrt(dx*dx + dy*dy + dz*dz);  // calculate the distance from "i" to "j"

            // if distance from "i" to "j" is in range for their atom types, 
            // set one of the atoms to which "i" is bonded to "j", and vice-versa.
			if (dist >= mindist[bond_index[i]][bond_index[j]] && 
                dist <= maxdist[bond_index[i]][bond_index[j]]) {  

                // bonded[x][5] is the current number of bonds that atom "x" has.
	            // Remember:   int bonded[MAX_ATOMS][6];	

				bonded[i][ bonded[i][5] ] = j;
				bonded[j][ bonded[j][5] ] = i;
				bonded[i][5] += 1;
				bonded[j][5] += 1;
				
			} //  dist is in bonding range for these atom types
		} //  j
 	} //  i

	return;
} // end of get bonds

/*----------------------------------------------------------------------------*/

void printbonds(const int natom, const int bonded[MAX_ATOMS][6], const char *message, const int B_print_all_bonds)
{
    register int i, j;
    pr(logFile, message);
    for (i = 0; i<natom; i++) {  // loop over atoms, "i", from 1 to natom
        pr(logFile, "DEBUG:  atom %d  bonded to ", i+1);
        if (B_print_all_bonds == 1) {
            for (j=0; j<6; j++) {  // loop over all the slots for this atom's "j" bonds
                pr(logFile, "  %d", bonded[i][j]+1);
            } // j
        } else {
            for (j=0; j<bonded[i][5]; j++) {  // loop over this atom's "j" bonds
                pr(logFile, "  %d", bonded[i][j]+1);
            } // j
        }
        pr(logFile, "\n");
    } // i
} // end of printbonds

/*----------------------------------------------------------------------------*/

void print_1_4_message(FILE *file, Boole B_include_1_4_interactions,  Real scale_1_4)
{
    if (B_include_1_4_interactions == FALSE) {
        pr(file, "1,4-interactions will be _ignored_ in the non-bonded internal energy calculation.\n\n");
    } else {
        pr(file, "1,4-interactions will be _included_ in the non-bonded internal energy calculation.\n\n");
        pr(file, "1,4-interaction energies will be will be scaled by a factor of %.2lf .\n\n", (double)scale_1_4);
        pr(file, "NOTE:  Computed internal energies will differ from the standard AutoDock free energy function.\n\n");
    }
    (void) fflush(file);
} // end of print_1_4_message

/* EOF */
