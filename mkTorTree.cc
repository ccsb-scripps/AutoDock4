/*

 $Id: mkTorTree.cc,v 1.13 2007/04/27 06:01:50 garrett Exp $

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
#include <ctype.h>
#include "mkTorTree.h"
#include "PDBQT_tokens.h"


extern int true_ligand_atoms;
extern FILE *logFile;
extern char  *programname;
    

void mkTorTree( int   atomnumber[ MAX_RECORDS ],
                char  Rec_line[ MAX_RECORDS ][ LINE_LEN ],
                int   nrecord,

                int   tlist[ MAX_TORS ][ MAX_ATOMS ],
                int   *P_ntor,
                int   *P_ntor_ligand,

                char  smFileName[ MAX_CHARS ],

                char  pdbaname[ MAX_ATOMS ][ 5 ],
                Boole *P_B_constrain,
                int   *P_atomC1,
                int   *P_atomC2,
                Real  *P_sqlower,
                Real  *P_squpper,
                int   *P_ntorsdof,
                int   ignore_inter[MAX_ATOMS])

{

    int   atomlast = 0;
    int   ii = 0;
    int   imax = 0;
    int   itor = 0;
    int   keyword_id = -1;
    int   nbranches = 0;
    int   ntor=0;
    int   tlistsort[ MAX_TORS ][ MAX_ATOMS ];
    int   found_new_res = 0;
    int   nres = 0;
    int   natoms_in_res = 0;
    int   natoms = 0;
    int   nrestor=0;
    int   found_first_res=0;

    register int   i = 0;
    register int   j = 0;
    register int   k = 0;

    char  error_message[ LINE_LEN ];
    char  rec5[ 5 ];

    Real lower = 0.;
    Real temp  = 0.;
    Real upper = 0.01;

#ifdef DEBUG
    int   oo = 0;
    char  C = 'A';
#endif /* DEBUG */

    /* ntor = *P_ntor; */

    /* TorsionTree torstree; */

    for (i = 0; i  < MAX_TORS;  i++ ) {
        for (j = 0;  j < MAX_ATOMS;  j++ ) {
            tlistsort[ i ][ j ] = 0;
        }
    }

    i = 0;
    j = 0;

    /* ________________________________________________________________
      |  Work out the torsion angle tree.                              |
      |________________________________________________________________|
      |  tlist contains information on torsion angles:                 |
      |                                                                |
      |  tlist[ i ][ ATM1 ]           = | atom numbers defining torsion|
      |  tlist[ i ][ ATM2 ]             |                              |
      |  tlist[ i ][ NUM_ATM_MOVED ]  = number of atoms to be rotated  |
      |  tlist[ i ][ 3 ] and on       = atom IDs to be rotated.        |
      |________________________________________________________________|
      | NOTE:  code does not explicitly check for keyword 'ROOT'.      |
      |________________________________________________________________|
    */

#ifdef DEBUG
    pr( logFile, "\n\n" );
    pr( logFile, "                                                  |Atoms|Total\n" );
    pr( logFile, "                                                  | 1| 2|Moved\n" );
    pr( logFile, "   Atom                                           |__|__|__|\n" );
    pr( logFile, "i  #  rec5 C atomlast nbranches j  ntor tlist[  ] [ 0| 1| 2| 3  4  5  6  7  8  9  10 11 12 13 14 15 ]\n" );
    pr( logFile, "__ __ ____ _ ________ _________ __ ____ _____________________________________________________________\n" );
#endif /* DEBUG */

    for (i = 0;  i < nrecord;  i++) {

        if ( (keyword_id = parse_PDBQT_line( Rec_line[ i ] )) == -1) {
            pr( logFile, "%s: Unrecognized keyword found while parsing PDBQT file, line:\n|%s|\n", programname, Rec_line[ i ] );
            continue;
        }

#ifdef DEBUG
        pr( logFile, "PDBQT-Line %d: %s", i+1, Rec_line[i] );
#endif /* DEBUG */

        switch( keyword_id ) {

    /*____________________________________________________________*/
            case PDBQ_REMARK:

                pr( logFile, "%s", Rec_line[ i ] );
                break;

    /*____________________________________________________________*/
            case PDBQ_NULL:

                break;

    /*____________________________________________________________*/
            case PDBQ_ATOM: 
            case PDBQ_HETATM:

             /* This is an ATOM or HETATM. */

             atomlast = atomnumber[ i ];

             if (found_new_res) {
                    /* We are in a residue. */
                    if (natoms_in_res < 2) {
                        /* flag the first two atoms in each new 
                         * residue to prevent them being
                         * included in the intermolecular
                         * energy calculation.  */
                         ignore_inter[natoms] = 1;
                    }
                    /* Keep counting the number of atoms in the residue. */
                    natoms_in_res++;
                } else {
                    /* We are not in a residue.
                     *
                     * "found_new_res" can only be reset to FALSE 
                     * if we encounter an "END_RES" record. 
                     * By default, found_new_res is set to FALSE. */
                     
                    /* reset the atom counter */
                    natoms_in_res = 0;
                }
                /* Increment atom counter for all atoms in PDBQT file */
                natoms++;

#ifdef DEBUG
                C = 'A';
                PrintDebugTors;
                PrintDebugTors2;
                pr( logFile, "]\n" );
#endif /* DEBUG */

                break;

    /*____________________________________________________________*/
            case PDBQ_BRANCH:
                if (found_first_res) {
                    sscanf(Rec_line[ i ],"%*s %d %*d", &nrestor );
                    tlist[ ntor ][ ATM1 ]= nrestor + true_ligand_atoms;
                } else {
                    sscanf(Rec_line[ i ],"%*s %d %*d", &tlist[ ntor ][ ATM1 ] );
                }
                tlist[ ntor ][ ATM2 ] = atomnumber[ i+1 ];
                --tlist[ ntor ][ ATM1 ];
                if ( tlist[ ntor ][ ATM2 ] == tlist[ ntor ][ ATM1 ]) {
                    prStr( error_message, "ERROR: line %d:\n%s\nThe two atoms defining torsion %d are the same!", (i+1), Rec_line[ i ], (ntor+1) );
                    stop( error_message );
                    exit( -1 );
                } /* endif */
                nbranches = 0;

#ifdef DEBUG
                C = 'B';
                PrintDebugTors;
                PrintDebugTors2;
                pr( logFile, "]\n" );
#endif /* DEBUG */

                for ( j = (i+2); j < nrecord; j++) {
                    for ( ii = 0; ii < 4; ii++ ) {
                        rec5[ii] = (char)tolower( (int)Rec_line[ j ][ ii ] );
                    }
                    rec5[4] = '\0';
#ifdef DEBUG
                    C = 'b';
                    PrintDebugTors;
                    PrintDebugTors2;
                    pr( logFile, "]\n" );
#endif /* DEBUG */

                    if (equal(rec5,"endb", 4) && (nbranches == 0))  break;
                    if (equal(rec5,"endb", 4) && (nbranches != 0))  --nbranches;
                    if (equal(rec5,"bran", 4))                      ++nbranches;
                    if (equal(rec5,"atom", 4) || equal(rec5,"heta", 4)) {
                        tlist[ ntor ][ tlist[ ntor ][ NUM_ATM_MOVED ] + 3 ] = atomnumber[ j ];
                        ++tlist[ ntor ][ NUM_ATM_MOVED ];
                    } /* endif */
                } /* j */
                ++ntor;

#ifdef DEBUG
                PrintDebugTors;
                PrintDebugTors2;
                pr( logFile, "]\n" );
#endif /* DEBUG */

                break;

    /*____________________________________________________________*/
            case PDBQ_TORS:

                tlist[ ntor ][ ATM2 ] = atomnumber[ i+1 ];
                tlist[ ntor ][ ATM1 ] = atomlast;
                if ( tlist[ ntor ][ ATM2 ] == tlist[ ntor ][ ATM1 ]) {
                    prStr( error_message, "ERROR: line %d:\n%s\nThe two atoms defining torsion %d are the same!", (i+1), Rec_line[ i ], (ntor+1) );
                    stop( error_message );
                    exit( -1 );
                }
                nbranches = 0;

#ifdef DEBUG
                C = 'T';
                PrintDebugTors;
                PrintDebugTors2;
                pr( logFile, "]\n" );
#endif /* DEBUG */

                for ( j = (i+2); j < nrecord; j++) {
                    for (ii = 0; ii < 4; ii++) {
                        rec5[ ii ] = (char)tolower( (int)Rec_line[ j ][ ii ] );
                    }
#ifdef DEBUG
                    C = 't';
                    PrintDebugTors;
                    PrintDebugTors2;
                    pr( logFile, "]\n" );
#endif /* DEBUG */

                    if (equal(rec5,"endb", 4) && (nbranches == 0))  break;
                    if (equal(rec5,"endb", 4) && (nbranches != 0))  --nbranches;
                    if (equal(rec5,"bran", 4))                      ++nbranches;
                    if (equal(rec5,"atom", 4) || equal(rec5,"heta", 4)) {
                        tlist[ ntor ][ tlist[ ntor ][ NUM_ATM_MOVED ] + 3 ] = atomnumber[ j ];
                        ++tlist[ ntor ][ NUM_ATM_MOVED ];
                    } /* endif */

#ifdef DEBUG
                    PrintDebugTors;
                    PrintDebugTors2;
                    pr( logFile, "]\n" );
#endif /* DEBUG */

                } /* j */
                ++ntor;

#ifdef DEBUG
                C = 't';
                PrintDebugTors;
                PrintDebugTors2;
                pr( logFile, "]\n" );
#endif /* DEBUG */

                break;

    /*____________________________________________________________*/
            case PDBQ_CONSTRAINT:

                sscanf(Rec_line[ i ],"%*s %d %d " FDFMT2, P_atomC1, P_atomC2, &lower, &upper);

                *P_B_constrain = TRUE;

                upper = fabs( (double)upper );
                lower = fabs( (double)lower );

                pr( logFile, "Constrain the distance between atom %d and atom %d to be within %.3f and %.3f Angstroms.\n\n", *P_atomC1, *P_atomC2, lower, upper);

                if (lower > upper) {
                    pr( logFile, "WARNING!  The lower bound was larger than the upper bound. I will switch these around.\n\n");
                    temp = upper;
                    upper = lower;
                    lower = temp;

                } else if (lower == upper) {
                    pr( logFile, "WARNING!  The lower bound is the same as the upper bound.\n\n");
                    upper += 0.01;
                }
                *P_sqlower = lower * lower;
                *P_squpper = upper * upper;
                break;

    /*____________________________________________________________*/
            case PDBQ_BEGIN_RES:
                found_new_res = 1;
                // if this is the first BEGIN_RES tag, then set the number of torsions in the ligand
                if (found_first_res == 0) {
                    *P_ntor_ligand = ntor;
                }
                found_first_res++;
                natoms_in_res = 0; /* reset number of atoms in this residue */
                break;

    /*____________________________________________________________*/
            case PDBQ_END_RES:
                found_new_res = 0;
                nres++;
                pr(logFile, "Residue number %d has %d moving atoms.\n\n", nres, natoms_in_res-2);
                break;

    /*____________________________________________________________*/
            case PDBQ_TORSDOF:
                sscanf(Rec_line[i], "%*s %d", P_ntorsdof);
                pr( logFile, "\nTORSDOF record detected: number of torsional degress of freedom has been set to %d.\n", *P_ntorsdof );
                break;

    /*____________________________________________________________*/
            default:
                break;
    /*____________________________________________________________*/
        } /* switch -- finished parsing this line of PDBQT file*/
    } /* i --- do next record in PDBQT file... */

    /*
    \   Sort Torsion list on number of atoms moved,
     \______________________________________________________________
    */

    if (ntor > MAX_TORS) {
        prStr( error_message, "ERROR: Too many torsions have been found (i.e. %d); maximum allowed is %d.\n Either: change the \"#define MAX_TORS\" line in constants.h\n Or:     edit \"%s\" to reduce the number of torsions defined.", (ntor+1), MAX_TORS, smFileName );
        stop( error_message );
        exit( -1 );
    } else {
        *P_ntor = ntor;
    }
    //if there are no flexible residues, still need to set P_ntor_ligand
    if (found_first_res == 0) {
         *P_ntor_ligand = ntor;
    }

#define check_atomnumber_ok( a )  (((a) >= 0) && ((a) < natoms))

    Boole B_atom_number_OK = TRUE;

    for (itor=0; itor<ntor; itor++ ) {
        B_atom_number_OK &= check_atomnumber_ok( tlist[ itor ][ ATM1 ] );
        B_atom_number_OK &= check_atomnumber_ok( tlist[ itor ][ ATM2 ] );
        if (B_atom_number_OK) for (int i=0;  i < tlist[ itor ][ NUM_ATM_MOVED ]; i++ ) {
                B_atom_number_OK &= check_atomnumber_ok( tlist[ itor ][ 3+i ] );
            }
        if (B_atom_number_OK) continue;
        pr( logFile, "%s: ERROR:  Torsion number %d between atom %d and atom %d has one or more atoms (out of %d atoms) that are out of range.\n\n", programname, itor+1, 1+tlist[itor][ATM1], 1+tlist[itor][ATM2], tlist[itor][NUM_ATM_MOVED] );
        pr( stderr, "%s: ERROR:  Torsion number %d between atom %d and atom %d has one or more atoms (out of %d atoms) that are out of range.\n\n", programname, itor+1, 1+tlist[itor][ATM1], 1+tlist[itor][ATM2], tlist[itor][NUM_ATM_MOVED] );
        exit(-1);
    }

    itor = 0;
    for (i=0;  i<MAX_ATOMS;  i++) {
        for ( j=0; j<ntor; j++ ) {
            if (tlist[ j ][ NUM_ATM_MOVED ] == i) {
                for (k=0;  k<MAX_ATOMS;  k++) {
                    tlistsort[ itor ][ k ] = tlist[ j ][ k ];
                }
                ++itor;
            }
        }
    }
    for ( i=0; i<ntor; i++ ) {
        for (j=0;  j<MAX_ATOMS;  j++) {
            tlist[ i ][ j ] = tlistsort[ i ][ j ];
        }
    }
    if (ntor > 0) {
        pr( logFile, "\n\nNumber of Rotatable Bonds in Small Molecule =\t%d torsions\n", ntor);
        pr( logFile, "\n\nTORSION TREE\n____________\n\nSorted in order of increasing number of atoms moved:\n\n" );
     
        pr( logFile, "Torsion                    #\n" );
        pr( logFile, " #  Atom1--Atom2 Moved List of Atoms Moved\n" );
        pr( logFile, "___ ____________ _____ ________________________________________________________\n");
        for ( j=0; j<ntor; j++ ) {
            pr( logFile, "%2d  %5s--%-5s  %3d  ", j+1, pdbaname[ tlist[ j ][ ATM1 ] ], pdbaname[ tlist[ j ][ ATM2 ] ], tlist[ j ][ NUM_ATM_MOVED ] );
            imax = tlist[ j ][ NUM_ATM_MOVED ] + 2;
            for ( i = 3; i <= imax; i++ ) {
                pr( logFile, "%s%c", pdbaname[ tlist[ j ][ i ]], (i<imax)?',':'.' );
            }
            pr( logFile, "\n" );
        }
        pr( logFile, "\n" );
    } else { 
        pr( logFile, "\n*** No Rotatable Bonds detected in Small Molecule. ***\n\n" );
    }
}
/* EOF */
