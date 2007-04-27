/*

 $Id: clmode.cc,v 1.7 2007/04/27 06:01:48 garrett Exp $

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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <sys/types.h>
#include <sys/times.h>
#include "clmode.h"


extern FILE *logFile;
extern char *programname;

void  clmode( int   num_atm_maps,
              Real clus_rms_tol,
              char  hostnm[MAX_CHARS],
              Clock jobStart,
              struct tms tms_jobStart,
              Boole write_all_clusmem,
              char  clusFN[MAX_CHARS],
              Real crdpdb[MAX_ATOMS][SPACE],
              Real sml_center[SPACE],
              Boole symmetry_flag,
              char  rms_ref_crds[MAX_CHARS] )

{
    FILE *clusFile;
    register int xyz = 0;
    int   anum = 0;
    char  atomstuff[MAX_ATOMS][MAX_CHARS];
    Real crdSave[MAX_RUNS][MAX_ATOMS][SPACE];
    Real econf[MAX_RUNS];
    Real eSave[2];
    Boole haveAtoms = FALSE;
    Boole haveTypes = FALSE;
    int   ii = 0;
    int   lastanum = -1;
    char  line[LINE_LEN];
    int   atomCounter = 0;
    int   natom = 0;
    int   natom_1 = -1;
    int   nconf = 0;
    int   confCounter = 0;
    int   ntype[MAX_ATOMS];
    char  pdbaname[MAX_ATOMS][5];
    Real q = 0.;
    char  rec5[5];
    int   nsaved = 0;
    char  anumStr[5];
    int   type[MAX_ATOMS];
    Real clu_rms[MAX_RUNS][MAX_RUNS];
    int   cluster[MAX_RUNS][MAX_RUNS];
    register int i = 0;
    register int j = 0;
    int   irunmax = -1;
    int   isort[MAX_RUNS];
    int   ncluster = 0;
    int   num_in_clu[MAX_RUNS];
    Real ref_crds[MAX_ATOMS][SPACE];
    int   ref_natoms = -1;
    Real ref_rms[MAX_RUNS];
    Boole haveEnergy = FALSE;
    ParameterEntry thisparm;

    for (j = 0; j < MAX_RUNS; j++) {
        num_in_clu[j] = 0;
        isort[j] = j;
        econf[j] = 0.;
    }
    /*
     * Open file containing coordinates to be clustered...
     */
    if ( openFile( clusFN , "r", &clusFile, jobStart, tms_jobStart, TRUE ) ) {
        pr( logFile, "Conformations to be clustered are in this file: \"%s\"\n\n", clusFN );
    }
    /*
     * Read in the conformations
     * All we need are the xyz's of each conformation,
     * and their Energies, plus the Run number/parent dlg file.
     */
    while ( fgets( line, LINE_LEN, clusFile) != NULL ) {

        pr( logFile, "INPUT-PDBQ: %s", line);

        for (ii = 0; ii < 4; ii++) { rec5[ii] = tolower( (int)line[ii] ); };

        if (( strindex( line, "USER    Total Interaction Energy of Complex") >= 0 )
         || ( strindex( line, "REMARK  Total Interaction Energy of Complex") >= 0 )) {
            /*
             * Read in the energy of this conformation;
             * This is preferred over "Final Docked Energy" because this is never
             * printed out rounded up as +4.25e+03, but always as +4246.45, e.g.:
             */
            if ( haveAtoms ) {
                econf[confCounter] = 0.;
                sscanf( line, "%*s %*s %*s %*s %*s %*s %*s " FDFMT, &econf[confCounter]);
                haveEnergy = TRUE;
            } else {
                /* ! haveAtoms
                 * We have not seen any atoms yet, so save this energy. 
                 */
                eSave[nsaved]=0.;
                sscanf( line, "%*s %*s %*s %*s %*s %*s %*s " FDFMT, &eSave[nsaved] );
                ++nsaved;
            }

        } else if ( (( strindex( line, "USER    Final Docked Energy") >= 0 ) 
                  || ( strindex( line, "REMARK  Final Docked Energy") >= 0 )) && ( ! haveEnergy ) ) {
            /*
             * Read in the energy of this conformation if we don't already
             * have an energy:
             */
            if ( haveAtoms ) {
                econf[confCounter] = 0.;
                sscanf( line, "%*s %*s %*s %*s %*s " FDFMT, &econf[confCounter]);
                haveEnergy = TRUE;
            } else {
                /* ! haveAtoms
                 * We have not seen any atoms yet, so save this energy. 
                 */
                eSave[nsaved]=0.;
                sscanf( line, "%*s %*s %*s %*s %*s " FDFMT, &eSave[nsaved] );
                ++nsaved;
            }

        } else if (equal( rec5,"atom", 4) || equal( rec5,"heta", 4)) {

            int serial;

            /* 
             * This line should contain coordinates, partial charge & 
             * atom type for one atom.
             * Let's save the coordinates for this atom, atomCounter.
             */
            readPDBQTLine( line, &serial, crdSave[confCounter][atomCounter], &q, &thisparm );

            if ( ! haveAtoms ) {
                /*
                 * We do not have any atoms for this conformation,
                 */
                sscanf( &line[6], "%s", anumStr );

                if ((anum = atoi(anumStr)) < lastanum) { /* initially, lastanum is -1, while anum is probably never -1, so this is false initially */
                    /*
                     * haveAtoms is FALSE, but this line begins with "atom" or
                     * "heta", so this must be the...
                     *
                     * Start of next conformation,
                     */
                    /* This is an atom line, so haveAtoms must be set to true: */
                    haveAtoms = TRUE;
                    /* We must also have read in the atom types: */
                    haveTypes = TRUE;
                    /* Transfer the saved energies to the econf arrays: */
                    econf[0] = eSave[0];
                    econf[1] = eSave[1];
                    /* Now we have the energy: */
                    haveEnergy = TRUE;
                    /* Increment the number of conformations */
                    ++confCounter;
                    for (xyz = 0;  xyz < SPACE;  xyz++) {
                        crdSave[confCounter][0][xyz] = crdSave[confCounter-1][atomCounter][xyz]; 
                    }
                    natom = atomCounter; /* number of atoms is set to the atom counter*/
                    natom_1 = natom - 1;  /* number of atoms minus 1, for 0-based counting */
                    atomCounter = 0; /* reset the counter "atomCounter" */
                } else {
                    /* 
                     * First of all, determine the atom types for all the atoms in the
                     * molecule: 
                     */
                    strncpy( atomstuff[atomCounter], line, (size_t)30 );
                    atomstuff[atomCounter][30] = '\0';
                    if ( ! haveTypes ) {
                        type[atomCounter] = -1;
                        sscanf( &line[12], "%s", pdbaname[atomCounter] );
                        /*
                         * Determine this atom's atom type:
                         */
                        type[atomCounter] = get_atom_type(pdbaname[atomCounter]);

                        if (type[atomCounter] == -1) {
                            pr( logFile, "\nNOTE: Atom number %d, using default atom type 1...\n\n", atomCounter+1);
                            type[atomCounter] = 1;
                        } else {
                            pr( logFile, "\nAtom number %d, recognized atom type = %d...\n\n", atomCounter+1, type[atomCounter]+1);
                        }
                        /* 
                         * Increment the number of atoms with this atom type:
                         */
                        ++ntype[ type[atomCounter] ];
                    }
                }
                /*
                 * Update the value of the last atom's serial number:
                 */
                lastanum = anum;
            } else if ( atomCounter == natom_1 ) {  /* initially, atomCounter=0, and natom_1= -1, so this is not true initially */
                /*
                 * We have all the atoms we expect for one molecule:
                 * Increment total number of conformations,
                 */
                ++confCounter;
                atomCounter = -1; /*  Pre-zero out the "atomCounter" counter... */
                haveEnergy = FALSE; /* we don't have energy yet for next conf. */
            }
            /* 
             * Just increment the number of atoms, atomCounter:
             */
            ++atomCounter;
            
        } /* This was an "atom" or "heta" line */
    } /* end while there is a new line. */

    irunmax = confCounter;
    nconf = confCounter;

    pr( logFile, "\nNumber of conformations found = %d\n", nconf );

    pr( logFile, "\nNumber of atoms per conformation = %d\n\n", natom );

    for (i=0; i<num_atm_maps; i++) {
        pr( logFile, "Number of atoms with type %d = %d\n", i+1, ntype[i]);
    }

    if (strncmp(rms_ref_crds,"unspecified filename",20) != 0) {
        /*
         * Read in reference structure, specified by the "rmsref" command...
         */
        if ((ref_natoms = getpdbcrds( rms_ref_crds, ref_crds)) == -1) {
     
            fprintf( logFile, "%s: Problems while reading \"%s\".\n", programname, rms_ref_crds);
            fprintf( logFile, "Will attempt to use the input PDBQ file coordinates as reference.\n");
     
        } else if (ref_natoms != natom) {
     
            pr( logFile, "%s: ERROR!  Wrong number of atoms in reference structure.\n", programname);
            pr( logFile, "Input PDBQ structure has %d atoms, but reference structure has %d atoms.\n\n", natom, ref_natoms);
            ref_natoms = -1;
     
        }
    }

    if (nconf <= 1) {

        pr( logFile, "\nSorry!  Unable to perform cluster analysis, because not enough structures were read in.\n");

    } else {

#ifdef DEBUG
     for (i = 0;  i < nconf;  i++) { 
        pr( logFile, "i=%-3d\tisort[i]=%-3d\teconf[isort[i]]=%+7.2f\n", 
                      i, isort[i], econf[isort[i]] );
    }
#endif /* DEBUG */

        pr( logFile, "\nSorting %d conformations by their energy.\n", irunmax);
        flushLog;

        sort_enrg( econf, isort, nconf );

#ifdef DEBUG
     for (i = 0;  i < nconf;  i++) { pr( logFile, "i=%-3d\tisort[i]=%-3d\teconf[isort[i]]=%+7.2f\n", i, isort[i], econf[isort[i]] ); }
#endif /* DEBUG */

        pr( logFile, "\nPerforming cluster analysis, using a cluster RMS tolerance of %.1f\n", clus_rms_tol );
        flushLog;

        ncluster = cluster_analysis( clus_rms_tol, cluster, num_in_clu, isort, 
                                     nconf, natom, type, crdSave, crdpdb, 
                                     sml_center, clu_rms, symmetry_flag,
                                     ref_crds, ref_natoms, ref_rms);

        pr( logFile, "\nOutputting structurally similar clusters, ranked in order of increasing energy.\n" );
        flushLog;

        prClusterHist( ncluster, irunmax, clus_rms_tol, num_in_clu, 
                       cluster, econf, clu_rms, ref_rms);

        bestpdb( ncluster, num_in_clu, cluster, econf, crdSave, 
                 atomstuff, natom, write_all_clusmem, ref_rms);

    }/*if we have more than 1 conformation... */

/*
 *  End cluster_mode and END PROGRAM...
 */
    success( hostnm, jobStart, tms_jobStart );

    exit((int)0);
}
/* EOF */
