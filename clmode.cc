/* clmode.cc */


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

void  clmode( char  atm_typ_str[ATOM_MAPS],
              int   num_atm_maps,
              float clus_rms_tol,
              char  hostnm[MAX_CHARS],
              Clock jobStart,
              struct tms tms_jobStart,
              Boole B_write_all_clusmem,
              char  clusFN[MAX_CHARS],
              float crdpdb[MAX_ATOMS][SPACE],
              float sml_center[SPACE],
              Boole B_symmetry_flag,
              char  rms_ref_crds[MAX_CHARS] )

{
    FILE *clusFile;
    register int xyz = 0;
    int   anum = 0;
    char  atomstuff[MAX_ATOMS][MAX_CHARS];
    float crdSave[MAX_RUNS][MAX_ATOMS][SPACE];
    float econf[MAX_RUNS];
    float eSave[2];
    Boole B_haveAtoms = FALSE;
    Boole B_haveTypes = FALSE;
    int   ii = 0;
    int   lastanum = -1;
    char  line[LINE_LEN];
    int   nat = 0;
    int   natom = 0;
    int   natom_1 = -1;
    int   nconf = 0;
    int   ntype[MAX_ATOMS];
    char  pdbaname[MAX_ATOMS][5];
    float q = 0.;
    char  rec5[5];
    int   ss = 0;
    char  anumStr[5];
    int   type[MAX_ATOMS];
    float clu_rms[MAX_RUNS][MAX_RUNS];
    int   cluster[MAX_RUNS][MAX_RUNS];
    register int i = 0;
    register int j = 0;
    int   irunmax = -1;
    int   isort[MAX_RUNS];
    int   ncluster = 0;
    int   num_in_clu[MAX_RUNS];
    float ref_crds[MAX_ATOMS][SPACE];
    int   ref_natoms = -1;
    float ref_rms[MAX_RUNS];

    for (j = 0; j < MAX_RUNS; j++) {
        num_in_clu[j] = 0;
        isort[j] = j;
        econf[j] = 0.;
    }
/*
**  Open file containing coordinates to be clustered...
*/
    if ( openFile( clusFN , "r", &clusFile, jobStart, tms_jobStart, TRUE ) ) {
        pr( logFile, "Conformations to be clustered are in this file: \"%s\"\n\n", clusFN );
    }
/*
    Read in the conformations
    All we need are the xyz's of each conformation,
    and their Energies, plus the Run number/parent dlg file.
*/
    while ( fgets( line, LINE_LEN, clusFile) != NULL ) {

        pr( logFile, "INPUT-PDBQ: %s", line);

        for (ii = 0; ii < 4; ii++) { 
            rec5[ii] = tolower( (int)line[ii] );
        }
        if (( strindex( line, "USER    Total Interaction Energy of Complex") >= 0 )
         || ( strindex( line, "REMARK  Total Interaction Energy of Complex") >= 0 )) {
/*
            Read in the energy of this conformation
*/
            if ( B_haveAtoms ) {
                econf[nconf] = 0.;
                sscanf( line, "%*s %*s %*s %*s %*s %*s %*s %f", &econf[nconf]);
            } else {
                eSave[ss]=0.;
                sscanf( line, "%*s %*s %*s %*s %*s %*s %*s %f", &eSave[ss] );
                ++ss;
            }

        } else if (( strindex( line, "USER    Final Docked Energy") >= 0 ) 
                || ( strindex( line, "REMARK  Final Docked Energy") >= 0 )) {
/*
            Read in the energy of this conformation
*/
            if ( B_haveAtoms ) {
                econf[nconf] = 0.;
                sscanf( line, "%*s %*s %*s %*s %*s %f", &econf[nconf]);
            } else {
                eSave[ss]=0.;
                sscanf( line, "%*s %*s %*s %*s %*s %f", &eSave[ss] );
                ++ss;
            }

        } else if (equal( rec5,"atom", 4) || equal( rec5,"heta", 4)) {

            readPDBQLine( line, crdSave[nconf][nat], &q );

            if ( ! B_haveAtoms ) {
                sscanf( &line[6], "%s", anumStr );
                if ( ( anum = atoi( anumStr )) < lastanum ) {
/*
                    Start of next conformation,
*/
                    B_haveAtoms = TRUE;
                    B_haveTypes = TRUE;
                    ++nconf;
                    for (xyz = 0;  xyz < SPACE;  xyz++) {
                        crdSave[nconf][0][xyz] = crdSave[nconf-1][nat][xyz]; 
                    }
                    econf[0] = eSave[0];
                    econf[1] = eSave[1];
                    natom = nat;
                    natom_1 = natom-1;
                    nat = 0;
                } else {
                    strncpy( atomstuff[nat], line, (size_t)30 );
                    atomstuff[nat][30] = '\0';
                    if ( ! B_haveTypes ) {
                        type[nat] = -1;
                        sscanf( &line[12], "%s", pdbaname[nat] );
                        type[nat] = get_atom_type( pdbaname[nat], atm_typ_str );
                        if ( type[nat] == -1 ) {
                            pr( logFile, "\nNOTE: Atom number %d, using default atom type 1...\n\n", nat+1);
                            type[nat] = 1;
                        } else {
                            pr( logFile, "\nAtom number %d, recognized atom type = %d...\n\n", nat+1, type[nat]+1);
                        }
                        ++ntype[ type[nat] ];
                    }
                }
                lastanum = anum;
            } else if ( nat == natom_1 ) {
/*
                Increment total number of conformations,
*/
                ++nconf;
                nat = -1; /*  Pre-zero out the "nat" counter... */
            }
            ++nat;
        } 
    } /* end while */

    irunmax = nconf;

    pr( logFile, "\nNumber of conformations found = %d\n", nconf );

    pr( logFile, "\nNumber of atoms per conformation = %d\n\n", natom );

    for (i=0; i<num_atm_maps; i++) {
        pr( logFile, "Number of atoms with type %d = %d\n", i+1, ntype[i]);
    }

    if (strncmp(rms_ref_crds,"unspecified filename",20) != 0) {
/*
        Read in reference structure, specified by the "rmsref" command...
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
                                     sml_center, clu_rms, B_symmetry_flag,
                                     ref_crds, ref_natoms, ref_rms);

        pr( logFile, "\nOutputting structurally similar clusters, ranked in order of increasing energy.\n" );
        flushLog;

        prClusterHist( ncluster, irunmax, clus_rms_tol, num_in_clu, 
                       cluster, econf, clu_rms, ref_rms);

        bestpdb( ncluster, num_in_clu, cluster, econf, crdSave, 
                 atomstuff, natom, B_write_all_clusmem, ref_rms);

    }/*if we have more than 1 conformation... */

/*
**  End cluster_mode and END PROGRAM...
*/
    success( hostnm, jobStart, tms_jobStart );

    exit((int)0);
}
/* EOF */
