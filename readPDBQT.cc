/*

 $Id: readPDBQT.cc,v 1.1 2005/03/11 02:11:31 garrett Exp $

*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <sys/times.h>
#include <sys/param.h>
#include <time.h>
#include <ctype.h> /* tolower */
#include "readPDBQT.h"
#include "pdbqtokens.h"
#include <search.h>
#include "structs.h"

/*----------------------------------------------------------------------------*/

extern int    debug;
extern int    parse_tors_mode;
extern FILE  *logFile;
extern char  *programname;
extern int    true_ligand_atoms;
extern int    oldpdbq;

/*----------------------------------------------------------------------------*/
Molecule readPDBQT( char  pdbq_line[ LINE_LEN ],

              int   num_atm_maps,

              int   *P_natom,
              FloatOrDouble crdpdb[ MAX_ATOMS ][ NTRN ],
              FloatOrDouble charge[ MAX_ATOMS ],
              Boole *P_B_haveCharges,
              int   map_index[ MAX_ATOMS ],        // was:   int type[ MAX_ATOMS ]
              int   bond_index[ MAX_ATOMS ],
              char  pdbaname[ MAX_ATOMS ][ 5 ],
              char  pdbqFileName[ MAX_CHARS ],
              char  atomstuff[ MAX_ATOMS ][ MAX_CHARS ],
              int   Htype,

              Boole *P_B_constrain,
              int   *P_atomC1,
              int   *P_atomC2,
              FloatOrDouble *P_sqlower,
              FloatOrDouble *P_squpper,

              int   *P_ntor1,
              int   *P_ntor,
              int   tlist[ MAX_TORS ][ MAX_ATOMS ],
              FloatOrDouble vt[ MAX_TORS ][ NTRN ],

              int   *P_Nnb,
              int   nonbondlist[MAX_NONBONDS][MAX_NBDATA],

              Clock jobStart,
              struct tms tms_jobStart,
              char  hostnm[ MAX_CHARS ],
              int   *P_ntorsdof,
              int   outlev,
              int   ignore_inter[MAX_ATOMS],
              int   B_include_1_4_interactions)

{
    FILE *pdbqFile;
    static char  dummy[ LINE_LEN ];
    static char  error_message[ LINE_LEN ];
    static char  message[ LINE_LEN ];
    static char  pdbq_record[ MAX_RECORDS ][ LINE_LEN ];
    static char  rec5[ 5 ];

    FloatOrDouble aq = 0.;
    FloatOrDouble lq = 0.;
    FloatOrDouble total_charge = 0.;
    FloatOrDouble total_charge_ligand = 0.;
    FloatOrDouble total_charge_residues = 0.;
    FloatOrDouble uq = 0.;

    static int atomnumber[ MAX_RECORDS ];
    int   ii = 0;
    int   iq = 0;
    int   natom = 0;
    static int   nbmatrix[ MAX_ATOMS ][ MAX_ATOMS ];
    int   nrecord = 0;
    int   ntor = 0;
    static int   ntype[ MAX_ATOMS ];
    static int   piece[ MAX_ATOMS ];
    int   found_begin_res = 0;  // 0 means we have not yet found a BEGIN_RES record...
    int   keyword_id = -1;
    int   nres = 0;
    int   npiece = 0;

    Boole B_has_conect_records = FALSE;

    int bonded[MAX_ATOMS][6];

    register int   i = 0;
    register int   j = 0;

    static FloatOrDouble QTOL = 0.005;

    Molecule mol;

    ParameterEntry thisparm;

    for (i=0; i<MAX_RECORDS; i++) {
        atomnumber[i] = 0;
    }

    for (j = 0;  j < MAX_ATOMS;  j++ ) {
        ntype[ j ] = 0;
        piece[ j ] = 0;
    }

    /*
    **  Attempt to open the ligand PDBQT file...
    */
    sscanf( pdbq_line, "%*s %s", pdbqFileName );
    if ( openFile( pdbqFileName,"r",&pdbqFile,jobStart,tms_jobStart,TRUE )) {
        pr( logFile, "Atomic coordinate, partial charge, PDBQT file = \"%s\"\n\n", pdbqFileName );
    }

    /*
    **  Count the number of records in the PDBQT file first...
    */
    nrecord = 0;
    while( fgets( dummy, LINE_LEN, pdbqFile ) != NULL ) {
        ++nrecord;
    }
    (void)fclose( pdbqFile );
    
    if ( nrecord > MAX_RECORDS ) {
        prStr( error_message, "ERROR: %d records read in, but only dimensioned for %d.\nChange \"MAX_RECORDS\" in \"constants.h\".", nrecord, MAX_RECORDS);
        stop( error_message );
        exit( -1 );
    } else {
        /*
         * Read in the input PDBQT file...
         */
        if ( openFile( pdbqFileName,"r", &pdbqFile, jobStart, tms_jobStart, TRUE) ) {
            pr( logFile, "\nINPUT PDBQT FILE:" );
            pr( logFile, "\n________________\n\n\n" );
            for (i=0; i<nrecord; i++) {
                if ( fgets( pdbq_record[ i ], LINE_LEN, pdbqFile ) != NULL ) {
                    pr( logFile, "INPUT-PDBQT: %s", pdbq_record[ i ] );
                }
            } /* i */
            pr( logFile, UnderLine );
        } /*if*/
        (void)fclose( pdbqFile );
    } /*if*/

    /*  Read in the ATOMs and HETATMs, 
     *  count them; 
     *  store the (x,y,z) coordinates...
     *  
     *  Also, check for any BEGIN_RES records, for receptor flexibility...
     */
    natom = 0;
    for (i = 0;  i < nrecord;  i++) {       // loop over all the lines in the ligand file
        strncpy( pdbq_line, pdbq_record[i], (size_t) LINE_LEN );
        // parse this ligand file's line
        keyword_id = parse_pdbq_line(pdbq_line);

        if ((keyword_id == PDBQ_ATOM) || (keyword_id == PDBQ_HETATM)) { // an atom line

            ENTRY item, *found_item; /*see hsearch(3C)*/
            ParameterEntry * found_parm;

            atomnumber[i] = natom;  // set the serial atomnumber[i] for this i-th atom

            // set up piece array
            // by reading in the records of the PDBQT file;  
            // each "piece" is a self-contained rigid entity.
            piece[natom] = npiece;
             
            // get the atom's coordinates, charge, and parameters
            readPDBQTLine( pdbq_line, crdpdb[natom], &charge[natom], &thisparm ); // sets "autogrid_type" in thisparm

            mol.crdpdb[natom][X] = crdpdb[natom][X];
            mol.crdpdb[natom][Y] = crdpdb[natom][Y];
            mol.crdpdb[natom][Z] = crdpdb[natom][Z];
            mol.crd[natom][X]    = crdpdb[natom][X];
            mol.crd[natom][Y]    = crdpdb[natom][Y];
            mol.crd[natom][Z]    = crdpdb[natom][Z];

            if (found_begin_res) {
                total_charge_residues += charge[natom];
            } else {
                total_charge_ligand += charge[natom];
            }
            *P_B_haveCharges = TRUE;

            strncpy( atomstuff[natom], pdbq_line, (size_t)30 );
            atomstuff[natom][30] = '\0';
            strcpy(mol.atomstr[natom], atomstuff[natom]);

            sscanf( &pdbq_line[ 12 ], "%s", pdbaname[natom] );
            map_index[natom] = -1;

            // "map_index" is used as an index into the "map" array,
            // to look up the correct energies in the current grid cell,
            // thus:  map[][][][map_index[natom]]

            // AutoGrid 4 Typing --/
            item.key = thisparm.autogrid_type;
            if ((found_item = hsearch (item, FIND)) != NULL) {
                found_parm = (ParameterEntry *)found_item->data;
                map_index[natom] = found_parm->map_index;
                bond_index[natom] = found_parm->bond_index;
                if (outlev > 0) {
                    (void) fprintf( logFile, "Found parameters for ligand atom %d, atom type \"%s\", grid map index = %d\n",
                            natom, found_parm->autogrid_type, found_parm->map_index );
                }
            } else {
                // We could not find this parameter -- return error here
                 prStr( message, "\n%s: *** WARNING!  Unknown ligand atom type \"%s\" found.  You should add  parameters for it to the parameter library first! ***\n\n", programname, total_charge_ligand);
                 pr_2x( stderr, logFile, message );
            }
            // /--AutoGrid4 
            
            if (map_index[natom] == -1) {
                pr( logFile, "%s: WARNING:  atom type not found, using the default, atom type = 1, instead.\n", programname);
                map_index[natom] = 0;  // we are 0-based internally, 1-based in printed output, for human's sake
            }

            ++ntype[ map_index[natom] ]; /* increment the number of atoms with this atomtype */
            ++natom; /* increment the number of atoms in PDBQT file */

        } else {
            ++npiece;
        }

        if (!found_begin_res) {
            /* No BEGIN_RES found yet. */
            /*
             * Keep updating "true_ligand_atoms" until we find a "BEGIN_RES",
             */
            true_ligand_atoms = natom; /* Set number of atoms in ligand now. */
            if (keyword_id == PDBQ_BEGIN_RES) {
                /* then a flexible receptor sidechain was found in the PDBQ file. */
                found_begin_res = 1; /* flag that we've found a BEGIN_RES record. */
                pr( logFile, "\nNumber of atoms in movable ligand = %d\n\n", true_ligand_atoms );
            }
        }

        if (keyword_id == PDBQ_BEGIN_RES) {
            nres++;  // increment the number of residues
        }

        if (keyword_id == PDBQ_CONECT) {
            // At least some of the atoms in the "ligand" may have their connectivity specified
            // so we could set up their bonded entries.  Future versions...
            B_has_conect_records = TRUE;
        }
        
    } /* i, next record in PDBQT file */

    pr( logFile, "\nNumber of atoms found in flexible receptor sidechains =\t%d atoms\n\n", natom - true_ligand_atoms);

    pr( logFile, "Total number of atoms found in PDBQT file =\t%d atoms\n\n", natom);

    pr( logFile, "Number of flexible receptor sidechains found =\t%d residues\n\n", nres);

    if (natom > MAX_ATOMS) {
        prStr( error_message, "ERROR: Too many atoms found (i.e. %d); maximum allowed is %d.\nChange the \"#define MAX_ATOMS\" line in \"constants.h\"\n.", natom, MAX_ATOMS );
        stop( error_message );
        exit( -1 );
    } else {
        *P_natom = natom;
        mol.natom = natom;
    }

    pr( logFile, "\nSummary of number of atoms of a given atom type:\n");
    pr( logFile, "------------------------------------------------\n\n");
    for (i=0; i<num_atm_maps; i++) {
        pr( logFile, "Number of atoms with atom type %d = %2d\n", i+1, ntype[ i ]);
    }

    pr( logFile, "\n\nSummary of total charge on ligand, residues and overall:\n");
    pr( logFile, "-----------------------------------------------------\n");
    /*
     * Check total charge on ligand
     */
    pr( logFile, "\nTotal charge on ligand =\t\t\t\t%+.3f e\n", total_charge_ligand );
    iq = (int) ( ( aq = fabs( total_charge_ligand ) ) + 0.5 );
    lq = iq - QTOL;
    uq = iq + QTOL;
    if ( ! ((aq >= lq) && (aq <= uq)) ) {
        prStr( message, "\n%s: *** WARNING!  Non-integral total charge (%.3f e) on ligand! ***\n\n", programname, total_charge_ligand);
        pr_2x( stderr, logFile, message );
    }

    /*
     * Check total charge on residues
     */
    pr( logFile, "\nTotal charge on residues =\t\t\t\t%+.3f e\n\n", total_charge_residues );
    iq = (int) ( ( aq = fabs( total_charge_residues ) ) + 0.5 );
    lq = iq - QTOL;
    uq = iq + QTOL;
    if ( ! ((aq >= lq) && (aq <= uq)) ) {
        prStr( message, "\n%s: *** WARNING!  Non-integral total charge (%.3f e) on residues! ***\n\n", programname, total_charge_residues);
        pr_2x( stderr, logFile, message );
    }

    /*
     * Check total charge on all PDBQT atoms
     */
    total_charge = total_charge_ligand + total_charge_residues;
    pr( logFile, "Total charge on all moving atoms (ligand + residues) =\t\t%+.3f e\n\n\n", total_charge);
    iq = (int) ( ( aq = fabs( total_charge) ) + 0.5 );
    lq = iq - QTOL;
    uq = iq + QTOL;
    if ( ! ((aq >= lq) && (aq <= uq)) ) {
        prStr( message, "\n%s: *** WARNING!  Non-integral total charge (%.3f e) on all moving atoms! ***\n\n", programname, total_charge);
        pr_2x( stderr, logFile, message );
    }

    /* 
     *  Work out where the torsions are; and what they move...
     *
     *  Also, detect which atoms we should ignore in the intermolecular energy
     *  calculation (ignore_inter[MAX_ATOMS] array)
     */
    mkTorTree( atomnumber, pdbq_record, nrecord, 
                tlist, &ntor, 
                pdbqFileName, pdbaname, 
                P_B_constrain, P_atomC1, P_atomC2,
                P_sqlower, P_squpper, P_ntorsdof,
                ignore_inter);

    *P_ntor  = ntor;
    *P_ntor1 = ntor - 1;

    mol.S.ntor = ntor;

    if (ntor > 0) {
        /*
        **  Create list of internal non-bonds to be used in internal energy calculation...
        */
        if (debug>0) {
            pr( logFile, "Finding bonds.\n\n");
        }

        // initialise the bonded array
        for (i=0; i<natom; i++) {
            for (j=0; j<5; j++) {
                bonded[i][j] = -1;
            } // j
            bonded[i][5] = 0;
        } // i

        if (debug>0) { 
            printbonds(natom, bonded, "\nDEBUG:  1. BEFORE getbonds, bonded[][] array is:\n\n", 1); 
        }

        getbonds( crdpdb, natom, bond_index, bonded);

        if (debug>0) {
            printbonds(natom, bonded, "\nDEBUG:  2. AFTER getbonds, bonded[][] array is:\n\n", 0);
            pr( logFile, "Detecting all non-bonds.\n\n");
        }

        nonbonds( crdpdb, nbmatrix, natom, bond_index, B_include_1_4_interactions, bonded);

        if (debug>0) {
            printbonds(natom, bonded, "\nDEBUG:  4. AFTER nonbonds, bonded[][] array is:\n\n", 0);
            pr( logFile, "Weeding out non-bonds in rigid parts of the torsion tree.\n\n");
        }

        weedbonds( natom, pdbaname, piece, ntor, tlist, nbmatrix, P_Nnb, nonbondlist, outlev, map_index );

        if (debug>0) {
            pr( logFile, "Calculating unit vectors for each torsion.\n\n");
        }

        torNorVec( crdpdb, ntor, tlist, vt );

        for ( i=0; i<MAX_TORS; i++) {
            mol.vt[i][X] = vt[i][X];
            mol.vt[i][Y] = vt[i][Y];
            mol.vt[i][Z] = vt[i][Z];
            for ( j=0; j<MAX_ATOMS; j++) {
                mol.tlist[i][j] = tlist[i][j];
            }
        }
    } else {
        fprintf( logFile, ">>> No torsions detected, so skipping \"nonbonds\", \"weedbonds\" and \"torNorVec\" <<<\n\n");
    }

    /*
    **  End program if just parsing torsions...
    */
    if (parse_tors_mode) {
        prStr( message, "\n\n *** PARSE TORSIONS MODE - Stopping here ***\n\n");
        success( hostnm, jobStart, tms_jobStart );
        exit( 0 );
    }

    return mol;
}


/*----------------------------------------------------------------------------*/
/* readPDBQTLine.cc */

void readPDBQTLine( char line[LINE_LEN],
                    FloatOrDouble crd[SPACE],
                    FloatOrDouble *ptr_q,
                    ParameterEntry *thisparm)

/*----------------------------------------------------------------------------*/
{
    char str[4][WORDLEN];

    // X, Y, Z coordinates
    (void) sscanf( &line[30], "%s %s %s", str[X], str[Y], str[Z] );

    crd[X] = atof( str[X] );
    crd[Y] = atof( str[Y] );
    crd[Z] = atof( str[Z] );

#ifdef DEBUG
    (void) fprintf(stderr, "readPDBQTLine: oldpdbq=%d\n%s", oldpdbq, line );
#endif /* DEBUG */

    if (oldpdbq) {
        (void) sscanf( &line[54], "%s", str[3] );
        // partial charge, q
        *ptr_q = atof( str[3] );
    } else {
        (void) sscanf( &line[70], "%s", str[3] );
        // partial charge, q
        *ptr_q = atof( str[3] );
        // atom type name
        (void) sscanf(&line[77], "%s", thisparm->autogrid_type );
    }

#ifdef DEBUG
    fprintf(stderr, "readPDBQTLine:  %.3f, %.3f, %.3f, %.3f\n", crd[X], crd[Y], crd[Z], *ptr_q );
#endif /* DEBUG */
}

/*----------------------------------------------------------------------------*/
/* EOF */
