/*

 $Id: readPDBQ.cc,v 1.6 2004/11/16 23:42:53 garrett Exp $

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
#include "readPDBQ.h"
#include "pdbqtokens.h"

/*----------------------------------------------------------------------------*/

extern int    debug;
extern int    parse_tors_mode;
extern FILE  *logFile;
extern char  *programname;
extern int    true_ligand_atoms;
extern int    oldpdbq;

/*----------------------------------------------------------------------------*/
Molecule readPDBQ( char  thisline[ LINE_LEN ],

              char  atm_typ_str[ ATOM_MAPS ],
              int   num_atm_maps,

              int   *P_natom,
              FloatOrDouble crdpdb[ MAX_ATOMS ][ NTRN ],
              FloatOrDouble charge[ MAX_ATOMS ],
              Boole *P_B_haveCharges,
              int   atmType[ MAX_ATOMS ],
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
              int   Nnbonds[ MAX_ATOMS ],
              int   nonbondlist[MAX_NONBONDS][4],

              Clock jobStart,
              struct tms tms_jobStart,
              char  hostnm[ MAX_CHARS ],
              int   *P_ntorsdof,
              int   outlev,
              int   ignore_inter[MAX_ATOMS])

{
    FILE *pdbqFile;
    static char  dummy[ LINE_LEN ];
    static char  error_message[ LINE_LEN ];
    static char  message[ LINE_LEN ];
    static char  record[ MAX_RECORDS ][ LINE_LEN ];
    static char  rec5[ 5 ];

    FloatOrDouble aq = 0.;
    FloatOrDouble lq = 0.;
    FloatOrDouble total_charge = 0.;
    FloatOrDouble total_charge_ligand = 0.;
    FloatOrDouble total_charge_residues = 0.;
    FloatOrDouble uq = 0.;

    static int Rec_atomnumber[ MAX_RECORDS ];
    int   ii = 0;
    int   iq = 0;
    int   iatom = 0;
    static int   nbmatrix_binary[ MAX_ATOMS ][ MAX_ATOMS ];
    int   nrecord = 0;
    int   ntor = 0;
    static int   ntype[ MAX_ATOMS ];
    static int   piece[ MAX_ATOMS ];
    int   found_begin_res = 0;
    int   keyword_id = -1;
    int   nres = 0;

    register int   i = 0;
    register int   j = 0;

    static FloatOrDouble QTOL = 0.005;

    Molecule mol;

    for (i=0; i<MAX_RECORDS; i++) {
        Rec_atomnumber[i] = 0;
    }

    for (j = 0;  j < MAX_ATOMS;  j++ ) {
        ntype[ j ] = 0;
        piece[ j ] = 0;
    }

    /*
    **  Attempt to open the ligand PDBQ file...
    */
    sscanf( thisline, "%*s %s", pdbqFileName );
    if ( openFile( pdbqFileName,"r",&pdbqFile,jobStart,tms_jobStart,TRUE )) {
        pr( logFile, "Atomic coordinate, partial charge, PDBQ file = \"%s\"\n\n", pdbqFileName );
    }

    /*
    **  Count the number of records in the PDBQ file first...
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
         * Read in the input PDBQ file...
         */
        if ( openFile( pdbqFileName,"r",&pdbqFile,jobStart,tms_jobStart,TRUE )) {
            pr( logFile, "\nINPUT PDBQ FILE:" );
            pr( logFile, "\n________________\n\n\n" );
            for (i=0; i<nrecord; i++) {
                if ( fgets( record[ i ], LINE_LEN, pdbqFile ) != NULL ) {
                    pr( logFile, "INPUT-PDBQ: %s", record[ i ] );
                }
            } /* i */
            pr( logFile, UnderLine );
        } /*if*/
        (void)fclose( pdbqFile );
    } /*if*/

    /*
     *  Read in the ATOMs and HETATMs, count them; store the (x,y,z) coordinates...
     *  
     *  Also, check for any BEGIN_RES records, for receptor flexibility...
     */
    iatom = 0;
    for (i = 0;  i < nrecord;  i++) {
        strncpy( thisline, record[ i ], (size_t)LINE_LEN );
        keyword_id = parse_pdbq_line(thisline);

        if ((keyword_id == PDBQ_ATOM) || (keyword_id == PDBQ_HETATM)) {
            Rec_atomnumber[ i ] = iatom;
             
            readPDBQLine( thisline, crdpdb[ iatom ], &charge[ iatom ] );
            mol.crdpdb[iatom][X] = crdpdb[iatom][X];
            mol.crdpdb[iatom][Y] = crdpdb[iatom][Y];
            mol.crdpdb[iatom][Z] = crdpdb[iatom][Z];
            mol.crd[iatom][X]    = crdpdb[iatom][X];
            mol.crd[iatom][Y]    = crdpdb[iatom][Y];
            mol.crd[iatom][Z]    = crdpdb[iatom][Z];
            if (found_begin_res) {
                total_charge_residues += charge[ iatom ];
            } else {
                total_charge_ligand += charge[ iatom ];
            }
            *P_B_haveCharges = TRUE;

            strncpy( atomstuff[ iatom ], thisline, (size_t)30 );
            atomstuff[iatom][30] = '\0';
            strcpy(mol.atomstr[iatom], atomstuff[iatom]);

            sscanf( &thisline[ 12 ], "%s", pdbaname[ iatom ] );
            atmType[ iatom ] = -1;
            atmType[ iatom ] = get_atom_type( pdbaname[ iatom ], atm_typ_str );
            if (atmType[ iatom ] == -1) {
                pr( logFile, "%s: atom type error, using the default, atom type = 1\n", programname);
                atmType[ iatom ] = 1;
            }
            ++ntype[ atmType[ iatom ] ]; /* count the number of atoms with this atomtype */
            ++iatom; /* count the number of atoms in PDBQ file */

        } else if (!found_begin_res) {
            /* No BEGIN_RES found yet. */
            /*
             * Keep updating "true_ligand_atoms" until we find a "BEGIN_RES",
             */
            true_ligand_atoms = iatom; /* Set number of atoms in ligand now. */
            if (keyword_id == PDBQ_BEGIN_RES) {
                /* then a flexible receptor sidechain was found in the PDBQ file. */
                found_begin_res = 1; /* flag that we've found a BEGIN_RES record. */
                pr( logFile, "\nNumber of atoms in movable ligand = %d\n\n", true_ligand_atoms );
            }
        }
        if (keyword_id == PDBQ_BEGIN_RES) {
            nres++;
        }
    } /* i, next record in PDBQ file */

    pr( logFile, "Number of atoms found in flexible receptor sidechains =\t%d atoms\n\n", iatom - true_ligand_atoms);

    pr( logFile, "Total number of atoms found in PDBQ file =\t%d atoms\n\n", iatom);

    pr( logFile, "Number of flexible receptor sidechains found =\t%d residues\n\n", nres);

    if (iatom > MAX_ATOMS) {
        prStr( error_message, "ERROR: Too many atoms found (i.e. %d); maximum allowed is %d.\nChange the \"#define MAX_ATOMS\" line in \"constants.h\"\n.", iatom, MAX_ATOMS );
        stop( error_message );
        exit( -1 );
    } else {
        *P_natom = iatom;
        mol.natom = iatom;
    }

    pr( logFile, "Summary of number of atoms with a given atom type:\n");
    pr( logFile, "--------------------------------------------------\n\n");
    for (i=0; i<num_atm_maps; i++) {
        pr( logFile, "Number of atoms with atom type %d = %d\n", i+1, ntype[ i ]);
    }

    pr( logFile, "\nSummary of total charge on ligand, residues and PDBQ:\n");
    pr( logFile, "-----------------------------------------------------\n\n");
    /*
     * Check total charge on ligand
     */
    pr( logFile, "\nTotal charge on ligand =\t\t%+.3f e\n", total_charge_ligand );
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
    pr( logFile, "\nTotal charge on residues =\t\t%+.3f e\n", total_charge_residues );
    iq = (int) ( ( aq = fabs( total_charge_residues ) ) + 0.5 );
    lq = iq - QTOL;
    uq = iq + QTOL;
    if ( ! ((aq >= lq) && (aq <= uq)) ) {
        prStr( message, "\n%s: *** WARNING!  Non-integral total charge (%.3f e) on residues! ***\n\n", programname, total_charge_residues);
        pr_2x( stderr, logFile, message );
    }

    /*
     * Check total charge on all PDBQ atoms
     */
    total_charge = total_charge_ligand + total_charge_residues;
    pr( logFile, "\nTotal charge on all PDBQ atoms (ligand + residues) =\t\t%+.3f e\n", total_charge);
    iq = (int) ( ( aq = fabs( total_charge) ) + 0.5 );
    lq = iq - QTOL;
    uq = iq + QTOL;
    if ( ! ((aq >= lq) && (aq <= uq)) ) {
        prStr( message, "\n%s: *** WARNING!  Non-integral total charge (%.3f e) on all PDBQ atoms! ***\n\n", programname, total_charge);
        pr_2x( stderr, logFile, message );
    }

    /* 
     *  Work out where the torsions are; and what they move...
     *
     *  Also, detect which atoms we should ignore in the intermolecular energy
     *  calculation (ignore_inter[MAX_ATOMS] array)
     */
    mkTorTree( Rec_atomnumber, record, nrecord, 
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
        **  Create list of internal non-bond distances to check...
        */
        nonbonds( crdpdb, nbmatrix_binary, iatom, Rec_atomnumber, nrecord, record, piece, Htype, atmType );
        weedbonds( iatom, pdbaname, piece, ntor, tlist, P_Nnb, Nnbonds, nbmatrix_binary, nonbondlist, outlev, atmType );
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
/* readPDBQLine.cc */

void readPDBQLine( char line[LINE_LEN],
               FloatOrDouble crd[SPACE],
               FloatOrDouble *ptr_q )

/*----------------------------------------------------------------------------*/
{
    char str[4][WORDLEN];

    sscanf( &line[30], "%s %s %s", str[X], str[Y], str[Z] );

    crd[X] = atof( str[X] );
    crd[Y] = atof( str[Y] );
    crd[Z] = atof( str[Z] );

#ifdef DEBUG
    fprintf(stderr, "readPDBQLine: oldpdbq=%d\n%s", oldpdbq, line );
#endif /* DEBUG */

    if (oldpdbq) {
        sscanf( &line[54], "%s", str[3] );
        *ptr_q = atof( str[3] );
    } else {
        sscanf( &line[70], "%s", str[3] );
        *ptr_q = atof( str[3] );
    }

#ifdef DEBUG
    fprintf(stderr, "readPDBQLine:  %.3f, %.3f, %.3f, %.3f\n", crd[X], crd[Y], crd[Z], *ptr_q );
#endif /* DEBUG */
}

/*----------------------------------------------------------------------------*/
/* EOF */
