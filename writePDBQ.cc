/*

 $Id: writePDBQ.cc,v 1.3 2004/02/12 05:50:49 garrett Exp $

*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/* writePDBQ.cc */

#ifndef WRITEPDBQSTATE
                      /* writePDBQ.cc */
#else /* WRITEPDBQSTATE */
                      /* writeStateOfPDBQ.cc */
#endif /* WRITEPDBQSTATE */

#ifndef WRITEPDBQSTATE
#include <math.h>
#endif /* not WRITEPDBQSTATE */

#include <stdio.h>
#include <string.h>
#include "structs.h"
#include "stateLibrary.h"
#include "writePDBQ.h"

#ifdef WRITEPDBQSTATE
    #include "printEnergies.h"
    #include "assert.h"
#endif /* WRITEPDBQSTATE */



extern int keepresnum;
extern FILE *logFile;

#ifndef WRITEMOLASPDBQFUNC

#ifndef WRITEPDBQSTATE

void writePDBQ(  int irun,

#else /* WRITEPDBQSTATE */

void writeStateOfPDBQ(  int   irun,

#endif /* WRITEPDBQSTATE */
                    char smFileName[MAX_CHARS],
                    char dpfFN[MAX_CHARS],
                    FloatOrDouble sml_center[SPACE],
#ifndef WRITEPDBQSTATE
                    State state,
#else /* WRITEPDBQSTATE */
                    State *Ptr_state,
#endif /* WRITEPDBQSTATE */
                    int ntor,
#ifndef WRITEPDBQSTATE
                    FloatOrDouble eintra,
                    FloatOrDouble einter,
#else /* WRITEPDBQSTATE */
                    FloatOrDouble *Ptr_eintra,
                    FloatOrDouble *Ptr_einter,
#endif /* WRITEPDBQSTATE */
                    int natom,
                    char atomstuff[MAX_ATOMS][MAX_CHARS],
                    FloatOrDouble crd[MAX_ATOMS][SPACE],
                    FloatOrDouble emap[MAX_ATOMS],
                    FloatOrDouble elec[MAX_ATOMS],
                    FloatOrDouble charge[MAX_ATOMS],
                    int ligand_is_inhibitor,
#ifndef WRITEPDBQSTATE
                    FloatOrDouble torsFreeEnergy,
                    int   outlev)
#else /* WRITEPDBQSTATE */
                    FloatOrDouble torsFreeEnergy,
                    FloatOrDouble vt[MAX_TORS][SPACE],
                    int   tlist[MAX_TORS][MAX_ATOMS],
                    FloatOrDouble crdpdb[MAX_ATOMS][SPACE],
                    int   nonbondlist[MAX_NONBONDS][4],
                    FloatOrDouble e_internal[NEINT][ATOM_MAPS][ATOM_MAPS],
                    int   type[MAX_ATOMS],
                    int   Nnb,
                    Boole B_calcIntElec,
                    FloatOrDouble q1q2[MAX_NONBONDS],
                    FloatOrDouble map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS],
                    FloatOrDouble inv_spacing,
                    FloatOrDouble xlo,
                    FloatOrDouble ylo,
                    FloatOrDouble zlo,
                    FloatOrDouble xhi,
                    FloatOrDouble yhi,
                    FloatOrDouble zhi,
                    Boole B_template,
                    FloatOrDouble template_energy[MAX_ATOMS],
                    FloatOrDouble template_stddev[MAX_ATOMS],
                    int   outlev)
#endif /* WRITEPDBQSTATE */

{
    int   i=0;
    FloatOrDouble MaxValue = 99.99;
#ifndef WRITEPDBQSTATE
    char rec14[14],rec9[9];
#else /* WRITEPDBQSTATE */
    char  AtmNamResNamNum[14],AtmNamResNam[9];
    Boole B_outside = FALSE;

    for (i=0; i<14; i++) {
        AtmNamResNamNum[i] = '\0';
    }

    for (i=0; i<9; i++) {
        AtmNamResNam[i] = '\0';
    }

// new stuff begins

    if ((outlev > -1) && (outlev < 3)) {                        // gmm 2001-11-01
        // in:
        // printState( logFile, (*Ptr_state), 2 );           // before 2001-11-01
        // the 2 is the level of detail: 2 is high, 0 is low
        printState( logFile, (*Ptr_state), outlev );            // gmm 2001-11-01
    } else if (outlev == -1) {                                  // gmm 2001-11-01
        printState( logFile, (*Ptr_state), 0 );                 // gmm 2001-11-01
    }                                                           // gmm 2001-11-01

    cnv_state_to_coords( (*Ptr_state),vt,tlist,ntor,crdpdb,crd,natom );
    B_outside = FALSE;
    for (i=0; (i<natom)&&(!B_outside); i++) {
        B_outside = is_out_grid( crd[i][0], crd[i][1], crd[i][2]);
    }
    if (!B_outside) {
        if (ntor > 0) {
            *Ptr_eintra = eintcal( nonbondlist, e_internal, crd, Nnb, B_calcIntElec, q1q2);
        } else {
            *Ptr_eintra = 0.0;
        }
        if (B_template) {
            *Ptr_einter = byatom_template_trilinterp( crd, charge, type, natom, map, inv_spacing, elec, emap, xlo, ylo, zlo,
                                               template_energy, template_stddev);
        } else {
            *Ptr_einter = trilinterp( crd, charge, type, natom, map, inv_spacing, elec, emap, xlo, ylo, zlo);
        }
    } else {
        *Ptr_eintra = *Ptr_einter = BIG;  // defined in constants.h
    }
    // new stuff ends

#endif /* WRITEPDBQSTATE */

    if (outlev > -1) {                                          // gmm 2001-11-01
    // output of coordinates                                    // gmm 2001-11-01

#ifndef WRITEPDBQSTATE
        pr( logFile, "\n" );
#endif /* not WRITEPDBQSTATE */
        pr( logFile, "DOCKED: MODEL     %4d\n", irun+1 );
        pr( logFile, "DOCKED: USER    Run = %d\n", irun+1 );
        pr( logFile, "DOCKED: USER    DPF = %s\n", dpfFN );

#ifndef WRITEPDBQSTATE
        printEnergies( einter, eintra, torsFreeEnergy, "DOCKED: USER    ", ligand_is_inhibitor);
#else /* WRITEPDBQSTATE */
        printEnergies( *Ptr_einter, *Ptr_eintra, torsFreeEnergy, "DOCKED: USER    ", ligand_is_inhibitor);
#endif /* WRITEPDBQSTATE */

        (void)fprintf( logFile, "DOCKED: USER    NEWDPF move %s\n", smFileName );
        (void)fprintf( logFile, "DOCKED: USER    NEWDPF about %f %f %f\n", sml_center[X],sml_center[Y],sml_center[Z]);
#ifndef WRITEPDBQSTATE
        (void)fprintf( logFile, "DOCKED: USER    NEWDPF tran0 %f %f %f\n", state.T.x, state.T.y, state.T.z );
        (void)fprintf( logFile, "DOCKED: USER    NEWDPF quat0 %f %f %f %f\n", state.Q.nx, state.Q.ny, state.Q.nz, Deg( state.Q.ang ) );
#else /* WRITEPDBQSTATE */
        (void)fprintf( logFile, "DOCKED: USER    NEWDPF tran0 %f %f %f\n", (*Ptr_state).T.x, (*Ptr_state).T.y, (*Ptr_state).T.z );
        (void)fprintf( logFile, "DOCKED: USER    NEWDPF quat0 %f %f %f %f\n", (*Ptr_state).Q.nx, (*Ptr_state).Q.ny, (*Ptr_state).Q.nz, Deg( (*Ptr_state).Q.ang ) );
#endif /* WRITEPDBQSTATE */
        if (ntor > 0) {
            (void)fprintf( logFile, "DOCKED: USER    NEWDPF ndihe %d\n", ntor );
            (void)fprintf( logFile, "DOCKED: USER    NEWDPF dihe0 " );
            for (i=0; i<ntor; i++) {
#ifndef WRITEPDBQSTATE
                (void)fprintf( logFile, "%.2f ", Deg(state.tor[i]) );
#else /* WRITEPDBQSTATE */
                (void)fprintf( logFile, "%.2f ", Deg((*Ptr_state).tor[i]) );
#endif /* WRITEPDBQSTATE */
            }/*i*/
            (void)fprintf( logFile, "\n" );
        }/*endif*/

        (void)fprintf( logFile, "DOCKED: USER                              x       y       z     vdW  Elec       q\n" );
        (void)fprintf( logFile, "DOCKED: USER                           _______ _______ _______ _____ _____    ______\n" );
        if (keepresnum > 0) {
            for (i = 0;  i < natom;  i++) {
#ifndef WRITEPDBQSTATE
                strncpy( rec14, &atomstuff[i][13], (size_t)13);
                (void)fprintf(logFile, FORMAT_PDBQ_ATOM_RESSTR, "DOCKED: ", i+1, rec14, crd[i][X], crd[i][Y], crd[i][Z], min(emap[i], MaxValue), min(elec[i], MaxValue), charge[i]);
                (void)fprintf(logFile, "\n");
#else /* WRITEPDBQSTATE */
                assert( i>=0 && i < natom);
                strncpy( AtmNamResNamNum, &atomstuff[i][13], (size_t)13);
                AtmNamResNamNum[13] = '\0'; // ensure string is NULL terminated!
                (void)fprintf(logFile, FORMAT_PDBQ_ATOM_RESSTR, "DOCKED: ", i+1, AtmNamResNamNum,      crd[i][X], crd[i][Y], crd[i][Z], min(emap[i], MaxValue), min(elec[i], MaxValue), charge[i]);
                (void)fprintf(logFile, "\n");
#endif /* WRITEPDBQSTATE */
            }/*i*/
        } else {
            for (i = 0;  i < natom;  i++) {
#ifndef WRITEPDBQSTATE
                strncpy( rec9, &atomstuff[i][13], (size_t)8);
                (void)fprintf(logFile, FORMAT_PDBQ_ATOM_RESNUM, "DOCKED: ", i+1, rec9,irun+1,crd[i][X], crd[i][Y], crd[i][Z], min(emap[i], MaxValue), min(elec[i], MaxValue), charge[i]);
                (void)fprintf(logFile, "\n");
#else /* WRITEPDBQSTATE */
                assert( i>=0 && i < natom);
                strncpy( AtmNamResNam, &atomstuff[i][13], (size_t)8);
                AtmNamResNam[8] = '\0'; // ensure string is NULL terminated!
                (void)fprintf(logFile, FORMAT_PDBQ_ATOM_RESNUM, "DOCKED: ", i+1, AtmNamResNam,irun+1,crd[i][X], crd[i][Y], crd[i][Z], min(emap[i], MaxValue), min(elec[i], MaxValue), charge[i]);
                (void)fprintf(logFile, "\n");
#endif /* WRITEPDBQSTATE */
            }/*i*/
        }/*endif*/
        (void)fprintf( logFile, "DOCKED: TER\n" );
        (void)fprintf( logFile, "DOCKED: ENDMDL\n" );
        // (void)fprintf( logFile, UnderLine ); // gmm 2001-11-01
        (void)fflush(  logFile );
        // output of coordinates // gmm 2001-11-01
    }
}

#else /* WRITEMOLASPDBQFUNC */

/***************************************************************/

void writeMolAsPDBQ(Molecule *mol, FILE *output)

/* write out the Molecule "mol" as a PDBQ formatted file to 
   the file pointed to by "output". */
{
    int i;
    char rec14[14],rec9[9];

    (void)fprintf( output, "\n" );
    (void)fprintf( output, "DIAGNOSTIC: MODEL     %4d\n", 1 );

    (void)fprintf( output, "DIAGNOSTIC: USER    NEWDPF tran0 %f %f %f\n", mol->S.T.x, mol->S.T.y, mol->S.T.z );
    (void)fprintf( output, "DIAGNOSTIC: USER    NEWDPF quat0 %f %f %f %f\n", mol->S.Q.nx, mol->S.Q.ny, mol->S.Q.nz, ( (mol->S.Q.ang) * 57.29577951 ) );
    if (mol->S.ntor > 0) {
        (void)fprintf( output, "DIAGNOSTIC: USER    NEWDPF ndihe %d\n", mol->S.ntor );
        (void)fprintf( output, "DIAGNOSTIC: USER    NEWDPF dihe0 " );
        for (i=0; i<mol->S.ntor; i++) {
            (void)fprintf( output, "%.2f ", ( (mol->S.tor[i]) * 57.29577951 ) );
        } 
        (void)fprintf( output, "\n" );
    } 

    (void)fprintf( output, "DIAGNOSTIC: USER                              x       y       z     vdW  Elec       q\n" );
    (void)fprintf( output, "DIAGNOSTIC: USER                           _______ _______ _______ _____ _____    ______\n" );
    if (keepresnum > 0) {
        for (i = 0;  i < mol->natom;  i++) {
            strncpy( rec14, &(mol->atomstr[i][13]), (size_t)13);
            (void)fprintf(output, FORMAT_PDBQ_ATOM_RESSTR, "DIAGNOSTIC: ", i+1, rec14,      mol->crd[i][0], mol->crd[i][1], mol->crd[i][2], 0.0, 0.0, 0.0);
            (void)fprintf(output, "\n");
        } 
    } else {
        for (i = 0;  i < mol->natom;  i++) {
            strncpy( rec9, &(mol->atomstr[i][13]), (size_t)8);
            (void)fprintf(output, FORMAT_PDBQ_ATOM_RESNUM, "DIAGNOSTIC: ", i+1, rec9,1, mol->crd[i][0],  mol->crd[i][1],  mol->crd[i][2], 0.0, 0.0, 0.0);
            (void)fprintf(output, "\n");
        } 
    } 
    (void)fprintf( output, "DIAGNOSTIC: TER\n" );
    (void)fprintf( output, "DIAGNOSTIC: ENDMDL\n" );
    (void)fprintf( output, "_______________________________________________________________________________\n\n" );
    (void)fflush(  output );
}

#endif /* WRITEMOLASPDBQFUNC */
/* EOF */
