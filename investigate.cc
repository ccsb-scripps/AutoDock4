/*

 $Id: investigate.cc,v 1.15 2007/04/27 06:01:49 garrett Exp $

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
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <sys/times.h>
#include <sys/param.h>
#include <time.h>
#include "structs.h"
#include "investigate.h"

#define RANDOM_MODE 1
#define CHANGE_MODE 2

extern FILE *logFile;
extern char *programname;


void investigate( int   Nnb,
                    Real charge[MAX_ATOMS],
                    Real abs_charge[MAX_ATOMS],
                    Real qsp_abs_charge[MAX_ATOMS],
                    Boole B_calcIntElec,
                    Real crd[MAX_ATOMS][SPACE],
                    Real crdpdb[MAX_ATOMS][SPACE],

                    EnergyTables *ptr_ad_energy_tables,

                    int   maxTests,
                    Real map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS],
                    int   natom,
                    NonbondParam *nonbondlist,
                    int   ntor,
                    int   outlev,
                    int   tlist[MAX_TORS][MAX_ATOMS],
                    int   type[MAX_ATOMS],
                    Real vt[MAX_TORS][SPACE],
                    Boole B_isGaussTorCon,
                    unsigned short US_torProfile[MAX_TORS][NTORDIVS],
                    Boole B_isTorConstrained[MAX_TORS],
                    Boole B_ShowTorE,
                    unsigned short US_TorE[MAX_TORS],
                    Real F_TorConRange[MAX_TORS][MAX_TOR_CON][2],
                    int   N_con[MAX_TORS],
                    Boole B_symmetry_flag,
                    char  FN_rms_ref_crds[MAX_CHARS],
                    int   OutputEveryNTests,
                    int   NumLocalTests,
                    Real trnStep,
                    Real torStep,
                    
                    int   ignore_inter[MAX_ATOMS],
                    
                    const Boole         B_include_1_4_interactions,
                    const Real scale_1_4,

                    const ParameterEntry parameterArray[MAX_MAPS],

                    const Real unbound_internal_FE,
                    GridMapSetInfo *info,
                    Boole B_use_non_bond_cutoff,
                    Boole B_have_flexible_residues)

{
    Boole B_outside = FALSE;

    int Itor = 0;
    register int Test = -1;
    int indx;
    int ref_natoms = -1;
    register int i = 0;
    //register int XYZ = 0;

    Real e = 0.;
    Real ref_crds[MAX_ATOMS][SPACE];
    Real rms;
    Real MaxRms = 20.0;
    Real RmsBinSize = 0.25;
    Real MinEnergyInRmsBin[NUMRMSBINS];
    int   NumberInRmsBin[NUMRMSBINS];
    int   NumberRandomInRmsBin[NUMRMSBINS];
    int   NumberChangeInRmsBin[NUMRMSBINS];
    int   RmsBinNum = 0;
    //int   NumMaxRms = 0;
    register int NumOutside = 0;
    register int LocalTest = 0;
    int   mode;

    time_t time_seed;

    State sNow; /* qtnNow, torNow */


/*  Initialize
*/
    for (i=0; i<NUMRMSBINS; i++) {
        MinEnergyInRmsBin[i] = BIG;
        NumberInRmsBin[i] = 0;
        NumberRandomInRmsBin[i] = 0;
        NumberChangeInRmsBin[i] = 0;
    }
    sNow.ntor = ntor;

/*  Read in reference coordinates
*/
    if (strncmp(FN_rms_ref_crds,"unspecified filename",20) != 0) {
        if ((ref_natoms = getpdbcrds( FN_rms_ref_crds, ref_crds)) == -1) {
            fprintf( logFile, "%s: Problems while reading \"%s\".\n", programname, FN_rms_ref_crds);
            fprintf( logFile, "Will attempt to use the input PDBQ file coordinates as reference.\n");
        } else if (ref_natoms != natom) {
            pr( logFile, "%s: ERROR!  Wrong number of atoms in reference structure.\n", programname);
            pr( logFile, "Input PDBQ structure has %d atoms, but reference structure has %d atoms.\n\n", natom, ref_natoms);
            ref_natoms = -1;
            exit(-1);
        }
    }

/* Begin investigating the force field,
   by recording the lowest energy for this rmsd
   from the crystal structure
   and from the minimized crystal structure.
*/

    pr( logFile, "\n\n\t\tBEGINNING INVESTIGATION OF FORCE FIELD\n");
    pr( logFile,     "\t\t______________________________________\n\n\n\n" );

/*  Initialize random number generator with a time-dependent seed...  */

    time_seed = time( &time_seed );
    seed_random( time_seed );

    for ( Test = 0;  Test < maxTests;  Test++ ) {

        for (LocalTest = 0; LocalTest < NumLocalTests; LocalTest++, Test++ ) {

            if (LocalTest == 0) {
                mode = RANDOM_MODE;
            } else {
                mode = CHANGE_MODE;
            }

            do { /* while (rms > MaxRms); */
                do { /* while (B_outside); */
                    if (mode == RANDOM_MODE) {
                        sNow = mkRandomState( ntor, F_TorConRange, N_con, info );
                        if (outlev > 2) {
                            fprintf(logFile, "mkRandomState:  ");
                            writeState(logFile, sNow);
                            fflush(logFile);
                        }
                    } else {
                        sNow = changeState( sNow, trnStep, torStep,
                                              ntor, F_TorConRange, N_con);
                        if (outlev > 2) {
                            fprintf(logFile, "changeState:  ");
                            writeState(logFile, sNow);
                            fflush(logFile);
                        }
                    }
                    cnv_state_to_coords( sNow, vt, tlist, ntor, crdpdb, crd, natom );
     
                    /* Check to see if any atom is outside the grid...  */
                    for (i = 0;  i < natom;  i++) {
                        B_outside= is_out_grid_info(crd[i][X], crd[i][Y], crd[i][Z]);
                        if ( B_outside ) {  /* Outside grid! */
                            ++NumOutside;
                            if (mode == CHANGE_MODE) {
                                /* changing pushes out of grid, so switch mode*/
                                mode = RANDOM_MODE;
                            }
                            break;/*...out of i*/
                        }/*outside*/
                    }/*for atoms i*/
                    /* If an atom is outside, do again */
                } while (B_outside);
                /* Now, ligand is inside grid */
                /* Calculate RMSD from reference structure */
                rms = getrms( crd, ref_crds, B_symmetry_flag, natom, type);
            } while (rms > MaxRms);
            /* Calculate Energy of System, */
            e = trilinterp( 0, natom, crd, charge, abs_charge, type, map, info, ALL_ATOMS_INSIDE_GRID, ignore_inter, NULL_ELEC, NULL_EVDW, NULL_ELEC_TOTAL, NULL_EVDW_TOTAL)
                    + eintcal( nonbondlist, ptr_ad_energy_tables, crd, Nnb, B_calcIntElec, B_include_1_4_interactions, scale_1_4, qsp_abs_charge, parameterArray, B_use_non_bond_cutoff, B_have_flexible_residues) - unbound_internal_FE;
            if (B_isGaussTorCon) {
                for (Itor = 0; Itor < ntor; Itor++) {
                    if (B_isTorConstrained[Itor] == 1) {
                        indx = RadiansToDivs( sNow.tor[Itor] );
                        if (B_ShowTorE) {
                            e += (Real)( US_TorE[Itor] 
                                          = US_torProfile[Itor][indx] );
                        } else {
                            e += (Real)US_torProfile[Itor][indx];
                        }
                    }
                }
            }
            /* Update minimum energy for this RMSD bin */
            RmsBinNum = (int)(rms / RmsBinSize);
            ++NumberInRmsBin[RmsBinNum];
            if (mode == RANDOM_MODE) {
                ++NumberRandomInRmsBin[RmsBinNum];
            } else if (mode == CHANGE_MODE) {
                ++NumberChangeInRmsBin[RmsBinNum];
            }
            if (e <= MinEnergyInRmsBin[RmsBinNum]) {
                MinEnergyInRmsBin[RmsBinNum] = e;
            }
            /* Output if it is time, */
            if (outlev > 0) {
                if ((Test+1)%OutputEveryNTests == 0) {
                    fprintf(logFile, "NumberOfTests= %d\n", Test+1);
                    fprintf(logFile, "-------------\n");
                    for (i=0; i<NUMRMSBINS; i++) {
                        fprintf(logFile, "%2d %5.2f-%5.2f:  %9.2f\t%7d\t%7d\t%7d\n", i+1, i*RmsBinSize, (i+1)*RmsBinSize, MinEnergyInRmsBin[i], NumberInRmsBin[i], NumberRandomInRmsBin[i], NumberChangeInRmsBin[i]);
                    }
                    fprintf(logFile, "\n");
                    fprintf(logFile, "NumOutside= %d\n", NumOutside);
                    fprintf(logFile, "\n");
                    fprintf(logFile, "\n");
                    fflush(logFile);
                }
            }

        } /*LocalTest*/
    } /* Loop over Test */
}
/* EOF */
