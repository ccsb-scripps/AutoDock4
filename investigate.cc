/*

 $Id: investigate.cc,v 1.7 2005/09/28 22:54:20 garrett Exp $

*/

/* investigate.cc */

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
                    FloatOrDouble charge[MAX_ATOMS],
                    FloatOrDouble abs_charge[MAX_ATOMS],
                    FloatOrDouble qsp_abs_charge[MAX_ATOMS],
                    Boole B_calcIntElec,
                    FloatOrDouble q1q2[MAX_NONBONDS],
                    FloatOrDouble crd[MAX_ATOMS][SPACE],
                    FloatOrDouble crdpdb[MAX_ATOMS][SPACE],

                    EnergyTables *ptr_ad_energy_tables,

                    int   maxTests,
                    FloatOrDouble map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS],
                    int   natom,
                    int   nonbondlist[MAX_NONBONDS][MAX_NBDATA],
                    int   ntor,
                    int   outlev,
                    int   tlist[MAX_TORS][MAX_ATOMS],
                    int   type[MAX_ATOMS],
                    FloatOrDouble vt[MAX_TORS][SPACE],
                    Boole B_isGaussTorCon,
                    unsigned short US_torProfile[MAX_TORS][NTORDIVS],
                    Boole B_isTorConstrained[MAX_TORS],
                    Boole B_ShowTorE,
                    unsigned short US_TorE[MAX_TORS],
                    FloatOrDouble F_TorConRange[MAX_TORS][MAX_TOR_CON][2],
                    int   N_con[MAX_TORS],
                    Boole B_symmetry_flag,
                    char  FN_rms_ref_crds[MAX_CHARS],
                    int   OutputEveryNTests,
                    int   NumLocalTests,
                    FloatOrDouble trnStep,
                    FloatOrDouble torStep,
                    
                    int   ignore_inter[MAX_ATOMS],
                    
                    const Boole         B_include_1_4_interactions,
                    const FloatOrDouble scale_1_4,

                    const ParameterEntry parameterArray[MAX_MAPS],

                    const FloatOrDouble unbound_internal_FE,
                    GridMapSetInfo *info )

{
    Boole B_outside = FALSE;

    int Itor = 0;
    register int Test = -1;
    int indx;
    int ref_natoms = -1;
    register int i = 0;
    //register int XYZ = 0;

    FloatOrDouble e = 0.;
    FloatOrDouble ref_crds[MAX_ATOMS][SPACE];
    FloatOrDouble rms;
    FloatOrDouble MaxRms = 20.0;
    FloatOrDouble RmsBinSize = 0.25;
    FloatOrDouble MinEnergyInRmsBin[NUMRMSBINS];
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
    pr( logFile,     "\t\t_____________________________________\n\n\n\n" );

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
            e = quicktrilinterp4( crd, charge, abs_charge, type, natom, map, ignore_inter, info) 
                    + eintcal( nonbondlist, ptr_ad_energy_tables, crd, Nnb, B_calcIntElec, q1q2, B_include_1_4_interactions, scale_1_4, qsp_abs_charge, parameterArray, unbound_internal_FE);
            if (B_isGaussTorCon) {
                for (Itor = 0; Itor < ntor; Itor++) {
                    if (B_isTorConstrained[Itor] == 1) {
                        indx = Rad2Div( sNow.tor[Itor] );
                        if (B_ShowTorE) {
                            e += (FloatOrDouble)( US_TorE[Itor] 
                                          = US_torProfile[Itor][indx] );
                        } else {
                            e += (FloatOrDouble)US_torProfile[Itor][indx];
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
