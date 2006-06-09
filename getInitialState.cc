/*

 $Id: getInitialState.cc,v 1.14 2006/06/09 01:54:16 garrett Exp $

*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/* getInitialState.cc */

#include <math.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/times.h>
#include <sys/param.h>
#include <time.h>
#include "getInitialState.h"
#include "trilinterp.h"
#include "calculateEnergies.h"


extern FILE *logFile;
extern char *programname;


void getInitialState(

            Real *Addr_e0total,
            Real e0max,

            State *sInit, /* was qtn0[QUAT] and tor0[MAX_TORS] */
            State *sMinm, /* was qtnMin[QUAT] and torMin[MAX_TORS] */
            State *sLast, /* was qtnLast[QUAT] and torLast[MAX_TORS] */

            Boole B_RandomTran0,
            Boole B_RandomQuat0,
            Boole B_RandomDihe0,

            Real charge[MAX_ATOMS],
            Real abs_charge[MAX_ATOMS],
            Real qsp_abs_charge[MAX_ATOMS],
            Real q1q2[MAX_NONBONDS],
            Real crd[MAX_ATOMS][SPACE],
            Real crdpdb[MAX_ATOMS][SPACE],
            char  atomstuff[MAX_ATOMS][MAX_CHARS],
            Real elec[MAX_ATOMS],
            Real emap[MAX_ATOMS],

            EnergyTables *ptr_ad_energy_tables,

            Boole B_calcIntElec,
            Real map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS],
            int   natom,
            int   Nnb,
            int   **nonbondlist,
            int   ntor,
            int   tlist[MAX_TORS][MAX_ATOMS],
            int   type[MAX_ATOMS],
            Real vt[MAX_TORS][SPACE],
            int   irun1,
            int   outlev,
            int   MaxRetries,

            Real torsFreeEnergy,

            int   ligand_is_inhibitor,

            int   ignore_inter[MAX_ATOMS],

            const Boole         B_include_1_4_interactions,
            const Real scale_1_4,

            const ParameterEntry parameterArray[MAX_MAPS],

            const Real unbound_internal_FE,

            GridMapSetInfo *info,
            Boole B_use_non_bond_cutoff,
            Boole B_have_flexible_residues
           )

{
    Real e0total = 0.;
    Real e0inter = 0.;
    Real e0intra = 0.;
    Real e0min = BIG;
    int   retries = 0;
    register int i = 0;
    Clock  initStart;
    Clock  initEnd;
    struct tms tms_initStart;
    struct tms tms_initEnd;
    EnergyBreakdown eb;


    /* Store time at start of initialization process... */

    initStart = times( &tms_initStart );

    retries = 0;
    do {
        /*
        ** while e0total, initial energy, is too high,
        */

        /*
        ** Initialize all state variables...
        */
        if (B_RandomTran0) {
            sInit->T.x = random_range( info->lo[X], info->hi[X] );
            sInit->T.y = random_range( info->lo[Y], info->hi[Y] );
            sInit->T.z = random_range( info->lo[Z], info->hi[Z] );
            if (outlev > 1) {
                pr( logFile, "Random initial translation,  tran0 %.3f %.3f %.3f\n", sInit->T.x, sInit->T.y, sInit->T.z);
            }
        }/*if*/
        if (B_RandomQuat0) {
            sInit->Q.nx  = random_range( -1., 1. );
            sInit->Q.ny  = random_range( -1., 1. );
            sInit->Q.nz  = random_range( -1., 1. );
            sInit->Q.ang = Rad( random_range( 0., 360.) );/*convert to radians*/

            mkUnitQuat( &(sInit->Q) );

            if (outlev > 1) {
                pr( logFile, "Random initial quaternion,  quat0 %.3f %.3f %.3f %.1f\n", sInit->Q.nx, sInit->Q.ny, sInit->Q.nz, Deg( sInit->Q.ang ) );
            }
        }/*if*/
        if ( B_RandomDihe0 && (ntor > 0) ) {
            if (outlev > 1) {
                pr( logFile, "Random initial torsions, ndihe = %d\ndihe0 = ", ntor);
            }
            sInit->ntor = ntor;
            for (i=0; i<ntor; i++) {
                sInit->tor[i] = random_range(-180.,180.);
                if (outlev > 1) {
                    pr( logFile, "%7.2f ", sInit->tor[i] ); /*in degrees*/
                }
                sInit->tor[i] = Rad( sInit->tor[i] ); /*now in radians*/
            }
            if (outlev > 1) {
                pr( logFile, "\n");
            }
        }/*if*/

        copyState( sMinm, *sInit );
        copyState( sLast, *sInit );

/* _________________________________________________________________________
**
** Initialize the automated docking simulation,
** _________________________________________________________________________
*/
        initautodock( atomstuff, crd, crdpdb, 
            natom, ntor, sInit, tlist, vt, outlev, info);
        
        e0inter = trilinterp( 0, natom, crd, charge, abs_charge, type, map, 
                    info, ALL_ATOMS_INSIDE_GRID, ignore_inter, elec, emap,
                    NULL_ELEC_TOTAL, NULL_EVDW_TOTAL);
        e0intra = eintcal( nonbondlist, ptr_ad_energy_tables, crd, Nnb, B_calcIntElec, q1q2, B_include_1_4_interactions, scale_1_4, qsp_abs_charge, parameterArray, B_use_non_bond_cutoff, B_have_flexible_residues) - unbound_internal_FE;
        e0total = e0inter + e0intra;

        if (e0total < e0min) {
            /*
            ** This energy is lower so update
            ** the initialization minimum-energy state variables,
            */
            e0min = e0total;
            copyState( sMinm, *sInit );
        }

        /*
        if ( (e0total > e0max) && ((!B_RandomTran0)||(!B_RandomQuat0)) ) {
            B_RandomTran0 = TRUE;
            B_RandomQuat0 = TRUE;
        }
        */

        ++retries;
        if ((retries > 0) && (retries < MaxRetries) &&
            (e0total > e0max) && (outlev > 1)) {

            pr(logFile, "Initial total energy, e0total = %.3f, too high!\n", e0total);
            pr(logFile, "Number of attempts = %d (run %d)\n\n", retries, irun1);
            pr(logFile, "Will try again...\n");

        } else if (retries >= MaxRetries) {

            pr( logFile, "Sorry, too many retries (%d).  Continuing...\n\nWill use the state with the lowest energy found, %.2f\n\n", MaxRetries, e0min);
            e0total = e0min;
            copyState( sInit, *sMinm );

            break;
        }
        fflush( logFile );
    } while ( e0total > e0max );

    cnv_state_to_coords( *sInit, vt, tlist, ntor, crdpdb, crd, natom );

    eb = calculateEnergies( natom, ntor, unbound_internal_FE, torsFreeEnergy, B_have_flexible_residues,
         crd, charge, abs_charge, type, map, info, SOME_ATOMS_OUTSIDE_GRID, 
         ignore_inter, elec, emap, NULL_ELEC_TOTAL, NULL_EVDW_TOTAL,
         nonbondlist, ptr_ad_energy_tables, Nnb, B_calcIntElec, q1q2, 
         B_include_1_4_interactions, scale_1_4, qsp_abs_charge, parameterArray, B_use_non_bond_cutoff );

    copyState( sMinm, *sInit );
    copyState( sLast, *sInit );

    prInitialState( &eb, natom, crd, atomstuff, type, emap, elec, charge, ligand_is_inhibitor, B_have_flexible_residues );

    initEnd = times( &tms_initEnd );
    pr(logFile, "Number of initialization attempts = %d (run %d)\n", retries, irun1);
    pr(logFile, "Time spent initializing: (Real, CPU, System): ");
    timesys( initEnd - initStart, &tms_initStart, &tms_initEnd );

    *Addr_e0total = eb.e_inter + eb.e_intra;
    fflush( logFile );

}
/* EOF */
