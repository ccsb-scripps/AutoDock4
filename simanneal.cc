/*

 $Id: simanneal.cc,v 1.2 2003/02/26 01:40:37 garrett Exp $

*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/* simanneal.cc */

#include <math.h>

    #include <stdlib.h>
    #include <stdio.h>
    #include <string.h>
    #include <sys/types.h>
    #include <sys/times.h>
    #include <sys/param.h>
    #include <time.h>
    #include "simanneal.h"
    #include "energy.h"



extern FILE *logFile;
extern char *programname;


void simanneal( int   *Addr_nconf,
		int   Nnb,
		FloatOrDouble WallEnergy,
		char  atomstuff[MAX_ATOMS][MAX_CHARS],
		FloatOrDouble charge[MAX_ATOMS],
		Boole B_calcIntElec,
		FloatOrDouble q1q2[MAX_NONBONDS],
		FloatOrDouble crd[MAX_ATOMS][SPACE],
		FloatOrDouble crdpdb[MAX_ATOMS][SPACE],
		char  FN_dpf[MAX_CHARS],
		FloatOrDouble e_internal[NEINT][ATOM_MAPS][ATOM_MAPS],
		FloatOrDouble econf[MAX_RUNS],
		Boole B_either,
		FloatOrDouble elec[MAX_ATOMS],
		FloatOrDouble emap[MAX_ATOMS],
		FloatOrDouble xhi,
		FloatOrDouble yhi,
		FloatOrDouble zhi,
		int   NcycMax,
		FloatOrDouble inv_spacing,
		int   irunmax,
		Clock jobStart,
		FloatOrDouble xlo,
		FloatOrDouble ylo,
		FloatOrDouble zlo,
		FloatOrDouble map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS],
		int   naccmax,
		int   natom,
		int   nonbondlist[MAX_NONBONDS][2],
		int   nrejmax,
		int   ntor1,
		int   ntor,
		int   outlev,

		State sInit, /* tor0, qtn0 */
		State sHist[MAX_RUNS], /* was qtnHist, torHist */

		FloatOrDouble qtwFac,
		Boole B_qtwReduc,
		FloatOrDouble qtwStep0,
		Boole B_selectmin,
		char  FN_ligand[MAX_CHARS],
		FloatOrDouble sml_center[SPACE],
		FloatOrDouble RT0,
		Boole B_RTChange,
		FloatOrDouble RTFac,
		struct tms tms_jobStart,
		int   tlist[MAX_TORS][MAX_ATOMS],
		FloatOrDouble torFac,
		Boole B_torReduc,
		FloatOrDouble torStep0,
		char  FN_trj[MAX_CHARS],
		int   trj_cyc_max,
		int   trj_cyc_min,
		int   trj_freq,
		FloatOrDouble trnFac,
		Boole B_trnReduc,
		FloatOrDouble trnStep0,
		int   type[MAX_ATOMS],
		FloatOrDouble vt[MAX_TORS][SPACE],
		Boole B_writeTrj,
		Boole B_constrain,
		int   atomC1,
		int   atomC2,
		FloatOrDouble sqlower,
		FloatOrDouble squpper,
		Boole B_linear_schedule,
		FloatOrDouble RTreduc,
		/*FloatOrDouble maxrad,*/
		Boole B_watch,
		char  FN_watch[MAX_CHARS],
		Boole B_isGaussTorCon,
		unsigned short US_torProfile[MAX_TORS][NTORDIVS],
		Boole B_isTorConstrained[MAX_TORS],
		Boole B_ShowTorE,
		unsigned short US_TorE[MAX_TORS],
		FloatOrDouble F_TorConRange[MAX_TORS][MAX_TOR_CON][2],
		int N_con[MAX_TORS],
		Boole B_RandomTran0,
		Boole B_RandomQuat0,
		Boole B_RandomDihe0,
		FloatOrDouble e0max,
		FloatOrDouble torsFreeEnergy,
		int   MaxRetries,
    int   ligand_is_inhibitor)

{
    char message[LINE_LEN];


    FILE *FP_trj;

    State sNow; /* qtnNow, torNow */
    State sChange; /* qtnChange, torChange */
    State sLast; /* qtnLast, torLast */
    State sMin; /* qtnMin, torMin */
    State sSave; /* qtnSave, torSave */

    FloatOrDouble d[SPACE];
    FloatOrDouble e = 0.;
    FloatOrDouble e0 = 0.;
    FloatOrDouble einter = 0.;
    FloatOrDouble eintra = 0.;
    FloatOrDouble eLast = 0.;
    FloatOrDouble eMin = BIG_ENERGY;
    FloatOrDouble etot = 0.0;
    FloatOrDouble inv_RT = 0.;
    FloatOrDouble qtwStep;
    FloatOrDouble rsqC1C2;
    FloatOrDouble RT = 616.;
    FloatOrDouble torTmp;
    FloatOrDouble torStep;
    FloatOrDouble trnStep;
    /* ** FloatOrDouble xloTrn; ** FloatOrDouble xhiTrn; ** FloatOrDouble yloTrn; ** FloatOrDouble yhiTrn; ** FloatOrDouble zloTrn; ** FloatOrDouble zhiTrn; ** FloatOrDouble lo[SPACE]; ** FloatOrDouble trnStepHi; ** FloatOrDouble qtwStepHi; ** FloatOrDouble torStepHi; */

    Boole B_inRange = FALSE;
    Boole B_outside = FALSE;
    Boole B_within_constraint = TRUE;

    int count = 0;
    int Itor = 0;
    int icycle = -1;
    int icycle1 = -1;
    int irun = -1;
    int irun1 = -1;
    int indx;

    register int i = 0;
    register int nacc = 0;
    register int nAcc = 0;
    register int nAccProb = 0;
    register int nedge = 0;
    register int nrej = 0;
    register int ntot = 1;
    register int xyz = 0;

    struct tms tms_cycStart;
    struct tms tms_cycEnd;
    struct tms tms_jobEnd;

    Clock cycStart;
    Clock cycEnd;
    Clock jobEnd;

    time_t time_seed;

    /* trnStepHi = HI_NRG_JUMP_FACTOR * trnStep; ** qtwStepHi = HI_NRG_JUMP_FACTOR * qtwStep; ** torStepHi = HI_NRG_JUMP_FACTOR * torStep; */
    /* lo[X] = xlo;  lo[Y] = ylo;  lo[Z] = zlo;*/
    /* xloTrn = xlo + maxrad; ** xhiTrn = xhi - maxrad; ** yloTrn = ylo + maxrad; ** yhiTrn = yhi - maxrad; ** zloTrn = zlo + maxrad; ** zhiTrn = zhi - maxrad; */

/* Open the trajectory file for writing, =====================================*/

    if ( B_writeTrj ) { 
        if ( (FP_trj = fopen(FN_trj, "w")) == NULL ) {
            prStr( message, "\n%s: can't create trajectory file %s\n", programname, FN_trj);
            pr_2x( stderr, logFile, message );
            prStr( message, "\n%s: Unsuccessful Completion.\n\n", programname);
            pr_2x( stderr, logFile, message );
            jobEnd = times( &tms_jobEnd );
            timesys( jobEnd - jobStart, &tms_jobStart, &tms_jobEnd );
            pr_2x( logFile, stderr, UnderLine );
            exit( -1 );
        }/*END PROGRAM*/
    }/*endif*/

/* Begin the automated docking simulation, ===================================*/

    pr( logFile, "\n\n\t\tBEGINNING MONTE CARLO SIMULATED ANNEALING\n");
    pr( logFile, "     \t\t_________________________________________\n\n\n\n" );


    for ( irun = 0;  irun < irunmax;  irun++ ) { /*===========================*/

        irun1 = 1 + irun;

	/*
	** Initialize random number generator with a time-dependent seed...
	*/
	time_seed = time( &time_seed );
	seed_random( time_seed );

	if (outlev > 0) {
	    pr(logFile, "\n\tINITIALIZING AUTOMATED DOCKING SIMULATION\n" );
	    pr(logFile, "\t_________________________________________\n\n" );
	    pr( logFile, "RUN %d...\n\n", irun1);
	    pr( logFile, "Time-dependent Seed:  %ld\n\n", time_seed);
	}

        if ( B_writeTrj ) {
            pr( FP_trj, "ntorsions %d\nrun %d\n", ntor, irun1 );
            fflush( FP_trj );
        }

	getInitialState( &e0, e0max,
			 &sInit, &sMin, &sLast, 
			 B_RandomTran0, B_RandomQuat0, B_RandomDihe0, 
			 charge, q1q2, crd, crdpdb, atomstuff,
			 elec, emap, e_internal, B_calcIntElec,
			 xhi, yhi, zhi, xlo, ylo, zlo,
			 inv_spacing, map, natom, Nnb, nonbondlist,
			 ntor, tlist, type, vt, irun1, outlev, MaxRetries,
			 torsFreeEnergy, ligand_is_inhibitor);

        RT = RT0;		/* Initialize the "annealing" temperature */
	if (RT <= APPROX_ZERO) { RT = 616.; }
	inv_RT = 1. / RT;

        eMin    = min(BIG_ENERGY, e0);
        eLast   = e0;
        trnStep = trnStep0;	/* translation*/
        qtwStep = qtwStep0;	/* quaternion angle, w*/
        torStep = torStep0;	/* torsion angles*/

        if (outlev > 0) {
	    pr( logFile, "\n\n\t\tBEGINNING SIMULATED ANNEALING");
	    pr( logFile, "\n\t\t_____________________________\n\n");
	    pr( logFile, "\n      \t      \tMinimum     Average     | Acc/    Accepted:    Rejected:     |          |  xyz-Translation  |        Time:        \n");
	      pr( logFile, "Run:  \tCycle:\tEnergy:     Energy:     |   /Rej: Down:  Up:   Total: Edge:  |   RT:    |   of Min.Energy   |  Real, CPU, System  \n" );
	      pr( logFile, "______\t______\t___________ ___________ | ______ ______ ______ ______ ______ | ________ | _________________ | ____________________\n" );
/*	                                  12345678901 12345678901   123456 123456 123456 123456 123456   12345678   12345 12345 12345
**	                 "%D /%D\T%D /%D\T%+11.2F     %+11.2F       %6.2F  %6D    %6D    %6D    %6D      %8.1F     %5.2F %5.2F %5.2F   ",
**	                 IRUN1, IRUNMAX, ICYCLE1, ICYCLEMAX, EmIN, ETOT/NTOT, QTNmIN[x], QTNmIN[y], QTNmIN[z], (NREJ!=0) ? (FLOAT)NACC/NREJ : -999., NACC, NREJ, NEDGE, rt
*/
        }
        fflush(logFile);
/*____________________________________________________________________________*/

        for ( icycle = 0;  icycle < NcycMax;  icycle++ ) {

            cycStart = times( &tms_cycStart );

            icycle1 = icycle + (ntot = 1);
            B_inRange = (icycle >= trj_cyc_min) && (icycle <= trj_cyc_max);
            if ( B_writeTrj && B_inRange ) {
		pr( FP_trj, "cycle %d\ntemp %f\n", icycle1, RT );
		fflush(  FP_trj  );
            }
            nAcc = nAccProb = nacc = nrej = nedge = 0;
	    etot = eLast;

/*____________________________________________________________________________*/

            do {
		/*
		** Do one Monte Carlo step,
                ** while the number of accepted steps < naccmax 
		** and number of rejected steps < nrejmax...
		*/

		mkNewState( &sNow, &sLast, &sChange,
                    vt, tlist, ntor, crd, crdpdb, natom,
		    trnStep,
		    /*qtwStep,*/
		    torStep,
		    F_TorConRange, N_con);

		if (B_constrain) {
		    for (xyz = 0;  xyz < SPACE;  xyz++) {
			d[xyz] = crd[atomC1][xyz] - crd[atomC2][xyz];
		    }
		    rsqC1C2 = sqhypotenuse(d[X],  d[Y], d[Z]);
		    if (! (B_within_constraint = (rsqC1C2 > sqlower) && 
						 (rsqC1C2 < squpper))) {
			copyState( &sNow, sLast );
		    }
		}/*if B_constrain*/

		if (B_within_constraint) {
		    /* 
		    ** Normally, this is true.
		    ** If the distance-constraint was set,
		    ** and was just violated, this is false.
		    */

		    for (i = 0;  i < natom;  i++) {
			B_outside= is_out_grid(crd[i][X], crd[i][Y], crd[i][Z]);
			if ( B_outside ) {
			    /*
			    ** Outside grid!
			    */
			    ++nedge;
			    ++nrej;
			    /*pr(logFile, "e"); / *###*/
			    /* etot += WallEnergy; ++ntot;*/
			    /*
			    ** Undo this move,
			    */
			    copyState( &sNow, sLast );
			    /*
			    ** Subtract just the translation;
			    ** hence "doubling back" from the edge.
			    */
			    sNow.T.x -= sChange.T.x;
			    sNow.T.y -= sChange.T.y;
			    sNow.T.z -= sChange.T.z;

			    if ( B_writeTrj && B_inRange && B_either && (++count == trj_freq)) {
				count = 0;
				e = eintra = WallEnergy;
				output_state( FP_trj, sNow, ntor, nacc+nrej, e,
				    eintra, (char)'e', B_watch, FN_watch, 
				    atomstuff, natom, crd);
			    }/*writeTrj*/

			    break;/*...out of i*/

			}/*outside*/
		    }/*for atoms i*/

		    if ( !B_outside ){ /*inside grid maps*/

 			/* Calculate Energy of System, =======================*/

			/*
			** MORE ACCURATE METHOD, (SLOWER):
			*/
			e = quicktrilinterp( crd, charge, type, natom, map, inv_spacing, xlo, ylo, zlo) + (eintra = eintcal( nonbondlist, e_internal, crd, type, Nnb, B_calcIntElec, q1q2));

			/*
			** LESS ACCURATE  METHOD (FASTER):
			** (alternative to above line).
			** 
			** Only calculate internal energy if 
			** trilinterp energy is low enough.
			** 
			** if ( (e = quicktrilinterp( crd, charge, type, natom, 
			** map, inv_spacing, xlo, ylo, zlo)) < 
			** ENERGY_CUTOFF) { 
			**   e += (eintra = eintcal( nonbondlist, e_internal,
			**   crd, type, Nnb, B_calcIntElec, q1q2 ));
			** }
			*/

			if (B_isGaussTorCon) {
			    /*** This looks wrong... for (Itor = 0; Itor <= ntor; Itor++) { ***/
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

			/* Apply the Metropolis energy test... ===============*/

			if (e <= eLast) {
			    /*
			    **  Accept this move immediately.
			    */
			    ++nacc;
			    ++nAcc;
			    etot += (eLast = e);
			    ++ntot;
			    copyState( &sLast, sNow );

			    /* pr(logFile, "A"); / *###*/
			    if (e < eMin) {
				/*
				** Update minimum-energy state variables,
				*/
				eMin = e;
				copyState( &sMin,  sNow );
			    }
			    if ( B_writeTrj && B_inRange && (++count == trj_freq)){
				count = 0;
				output_state( FP_trj, sNow, ntor, nacc+nrej, e,
				    eintra, (char)'A', B_watch, FN_watch, 
				    atomstuff, natom, crd);
			    }/*write trajectory*/

			} else {

                            /* Probabilistic move. ===========================*/

			    if (exp((double)((eLast-e)*inv_RT))<local_random()){
				/*
				** Failed the probability test. 
				** Reject this move.
				*/
				++nrej;
				/* pr(logFile, "R"); / *###*/
				copyState( &sNow, sLast );

				if ( B_writeTrj && B_inRange && B_either && 
				    (++count == trj_freq)) {
				    count = 0;
				    output_state( FP_trj, sNow, ntor, nacc+nrej,
				        e, eintra, (char)'R', B_watch, FN_watch,
					atomstuff, natom, crd);
				}/*write trajectory*/

			    } else {
				/*
				** Passed the probability test.  
				** Accept this move.
				** A chance to escape a local minimum...
				*/
				++nacc;
				++nAccProb;
				etot += e;
				++ntot;
				/*pr(logFile, "a"); / *###*/

				eLast = e;
				copyState( &sLast, sNow );
				if ( B_writeTrj && B_inRange && (++count == trj_freq))  {
				    count = 0;
				    output_state(FP_trj, sNow, ntor, nacc+nrej,
					e, eintra, (char)'a', B_watch, FN_watch,
					atomstuff, natom, crd);
				}/*write trajectory*/
			    }/*passed Monte Carlo probablility test*/
			}/*e > eLast, Probabilistic move...*/
		    }/*inside grid maps*/
		}/*within_constraint*/

            } while ( nacc < naccmax  &&  nrej < nrejmax );
/*____________________________________________________________________________*/

	    if ((nacc == 0) && (nedge == nrejmax)) {
		/*
		**  Clear indication that ligand got stuck on an edge...
		*/
		pr(logFile, "\n\n>>> Ligand appears to be stuck on an edge: forced to re-initialize. <<<\n\n");
		--icycle;
		getInitialState( &e0, e0max,
			 &sInit, &sMin, &sLast, 
			 TRUE, TRUE, TRUE, 
			 charge, q1q2, crd, crdpdb, atomstuff,
			 elec, emap, e_internal, B_calcIntElec,
			 xhi, yhi, zhi, xlo, ylo, zlo,
			 inv_spacing, map, natom, Nnb, nonbondlist,
			 ntor, tlist, type, vt, irun1, outlev, MaxRetries,
			 torsFreeEnergy, ligand_is_inhibitor);

	    } else {

		if ( B_trnReduc )  trnStep *= trnFac;       
		if ( B_qtwReduc )  qtwStep *= qtwFac;
		if ( B_torReduc )  torStep *= torFac;
		/*
		**  Output-level dependent diagnostics...
		*/
		if (outlev > 0) {
		    /*pr(logFile, "\n"); / *###*/
		    pr( logFile, "%d /%d\t%d /%d\t%+11.2f %+11.2f   %6.2f %6d %6d %6d %6d   %8.1f   %5.2f %5.2f %5.2f   ", irun1, irunmax, icycle1, NcycMax, eMin, etot/ntot, (nrej!=0) ? (FloatOrDouble)nacc/nrej : 999.99, nAcc, nAccProb, nrej, nedge, RT, sMin.T.x, sMin.T.y, sMin.T.z );
		    cycEnd = times( &tms_cycEnd );
		    timesys( cycEnd - cycStart, &tms_cycStart, &tms_cycEnd );
		    if (outlev > 1) {
			pr( logFile, "\tEnergy:   \tState:\n\t__________\t____________________________________________________________\nMinimum\t%+6.2f\t(%+.2f,%+.2f,%+.2f), q = [w,(x,y,z)] = [%5.1f deg, (%+.2f,%+.2f,%+.2f)],\n", eMin, sMin.T.x, sMin.T.y, sMin.T.z, Deg(sMin.Q.ang) , sMin.Q.nx, sMin.Q.ny, sMin.Q.nz );
			pr( logFile, "\nLast\t%+6.2f\t(%+.2f,%+.2f,%+.2f), q = [w,(x,y,z)] = [%5.1f deg, (%+.2f,%+.2f,%+.2f)],\n", eLast, sLast.T.x, sLast.T.y, sLast.T.z, Deg(sLast.Q.ang) , sLast.Q.nx, sLast.Q.ny, sLast.Q.nz );
			if (ntor > 0) {
			    pr( logFile, "Minimum:\t(" );
			    for (i=0; i<ntor; i++) {
				pr( logFile, "%.1f%s ", Deg(sMin.tor[i]), (i < ntor1)?",":" deg)" );
			    }
			    pr( logFile, "\nLast:\t(" );
			    for (i=0; i<ntor; i++) {
				pr( logFile, "%.1f%s ", Deg(sLast.tor[i]), (i < ntor1)?",":" deg)" );
			    }
			    pr( logFile, "\n" );
			}
			if ( B_trnReduc )
			    pr( logFile, "\nTranslation step size reduced; now =\t\t +/- %.2f A\n", trnFac);
			if ( B_qtwReduc )
			    pr( logFile, "\nQuaternion Rotation step size reduced; now =\t +/- %.2f deg\n", qtwFac);
			if ( B_torReduc )
			    pr( logFile, "\nTorsion step size reduced; now =\t\t +/- %.2f deg\n", torFac);
		    }/*outlev > 1*/
		    flushLog;
		}/*outlev > 0*/
		/*
		** Reduce temperature,
		*/
		if ( B_linear_schedule ) {
		    RT -= RTreduc;

		    if (RT <= APPROX_ZERO) {
			inv_RT = 1. / APPROX_ZERO;
		    } else {
			inv_RT = 1. / RT;
		    }
		} else if ( B_RTChange ) {
		    inv_RT = 1./( RT *= RTFac );
		}
		/*
		** Start next cycle at minimum state?
		*/
		if ( B_selectmin ) {
		    eLast = eMin;
		    copyState( &sLast, sMin );
		}
	    } /* Prepares for next cycle */

        } /* icycle */
/*____________________________________________________________________________*/

        if ( B_selectmin ) {
	    copyState( &sSave, sMin );
        } else {
	    copyState( &sSave, sLast );
        }

	sSave.Q.ang = WrpRad( ModRad( sSave.Q.ang ) );

        for (i=0; i<ntor; i++) {
            sSave.tor[i] = WrpRad( ModRad( sSave.tor[i] ) );
        }
   
        pr( logFile, UnderLine );
        pr( logFile, "\n\n\t\tFINAL DOCKED STATE\n" );
        pr( logFile,     "\t\t__________________\n\n\n" );
        pr( logFile, "Run Number %d, \n\nFinal Energy = %+.2f\n", irun1, eLast);

        pr( logFile, "Final Translation = %.2f, %.2f, %.2f\n", sSave.T.x, sSave.T.y, sSave.T.z );
        pr( logFile, "Final Quaternion Rotation Angle = %5.1f deg\n", Deg(sSave.Q.ang) );
        pr( logFile, "Final Quaternion Unit Vector = ( %+.2f, %+.2f, %+.2f )\n", sSave.Q.nx, sSave.Q.ny, sSave.Q.nz );

	copyState( &sHist[ *Addr_nconf ], sSave );

        if (ntor > 0) {
            pr( logFile, "Final Torsions:\n" );
            for (i=0; i<ntor; i++) {
		torTmp = Deg( sSave.tor[i] );
		torTmp = ModDeg( torTmp );
		torTmp = WrpDeg( torTmp );
                sHist[ *Addr_nconf ].tor[i] = sSave.tor[i];
                pr( logFile, "          %2d = %7.2f deg", i+1, torTmp);
		if ((B_isTorConstrained[i] == 1) && B_ShowTorE) {
		    pr(logFile, ", Energetic penalty = %uhd\n", US_TorE[i]);
		} else {
		    pr(logFile, "\n");
		}
            }
        }

	cnv_state_to_coords( sSave, vt, tlist, ntor, crdpdb, crd, natom );

	if (ntor > 0) {
	    eintra = eintcal( nonbondlist, e_internal, crd, type, Nnb, 
		B_calcIntElec, q1q2);
	} else {
	    eintra = 0.0;
	}
        einter = trilinterp( crd, charge, type, natom, map, inv_spacing, 
	    elec, emap, xlo, ylo, zlo);

	writePDBQ( irun, FN_ligand, FN_dpf, sml_center, sSave, ntor,
	  eintra, einter, natom, atomstuff, crd, emap, elec, charge,
	  ligand_is_inhibitor, torsFreeEnergy, outlev);

        econf[(*Addr_nconf)] = eLast;
        ++(*Addr_nconf);

    } /* Loop over runs ======================================================*/

    if ( B_writeTrj ) {
    	fclose( FP_trj );
    }
}
/* EOF */
