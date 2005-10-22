/*

 $Id: writeStateOfPDBQ.cc,v 1.6 2005/10/22 04:40:02 garrett Exp $

*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/* writePDBQState.cc */

#include <math.h>
#include <stdio.h>
#include <string.h>
#include "assert.h"
#include "writeStateOfPDBQ.h"

extern int      keepresnum;
extern FILE    *logFile;

extern int      write_stateFile;
extern FILE    *stateFile;

void
writeStateOfPDBQ(int irun, FourByteLong seed[2],

		 char smFileName[MAX_CHARS],
		 char dpfFN[MAX_CHARS],
		 FloatOrDouble sml_center[SPACE],
		 State * Ptr_state,
		 int ntor,
		 FloatOrDouble * Ptr_eintra,
		 FloatOrDouble * Ptr_einter,
		 int natom,
		 char atomstuff[MAX_ATOMS][MAX_CHARS],
		 FloatOrDouble crd[MAX_ATOMS][SPACE],
		 FloatOrDouble emap[MAX_ATOMS],
		 FloatOrDouble elec[MAX_ATOMS],
		 FloatOrDouble charge[MAX_ATOMS],
		 FloatOrDouble abs_charge[MAX_ATOMS],
		 FloatOrDouble qsp_abs_charge[MAX_ATOMS],
		 int ligand_is_inhibitor,
		 FloatOrDouble torsFreeEnergy,
		 FloatOrDouble vt[MAX_TORS][SPACE],
		 int tlist[MAX_TORS][MAX_ATOMS],
		 FloatOrDouble crdpdb[MAX_ATOMS][SPACE],
		 int nonbondlist[MAX_NONBONDS][MAX_NBDATA],
         EnergyTables *ptr_ad_energy_tables,
		 int type[MAX_ATOMS],
		 int Nnb,
		 Boole B_calcIntElec,
		 FloatOrDouble q1q2[MAX_NONBONDS],
         FloatOrDouble map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS],
		 Boole B_template,
		 FloatOrDouble template_energy[MAX_ATOMS],
		 FloatOrDouble template_stddev[MAX_ATOMS],
		 int outlev,
		 int ignore_inter[MAX_ATOMS],
		 const Boole B_include_1_4_interactions,
		 const FloatOrDouble scale_1_4,
		 const ParameterEntry parameterArray[MAX_MAPS],
		 const FloatOrDouble unbound_internal_FE,

         GridMapSetInfo *info)

{
	int             i = 0;

	FloatOrDouble   emap_total = 0.0L;
	FloatOrDouble   elec_total = 0.0L;
	FloatOrDouble   MaxValue = 99.99L;

	char            AtmNamResNamNum[14], AtmNamResNam[9];

	Boole           B_outside = FALSE;


	for (i = 0; i < 14; i++) {
		AtmNamResNamNum[i] = '\0';
	}

	for (i = 0; i < 9; i++) {
		AtmNamResNam[i] = '\0';
	}

	if ((outlev > -1) && (outlev < 3)) {
        //the 2 is the level of detail:2 is high, 0 is low
        printState(logFile, (*Ptr_state), outlev);
	} else if (outlev == -1) {
		printState(logFile, (*Ptr_state), 0);
	}
	cnv_state_to_coords((*Ptr_state), vt, tlist, ntor, crdpdb, crd, natom);

	B_outside = FALSE;
	for (i = 0; (i < natom) && (!B_outside); i++) {
		B_outside = is_out_grid_info(crd[i][0], crd[i][1], crd[i][2]);
	}
	if (!B_outside) {
		if (ntor > 0) {
			*Ptr_eintra = eintcal(nonbondlist, ptr_ad_energy_tables, crd, Nnb, B_calcIntElec, q1q2, B_include_1_4_interactions, scale_1_4, qsp_abs_charge, parameterArray, unbound_internal_FE);
		} else {
			*Ptr_eintra = 0.0;
		}
		if (B_template) {
			*Ptr_einter = byatom_template_trilinterp(crd, charge, abs_charge, type, natom, map, elec, emap, 
					  template_energy, template_stddev, info );
		} else {
			*Ptr_einter = trilinterp4(crd, charge, abs_charge, type, natom, map, elec, emap, ignore_inter, info );
		}

	} else {
		*Ptr_eintra = *Ptr_einter = BIG;
		//BIG is defined in constants.h
	}

	//Sum the non - bonded energies(emap[i]) and the electrostatic energies(elec[i])
	                emap_total = 0.0L;
	elec_total = 0.0L;
	for (i = 0; i < natom; i++) {
		emap_total += emap[i];
		elec_total += elec[i];
	}


	if (outlev > -1) {
		//output of coordinates(gmm 2001 - 11 - 01)
        pr( logFile, "DOCKED: MODEL     %4d\n", irun+1 );
        pr( logFile, "DOCKED: USER    Run = %d\n", irun+1 );
        pr( logFile, "DOCKED: USER    DPF = %s\n", dpfFN );
        pr( logFile, "DOCKED: USER  \n" );
        
        printEnergies(*Ptr_einter, *Ptr_eintra, torsFreeEnergy, "DOCKED: USER    ", ligand_is_inhibitor, emap_total, elec_total, unbound_internal_FE);

		if (write_stateFile) {
			pr(stateFile, "\n");
			pr(stateFile, "\t<run id=\"%4d\">\n", irun + 1);
			pr(stateFile, "\t\t<seed>%ld %ld</seed>\n", seed[0], seed[1]);
			pr(stateFile, "\t\t<dpf>%s</dpf>\n", dpfFN);
			printStateEnergies(*Ptr_einter, *Ptr_eintra, torsFreeEnergy, "DOCKED: USER    ", ligand_is_inhibitor, unbound_internal_FE);
		}
		(void) fprintf(logFile, "DOCKED: USER    NEWDPF move %s\n", smFileName);
		(void) fprintf(logFile, "DOCKED: USER    NEWDPF about %f %f %f\n", sml_center[X], sml_center[Y], sml_center[Z]);
		(void) fprintf(logFile, "DOCKED: USER    NEWDPF tran0 %f %f %f\n", (*Ptr_state).T.x, (*Ptr_state).T.y, (*Ptr_state).T.z);
		(void) fprintf(logFile, "DOCKED: USER    NEWDPF quat0 %f %f %f %f\n", (*Ptr_state).Q.nx, (*Ptr_state).Q.ny, (*Ptr_state).Q.nz, Deg(WrpRad(ModRad((*Ptr_state).Q.ang))));
		if (ntor > 0) {
			(void) fprintf(logFile, "DOCKED: USER    NEWDPF ndihe %d\n", ntor);
			(void) fprintf(logFile, "DOCKED: USER    NEWDPF dihe0 ");
			for (i = 0; i < ntor; i++) {
				(void) fprintf(logFile, "%.2f ", Deg((*Ptr_state).tor[i]));
			}
			(void) fprintf(logFile, "\n");

		}		/* write state file */
		if (write_stateFile) {
			pr(stateFile, "\t\t<move>%s</move>\n", smFileName);
			pr(stateFile, "\t\t<about>%f %f %f</about>\n", sml_center[X], sml_center[Y], sml_center[Z]);

			pr(stateFile, "\t\t<tran0>%f %f %f</tran0>\n", (*Ptr_state).T.x, (*Ptr_state).T.y, (*Ptr_state).T.z);
			pr(stateFile, "\t\t<quat0>%f %f %f %f</quat0>\n", (*Ptr_state).Q.nx, (*Ptr_state).Q.ny, (*Ptr_state).Q.nz, Deg(WrpRad(ModRad((*Ptr_state).Q.ang))));
			if (ntor > 0) {
				pr(stateFile, "\t\t<ndihe>%d</ndihe>\n", ntor);
				pr(stateFile, "\t\t<dihe0>");
				for (i = 0; i < ntor; i++) {
					(void) fprintf(stateFile, "%.2f ", Deg(WrpRad(ModRad((*Ptr_state).tor[i]))));
				}
				(void) fprintf(stateFile, "\n");
				pr(stateFile, "</dihe0>\n");
			}
			pr(stateFile, "\t</run>\n");
		}		/* end write state file */
		(void) fprintf(logFile, "DOCKED: USER                              x       y       z     vdW  Elec       q\n");
		(void) fprintf(logFile, "DOCKED: USER                           _______ _______ _______ _____ _____    ______\n");
		if (keepresnum > 0) {
			for (i = 0; i < natom; i++) {
				assert(i >= 0 && i < natom);
				strncpy(AtmNamResNamNum, &atomstuff[i][13], (size_t) 13);
				AtmNamResNamNum[13] = '\0';
				(void) fprintf(logFile, FORMAT_PDBQ_ATOM_RESSTR, "DOCKED: ", i + 1, AtmNamResNamNum, crd[i][X], crd[i][Y], crd[i][Z], min(emap[i], MaxValue), min(elec[i], MaxValue), charge[i]);
				(void) fprintf(logFile, "\n");
			}
		} else {
			for (i = 0; i < natom; i++) {
				assert(i >= 0 && i < natom);
				strncpy(AtmNamResNam, &atomstuff[i][13], (size_t) 8);
				AtmNamResNam[8] = '\0';
				(void) fprintf(logFile, FORMAT_PDBQ_ATOM_RESNUM, "DOCKED: ", i + 1, AtmNamResNam, irun + 1, crd[i][X], crd[i][Y], crd[i][Z], min(emap[i], MaxValue), min(elec[i], MaxValue), charge[i]);
				(void) fprintf(logFile, "\n");
			}
		}
		(void) fprintf(logFile, "DOCKED: TER\n");
		(void) fprintf(logFile, "DOCKED: ENDMDL\n");
		//(void) fprintf(logFile, UnderLine);
		(void) fflush(logFile);
	}
}

void write_emap_elec( State *Ptr_state,
                     int ntor,
                     int natom,
                     FloatOrDouble crd[MAX_ATOMS][SPACE],
                     FloatOrDouble charge[MAX_ATOMS],
                     FloatOrDouble abs_charge[MAX_ATOMS],
                     FloatOrDouble qsp_abs_charge[MAX_ATOMS],
                     FloatOrDouble vt[MAX_TORS][SPACE],
                     int tlist[MAX_TORS][MAX_ATOMS],
                     FloatOrDouble crdpdb[MAX_ATOMS][SPACE],
                     int nonbondlist[MAX_NONBONDS][MAX_NBDATA],
                     EnergyTables *ptr_ad_energy_tables,
                     int type[MAX_ATOMS],
                     int Nnb,
                     Boole B_calcIntElec,
                     FloatOrDouble q1q2[MAX_NONBONDS],
                     FloatOrDouble map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS],
                     Boole B_template,
                     FloatOrDouble template_energy[MAX_ATOMS],
                     FloatOrDouble template_stddev[MAX_ATOMS],
                     int ignore_inter[MAX_ATOMS],
                     const Boole B_include_1_4_interactions,
                     const FloatOrDouble scale_1_4,
                     const ParameterEntry parameterArray[MAX_MAPS],
                     const FloatOrDouble unbound_internal_FE,
                     GridMapSetInfo *info )

{
	int           i = 0;
	FloatOrDouble emap_total = 0.0L;
	FloatOrDouble elec_total = 0.0L;
	FloatOrDouble emap[MAX_ATOMS];
	FloatOrDouble elec[MAX_ATOMS];


	cnv_state_to_coords( (*Ptr_state), vt, tlist, ntor, crdpdb, crd, natom);

	Boole B_outside = FALSE;
	for (i = 0; (i < natom) && (!B_outside); i++) {
		B_outside = is_out_grid_info(crd[i][0], crd[i][1], crd[i][2]);
	}
	if (!B_outside) {
		if (ntor > 0) {
			(void) eintcal(nonbondlist, ptr_ad_energy_tables, crd, Nnb, B_calcIntElec, q1q2, 
                           B_include_1_4_interactions, scale_1_4, qsp_abs_charge, 
                           parameterArray, unbound_internal_FE);
		}
		if (B_template) {
			(void) byatom_template_trilinterp(crd, charge, abs_charge, type, natom, map, 
                                              elec, emap, template_energy, template_stddev, info );
		} else {
			(void) trilinterp4(crd, charge, abs_charge, type, natom, map, elec, emap, ignore_inter, info );
		}
	}

	// Sum the non-bonded energies (emap[i]) and the electrostatic energies (elec[i])
    emap_total = 0.0L;
	elec_total = 0.0L;
	for (i = 0; i < natom; i++) {
		emap_total += emap[i];
		elec_total += elec[i];
	}
    pr( logFile, "%9.2e %9.2e", emap_total, elec_total );
}
/* EOF */
