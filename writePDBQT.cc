/*

 $Id: writePDBQT.cc,v 1.1 2006/04/12 03:46:36 garrett Exp $

*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>
#include <stdio.h>
#include <string.h>
#include "assert.h"
#include "writePDBQT.h"
#include "parse_PDBQT_line.h"

extern int keepresnum;
extern FILE *logFile;
extern int write_stateFile;
extern FILE *stateFile;

void
writePDBQT(int irun, FourByteLong seed[2],

		 char smFileName[MAX_CHARS],
		 char dpfFN[MAX_CHARS],
		 FloatOrDouble sml_center[SPACE],
		 State state,
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
		 int **nonbondlist,
         EnergyTables *ptr_ad_energy_tables,
		 int type[MAX_ATOMS],  // aka 'map_index' in 'ParameterEntry' structures
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

         GridMapSetInfo *info,
         int state_type,  // 0 means the state is unbound, 1 means the state is docked
         char PDBQT_record[MAX_RECORDS][LINE_LEN]
         )

{
	int             i = 0;
	FloatOrDouble   emap_total = 0.0L;
	FloatOrDouble   elec_total = 0.0L;
	FloatOrDouble   MaxValue = 99.99L;
	char            AtmNamResNamNum[14], AtmNamResNam[9];
    char            state_type_string[MAX_CHARS];
    char            state_type_prefix_string[MAX_CHARS];
    char            state_type_prefix_USER_string[MAX_CHARS];
	Boole           B_outside = FALSE;

    // Initialize various character strings
    if (state_type == 0) {
        strcpy(state_type_string, "UNBOUND");
        strcpy(state_type_prefix_string, "UNBOUND: ");
        strcpy(state_type_prefix_USER_string, "UNBOUND: USER    ");
    } else if (state_type == 1) {
        strcpy(state_type_string, "DOCKED");
        strcpy(state_type_prefix_string, "DOCKED: ");
        strcpy(state_type_prefix_USER_string, "DOCKED: USER    ");
    }
	for (i = 0; i < 14; i++) { AtmNamResNamNum[i] = '\0'; }
	for (i = 0; i < 9; i++) { AtmNamResNam[i] = '\0'; }

    // Write out the state variables
	if ((outlev > -1) && (outlev < 3)) {
        // "outlev" is the level of detail: 2 is high, 0 is low
        pr(logFile,"State:\t");
        printState(logFile, state, outlev);
        pr(logFile,"\n\n");
	} else if (outlev == -1) {
		printState(logFile, state, 0);
	}

    // Convert state variables to x,y,z-coordinates
	cnv_state_to_coords(state, vt, tlist, ntor, crdpdb, crd, natom);

    // Calculate the internal energy
    if (ntor > 0) {
        *Ptr_eintra = eintcal(nonbondlist, ptr_ad_energy_tables, crd, Nnb, B_calcIntElec, q1q2, B_include_1_4_interactions, scale_1_4, qsp_abs_charge, parameterArray) - unbound_internal_FE;
    } else {
        *Ptr_eintra = 0.0;
    }

    // Only for DOCKED states, not for UNBOUND states
    if (state_type == 1) {
        // Calculate the intermolecular energy
        B_outside = FALSE;
        for (i = 0; (i < natom) && (!B_outside); i++) {
            B_outside = is_out_grid_info(crd[i][0], crd[i][1], crd[i][2]);
        }
        if (!B_outside) {
            if (!B_template) {
                *Ptr_einter = trilinterp( crd, charge, abs_charge, type, natom, map, 
                        info, ALL_ATOMS_INSIDE_GRID, ignore_inter, elec, emap, &elec_total, &emap_total);
            } else {
                *Ptr_einter = template_trilinterp( crd, charge, abs_charge, type, natom, map, 
                        info, ALL_ATOMS_INSIDE_GRID, 
                        NULL_IGNORE_INTERMOL, template_energy, template_stddev, elec /* set */ , emap /* set */,  &elec_total, &emap_total);
            }
        } else {
            *Ptr_einter = BIG; // BIG is defined in constants.h
        }
    }

	if (outlev > -1) {
		// output of coordinates
        pr( logFile, "%s: MODEL     %4d\n", state_type_string, irun+1 );
        pr( logFile, "%s: USER    Run = %d\n", state_type_string, irun+1 );
        pr( logFile, "%s: USER    DPF = %s\n", state_type_string, dpfFN );
        pr( logFile, "%s: USER  \n", state_type_string );
        
        printEnergies(*Ptr_einter, *Ptr_eintra, torsFreeEnergy, state_type_prefix_USER_string, ligand_is_inhibitor, emap_total, elec_total, unbound_internal_FE);

        // Write part of the "XML" state file
		if (write_stateFile) {
			pr(stateFile, "\n");
			pr(stateFile, "\t<run id=\"%4d\">\n", irun + 1);
			pr(stateFile, "\t\t<seed>%ld %ld</seed>\n", seed[0], seed[1]);
			pr(stateFile, "\t\t<dpf>%s</dpf>\n", dpfFN);
            printStateEnergies(*Ptr_einter, *Ptr_eintra, torsFreeEnergy, state_type_prefix_USER_string, ligand_is_inhibitor, unbound_internal_FE);
		} // End write state file

		(void) fprintf(logFile, "%s: USER    NEWDPF move %s\n", state_type_string, smFileName);
		(void) fprintf(logFile, "%s: USER    NEWDPF about %f %f %f\n", state_type_string, sml_center[X], sml_center[Y], sml_center[Z]);
		(void) fprintf(logFile, "%s: USER    NEWDPF tran0 %f %f %f\n", state_type_string, state.T.x, state.T.y, state.T.z);
		(void) fprintf(logFile, "%s: USER    NEWDPF quat0 %f %f %f %f\n", state_type_string, state.Q.nx, state.Q.ny, state.Q.nz, Deg(WrpRad(ModRad(state.Q.ang))));
		if (ntor > 0) {
			(void) fprintf(logFile, "%s: USER    NEWDPF ndihe %d\n", state_type_string, ntor);
			(void) fprintf(logFile, "%s: USER    NEWDPF dihe0 ", state_type_string);
			for (i = 0; i < ntor; i++) {
				(void) fprintf(logFile, "%.2f ", Deg(state.tor[i]));
			}
			(void) fprintf(logFile, "\n");

		}
        
        // Write remaining part of the "XML" state file
		if (write_stateFile) {
			pr(stateFile, "\t\t<move>%s</move>\n", smFileName);
			pr(stateFile, "\t\t<about>%f %f %f</about>\n", sml_center[X], sml_center[Y], sml_center[Z]);

			pr(stateFile, "\t\t<tran0>%f %f %f</tran0>\n", state.T.x, state.T.y, state.T.z);
			pr(stateFile, "\t\t<quat0>%f %f %f %f</quat0>\n", state.Q.nx, state.Q.ny, state.Q.nz, Deg(WrpRad(ModRad(state.Q.ang))));
			if (ntor > 0) {
				pr(stateFile, "\t\t<ndihe>%d</ndihe>\n", ntor);
				pr(stateFile, "\t\t<dihe0>");
				for (i = 0; i < ntor; i++) {
					(void) fprintf(stateFile, "%.2f ", Deg(WrpRad(ModRad(state.tor[i]))));
				}
				(void) fprintf(stateFile, "\n");
				pr(stateFile, "</dihe0>\n");
			}
			pr(stateFile, "\t</run>\n");
		} // End write state file

        (void) fprintf(logFile, "%s: USER  \n", state_type_string);

        // Count the number of non-NULL records in the PDBQT file
        int nrecord=0;
        int r=0;
        for (r=0; PDBQT_record[r][0] != '\0'; r++) { }
        nrecord=r;

        int keyword_id = -1;
        int print_header = FALSE;
        // Zero the atom counter,
        i=0;
        for (r=0; r<nrecord; r++) {
            // If this record is neither an ATOM nor a HETATM then print it,
            // else print the new coordinates of this atom.
            keyword_id = parse_PDBQT_line(PDBQT_record[r]);
            if (keyword_id == PDBQ_ROOT) {
                // Print the header just before we print out the ROOT record
                print_header = TRUE;
            }
            if ((keyword_id == PDBQ_ATOM) || (keyword_id == PDBQ_HETATM)) {
                assert(i >= 0 && i < natom);
                if (keepresnum > 0) {
                    // Retain the original Residue Numbering
                    strncpy(AtmNamResNamNum, &atomstuff[i][13], (size_t) 13);
                    AtmNamResNamNum[13] = '\0';
                    (void) fprintf(logFile, FORMAT_PDBQT_ATOM_RESSTR, state_type_prefix_string, i + 1, AtmNamResNamNum, crd[i][X], crd[i][Y], crd[i][Z], min(emap[i], MaxValue), min(elec[i], MaxValue), charge[i], parameterArray[type[i]].autogrid_type );
                } else {
                    // Change the residue number to the run number
                    strncpy(AtmNamResNam, &atomstuff[i][13], (size_t) 8);
                    AtmNamResNam[8] = '\0';
                    (void) fprintf(logFile, FORMAT_PDBQT_ATOM_RESNUM, state_type_prefix_string, i + 1, AtmNamResNam, irun + 1, crd[i][X], crd[i][Y], crd[i][Z], min(emap[i], MaxValue), min(elec[i], MaxValue), charge[i], parameterArray[type[i]].autogrid_type);
                }
                (void) fprintf(logFile, "\n");
                // Increment the atom counter
                i++;
            } else {
                if (print_header) {
                    (void) fprintf(logFile, "%s: USER                              x       y       z     vdW  Elec       q    Type\n", state_type_string);
                    (void) fprintf(logFile, "%s: USER                           _______ _______ _______ _____ _____    ______ ____\n", state_type_string);
                    // Make sure we only print the header once
                    print_header = FALSE;
                }
                (void) fprintf(logFile, "%s%s", state_type_prefix_string, PDBQT_record[r]);
            }
        } // r

		(void) fprintf(logFile, "%s: TER\n", state_type_string);
		(void) fprintf(logFile, "%s: ENDMDL\n", state_type_string);
		//(void) fprintf(logFile, UnderLine);
		(void) fflush(logFile);
	}
}

/* EOF */
