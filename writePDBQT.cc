/*

 $Id: writePDBQT.cc,v 1.14 2007/04/27 06:01:52 garrett Exp $

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
#include <stdio.h>
#include <string.h>
#include "assert.h"
#include "writePDBQT.h"
#include "parse_PDBQT_line.h"
#include "calculateEnergies.h"

extern int keepresnum;
extern FILE *logFile;
extern int write_stateFile;
extern FILE *stateFile;
extern int true_ligand_atoms;
extern int Nnb_array[3];
extern Real nb_group_energy[3];

void
writePDBQT(int irun, FourByteLong seed[2],

		 char smFileName[MAX_CHARS],
		 char dpfFN[MAX_CHARS],
		 Real sml_center[SPACE],
		 State state,
		 int ntor,
		 Real * Ptr_eintra,
		 Real * Ptr_einter,
		 int natom,
		 char atomstuff[MAX_ATOMS][MAX_CHARS],
		 Real crd[MAX_ATOMS][SPACE],
		 Real emap[MAX_ATOMS],
		 Real elec[MAX_ATOMS],
		 Real charge[MAX_ATOMS],
		 Real abs_charge[MAX_ATOMS],
		 Real qsp_abs_charge[MAX_ATOMS],
		 int ligand_is_inhibitor,
		 Real torsFreeEnergy,
		 Real vt[MAX_TORS][SPACE],
		 int tlist[MAX_TORS][MAX_ATOMS],
		 Real crdpdb[MAX_ATOMS][SPACE],
		 NonbondParam *nonbondlist,
         EnergyTables *ptr_ad_energy_tables,
		 int type[MAX_ATOMS],  // aka 'map_index' in 'ParameterEntry' structures
		 int Nnb,
		 Boole B_calcIntElec,
         Real map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS],
		 int outlev,
		 int ignore_inter[MAX_ATOMS],
		 const Boole B_include_1_4_interactions,
		 const Real scale_1_4,
		 const ParameterEntry parameterArray[MAX_MAPS],
		 const Real unbound_internal_FE,

         GridMapSetInfo *info,
         int state_type,  // 0 means the state is unbound, 1 means the state is docked
         char PDBQT_record[MAX_RECORDS][LINE_LEN],
         Boole B_use_non_bond_cutoff,
         Boole B_have_flexible_residues
         )

{
	int i = 0;
    EnergyBreakdown eb;

	Real emap_total = 0.0L;
	Real elec_total = 0.0L;
	Real MaxValue = 99.99L;
	Real MinValue = -99.99L;

    Real e_inter_moving_fixed      = 0.0L;  // (1)  // trilinterp( 0, true_ligand_atoms, ...)
    Real e_intra_moving_fixed_rec  = 0.0L;  // (2)  // trilinterp( true_ligand_atoms, natom, ...)
    Real e_intra_moving_moving_lig = 0.0L;  // (3)  // eintcal( 0, nb_array[0], ...)            // nb_group_energy[INTRA_LIGAND]
    Real e_inter_moving_moving     = 0.0L;  // (4)  // eintcal( nb_array[0], nb_array[1], ...)  // nb_group_energy[INTER]
    Real e_intra_moving_moving_rec = 0.0L;  // (5)  // eintcal( nb_array[1], nb_array[2], ...)  // nb_group_energy[INTRA_RECEPTOR]

    Real e_inter = 0.0;      // total    intermolecular energy = (1) + (4)
    Real e_intra_lig = 0.0;  // ligand   intramolecular energy = (3)
    Real e_intra_rec = 0.0;  // receptor intramolecular energy = (2) + (5)

    e_inter     = e_inter_moving_fixed + e_inter_moving_moving;          // total    intermolecular energy = (1) + (4)
    e_intra_lig = e_intra_moving_moving_lig;                             // ligand   intramolecular energy = (3)
    e_intra_rec = e_intra_moving_fixed_rec + e_intra_moving_moving_rec;  // receptor intramolecular energy = (2) + (5)

	char AtmNamResNamNum[14], AtmNamResNam[9];
    char state_type_string[MAX_CHARS];
    char state_type_prefix_string[MAX_CHARS];
    char state_type_prefix_USER_string[MAX_CHARS];
	Boole B_outside = FALSE;
    Real this_emap = 0.;
    Real this_elec = 0.;

    // Initialise various character strings
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

    initialise_energy_breakdown( &eb, torsFreeEnergy, unbound_internal_FE );

    // Write out the state variables
	if ((outlev > -1) && (outlev < 3)) {
        // "outlev" is the level of detail: 2 is high, 0 is low
        pr(logFile,"State:\t");
        printState(logFile, state, outlev);
        pr(logFile,"\n\n");
	} else if (outlev < 0) {
		printState(logFile, state, 0);
	}

    // Convert state variables to x,y,z-coordinates
	cnv_state_to_coords( state, vt, tlist, ntor, crdpdb, crd, natom );

    // Check if any atoms are outside the grid box
    B_outside = FALSE;
    for (i = 0; (i < natom) && (!B_outside); i++) {
        B_outside = is_out_grid_info(crd[i][0], crd[i][1], crd[i][2]);
    }

    // Calculate the energy breakdown
    eb = calculateEnergies( natom, ntor, unbound_internal_FE, torsFreeEnergy, B_have_flexible_residues,
         crd, charge, abs_charge, type, map, info, B_outside, 
         ignore_inter, elec, emap, &elec_total, &emap_total,
         nonbondlist, ptr_ad_energy_tables, Nnb, B_calcIntElec,
         B_include_1_4_interactions, scale_1_4, qsp_abs_charge, parameterArray, B_use_non_bond_cutoff );

    // Set the total intramolecular energy (sum of intramolecular energies of ligand and of protein)
    if (ntor > 0) {
        // Add the intramolecular energy of the receptor, for the (moving, fixed) atom pairs // (2)
        *Ptr_eintra = nb_group_energy[INTRA_LIGAND] + nb_group_energy[INTRA_RECEPTOR] + eb.e_intra_moving_fixed_rec;
    } else {
        *Ptr_eintra = 0.0;
    }

    // Set the total intermolecular energy
    if (state_type == 1) {
        // DOCKED
        // Set *Ptr_einter, the intermolecular energy, only for DOCKED states, not for UNBOUND states
        *Ptr_einter = eb.e_inter;
    } else {
        // UNBOUND
        // "intermolecular" energy is meaningless for unbound state, so set this to zero
        *Ptr_einter = 0.0;
        eb.e_inter = 0.0;
        emap_total = 0.0;
        elec_total = 0.0;
        eb.e_inter_moving_fixed = 0.0;
        eb.e_inter_moving_moving = 0.0;
    }

	if (outlev > -1) {
		// output of coordinates
        pr( logFile, "%s: MODEL     %4d\n", state_type_string, irun+1 );
        pr( logFile, "%s: USER    Run = %d\n", state_type_string, irun+1 );
        pr( logFile, "%s: USER    DPF = %s\n", state_type_string, dpfFN );
        pr( logFile, "%s: USER  \n", state_type_string );
        
        printEnergies( &eb, state_type_prefix_USER_string, ligand_is_inhibitor, emap_total, elec_total, B_have_flexible_residues );

        // Write part of the "XML" state file
		if (write_stateFile) {
			pr(stateFile, "\n");
			pr(stateFile, "\t<run id=\"%4d\">\n", irun + 1);
			pr(stateFile, "\t\t<seed>%ld %ld</seed>\n", seed[0], seed[1]);
			pr(stateFile, "\t\t<dpf>%s</dpf>\n", dpfFN);
            printStateEnergies( &eb, state_type_prefix_USER_string, ligand_is_inhibitor );
		} // End write state file

		(void) fprintf(logFile, "%s: USER    NEWDPF move %s\n", state_type_string, smFileName);
		(void) fprintf(logFile, "%s: USER    NEWDPF about %f %f %f\n", state_type_string, sml_center[X], sml_center[Y], sml_center[Z]);
		(void) fprintf(logFile, "%s: USER    NEWDPF tran0 %f %f %f\n", state_type_string, state.T.x, state.T.y, state.T.z);
		(void) fprintf(logFile, "%s: USER    NEWDPF quaternion0 %f %f %f %f\n", state_type_string, state.Q.x, state.Q.y, state.Q.z, state.Q.w);
        state.Q = convertQuatToRot( state.Q );
		(void) fprintf(logFile, "%s: USER    NEWDPF axisangle0 %f %f %f %f\n", state_type_string, state.Q.nx, state.Q.ny, state.Q.nz, RadiansToDegrees(WrpRad(ModRad(state.Q.ang))));
		(void) fprintf(logFile, "%s: USER    NEWDPF quat0 %f %f %f %f # deprecated\n", state_type_string, state.Q.nx, state.Q.ny, state.Q.nz, RadiansToDegrees(WrpRad(ModRad(state.Q.ang))));
		if (ntor > 0) {
            // ndihe is deprecated; uses the number of torsions in the PDBQT's torsion tree
			// (void) fprintf(logFile, "%s: USER    NEWDPF ndihe %d\n", state_type_string, ntor);
			(void) fprintf(logFile, "%s: USER    NEWDPF dihe0 ", state_type_string);
			for (i = 0; i < ntor; i++) {
				(void) fprintf(logFile, "%.2f ", RadiansToDegrees(WrpRad(ModRad(state.tor[i]))));
			}
			(void) fprintf(logFile, "\n");

		}
        
        // Write remaining part of the "XML" state file
		if (write_stateFile) {
			pr(stateFile, "\t\t<move>%s</move>\n", smFileName);
			pr(stateFile, "\t\t<about>%f %f %f</about>\n", sml_center[X], sml_center[Y], sml_center[Z]);

			pr(stateFile, "\t\t<tran0>%f %f %f</tran0>\n", state.T.x, state.T.y, state.T.z);
			pr(stateFile, "\t\t<quaternion0>%f %f %f %f</quaternion0>\n", state.Q.x, state.Q.y, state.Q.z, state.Q.w);
            state.Q = convertQuatToRot( state.Q );
			pr(stateFile, "\t\t<quat0>%f %f %f %f</quat0>\n", state.Q.nx, state.Q.ny, state.Q.nz, RadiansToDegrees(WrpRad(ModRad(state.Q.ang))));
			pr(stateFile, "\t\t<axisangle0>%f %f %f %f</axisangle0>\n", state.Q.nx, state.Q.ny, state.Q.nz, RadiansToDegrees(WrpRad(ModRad(state.Q.ang))));
			if (ntor > 0) {
				pr(stateFile, "\t\t<ndihe>%d</ndihe>\n", ntor);
				pr(stateFile, "\t\t<dihe0>");
				for (i = 0; i < ntor; i++) {
					(void) fprintf(stateFile, "%.2f ", RadiansToDegrees(WrpRad(ModRad(state.tor[i]))));
				}
				(void) fprintf(stateFile, "\n");
				pr(stateFile, "</dihe0>\n");
			}
			pr(stateFile, "\t</run>\n");
		} // End write state file

        (void) fprintf(logFile, "%s: USER  \n", state_type_string);

        // Count the number of non-NULL records in the PDBQT file
        int nrecord = 0;
        int r = 0;
        for (r = 0; PDBQT_record[r][0] != '\0'; r++) { }
        nrecord = r;

        int keyword_id = -1;
        int print_header = FALSE;
        // Zero the atom counter,
        i = 0;
        for (r = 0; r < nrecord; r++) {
            // If this record is neither an ATOM nor a HETATM then print it,
            // else print the new coordinates of this atom.
            keyword_id = parse_PDBQT_line(PDBQT_record[r]);
            if (keyword_id == PDBQ_ROOT) {
                // Print the header just before we print out the ROOT record
                print_header = TRUE;
            }
            if ((keyword_id == PDBQ_ATOM) || (keyword_id == PDBQ_HETATM)) {
                assert(i >= 0 && i < natom);
                // If the state_type is unbound, then ignore the per-atom intermolecular
                // emap and elec values; set these to 0.
                if (state_type == 1) {
                    // DOCKED
                    this_emap = (emap[i] >= 0.) ? min(emap[i], MaxValue) : max(emap[i], MinValue);
                    this_elec = (elec[i] >= 0.) ? min(elec[i], MaxValue) : max(elec[i], MinValue);
                } else {
                    // UNBOUND
                    this_emap = 0.;
                    this_elec = 0.;
                }
                if (keepresnum > 0) {
                    // Retain the original Residue Numbering
                    strncpy(AtmNamResNamNum, &atomstuff[i][13], (size_t) 13);
                    AtmNamResNamNum[13] = '\0';
                    (void) fprintf(logFile, FORMAT_PDBQT_ATOM_RESSTR, state_type_prefix_string, 
                                   i + 1, AtmNamResNamNum, crd[i][X], crd[i][Y], crd[i][Z], 
                                   this_emap, this_elec,
                                   charge[i], parameterArray[type[i]].autogrid_type );
                } else {
                    // Change the residue number to the run number
                    strncpy(AtmNamResNam, &atomstuff[i][13], (size_t) 8);
                    AtmNamResNam[8] = '\0';
                    (void) fprintf(logFile, FORMAT_PDBQT_ATOM_RESNUM, state_type_prefix_string, 
                                   i + 1, AtmNamResNam, irun + 1, crd[i][X], crd[i][Y], crd[i][Z], 
                                   this_emap, this_elec,
                                   charge[i], parameterArray[type[i]].autogrid_type);
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
    } // outlev > -1
} // writePDBQT()

void print_PDBQT( FILE *logFile, 
                  const int true_ligand_atoms,
                  const char atomstuff[MAX_ATOMS][MAX_CHARS],
                  const Real crdpdb[MAX_ATOMS][SPACE],
                  const Real charge[MAX_ATOMS],
                  const ParameterEntry parameterArray[MAX_MAPS],
                  const int type[MAX_ATOMS],
                  const char prefix[MAX_CHARS] )
{ // Print out the coordinates
    register int i=0;
    char AtmNamResNamNum[14];
    for (i=0; i<true_ligand_atoms; i++) {
        strncpy( AtmNamResNamNum, &atomstuff[i][13], (size_t) 13 );
        AtmNamResNamNum[13] = '\0';
        (void) fprintf( logFile, FORMAT_PDBQT_ATOM_RESSTR, prefix, 
                        i + 1, AtmNamResNamNum, crdpdb[i][X], crdpdb[i][Y], crdpdb[i][Z], 
                        1., 0.,
                        charge[i], parameterArray[type[i]].autogrid_type );
        (void) fprintf( logFile, "\n" ); 
    }
    pr( logFile, "\n\n" );
} // end Print out the coordinates

/* EOF */
