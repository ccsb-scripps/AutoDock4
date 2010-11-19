/*

 $Id: readPDBQT.cc,v 1.16.2.1 2010/11/19 20:09:29 rhuey Exp $

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
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <sys/times.h>
#include <sys/param.h>
#include <time.h>
#include <ctype.h>		/* tolower */
#include <search.h>
#include "readPDBQT.h"
#include "PDBQT_tokens.h"
#include "structs.h"
#include "atom_parameter_manager.h"

/*----------------------------------------------------------------------------*/

extern int      debug;
extern int      parse_tors_mode;
extern FILE    *logFile;
extern char    *programname;
extern int      true_ligand_atoms;
extern int		ntor_lig[MAX_LIGANDS];  // number of torsion in each true ligand
extern int      ntor_res[MAX_LIGANDS];
extern int		total_ntor_res;
extern int      gene_index_lig[MAX_LIGANDS][2];  //gene num start_point & end_point of a ligand.
/*----------------------------------------------------------------------------*/
Molecule readPDBQT(char input_line[LINE_LEN],
                    int num_atom_maps,
                    int *P_natom,
                    Real crdpdb[MAX_ATOMS][NTRN],
                    Real crdreo[MAX_ATOMS][NTRN],
                    Real charge[MAX_ATOMS],
                    Boole * P_B_haveCharges,
                    int map_index[MAX_ATOMS], //was:int type[MAX_ATOMS]
                    int bond_index[MAX_ATOMS],
                    char pdbaname[MAX_ATOMS][5],
                    // added for multiple ligands  - Huameng 09/27/07
					int nlig,                // number of actual ligands       
					int natom_in_lig[MAX_LIGANDS],  // atom_number in each ligand
                    char FN_ligands[MAX_LIGANDS][MAX_CHARS],     // array of ligand file names
                    char FN_flexres[MAX_CHARS],
                    Boole B_have_flexible_residues,
                    char atomstuff[MAX_ATOMS][MAX_CHARS],
                    int Htype,
                    Boole * P_B_constrain,
                    int *P_atomC1,
                    int *P_atomC2,
                    Real *P_sqlower,
                    Real *P_squpper,
                    int *P_ntor1,
                    int *P_ntor,
                    int *P_ntor_ligand,   // the number of torsions in the ligand (excluding the flexible residues in receptor)
                    int tlist[MAX_TORS][MAX_ATOMS],
                    Real vt[MAX_TORS][NTRN],
                    int *P_Nnb,
                    NonbondParam *nonbondlist,
                    Clock jobStart,
                    struct tms tms_jobStart,
                    char hostnm[MAX_CHARS],
                    int *P_ntorsdof,
                    int outlev,
                    int ignore_inter[MAX_ATOMS],
                    int B_include_1_4_interactions,
                    Atom atoms[MAX_ATOMS],
                    char PDBQT_record[MAX_RECORDS][LINE_LEN]
                    )

{
	FILE           *FP_ligand;
	FILE           *FP_flexres;
	static char     dummy[LINE_LEN];
	static char     error_message[LINE_LEN];
	static char     message[LINE_LEN];

	Real   aq = 0.;
	Real   lq = 0.;
	Real   total_charge_ligand = 0.;
	Real   uq = 0.;

	static int      atomnumber[MAX_RECORDS];
	int             iq = 0;
	int             natom = 0;
	static int      nbmatrix[MAX_ATOMS][MAX_ATOMS];
	int             nrecord = 0;
    int             nligand_record = 0;
	int             ntor = 0;
	static int      ntype[MAX_ATOMS];
	static int      rigid_piece[MAX_ATOMS];
	int             found_begin_res = 0; // found_begin_res == 0 means we have not yet found a BEGIN_RES record...
   
    int             keyword_id = -1;
	int             nres = 0;
	// added -Huameng 09/27/07
	int 			found_begin_lig = 0; // means we have not yet found a BEGIN_LIGAND record... 
	int             ilig = 0;
	int             ntor_true_ligand = 0;
	int             nrigid_piece = 0;

	Boole           B_has_conect_records = FALSE;
	Boole           B_is_in_branch = FALSE;

	int             bonded[MAX_ATOMS][6];

	register int    i = 0;
	register int    j = 0;

	static Real QTOL = 0.005;

	Molecule        mol;

	ParameterEntry  this_parameter_entry;

	for (i = 0; i < MAX_RECORDS; i++) {
		atomnumber[i] = 0;
	}

	for (j = 0; j < MAX_ATOMS; j++) {
		ntype[j] = 0;
		rigid_piece[j] = 0;
	}

    //  Attempt to open the ligand PDBQT file...
    
    // -Huameng 09/24/07  open multiple ligand files
	//sscanf(input_line, "%*s %s", FN_ligand);
	//if (openFile(FN_ligand, "r", &FP_ligand, jobStart, tms_jobStart, TRUE)) {
	//	pr(logFile, "Ligand PDBQT file = \"%s\"\n\n", FN_ligand);
	//}
	
	// Count the number of records in the Ligand PDBQT file first...
	for(int i=0; i < nlig;  i++) {
		if (openFile(FN_ligands[i], "r", &FP_ligand, jobStart, tms_jobStart, TRUE)) {
			pr(logFile, "Ligand PDBQT file=\"%s\"\n", FN_ligands[i]);
		}	
		while (fgets(dummy, LINE_LEN, FP_ligand) != NULL) {
			++nrecord;
		}
		(void) fclose(FP_ligand);
	}
	
	// now open flexible residue pdbqt file and count records
    if (B_have_flexible_residues) {
        //  Attempt to open the flexible residue PDBQT file...
        if (openFile(FN_flexres, "r", &FP_flexres, jobStart, tms_jobStart, TRUE)) {
            pr(logFile, "Flexible Residues PDBQT file=\"%s\"\n", FN_flexres);
        }
        //  Count the number of records in the Flexible Residue PDBQT file next, 
        //  but don't reset the counter...
        while (fgets(dummy, LINE_LEN, FP_flexres) != NULL) {
            ++nrecord;
        }
        (void) fclose(FP_flexres);
    } else {
    	pr(logFile, "No Flexible Residues specified in dpf file\n");
    }

    // Set the nrecord-th entry of PDBQT_record to NULL, aka '\0'
    (void) strcpy(PDBQT_record[nrecord], "\0");

    // Read in the PDBQT file(s) if there are not too many records
	if (nrecord > MAX_RECORDS) {
		prStr(error_message, "ERROR: %d records read in, but only dimensioned for %d.\nChange \"MAX_RECORDS\" in \"constants.h\".", nrecord, MAX_RECORDS);
		stop(error_message);
		exit(-1);
	} else {
        // Read in the multiple input Ligand PDBQT files...
        // -Huameng 09/24/07 
        i = 0;     
        for(int k=0; k < nlig;  k++) {
			if (openFile(FN_ligands[k], "r", &FP_ligand, jobStart, tms_jobStart, TRUE)) {
				pr(logFile,   "INPUT LIGAND PDBQT FILE %s:", FN_ligands[k]);
				pr(logFile, "\n________________________\n\n");
				// add BEGIN_LIG to PDBQT_record[0]
				sprintf(PDBQT_record[i], "BEGIN_LIG %d\n", k);
				i++;				
				pr(logFile, "%s", PDBQT_record[i]);
				while (fgets(PDBQT_record[i], LINE_LEN, FP_ligand) != NULL) {
					pr(logFile, "INPUT-LIGAND-PDBQT: %s", PDBQT_record[i]);
					// get torions for each ligand  -Huameng 11/25/2007
					if(strstr(PDBQT_record[i], "active torsions:") != NULL) {
						sscanf(PDBQT_record[i], "%*s %d", &ntor_lig[k]);
						pr(logFile, "LIGAND %d, torsions=%d\n", k, ntor_lig[k]);
						ntor_true_ligand += ntor_lig[k];
					}
					i++;
				}				
				pr(logFile, "\n%s", UnderLine);
			} // if
			(void) fclose(FP_ligand);
			
			// add END_LIG to PDBQT_record at the end of each ligand
			sprintf(PDBQT_record[i], "END_LIG %d\n", k);
			i++;					
        } // nlig, k
        
        // end reading true ligand PDBQT rec        
        nligand_record = i;
        
        // Read in the input Flexible Residues PDBQT file... 
        // continue i...
        if (B_have_flexible_residues) {       	
        	int res = 0;           
            if (openFile(FN_flexres, "r", &FP_flexres, jobStart, tms_jobStart, TRUE)) {
                pr(logFile,   "INPUT FLEXIBLE RESIDUES PDBQT FILE:");
                pr(logFile, "\n___________________________________\n\n");
                
                while(fgets(PDBQT_record[i], LINE_LEN, FP_flexres) != NULL) {
                    pr(logFile, "INPUT-FLEXRES-PDBQT: %s", PDBQT_record[i]);
                    
                    // get torions for each felx. residue  -Huameng 11/25/2007
					if(strstr(PDBQT_record[i], "active torsions:") != NULL) {
						sscanf(PDBQT_record[i], "%*s %d", &ntor_res[res]);							
						//pr(logFile, "flex residue %d, torsions=%d\n", res, ntor_res[res]);							
						total_ntor_res += ntor_res[res];
						res++;
					}   
					i++;					                    
                } // while
             
                pr(logFile, UnderLine);
            } // if
            (void) fclose(FP_flexres);
        }
        
        // have to set nrecord to i since we added BEGIN_LIG, END_LIG records
        nrecord = i;
        
        pr(logFile, "Records (ligand + res) in PDBQT_record[]=%d, nligand_record=%d, total_records=%d\n", 
        			i, nligand_record, nrecord);
        pr(logFile, "Total torsions in flexible residues = %d.\n", total_ntor_res);
        
        // handle multi-ligand -huameng 12/05/2007
        // set values for gene_index_lig[MAX_LIGANDS][2]
        int start_idx = 0;
        for(i =0; i < nlig; i++) {
        	
         	gene_index_lig[i][0] = start_idx;        	
         	// rotation & quarterion part        	
         	gene_index_lig[i][1] = start_idx + 7;
         	
         	start_idx = gene_index_lig[i][1];
         	pr(logFile, "gene_index_lig[i][0] = %d,  gene_index_lig[i][1] = %d\n", 
         				 gene_index_lig[i][0], gene_index_lig[i][1]);
         }
        
      	 //pr(logFile, "Residue torsion gene_index_lig[nlig][0] = %d,  gene_index_lig[nlig][1] = %d\n", 
         //			gene_index_lig[nlig][0], gene_index_lig[nlig][1]);
         			
	} // if (read records in PDBQT files)

	// Count the ATOMs and HETATMs; store the (x,y,z) coordinates...
	// Also, check for any BEGIN_RES records, for receptor flexibility...
	pr(logFile, "\nDetermining Atom Types and Parameters for the Moving Atoms\n");
	pr(logFile,   "__________________________________________________________\n\n");
	natom = 0;

    // Loop over all the lines in either the ligand or the "reconstructed" combined-ligand-flexible-residues file
	for (i = 0; i < nrecord; i++) {
		strncpy(input_line, PDBQT_record[i], (size_t) LINE_LEN);
		pr(logFile, "PDBQT_record [%d]: %s", i, PDBQT_record[i]);
		// Parse this line in the ligand file
		keyword_id = parse_PDBQT_line(input_line);

        switch ( keyword_id ) {
            case PDBQ_ATOM:
            case PDBQ_HETATM:
            	// added for multiple ligands -huameng 09/27/07
                // count number of atoms in each moving ligand
                if( found_begin_lig == 0 ) {
                	// Flag this as an error
                    // Incorrectly nested BRANCH/ENDBRANCH records
                    pr( logFile, "%s: ERROR: BEGIN_LIG record must be given before any ATMO/HEATM; see line %d iligand=%d.\n\n", programname, i+1, ilig);
                    pr( stderr, "%s: ERROR: BEGIN_LIG record must be given before any ATMO/HEATM; see line %d iligand=%d.\n\n", programname, i+1, ilig);
                    exit( -1 );
                }
                if ( ! B_is_in_branch ) {
                    // Flag this as an error
                    // Incorrectly nested BRANCH/ENDBRANCH records
                    pr( logFile, "%s: ERROR:  All ATOM and HETATM records must be given before any nested BRANCHes; see line %d in PDBQT file \"%s\".\n\n", programname, i+1, FN_ligands[ilig]);
                    pr( stderr, "%s: ERROR:  All ATOM and HETATM records must be given before any nested BRANCHes; see line %d in PDBQT file \"%s\".\n\n", programname, i+1, FN_ligands[ilig]);
                    exit( -1 );
                }
                
                // Check that the line is at least 78 characters long
                if (strlen(input_line) < 78) {
                    pr(logFile, "%s: FATAL ERROR: line %d is too short!\n", programname, i+1);
                    pr(logFile, "%s: FATAL ERROR: line \"%s\".\n", programname, input_line);
                    pr(stderr, "%s: FATAL ERROR: line %d is too short!\n", programname, i+1);
                    pr(stderr, "%s: FATAL ERROR: line \"%s\".\n", programname, input_line);
                    exit(-1);
                }

                ParameterEntry * found_parm;
                int serial;

                // Set up rigid_piece array by reading in the records of the PDBQT file;
                // each "rigid_piece" is a self-contained rigid entity.
                rigid_piece[natom] = nrigid_piece;

                // Read the coordinates and store them in crdpdb[],
                // read the partial atomic charge and store it in charge[],
                // and read the parameters of this atom and store them in this_parameter_entry
                // set the "autogrid_type" in this_parameter_entry
                readPDBQTLine(input_line, &serial, crdpdb[natom], &charge[natom], &this_parameter_entry);

                /*
                // Verify the serial number for this atom
                if ( serial != (natom + 1) ) {
                    pr( logFile, "%s: ERROR:  ATOM and HETATM records must be numbered sequentially from 1.  See line %d in PDBQT file \"%s\".\n\n", programname, i+1, FN_ligand);
                    pr( stderr, "%s: ERROR:  ATOM and HETATM records must be numbered sequentially from 1.  See line %d in PDBQT file \"%s\".\n\n", programname, i+1, FN_ligand);
                    exit( -1 );
                }
                */

                // Set the serial atomnumber[i] for this atom
                atomnumber[i] = natom;

                for (j = 0; j < NTRN; j++) {
                    mol.crdpdb[natom][j] = crdpdb[natom][j];
                    mol.crd[natom][j] = crdpdb[natom][j];
                    // crdreo[natom][j] = crdpdb[natom][j];
                }

                if (!found_begin_res) {
                    // Only accumulate charges on the ligand...
                    total_charge_ligand += charge[natom];
                }
                *P_B_haveCharges = TRUE;

                strncpy(atomstuff[natom], input_line, (size_t) 30);
                atomstuff[natom][30] = '\0';
                strcpy(mol.atomstr[natom], atomstuff[natom]);

                sscanf(&input_line[12], "%s", pdbaname[natom]);

                // "map_index" is used as an index into the AutoGrid "map" array to look up 
                // the correct energies in the current grid cell, thus:	map[][][][map_index[natom]]
                map_index[natom] = -1;

                // "apm_find" is the new AutoDock 4 atom typing mechanism
                found_parm = apm_find(this_parameter_entry.autogrid_type);
                if (found_parm != NULL) {
                    map_index[natom] = found_parm->map_index;
                    bond_index[natom] = found_parm->bond_index;
                    if (outlev > 0) {
                        (void) fprintf(logFile, "Found parameters for atom %d, atom type \"%s\", grid map index = %d\n",
                                   natom + 1, found_parm->autogrid_type, found_parm->map_index);
                    }
                } else {
                    // We could not find this parameter -- return an error
                    prStr(message, "\n%s: *** WARNING!  Unknown atom type \"%s\" found.  You should add parameters for it to the parameter library first! ***\n\n", programname, this_parameter_entry.autogrid_type);
                    pr_2x(stderr, logFile, message);
                }

                if (map_index[natom] == -1) {
                    prStr(message, "%s: WARNING: the atom type (%s) of atom number %d could not be found;\n\tcheck that this atom type is listed after the \"ligand_types\" keyword in the DPF,\n\tand make sure to add a \"map\" keyword to the DPF for this atom type.\n\tNote that AutoDock will use the default atom type = 1, instead.\n\n", programname, this_parameter_entry.autogrid_type, natom+1);
                    pr_2x(stderr, logFile, message);
                    map_index[natom] = 0; // we are 0-based internally, 1-based in printed output
                }

                // Increment the number of atoms having this atomtype
                ++ntype[map_index[natom]];

                // Increment the number of atoms found in PDBQT file
                ++natom;
                
                
                // added for multiple ligands -huameng 09/27/07
                // count atoms in each moving ligand only 
                // index = nligand -1 
                if(!found_begin_res) {            
                	++natom_in_lig[ilig - 1];
                }
                
                break;

            case PDBQ_ROOT:
            case PDBQ_BRANCH:
                B_is_in_branch = TRUE;
                ++nrigid_piece;
                break;

            case PDBQ_ENDROOT:
            case PDBQ_ENDBRANCH:
                B_is_in_branch = FALSE;
                break;

            case PDBQ_NULL:
            case PDBQ_REMARK:
            case PDBQ_TORS:
            case PDBQ_ENDTORS:
            case PDBQ_TORSDOF:
            case PDBQ_CONSTRAINT:                    
            case PDBQ_END_RES:
                break;
			
            case PDBQ_CONECT:
                // At least some of the atoms in the "ligand" may have their connectivity specified
                // so we could set up their bonded entries. For future versions...
                B_has_conect_records = TRUE;
                break;              
			// handle multi-ligands  -Huameng 09/26/07
			case PDBQ_END_LIG:
				found_begin_lig = 0;
				break;
			case PDBQ_BEGIN_LIG:					
				if (!found_begin_lig) {					
					// Flag that we've found a BEGIN_LIG record.
					found_begin_lig = 1;
				}
				ilig++;						
                break;
            case PDBQ_BEGIN_RES:
            	// end reading of true ligands
            	found_begin_lig = ilig;
            	
                if (!found_begin_res) {
                    // Then a flexible receptor sidechain was found in the PDBQ file.
                    // Flag that we've found a BEGIN_RES record.
                    found_begin_res = 1;
                    pr(logFile, "\nNumber of atoms in movable ligands = %d\n\n", true_ligand_atoms);
                }
                // Increment the number of residues
                nres++;
                break;

            case PDBQ_UNRECOGNIZED:
            default:
                pr(logFile, "%s: WARNING: Unrecognized PDBQT record type in line:\n", programname );
                pr(logFile, "%s: WARNING: %s\n", programname, input_line );
                break;

        } // end switch( keyword_id )

		if (!found_begin_res) {
			// No BEGIN_RES tag has been found yet.
            // Keep updating "true_ligand_atoms" until we find a "BEGIN_RES".
            // "true_ligand_atoms" is the number of atoms in the moving ligand, 
            // and excludes all atoms in the flexible sidechain residues of the receptor.
			true_ligand_atoms = natom;
        }

	} // i, next record in PDBQT file 

	pr(logFile, "\nNumber of atoms found in flexible receptor sidechains (\"residues\") =\t%d atoms\n\n", natom - true_ligand_atoms);

	pr(logFile, "Total number of atoms found = %d atoms\t ligand atoms=%d\n\n", natom, true_ligand_atoms);
	pr(logFile, "Natom in ligand %d = %d\n\n", ilig, natom_in_lig[ilig - 1]);
	pr(logFile, "Number of actual ligands read in = %d ligands\n\n", nlig);
	pr(logFile, "Number of flexible residues in the receptor =\t%d residues\n\n", nres);

	if (natom > MAX_ATOMS) {
		prStr(error_message, "ERROR: Too many atoms found (i.e. %d); maximum allowed is %d.\nChange the \"#define MAX_ATOMS\" line in \"constants.h\"\n.", natom, MAX_ATOMS);
		stop(error_message);
		exit(-1);
	} else {
		*P_natom = natom;
		mol.natom = natom;
	}

	pr(logFile, "\nSummary of number of atoms of a given atom type:\n");
	pr(logFile, "------------------------------------------------\n\n");
	for (i = 0; i < num_atom_maps; i++) {
		pr(logFile, "Number of atoms with atom type %d = %2d\n", i + 1, ntype[i]);
	}

	pr(logFile, "\n\nSummary of total charge on ligand, residues and overall:\n");
	pr(logFile, "-------------------------------------------------------\n");

    // Check total charge on ligand
	pr(logFile, "\nTotal charge on ligand                               =\t%+.3f e\n", total_charge_ligand);
	iq = (int) ((aq = fabs(total_charge_ligand)) + 0.5);
	lq = iq - QTOL;
	uq = iq + QTOL;
	if (!((aq >= lq) && (aq <= uq))) {
		prStr(message, "\n%s: *** WARNING!  Non-integral total charge (%.3f e) on ligand! ***\n\n", programname, total_charge_ligand);
		pr_2x(stderr, logFile, message);
	}

	/*
	 * Work out where the torsions are; and what they move...
	 *
	 * Also, detect which atoms we should ignore in the
	 * intermolecular energy calculation (ignore_inter[MAX_ATOMS]
	 * array)
	 */
	mkTorTree(atomnumber, PDBQT_record, nrecord, tlist, &ntor, P_ntor_ligand, FN_ligands[0], pdbaname,
              P_B_constrain, P_atomC1, P_atomC2, P_sqlower, P_squpper, P_ntorsdof, ignore_inter);

	*P_ntor = ntor;
	*P_ntor1 = ntor - 1;
	
	mol.S.ntor = ntor;
	
	// for multi-ligand even ntor = 0, we still need build a list of bonds.
	// Huameng to fix a bug 10/05/2010
	if (nlig > 0 || ntor > 0) {
        //  Create a list of internal non-bonds to be used in internal energy calculation...
		if (debug > 0) {
			pr(logFile, "Finding bonds.\n\n");
		}
		// Initialise the bonded array
        for (i = 0; i < natom; i++) {
			for (j = 0; j < 5; j++) {
				bonded[i][j] = -1;
			} // j
            bonded[i][5] = 0;
		} // i
		
		if (debug > 0) {
			printbonds(natom, bonded, "\nDEBUG:  1. BEFORE getbonds, bonded[][] array is:\n\n", 1);
		}
        // get all the bonds in multiple ligands
        // -modified for multiple ligands  -Huameng 09/26/07
		//getbonds(crdpdb, 0, true_ligand_atoms, bond_index, bonded);
		int from_atom = 0;
		int to_atom = 0;
		for(i = 0; i < nlig; i++) {
						
    		to_atom = from_atom + natom_in_lig[i];
	    
			getbonds(crdpdb, from_atom, to_atom, bond_index, bonded);			
			//set the start atom for next ligand
			from_atom = to_atom;
		}
		
		// make sure the last to_atom equals to true_ligand_atoms 
		pr(logFile, "After getbonds of ligands to_atoms=%d, true_ligand_atoms=%d\n\n", to_atom, true_ligand_atoms);
		if(to_atom != true_ligand_atoms ) {
			pr(stderr, "Error: after getbonds of multiple ligands to_atoms=%d, all_ligand_atoms=%d\n\n", to_atom, true_ligand_atoms);
			exit(-1);
		}
		
        if (B_have_flexible_residues) {
            // find all the bonds in the receptor
            getbonds(crdpdb, true_ligand_atoms, natom, bond_index, bonded);
        }
        			
		if (debug > 0) {
			printbonds(natom, bonded, "\nDEBUG:  2. AFTER getbonds, bonded[][] array is:\n\n", 0);
			pr(logFile, "Detecting all non-bonds.\n\n");
		}
		
		// get non bonds
		nonbonds(crdpdb, nbmatrix, natom, bond_index, B_include_1_4_interactions, bonded);

		if (debug > 0) {
			printbonds(natom, bonded, "\nDEBUG:  4. AFTER nonbonds, bonded[][] array is:\n\n", 0);
			pr(logFile, "Weeding out non-bonds in rigid parts of the torsion tree.\n\n");
		}
		// add natom_in_lig array -Huameng 09/28/07
		weedbonds(natom, pdbaname, rigid_piece, ntor, tlist, nbmatrix, P_Nnb, nonbondlist, outlev, map_index,
				  natom_in_lig, nlig);

		print_nonbonds(natom, pdbaname, rigid_piece, ntor, tlist, nbmatrix, *P_Nnb, nonbondlist, outlev, map_index);

        // Update the unit vectors for the torsion rotations
        if(ntor > 0){
        	update_torsion_vectors( crdpdb, ntor, tlist, vt, &mol, debug );
        }
		flushLog;

	} else {
		fprintf(logFile, ">>> No torsions detected or single ligand docking, so skipping \"nonbonds\", \"weedbonds\" and \"torNorVec\" <<<\n\n");
	}

    //  End program if just parsing torsions...
	if (parse_tors_mode) {
		prStr(message, "\n\n *** PARSE TORSIONS MODE - Stopping here ***\n\n");
		fprintf(logFile, message);
		success(hostnm, jobStart, tms_jobStart);
		exit(0);
	}
	return mol;
}


/*----------------------------------------------------------------------------*/
/* readPDBQTLine.cc */

void
readPDBQTLine( char line[LINE_LEN],
               int  *ptr_serial,
               Real crd[SPACE],
               Real *ptr_q,
               ParameterEntry *this_parameter_entry )
/*----------------------------------------------------------------------------*/
{
    char char8[9];
    char char6[7];
    char char5[6];
    char char2[3];
	static char message[LINE_LEN];

    // Initialise char5
    (void) strcpy( char5, "    0" );
    char5[5] = '\0';

    // Initialise char8
    (void) strcpy( char8, "   0.000" );
    char8[8] = '\0';

    // Initialise char6
    (void) strcpy( char6, "  0.00" );
    char6[6] = '\0';

    // Initialise char2
    (void) strcpy( char2, "C " );
    char2[2] = '\0';

#define check_sscanf( str, fmt, val, fieldname )  if (1 != sscanf( str, fmt, val ))  {\
    sprintf(message, "\n%s: WARNING! Could not read " fieldname " in PDBQT line \"%s\".\n", programname, line );\
    pr_2x(stderr, logFile, message); }

    // Read in the serial number of this atom
    (void) strncpy( char5, &line[6], (size_t)5 );
    char5[5] = '\0';
    check_sscanf( char5, "%d", ptr_serial, "serial number" );

	// Read in the X, Y, Z coordinates
    (void) strncpy( char8, &line[30], (size_t)8 );
    char8[8] = '\0';
    check_sscanf( char8, FDFMT, &crd[X], "x-coordinate" );

    (void) strncpy( char8, &line[38], (size_t)8 );
    char8[8] = '\0';
    check_sscanf( char8, FDFMT, &crd[Y], "y-coordinate" );

    (void) strncpy( char8, &line[46], (size_t)8 );
    char8[8] = '\0';
    check_sscanf( char8, FDFMT, &crd[Z], "z-coordinate" );

#ifdef DEBUG
	(void) fprintf(stderr, "readPDBQTLine: %s", line);
#endif				/* DEBUG */

    // partial charge, q
    (void) strncpy( char6, &line[70], (size_t)6 );
    char6[6] = '\0';
    check_sscanf( char6, FDFMT, ptr_q, "partial charge" );

    // atom type name
    (void) strncpy( char2, &line[77], (size_t)2 );
    char2[2] = '\0';
    check_sscanf( char2, "%s", this_parameter_entry->autogrid_type, "atom type" );

#undef check_sscanf

#ifdef DEBUG
	fprintf(stderr, "readPDBQTLine:  %d, %.3f, %.3f, %.3f, %.3f, %s\n", *ptr_serial, crd[X], crd[Y], crd[Z], *ptr_q,
    this_parameter_entry->autogrid_type);
#endif				/* DEBUG */
}

/* EOF */
