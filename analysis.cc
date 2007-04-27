/*

 $Id: analysis.cc,v 1.27 2007/04/27 06:01:47 garrett Exp $

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
#include <sys/types.h>
#include "constants.h"
#include "structs.h"
#include "getpdbcrds.h"
#include "stateLibrary.h"
#include "cnv_state_to_coords.h"
#include "sort_enrg.h"
#include "cluster_analysis.h"
#include "prClusterHist.h"
#include "getrms.h"
#include "calculateEnergies.h"
#include "print_rem.h"
#include "strindex.h"
#include "print_avsfld.h"
#include "printEnergies.h"
#include "analysis.h"
#include "eintcal.h"

extern FILE  *logFile;
extern int   keepresnum;
extern char  dock_param_fn[];
extern char  *programname;
extern int   true_ligand_atoms;

void analysis( int   Nnb, 
               char  atomstuff[MAX_ATOMS][MAX_CHARS], 
               Real charge[MAX_ATOMS], 
               Real abs_charge[MAX_ATOMS], 
               Real qsp_abs_charge[MAX_ATOMS], 
               Boole B_calcIntElec,
               Real clus_rms_tol, 
               Real crdpdb[MAX_ATOMS][SPACE], 

               const EnergyTables *ptr_ad_energy_tables,

               Real  map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS], 
               Real  econf[MAX_RUNS], 
               int   irunmax, 
               int   natom, 
               NonbondParam *nonbondlist, 
               int   nconf, 
               int   ntor, 
               State hist[MAX_RUNS], 
               char  smFileName[MAX_CHARS], 
               Real  sml_center[SPACE],
               Boole B_symmetry_flag, 
               int   tlist[MAX_TORS][MAX_ATOMS], 
               int   type[MAX_ATOMS], 
               Real  vt[MAX_TORS][SPACE],
               char  FN_rms_ref_crds[MAX_CHARS],
               Real  torsFreeEnergy,
               Boole B_write_all_clusmem,
               int   ligand_is_inhibitor,
               int   outlev,
			   int   ignore_inter[MAX_ATOMS],
               const Boole   B_include_1_4_interactions,
               const Real scale_1_4,

               const ParameterEntry parameterArray[MAX_MAPS],
               const Real unbound_internal_FE,

               GridMapSetInfo *info,
               Boole B_use_non_bond_cutoff,
               Boole B_have_flexible_residues,
               Boole B_rms_atoms_ligand_only

              )

{
    /* register int   imol = 0; */
    static char  filename[MAX_CHARS];
    static char  label[MAX_CHARS];
    static char  rec14[14];
    static char  rec9[9];

    static Real clu_rms[MAX_RUNS][MAX_RUNS];
    static Real crdSave[MAX_RUNS][MAX_ATOMS][SPACE];
    static Real crd[MAX_ATOMS][SPACE];
    static Real elec[MAX_ATOMS];
    static Real elec_total;
    static Real emap[MAX_ATOMS];
    static Real emap_total;
    // Real lo[3];
    static Real ref_crds[MAX_ATOMS][SPACE];
    static Real ref_rms[MAX_RUNS];
    Real torDeg = 0.;
    Real modtorDeg = 0.;
    Real MaxValue = 99.99;

    int   c = 0;
    int   c1 = 0;
    static int   cluster[MAX_RUNS][MAX_RUNS];
    int   i1=1;
    int   indpf = 0;
    static int   isort[MAX_RUNS];
    int   ncluster = 1; int   num_in_clu[MAX_RUNS];
    int   off[VECLENMAX];
    int   ref_natoms = -1;
    int   veclen = 0;
    int   kmax = 0;
    int   n_rms_atoms = 0;

    register int   i = 0;
    register int   j = 0;
    register int   k = 0;
    register int   t = 0;

    State save;


    pr( logFile, "\n\t\tCLUSTER ANALYSIS OF CONFORMATIONS\n\t\t_________________________________\n\nNumber of conformations = %d\n", nconf );

    // Initialise these arrays
    for (j = 0; j < MAX_RUNS; j++) {
        num_in_clu[j] = 0;
        isort[j] = j;
        for (i = 0; i < MAX_RUNS; i++) {
            cluster[j][i] = 0;
        }
    }

    // Set the number of atoms to cluster on
    pr( logFile, "\n");
    if (B_rms_atoms_ligand_only == TRUE) {
        // use only the ligand atoms for clustering
        n_rms_atoms = true_ligand_atoms;
        pr( logFile, "RMSD cluster analysis will be performed using the ligand atoms only (%d / %d total atoms).\n", n_rms_atoms, natom);
    } else {
        // use all the moving atoms in the receptor plus the ligand atoms for clustering
        n_rms_atoms = natom;
        pr( logFile, "RMSD cluster analysis will be performed using all the moving atoms in the receptor\n plus the ligand atoms (%d / %d total atoms).\n", n_rms_atoms, natom);
    }

    // Read in reference coordinates...
    if (strncmp(FN_rms_ref_crds,"unspecified filename",20) != 0) {
        if ((ref_natoms = getpdbcrds( FN_rms_ref_crds, ref_crds)) == -1) {
            fprintf( logFile, "%s: Problems while reading the ligand reference coordinates file \"%s\".\n", programname, FN_rms_ref_crds);
            fprintf( logFile, "%s: Will attempt to use the input ligand PDBQT coordinates as reference instead.\n", programname);
        } else if (ref_natoms != natom) {
            // intention is to compare the number of reference atoms with the number of atoms we are comparing
            // if receptor is flexible, natom will include both the receptor and ligand atoms, but if the
            // receptor is rigid, natom will be equal to true_ligand_atoms.
            pr( logFile, "%s: ERROR!  Wrong number of atoms in reference structure.\n", programname);
            pr( logFile, "%s:         The reference structure should consist of only the ligand atoms.\n", programname);
            pr( logFile, "%s:         Input ligand PDBQT structure has %d atoms, but reference structure has %d atoms.\n\n", programname, true_ligand_atoms, ref_natoms);
            ref_natoms = -1;
        }
    }

    // Generate coordinates for each final transformation,
    for ( k=0; k<nconf; k++ ) {

        /* fprintf( logFile, "\n\nState hist[%d].\n", k); */
        if (outlev > -1) {
            printState( logFile, hist[k], 2 );
        }

        /* fprintf( logFile, "\nCopying state %d.\n", k); */
        copyState( &save, hist[k] );

        /* fprintf( logFile, "Converting state %d to coordinates.\n", k); */
        cnv_state_to_coords( save, vt, tlist, ntor, crdpdb, crd, natom);

        /* fprintf( logFile, "Saving coordinates of state %d.\n", k); */

        /* Save coordinates in crdSave array...  */
        (void)memcpy(crdSave[k], crd, natom*3*sizeof(Real));
    } /*k*/

    flushLog;

    // Sort conformations by energy and perform cluster analysis,
    if (nconf > 1) {
        sort_enrg( econf, isort, nconf );

        // NOTE:  We are clustering on only the ligand atoms, regardless
        // of flexibility in the receptor sidechains...
        ncluster = cluster_analysis( clus_rms_tol, cluster, num_in_clu, isort, 
                    nconf, n_rms_atoms, type, crdSave, crdpdb, 
                    sml_center, clu_rms, B_symmetry_flag,
                    ref_crds, ref_natoms, ref_rms);

        pr( logFile, "\nOutputting structurally similar clusters, ranked in order of increasing energy.\n" );
        flushLog;

        prClusterHist( ncluster, irunmax, clus_rms_tol,num_in_clu, cluster, econf, clu_rms, ref_rms);

        if (outlev > -1) {
            pr( logFile, "\n\tLOWEST ENERGY DOCKED CONFORMATION from EACH CLUSTER");
            pr( logFile, "\n\t___________________________________________________\n\n\n" );

            if (keepresnum > 0 ) {
                pr( logFile, "\nKeeping original residue number (specified in the input PDBQ file) for outputting.\n\n");
            } else {
                pr( logFile, "\nResidue number will be set to the conformation's cluster rank.\n\n");
            }
        }
    } else {
        pr( logFile, "\nSorry!  Unable to perform cluster analysis, because not enough conformations were generated.\n\n\n" );

        ncluster = 1;
        ref_rms[0] = getrms( crd, ref_crds, B_symmetry_flag, n_rms_atoms, type);
        clu_rms[0][0] = 0.;
        num_in_clu[0] = 1;
        cluster[0][0] = 0;
    }
    flushLog;

    // For each cluster, i
    for (i = 0;  i < ncluster;  i++) {
        i1 = i + 1;

        // c = cluster[i][0];
        if (B_write_all_clusmem) {
            kmax = num_in_clu[i];
        } else {
            kmax = 1;	/* write lowest-energy only */
        }

        // For each member, k, of this cluster
        for (k = 0;  k < kmax;  k++) {
            c = cluster[i][k];
            c1 = c + 1;

            (void)memcpy(crd, crdSave[c], natom*3*sizeof(Real));

            Boole B_outside=FALSE;
            register int ia=0;
            for (ia=0; (ia<natom)&&(!B_outside); ia++) {
                B_outside = is_out_grid_info(crd[ia][0], crd[ia][1], crd[ia][2]);
            }

            EnergyBreakdown eb;

            eb = calculateEnergies( natom, ntor, unbound_internal_FE, torsFreeEnergy, B_have_flexible_residues,
                 crd, charge, abs_charge, type, map, info, B_outside?SOME_ATOMS_OUTSIDE_GRID:ALL_ATOMS_INSIDE_GRID,
                 ignore_inter, elec, emap, &elec_total, &emap_total,
                 nonbondlist, ptr_ad_energy_tables, Nnb, B_calcIntElec,
                 B_include_1_4_interactions, scale_1_4, qsp_abs_charge, parameterArray, B_use_non_bond_cutoff );
     
            print_rem( logFile, i1, num_in_clu[i], c1, ref_rms[c]);

            printEnergies( &eb, "USER    ", ligand_is_inhibitor, emap_total, elec_total, B_have_flexible_residues);
     
            pr( logFile, "USER  \n");
            pr( logFile, "USER    DPF = %s\n", dock_param_fn);
            pr( logFile, "USER    NEWDPF move\t%s\n", smFileName );
            pr( logFile, "USER    NEWDPF about\t%f %f %f\n", sml_center[X],sml_center[Y],sml_center[Z]);
            pr( logFile, "USER    NEWDPF tran0\t%f %f %f\n", hist[c].T.x, hist[c].T.y, hist[c].T.z );
            pr( logFile, "USER    NEWDPF axisangle0\t%f %f %f %f\n", hist[c].Q.nx, hist[c].Q.ny, hist[c].Q.nz, RadiansToDegrees(hist[c].Q.ang) );
            pr( logFile, "USER    NEWDPF quaternion0\t%f %f %f %f\n", hist[c].Q.x, hist[c].Q.y, hist[c].Q.z, hist[c].Q.w );
            if (ntor > 0) {
                pr( logFile, "USER    NEWDPF ndihe\t%d\n", hist[c].ntor );
                pr( logFile, "USER    NEWDPF dihe0\t" );
                flushLog;
                for ( t = 0;  t < hist[c].ntor;  t++ ) {
                    torDeg = RadiansToDegrees(hist[c].tor[t]);
                    modtorDeg = ModDeg(torDeg);
                    pr( logFile, "%.2f ", WrpDeg(modtorDeg) );
                }/*t*/
                pr( logFile, "\n" );
            }/*if*/
            pr( logFile, "USER  \n");
            flushLog;
     
            if (outlev > -11) {
                if (keepresnum > 0) {
                    // Log File PDBQ coordinates [
                    pr( logFile, "USER                              x       y       z    vdW   Elec        q     RMS \n" );
                    // TODO output the ROOT, ENDROOT, BRANCH, ENDBRANCH, TORS records...
                    for (j = 0;  j < natom;  j++) {
                        strncpy( rec14, &atomstuff[j][13], (size_t)13);
                        rec14[13]='\0';
                        pr(logFile, FORMAT_PDBQ_ATOM_RESSTR, "", j+1, rec14, crd[j][X], crd[j][Y], crd[j][Z], min(emap[j], MaxValue), min(elec[j], MaxValue), charge[j]);
                        pr(logFile," %6.3f\n", ref_rms[c]); 
                    }
                    //]
                } else {
                    // Log File PDBQ coordinates [
                    // TODO output the ROOT, ENDROOT, BRANCH, ENDBRANCH, TORS records...
                    pr( logFile, "USER                 Rank         x       y       z    vdW   Elec        q     RMS \n");
                    for (j = 0;  j < natom;  j++) {
                        strncpy( rec9, &atomstuff[j][13], (size_t)8);
                        rec9[8]='\0';
                        pr(logFile, FORMAT_PDBQ_ATOM_RESNUM, "", j+1, rec9, i1, crd[j][X], crd[j][Y], crd[j][Z], min(emap[j], MaxValue), min(elec[j], MaxValue), charge[j]);
                        pr(logFile," %6.3f\n", ref_rms[c]); 
                    }/*j*/
                    //]
                }/*if*/
            }
            pr( logFile, "TER\n" );
            pr( logFile, "ENDMDL\n" );
            // End of outputting coordinates of this "MODEL"...
#ifdef EINTCALPRINT
            // Print detailed breakdown of internal energies of all non-bonds
            (void) eintcalPrint(nonbondlist, ptr_ad_energy_tables, crd, Nnb, B_calcIntElec, B_include_1_4_interactions, scale_1_4, qsp_abs_charge, parameterArray, B_use_non_bond_cutoff, B_have_flexible_residues);
#endif
            flushLog;
        } /*k*/
    } /*i   (Next cluster.) */
    pr( logFile, "\n\n" );

    // AVS Field file [
    off[0]=5;
    off[1]=6;
    off[2]=7;
    off[3]=8;
    off[4]=9;
    off[5]=10;
    if (keepresnum > 0) {
        off[6]=11;
        veclen = 7;
        strcpy( label, "x y z vdW Elec q RMS\0" );
    } else {
        off[6]=4;
        off[7]=11;
        veclen = 8;
        strcpy( label, "x y z vdW Elec q Rank RMS\0" );
    }
    indpf = strindex( dock_param_fn, ".dpf" );
    strncpy( filename, dock_param_fn, (size_t)indpf );
    filename[ indpf ] = '\0';
    strcat( filename, ".dlg.pdb\0" );

    print_avsfld( logFile, veclen, natom, ncluster, off, 12, label, filename );
    //]
}
/* EOF */
