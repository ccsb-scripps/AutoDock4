/*

 $Id: analysis.cc,v 1.15 2006/02/11 04:49:53 garrett Exp $

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
#include "eintcal.h"
#include "trilinterp.h"
#include "print_rem.h"
#include "strindex.h"
#include "print_avsfld.h"
#include "printEnergies.h"
#include "analysis.h"

extern FILE *logFile;
extern int   keepresnum;
extern char  dock_param_fn[];
extern char  *programname;

void analysis( int   Nnb, 
               char  atomstuff[MAX_ATOMS][MAX_CHARS], 
               FloatOrDouble charge[MAX_ATOMS], 
               FloatOrDouble abs_charge[MAX_ATOMS], 
               FloatOrDouble qsp_abs_charge[MAX_ATOMS], 
               Boole B_calcIntElec,
               FloatOrDouble q1q2[MAX_NONBONDS],
               FloatOrDouble clus_rms_tol, 
               FloatOrDouble crdpdb[MAX_ATOMS][SPACE], 

               EnergyTables *ptr_ad_energy_tables,

               FloatOrDouble map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS], 
               FloatOrDouble econf[MAX_RUNS], 
               int   irunmax, 
               int   natom, 
               int   nonbondlist[MAX_NONBONDS][MAX_NBDATA], 
               int   nconf, 
               int   ntor, 
               State hist[MAX_RUNS], 
               char  smFileName[MAX_CHARS], 
               FloatOrDouble sml_center[SPACE],
               Boole B_symmetry_flag, 
               int   tlist[MAX_TORS][MAX_ATOMS], 
               int   type[MAX_ATOMS], 
               FloatOrDouble vt[MAX_TORS][SPACE],
               char  FN_rms_ref_crds[MAX_CHARS],
               FloatOrDouble torsFreeEnergy,
               Boole B_write_all_clusmem,
               int ligand_is_inhibitor,
               Boole B_template,
               FloatOrDouble template_energy[MAX_ATOMS],
               FloatOrDouble template_stddev[MAX_ATOMS],
               int           outlev,
			   int           ignore_inter[MAX_ATOMS],
               const Boole   B_include_1_4_interactions,
               const FloatOrDouble scale_1_4,

               const ParameterEntry parameterArray[MAX_MAPS],
               const FloatOrDouble unbound_internal_FE,

               GridMapSetInfo *info
              )

{
    /* register int   imol = 0; */
    static char  filename[MAX_CHARS];
    static char  label[MAX_CHARS];
    static char  rec14[14];
    static char  rec9[9];

    static FloatOrDouble clu_rms[MAX_RUNS][MAX_RUNS];
    static FloatOrDouble crdSave[MAX_RUNS][MAX_ATOMS][SPACE];
    static FloatOrDouble crd[MAX_ATOMS][SPACE];
    FloatOrDouble einter = 0.;
    FloatOrDouble eintra = 0.;
    static FloatOrDouble elec[MAX_ATOMS];
    static FloatOrDouble elec_total;
    static FloatOrDouble emap[MAX_ATOMS];
    static FloatOrDouble emap_total;
    // FloatOrDouble lo[3];
    static FloatOrDouble ref_crds[MAX_ATOMS][SPACE];
    static FloatOrDouble ref_rms[MAX_RUNS];
    FloatOrDouble torDeg = 0.;
    FloatOrDouble modtorDeg = 0.;
    FloatOrDouble MaxValue = 99.99;

    int   c = 0;
    int   c1 = 0;
    static int   cluster[MAX_RUNS][MAX_RUNS];
    int   i1=1;
    int   indpf = 0;
    static int   isort[MAX_RUNS];
    int   ncluster = 1;
    int   num_in_clu[MAX_RUNS];
    int   off[VECLENMAX];
    int   ref_natoms = -1;
    int   veclen = 0;
    int   kmax = 0;

    register int   i = 0;
    register int   j = 0;
    register int   k = 0;
    register int   t = 0;

    State save;

    pr( logFile, "\n\t\tCLUSTER ANALYSIS OF CONFORMATIONS\n\t\t_________________________________\n\nNumber of conformations = %d\n", nconf );

    for (j = 0; j < MAX_RUNS; j++) {
        num_in_clu[j] = 0;
        isort[j] = j;
        for (i = 0; i < MAX_RUNS; i++) {
            cluster[j][i] = 0;
        }
    }

    // Read in reference coordinates...
    if (strncmp(FN_rms_ref_crds,"unspecified filename",20) != 0) {
        if ((ref_natoms = getpdbcrds( FN_rms_ref_crds, ref_crds)) == -1) {
            fprintf( logFile, "%s: Problems while reading reference coordinates file \"%s\".\n", programname, FN_rms_ref_crds);
            fprintf( logFile, "Will attempt to use the input PDBQ file coordinates as reference instead.\n");
        } else if (ref_natoms != natom) {
            pr( logFile, "%s: ERROR!  Wrong number of atoms in reference structure.\n", programname);
            pr( logFile, "Input PDBQ structure has %d atoms, but reference structure has %d atoms.\n\n", natom, ref_natoms);
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
        (void)memcpy(crdSave[k], crd, natom*3*sizeof(FloatOrDouble));
    } /*k*/

    flushLog;

    // Sort conformations by energy and perform cluster analysis,
    if (nconf > 1) {
        sort_enrg( econf, isort, nconf );

        ncluster = cluster_analysis( clus_rms_tol, cluster, num_in_clu, isort, 
                    nconf, natom, type, crdSave, crdpdb, 
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
        ref_rms[0] = getrms( crd, ref_crds, B_symmetry_flag, natom, type);
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

            (void)memcpy(crd, crdSave[c], natom*3*sizeof(FloatOrDouble));
     
            if (ntor > 0) {
                eintra = eintcal( nonbondlist, ptr_ad_energy_tables, crd, Nnb, B_calcIntElec, q1q2, B_include_1_4_interactions, scale_1_4, qsp_abs_charge, parameterArray, unbound_internal_FE);
            } else {
                // eintra = torsFreeEnergy;
                eintra = 0.0;
            }
            if (!B_template) {
                // we assume that some atoms might be outside the grid -
                // we do not know to be honest, but this is a safer usage of trilinterp
                 einter = trilinterp( crd, charge, abs_charge, type, natom, map, 
                    info, SOME_ATOMS_OUTSIDE_GRID, ignore_inter, 
                    elec, emap, &elec_total, &emap_total);
            } else {
                 einter = template_trilinterp( crd, charge, abs_charge, type, natom, map, 
                        info, SOME_ATOMS_OUTSIDE_GRID, ignore_inter, 
                        template_energy, template_stddev, 
                        elec /* set */ , emap /* set */, &elec_total, &emap_total );
            }

            print_rem( logFile, i1, num_in_clu[i], c1, ref_rms[c]);
            printEnergies( einter, eintra, torsFreeEnergy, "USER    ", ligand_is_inhibitor, emap_total, elec_total, unbound_internal_FE );
     
            pr( logFile, "USER  \n");
            pr( logFile, "USER    DPF = %s\n", dock_param_fn);
            pr( logFile, "USER    NEWDPF move\t%s\n", smFileName );
            pr( logFile, "USER    NEWDPF about\t%f %f %f\n", sml_center[X],sml_center[Y],sml_center[Z]);
            pr( logFile, "USER    NEWDPF tran0\t%f %f %f\n", hist[c].T.x, hist[c].T.y, hist[c].T.z );
            pr( logFile, "USER    NEWDPF quat0\t%f %f %f %f\n", hist[c].Q.nx, hist[c].Q.ny, hist[c].Q.nz, Deg(hist[c].Q.ang) );
            if (ntor > 0) {
                pr( logFile, "USER    NEWDPF ndihe\t%d\n", hist[c].ntor );
                pr( logFile, "USER    NEWDPF dihe0\t" );
                flushLog;
                for ( t = 0;  t < hist[c].ntor;  t++ ) {
                    torDeg = Deg(hist[c].tor[t]);
                    modtorDeg = ModDeg(torDeg);
                    pr( logFile, "%.2f ", Wrp(modtorDeg) );
                }/*t*/
                pr( logFile, "\n" );
            }/*if*/
            pr( logFile, "USER  \n");
            flushLog;
     
            if (keepresnum > 0) {
                if (outlev > -11) {
                    // Log File PDBQ coordinates [
                    if (!B_template) {
                        pr( logFile, "USER                              x       y       z    vdW   Elec        q     RMS \n" );
                    } else {
                        pr( logFile, "USER                              x       y       z    vdW   Template    q     RMS \n" );
                    }
                    for (j = 0;  j < natom;  j++) {
                        strncpy( rec14, &atomstuff[j][13], (size_t)13);
                        rec14[13]='\0';
                        pr(logFile, FORMAT_PDBQ_ATOM_RESSTR, "", j+1, rec14, crd[j][X], crd[j][Y], crd[j][Z], min(emap[j], MaxValue), min(elec[j], MaxValue), charge[j]);
                        pr(logFile," %6.3f\n", ref_rms[c]); 
                    }
                    //]
                }
            } else {
                if (outlev > -11) {
                    // Log File PDBQ coordinates [
                    if (!B_template) {
                        pr( logFile, "USER                 Rank         x       y       z    vdW   Elec        q     RMS \n");
                    } else {
                        pr( logFile, "USER                 Rank         x       y       z    vdW   Template    q     RMS \n" );
                    }
                    for (j = 0;  j < natom;  j++) {
                        strncpy( rec9, &atomstuff[j][13], (size_t)8);
                        rec9[8]='\0';
                        pr(logFile, FORMAT_PDBQ_ATOM_RESNUM, "", j+1, rec9, i1, crd[j][X], crd[j][Y], crd[j][Z], min(emap[j], MaxValue), min(elec[j], MaxValue), charge[j]);
                        pr(logFile," %6.3f\n", ref_rms[c]); 
                    }/*j*/
                    //]
                }
            }/*if*/
            pr( logFile, "TER\n" );
            pr( logFile, "ENDMDL\n" );
            // End of outputting coordinates of this "MODEL"...
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
        if (!B_template) {
            strcpy( label, "x y z vdW Elec q RMS\0" );
        } else {
            strcpy( label, "x y z vdW Template q RMS\0" );
        }
    } else {
        off[6]=4;
        off[7]=11;
        veclen = 8;
        if (!B_template) {
            strcpy( label, "x y z vdW Elec q Rank RMS\0" );
        } else {
            strcpy( label, "x y z vdW Template q Rank RMS\0" );
        }
    }
    indpf = strindex( dock_param_fn, ".dpf" );
    strncpy( filename, dock_param_fn, (size_t)indpf );
    filename[ indpf ] = '\0';
    strcat( filename, ".dlg.pdb\0" );

    print_avsfld( logFile, veclen, natom, ncluster, off, 12, label, filename );
    //]
}
/* EOF */
