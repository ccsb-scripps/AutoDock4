/* analysis.cc */

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
               float charge[MAX_ATOMS], 
               Boole B_calcIntElec,
               float q1q2[MAX_NONBONDS],
               float clus_rms_tol, 
               float crdpdb[MAX_ATOMS][SPACE], 
               float e_internal[NEINT][ATOM_MAPS][ATOM_MAPS], 
               float inv_spacing, 
               float map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS], 
               float econf[MAX_RUNS], 
               int   irunmax, 
               float xlo, 
               float ylo, 
               float zlo, 
               int   natom, 
               int   nonbondlist[MAX_NONBONDS][2], 
               int   nconf, 
               int   ntor, 
               State hist[MAX_RUNS], 
               char  smFileName[MAX_CHARS], 
               float sml_center[SPACE],
               Boole B_symmetry_flag, 
               int   tlist[MAX_TORS][MAX_ATOMS], 
               int   type[MAX_ATOMS], 
               float vt[MAX_TORS][SPACE],
               char  FN_rms_ref_crds[MAX_CHARS],
               float torsFreeEnergy,
               Boole B_write_all_clusmem,
               int ligand_is_inhibitor,
               Boole B_template,
               float template_energy[MAX_ATOMS],
               float template_stddev[MAX_ATOMS])

{
    /* register int   imol = 0; */
    char  filename[MAX_CHARS];
    char  label[MAX_CHARS];
    char  rec14[14];
    char  rec9[9];

    float clu_rms[MAX_RUNS][MAX_RUNS];
    float crdSave[MAX_RUNS][MAX_ATOMS][SPACE];
    float crd[MAX_ATOMS][SPACE];
    float einter = 0.;
    float eintra = 0.;
    float elec[MAX_ATOMS];
    float emap[MAX_ATOMS];
    float ref_crds[MAX_ATOMS][SPACE];
    float ref_rms[MAX_RUNS];
    float torDeg = 0.;
    float modtorDeg = 0.;
    float MaxValue = 99.99;

    int   c = 0;
    int   c1 = 0;
    int   cluster[MAX_RUNS][MAX_RUNS];
    int   i1=1;
    int   indpf = 0;
    int   isort[MAX_RUNS];
    int   ncluster = 1;
    int   num_in_clu[MAX_RUNS];
    int   off[VECLENMAX];
    int   ref_natoms = -1;
    int   veclen = 0;
    int   kmax = 0;

    // register int   XYZ = 0;
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

    if (strncmp(FN_rms_ref_crds,"unspecified filename",20) != 0) {
        /*
        Read in reference coordinates...
        */
        if ((ref_natoms = getpdbcrds( FN_rms_ref_crds, ref_crds)) == -1) {
            fprintf( logFile, "%s: Problems while reading reference coordinates file \"%s\".\n", programname, FN_rms_ref_crds);
            fprintf( logFile, "Will attempt to use the input PDBQ file coordinates as reference instead.\n");
        } else if (ref_natoms != natom) {
            pr( logFile, "%s: ERROR!  Wrong number of atoms in reference structure.\n", programname);
            pr( logFile, "Input PDBQ structure has %d atoms, but reference structure has %d atoms.\n\n", natom, ref_natoms);
            ref_natoms = -1;
        }
    }

    /*
    Generate coordinates for each final transformation,
    */
    for ( k=0; k<nconf; k++ ) {

        /* fprintf( logFile, "\n\nState hist[%d].\n", k); */
        printState( logFile, hist[k], 2 );

        /* fprintf( logFile, "\nCopying state %d.\n", k); */
        copyState( &save, hist[k] );

        /* fprintf( logFile, "Converting state %d to coordinates.\n", k); */
        cnv_state_to_coords( save, vt, tlist, ntor, crdpdb, crd, natom);

        /* fprintf( logFile, "Saving coordinates of state %d.\n", k); */

        /* Save coordinates in crdSave array...  */
        (void)memcpy(crdSave[k], crd, natom*3*sizeof(float));
        // for (j = 0;  j < natom;  j++) {
            // for (XYZ = 0;  XYZ < SPACE;  XYZ++) {
                // crdSave[k][j][XYZ] = crd[j][XYZ]; // } /* XYZ */ // } /*j*/
    } /*k*/

    flushLog;

    if (nconf > 1) {
        sort_enrg( econf, isort, nconf );

        ncluster = cluster_analysis( clus_rms_tol, cluster, 
                    num_in_clu, isort, 
                    nconf, natom, type, crdSave, crdpdb, 
                    sml_center, clu_rms, B_symmetry_flag,
                    ref_crds, ref_natoms, ref_rms);

        pr( logFile, "\nOutputting structurally similar clusters, ranked in order of increasing energy.\n" );
        flushLog;

        prClusterHist( ncluster, irunmax, clus_rms_tol,num_in_clu, 
                       cluster, econf, clu_rms, ref_rms);

        pr( logFile, "\n\tLOWEST ENERGY DOCKED CONFORMATION from EACH CLUSTER");
        pr( logFile, "\n\t___________________________________________________\n\n\n" );

        if (keepresnum > 0 ) {
            pr( logFile, "\nKeeping original residue number (specified in the input PDBQ file) for outputting.\n\n");
        } else {
            pr( logFile, "\nResidue number will be set to the conformation's cluster rank.\n\n");
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

    for (i = 0;  i < ncluster;  i++) { // each cluster, i
        i1 = i + 1;

        // c = cluster[i][0];
        if (B_write_all_clusmem) {
            kmax = num_in_clu[i];
        } else {
            kmax = 1;	/* write lowest-energy only */
        }

        for (k = 0;  k < kmax;  k++) { // each member of this cluster, k
            c = cluster[i][k];
            c1 = c + 1;

            (void)memcpy(crd, crdSave[c], natom*3*sizeof(float));
            // for (j = 0;  j < natom;  j++) {
                // for (XYZ = 0;  XYZ < SPACE;  XYZ++) {
                    // crd[j][XYZ] = crdSave[c][j][XYZ]; // } // }/*j*/
     
            if (ntor > 0) {
                // eintra = eintcal( nonbondlist,e_internal,crd,type,Nnb,B_calcIntElec,q1q2 ) + torsFreeEnergy;
                eintra = eintcal( nonbondlist,e_internal,crd,type,Nnb,B_calcIntElec,q1q2 );
            } else {
                // eintra = torsFreeEnergy;
                eintra = 0.0;
            }
            if (!B_template) {
                 einter = trilinterp( crd, charge, type, natom, map, inv_spacing, elec, emap, xlo, ylo, zlo );
            } else {
                 einter = byatom_template_trilinterp( crd, charge, type, natom, map, inv_spacing, elec, emap, xlo, ylo, zlo,
                                                      template_energy, template_stddev);
            }
     

            print_rem( logFile, i1, num_in_clu[i], c1, ref_rms[c]);
            printEnergies( einter, eintra, torsFreeEnergy, "USER    ", ligand_is_inhibitor );
     
            pr( logFile, "USER\n");
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
            pr( logFile, "USER\n");
            flushLog;
            off[0]=5; off[1]=6; off[2]=7; off[3]=8; off[4]=9; off[5]=10;
     
            if (keepresnum > 0) {
                if (!B_template) {
                    pr( logFile, "USER                              x       y       z    vdW   Elec        q     RMS \n" );
                } else {
                    pr( logFile, "USER                              x       y       z    vdW   Template    q     RMS \n" );
                }
                /*
                 * 123456
                 * +99.99
                 */
                for (j = 0;  j < natom;  j++) {
                    strncpy( rec14, &atomstuff[j][13], (size_t)13);
                    rec14[13]='\0';
                    pr(logFile, FORMAT_PDBQ_ATOM_RESSTR, "", j+1, rec14, crd[j][X], crd[j][Y], crd[j][Z], min(emap[j], MaxValue), min(elec[j], MaxValue), charge[j]);
                    pr(logFile," %6.3f\n", ref_rms[c]); 
                    /* 
                    pr( logFile, "ATOM  %5d  %13s    %8.3f%8.3f%8.3f%+6.2f%+6.2f    %+6.3f %6.3f\n", j+1, rec14, crd[j][X], crd[j][Y], crd[j][Z], min(emap[j], MaxValue), min(elec[j], MaxValue), charge[j], ref_rms[c]);
                    */
                }
                if (!B_template) {
                    strcpy( label, "x y z vdW Elec q RMS\0" );
                } else {
                    strcpy( label, "x y z vdW Template q RMS\0" );
                }
                veclen = 7;
                off[6]=11;
            } else {
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
                    /*
                    pr( logFile, "ATOM  %5d  %8s%5d    %8.3f%8.3f%8.3f%+6.2f%+6.2f    %+6.3f %6.3f\n",
                    j+1, rec9, i1, crd[j][X], crd[j][Y], crd[j][Z], min(emap[j], MaxValue), min(elec[j], MaxValue), charge[j], ref_rms[c]);
                    */
                }/*j*/
                if (!B_template) {
                    strcpy( label, "x y z vdW Elec q Rank RMS\0" );
                } else {
                    strcpy( label, "x y z vdW Template q Rank RMS\0" );
                }
                veclen = 8;
                off[6]=4;
                off[7]=11;
            }/*if*/
            pr( logFile, "TER\n" );
            pr( logFile, "ENDMDL\n" );
            flushLog;
        } /*k*/
    } /*i   (Next cluster.) */
    pr( logFile, "\n\n" );

    indpf = strindex( dock_param_fn, ".dpf" );
    strncpy( filename, dock_param_fn, (size_t)indpf );
    filename[ indpf ] = '\0';
    strcat( filename, ".dlg.pdb\0" );

    print_avsfld( logFile, veclen, natom, ncluster, off, 12, label, filename );

}
/* EOF */
