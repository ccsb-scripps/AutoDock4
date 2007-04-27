/*

 $Id: bestpdb.cc,v 1.4 2007/04/27 06:01:47 garrett Exp $

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

/* bestpdb.cc */

#include <stdio.h>
#include <string.h>
#include "constants.h"
#include "print_rem.h"
#include "strindex.h"
#include "print_avsfld.h"
#include "bestpdb.h"


extern FILE *logFile;
extern int keepresnum;
extern char dock_param_fn[];

void bestpdb( int ncluster, 
	      int num_in_clu[MAX_RUNS],
	      int cluster[MAX_RUNS][MAX_RUNS],
	      Real econf[MAX_RUNS],
	      Real crd[MAX_RUNS][MAX_ATOMS][SPACE],
	      char atomstuff[MAX_ATOMS][MAX_CHARS],
	      int natom,
	      Boole B_write_all_clusmem,
	      Real ref_rms[MAX_RUNS])

{
    register int  i=0,
	          j=0,
	          k=0,
		  confnum=0;

    int           c = 0,
		  kmax = 0,
		  /* imol = 0, */
		  indpf = 0,
		  off[7],
		  nframes = 0,
	          stride = 0,
		  c1 = 1,
		  i1 = 1;
		  

    char          rec13[14],
		  filnm[MAX_CHARS],
		  label[MAX_CHARS],
		  rec8[9];

    pr( logFile, "\n\tLOWEST ENERGY DOCKED CONFORMATION from EACH CLUSTER");
    pr( logFile, "\n\t___________________________________________________\n\n\n" );

    if (keepresnum > 0 ) {
	pr( logFile, "\nKeeping original residue number (specified in the input PDBQ file) for outputting.\n\n");
    } else {
	pr( logFile, "\nResidue number will be the conformation's rank.\n\n");
    }

    for (i = 0;  i < ncluster;  i++) {
        i1 = i + 1;

        if (B_write_all_clusmem) {
            kmax = num_in_clu[i];
        } else {
            kmax = 1;	/* write lowest-energy only */
        }

        for ( k = 0; k < kmax; k++ ) {
            c = cluster[i][k];
            c1 = c + 1;

            fprintf( logFile, "USER    DPF = %s\n", dock_param_fn);
            fprintf( logFile, "USER    Conformation Number = %d\n", ++confnum);

            print_rem(logFile, i1, num_in_clu[i], c1, ref_rms[c]);

            if (keepresnum > 0) {
                fprintf( logFile, "USER                              x       y       z   Rank Run  Energy    RMS\n");
                for (j = 0;  j < natom;  j++) {
                    strncpy( rec13, &atomstuff[j][13], (size_t)13);
                    rec13[13]='\0';
                    fprintf( logFile, FORMAT_PDBQ_ATOM_RANKRUN_STR, j+1, rec13, crd[c][j][X], crd[c][j][Y], crd[c][j][Z], i1, c1, econf[c], ref_rms[c] );
                } /* j */
            } else {
                fprintf( logFile, "USER                   Rank       x       y       z    Run   Energy    RMS\n");
                for (j = 0;  j < natom;  j++) {
                    strncpy( rec8, &atomstuff[j][13], (size_t)8);
                    rec8[8]='\0';
                    fprintf( logFile, FORMAT_PDBQ_ATOM_RUN_NUM, j+1, rec8, i1, crd[c][j][X], crd[c][j][Y], crd[c][j][Z], c1, econf[c], ref_rms[c] );
                } /* j */
            }
            fprintf( logFile, "TER\n" );
            fprintf( logFile, "ENDMDL\n" );
            fflush( logFile );

            nframes++;
        } /* for k */

    } /* for i */

    fprintf( logFile, "\n" );

    strcpy(label, "x y z Rank Run Energy RMS\0" );
    if (keepresnum > 0) {
        off[0]=5; off[1]=6; off[2]=7; off[3]=8; off[4]=9; off[5]=10; off[6]=11;
        stride=12;
    } else {
        off[0]=5; off[1]=6; off[2]=7; off[3]=4; off[4]=8; off[5]=9; off[6]=10;
        stride=11;
    } /* if */
     
    indpf = strindex( dock_param_fn, ".dpf" );
    strncpy( filnm, dock_param_fn, (size_t)indpf );
    filnm[ indpf ] = '\0';
    strcat( filnm, ".dlg.pdb\0" );
     
    print_avsfld( logFile, 7, natom, nframes, off, stride, label, filnm );
}
/* EOF */
