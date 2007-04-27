/*

 $Id: output_state.cc,v 1.7 2007/04/27 06:01:50 garrett Exp $

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

#include <stdio.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "structs.h"
#include "output_state.h"

/* LOCK_SH 1       shared lock */
/* LOCK_EX 2       exclusive lock */
/* LOCK_NB 4       don't block when locking */
/* LOCK_UN 8       unlock */
#define PERMS 0666        /* hexadecimal permissions for watch-file */

/*----------------------------------------------------------------------------*/
void output_state( FILE *fp,
		   State S,
                   int ntor,
                   int istep,
                   Real energy,
                   Real eint,
                   char lastmove,
                   Boole B_watch,
                   char FN_watch[MAX_CHARS],
                   char atomstuff[MAX_ATOMS][MAX_CHARS],
                   int natom,
                   Real crd[MAX_ATOMS][SPACE])
/*----------------------------------------------------------------------------*/
{
    int i;
	/*int lockf_status;*/
	
//#ifndef __ppc__
#if defined( __ppc__ ) || defined( __CYGWIN__ )
    // F_LOCK is not supported
#else
    int FD_watch;
    FILE *FP_watch;
#endif

    fprintf(fp, "state %d %c %f %f  %lf %lf %lf  %lf %lf %lf %lf\n",
        istep, lastmove, energy, eint, S.T.x, S.T.y, S.T.z,
        S.Q.nx, S.Q.ny, S.Q.nz, RadiansToDegrees( S.Q.ang ) );

    for (i=0; i<ntor; i++) {
        fprintf(fp, "%f\n", RadiansToDegrees( S.tor[i]) );
    }

/* >>>> NOW USES lockf !!!! <<<< */

// #ifndef __ppc__
#if defined( __ppc__ ) || defined( __CYGWIN__ )
	// F_LOCK is not supported
#else
    if (B_watch) {
        if ((FD_watch = creat( FN_watch, PERMS )) != -1) {;
            /* creates new file, or re-write old one */

            if ((FP_watch = fdopen( FD_watch, "w")) != NULL ) {
                /*lockf_status = lockf( FD_watch, F_LOCK, 0 ); */
                (void) lockf( FD_watch, F_LOCK, 0 ); 

                for (i = 0;  i < natom;  i++) {
		    fprintf( FP_watch, "%30s%8.3f%8.3f%8.3f\n",
		    atomstuff[i], crd[i][X], crd[i][Y], crd[i][Z]);
                }
                fclose( FP_watch ); /*lockf_status=lockf(FD_watch,F_ULOCK,0);*/
            }
        }
    }
#endif

    return;
}
/* EOF */
