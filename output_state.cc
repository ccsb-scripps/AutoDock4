/*

 $Id: output_state.cc,v 1.2 2003/02/26 01:24:26 garrett Exp $

*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/* output_state.cc */

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
                   FloatOrDouble energy,
                   FloatOrDouble eint,
                   char lastmove,
                   Boole B_watch,
                   char FN_watch[MAX_CHARS],
                   char atomstuff[MAX_ATOMS][MAX_CHARS],
                   int natom,
                   FloatOrDouble crd[MAX_ATOMS][SPACE])
/*----------------------------------------------------------------------------*/
{
    int i;
	/*int lockf_status;*/
	int FD_watch;
    FILE *FP_watch;

    fprintf(fp, "state %d %c %f %f  %lf %lf %lf  %lf %lf %lf %lf\n",
        istep, lastmove, energy, eint, S.T.x, S.T.y, S.T.z,
        S.Q.nx, S.Q.ny, S.Q.nz, Deg( S.Q.ang ) );

    for (i=0; i<ntor; i++) {
        fprintf(fp, "%f\n", Deg( S.tor[i]) );
    }

/* >>>> NOW USES lockf !!!! <<<< */

#ifndef __ppc__
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
