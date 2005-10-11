/*

 $Id: output_state.cc,v 1.3.8.1 2005/10/11 00:12:16 alther Exp $

*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/* output_state.cc */

#include <stdio.h>
#include <fcntl.h>

#ifdef _WIN32
   #include <io.h>
   #include <sys/locking.h>
#else
   #include <unistd.h>
#endif

#include <sys/stat.h>
#include <sys/types.h>
#include "structs.h"
#include "output_state.h"

/* LOCK_SH 1       shared lock */
/* LOCK_EX 2       exclusive lock */
/* LOCK_NB 4       don't block when locking */
/* LOCK_UN 8       unlock */

/* hexadecimal permissions for watch-file */
#ifdef _WIN32
   #define PERMS (_S_IREAD | _S_IWRITE)
#else
   #define PERMS 0666
#endif

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
	
#ifndef __ppc__
	int FD_watch;
    FILE *FP_watch;
#endif

    fprintf(fp, "state %d %c %f %f  %lf %lf %lf  %lf %lf %lf %lf\n",
        istep, lastmove, energy, eint, S.T.x, S.T.y, S.T.z,
        S.Q.nx, S.Q.ny, S.Q.nz, Deg( S.Q.ang ) );

    for (i=0; i<ntor; i++) {
        fprintf(fp, "%f\n", Deg( S.tor[i]) );
    }

/* >>>> NOW USES lockf !!!! <<<< */

#ifndef __ppc__
    if (B_watch) {

       /* creates new file, or re-write old one */
#ifdef _WIN32
       if ((FD_watch = _creat( FN_watch, PERMS )) != -1) {
#else
       if ((FD_watch = creat( FN_watch, PERMS )) != -1) {
#endif

            if ((FP_watch = fdopen( FD_watch, "w")) != NULL ) {
                /*lockf_status = lockf( FD_watch, F_LOCK, 0 ); */
#ifdef _WIN32
                (void) _locking(FD_watch, _LK_LOCK, 0);
#else
                (void) lockf( FD_watch, F_LOCK, 0 );
#endif

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
