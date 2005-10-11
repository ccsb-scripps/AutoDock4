/*

 $Id: input_state.cc,v 1.2.8.1 2005/10/11 00:05:46 alther Exp $

*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/* input_state.cc */

#include <stdio.h>
#include <string.h>
#include "input_state.h"
#include "qmultiply.h"

#define LINELEN 132

int input_state( State *S,
		 FILE  *fp,
		 char  line[LINE_LEN],
		 int   ntor,
		 int   *p_istep,
		 FloatOrDouble *p_energy,
		 FloatOrDouble *p_eint,
		 char  *p_lastmove )
{
    int i, istep, status;
    FloatOrDouble energy, eint;
    char lastmove;
    char myline[LINELEN];

#ifdef DEBUG
    fprintf(stderr, "line=|%s|\n", line);
#endif /* DEBUG */

    #ifdef USE_DOUBLE
        status = sscanf(line, "%*s %d %1s %lf %lf %lf %lf %lf %lf %lf %lf %lf", &istep, &lastmove, &energy, &eint,  &(S->T.x), &(S->T.y), &(S->T.z),  &(S->Q.nx), &(S->Q.ny), &(S->Q.nz),  &(S->Q.ang) );
    #else
        status = sscanf(line, "%*s %d %1s %f %f %lf %lf %lf %lf %lf %lf %lf", &istep, &lastmove, &energy, &eint,  &(S->T.x), &(S->T.y), &(S->T.z),  &(S->Q.nx), &(S->Q.ny), &(S->Q.nz),  &(S->Q.ang) );
    #endif

    if (status != 0) {
	S->Q.ang = Rad( S->Q.ang );
	mkUnitQuat( &(S->Q) );

        *p_istep = istep;
	*p_energy = energy;
	*p_eint = eint;
	*p_lastmove = lastmove;

        for (i=0; i<ntor; i++) {
	    (void) fgets(myline, LINELEN, fp);
            sscanf(myline, "%lf", &(S->tor[i]) ); /* input torsions are in degrees */
            S->tor[i] = Rad( S->tor[i] );    /* now in radians */
        }
    }
    return( status );
}
/* EOF */
