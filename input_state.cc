/* input_state.cc */

#include <math.h>


    #include <stdio.h>
    #include <string.h>
    #include "input_state.h"

#define LINELEN 132

int input_state( State *S,
		 FILE  *fp,
		 char  line[LINE_LEN],
		 int   ntor,
		 int   *p_istep,
		 float *p_energy,
		 float *p_eint,
		 char  *p_lastmove )
{
    int i, istep, status;
    float energy, eint;
    char lastmove;
    char myline[LINELEN];

#ifdef DEBUG
    fprintf(stderr, "line=|%s|\n", line);
#endif /* DEBUG */

    status = sscanf(line, "%*s %d %1s %f %f %lf %lf %lf %lf %lf %lf %lf", &istep, &lastmove, &energy, &eint,  &(S->T.x), &(S->T.y), &(S->T.z),  &(S->Q.nx), &(S->Q.ny), &(S->Q.nz),  &(S->Q.ang) );

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
