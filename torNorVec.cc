/* torNorVec.cc */

#include <math.h>

    #include <stdio.h>
    #include <stdlib.h>
    #include "torNorVec.h"
    #ifdef DEBUG
	#include <ctype.h>
    #endif /* DEBUG */


extern FILE *logFile;


void torNorVec( float crdpdb[MAX_ATOMS][SPACE],
		int ntor,
		int tlist[MAX_TORS][MAX_ATOMS],
		float vt[MAX_TORS][SPACE] )
{

    register int xyz = 0;
    register int j = 0;

    float magVec = 0.;
    float imagVec = 0.;
    float v[SPACE];

    char error_message[LINE_LEN];

    /*_____________________________________________________________
      | Calculate normal vectors of torsion bonds,                 |
      |____________________________________________________________|
      |      Note: torsions must be rotated from leaves to root    |
      |       or this pre-calculation of normal vectors will fail. |
      |____________________________________________________________| */

    pr( logFile, "\nCalculating normalized torsion vectors.\n\n" );

    for (j=0; j<ntor; j++) {
	for (xyz = 0;  xyz < SPACE;  xyz++) {
	    v[xyz] = crdpdb[ tlist[j][ATM2] ][xyz] - crdpdb[ tlist[j][ATM1] ][xyz];

#ifdef DEBUG
		pr( logFile, "\n__norm__ Torsion %d, crdpdb[ %d,ATM2 ][%c] = %.3f, crdpdb[ %d,ATM1 ][%c] = %.3f\n",
		    j, tlist[j][ATM2], toascii(88+xyz), crdpdb[tlist[j][ATM2]][xyz], 
		    tlist[j][ATM1], toascii(88+xyz), crdpdb[tlist[j][ATM1]][xyz] );
		flushLog; 
#endif /* DEBUG */

	} /* xyz */

	magVec = hypotenuse( v[X], v[Y], v[Z] );  /* Magnitude of vector v[xyz] */

	if (magVec == 0.) {
	    prStr( error_message, "Torsion %d, normal vector, magVec, is 0; imminent division by zero caught.", j );
	    stop( error_message );
	    exit( -1 );
	}

	imagVec = 1. / magVec;

	for (xyz = 0;  xyz < SPACE;  xyz++) {
	    vt[j][xyz] = v[xyz] * imagVec;  /* Normalize. */
	} /* xyz */

    } /* j */
    flushLog;
}
/* EOF */
