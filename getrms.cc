/*

 $Id: getrms.cc,v 1.4 2006/07/15 02:38:44 garrett Exp $

*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/* getrms.cc */

#include <math.h>

    #include "getrms.h"



Real getrms ( Real Crd[MAX_ATOMS][SPACE], 
               Real CrdRef[MAX_ATOMS][SPACE], 
               Boole B_symmetry_flag, 
               int natom, 
               int type[MAX_ATOMS] )

{
    double sqrSum, sqrMin, dc[SPACE];
    register int i, j, xyz;

    sqrSum = 0.;

    if (B_symmetry_flag) {
        for (i = 0;  i < natom;  i++) {
            sqrMin = BIG;
            for (j = 0;  j < natom;  j++) {                
                if (type[i] == type[j]) {
                    for (xyz = 0;  xyz < SPACE;  xyz++) {
                        dc[xyz]= Crd[i][xyz] - CrdRef[j][xyz];
                    } /* xyz */
                    sqrMin = min( sqhypotenuse(dc[X], dc[Y], dc[Z]), sqrMin );
                }
            } /*  next j  */
            sqrSum += sqrMin;
        } /*  next i  */
    } else {
        for (i = 0;  i < natom;  i++) {
            for (xyz = 0;  xyz < SPACE;  xyz++) {
                dc[xyz]= Crd[i][xyz] - CrdRef[i][xyz];
            } /* xyz */
            sqrSum += sqhypotenuse( dc[X], dc[Y], dc[Z] );
        } /*  next i  */
    }

    return ( sqrt( sqrSum / (double)natom )  );
}
/* EOF */
