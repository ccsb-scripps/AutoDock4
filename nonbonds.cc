/*

 $Id: nonbonds.cc,v 1.2 2003/02/26 01:22:11 garrett Exp $

*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/* nonbonds.cc */

#include <math.h>

    #include <stdio.h>
    #include <string.h>
    #include <ctype.h>
    #include "nonbonds.h"

#ifdef DEBUG
extern  FILE    *logFile;
#endif /* DEBUG */

void nonbonds( FloatOrDouble crdpdb[MAX_ATOMS][SPACE],
               int nbmatrix_binary[MAX_ATOMS][MAX_ATOMS],
               int natom,
               int atomnumber[MAX_RECORDS], 
               int nrecord,
               char record[MAX_RECORDS][LINE_LEN],
               int piece[MAX_ATOMS],
               int Htype,
               int type[MAX_ATOMS] )

{
    char rec5[5];
    static FloatOrDouble d12max  = 2.101;
    static FloatOrDouble d12min  = 1.000;
    static FloatOrDouble d13max  = 2.650;
    static FloatOrDouble d13min  = 2.100;
    static FloatOrDouble d14max  = 3.970;
    static FloatOrDouble d14min  = 2.651;
    static FloatOrDouble dH12max = 1.437;
    static FloatOrDouble dH12min = 0.900;
    static FloatOrDouble dH13max = 2.200;
    static FloatOrDouble dH13min = 1.700;
    static FloatOrDouble dH14max = 3.310;
    static FloatOrDouble dH14min = 2.201;
    FloatOrDouble rij = 0.;
    FloatOrDouble threshold = 4.0;

    int ii = 0;
    int npiece = 0;
    register int xyz = 0;
    register int i = 0;
    register int j = 0;

    FloatOrDouble d[SPACE];

/*
** 
** Create list of internal non-bond distances to check...
** 
*/

#ifdef DEBUG
    pr( logFile, "\n__A__              " );
    for (j = 0;  j < natom;  j++) {
	/*pr( logFile, " %c %-2d", atm_typ_str[type[j]], j);*/
	pr( logFile, " ? %-2d", j);
    }
    pr( logFile, "\n__A__              " );
    for (j = 0;  j < natom;  j++) {
	pr( logFile, " ____" );
    }
#endif /* DEBUG */

    for (i = 0;  i < natom;  i++) {

#ifdef DEBUG
        pr( logFile, "\n__A__ d[? %-2d][ ] = ",i);
        /* pr( logFile, "\n__A__ d[%c %-2d][ ] = ",atm_typ_str[type[i]],i); */
#endif /* DEBUG */

        for (j = 0;  j < natom;  j++) {
            nbmatrix_binary[i][j] = 1;    /*** MUST KEEP THIS LINE! ***/
            for (xyz = 0;  xyz < SPACE;  xyz++) {
                d[xyz] = crdpdb[i][xyz] - crdpdb[j][xyz] ;
            }/*xyz*/
            rij = sqhypotenuse( d[X], d[Y], d[Z] );

#ifdef DEBUG
            pr( logFile, (rij < 100.)?"|%4.2f":"|%4.1f", sqrt(rij) );
#endif /* DEBUG */

            if (rij < threshold) {
                nbmatrix_binary[i][j] = nbmatrix_binary[j][i] = 0;
            }/*if*/
        }/*j*/

#ifdef DEBUG
            pr( logFile, "\n__A__              " );
#endif /* DEBUG */

        for (j = 0;  j < natom;  j++) {
            for (xyz = 0;  xyz < SPACE;  xyz++) {
                d[xyz] = crdpdb[i][xyz] - crdpdb[j][xyz];
            }/*xyz*/
            rij = hypotenuse( d[X], d[Y], d[Z] );
            if ((type[i] != Htype) && (type[j] != Htype)) {
                if (rij < d12min) {

#ifdef DEBUG
                    pr( logFile, "|****" );
#endif /* DEBUG */

                    /* sep = 0; */
                    nbmatrix_binary[i][j] = nbmatrix_binary[j][i] = 0;
                } else if ((rij >= d12min) && (rij <= d12max)) {

#ifdef DEBUG
                    pr( logFile, "|_12_" );
#endif /* DEBUG */

                    /* sep = 2; */
                    nbmatrix_binary[i][j] = nbmatrix_binary[j][i] = 0;
                } else if ((rij >= d13min) && (rij <= d13max)) {

#ifdef DEBUG
                    pr( logFile, "|_13_" );
#endif /* DEBUG */

                    /* sep = 3; */
                    nbmatrix_binary[i][j] = nbmatrix_binary[j][i] = 0;
                } else if ((rij >= d14min) && (rij <= d14max)) {

#ifdef DEBUG
                    pr( logFile, "|_14_" );
#endif /* DEBUG */

                    /* sep = 4; */
                } else if (rij>d14max) {

#ifdef DEBUG
                    pr( logFile, "|_15+" );
#endif /* DEBUG */

                    /* sep = 5; */
                } else {

#ifdef DEBUG
                    pr( logFile, "|_?__" );
#endif /* DEBUG */

                    /* sep = 9; */
                    nbmatrix_binary[i][j] = nbmatrix_binary[j][i] = 0;
                }
            } else {
                if (rij < dH12min) {

#ifdef DEBUG
                    pr( logFile, "|H***" );
#endif /* DEBUG */

                    /* sep = 0; */
                    nbmatrix_binary[i][j] = nbmatrix_binary[j][i] = 0;
                } else if ((rij >= dH12min) && (rij <= dH12max)) {

#ifdef DEBUG
                    pr( logFile, "|H12_" );
#endif /* DEBUG */

                    /* sep = 2; */
                    nbmatrix_binary[i][j] = nbmatrix_binary[j][i] = 0;
                } else if ((rij >= dH13min) && (rij <= dH13max)) {

#ifdef DEBUG
                    pr( logFile, "|H13_" );
#endif /* DEBUG */

                    /* sep = 3; */
                    nbmatrix_binary[i][j] = nbmatrix_binary[j][i] = 0;
                } else if ((rij >= dH14min) && (rij <= dH14max)) {
                    /* sep = 4; */

#ifdef DEBUG
                    pr( logFile, "|H14_" );
#endif /* DEBUG */

                } else if (rij > dH14max) {
                    /* sep = 5; */

#ifdef DEBUG
                    pr( logFile, "|H15+" );
#endif /* DEBUG */

                } else {
                    /* sep = 9; */

#ifdef DEBUG
                    pr( logFile, "| ?  " );
#endif /* DEBUG */

                    nbmatrix_binary[i][j] = nbmatrix_binary[j][i] = 0;
                }/*endif*/
            }/*endif*/
        }/*j*/
    }/*i*/

#ifdef DEBUG
    pr( logFile, "\n\n" );
#endif /* DEBUG */

    for (j = 0;  j < nrecord;  j++) {
        i = atomnumber[j];
        for (ii = 0; ii < 4; ii++)  {
	    rec5[ii] = (char)tolower( (int)record[j][ii] ); 
	}
        if ( equal(rec5,"atom", 4) || equal(rec5,"heta", 4)) {
            piece[i] = npiece;
        } else {
            ++npiece;
        }/*if*/
    }/*j*/
}
/* EOF */
