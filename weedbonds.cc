/* weedbonds.cc */

    #include <stdio.h>
    #include <stdlib.h>
    #include "weedbonds.h"


extern FILE *logFile;


void weedbonds( int natom,
		char pdbaname[MAX_ATOMS][5],
		int piece[MAX_ATOMS],
		int ntor,
		int tlist[MAX_TORS][MAX_ATOMS],
		int *Addr_Nnb,
		int Nnbonds[MAX_ATOMS],
		int nbmatrix_binary[MAX_ATOMS][MAX_ATOMS],
		int nonbondlist[MAX_NONBONDS][2] )

{
    static int OUTNUMATM = 10;

    int a11=0;
    int a12=0;
    int a21=0;
    int a22=0;
    int p11 = 0;
    int p12 = 0;
    int p21 = 0;
    int p22 = 0;
    int i_atmnum = 0;
    int j_atmnum = 0;
    int n_a = 0;
    int p = 0;
    int repflag = FALSE;
    int Nnb = 0;

    register int i = 0;
    register int j = 0;
    register int k = 0;

    char error_message[LINE_LEN];

/*___________________________________________________________________________
|    ENDBRANCH---TORS---BRANCH---R O O T---BRANCH---ENDBRANCH                |
|                                  /              \                          |
|                                BRANCH            BRANCH--TORS---ENDBRANCH  |
|                                /                  \                        |
|                               ENDBRANCH            ENDBRANCH               |
|____________________________________________________________________________|
|  Eliminate all rigidly bonded atoms:                                       |
|                                     * those atoms which are at the ROOT;   |
|                                     * atoms between TORS and BRANCH;       |
|                                     * atoms between BRANCH and ENDBRANCH.  |
|  This is necessary for internal energy calculations.                       |
|____________________________________________________________________________|
| Weed out bonds in rigid pieces,                                            |
|____________________________________________________________________________|
*/
    for (j = 0;  j < natom;  j++) {
	for (i = 0;  i < natom;  i++) {
	    if (piece[i] == piece[j]) {
		nbmatrix_binary[j][i] = nbmatrix_binary[i][j] = 0;
	    } /* endif */
	} /*  i  */
    } /*j*/
    /* 
    \   Weed out bonds across torsions,
     \______________________________________________________________
    */
    for (i=0; i<ntor; i++) {
	nbmatrix_binary[ tlist[i][ATM2] ][ tlist[i][ATM1] ] = 0;

#ifdef DEBUG
/**/	    pr( logFile, "__2__ i=%-2d: nbmatrix_binary[tlist[%-2d][ATM2]=%-2d][tlist[%-2d][ATM1]=%-2d]=%d\n",i,i,tlist[i][ATM2],i,tlist[i][ATM1],nbmatrix_binary[tlist[i][ATM2]][tlist[i][ATM1]]);
#endif /* DEBUG */

    } /* i */
    /*
    \  Weed out bonds from atoms directly connected to rigid pieces,
     \______________________________________________________________
    */
    for (i=0; i<ntor; i++) {

	a11 = tlist[i][ATM1];
	a21 = tlist[i][ATM2];
	p11 = piece[a11];
	p21 = piece[a21];

	for (j=0; j<ntor; j++) {

	    a12 = tlist[j][ATM1];
	    a22 = tlist[j][ATM2];
	    p12 = piece[a12];
	    p22 = piece[a22];

	    if (p11 == p12)  { nbmatrix_binary[ a22 ][ a21 ] = nbmatrix_binary[ a21 ][ a22 ] = 0; }
	    if (p11 == p22)  { nbmatrix_binary[ a12 ][ a21 ] = nbmatrix_binary[ a21 ][ a12 ] = 0; }
	    if (p21 == p12)  { nbmatrix_binary[ a22 ][ a11 ] = nbmatrix_binary[ a11 ][ a22 ] = 0; }
	    if (p21 == p22)  { nbmatrix_binary[ a12 ][ a11 ] = nbmatrix_binary[ a11 ][ a12 ] = 0; }

	} /* j */

	for (k = 0;  k < natom;  k++) {

	    p = piece[k];
	    if (p11 == p)  { nbmatrix_binary[ k ][ a21 ] = nbmatrix_binary[ a21 ][ k ] = 0; }
	    if (p21 == p)  { nbmatrix_binary[ k ][ a11 ] = nbmatrix_binary[ a11 ][ k ] = 0; }

	} /* k */
    } /* i */

    for ( i = 0;  i < (natom-1);  i++ ) {

#ifdef DEBUG
/**/	pr( logFile,"> i = %d (CHECKING SYMMETRY OF NBMATRIX)\n\n", i);
#endif /* DEBUG */

	for ( j = i+1;  j < natom;  j++ ) {

#ifdef DEBUG
/**/	    pr( logFile,"< i,j = %2d,%2d  nbmatrix_binary[][] = %d/%d\n", i, j, nbmatrix_binary[i][j], nbmatrix_binary[j][i] );
#endif /* DEBUG */

	    if ((nbmatrix_binary[i][j] == 1) && (nbmatrix_binary[j][i] == 1)) {
		nonbondlist[Nnb][ATM1] = i;
		nonbondlist[Nnb][ATM2] = j;

#ifdef DEBUG
/**/		pr( logFile,"< nonbondlist[%2d][0,1] = %2d,%2d\n", Nnb, nonbondlist[Nnb][ATM1], nonbondlist[Nnb][ATM2] );
#endif /* DEBUG */

		++Nnb;
	    } else if ( (nbmatrix_binary[i][j] == 1) && (nbmatrix_binary[j][i] == 0)
		     || (nbmatrix_binary[i][j] == 0) && (nbmatrix_binary[j][i] == 1) ) {
		pr( logFile, "ERROR: ASSYMMETRY detected in Non-Bond Matrix at %d,%d\n", i,j);
	    }
	} /* j */
    } /* i */

    if (Nnb > MAX_NONBONDS) {
	prStr( error_message, "too many non-bonded interactions (%d) in small molecule\n\t(increase MAX_NONBONDS from %d).", Nnb, MAX_NONBONDS );
	stop( error_message );
	exit( -1 );
    } else {
	*Addr_Nnb = Nnb;
    }

    if (ntor > 0) {
	pr( logFile, "\n\nMatrix of Non-Bonded Interactions:\n" );
	pr( logFile, "__________________________________\n\n" );
	pr( logFile, "Key:\nX = non-bonded interaction\n" );
	pr( logFile, "_ = 1,2 or 1,3 interaction\n\n" );
	pr( logFile, "\nAtom: ID: " ); for (j = 0;  j < natom;  j++)   fprintf( logFile, "%2d", (1+j) );
	pr( logFile, "\n_____ ___ " ); for (j = 0;  j < natom;  j++)   fprintf( logFile, "__" );
	for (j = 0;  j < natom;  j++) {
	    pr( logFile, "\n%4s  %2d  ", pdbaname[j], 1+j );
	    for (i = 0;  i < natom;  i++) {
		pr( logFile, "|%c", (nbmatrix_binary[j][i])?'X':'_' );
	    } /* i */
	} /* j */
	pr( logFile, "\n\n" );
	flushLog;
    } /*  endif  */

#ifdef DEBUG
/**/	for (i = 0;  i < Nnb;  i++) {
/**/	    pr( logFile,"> nonbondlist[%2d][0,1] = %2d,%2d\n", i,nonbondlist[i][ATM1],nonbondlist[i][ATM2] );
/**/	} /*  i  */
#endif /* DEBUG */

    for (i = 0;  i < Nnb;  i++) {

	i_atmnum = nonbondlist[i][ATM1];
	j_atmnum = nonbondlist[i][ATM2];

#ifdef DEBUG
/**/	    pr( logFile, "* Assigning nbmatrix_binary[%d][%d] to %d *\n", Nnbonds[i_atmnum], i_atmnum, j_atmnum );
#endif /* DEBUG */

	nbmatrix_binary[Nnbonds[i_atmnum]][i_atmnum] = j_atmnum;
	
	/* NOTE: re-utilizes the nbmatrix_binary array;
	\        nbmatrix_binary is `corrupted' after this...
	 \ 
	  \ Normally, nbmatrix_binary contains 0s and 1s, but this
	   \ assigns atom-numbers, which can be >1.
	  */

#ifdef DEBUG
/**/	    pr( logFile, ">>--> i = %d, j_atmnum = %d, nbmatrix_binary[ %d ][ %d ] = %d\n",
/**/	    i, j_atmnum, Nnbonds[i_atmnum],i_atmnum, nbmatrix_binary[Nnbonds[i_atmnum]][i_atmnum] );
#endif /* DEBUG */

	++Nnbonds[i_atmnum];
    } /*  i  */

#ifdef DEBUG
/**/    if (ntor > 0) {
/**/	    pr( logFile, "\nCORRUPTED Matrix of Non-Bonded Interactions:\n" );
/**/	    pr( logFile, "____________________________________________\n\n" );
/**/	    pr( logFile, "\nAtom: ID: " ); for (j = 0;  j < natom;  j++)   fprintf( logFile, "%2d", j );
/**/	    pr( logFile, "\n_____ ___ " ); for (j = 0;  j < natom;  j++)   fprintf( logFile, "__" );
/**/	    for (i = 0;  i < natom;  i++) {
/**/	        pr( logFile, "\n%4s  %2d  ", pdbaname[i], i );
/**/	        for (j = 0;  j < natom;  j++) {
/**/	    	    pr( logFile, "%2d", nbmatrix_binary[j][i] );
/**/	        } /* j */
/**/	    } /* i */
/**/	    pr( logFile, "\n\n" );
/**/	    flushLog;
/**/    } /* endif */
#endif /* DEBUG */

    if (ntor > 0) {
	pr( logFile, "\n\nList of Internal Non-Bonded Interactions:\n" );
	pr( logFile, "_________________________________________\n\n" );
	pr( logFile, "First Second Non-Bonded\n" );
	pr( logFile, "Atom  Atom(s)\n" );
	pr( logFile, "_____ " );
    } /* endif */

    for (i = 0;  i < MAX_ATOMS;  i++) {
	if (Nnbonds[i] != 0) {
	    for ( j = 0;  j < OUTNUMATM;  j++ ) {
		pr( logFile, "_____" );
	    } /* j */
	    break;
	} /* endif */
    } /* i */

    if (ntor > 0) {
	pr( logFile, "_\n" );
    }
    flushLog;

    for (i = 0;  i < MAX_ATOMS;  i++) {
	if (Nnbonds[i] != 0) {
	    repflag = FALSE;
	    n_a = 0;
	    pr( logFile, "%4s: ", pdbaname[i] );
	    for ( j = 0; j < Nnbonds[i]; j++) {
		++n_a; 
		if (n_a >= OUTNUMATM) {
		    n_a = 0;
		    pr( logFile, "\n      " );
		}
		pr( logFile, "%s", pdbaname[ nbmatrix_binary[j][i] ] );
		for ( k = 0; k < j; k++ ) {
		    if (nbmatrix_binary[k][i] == nbmatrix_binary[j][i]) {
			repflag = TRUE;
			break;
		    }
		} /*  k  */
		if (repflag&&(k != j)) {
		    pr( logFile,"(ERROR! %d & %d)", j, k );
		}
		pr( logFile, "%s", (j == (Nnbonds[i]-1))?".":", " );
	    } /*  j  */
	    pr( logFile,"\n\n");
	} /*  endif  */
    } /*  i  */

    pr( logFile, "\nInternal Non-bonded Interactions before,\t%d\n", (natom+1)*natom/2);
    pr( logFile, "                       and after weeding =\t%d\n\n", Nnb);
}
/* EOF */
