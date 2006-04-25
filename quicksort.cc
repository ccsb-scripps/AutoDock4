/*

 $Id: quicksort.cc,v 1.3 2006/04/25 22:33:01 garrett Exp $

*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/* quicksort.cc */

    #include "quicksort.h"
#ifdef DEBUG
#include <stdio.h>
#include <string.h>
#endif /* DEBUG */


/* 
\  Based on C.A.R.Hoare's original
 \ algorithm of 1962...
*/

#ifdef DEBUG
extern FILE *logFile;
#endif /* DEBUG */

void quicksort( Real e[], 
		int isort[],
		int left,
		int right )

{
    int i, last;

#ifdef DEBUG
            char array[101];
#endif  /* DEBUG */

    if (left >= right) {
	return;
    }

#ifdef DEBUG
            strncpy( array, "----------------------------------------------------------------------------------------------------", (size_t)100 );
            array[100] = '\0';
            for (i=left+1; i<right; i++) array[i]=' ';
            array[left] = 'L';
            array[right] = 'R';
            fprintf( logFile, "%s\n", array );
#endif  /* DEBUG */

    swap( isort, left, (left+right)/2 );

    last = left;

    for (i = left+1; i <= right; i++) {
	if (e[isort[i]] < e[isort[left]]) {
	    swap( isort, ++last, i );
	}
    }

    swap( isort, (int)left, (int)last );

    quicksort( e, isort, left,  last-1 );
    quicksort( e, isort, last+1, right  );
}
/* EOF */
