/*

 $Id: swap.cc,v 1.2 2003/02/26 01:46:44 garrett Exp $

*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/* swap.cc */

    #include "swap.h"
    #ifdef DEBUG
	#include <stdio.h>
	#include <string.h>
    #endif /* DEBUG */


#ifdef DEBUG
extern FILE *logFile;
#endif /* DEBUG */

void swap ( int v[],
	    int i, 
	    int j )

{
    int temp;

#ifdef DEBUG
    int k;
    char array[101];
    strncpy( array, "----------------------------------------------------------------------------------------------------", (size_t)100 );
    array[100] = '\0';
    for (k=i+1; k<j; k++) array[k]=' ';
    array[i] = '<';
    array[j] = '>';
    fprintf( logFile, "%s", array );
    fprintf( logFile, " swapping %d & %d.\n", i, j);
#endif /* DEBUG */

    temp = v[i];
    v[i] = v[j];
    v[j] = temp;
}
/* EOF */
