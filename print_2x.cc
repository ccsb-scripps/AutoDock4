/*

 $Id: print_2x.cc,v 1.2 2003/02/26 01:27:35 garrett Exp $

*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/* print_2x.cc */

    #include <stdio.h>
    #include "print_2x.h"

void
print_2x( FILE *stream1,
	  FILE *stream2,
	  char *string )

{
	fprintf( stream1, string );
	fprintf( stream2, string );
}
/* EOF */
