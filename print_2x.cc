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
