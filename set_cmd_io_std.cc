/* set_cmd_io_std.cc */

    #include <stdio.h>
    #include "set_cmd_io_std.h"


extern FILE *command_in_fp;
extern FILE *command_out_fp;

/*----------------------------------------------------------------------------*/

void set_cmd_io_std( void )

/*----------------------------------------------------------------------------*/
{
    command_in_fp = stdin;
    command_out_fp = stdout;
}
/* EOF */
