/*

 $Id: set_cmd_io_std.cc,v 1.3 2004/11/16 23:42:53 garrett Exp $

*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

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
