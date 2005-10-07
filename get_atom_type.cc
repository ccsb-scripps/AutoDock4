/*

 $Id: get_atom_type.cc,v 1.4 2005/08/15 23:00:09 garrett Exp $

*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/* get_atom_type.cc */

#include <stdio.h>
#include <string.h>
#include "get_atom_type.h"
#include "structs.h"
#include "autocomm.h"
#include "print_2x.h"
#include "atom_parameter_manager.h"

extern char *programname;
extern FILE *logFile;
extern int debug;


int get_atom_type( char aname[4] )
{
    ParameterEntry * found_parm;
    ParameterEntry thisparm;
    int map_index = -1;
    int bond_index = -1;
    char message[MAX_CHARS];

    // "map_index" is used as an index into the "map" array,
    // to look up the correct energies in the current grid cell,
    // thus:  map[][][][map_index]

    // AutoGrid 4 Typing
    strcpy(thisparm.autogrid_type, aname);
    found_parm = apm_find(thisparm.autogrid_type);
    if (found_parm != NULL) {
        map_index = found_parm->map_index;
        bond_index = found_parm->bond_index;
        if (debug > 0) {
            (void) fprintf( logFile, "Found parameters for ligand atom named %s, atom type \"%s\", grid map index = %d\n",
                    aname, found_parm->autogrid_type, found_parm->map_index );
        }

    } else {
        // We could not find this parameter -- return error here
         prStr( message, "\n%s: *** WARNING!  Unknown ligand atom type \"%s\" found.  You should add  parameters for it to the parameter library first! ***\n\n", programname, aname);
         pr_2x( stderr, logFile, message );

    }
    // /--AutoGrid4 
    
    if (map_index == -1) {
        pr( logFile, "%s: WARNING:  atom type not found, using the default, atom type = 1, instead.\n", programname);
        map_index = 0;  // we are 0-based internally, 1-based in printed output, for human's sake
        bond_index = 0;  // we are 0-based internally, 1-based in printed output, for human's sake
    }
    return (map_index);
}


/* EOF */
