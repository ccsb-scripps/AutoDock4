/*

 $Id: get_atom_type.cc,v 1.3 2005/03/11 02:11:30 garrett Exp $

*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/* get_atom_type.cc */

#include <stdio.h>
#include <search.h>
#include "get_atom_type.h"
#include "structs.h"
#include "autocomm.h"
#include "print_2x.h"

extern char *programname;
extern FILE *logFile;
extern int debug;


int get_atom_type( char aname[4] )
{
    ENTRY item, *found_item; /*see hsearch(3C)*/
    ParameterEntry * found_parm;
    ParameterEntry thisparm;
    int map_index = -1;
    int bond_index = -1;
    char message[MAX_CHARS];


    // "map_index" is used as an index into the "map" array,
    // to look up the correct energies in the current grid cell,
    // thus:  map[][][][map_index]

    // AutoGrid 4 Typing --/
    item.key = thisparm.autogrid_type;
    if ((found_item = hsearch (item, FIND)) != NULL) {
        found_parm = (ParameterEntry *)found_item->data;
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
