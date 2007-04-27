/*

 $Id: get_atom_type.cc,v 1.5 2007/04/27 06:01:48 garrett Exp $

 AutoDock 

 Copyright (C) 1989-2007,  Garrett M. Morris, David S. Goodsell, Ruth Huey, Arthur J. Olson, 
 All Rights Reserved.

 AutoDock is a Trade Mark of The Scripps Research Institute.

 This program is free software; you can redistribute it and/or
 modify it under the terms of the GNU General Public License
 as published by the Free Software Foundation; either version 2
 of the License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

 */


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

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
