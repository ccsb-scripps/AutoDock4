
#define ENERGYPENALTY 500.0	/* Energy factor which is multiplied by distance
				   from centre of grid, to penalize atoms
				   outside grid */

#ifndef QUICKTRILINTERP
#define QUICKTRILINTERP
#include "constants.h"
float  quicktrilinterp( CONST_FLOAT tcoord[MAX_ATOMS][SPACE], 
            CONST_FLOAT charge[MAX_ATOMS], 
            CONST_INT   type[MAX_ATOMS], 
            CONST_INT   natom, 
            CONST_FLOAT map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS],
            CONST_FLOAT inv_spacing, 
            CONST_FLOAT xlo, 
            CONST_FLOAT ylo, 
            CONST_FLOAT zlo );
#endif

#ifndef TRILINTERP
#define TRILINTERP
#include "constants.h"
float  trilinterp(CONST_FLOAT tcoord[MAX_ATOMS][SPACE], 
            CONST_FLOAT charge[MAX_ATOMS], 
            CONST_INT   type[MAX_ATOMS], 
            CONST_INT   natom, 
            CONST_FLOAT map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS],
            CONST_FLOAT inv_spacing, 
                  float elec[MAX_ATOMS],
                  float evdW[MAX_ATOMS],
            CONST_FLOAT xlo, 
            CONST_FLOAT ylo, 
            CONST_FLOAT zlo );
#endif

#ifndef OUTSIDETRILINTERP
#define OUTSIDETRILINTERP
#include "constants.h"
float  outsidetrilinterp( CONST_FLOAT tcoord[MAX_ATOMS][SPACE], 
            CONST_FLOAT charge[MAX_ATOMS], 
            CONST_INT   type[MAX_ATOMS], 
            CONST_INT   natom, 
            CONST_FLOAT map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS],
            CONST_FLOAT inv_spacing, 
	    // float elec[MAX_ATOMS],
	    // float emap[MAX_ATOMS],
            CONST_FLOAT xlo, 
            CONST_FLOAT ylo, 
            CONST_FLOAT zlo,
            CONST_FLOAT xhi, 
            CONST_FLOAT yhi, 
            CONST_FLOAT zhi,
            CONST_FLOAT xcen, 
            CONST_FLOAT ycen, 
            CONST_FLOAT zcen );
#endif

#ifndef TEMPLATETRILINTERP
#define TEMPLATETRILINTERP
#include "constants.h"
float  template_trilinterp( CONST_FLOAT tcoord[MAX_ATOMS][SPACE], 
            CONST_FLOAT charge[MAX_ATOMS], 
            CONST_INT   type[MAX_ATOMS], 
            CONST_INT   natom, 
            CONST_FLOAT map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS],
            CONST_FLOAT inv_spacing, 
            CONST_FLOAT xlo, 
            CONST_FLOAT ylo, 
            CONST_FLOAT zlo,
            CONST_FLOAT template_energy[MAX_ATOMS],
            CONST_FLOAT template_stddev[MAX_ATOMS]);
#endif

#ifndef OUTSIDETEMPLTRILINTERP
#define OUTSIDETEMPLTRILINTERP
#include "constants.h"
float  outside_templ_trilinterp( CONST_FLOAT tcoord[MAX_ATOMS][SPACE], 
            CONST_FLOAT charge[MAX_ATOMS], 
            CONST_INT   type[MAX_ATOMS], 
            CONST_INT   natom, 
            CONST_FLOAT map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS],
            CONST_FLOAT inv_spacing, 
            CONST_FLOAT xlo, 
            CONST_FLOAT ylo, 
            CONST_FLOAT zlo,
            CONST_FLOAT xhi, 
            CONST_FLOAT yhi, 
            CONST_FLOAT zhi,
            CONST_FLOAT xcen, 
            CONST_FLOAT ycen, 
            CONST_FLOAT zcen,
            CONST_FLOAT template_energy[MAX_ATOMS],
            CONST_FLOAT template_stddev[MAX_ATOMS]);
#endif

#ifndef BYATOM_TEMPLATE_TRILINTERP
#define BYATOM_TEMPLATE_TRILINTERP
#include "constants.h"
float  byatom_template_trilinterp( CONST_FLOAT tcoord[MAX_ATOMS][SPACE], 
                                   CONST_FLOAT charge[MAX_ATOMS], 
                                   CONST_INT   type[MAX_ATOMS], 
                                   CONST_INT   natom, 
                                   CONST_FLOAT map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS],
                                   CONST_FLOAT inv_spacing, 
		                           float elec[MAX_ATOMS], 
		                           float emap[MAX_ATOMS], 
                                   CONST_FLOAT xlo, 
                                   CONST_FLOAT ylo, 
                                   CONST_FLOAT zlo,
                                   CONST_FLOAT template_energy[MAX_ATOMS],
                                   CONST_FLOAT template_stddev[MAX_ATOMS]);
#endif
