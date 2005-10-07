
#define ENERGYPENALTY 500.0	/* Energy factor which is multiplied by distance
				   from centre of grid, to penalize atoms
				   outside grid */
#ifndef TRILINTERP
#define TRILINTERP

#include "constants.h"
#include "structs.h"

FloatOrDouble  trilinterp(CONST_FLOAT tcoord[MAX_ATOMS][SPACE], 
            CONST_FLOAT charge[MAX_ATOMS], 
            CONST_FLOAT abs_charge[MAX_ATOMS], 
            CONST_INT   type[MAX_ATOMS], 
            CONST_INT   natom, 
            CONST_FLOAT map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS],
                  FloatOrDouble elec[MAX_ATOMS],
                  FloatOrDouble evdW[MAX_ATOMS],
            GridMapSetInfo *info );
#endif

#ifndef QUICKTRILINTERP
#define QUICKTRILINTERP
#include "constants.h"
FloatOrDouble  quicktrilinterp( CONST_FLOAT tcoord[MAX_ATOMS][SPACE], 
            CONST_FLOAT charge[MAX_ATOMS], 
            CONST_FLOAT abs_charge[MAX_ATOMS], 
            CONST_INT   type[MAX_ATOMS], 
            CONST_INT   natom, 
            CONST_FLOAT map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS],
            GridMapSetInfo *info );
#endif


#ifndef OUTSIDETRILINTERP
#define OUTSIDETRILINTERP
#include "constants.h"
FloatOrDouble  outsidetrilinterp( CONST_FLOAT tcoord[MAX_ATOMS][SPACE], 
            CONST_FLOAT charge[MAX_ATOMS], 
            CONST_FLOAT abs_charge[MAX_ATOMS], 
            CONST_INT   type[MAX_ATOMS], 
            CONST_INT   natom, 
            CONST_FLOAT map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS],
	    // FloatOrDouble elec[MAX_ATOMS],
	    // FloatOrDouble emap[MAX_ATOMS],
            GridMapSetInfo *info );
#endif

#ifndef OUTSIDETRILINTERPBYATOM
#define OUTSIDETRILINTERPBYATOM
#include "constants.h"
FloatOrDouble  outsidetrilinterpbyatom( CONST_FLOAT tcoord[MAX_ATOMS][SPACE], 
            CONST_FLOAT charge[MAX_ATOMS], 
            CONST_FLOAT abs_charge[MAX_ATOMS], 
            CONST_INT   type[MAX_ATOMS], 
            CONST_INT   natom, 
            CONST_FLOAT map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS],
	        FloatOrDouble elec[MAX_ATOMS],
	        FloatOrDouble emap[MAX_ATOMS],
            GridMapSetInfo *info );
#endif

#ifndef TEMPLATETRILINTERP
#define TEMPLATETRILINTERP
#include "constants.h"
FloatOrDouble  template_trilinterp( CONST_FLOAT tcoord[MAX_ATOMS][SPACE], 
            CONST_FLOAT charge[MAX_ATOMS], 
            CONST_FLOAT abs_charge[MAX_ATOMS], 
            CONST_INT   type[MAX_ATOMS], 
            CONST_INT   natom, 
            CONST_FLOAT map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS],
            CONST_FLOAT template_energy[MAX_ATOMS],
            CONST_FLOAT template_stddev[MAX_ATOMS],
            GridMapSetInfo *info );
#endif

#ifndef OUTSIDETEMPLTRILINTERP
#define OUTSIDETEMPLTRILINTERP
#include "constants.h"
FloatOrDouble  outside_templ_trilinterp( CONST_FLOAT tcoord[MAX_ATOMS][SPACE], 
            CONST_FLOAT charge[MAX_ATOMS], 
            CONST_FLOAT abs_charge[MAX_ATOMS], 
            CONST_INT   type[MAX_ATOMS], 
            CONST_INT   natom, 
            CONST_FLOAT map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS],
            CONST_FLOAT template_energy[MAX_ATOMS],
            CONST_FLOAT template_stddev[MAX_ATOMS],
            GridMapSetInfo *info );
#endif

#ifndef BYATOM_TEMPLATE_TRILINTERP
#define BYATOM_TEMPLATE_TRILINTERP
#include "constants.h"
FloatOrDouble  byatom_template_trilinterp( CONST_FLOAT tcoord[MAX_ATOMS][SPACE], 
                                   CONST_FLOAT charge[MAX_ATOMS], 
                                   CONST_FLOAT abs_charge[MAX_ATOMS], 
                                   CONST_INT   type[MAX_ATOMS], 
                                   CONST_INT   natom, 
                                   CONST_FLOAT map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS],
		                           FloatOrDouble elec[MAX_ATOMS], 
		                           FloatOrDouble emap[MAX_ATOMS], 
                                   CONST_FLOAT template_energy[MAX_ATOMS],
                                   CONST_FLOAT template_stddev[MAX_ATOMS],
                                   GridMapSetInfo *info );
#endif

#ifndef TRILINTERP4
#define TRILINTERP4
#include "constants.h"
FloatOrDouble  trilinterp4(CONST_FLOAT tcoord[MAX_ATOMS][SPACE], 
            CONST_FLOAT charge[MAX_ATOMS], 
            CONST_FLOAT abs_charge[MAX_ATOMS], 
            CONST_INT   type[MAX_ATOMS], 
            CONST_INT   natom, 
            CONST_FLOAT map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS],
                  FloatOrDouble elec[MAX_ATOMS],
                  FloatOrDouble evdW[MAX_ATOMS],
            int ignore_inter[MAX_ATOMS],
            GridMapSetInfo *info );
#endif

#ifndef OUTSIDETRILINTERP4BYATOM
#define OUTSIDETRILINTERP4BYATOM
#include "constants.h"
FloatOrDouble  outsidetrilinterp4byatom( CONST_FLOAT tcoord[MAX_ATOMS][SPACE], 
            CONST_FLOAT charge[MAX_ATOMS], 
            CONST_FLOAT abs_charge[MAX_ATOMS], 
            CONST_INT   type[MAX_ATOMS], 
            CONST_INT   natom, 
            CONST_FLOAT map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS],
	    FloatOrDouble elec[MAX_ATOMS],
	    FloatOrDouble emap[MAX_ATOMS],
            int ignore_inter[MAX_ATOMS],
            GridMapSetInfo *info );
#endif

#ifndef OUTSIDETRILINTERP4
#define OUTSIDETRILINTERP4
#include "constants.h"
FloatOrDouble  outsidetrilinterp4( CONST_FLOAT tcoord[MAX_ATOMS][SPACE], 
            CONST_FLOAT charge[MAX_ATOMS], 
            CONST_FLOAT abs_charge[MAX_ATOMS], 
            CONST_INT   type[MAX_ATOMS], 
            CONST_INT   natom, 
            CONST_FLOAT map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS],
            int ignore_inter[MAX_ATOMS],
            GridMapSetInfo *info );
#endif

#ifndef QUICKTRILINTERP4
#define QUICKTRILINTERP4
#include "constants.h"
FloatOrDouble  quicktrilinterp4( CONST_FLOAT tcoord[MAX_ATOMS][SPACE], 
            CONST_FLOAT charge[MAX_ATOMS], 
            CONST_FLOAT abs_charge[MAX_ATOMS], 
            CONST_INT   type[MAX_ATOMS], 
            CONST_INT   natom, 
            CONST_FLOAT map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS],
            int ignore_inter[MAX_ATOMS],
            GridMapSetInfo *info );
#endif

