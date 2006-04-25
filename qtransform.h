
#ifndef QTRANSFORM
#define QTRANSFORM
#include "constants.h"
#include "structs.h"
void qtransform( const Coord T,
	 	 const Quat  q,
                 Real tcoord[MAX_ATOMS][SPACE],
		 const int   natom);
#endif
