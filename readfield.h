
#ifndef READFIELD
#define READFIELD
#include "constants.h"
#include "openfile.h"
#include "stop.h"
void    readfield( FloatOrDouble *P_inv_spacing, 
                FloatOrDouble *P_spacing, 
                char  gdfldFileName[MAX_CHARS], 
                char  gpfFileName[MAX_CHARS], 
                int   gridpts1[SPACE], 
                int   gridpts[SPACE], 
		FloatOrDouble *xhi,
		FloatOrDouble *yhi,
		FloatOrDouble *zhi,
                Clock jobStart, 
                char  line[LINE_LEN], 
                FloatOrDouble *xlo, 
                FloatOrDouble *ylo, 
                FloatOrDouble *zlo, 
                char  macromolFileName[MAX_CHARS], 
                FloatOrDouble maP_center[SPACE], 
		struct tms tms_jobStart );
#endif
