
#ifndef READFIELD
#define READFIELD
#include "constants.h"
#include "openfile.h"
#include "stop.h"
void    readfield( float *P_inv_spacing, 
                float *P_spacing, 
                char  gdfldFileName[MAX_CHARS], 
                char  gpfFileName[MAX_CHARS], 
                int   gridpts1[SPACE], 
                int   gridpts[SPACE], 
		float *xhi,
		float *yhi,
		float *zhi,
                Clock jobStart, 
                char  line[LINE_LEN], 
                float *xlo, 
                float *ylo, 
                float *zlo, 
                char  macromolFileName[MAX_CHARS], 
                float maP_center[SPACE], 
		struct tms tms_jobStart );
#endif
