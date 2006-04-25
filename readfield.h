
#ifndef READFIELD
#define READFIELD
#include "constants.h"
#include "openfile.h"
#include "stop.h"
#include "structs.h"
/*
void    readfield( Real *P_inv_spacing, 
                Real *P_spacing, 
                char  gdfldFileName[MAX_CHARS], 
                char  gpfFileName[MAX_CHARS], 
                int   gridpts1[SPACE], 
                int   gridpts[SPACE], 
		Real *xhi,
		Real *yhi,
		Real *zhi,
                Clock jobStart, 
                char  line[LINE_LEN], 
                Real *xlo, 
                Real *ylo, 
                Real *zlo, 
                char  macromolFileName[MAX_CHARS], 
                Real maP_center[SPACE], 
		struct tms tms_jobStart );
*/

void readfield( GridMapSetInfo *info, // *ptr_map_set_info
                char line[LINE_LEN],
                Clock jobStart,
                struct tms tms_jobStart );


#endif
