
#ifndef MAPC2F
#define MAPC2F
#include "constants.h"
#include "openfile.h"
#include "warn_bad_file.h"
#include "strindex.h"
#include "print_2x.h"
#include "check_header_line.h"
#include "warn_bad_file.h"
#include "check_header_float.h"
#include "check_header_int.h"
#include "timesys.h"
float   mapc2f( char C_mapValue );
#endif

#ifndef READMAP
#define READMAP
#include "constants.h"
#include "openfile.h"
#include "warn_bad_file.h"
#include "strindex.h"
#include "print_2x.h"
#include "check_header_line.h"
#include "warn_bad_file.h"
#include "check_header_float.h"
#include "check_header_int.h"
#include "timesys.h"
void    readmap( Boole *P_B_HaveMap, 
                int   *P_Imap,
                int   *P_NumAtmMaps, 
                float *P_ExtSpacing, 
                char  AtmTypStr[ATOM_MAPS], 
                char  ExtFldFileName[MAX_CHARS], 
                int   ExtGridPts1[SPACE], 
                int   ExtGridPts[SPACE], 
                Clock jobStart, 
                char  line[LINE_LEN], 
                char  ExtMacromolFileName[MAX_CHARS], 
                float map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS], 
                float MapCenter[SPACE], 
                float MapMax[MAX_MAPS], 
                float MapMin[MAX_MAPS], 
		struct tms tmsJobStart,
		Boole B_charMap);
#endif
