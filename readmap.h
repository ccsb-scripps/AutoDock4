
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
FloatOrDouble   mapc2f( char C_mapValue );
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
                FloatOrDouble *P_ExtSpacing, 
                char  AtmTypStr[ATOM_MAPS], 
                char  ExtFldFileName[MAX_CHARS], 
                int   ExtGridPts1[SPACE], 
                int   ExtGridPts[SPACE], 
                Clock jobStart, 
                char  line[LINE_LEN], 
                char  ExtMacromolFileName[MAX_CHARS], 
                FloatOrDouble map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS], 
                FloatOrDouble MapCenter[SPACE], 
                FloatOrDouble MapMax[MAX_MAPS], 
                FloatOrDouble MapMin[MAX_MAPS], 
		struct tms tmsJobStart,
		Boole B_charMap);
#endif
