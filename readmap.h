
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
Real   mapc2f( char C_mapValue );
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
#include "structs.h"

void readmap( char line[LINE_LEN],
             int outlev,

             Clock jobStart,
             struct tms tmsJobStart,
        
             Boole B_charMap,

             Boole *P_B_HaveMap, 
             int *P_imap, 
             
             GridMapSetInfo *info,
             Real map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS]
             // double *maps 
             );

#endif
