#ifndef DPFTYPES
#define DPFTYPES

#include "constants.h"
#include "stop.h"

void dpftypes( int   *P_Htype, 
                  int   *P_num_all_maps,
                  int   *P_num_atm_maps, 
                  char  atm_tyP_str[ATOM_MAPS], 
                  char  line[LINE_LEN] );
#endif
