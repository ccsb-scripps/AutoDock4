#ifndef MAIN
#define MAIN

#include "constants.h"
#include "setflags.h"
#include "stop.h"
#include "banner.h"
#include "printdate.h"
#include "parse_dpf_line.h"
#include "print_2x.h"
#include "strindex.h"
#include "stateLibrary.h"
#include "timesyshms.h"
#include "cnv_state_to_coords.h"
#include "dpftypes.h"
#include "readfield.h"
#include "readmap.h"
#include "cmdmode.h"
#include "success.h"
#include "readPDBQ.h"
#include "intnbtable.h"
#include "nbe.h"
#include "clmode.h"
#include "simanneal.h"
#include "eintcal.h"
#include "writePDBQ.h"
#include "evaluate_energy.h"
#include "analysis.h"
#include "printEnergies.h"
#include "investigate.h"
#include "qmultiply.h"

int  main( int  argc, char **argv, char **envp);

#endif

