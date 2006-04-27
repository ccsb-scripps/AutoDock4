#ifndef MAIN
#define MAIN

#include "analysis.h"
#include "assert.h"
#include "atom_parameter_manager.h"
#include "autoglobal.h"
#include "banner.h"
#include "clmode.h"
#include "cmdmode.h"
#include "cnv_state_to_coords.h"
#include "constants.h"
#include "eintcal.h"
#include "evaluate_energy.h"
#include "intnbtable.h"
#include "investigate.h"
#include "nbe.h"
#include "parse_dpf_line.h"
#include "parse_param_line.h"
#include "parsetypes.h"
#include "printEnergies.h"
#include "print_2x.h"
#include "printdate.h"
#include "qmultiply.h"
#include "readPDBQT.h"
#include "readfield.h"
#include "readmap.h"
#include "read_parameter_library.h"
#include "setflags.h"
#include "simanneal.h"
#include "stateLibrary.h"
#include "stop.h"
#include "strindex.h"
#include "structs.h"
#include "success.h"
#include "timesyshms.h"
#include "writePDBQT.h"

#define UNBOUND 0
#define DOCKED 1

//int  main( int  argc, char **argv, char **envp);
int main (int argc, char * const argv[], char * const envp[]);

#endif

