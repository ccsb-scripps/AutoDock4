
#ifndef CLMODE
#define CLMODE
#include "constants.h"
#include "strindex.h"
#include "readPDBQT.h"
#include "get_atom_type.h"
#include "getpdbcrds.h"
#include "sort_enrg.h"
#include "cluster_analysis.h"
#include "prClusterHist.h"
#include "bestpdb.h"
#include "success.h"
#include "qmultiply.h"
#include "openfile.h"
void  clmode( int   num_atm_maps, 
              FloatOrDouble clus_rms_tol, 
              char  hostnm[MAX_CHARS], 
              Clock jobStart,
              struct tms tms_jobStart, 
              Boole B_write_all_clusmem, 
              char  clusFN[MAX_CHARS], 
              FloatOrDouble crdpdb[MAX_ATOMS][SPACE], 
              FloatOrDouble sml_center[SPACE], 
              Boole B_symmetry_flag,
              char  rms_ref_crds[MAX_CHARS] );
#endif
