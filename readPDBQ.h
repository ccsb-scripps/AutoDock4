#ifndef READPDBQ
#define READPDBQ

#include "structs.h"
#include "constants.h"
#include "stop.h"
#include "get_atom_type.h"
#include "print_2x.h"
#include "mkTorTree.h"
#include "nonbonds.h"
#include "weedbonds.h"
#include "torNorVec.h"
#include "success.h"
#include "openfile.h"

#include "constants.h"
void  readPDBQLine( char  line[LINE_LEN],
                float crd[SPACE], 
                float *P_q );

Molecule readPDBQ( char  line[LINE_LEN],

	      char  atm_typ_str[ATOM_MAPS],
	      int   num_atm_maps,

	      int   *P_natom,
	      float crdpdb[MAX_ATOMS][NTRN],
	      float charge[MAX_ATOMS],
	      Boole *P_B_haveCharges,
	      int   type[MAX_ATOMS],
	      char  pdbaname[MAX_ATOMS][5],
	      char  pdbqFileName[MAX_CHARS],
	      char  atomstuff[MAX_ATOMS][MAX_CHARS],
	      int   Htype,

	      Boole *P_B_constrain,
	      int   *P_atomC1,
	      int   *P_atomC2,
	      float *P_sqlower,
	      float *P_squpper,

	      int   *P_ntor1,
	      int   *P_ntor,
	      int   tortree[MAX_TORS][MAX_ATOMS],
	      float vt[MAX_TORS][NTRN],

	      int   *P_Nnb,
	      int   Nnbonds[MAX_ATOMS],
	      int   nonbondlist[MAX_NONBONDS][2],

	      Clock jobStart,
	      struct tms tms_jobStart,
	      char  hostnm[MAX_CHARS],
        int   *P_ntorsdof );
#endif
