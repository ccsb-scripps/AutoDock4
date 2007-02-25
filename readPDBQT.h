#ifndef READPDBQT
#define READPDBQT

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
void  readPDBQTLine( char  line[LINE_LEN],
                     Real crd[SPACE], 
                     Real *P_q,
                     ParameterEntry *thisparm);

Molecule readPDBQT( char  line[LINE_LEN],

              int   num_atm_maps,

              int   *P_natom,
              Real crdpdb[MAX_ATOMS][NTRN],
              Real crdreo[MAX_ATOMS][NTRN],
              Real charge[MAX_ATOMS],
              Boole *P_B_haveCharges,
              int   type[MAX_ATOMS],
              int   bondtype[MAX_ATOMS],
              char  pdbaname[MAX_ATOMS][5],

              char  pdbqFileName[MAX_CHARS],
              char FN_flexres[MAX_CHARS],
              Boole B_have_flexible_residues,

              char  atomstuff[MAX_ATOMS][MAX_CHARS],
              int   Htype,

              Boole *P_B_constrain,
              int   *P_atomC1,
              int   *P_atomC2,
              Real *P_sqlower,
              Real *P_squpper,

              int   *P_ntor1,
              int   *P_ntor,
              int *P_ntor_ligand,
              int   tortree[MAX_TORS][MAX_ATOMS],
              Real vt[MAX_TORS][NTRN],

              int   *P_Nnb,
              NonbondParam *nonbondlist,

              Clock jobStart,
              struct tms tms_jobStart,
              char  hostnm[MAX_CHARS],

              int   *P_ntorsdof,
              int   outlev,

              int   ignore_inter[MAX_ATOMS],
              
              int   B_include_1_4_interactions,
              
              Atom  atoms[MAX_ATOMS],
              char  PDBQT_record[MAX_RECORDS][LINE_LEN]
              );
#endif
