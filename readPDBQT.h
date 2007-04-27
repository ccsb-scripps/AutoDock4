/*

 $Id: readPDBQT.h,v 1.10 2007/04/27 06:01:50 garrett Exp $

 AutoDock 

 Copyright (C) 1989-2007,  Garrett M. Morris, David S. Goodsell, Ruth Huey, Arthur J. Olson, 
 All Rights Reserved.

 AutoDock is a Trade Mark of The Scripps Research Institute.

 This program is free software; you can redistribute it and/or
 modify it under the terms of the GNU General Public License
 as published by the Free Software Foundation; either version 2
 of the License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

 */

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

void  readPDBQTLine( char line[LINE_LEN],
                     int  *ptr_serial,
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
