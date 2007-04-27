/*

 $Id: mkTorTree.h,v 1.7 2007/04/27 06:01:50 garrett Exp $

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

#ifndef MKTORTREE
#define MKTORTREE

#include "constants.h"
#include "parse_PDBQT_line.h"
#include "stop.h"

void  mkTorTree(int   atomnumber[MAX_RECORDS],
                char  record[MAX_RECORDS][LINE_LEN],
                int   nrecord,
                int   tlist[MAX_TORS][MAX_ATOMS],
                int   *P_ntor,
                int   *P_ntor_ligand,
                char  smFileName[MAX_CHARS],
                char  pdbaname[MAX_ATOMS][5],
                Boole *P_B_constrain,
                int   *P_atomC1,
                int   *P_atomC2,
                Real *P_sqlower,
                Real *P_squpper,
                int   *P_ntorsdof,
                int   ignore_inter[MAX_ATOMS]);
#endif
