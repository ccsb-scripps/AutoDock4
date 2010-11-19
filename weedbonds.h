/*

 $Id: weedbonds.h,v 1.9.2.1 2010/11/19 20:09:29 rhuey Exp $

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

#ifndef WEEDBONDS
#define WEEDBONDS

#include "constants.h"
#include "stop.h"
#include "structs.h"

void  weedbonds( int   natom,
                 char  pdbaname[MAX_ATOMS][5],
                 int   piece[MAX_ATOMS],
                 int   ntor,
                 int   tlist[MAX_TORS][MAX_ATOMS],
                 int   nbmatrix_binary[MAX_ATOMS][MAX_ATOMS],
                 int   *P_Nnb,
                 NonbondParam *nonbondlist,
                 int   outlev,
				 int   type[MAX_ATOMS],
				 int   natom_in_ligand[MAX_LIGANDS],
				 int   nligand);
#endif


#ifndef PRINT_NONBONDS
#define PRINT_NONBONDS

#include "constants.h"
#include "stop.h"
#include "structs.h"

void print_nonbonds(
                int natom,
                char pdbaname[MAX_ATOMS][5],
                int piece[MAX_ATOMS],
                int ntor,
                int tlist[MAX_TORS][MAX_ATOMS],
                int nbmatrix[MAX_ATOMS][MAX_ATOMS],
                int Nnb,
                NonbondParam *nonbondlist,
                int outlev,
                int type[MAX_ATOMS]);
#endif

