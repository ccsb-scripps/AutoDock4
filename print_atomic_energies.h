/*

 $Id: print_atomic_energies.h,v 1.5 2009/05/08 23:02:15 rhuey Exp $

 AutoDock 

Copyright (C) 2009 The Scripps Research Institute. All rights reserved.

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

#ifndef PRINT_ATOMIC_ENERGIES
#define PRINT_ATOMIC_ENERGIES

#include "constants.h"

void  print_atomic_energies(int   natom,
                            char  atomstuff[MAX_ATOMS][MAX_CHARS],
                            int   type[MAX_ATOMS],
                            Real emap[MAX_ATOMS],
                            Real elec[MAX_ATOMS],
                            Real charge[MAX_ATOMS] );
#endif
