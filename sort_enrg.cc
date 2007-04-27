/*

 $Id: sort_enrg.cc,v 1.5 2007/04/27 06:01:51 garrett Exp $

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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "sort_enrg.h"

extern FILE *logFile;

void sort_enrg( Real econf[MAX_RUNS],
                int isort[MAX_RUNS],
		int nconf )

{
/*__________________________________________________________________________
 | Sort conformations on energy                                             |
 |__________________________________________________________________________|
 | Searches through all conformations;  puts in isort[0] the index of the   |
 | lowest energy, in isort[1] the next lowest energy's index, and so on.    |
 |__________________________________________________________________________|
 | WARNING: Fails if any 2 or more econf[] energies are equal.              |
 |__________________________________________________________________________|*/

    quicksort( econf, isort, 0, nconf-1 );
}
/* EOF */
