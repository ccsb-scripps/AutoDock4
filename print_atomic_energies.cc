/*

 $Id: print_atomic_energies.cc,v 1.4 2007/04/27 06:01:50 garrett Exp $

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

#include <stdio.h>
#include <string.h>
#include "print_atomic_energies.h"


extern FILE *logFile;

/*----------------------------------------------------------------------------*/
void print_atomic_energies( int natom, 
			    char atomstuff[MAX_ATOMS][MAX_CHARS],
			    int type[MAX_ATOMS],
			    Real emap[MAX_ATOMS],
			    Real elec[MAX_ATOMS],
			    Real charge[MAX_ATOMS] )

/*----------------------------------------------------------------------------*/
{
    char rec[16];
    register int i;

    fprintf( logFile, 
     "Small Molecule  Atom  Non-bonded   Electrostatic  Partial\n" );
    fprintf( logFile, 
     "Atom & Residue  Type    Energy        Energy      Charge\n" );
    fprintf( logFile, 
     "______________  ____  ___________  _____________  ______\n" );

    for (i = 0;  i < natom;  i++) {
        strncpy( rec, &atomstuff[i][6], (size_t)14 );
        fprintf( logFile, "%.14s   %1d  %+11.2f  %+11.2f      %+6.3f\n", rec, (type[i]+1), emap[i], elec[i], charge[i] );
    }
}
/* EOF */
