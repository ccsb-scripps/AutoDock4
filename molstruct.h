/*

 $Id: molstruct.h,v 1.9 2009/05/08 23:02:14 rhuey Exp $

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

/* structure for molecules; side-chains, etc. */

typedef struct molecule {

	int   natom;
	double untrnfm_crdpdb[ MAX_ATOMS ][ NTRN ];
	double trnsfmd_crdpdb[ MAX_ATOMS ][ NTRN ];
	double charge[ MAX_ATOMS ];
	int   type[ MAX_ATOMS ];
	char  pdbaname[ MAX_ATOMS ][ 5 ];
	char  atomstuff[ MAX_ATOMS ][ MAX_CHARS ];

	Boole B_haveCharges;
	char  pdbqFileName[ PATH_MAX ];
	int   Htype;
	Boole B_constrain;
	int   atomC1;
	int   atomC2;
	double sqlower;
	double squpper;

	int   ntor1;
	int   ntor;
	int   tlist[ MAX_TORS ][ MAX_ATOMS ];
	double vt[ MAX_TORS ][ NTRN ];

	int   Nnb;
	int   Nnbonds[ MAX_ATOMS ];
	NonbondParam *nonbondlist;
} Molecule;
