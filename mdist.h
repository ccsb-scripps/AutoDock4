/*

 $Id: mdist.h,v 1.6 2007/04/27 06:01:49 garrett Exp $

 AutoDock 

 Copyright (C) 1989-2007,  Peter Reilly, Garrett M. Morris, David S. Goodsell, Ruth Huey, Arthur J. Olson, 
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

#include "autocomm.h"

void mdist();

enum {C,N,O,H,XX,P,S};  // see "bond_index" in the "AD4_parameters.dat" file.
//    0 1 2 3 4  5 6
#define NUM_ENUM_ATOMTYPES 7 // this should be the length of the enumerated atom types above
	
double mindist[NUM_ENUM_ATOMTYPES][NUM_ENUM_ATOMTYPES];
double maxdist[NUM_ENUM_ATOMTYPES][NUM_ENUM_ATOMTYPES];

void mdist() {

	register int i,j;

    // Zero all the mindist and maxdist elements.
	for (i=0; i<   NUM_ENUM_ATOMTYPES; i++) {
		for (j=0; j<   NUM_ENUM_ATOMTYPES; j++) {
			mindist[i][j] = 0.0L;
			maxdist[i][j] = 0.0L;
		}
	}

    // Set all the mindist and maxdist elements to the defaults for AutoDock versions 1 - 3...
	for (i=0; i<   NUM_ENUM_ATOMTYPES; i++) {
		for (j=0; j<   NUM_ENUM_ATOMTYPES; j++) {
			mindist[i][j] = 0.9;
			maxdist[i][j] = 2.1;
		}
	}

    /*
     * These values were contributed by Peter Reilly et al.:
     *
    mindist[C][C] = 1.32;
    maxdist[C][C] = 1.545;
    mindist[C][N] = 1.32;
    maxdist[C][N] = 1.39;
    mindist[C][O] = 1.20;
    maxdist[C][O] = 1.43;
    mindist[C][XX] = 1.07;
    maxdist[C][XX] = 1.15;
    mindist[C][S] = 1.80;
    maxdist[C][S] = 1.84;
    mindist[N][H] = 0.99;
    maxdist[N][H] = 1.03;
    mindist[N][XX] = 0.99;
    maxdist[N][XX] = 1.03;
    mindist[O][H] = 0.94;
    maxdist[O][H] = 0.98;
    mindist[O][P] = 1.47;
    maxdist[O][P] = 1.63;
    mindist[H][S] = 1.316;
    maxdist[H][S] = 1.356;
    mindist[XX][S] = 1.316;
    maxdist[XX][S] = 1.356;
    mindist[S][S] = 2.018;
    maxdist[S][S] = 2.058;
     */
    mindist[C][H] = 1.07;
    maxdist[C][H] = 1.15;
    mindist[H][C] = mindist[C][H];
    maxdist[H][C] = maxdist[C][H];

    mindist[N][H] = 0.99;
    maxdist[N][H] = 1.10;
    mindist[H][N] = mindist[N][H];
    maxdist[H][N] = maxdist[N][H];

    mindist[O][H] = 0.94;
    maxdist[O][H] = 1.10;
    mindist[H][O] = mindist[O][H];
    maxdist[H][O] = maxdist[O][H];

    mindist[S][H] = 1.316;
    maxdist[S][H] = 1.356;
    mindist[H][S] = mindist[S][H];
    maxdist[H][S] = maxdist[S][H];

    mindist[P][H] = 1.35;
    maxdist[P][H] = 1.40;
    mindist[H][P] = mindist[P][H];
    maxdist[H][P] = maxdist[P][H];

    mindist[N][O] = 1.11;  // N=O is ~ 1.21 Å, minus 0.1Å error
    maxdist[N][O] = 1.50;  // N-O is ~ 1.40 Å, plus 0.1 Å error
    mindist[O][N] = mindist[N][O];  // N=O is ~ 1.21 Å, minus 0.1Å error
    maxdist[O][N] = maxdist[N][O];  // N-O is ~ 1.40 Å, plus 0.1 Å error
}; 
