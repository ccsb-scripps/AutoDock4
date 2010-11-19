/*

 $Id: getrms.cc,v 1.5 2007/04/27 06:01:48 garrett Exp $

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

#include <math.h>
#include "getrms.h"


Real getrms ( Real Crd[MAX_ATOMS][SPACE], 
               Real CrdRef[MAX_ATOMS][SPACE], 
               Boole B_symmetry_flag, 
               int natom, 
               int type[MAX_ATOMS] )

{
    double sqrSum, sqrMin, dc[SPACE];
    register int i, j, xyz;

    sqrSum = 0.;

    if (B_symmetry_flag) {
        for (i = 0;  i < natom;  i++) {
            sqrMin = BIG;
            for (j = 0;  j < natom;  j++) {                
                if (type[i] == type[j]) {
                    for (xyz = 0;  xyz < SPACE;  xyz++) {
                        dc[xyz]= Crd[i][xyz] - CrdRef[j][xyz];
                    } /* xyz */
                    sqrMin = min( sqhypotenuse(dc[X], dc[Y], dc[Z]), sqrMin );
                }
            } /*  next j  */
            sqrSum += sqrMin;
        } /*  next i  */
    } else {
        for (i = 0;  i < natom;  i++) {
            for (xyz = 0;  xyz < SPACE;  xyz++) {
                dc[xyz]= Crd[i][xyz] - CrdRef[i][xyz];
            } /* xyz */
            sqrSum += sqhypotenuse( dc[X], dc[Y], dc[Z] );
        } /*  next i  */
    }

    return ( sqrt( sqrSum / (double)natom )  );
}
/* EOF */
