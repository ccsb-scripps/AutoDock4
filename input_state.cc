/*

 $Id: input_state.cc,v 1.7 2009/05/08 23:02:13 rhuey Exp $

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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>
#include <stdio.h>
#include <string.h>
#include "input_state.h"

#define LINELEN 132

int input_state( State *S,
		 FILE  *fp,
		 char  line[LINE_LEN],
		 int   ntor,
		 int   *p_istep,
		 Real *p_energy,
		 Real *p_eint,
		 char  *p_lastmove )
{
    int i, istep, status;
    Real energy, eint;
    char lastmove;
    char myline[LINELEN];

#ifdef DEBUG
    fprintf(stderr, "line=|%s|\n", line);
#endif /* DEBUG */

    status = sscanf(line, "%*s %d %1s " FDFMT " " FDFMT " %lf %lf %lf %lf %lf %lf %lf", &istep, &lastmove, &energy, &eint,  &(S->T.x), &(S->T.y), &(S->T.z),  &(S->Q.nx), &(S->Q.ny), &(S->Q.nz),  &(S->Q.ang) );

    if (status != 0) {
	S->Q.ang = DegreesToRadians( S->Q.ang );
	mkUnitQuat( &(S->Q) );

    *p_istep = istep;
	*p_energy = energy;
	*p_eint = eint;
	*p_lastmove = lastmove;

        for (i=0; i<ntor; i++) {
	    (void) fgets(myline, LINELEN, fp);
            sscanf(myline, "%lf", &(S->tor[i]) ); /* input torsions are in degrees */
            S->tor[i] = DegreesToRadians( S->tor[i] );    /* now in radians */
        }
    }
    return( status );
}
/* EOF */
