/*

 $Id: input_state.cc,v 1.6.2.1 2010/11/19 20:09:29 rhuey Exp $

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
    int i, istep, status = 0;
    Real energy, eint;
    char lastmove;
    char myline[LINELEN];

#ifdef DEBUG
    fprintf(stderr, "line=|%s|\n", line);
#endif /* DEBUG */
	// Handle multi-ligand -Huameng 11/08/2007
    for(i = 0; i < S->nlig; i++) {  
	    status = sscanf(line, "%*s %d %1s " FDFMT " " FDFMT " %lf %lf %lf %lf %lf %lf %lf", &istep, &lastmove, &energy, &eint,  &(S->T[i].x), &(S->T[i].y), &(S->T[i].z),  &(S->Q[i].nx), &(S->Q[i].ny), &(S->Q[i].nz),  &(S->Q[i].ang) );	
	    
	    if (status != 0) {
		S->Q[i].ang = DegreesToRadians( S->Q[i].ang );
		mkUnitQuat( &(S->Q[i]) );
    	}
    } // for nlig
    
    *p_istep = istep;
	*p_energy = energy;
	*p_eint = eint;
	*p_lastmove = lastmove;

    for (i=0; i<ntor; i++) {
    	(void) fgets(myline, LINELEN, fp);
        sscanf(myline, "%lf", &(S->tor[i]) ); /* input torsions are in degrees */
        S->tor[i] = DegreesToRadians( S->tor[i] );    /* now in radians */
    }
    
    return( status );
}
/* EOF */
