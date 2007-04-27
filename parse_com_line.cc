/*

 $Id: parse_com_line.cc,v 1.3 2007/04/27 06:01:50 garrett Exp $

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
#include <ctype.h>
#include "parse_com_line.h"
#include "cmdtokens.h"


int parse_com_line( char line[LINE_LEN] )

/******************************************************************************/
/*      Name: parse_com_line                                                  */
/*  Function: Parse the command line                                          */
/* Copyright: (C) 1994, TSRI                                                  */
/*----------------------------------------------------------------------------*/
/*    Author: Garrett Morris, The Scripps Research Institute                  */
/*      Date: 16/01/93                                                        */
/*----------------------------------------------------------------------------*/
/*    Inputs: line                                                            */
/*   Returns: integer token describing the command found.                     */
/*   Globals: none.                                                           */
/*----------------------------------------------------------------------------*/
/* Modification Record                                                        */
/* Date     Inits   Comments                                                  */
/* 06/09/95 RSH     Fixed bug                                                 */
/* 16/01/93 GMM     Entered code.                                             */
/******************************************************************************/

{
    int i, token = COM_NULL ;
    char c[4];

    for (i=0; i<4; i++)
        c[i] = (char)tolower( (int)line[i] );

    if ((c[0]=='\n')||(c[0]=='\0')) {
        token = COM_NULL;
    } else if ( (strncmp(c,"stop",4)==0) ||
                (strncmp(c,"exit",4)==0) ||
                (strncmp(c,"quit",4)==0)    ) {
        token = COM_STOP;
    } else if (strncmp(c,"eval",4)==0) {
        token = COM_EVAL;
    } else if (strncmp(c,"outc",4)==0) {
        token = COM_OUTC;
    } else if (strncmp(c,"oute",4)==0) {
        token = COM_OUTE;
    } else if (strncmp(c,"traj",4)==0) {
        token = COM_TRJ;
    } else if ( (strncmp(c,"epdb",4)==0) ||
                (strncmp(c,"ener",4)==0)    ) {
        token = COM_EPDB;
    }

    return(token);
}
/* EOF */
