/*

 $Id: parse_com_line.cc,v 1.2 2003/02/26 01:24:34 garrett Exp $

*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/* parse_com_line.cc */

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
