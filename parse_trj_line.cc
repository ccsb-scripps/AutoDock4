/* parse_trj_line.cc */

    #include <stdio.h>
    #include <string.h>
    #include <ctype.h>
    #include "parse_trj_line.h"
    #include "trjtokens.h"


int parse_trj_line( char line[LINE_LEN] )

/******************************************************************************/
/*      Name: parse_trj_line                                                  */
/*  Function: Parse the trajectory file line                                  */
/* Copyright: (C) 1994, TSRI                                                  */
/*----------------------------------------------------------------------------*/
/*    Author: Garrett Morris, The Scripps Research Institute                  */
/*      Date: 26/01/93                                                        */
/*----------------------------------------------------------------------------*/
/*    Inputs: line                                                            */
/*   Returns: integer token describing the keyword found.                     */
/*   Globals: none.                                                           */
/*----------------------------------------------------------------------------*/
/* Modification Record                                                        */
/* Date     Inits   Comments                                                  */
/* 26/01/93 GMM     Entered code.                                             */
/******************************************************************************/

{
    int l, i, token = -1 ;	       /* return -1 if nothing is recognized. */
    char c[LINE_LEN];

    l = strindex(line, " ");

    if (l == -1)
        l = strlen(line);
    
    for (i=0; i<l; i++)
        c[i] = (char)tolower( (int)line[i] );

    if ((c[0]=='\n')||(c[0]=='\0')) {
        token = TRJ_NULL;
    } else if (strncmp(c,"ntorsions",9)==0) {
        token = TRJ_NTOR;
    } else if (strncmp(c,"run",3)==0) {
        token = TRJ_RUN;
    } else if (strncmp(c,"cycle",5)==0) {
        token = TRJ_CYCLE;
    } else if (strncmp(c,"state",5)==0) {
        token = TRJ_STATE;
    } else if (strncmp(c,"temp",4)==0) {
        token = TRJ_TEMP;
    }

    return(token);
}
/* EOF */
