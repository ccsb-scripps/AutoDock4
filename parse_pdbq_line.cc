/* parse_pdbq_line.cc */


    #include <stdio.h>
    #include <string.h>
    #include <ctype.h>
    #include "parse_pdbq_line.h"
    #include "pdbqtokens.h"


extern FILE *logFile;


int parse_pdbq_line( char line[LINE_LEN] )

/******************************************************************************/
/*      Name: parse_pdbq_line                                                 */
/*  Function: Parse the PDBQ file line.                                       */
/* Copyright: (C) 1994, TSRI                                                  */
/*----------------------------------------------------------------------------*/
/*    Author: Garrett Morris, The Scripps Research Institute                  */
/*      Date: 11/06/93                                                        */
/*----------------------------------------------------------------------------*/
/*    Inputs: line                                                            */
/*   Returns: integer token describing the keyword found.                     */
/*   Globals: none.                                                           */
/*----------------------------------------------------------------------------*/
/* Modification Record                                                        */
/* Date     Inits   Comments                                                  */
/* 11/06/93 GMM     Entered code.                                             */
/******************************************************************************/

{
    int l, i, token = PDBQ_UNRECOGNIZED;
    char c[LINE_LEN];

    if ((l = strindex(line, " ")) == -1) {
        l = strlen(line);
    }
    for (i=0; i<l && i<(LINE_LEN-1); i++) {
        c[i] = (char)tolower( (int)line[i] );
    }
    c[i] = '\0';

    if ((c[0]=='\n')||(c[0]=='\0')) {
	token = PDBQ_NULL;
    } else if (l >= 4) {
	if ((strncmp(c,"rema",4)==0) || (strncmp(c,"user",4)==0)) {
	    token = PDBQ_REMARK;
	} else if (strncmp(c,"root",4)==0) {
	    token = PDBQ_ROOT;
	} else if (strncmp(c,"endr",4)==0) {
	    token = PDBQ_ENDROOT;
	} else if (strncmp(c,"atom",4)==0) {
	    token = PDBQ_ATOM;
	} else if (strncmp(c,"heta",4)==0) {
	    token = PDBQ_HETATM;
	} else if ((strncmp(c,"tors",4)==0) && (strncmp(c,"torsdof",7)!=0)) {
	    token = PDBQ_TORS;
	} else if ((strlen(c) >= 7) && (strncmp(c,"torsdof",7)==0)) {
	    token = PDBQ_TORSDOF;
	} else if (strncmp(c,"tdof",4)==0) {
	    token = PDBQ_TORSDOF;
	} else if (strncmp(c,"endt",4)==0) {
	    token = PDBQ_ENDTORS;
	} else if (strncmp(c,"bran",4)==0) {
	    token = PDBQ_BRANCH;
	} else if (strncmp(c,"endb",4)==0) {
	    token = PDBQ_ENDBRANCH;
	} else if (strncmp(c,"cons",4)==0) {
	    token = PDBQ_CONSTRAINT;
	}
    } else {
	token = PDBQ_NULL;
    }

    return( token );
}
/* EOF */
