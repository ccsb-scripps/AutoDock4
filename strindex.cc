/*

 $Id: strindex.cc,v 1.3 2005/09/28 22:54:21 garrett Exp $

*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif



#include "strindex.h"


int strindex( char s[], char t[] )

/******************************************************************************/
/*      Name: strindex                                                        */
/*  Function: return index of t in s, -1 if none.                             */
/* Copyright: (C) Garrett Morris, TSRI                                        */
/*----------------------------------------------------------------------------*/
/*    Author: Garrett Morris, The Scripps Research Institute                  */
/*      Date: 11/30/92                                                        */
/*----------------------------------------------------------------------------*/
/*    Inputs: s,t                                                             */
/*   Returns: i, index of t in s, -1 if not found.                            */
/*   Globals: none.                                                           */
/*----------------------------------------------------------------------------*/
/* Modification Record                                                        */
/* Date     Inits   Comments                                                  */
/* 11/30/92 GMM     Entered code.                                             */
/******************************************************************************/

{
    register int i,c1,c2;

    for (i=0; s[i] != '\0'; i++) {
        for (c1=i, c2=0; s[c1]==t[c2] && t[c2]!='\0'; c1++, c2++)
            ;
        if (c2>0 && t[c2]=='\0')
            return( i );
    }
    return( -1 );
}
/* EOF */
