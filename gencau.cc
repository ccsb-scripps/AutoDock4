/***** RCS INFO **********************************************************

$Id: gencau.cc,v 1.1.1.1 2001/08/13 22:05:53 gillet Exp $
$Source: /Users/mp/facil/autodock/git-luna/autodock-cvstar/gencau.cc,v $
$Log: gencau.cc,v $
Revision 1.1.1.1  2001/08/13 22:05:53  gillet
 import initial of autodock sources

// Revision 3.0  1996/03/11  05:40:00  halliday
// The function definition for the GA/LS hybrid.
//

****** RCS INFO *********************************************************/

/* gencau.c
 *
 */

#include <math.h>
#include "ranlib.h"


float rcauchy(float alpha, float beta)
/*
**********************************************************************
     float rcauchy(float alpha, float beta)
         GENerate random deviate from a CAUchy distribution
                              Function
     Generates a single random deviate from a Cauchy distribution
     with mean, alpha, and parameter, beta.
                              Arguments
**********************************************************************
*/
{
static float gencau;

    gencau = beta*scauchy1()+alpha;
    return gencau;
}

