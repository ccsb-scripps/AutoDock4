/*

 $Id: gencau.cc,v 1.2 2003/02/26 01:02:47 garrett Exp $

*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/***** RCS INFO **********************************************************

$Id: gencau.cc,v 1.2 2003/02/26 01:02:47 garrett Exp $
$Source: /Users/mp/facil/autodock/git-luna/autodock-cvstar/gencau.cc,v $
$Log: gencau.cc,v $
Revision 1.2  2003/02/26 01:02:47  garrett

General changes to generate identical output on each hardware/os
platform (as much as possible):

    -DUSE_DOUBLE compile switch introduced. Mixed float/double
    precision would be used otherwise, as before the "pre-integration-305"
    stage.  This switch controls the new typedef "FloatOrDouble", which is
    "double" when -DUSE_DOUBLE is used as a compile-time flag, otherwise this
    is typedef'ed to "float". (GMM, WML)

    Prior to switching all floats to doubles, the AutoDock docking
    log files often had differences in energies, coordinates or state
    variables in the sixth (or so) decimal place of floats.

Lindy added the HAVE_CONFIG_H compile-time switch with "config.h"
in preparation for "autoconf" and "configure"... (WML)

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


FloatOrDouble rcauchy(FloatOrDouble alpha, FloatOrDouble beta)
/*
**********************************************************************
     FloatOrDouble rcauchy(FloatOrDouble alpha, FloatOrDouble beta)
         GENerate random deviate from a CAUchy distribution
                              Function
     Generates a single random deviate from a Cauchy distribution
     with mean, alpha, and parameter, beta.
                              Arguments
**********************************************************************
*/
{
static FloatOrDouble gencau;

    gencau = beta*scauchy1()+alpha;
    return gencau;
}

