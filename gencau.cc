/*

 $Id: gencau.cc,v 1.3 2003/02/27 02:11:08 lindy Exp $

*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


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

