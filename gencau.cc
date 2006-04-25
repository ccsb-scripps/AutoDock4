/*

 $Id: gencau.cc,v 1.4 2006/04/25 22:32:11 garrett Exp $

*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include <math.h>
#include "ranlib.h"


Real rcauchy(Real alpha, Real beta)
/*
**********************************************************************
     Real rcauchy(Real alpha, Real beta)
         GENerate random deviate from a CAUchy distribution
                              Function
     Generates a single random deviate from a Cauchy distribution
     with mean, alpha, and parameter, beta.
                              Arguments
**********************************************************************
*/
{
static Real gencau;

    gencau = beta*scauchy1()+alpha;
    return gencau;
}

