/*

 $Id: gencau.cc,v 1.6 2006/11/03 03:37:19 garrett Exp $

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

    // scauchy2 is faster than scauchy1
    // according to testing involving looping over the scauchy variants and timing each
    gencau = beta*scauchy2()+alpha;
    return gencau;
}

