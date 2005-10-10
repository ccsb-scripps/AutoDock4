/*

 $Id: qtransform.cc,v 1.3.8.1 2005/10/10 16:39:32 alther Exp $

*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/* qtransform.cc */

#include "qtransform.h"
#include <stdio.h>
#include <string.h>


void qtransform( const Coord T,
                 const Quat  q,
                 FloatOrDouble tcoord[MAX_ATOMS][SPACE],
                 const int   natom)

/******************************************************************************/
/*      Name: qtransform                                                      */
/*  Function: Accelerated quaternion transformation                           */
/*            Performs both a rigid-body translation and rotation.            */
/*            Assumes quaternion vector q.nx,q.ny,q.nz is normalized outside  */
/*            this routine.                                                   */
/* Copyright: (C) 1995, TSRI                                                  */
/*----------------------------------------------------------------------------*/
/*   Authors: Garrett M. Morris, The Scripps Research Institute               */
/*            David Goodsell, UCLA                                            */
/*      Date: 11/23/94                                                        */
/*----------------------------------------------------------------------------*/
/*    Inputs: q, tcoord, natom                                                */
/*            q[] contains (xt,yt,zt,x,y,z,w(radians))                        */
/*   Returns: tcoord                                                          */
/*   Globals: QUAT, MAX_ATOMS                                                 */
/*----------------------------------------------------------------------------*/
/* Modification Record                                                        */
/* Date     Inits   Comments                                                  */
/* 05/15/92 GMM     Translated into C.                                        */
/* 11/23/92 GMM     Introduced (15+,9*) version, replacing (12+,36*) version. */
/* 11/23/92 GMM     Introduced h and inv_qmag, since * is faster than /.      */
/******************************************************************************/
{
    register int a;
    Coord  tmp;
    double w, x, y, z;
    double tx, ty, tz;
    double omtxx;
    double twx, txy, txz;
    double twy, tyy, tyz;
    double twz, tzz;
    double r11, r12, r13, r21, r22, r23, r31, r32, r33;

    w = q.w;
    x = q.x;
    y = q.y;
    z = q.z;


/*  12 adds, 36 multiplies...     */
/*    r11 = 1. - 2.*y*y - 2.*z*z; */
/*    r12 =      2.*x*y + 2.*w*z; */
/*    r13 =      2.*x*z - 2.*w*y; */
/*    r21 =      2.*x*y - 2.*w*z; */
/*    r22 = 1. - 2.*x*x - 2.*z*z; */
/*    r23 =      2.*y*z + 2.*w*x; */
/*    r31 =      2.*x*z + 2.*w*y; */
/*    r32 =      2.*y*z - 2.*w*x; */
/*    r33 = 1. - 2.*x*x - 2.*y*y; */

/*  15 adds, 9 multiplies...      */
    tx  = x+x;
    ty  = y+y;
    tz  = z+z;

    twx = w*tx;
    omtxx = 1. - x*tx;
    txy = y*tx;
    txz = z*tx;

    twy = w*ty;
    tyy = y*ty;
    tyz = z*ty;

    twz = w*tz;
    tzz = z*tz;

    r11 = 1. - tyy - tzz;
    r12 =      txy + twz;
    r13 =      txz - twy;
    r21 =      txy - twz;
    r22 = omtxx    - tzz;
    r23 =      tyz + twx;
    r31 =      txz + twy;
    r32 =      tyz - twx;
    r33 = omtxx    - tyy;

    for (a = 0;  a < natom;  a++) {
        tmp.x = ((double)tcoord[a][X])*r11 + ((double)tcoord[a][Y])*r21 + ((double)tcoord[a][Z])*r31 + T.x;
        tmp.y = ((double)tcoord[a][X])*r12 + ((double)tcoord[a][Y])*r22 + ((double)tcoord[a][Z])*r32 + T.y;
        tmp.z = ((double)tcoord[a][X])*r13 + ((double)tcoord[a][Y])*r23 + ((double)tcoord[a][Z])*r33 + T.z;
        tcoord[a][X] = tmp.x;
        tcoord[a][Y] = tmp.y;
        tcoord[a][Z] = tmp.z;
    }
}
/* EOF */
