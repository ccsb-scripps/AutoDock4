/* qmultiply.cc */

#include <math.h>

    #include <stdio.h>
    #include <string.h>
    #include "qmultiply.h"


void qmultiply( Quat *q,
                Quat *ql,
                Quat *qr )

/******************************************************************************/
/*      Name: qultiply                                                        */
/*  Function: Quaternion Multiplication (Accelerated)                         */
/*            [q]  =  [ql] [qr]                                               */
/*            [s1,v1][s2,v2] = [(s1*s2 - v1.v2), (s1*v2 + s2*v1 + v1^v2)]     */
/*                ~~     ~~              ~~ ~~       ~~      ~~   ~~ ~~       */
/*            Quaternion vector q[3,4,5] need not be normalized outside this  */
/*            routine.                                                        */
/*            q[] contains (xt,yt,zt, x,y,z,w(radians))                       */
/*            if ql[QW] is zero, qr is new vector                             */
/*            if qr[QW] is zero, ql is new vector                             */
/* Copyright: (C) 1994, TSRI                                                  */
/*----------------------------------------------------------------------------*/
/*   Authors: Garrett M. Morris, The Scripps Research Institute.              */
/*            David Goodsell, TSRI                                            */
/*      Date: 12/03/92                                                        */
/*----------------------------------------------------------------------------*/
/*    Inputs: ql = last quaternion                                            */
/*            qr = (random) rotation to be applied to ql                      */
/*   Returns: q  = resultant quaternion                                       */
/*   Globals: none.                                                           */
/*----------------------------------------------------------------------------*/
/* Modification Record                                                        */
/* Date     Inits   Comments                                                  */
/* 05/15/92 GMM     Translated into C                                         */
/* 12/03/92 GMM     Changed '/2.' to '*0.5'; introduced 'hqwl' and 'hqwr'.    */
/* 12/03/92 GMM     Replaced rqtot by inv_qag; was '/rqtot', now '*inv_qag'   */
/******************************************************************************/
{ 
    if (ql->w == 0.) {
        q->x = qr->x;
        q->y = qr->y;
        q->z = qr->z;
        q->w = qr->w;
        return;
    }

    if (qr->w == 0.) {
        q->x = ql->x;
        q->y = ql->y;
        q->z = ql->z;
        q->w = ql->w;
        return;
    }
    
    q->x = (double) (ql->w*qr->x + ql->x*qr->w + ql->y*qr->z - ql->z*qr->y);
    q->y = (double) (ql->w*qr->y + ql->y*qr->w + ql->z*qr->x - ql->x*qr->z);
    q->z = (double) (ql->w*qr->z + ql->z*qr->w + ql->x*qr->y - ql->y*qr->x);
    q->w = (double) (ql->w*qr->w - ql->x*qr->x - ql->y*qr->y - ql->z*qr->z);

    // q->qmag  = hypotenuse4( q->x,  q->y,  q->z,  q->w  );

/* make sure you put in the conversion from x to nx, w to ang */
}

void mkUnitQuat( Quat *q )
{	
    double inv_nmag, hqang, s;
	     
    inv_nmag = 1. / hypotenuse( q->nx, q->ny, q->nz );
    q->nx *= inv_nmag;       /* Normalize q */
    q->ny *= inv_nmag;       /* Normalize q */
    q->nz *= inv_nmag;       /* Normalize q */
      
    hqang = 0.5 * q->ang;
    s     = sin( hqang );
    
    q->w  = cos( hqang );
    q->x  = s * q->nx;
    q->y  = s * q->ny;
    q->z  = s * q->nz;
    
    /* q->qmag = hypotenuse4( q->x,  q->y,  q->z,  q->w  ); */
}

void printQuat( Quat q )
{
    printf("Quat(x,y,z,w)= %.2f %.2f %.2f %.2f\n", q.x, q.y, q.z, q.w);
}

/* EOF */
