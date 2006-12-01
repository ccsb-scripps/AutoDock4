/*

 $Id: qmultiply.cc,v 1.6 2006/12/01 02:23:51 garrett Exp $

*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/* qmultiply.cc */

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include "qmultiply.h"

#ifdef DEBUG_MUTATION
extern  FILE    *logFile;
#endif


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
    // essentially, convertQuatToRot( Quat *q )
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
} // mkUnitQuat( Quat *q )

void printQuat( FILE *fp, Quat q )
{
    (void) fprintf( fp, "Quat(x,y,z,w)=      %5.2f %5.2f %5.2f %5.2f\n", q.x, q.y, q.z, q.w);
    (void) fprintf( fp, "Mag(Quat(x,y,z,w))= %5.2f\n", sqrt(q.x*q.x + q.y*q.y + q.z*q.z + q.w*q.w) );
    (void) fprintf( fp, "Quat(nx,ny,nz,ang)= %5.2f %5.2f %5.2f %5.2f\n", q.nx, q.ny, q.nz, q.ang);
    (void) fprintf( fp, "Mag(nx,ny,nz)=      %5.2f\n", sqrt(q.nx*q.nx + q.ny*q.ny + q.nz*q.nz) );
} // printQuat( Quat q )

Quat normQuat( Quat q )
    // Normalise the 4D quaternion, x,y,z,w
{
    double mag4 = hypotenuse4( q.x, q.y, q.z, q.w );
    if (mag4 > APPROX_ZERO) {
        double inv_mag4 = 1. / mag4;
        q.x *= inv_mag4;
        q.y *= inv_mag4;
        q.z *= inv_mag4;
        q.w *= inv_mag4;
    }
    return q;
}

Quat normRot( Quat q )
    // Normalise the 3D rotation axis or vector nx,ny,nz
{
    double mag3 = hypotenuse( q.nx, q.ny, q.nz );
    if (mag3 > APPROX_ZERO) {
        double inv_mag3 = 1. / mag3;
        q.nx *= inv_mag3;
        q.ny *= inv_mag3;
        q.nz *= inv_mag3;
    }
    return q;
}

Quat convertQuatToRot( Quat q )
    // Update the (nx,ny,nz,ang) components of the quaternion q, 
    // to correspond to the (x,y,z,w) components.
{
#ifdef DEBUG_MUTATION
    fprintf( logFile, "q.w = %.3f\n", q.w );
#endif
    assert( fabs( q.w ) <= 1.0 );
    register double angle = 2 * acos( q.w );
    register double inv_sin_half_angle = 1 / sin( angle / 2 );
    Quat retval;

    retval.nx = q.x * inv_sin_half_angle;
    retval.ny = q.y * inv_sin_half_angle;
    retval.nz = q.z * inv_sin_half_angle;

    retval = normRot( retval );

    if (angle > PI)  angle -= TWOPI;  // by convention, angles should be in the range -PI to +PI.
    retval.ang = angle;

    retval.x = q.x;
    retval.y = q.y;
    retval.z = q.z;
    retval.w = q.w;

    return retval;
} // convertQuatToRot( Quat q )

Quat convertRotToQuat( Quat q )
    // Normalize the rotation-about-axis vector 
    // and convert the rotation-about-axis components (nx,ny,nz,ang)
    // to the corresponding quaternion components (x,y,z,w)
{	
    double hqang, s;
    Quat retval;

    retval.nx = q.nx;
    retval.ny = q.ny;
    retval.nz = q.nz;
    retval = normRot( retval );

    retval.ang = q.ang;
      
    hqang = 0.5 * q.ang;
    s = sin( hqang );
    
    retval.x = s * q.nx;
    retval.y = s * q.ny;
    retval.z = s * q.nz;
    retval.w = cos( hqang );
    
    /* q.qmag = hypotenuse4( q.x,  q.y,  q.z,  q.w  ); */
    return retval;
} // Quat convertRotToQuat( Quat q )

Quat uniformQuat( void )
    // Generate a uniformly-distributed random quaternion
{
    double x0, r1, r2, t1, t2;  // for uniformly distributed quaternion calculation
    Quat q;

    /*
    **  This should produce a uniformly distributed quaternion, according to
    **  Shoemake, Graphics Gems III.6, pp.124-132, "Uniform Random Rotations",
    **  published by Academic Press, Inc., (1992)
    */
    t1 = genunf(0., TWOPI);
    // q.x = sin( t1 ) * (  r1 = ( (genunf(0., 1.) < 0.5) ?  (-1.) : (+1.) ) * sqrt( 1. - (x0 = genunf(0., 1.)) )  );  // random sign version
    q.x = sin( t1 ) * (  r1 = sqrt( 1. - (x0 = genunf(0., 1.)) )  );  // strict Shoemake version
    q.y = cos( t1 ) * r1;
    t2 = genunf(0., TWOPI);
    // q.z = sin( t2 ) * (  r2 = ( (genunf(0., 1.) < 0.5) ?  (-1.) : (+1.) ) * sqrt( x0 )  );  // random sign version
    q.z = sin( t2 ) * (  r2 = sqrt( x0 )  );  // strict Shoemake version
    q.w = cos( t2 ) * r2;

    return q;
}

void unitQuat2rotation( Quat *q )
    // Convert from a unit quaternion to a rotation about an unit 3D-vector
{
    double inv_sin_half_ang;

    q->ang = 2. * acos( q->w );
    inv_sin_half_ang = 1. / sin( 0.5 * q->ang );
    q->nx  = q->x * inv_sin_half_ang; 
    q->ny  = q->y * inv_sin_half_ang; 
    q->nz  = q->z * inv_sin_half_ang; 
    
    return;
}


/* EOF */
