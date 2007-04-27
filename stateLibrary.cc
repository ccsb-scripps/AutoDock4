/*

 $Id: stateLibrary.cc,v 1.14 2007/04/27 06:01:51 garrett Exp $

 AutoDock 

 Copyright (C) 1989-2007,  Garrett M. Morris, David S. Goodsell, Ruth Huey, Arthur J. Olson, 
 All Rights Reserved.

 AutoDock is a Trade Mark of The Scripps Research Institute.

 This program is free software; you can redistribute it and/or
 modify it under the terms of the GNU General Public License
 as published by the Free Software Foundation; either version 2
 of the License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>
#include <stdio.h>
#include "stateLibrary.h"
#include "qmultiply.h"

extern FILE *logFile;

void initialiseState( State *S )
{
    register int i;
    S->T.x = 0.0;
    S->T.y = 0.0;
    S->T.z = 0.0;
    initialiseQuat( &(S->Q) );
    S->ntor = 0;
    for (i = 0; i  < MAX_TORS;  i++ ) {
        S->tor[i] = 0.0;
    }
}

void initialiseQuat( Quat *Q )
{
    Q->nx = 1.0;
    Q->ny = 0.0;
    Q->nz = 0.0;
    Q->ang = 0.0;
    Q->x = 0.0;
    Q->y = 0.0;
    Q->z = 0.0;
    Q->w = 1.0;
    Q->qmag = 1.0;
}

void copyState( State *D,  /* Destination -- copy to here */
                    State  S ) /* Source      -- copy this.   */
{
    register int i;
        
    /*
    D->T.x    = S.T.x;
    D->T.y    = S.T.y;
    D->T.z    = S.T.z;
    */
    D->T = S.T;
    
    /*
    D->Q.nx   = S.Q.nx;
    D->Q.ny   = S.Q.ny;
    D->Q.nz   = S.Q.nz;
    D->Q.ang  = S.Q.ang;
    D->Q.x    = S.Q.x;
    D->Q.y    = S.Q.y;
    D->Q.z    = S.Q.z;
    D->Q.w    = S.Q.w;
    D->Q.qmag = S.Q.qmag;
    */
    D->Q = S.Q;
 
    D->ntor   = S.ntor;
 
    for ( i=0; i < S.ntor; i++ ) {
            D->tor[i] = S.tor[i];
    }

    D->hasEnergy = S.hasEnergy;

    /*
    D->e.total = S.e.total;
    D->e.intra = S.e.intra;
    D->e.inter = S.e.inter;
    D->e.FE = S.e.FE;
    */
    D->e = S.e;
}

void printState( FILE *fp, 
                 State S, 
                 int detail )
{
    register int i;
    Real torDegTmp;

    switch( detail ) {
        case 0:
        case 1:
            writeState(fp,S);
            break;

        case 2:
        default:
            (void)fprintf( fp, "\nSTATE VARIABLES:\n________________\n\n" );
            (void)fprintf( fp, "Translation x,y,z         = %.3f %.3f %.3f\n", S.T.x, S.T.y, S.T.z );
            (void)fprintf( fp, "Quaternion x,y,z,w        = %.3f %.3f %.3f %.3f\n", S.Q.x, S.Q.y, S.Q.z, S.Q.w );
            S.Q = convertQuatToRot( S.Q );
            S.Q.ang = WrpRad( ModRad( S.Q.ang ));
            (void)fprintf( fp, "Axis-Angle nx,ny,nz,angle = %.3f %.3f %.3f %.3f\n", S.Q.nx, S.Q.ny, S.Q.nz, RadiansToDegrees(S.Q.ang) );
            //(void)fprintf( fp, "Quaternion qmag           = %.3f\n", S.Q.qmag );
            (void)fprintf( fp, "Number of Torsions        = %d\n", S.ntor );
            if (S.ntor > 0) {
                (void)fprintf( fp, "Torsions (degrees)        =");
                for (i=0; i<S.ntor; i++) {
                    S.tor[i] = WrpRad( ModRad( S.tor[i] ) );
                }
                for (i=0; i<S.ntor; i++) {
                    torDegTmp = RadiansToDegrees( S.tor[i] );
                    torDegTmp = ModDeg( torDegTmp );
                    torDegTmp = WrpDeg( torDegTmp );
                    pr( fp, " %.2f", torDegTmp );
                    //if ((B_isTorConstrained[i] == 1) && B_ShowTorE) {
                        //pr( fp, ", Energetic penalty = %uhd\n", US_TorE[i]);
                    //} else {
                        //pr( fp, "\n");
                    //}
                }
            }
            (void)fprintf( fp, "\n\n");
            break;

        case 3:
            // Writes only the translation component of the state
            (void)fprintf( fp, "%.3f %.3f %.3f", S.T.x, S.T.y, S.T.z );
            break;
    }
}

void writeState( FILE *fp, State S )
{
    register int i;
    Real torDegTmp;

    //    (void)fprintf( fp, "State= " );

    // Write translation.
    (void)fprintf( fp, "%.3f %.3f %.3f  ", S.T.x, S.T.y, S.T.z );

    // Convert quaternion to axis-angle.
    S.Q = convertQuatToRot( S.Q );

    // Write axis-angle.
    S.Q.ang = WrpRad( ModRad( S.Q.ang ));
    (void)fprintf( fp, "%.3f %.3f %.3f %.3f  ", S.Q.nx, S.Q.ny, S.Q.nz, RadiansToDegrees(S.Q.ang) );
    
    // Write torsion angles.
    if (S.ntor > 0) {
        for (i=0; i<S.ntor; i++) {
            S.tor[i] = WrpRad( ModRad( S.tor[i] ) );
        }
        for (i=0; i<S.ntor; i++) {
            torDegTmp = RadiansToDegrees( S.tor[i] );
            torDegTmp = ModDeg( torDegTmp );
            torDegTmp = WrpDeg( torDegTmp );
            // Commented out next line to make format more consistent, now all
            // numbers are space-delimited.
            //pr( fp, " %.2f%c", torDegTmp, (i==(S.ntor-1) ? '.' : ','));
            pr( fp, " %.2f", torDegTmp );
        }
    }
    // Leave fp on this line for energies which follow....
    //    (void)fprintf( fp, "\n");
}

int checkState( const State *D )
{
    register int i;
    int retval = 1;
    double magnitude_q = 0.;
        
    if (ISNAN(D->T.x)) {
        (void)fprintf(logFile,"checkState: (NaN) detected in x translation\n");
        retval = 0;
    }
    if (ISNAN(D->T.y)) {
        (void)fprintf(logFile,"checkState: (NaN) detected in y translation\n");
        retval = 0;
    }
    if (ISNAN(D->T.z)) {
        (void)fprintf(logFile,"checkState: (NaN) detected in z translation\n");
        retval = 0;
    
    }

    if (ISNAN(D->Q.x)) {
        (void)fprintf(logFile,"checkState: (NaN) detected in x quaternion\n");
        retval = 0;
    }
    if (ISNAN(D->Q.y)) {
        (void)fprintf(logFile,"checkState: (NaN) detected in y quaternion\n");
        retval = 0;
    }
    if (ISNAN(D->Q.z)) {
        (void)fprintf(logFile,"checkState: (NaN) detected in z quaternion\n");
        retval = 0;
    }
    if (ISNAN(D->Q.w)) {
        (void)fprintf(logFile,"checkState: (NaN) detected in w quaternion\n");
        retval = 0;
    }
    magnitude_q = hypotenuse4(D->Q.x,  D->Q.y,  D->Q.z,  D->Q.w);
    if (ISNAN(magnitude_q)) {
        (void)fprintf(logFile,"checkState: After computing the magnitude of quaternion, (NaN) was detected\n");
        retval = 0;
    }
 
    for ( i=0; i < D->ntor; i++ ) {
            if (ISNAN(D->tor[i])) {
                (void)fprintf(logFile,"checkState: (NaN) detected in torsion %d\n",i+1);
                retval = 0;
            }
    }

    if (ISNAN(D->Q.nx)) {
        (void)fprintf(logFile,"checkState: (NaN) detected in nx component of axis-angle\n");
        retval = 0;
    }
    if (ISNAN(D->Q.ny)) {
        (void)fprintf(logFile,"checkState: (NaN) detected in ny component of axis-angle\n");
        retval = 0;
    }
    if (ISNAN(D->Q.nz)) {
        (void)fprintf(logFile,"checkState: (NaN) detected in nz component of axis-angle\n");
        retval = 0;
    }
    if (ISNAN(D->Q.ang)) {
        (void)fprintf(logFile,"checkState: (NaN) detected in ang component of axis-angle\n");
        retval = 0;
    }


    return(retval);
}

Molecule copyStateToMolecule(State *S, Molecule *mol) /* S is the source */
{
    register int i;
    mol->S.T.x = S->T.x;
    mol->S.T.y = S->T.y;
    mol->S.T.z = S->T.z;
    mol->S.Q.nx = S->Q.nx;
    mol->S.Q.ny = S->Q.ny;
    mol->S.Q.nz = S->Q.nz;
    mol->S.Q.ang = S->Q.ang;
    mol->S.Q.x = S->Q.x;
    mol->S.Q.y = S->Q.y;
    mol->S.Q.z = S->Q.z;
    mol->S.Q.w = S->Q.w;
    mol->S.Q.qmag = S->Q.qmag;
    mol->S.ntor = S->ntor;
    for (i = 0; i  < MAX_TORS;  i++ ) {
        mol->S.tor[i] = S->tor[i];
    }
    return *mol;
}
/* EOF */
