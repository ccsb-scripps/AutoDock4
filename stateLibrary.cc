/*

 $Id: stateLibrary.cc,v 1.8 2005/10/06 22:52:49 lindy Exp $

*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/* stateLibrary.cc */

#include <math.h>

#include <stdio.h>
#include "stateLibrary.h"

extern FILE *logFile;

void initialiseState( State *S )
{
    register int i;
    S->T.x = 0.0;
    S->T.y = 0.0;
    S->T.z = 0.0;
    S->Q.nx = 1.0;
    S->Q.ny = 0.0;
    S->Q.nz = 0.0;
    S->Q.ang = 0.0;
    S->Q.x = 1.0;
    S->Q.y = 0.0;
    S->Q.z = 0.0;
    S->Q.w = 0.0;
    S->Q.qmag = 1.0;
    S->ntor = 0;
    for (i = 0; i  < MAX_TORS;  i++ ) {
        S->tor[i] = 0.0;
    }
}

void copyState( State *D,  /* Destination -- copy to here */
                    State  S ) /* Source      -- copy this.   */
{
    register int i;
        
    D->T.x    = S.T.x;
    D->T.y    = S.T.y;
    D->T.z    = S.T.z;
    
    D->Q.nx   = S.Q.nx;
    D->Q.ny   = S.Q.ny;
    D->Q.nz   = S.Q.nz;
    D->Q.ang  = S.Q.ang;
    D->Q.x    = S.Q.x;
    D->Q.y    = S.Q.y;
    D->Q.z    = S.Q.z;
    D->Q.w    = S.Q.w;
    D->Q.qmag = S.Q.qmag;
 
    D->ntor   = S.ntor;
 
    for ( i=0; i < S.ntor; i++ ) {
            D->tor[i] = S.tor[i];
    }
}

void printState( FILE *fp, 
                 State S, 
                 int detail )
{
    register int i;
    FloatOrDouble torDegTmp;

    switch( detail ) {
        case 0:
            writeState(fp,S);
            break;

        case 2:
        default:
            (void)fprintf( fp, "\nSTATE VARIABLES:\n________________\n\n" );
            (void)fprintf( fp, "Translation x,y,z         = %.3f %.3f %.3f\n", S.T.x, S.T.y, S.T.z );
            S.Q.ang = WrpRad( ModRad( S.Q.ang ));
            (void)fprintf( fp, "Quaternion nx,ny,nz,angle = %.3f %.3f %.3f %.3f\n", S.Q.nx, S.Q.ny, S.Q.nz, Deg(S.Q.ang) );
            (void)fprintf( fp, "Quaternion x,y,z,w        = %.3f %.3f %.3f %.3f\n", S.Q.x, S.Q.y, S.Q.z, S.Q.w );
            //(void)fprintf( fp, "Quaternion qmag           = %.3f\n", S.Q.qmag );
            (void)fprintf( fp, "Number of Torsions        = %d\n", S.ntor );
            if (S.ntor > 0) {
                (void)fprintf( fp, "Torsions (degrees)        =");
                for (i=0; i<S.ntor; i++) {
                    S.tor[i] = WrpRad( ModRad( S.tor[i] ) );
                }
                for (i=0; i<S.ntor; i++) {
                    torDegTmp = Deg( S.tor[i] );
                    torDegTmp = ModDeg( torDegTmp );
                    torDegTmp = WrpDeg( torDegTmp );
                    // Commented out next line to make format more consistent, now all
                    // numbers are space-delimited.
                    //pr( fp, " %.2f%c", torDegTmp, (i==(S.ntor-1) ? '.' : ',')); 
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
    FloatOrDouble torDegTmp;

    //    (void)fprintf( fp, "State= " );

    // Write translation.
    (void)fprintf( fp, "%.3f %.3f %.3f  ", S.T.x, S.T.y, S.T.z );

    // Write quaternion.
    S.Q.ang = WrpRad( ModRad( S.Q.ang ));
    (void)fprintf( fp, "%.3f %.3f %.3f %.3f  ", S.Q.nx, S.Q.ny, S.Q.nz,
		   Deg(S.Q.ang) );
    
    // Write torsion angles.
    if (S.ntor > 0) {
        for (i=0; i<S.ntor; i++) {
            S.tor[i] = WrpRad( ModRad( S.tor[i] ) );
        }
        for (i=0; i<S.ntor; i++) {
            torDegTmp = Deg( S.tor[i] );
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

int checkState(State *D)
{
    register int i;
    int retval = 1;
        
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
    if (ISNAN(D->Q.nx)) {
        (void)fprintf(logFile,"checkState: (NaN) detected in nx quaternion\n");
        retval = 0;
    }
    if (ISNAN(D->Q.ny)) {
        (void)fprintf(logFile,"checkState: (NaN) detected in ny quaternion\n");
        retval = 0;
    }
    if (ISNAN(D->Q.nz)) {
        (void)fprintf(logFile,"checkState: (NaN) detected in nz quaternion\n");
        retval = 0;
    }
    if (ISNAN(D->Q.ang)) {
        (void)fprintf(logFile,"checkState: (NaN) detected in quaternion angle\n");
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
    D->Q.qmag = hypotenuse4(D->Q.x,  D->Q.y,  D->Q.z,  D->Q.w);
    if (ISNAN(D->Q.qmag)) {
        (void)fprintf(logFile,"checkState: (NaN) detected in magnitude of quaternion\n");
        retval = 0;
    }
 
    for ( i=0; i < D->ntor; i++ ) {
            if (ISNAN(D->tor[i])) {
                (void)fprintf(logFile,"checkState: (NaN) detected in torsion %d\n",i+1);
                retval = 0;
            }
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
