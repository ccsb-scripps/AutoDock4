/*

 $Id: trilinterp.cc,v 1.4 2004/11/16 23:42:54 garrett Exp $

*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/* trilinterp.cc */

#include <math.h>
#include "trilinterp.h"


/* linear interpolation from l (when a=0) to h (when a=1)*/
/* (equal to (a*h)+((1-a)*l) )*/
#define LERP(a,l,h)	((l)+(((h)-(l))*(a)))

extern int ElecMap;

#ifdef DEBUG
#include <stdio.h>
extern FILE *logFile;
#endif

FloatOrDouble trilinterp( CONST_FLOAT tcoord[MAX_ATOMS][SPACE], 
          CONST_FLOAT charge[MAX_ATOMS], 
          CONST_INT   type[MAX_ATOMS], 
          CONST_INT   total_atoms, 
          CONST_FLOAT map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS],
          CONST_FLOAT inv_spacing, 
          FloatOrDouble elec[MAX_ATOMS], 
          FloatOrDouble emap[MAX_ATOMS], 
          CONST_FLOAT xlo,	/**/
          CONST_FLOAT ylo,	/*   FloatOrDouble lo[SPACE] ) SLOWER */
          CONST_FLOAT zlo )	/**/

/*
** FloatOrDouble tcoord[MAX_ATOMS][SPACE];	temporary coordinates
** FloatOrDouble charge[MAX_ATOMS];		partial atomic charges
** int   type[MAX_ATOMS];		atom type of each atom
** int   total_atoms;			number of atoms
** FloatOrDouble map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS];  
** 					intermolecular interaction energies
** FloatOrDouble inv_spacing;  			= 1/(grid point spacing, in Angstroms) 
** FloatOrDouble elec[MAX_ATOMS];  		electrostatic energies, atom by atom
** FloatOrDouble emap[MAX_ATOMS];  		intermolecular energies
** FloatOrDouble xlo,ylo,zlo;			minimum coordinates in x,y,z
*/

/* { */

/******************************************************************************/
/*      Name: trilinterp                                                      */
/*  Function: Trilinear interpolation of interaction energies from map[]      */
/*            using the coordinates in tcoord[].                              */
/* Copyright: (C) 1994, TSRI                                                  */
/*----------------------------------------------------------------------------*/
/*   Authors: Garrett M. Morris, TSRI, Accelerated C version 2.2              */
/*            David Goodsell, UCLA, Original FORTRAN version 1.0              */
/*      Date: 10/06/94                                                        */
/*----------------------------------------------------------------------------*/
/*    Inputs: tcoord, charge, type, total_atoms, map, inv_spacing, lo         */
/*   Returns: total energy                                                    */
/*   Globals: MAX_ATOMS, SPACE, MAX_ATOMS, MAX_GRID_PTS, MAX_MAPS.            */
/*----------------------------------------------------------------------------*/
/* Modification Record                                                        */
/* Date     Inits   Comments                                                  */
/* 05/05/91 GMM     Translated into C.                                        */
/* 01/03/94 GMM     Optimized code by examining 'cc -S trilinterp.c' output.  */
/* 10/06/94 GMM     Optional 10% gain in speed, using nearest point, not      */
/*                  trilinear interpolation. Compile with -DMINPOINT flag.    */
/******************************************************************************/

{
    register double emaptotal, electotal;
    register double u,   v,   w;
    register double p0u, p0v, p0w;
    register double p1u, p1v, p1w;
    register int AtomType;        /* atom type */
    register int u0,  v0,  w0;
    register int u1,  v1,  w1;
    register int i;               /* i-th atom */

#ifdef MINPOINT
    register int x,y,z;                                                    /*MINPOINT*/
#else
    register double e, m; 
#endif

    emaptotal = electotal = 0.;

    for (i=0; i<total_atoms; i++) {

        AtomType = type[i];

        u1  = (u0 = (int) (u = ((double)tcoord[i][X]-(double)xlo) * (double)inv_spacing)) + 1;
        p1u = 1. - (p0u = u - (double) u0);

        v1  = (v0 = (int) (v = ((double)tcoord[i][Y]-(double)ylo) * (double)inv_spacing)) + 1;
        p1v = 1. - (p0v = v - (double) v0);

        w1  = (w0 = (int) (w = ((double)tcoord[i][Z]-(double)zlo) * (double)inv_spacing)) + 1;
        p1w = 1. - (p0w = w - (double) w0);

#ifdef MINPOINT
    x = (p0u < p1u)? u0 : u1;				    /*MINPOINT*/
    y = (p0v < p1v)? v0 : v1;				    /*MINPOINT*/
    z = (p0w < p1w)? w0 : w1;				    /*MINPOINT*/
        						    /*MINPOINT*/
        electotal += (elec[i] = map[z][y][x][ElecMap] * charge[i]); /*MINPOINT*/
        emaptotal += (emap[i] = map[z][y][x][AtomType]); 	    /*MINPOINT*/
#else
        e = m = 0.;

        e += p1u * p1v * p1w * map[ w0 ][ v0 ][ u0 ][ElecMap];
        m += p1u * p1v * p1w * map[ w0 ][ v0 ][ u0 ][AtomType];

        m += p0u * p1v * p1w * map[ w0 ][ v0 ][ u1 ][AtomType];
        e += p0u * p1v * p1w * map[ w0 ][ v0 ][ u1 ][ElecMap];

        e += p1u * p0v * p1w * map[ w0 ][ v1 ][ u0 ][ElecMap];
        m += p1u * p0v * p1w * map[ w0 ][ v1 ][ u0 ][AtomType];

        m += p0u * p0v * p1w * map[ w0 ][ v1 ][ u1 ][AtomType];
        e += p0u * p0v * p1w * map[ w0 ][ v1 ][ u1 ][ElecMap];

        e += p1u * p1v * p0w * map[ w1 ][ v0 ][ u0 ][ElecMap];
        m += p1u * p1v * p0w * map[ w1 ][ v0 ][ u0 ][AtomType];

        m += p0u * p1v * p0w * map[ w1 ][ v0 ][ u1 ][AtomType];
        e += p0u * p1v * p0w * map[ w1 ][ v0 ][ u1 ][ElecMap];

        e += p1u * p0v * p0w * map[ w1 ][ v1 ][ u0 ][ElecMap];
        m += p1u * p0v * p0w * map[ w1 ][ v1 ][ u0 ][AtomType];

        m += p0u * p0v * p0w * map[ w1 ][ v1 ][ u1 ][AtomType];
        e += p0u * p0v * p0w * map[ w1 ][ v1 ][ u1 ][ElecMap];

        electotal += (elec[i] = e * charge[i]);
        emaptotal += (emap[i] = m); 

#endif /* not MINPOINT */

    }/*for  0 <= i < total_atoms*/

    return( (FloatOrDouble)(electotal + emaptotal) );
}

/* } */

/*----------------------------------------------------------------------------*/
/* quicktrilinterp.c */

/* { */

FloatOrDouble quicktrilinterp(    CONST_FLOAT tcoord[MAX_ATOMS][SPACE], 
                          CONST_FLOAT charge[MAX_ATOMS], 
                          CONST_INT   type[MAX_ATOMS], 
                          CONST_INT   total_atoms, 
                          CONST_FLOAT map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS], 
                          CONST_FLOAT inv_spacing, 
                          CONST_FLOAT xlo,
                          CONST_FLOAT ylo,
                          CONST_FLOAT zlo )

{
    register double etotal;
    register double u,   v,   w;
    register double p0u, p0v, p0w;
    register double p1u, p1v, p1w;
    register int AtomType;
    register int u0,  v0,  w0;
    register int u1,  v1,  w1;
    register int i;                /* i-th atom */

#ifdef MINPOINT
    register int x,y,z;
#else
    register double e, m; 
#endif

    etotal = 0.;

    for (i=0; i<total_atoms; i++) {

#ifdef DEBUG
    // gmm  19-FEB-2003
    (void)fprintf(logFile, "trilinterp.cc/quicktrilinterp(...)  atom i= %d\n", i);
#endif /* DEBUG */

        AtomType = type[i];

        u1  = (u0 = (int) (u = ((double)tcoord[i][X]-(double)xlo) * (double)inv_spacing)) + 1;
        p1u = 1. - (p0u = u - (double) u0);

        v1  = (v0 = (int) (v = ((double)tcoord[i][Y]-(double)ylo) * (double)inv_spacing)) + 1;
        p1v = 1. - (p0v = v - (double) v0);

        w1  = (w0 = (int) (w = ((double)tcoord[i][Z]-(double)zlo) * (double)inv_spacing)) + 1;
        p1w = 1. - (p0w = w - (double) w0);

#ifdef DEBUG
    // gmm  19-FEB-2003
    (void)fprintf(logFile, "trilinterp.cc/quicktrilinterp(...)\ttcoord[i][X]= %.9f\n\t\t\t\t\ttcoord[i][Y]= %.9f\n\t\t\t\t\ttcoord[i][Z]= %.9f\n", tcoord[i][X], tcoord[i][Y], tcoord[i][Z]);
    (void)fprintf(logFile, "trilinterp.cc/quicktrilinterp(...)\txlo= %.9f\n\t\t\t\t\tylo= %.9f\n\t\t\t\t\tzlo= %.9f\n", xlo, ylo, zlo);
    (void)fprintf(logFile, "trilinterp.cc/quicktrilinterp(...)\tinv_spacing= %.9f\n", inv_spacing);

    (void)fprintf(logFile, "trilinterp.cc/quicktrilinterp(...)\tu= %.9lf\n\t\t\t\t\tu0= %d\n\t\t\t\t\tu1= %d\n\t\t\t\t\tp0u= %.9lf\n\t\t\t\t\tp1u= %.9lf\n", u, u0, u1, p0u, p1u);
    (void)fprintf(logFile, "trilinterp.cc/quicktrilinterp(...)\tv= %.9lf\n\t\t\t\t\tv0= %d\n\t\t\t\t\tv1= %d\n\t\t\t\t\tp0v= %.9lf\n\t\t\t\t\tp1v= %.9lf\n", v, v0, v1, p0v, p1v);
    (void)fprintf(logFile, "trilinterp.cc/quicktrilinterp(...)\tw= %.9lf\n\t\t\t\t\tw0= %d\n\t\t\t\t\tw1= %d\n\t\t\t\t\tp0w= %.9lf\n\t\t\t\t\tp1w= %.9lf\n", w, w0, w1, p0w, p1w);
#endif /* DEBUG */

#ifdef MINPOINT
        x = (p0u < p1u)? u0 : u1;
        y = (p0v < p1v)? v0 : v1;
        z = (p0w < p1w)? w0 : w1;

        etotal += map[z][y][x][ElecMap] * charge[i] + map[z][y][x][AtomType]; 
#else
        // e = m = 0.;

        e = p1u * p1v * p1w * map[ w0 ][ v0 ][ u0 ][ElecMap];
        m = p1u * p1v * p1w * map[ w0 ][ v0 ][ u0 ][AtomType];

#ifdef DEBUG
    // gmm  19-FEB-2003
    (void)fprintf(logFile, "trilinterp.cc/quicktrilinterp(...)  1  map[ w0 ][v0 ][ u0 ][ElecMap]= %.9f\n", map[ w0 ][ v0 ][ u0 ][ElecMap]);
    (void)fprintf(logFile, "trilinterp.cc/quicktrilinterp(...)  1  map[ w0 ][v0 ][ u0 ][AtomType]= %.9f\n", map[ w0 ][ v0 ][ u0 ][AtomType]);
    (void)fprintf(logFile, "trilinterp.cc/quicktrilinterp(...)  1  p1u * p1v * p1w= %.9lf\n", p1u * p1v * p1w);
    (void)fprintf(logFile, "trilinterp.cc/quicktrilinterp(...)  1  e= %.9lf\n\t\t\t\t\t m= %.9lf\n", e, m);
#endif /* DEBUG */

        m += p0u * p1v * p1w * map[ w0 ][ v0 ][ u1 ][AtomType];
        e += p0u * p1v * p1w * map[ w0 ][ v0 ][ u1 ][ElecMap];

#ifdef DEBUG // gmm  19-FEB-2003
    (void)fprintf(logFile, "trilinterp.cc/quicktrilinterp(...)  2  e= %.9lf\n\t\t\t\t m= %.9lf\n", e, m);
#endif /* DEBUG */

        e += p1u * p0v * p1w * map[ w0 ][ v1 ][ u0 ][ElecMap];
        m += p1u * p0v * p1w * map[ w0 ][ v1 ][ u0 ][AtomType];

#ifdef DEBUG // gmm  19-FEB-2003
    (void)fprintf(logFile, "trilinterp.cc/quicktrilinterp(...)  3  e= %.9lf\n\t\t\t\t m= %.9lf\n", e, m);
#endif /* DEBUG */

        m += p0u * p0v * p1w * map[ w0 ][ v1 ][ u1 ][AtomType];
        e += p0u * p0v * p1w * map[ w0 ][ v1 ][ u1 ][ElecMap];

#ifdef DEBUG // gmm  19-FEB-2003
    (void)fprintf(logFile, "trilinterp.cc/quicktrilinterp(...)  4  e= %.9lf\n\t\t\t\t m= %.9lf\n", e, m);
#endif /* DEBUG */

        e += p1u * p1v * p0w * map[ w1 ][ v0 ][ u0 ][ElecMap];
        m += p1u * p1v * p0w * map[ w1 ][ v0 ][ u0 ][AtomType];

#ifdef DEBUG // gmm  19-FEB-2003
    (void)fprintf(logFile, "trilinterp.cc/quicktrilinterp(...)  5  e= %.9lf\n\t\t\t\t m= %.9lf\n", e, m);
#endif /* DEBUG */

        m += p0u * p1v * p0w * map[ w1 ][ v0 ][ u1 ][AtomType];
        e += p0u * p1v * p0w * map[ w1 ][ v0 ][ u1 ][ElecMap];
    
#ifdef DEBUG // gmm  19-FEB-2003
    (void)fprintf(logFile, "trilinterp.cc/quicktrilinterp(...)  6  e= %.9lf\n\t\t\t\t m= %.9lf\n", e, m);
#endif /* DEBUG */
    
        e += p1u * p0v * p0w * map[ w1 ][ v1 ][ u0 ][ElecMap];
        m += p1u * p0v * p0w * map[ w1 ][ v1 ][ u0 ][AtomType];

#ifdef DEBUG // gmm  19-FEB-2003
    (void)fprintf(logFile, "trilinterp.cc/quicktrilinterp(...)  7  e= %.9lf\n\t\t\t\t m= %.9lf\n", e, m);
#endif /* DEBUG */

        m += p0u * p0v * p0w * map[ w1 ][ v1 ][ u1 ][AtomType];
        e += p0u * p0v * p0w * map[ w1 ][ v1 ][ u1 ][ElecMap];

#ifdef DEBUG // gmm  19-FEB-2003
    (void)fprintf(logFile, "trilinterp.cc/quicktrilinterp(...)  8  e= %.9lf\n\t\t\t\t m= %.9lf\n", e, m);
#endif /* DEBUG */

        etotal += e * charge[i] + m; 

#ifdef DEBUG // gmm  19-FEB-2003
    (void)fprintf(logFile, "trilinterp.cc/quicktrilinterp(...)  9  etotal= %.9lf\n", etotal);
    (void)fprintf(logFile, "trilinterp.cc/quicktrilinterp(...)  9  (FloatOrDouble)etotal= %.5f\n", (FloatOrDouble)etotal);
#endif /* DEBUG */

#endif /* not MINPOINT */

    }/*for  0 <= i < total_atoms*/

    return( (FloatOrDouble)etotal );
}
/* } */

/*----------------------------------------------------------------------------*/

/* { */

FloatOrDouble outsidetrilinterp(CONST_FLOAT tcoord[MAX_ATOMS][SPACE], 
          	CONST_FLOAT charge[MAX_ATOMS], 
          	CONST_INT   type[MAX_ATOMS], 
          	CONST_INT   total_atoms, 
          	CONST_FLOAT map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS], 
          	CONST_FLOAT inv_spacing, 
                // FloatOrDouble elec[MAX_ATOMS],
        	      // FloatOrDouble emap[MAX_ATOMS],
          	CONST_FLOAT xlo,
          	CONST_FLOAT ylo,
          	CONST_FLOAT zlo,
          	CONST_FLOAT xhi,
          	CONST_FLOAT yhi,
          	CONST_FLOAT zhi,
          	CONST_FLOAT xcen,
          	CONST_FLOAT ycen,
          	CONST_FLOAT zcen )

{
    double	 etotal, epenalty;
    // CONST_FLOAT	 x, y, z; // Xcode-gmm was double 2004-02-03
    double	 x, y, z; // Xcode-gmm was double 2004-02-03
    double	 u,   v,   w;
    double	 p0u, p0v, p0w;
    double	 p1u, p1v, p1w;

    int		 AtomType;
    int		 u0,  v0,  w0;
    int		 u1,  v1,  w1;

    register int i;		/* i-th atom */

#ifdef MINPOINT
    int		 x,y,z;
#else
    double 	 e, m; 
#endif

    etotal = 0.;

    for (i=0; i<total_atoms; i++) {
        x = tcoord[i][X];
        y = tcoord[i][Y];
        z = tcoord[i][Z];

        if (is_out_grid(x,y,z)) {
            x -= xcen;
            y -= ycen;
            z -= zcen;
            // sqhypotenuse(x,y,z) is the square of the distance 
            // from grid's centre to atom
            epenalty = sqhypotenuse(x,y,z) * ENERGYPENALTY;
            // etotal += (elec[i] = epenalty) + (emap[i] = epenalty);
            etotal += (epenalty) + (epenalty);
        } else {
            AtomType = type[i];
     
            u1  = (u0 = (int) (u = ((double)tcoord[i][X]-(double)xlo) * (double)inv_spacing)) + 1;
            p1u = 1. - (p0u = u - (double) u0);
     
            v1  = (v0 = (int) (v = ((double)tcoord[i][Y]-(double)ylo) * (double)inv_spacing)) + 1;
            p1v = 1. - (p0v = v - (double) v0);
     
            w1  = (w0 = (int) (w = ((double)tcoord[i][Z]-(double)zlo) * (double)inv_spacing)) + 1;
            p1w = 1. - (p0w = w - (double) w0);
     
#ifdef MINPOINT
            x = (p0u < p1u)? u0 : u1;
            y = (p0v < p1v)? v0 : v1;
            z = (p0w < p1w)? w0 : w1;
     
            // etotal += (elec[i] = map[z][y][x][ElecMap] * charge[i]) + 
                      // (emap[i] = map[z][y][x][AtomType]); 
            etotal += (map[z][y][x][ElecMap] * charge[i]) + 
                      (map[i] = map[z][y][x][AtomType]); 

#else
            // e = m = 0.;
     
            e = p1u * p1v * p1w * map[ w0 ][ v0 ][ u0 ][ElecMap];
            m = p1u * p1v * p1w * map[ w0 ][ v0 ][ u0 ][AtomType];
     
            m += p0u * p1v * p1w * map[ w0 ][ v0 ][ u1 ][AtomType];
            e += p0u * p1v * p1w * map[ w0 ][ v0 ][ u1 ][ElecMap];
     
            e += p1u * p0v * p1w * map[ w0 ][ v1 ][ u0 ][ElecMap];
            m += p1u * p0v * p1w * map[ w0 ][ v1 ][ u0 ][AtomType];
     
            m += p0u * p0v * p1w * map[ w0 ][ v1 ][ u1 ][AtomType];
            e += p0u * p0v * p1w * map[ w0 ][ v1 ][ u1 ][ElecMap];

            e += p1u * p1v * p0w * map[ w1 ][ v0 ][ u0 ][ElecMap];
            m += p1u * p1v * p0w * map[ w1 ][ v0 ][ u0 ][AtomType];
     
            m += p1u * p0v * p0w * map[ w1 ][ v1 ][ u0 ][AtomType];
            e += p1u * p0v * p0w * map[ w1 ][ v1 ][ u0 ][ElecMap];
     
            e += p0u * p1v * p0w * map[ w1 ][ v0 ][ u1 ][ElecMap];
            m += p0u * p1v * p0w * map[ w1 ][ v0 ][ u1 ][AtomType];
     
            m += p0u * p0v * p0w * map[ w1 ][ v1 ][ u1 ][AtomType];
            e += p0u * p0v * p0w * map[ w1 ][ v1 ][ u1 ][ElecMap];

            // etotal += (elec[i] = e * charge[i]) + (emap[i] = m); 
            etotal += (e * charge[i]) + (m); 

#endif /* not MINPOINT */

        } /* inside grid */
    }/*for  0 <= i < total_atoms*/

    return( (FloatOrDouble)etotal );
}

// added 07-nov-2001

/*----------------------------------------------------------------------------*/

FloatOrDouble outsidetrilinterpbyatom(CONST_FLOAT tcoord[MAX_ATOMS][SPACE], 
                          CONST_FLOAT charge[MAX_ATOMS], 
                          CONST_INT   type[MAX_ATOMS], 
                          CONST_INT   total_atoms, 
                          CONST_FLOAT map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS], 
                          CONST_FLOAT inv_spacing, 
                          FloatOrDouble elec[MAX_ATOMS],
                          FloatOrDouble emap[MAX_ATOMS],
                          CONST_FLOAT xlo,
                          CONST_FLOAT ylo,
                          CONST_FLOAT zlo,
                          CONST_FLOAT xhi,
                          CONST_FLOAT yhi,
                          CONST_FLOAT zhi,
                          CONST_FLOAT xcen,
                          CONST_FLOAT ycen,
                          CONST_FLOAT zcen )

{
    register double etotal, epenalty;
    register double x, y, z;
    register double u,   v,   w;
    register double p0u, p0v, p0w;
    register double p1u, p1v, p1w;
    register int AtomType;
    register int u0,  v0,  w0;
    register int u1,  v1,  w1;
    register int i;                /* i-th atom */

#ifdef MINPOINT
    register int x,y,z;
#else
    register double e, m; 
#endif

    etotal = 0.;

    for (i=0; i<total_atoms; i++) {
        x = tcoord[i][X];
        y = tcoord[i][Y];
        z = tcoord[i][Z];

        if (is_out_grid(x,y,z)) {
            x -= xcen;
            y -= ycen;
            z -= zcen;
            // sqhypotenuse(x,y,z) is the square of the distance 
            // from grid's centre to atom
            epenalty = sqhypotenuse(x,y,z) * ENERGYPENALTY;
            etotal += (elec[i] = epenalty) + (emap[i] = epenalty);
        } else {
            AtomType = type[i];
     
            u1  = (u0 = (int) (u = ((double)tcoord[i][X]-(double)xlo) * (double)inv_spacing)) + 1;
            p1u = 1. - (p0u = u - (double) u0);
     
            v1  = (v0 = (int) (v = ((double)tcoord[i][Y]-(double)ylo) * (double)inv_spacing)) + 1;
            p1v = 1. - (p0v = v - (double) v0);
     
            w1  = (w0 = (int) (w = ((double)tcoord[i][Z]-(double)zlo) * (double)inv_spacing)) + 1;
            p1w = 1. - (p0w = w - (double) w0);
     
#ifdef MINPOINT
            x = (p0u < p1u)? u0 : u1;
            y = (p0v < p1v)? v0 : v1;
            z = (p0w < p1w)? w0 : w1;
     
            etotal += (elec[i] = map[z][y][x][ElecMap] * charge[i]) + 
                      (emap[i] = map[z][y][x][AtomType]); 

#else
            e = m = 0.;
     
            e += p1u * p1v * p1w * map[ w0 ][ v0 ][ u0 ][ElecMap];
            m += p1u * p1v * p1w * map[ w0 ][ v0 ][ u0 ][AtomType];
     
            m += p0u * p1v * p1w * map[ w0 ][ v0 ][ u1 ][AtomType];
            e += p0u * p1v * p1w * map[ w0 ][ v0 ][ u1 ][ElecMap];
     
            e += p1u * p0v * p1w * map[ w0 ][ v1 ][ u0 ][ElecMap];
            m += p1u * p0v * p1w * map[ w0 ][ v1 ][ u0 ][AtomType];
     
            m += p1u * p1v * p0w * map[ w1 ][ v0 ][ u0 ][AtomType];
            e += p1u * p1v * p0w * map[ w1 ][ v0 ][ u0 ][ElecMap];
     
            e += p0u * p0v * p1w * map[ w0 ][ v1 ][ u1 ][ElecMap];
            m += p0u * p0v * p1w * map[ w0 ][ v1 ][ u1 ][AtomType];
     
            m += p1u * p0v * p0w * map[ w1 ][ v1 ][ u0 ][AtomType];
            e += p1u * p0v * p0w * map[ w1 ][ v1 ][ u0 ][ElecMap];
     
            e += p0u * p1v * p0w * map[ w1 ][ v0 ][ u1 ][ElecMap];
            m += p0u * p1v * p0w * map[ w1 ][ v0 ][ u1 ][AtomType];
     
            m += p0u * p0v * p0w * map[ w1 ][ v1 ][ u1 ][AtomType];
            e += p0u * p0v * p0w * map[ w1 ][ v1 ][ u1 ][ElecMap];

            etotal += (elec[i] = e * charge[i]) + (emap[i] = m); 

#endif /* not MINPOINT */

        } /* inside grid */
    }/*for  0 <= i < total_atoms*/

    return( (FloatOrDouble)etotal );
} /* outsidetrilinterpbyatom */
/* } */

// added 15-jan-2001

/* { */
/*----------------------------------------------------------------------------*/

FloatOrDouble template_trilinterp( CONST_FLOAT tcoord[MAX_ATOMS][SPACE], 
                             CONST_FLOAT charge[MAX_ATOMS], 
                             CONST_INT   type[MAX_ATOMS], 
                             CONST_INT   total_atoms, 
                             CONST_FLOAT map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS], 
                             CONST_FLOAT inv_spacing, 
                             CONST_FLOAT xlo,
                             CONST_FLOAT ylo,
                             CONST_FLOAT zlo,
                           CONST_FLOAT template_energy[MAX_ATOMS],
                           CONST_FLOAT template_stddev[MAX_ATOMS])

{
    double	 etotal;
    double	 u,   v,   w;
    double	 p0u, p0v, p0w;
    double	 p1u, p1v, p1w;
    double diff;

    int		 AtomType;
    int		 u0,  v0,  w0;
    int		 u1,  v1,  w1;

    register int i;		/* i-th atom */

#ifdef MINPOINT
    int		 x,y,z;
#else
    double 	 e;  // electrostatic energy
    double 	 m; // affinity map
#endif

    etotal = 0.;

    for (i=0; i<total_atoms; i++) {

        AtomType = type[i];

        u1  = (u0 = (int) (u = ((double)tcoord[i][X]-(double)xlo) * (double)inv_spacing)) + 1;
        p1u = 1. - (p0u = u - (double) u0);

        v1  = (v0 = (int) (v = ((double)tcoord[i][Y]-(double)ylo) * (double)inv_spacing)) + 1;
        p1v = 1. - (p0v = v - (double) v0);

        w1  = (w0 = (int) (w = ((double)tcoord[i][Z]-(double)zlo) * (double)inv_spacing)) + 1;
        p1w = 1. - (p0w = w - (double) w0);

#ifdef MINPOINT
        x = (p0u < p1u)? u0 : u1;
        y = (p0v < p1v)? v0 : v1;
        z = (p0w < p1w)? w0 : w1;

        diff = (map[z][y][x][AtomType] + map[z][y][x][ElecMap] * charge[i]  - template_energy[i]) / template_stddev[i]; 
        etotal +=  + diff * diff;
#else
        e = m = 0.;

        e += p1u * p1v * p1w * map[ w0 ][ v0 ][ u0 ][ElecMap];
        m += p1u * p1v * p1w * map[ w0 ][ v0 ][ u0 ][AtomType];

        m += p0u * p1v * p1w * map[ w0 ][ v0 ][ u1 ][AtomType];
        e += p0u * p1v * p1w * map[ w0 ][ v0 ][ u1 ][ElecMap];

        e += p1u * p0v * p1w * map[ w0 ][ v1 ][ u0 ][ElecMap];
        m += p1u * p0v * p1w * map[ w0 ][ v1 ][ u0 ][AtomType];

        m += p1u * p1v * p0w * map[ w1 ][ v0 ][ u0 ][AtomType];
        e += p1u * p1v * p0w * map[ w1 ][ v0 ][ u0 ][ElecMap];

        e += p0u * p0v * p1w * map[ w0 ][ v1 ][ u1 ][ElecMap];
        m += p0u * p0v * p1w * map[ w0 ][ v1 ][ u1 ][AtomType];

        m += p1u * p0v * p0w * map[ w1 ][ v1 ][ u0 ][AtomType];
        e += p1u * p0v * p0w * map[ w1 ][ v1 ][ u0 ][ElecMap];

        e += p0u * p1v * p0w * map[ w1 ][ v0 ][ u1 ][ElecMap];
        m += p0u * p1v * p0w * map[ w1 ][ v0 ][ u1 ][AtomType];

        m += p0u * p0v * p0w * map[ w1 ][ v1 ][ u1 ][AtomType];
        e += p0u * p0v * p0w * map[ w1 ][ v1 ][ u1 ][ElecMap];

        diff = (m + e * charge[i] - template_energy[i]) / template_stddev[i]; 
        etotal += diff * diff;

#endif /* not MINPOINT */

    }/*for  0 <= i < total_atoms*/

    return((FloatOrDouble)sqrt(etotal / (FloatOrDouble)total_atoms ));
} // FloatOrDouble template_trilinterp( CONST_FLOAT tcoord[MAX_ATOMS][SPACE], 

/*----------------------------------------------------------------------------*/

FloatOrDouble outside_templ_trilinterp(CONST_FLOAT tcoord[MAX_ATOMS][SPACE], 
          	CONST_FLOAT charge[MAX_ATOMS], 
          	CONST_INT   type[MAX_ATOMS], 
          	CONST_INT   total_atoms, 
          	CONST_FLOAT map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS], 
          	CONST_FLOAT inv_spacing, 
          	CONST_FLOAT xlo,
          	CONST_FLOAT ylo,
          	CONST_FLOAT zlo,
          	CONST_FLOAT xhi,
          	CONST_FLOAT yhi,
          	CONST_FLOAT zhi,
          	CONST_FLOAT xcen,
          	CONST_FLOAT ycen,
          	CONST_FLOAT zcen,
            CONST_FLOAT template_energy[MAX_ATOMS],
            CONST_FLOAT template_stddev[MAX_ATOMS])

{
    double	 etotal, epenalty;
    double	 x, y, z;
    double	 u,   v,   w;
    double	 p0u, p0v, p0w;
    double	 p1u, p1v, p1w;
    double   diff;

    int		 AtomType;
    int		 u0,  v0,  w0;
    int		 u1,  v1,  w1;

    register int i;                /* i-th atom */

#ifdef MINPOINT
    int		 x,y,z;
#else
    double 	 e, m; 
#endif

    etotal = 0.;

    for (i=0; i<total_atoms; i++) {
        x = tcoord[i][X];
        y = tcoord[i][Y];
        z = tcoord[i][Z];

        if (is_out_grid(x,y,z)) {
            x -= xcen;
            y -= ycen;
            z -= zcen;
            // sqhypotenuse(x,y,z) is the square of the distance 
            // from grid's centre to atom
            epenalty = sqhypotenuse(x,y,z) * ENERGYPENALTY;
            // etotal += (elec[i] = epenalty) + (emap[i] = epenalty);
            etotal += (epenalty) + (epenalty);
        } else {
            AtomType = type[i];
     
            u1  = (u0 = (int) (u = ((double)tcoord[i][X]-(double)xlo) * (double)inv_spacing)) + 1;
            p1u = 1. - (p0u = u - (double) u0);
     
            v1  = (v0 = (int) (v = ((double)tcoord[i][Y]-(double)ylo) * (double)inv_spacing)) + 1;
            p1v = 1. - (p0v = v - (double) v0);
     
            w1  = (w0 = (int) (w = ((double)tcoord[i][Z]-(double)zlo) * (double)inv_spacing)) + 1;
            p1w = 1. - (p0w = w - (double) w0);
     
#ifdef MINPOINT
            x = (p0u < p1u)? u0 : u1;
            y = (p0v < p1v)? v0 : v1;
            z = (p0w < p1w)? w0 : w1;
     
        diff = (map[z][y][x][AtomType] + map[z][y][x][ElecMap] * charge[i] - template_energy[i]) / template_stddev[i];
        etotal += diff * diff;
#else
            e = m = 0.;
     
            e += p1u * p1v * p1w * map[ w0 ][ v0 ][ u0 ][ElecMap];
            m += p1u * p1v * p1w * map[ w0 ][ v0 ][ u0 ][AtomType];
     
            m += p0u * p1v * p1w * map[ w0 ][ v0 ][ u1 ][AtomType];
            e += p0u * p1v * p1w * map[ w0 ][ v0 ][ u1 ][ElecMap];
     
            e += p1u * p0v * p1w * map[ w0 ][ v1 ][ u0 ][ElecMap];
            m += p1u * p0v * p1w * map[ w0 ][ v1 ][ u0 ][AtomType];
     
            m += p1u * p1v * p0w * map[ w1 ][ v0 ][ u0 ][AtomType];
            e += p1u * p1v * p0w * map[ w1 ][ v0 ][ u0 ][ElecMap];
     
            e += p0u * p0v * p1w * map[ w0 ][ v1 ][ u1 ][ElecMap];
            m += p0u * p0v * p1w * map[ w0 ][ v1 ][ u1 ][AtomType];
     
            m += p1u * p0v * p0w * map[ w1 ][ v1 ][ u0 ][AtomType];
            e += p1u * p0v * p0w * map[ w1 ][ v1 ][ u0 ][ElecMap];
     
            e += p0u * p1v * p0w * map[ w1 ][ v0 ][ u1 ][ElecMap];
            m += p0u * p1v * p0w * map[ w1 ][ v0 ][ u1 ][AtomType];
     
            m += p0u * p0v * p0w * map[ w1 ][ v1 ][ u1 ][AtomType];
            e += p0u * p0v * p0w * map[ w1 ][ v1 ][ u1 ][ElecMap];

        diff = (m + e * charge[i] - template_energy[i]) / template_stddev[i];
            etotal += diff * diff; 

#endif /* not MINPOINT */

        } /* inside grid */
    }/*for  0 <= i < total_atoms*/

    return( (FloatOrDouble)sqrt(etotal / (FloatOrDouble)total_atoms));
}

/*----------------------------------------------------------------------------*/

FloatOrDouble byatom_template_trilinterp( CONST_FLOAT tcoord[MAX_ATOMS][SPACE], 
          	                      CONST_FLOAT charge[MAX_ATOMS], 
          	                      CONST_INT   type[MAX_ATOMS], 
          	                      CONST_INT   total_atoms, 
          	                      CONST_FLOAT map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS], 
          	                      CONST_FLOAT inv_spacing, 
                                  FloatOrDouble elec[MAX_ATOMS], 
                                  FloatOrDouble emap[MAX_ATOMS], 
          	                      CONST_FLOAT xlo,
          	                      CONST_FLOAT ylo,
          	                      CONST_FLOAT zlo,
                                  CONST_FLOAT template_energy[MAX_ATOMS],
                                  CONST_FLOAT template_stddev[MAX_ATOMS])
{
    double	 etotal;
    double	 u,   v,   w;
    double	 p0u, p0v, p0w;
    double	 p1u, p1v, p1w;
    double diff;

    int		 AtomType;
    int		 u0,  v0,  w0;
    int		 u1,  v1,  w1;

    register int i;		/* i-th atom */

#ifdef MINPOINT
    int		 x,y,z;
#else
    double 	 e;  // electrostatic energy
    double 	 m; // affinity map
#endif

    etotal = 0.;

    for (i=0; i<total_atoms; i++) {

        AtomType = type[i];

        u1  = (u0 = (int) (u = ((double)tcoord[i][X]-(double)xlo) * (double)inv_spacing)) + 1;
        p1u = 1. - (p0u = u - (double) u0);

        v1  = (v0 = (int) (v = ((double)tcoord[i][Y]-(double)ylo) * (double)inv_spacing)) + 1;
        p1v = 1. - (p0v = v - (double) v0);

        w1  = (w0 = (int) (w = ((double)tcoord[i][Z]-(double)zlo) * (double)inv_spacing)) + 1;
        p1w = 1. - (p0w = w - (double) w0);

#ifdef MINPOINT
        x = (p0u < p1u)? u0 : u1;
        y = (p0v < p1v)? v0 : v1;
        z = (p0w < p1w)? w0 : w1;

        diff = ((emap[i] = map[z][y][x][AtomType]) + (elec[i] = map[z][y][x][ElecMap] * charge[i]) - template_energy[i]) / template_stddev[i]; 
        etotal += diff * diff;
#else
        e = m = 0.;

        e += p1u * p1v * p1w * map[ w0 ][ v0 ][ u0 ][ElecMap];
        m += p1u * p1v * p1w * map[ w0 ][ v0 ][ u0 ][AtomType];

        m += p0u * p1v * p1w * map[ w0 ][ v0 ][ u1 ][AtomType];
        e += p0u * p1v * p1w * map[ w0 ][ v0 ][ u1 ][ElecMap];

        e += p1u * p0v * p1w * map[ w0 ][ v1 ][ u0 ][ElecMap];
        m += p1u * p0v * p1w * map[ w0 ][ v1 ][ u0 ][AtomType];

        m += p1u * p1v * p0w * map[ w1 ][ v0 ][ u0 ][AtomType];
        e += p1u * p1v * p0w * map[ w1 ][ v0 ][ u0 ][ElecMap];

        e += p0u * p0v * p1w * map[ w0 ][ v1 ][ u1 ][ElecMap];
        m += p0u * p0v * p1w * map[ w0 ][ v1 ][ u1 ][AtomType];

        m += p1u * p0v * p0w * map[ w1 ][ v1 ][ u0 ][AtomType];
        e += p1u * p0v * p0w * map[ w1 ][ v1 ][ u0 ][ElecMap];

        e += p0u * p1v * p0w * map[ w1 ][ v0 ][ u1 ][ElecMap];
        m += p0u * p1v * p0w * map[ w1 ][ v0 ][ u1 ][AtomType];

        m += p0u * p0v * p0w * map[ w1 ][ v1 ][ u1 ][AtomType];
        e += p0u * p0v * p0w * map[ w1 ][ v1 ][ u1 ][ElecMap];

        diff = ((emap[i] = m) + (elec[i] = e * charge[i]) - template_energy[i]) / template_stddev[i];
        etotal += diff * diff;

#endif /* not MINPOINT */

    }/*for  0 <= i < total_atoms*/

    return( (FloatOrDouble)sqrt(etotal / (FloatOrDouble)total_atoms));
}

//#------------------------------------------------------------------

/* { */

FloatOrDouble trilinterp4( CONST_FLOAT tcoord[MAX_ATOMS][SPACE], 
                  CONST_FLOAT charge[MAX_ATOMS], 
                  CONST_INT   type[MAX_ATOMS], 
                  CONST_INT   total_atoms, 
                  CONST_FLOAT map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS],
                  CONST_FLOAT inv_spacing, 
                  FloatOrDouble elec[MAX_ATOMS], 
                  FloatOrDouble emap[MAX_ATOMS], 
                  CONST_FLOAT xlo,
                  CONST_FLOAT ylo,
                  CONST_FLOAT zlo,
                  int ignore_inter[MAX_ATOMS])

{
    register double emaptotal, electotal;
    register double u,   v,   w;
    register double p0u, p0v, p0w;
    register double p1u, p1v, p1w;
    register int AtomType;        /* atom type */
    register int u0,  v0,  w0;
    register int u1,  v1,  w1;
    register int i;               /* i-th atom */

#ifdef MINPOINT
    register int x,y,z;                                                    /*MINPOINT*/
#else
    register double e, m; 
#endif

    emaptotal = electotal = 0.;

    for (i=0; i<total_atoms; i++) {

        if (!ignore_inter[i]) {

            AtomType = type[i];

            u1  = (u0 = (int) (u = (tcoord[i][X]-xlo) * inv_spacing)) + 1;
            p1u = 1. - (p0u = u - (double) u0);

            v1  = (v0 = (int) (v = (tcoord[i][Y]-ylo) * inv_spacing)) + 1;
            p1v = 1. - (p0v = v - (double) v0);

            w1  = (w0 = (int) (w = (tcoord[i][Z]-zlo) * inv_spacing)) + 1;
            p1w = 1. - (p0w = w - (double) w0);

    #ifdef MINPOINT
            x = (p0u < p1u)? u0 : u1;                                    /*MINPOINT*/
            y = (p0v < p1v)? v0 : v1;                                    /*MINPOINT*/
            z = (p0w < p1w)? w0 : w1;                                    /*MINPOINT*/
                                                                        /*MINPOINT*/
            electotal += (elec[i] = map[z][y][x][ElecMap] * charge[i]); /*MINPOINT*/
            emaptotal += (emap[i] = map[z][y][x][AtomType]);             /*MINPOINT*/
    #else
            e = m = 0.;

            e += p1u * p1v * p1w * map[ w0 ][ v0 ][ u0 ][ElecMap];
            m += p1u * p1v * p1w * map[ w0 ][ v0 ][ u0 ][AtomType];

            m += p0u * p1v * p1w * map[ w0 ][ v0 ][ u1 ][AtomType];
            e += p0u * p1v * p1w * map[ w0 ][ v0 ][ u1 ][ElecMap];

            e += p1u * p0v * p1w * map[ w0 ][ v1 ][ u0 ][ElecMap];
            m += p1u * p0v * p1w * map[ w0 ][ v1 ][ u0 ][AtomType];

            m += p1u * p1v * p0w * map[ w1 ][ v0 ][ u0 ][AtomType];
            e += p1u * p1v * p0w * map[ w1 ][ v0 ][ u0 ][ElecMap];

            e += p0u * p0v * p1w * map[ w0 ][ v1 ][ u1 ][ElecMap];
            m += p0u * p0v * p1w * map[ w0 ][ v1 ][ u1 ][AtomType];

            m += p1u * p0v * p0w * map[ w1 ][ v1 ][ u0 ][AtomType];
            e += p1u * p0v * p0w * map[ w1 ][ v1 ][ u0 ][ElecMap];

            e += p0u * p1v * p0w * map[ w1 ][ v0 ][ u1 ][ElecMap];
            m += p0u * p1v * p0w * map[ w1 ][ v0 ][ u1 ][AtomType];

            m += p0u * p0v * p0w * map[ w1 ][ v1 ][ u1 ][AtomType];
            e += p0u * p0v * p0w * map[ w1 ][ v1 ][ u1 ][ElecMap];

            electotal += (elec[i] = e * charge[i]);
            emaptotal += (emap[i] = m); 

    #endif /* not MINPOINT */

        } else {
            elec[i] = 0.0;
            emap[i] = 0.0;
        }/* if (ignore_inter[i]) */

    }/*for  0 <= i < total_atoms*/

    return( (FloatOrDouble)(electotal + emaptotal) );
}

/* } */

//#------------------------------------------------------------------

/* { */

FloatOrDouble outsidetrilinterp4(CONST_FLOAT tcoord[MAX_ATOMS][SPACE],
                          CONST_FLOAT charge[MAX_ATOMS], 
                          CONST_INT   type[MAX_ATOMS], 
                          CONST_INT   total_atoms, 
                          CONST_FLOAT map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS], 
                          CONST_FLOAT inv_spacing, 
                          CONST_FLOAT xlo,
                          CONST_FLOAT ylo,
                          CONST_FLOAT zlo,
                          CONST_FLOAT xhi,
                          CONST_FLOAT yhi,
                          CONST_FLOAT zhi,
                          CONST_FLOAT xcen,
                          CONST_FLOAT ycen,
                          CONST_FLOAT zcen,
                          int ignore_inter[MAX_ATOMS])

{
    register double etotal, epenalty;
    register double x, y, z;
    register double u,   v,   w;
    register double p0u, p0v, p0w;
    register double p1u, p1v, p1w;
    register int AtomType;
    register int u0,  v0,  w0;
    register int u1,  v1,  w1;
    register int i;                /* i-th atom */

#ifdef MINPOINT
    register int x,y,z;
#else
    register double e, m; 
#endif

    etotal = 0.;

    for (i=0; i<total_atoms; i++) {

        if (!ignore_inter[i]) {

            x = tcoord[i][X];
            y = tcoord[i][Y];
            z = tcoord[i][Z];

            if (is_out_grid(x,y,z)) {
                x -= xcen;
                y -= ycen;
                z -= zcen;
                // sqhypotenuse(x,y,z) is the square of the distance 
                // from grid's centre to atom
                epenalty = sqhypotenuse(x,y,z) * ENERGYPENALTY;
                etotal += epenalty + epenalty;
            } else {
                AtomType = type[i];
         
                u1  = (u0 = (int) (u = (tcoord[i][X]-xlo) * inv_spacing)) + 1;
                p1u = 1. - (p0u = u - (double) u0);
         
                v1  = (v0 = (int) (v = (tcoord[i][Y]-ylo) * inv_spacing)) + 1;
                p1v = 1. - (p0v = v - (double) v0);
         
                w1  = (w0 = (int) (w = (tcoord[i][Z]-zlo) * inv_spacing)) + 1;
                p1w = 1. - (p0w = w - (double) w0);
         
    #ifdef MINPOINT
                x = (p0u < p1u)? u0 : u1;
                y = (p0v < p1v)? v0 : v1;
                z = (p0w < p1w)? w0 : w1;
         
                etotal += map[z][y][x][ElecMap] * charge[i] + map[z][y][x][AtomType]; 

    #else
                e = m = 0.;
         
                e += p1u * p1v * p1w * map[ w0 ][ v0 ][ u0 ][ElecMap];
                m += p1u * p1v * p1w * map[ w0 ][ v0 ][ u0 ][AtomType];
         
                m += p0u * p1v * p1w * map[ w0 ][ v0 ][ u1 ][AtomType];
                e += p0u * p1v * p1w * map[ w0 ][ v0 ][ u1 ][ElecMap];
         
                e += p1u * p0v * p1w * map[ w0 ][ v1 ][ u0 ][ElecMap];
                m += p1u * p0v * p1w * map[ w0 ][ v1 ][ u0 ][AtomType];
         
                m += p1u * p1v * p0w * map[ w1 ][ v0 ][ u0 ][AtomType];
                e += p1u * p1v * p0w * map[ w1 ][ v0 ][ u0 ][ElecMap];
         
                e += p0u * p0v * p1w * map[ w0 ][ v1 ][ u1 ][ElecMap];
                m += p0u * p0v * p1w * map[ w0 ][ v1 ][ u1 ][AtomType];
         
                m += p1u * p0v * p0w * map[ w1 ][ v1 ][ u0 ][AtomType];
                e += p1u * p0v * p0w * map[ w1 ][ v1 ][ u0 ][ElecMap];
         
                e += p0u * p1v * p0w * map[ w1 ][ v0 ][ u1 ][ElecMap];
                m += p0u * p1v * p0w * map[ w1 ][ v0 ][ u1 ][AtomType];
         
                m += p0u * p0v * p0w * map[ w1 ][ v1 ][ u1 ][AtomType];
                e += p0u * p0v * p0w * map[ w1 ][ v1 ][ u1 ][ElecMap];

                etotal += e * charge[i] + m; 

    #endif /* not MINPOINT */

            } /* inside grid */

        }

    }/*for  0 <= i < total_atoms*/

    return( (FloatOrDouble)etotal );
}

/* } */

//#------------------------------------------------------------------

/* { */

FloatOrDouble outsidetrilinterp4byatom(CONST_FLOAT tcoord[MAX_ATOMS][SPACE],
                          CONST_FLOAT charge[MAX_ATOMS], 
                          CONST_INT   type[MAX_ATOMS], 
                          CONST_INT   total_atoms, 
                          CONST_FLOAT map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS], 
                          CONST_FLOAT inv_spacing, 
                        FloatOrDouble elec[MAX_ATOMS],
                        FloatOrDouble emap[MAX_ATOMS],
                          CONST_FLOAT xlo,
                          CONST_FLOAT ylo,
                          CONST_FLOAT zlo,
                          CONST_FLOAT xhi,
                          CONST_FLOAT yhi,
                          CONST_FLOAT zhi,
                          CONST_FLOAT xcen,
                          CONST_FLOAT ycen,
                          CONST_FLOAT zcen,
                          int ignore_inter[MAX_ATOMS])

{
    register double etotal, epenalty;
    register double x, y, z;
    register double u,   v,   w;
    register double p0u, p0v, p0w;
    register double p1u, p1v, p1w;
    register int AtomType;
    register int u0,  v0,  w0;
    register int u1,  v1,  w1;
    register int i;                /* i-th atom */

#ifdef MINPOINT
    register int x,y,z;
#else
    register double e, m; 
#endif

    etotal = 0.;

    for (i=0; i<total_atoms; i++) {

        if (!ignore_inter[i]) {

            x = tcoord[i][X];
            y = tcoord[i][Y];
            z = tcoord[i][Z];

            if (is_out_grid(x,y,z)) {
                x -= xcen;
                y -= ycen;
                z -= zcen;
                // sqhypotenuse(x,y,z) is the square of the distance 
                // from grid's centre to atom
                epenalty = sqhypotenuse(x,y,z) * ENERGYPENALTY;
                etotal += (elec[i] = epenalty) + (emap[i] = epenalty);
            } else {
                AtomType = type[i];
         
                u1  = (u0 = (int) (u = (tcoord[i][X]-xlo) * inv_spacing)) + 1;
                p1u = 1. - (p0u = u - (double) u0);
         
                v1  = (v0 = (int) (v = (tcoord[i][Y]-ylo) * inv_spacing)) + 1;
                p1v = 1. - (p0v = v - (double) v0);
         
                w1  = (w0 = (int) (w = (tcoord[i][Z]-zlo) * inv_spacing)) + 1;
                p1w = 1. - (p0w = w - (double) w0);
         
    #ifdef MINPOINT
                x = (p0u < p1u)? u0 : u1;
                y = (p0v < p1v)? v0 : v1;
                z = (p0w < p1w)? w0 : w1;
         
                etotal += (elec[i] = map[z][y][x][ElecMap] * charge[i]) + 
                          (emap[i] = map[z][y][x][AtomType]); 

    #else
                e = m = 0.;
         
                e += p1u * p1v * p1w * map[ w0 ][ v0 ][ u0 ][ElecMap];
                m += p1u * p1v * p1w * map[ w0 ][ v0 ][ u0 ][AtomType];
         
                m += p0u * p1v * p1w * map[ w0 ][ v0 ][ u1 ][AtomType];
                e += p0u * p1v * p1w * map[ w0 ][ v0 ][ u1 ][ElecMap];
         
                e += p1u * p0v * p1w * map[ w0 ][ v1 ][ u0 ][ElecMap];
                m += p1u * p0v * p1w * map[ w0 ][ v1 ][ u0 ][AtomType];
         
                m += p1u * p1v * p0w * map[ w1 ][ v0 ][ u0 ][AtomType];
                e += p1u * p1v * p0w * map[ w1 ][ v0 ][ u0 ][ElecMap];
         
                e += p0u * p0v * p1w * map[ w0 ][ v1 ][ u1 ][ElecMap];
                m += p0u * p0v * p1w * map[ w0 ][ v1 ][ u1 ][AtomType];
         
                m += p1u * p0v * p0w * map[ w1 ][ v1 ][ u0 ][AtomType];
                e += p1u * p0v * p0w * map[ w1 ][ v1 ][ u0 ][ElecMap];
         
                e += p0u * p1v * p0w * map[ w1 ][ v0 ][ u1 ][ElecMap];
                m += p0u * p1v * p0w * map[ w1 ][ v0 ][ u1 ][AtomType];
         
                m += p0u * p0v * p0w * map[ w1 ][ v1 ][ u1 ][AtomType];
                e += p0u * p0v * p0w * map[ w1 ][ v1 ][ u1 ][ElecMap];

                etotal += (elec[i] = e * charge[i]) + (emap[i] = m); 

    #endif /* not MINPOINT */

            } /* inside grid */
        } else { 
            /* if (ignore_inter[i]) */
            elec[i] = 0.0;
            emap[i] = 0.0;
        }
    }/*for  0 <= i < total_atoms*/

    return( (FloatOrDouble)etotal );
}

/* } */

//#------------------------------------------------------------------

/* { */

FloatOrDouble quicktrilinterp4(    CONST_FLOAT tcoord[MAX_ATOMS][SPACE], 
                          CONST_FLOAT charge[MAX_ATOMS], 
                          CONST_INT   type[MAX_ATOMS], 
                          CONST_INT   total_atoms, 
                          CONST_FLOAT map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS], 
                          CONST_FLOAT inv_spacing, 
                          CONST_FLOAT xlo,
                          CONST_FLOAT ylo,
                          CONST_FLOAT zlo,
                          int ignore_inter[MAX_ATOMS])

{
    register double etotal;
    register double u,   v,   w;
    register double p0u, p0v, p0w;
    register double p1u, p1v, p1w;
    register int AtomType;
    register int u0,  v0,  w0;
    register int u1,  v1,  w1;
    register int i;                /* i-th atom */

#ifdef MINPOINT
    register int x,y,z;
#else
    register double e, m; 
#endif

    etotal = 0.;

    for (i=0; i<total_atoms; i++) {

        if (!ignore_inter[i]) {

            AtomType = type[i];

            u1  = (u0 = (int) (u = (tcoord[i][X]-xlo) * inv_spacing)) + 1;
            p1u = 1. - (p0u = u - (double) u0);

            v1  = (v0 = (int) (v = (tcoord[i][Y]-ylo) * inv_spacing)) + 1;
            p1v = 1. - (p0v = v - (double) v0);

            w1  = (w0 = (int) (w = (tcoord[i][Z]-zlo) * inv_spacing)) + 1;
            p1w = 1. - (p0w = w - (double) w0);

    #ifdef MINPOINT
            x = (p0u < p1u)? u0 : u1;
            y = (p0v < p1v)? v0 : v1;
            z = (p0w < p1w)? w0 : w1;

            etotal += map[z][y][x][ElecMap] * charge[i] + map[z][y][x][AtomType]; 
    #else
            e = m = 0.;

            e += p1u * p1v * p1w * map[ w0 ][ v0 ][ u0 ][ElecMap];
            m += p1u * p1v * p1w * map[ w0 ][ v0 ][ u0 ][AtomType];

            m += p0u * p1v * p1w * map[ w0 ][ v0 ][ u1 ][AtomType];
            e += p0u * p1v * p1w * map[ w0 ][ v0 ][ u1 ][ElecMap];

            e += p1u * p0v * p1w * map[ w0 ][ v1 ][ u0 ][ElecMap];
            m += p1u * p0v * p1w * map[ w0 ][ v1 ][ u0 ][AtomType];

            m += p1u * p1v * p0w * map[ w1 ][ v0 ][ u0 ][AtomType];
            e += p1u * p1v * p0w * map[ w1 ][ v0 ][ u0 ][ElecMap];

            e += p0u * p0v * p1w * map[ w0 ][ v1 ][ u1 ][ElecMap];
            m += p0u * p0v * p1w * map[ w0 ][ v1 ][ u1 ][AtomType];

            m += p1u * p0v * p0w * map[ w1 ][ v1 ][ u0 ][AtomType];
            e += p1u * p0v * p0w * map[ w1 ][ v1 ][ u0 ][ElecMap];

            e += p0u * p1v * p0w * map[ w1 ][ v0 ][ u1 ][ElecMap];
            m += p0u * p1v * p0w * map[ w1 ][ v0 ][ u1 ][AtomType];

            m += p0u * p0v * p0w * map[ w1 ][ v1 ][ u1 ][AtomType];
            e += p0u * p0v * p0w * map[ w1 ][ v1 ][ u1 ][ElecMap];

            etotal += e * charge[i] + m; 

    #endif /* not MINPOINT */
        }

    }/*for  0 <= i < total_atoms*/

    return( (FloatOrDouble)etotal );
}
/* } */

/* EOF */
