/* trilinterp.cc */

#include <math.h>
#include "trilinterp.h"


/* linear interpolation from l (when a=0) to h (when a=1)*/
/* (equal to (a*h)+((1-a)*l) )*/
#define LERP(a,l,h)	((l)+(((h)-(l))*(a)))

extern int ElecMap;


float trilinterp( CONST_FLOAT tcoord[MAX_ATOMS][SPACE], 
		  CONST_FLOAT charge[MAX_ATOMS], 
		  CONST_INT   type[MAX_ATOMS], 
		  CONST_INT   total_atoms, 
		  CONST_FLOAT map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS],
		  CONST_FLOAT inv_spacing, 
		  float elec[MAX_ATOMS], 
		  float emap[MAX_ATOMS], 
		  CONST_FLOAT xlo,	/**/
		  CONST_FLOAT ylo,	/*   float lo[SPACE] ) SLOWER */
		  CONST_FLOAT zlo )	/**/

/*
** float tcoord[MAX_ATOMS][SPACE];	temporary coordinates
** float charge[MAX_ATOMS];		partial atomic charges
** int   type[MAX_ATOMS];		atom type of each atom
** int   total_atoms;			number of atoms
** float map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS];  
** 					intermolecular interaction energies
** float inv_spacing;  			= 1/(grid point spacing, in Angstroms) 
** float elec[MAX_ATOMS];  		electrostatic energies, atom by atom
** float emap[MAX_ATOMS];  		intermolecular energies
** float xlo,ylo,zlo;			minimum coordinates in x,y,z
*/

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
    double	 emaptotal, electotal;
    double	 u,   v,   w;
    double	 p0u, p0v, p0w;
    double	 p1u, p1v, p1w;

    int		 AtomType;	/* atom type */
    int		 u0,  v0,  w0;
    int		 u1,  v1,  w1;

    register int i;		/* i-th atom */

#ifdef MINPOINT
    int		 x,y,z;						    /*MINPOINT*/
#else
    double 	 e, m; 
#endif

    emaptotal = electotal = 0.;

    for (i=0; i<total_atoms; i++) {

        AtomType = type[i];

        u1  = (u0 = (int) (u = (tcoord[i][X]-xlo) * inv_spacing)) + 1;
        p1u = 1. - (p0u = u - (double) u0);

        v1  = (v0 = (int) (v = (tcoord[i][Y]-ylo) * inv_spacing)) + 1;
        p1v = 1. - (p0v = v - (double) v0);

        w1  = (w0 = (int) (w = (tcoord[i][Z]-zlo) * inv_spacing)) + 1;
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

    return( (float)(electotal + emaptotal) );
}
/* End of Function */

/*----------------------------------------------------------------------------*/
/* quicktrilinterp.c */

float quicktrilinterp(	CONST_FLOAT tcoord[MAX_ATOMS][SPACE], 
		  	CONST_FLOAT charge[MAX_ATOMS], 
		  	CONST_INT   type[MAX_ATOMS], 
		  	CONST_INT   total_atoms, 
		  	CONST_FLOAT map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS], 
		  	CONST_FLOAT inv_spacing, 
		  	CONST_FLOAT xlo,
		  	CONST_FLOAT ylo,
		  	CONST_FLOAT zlo )

{
    double	 etotal;
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

    }/*for  0 <= i < total_atoms*/

    return( (float)etotal );
}
/* EOF */

/*----------------------------------------------------------------------------*/

float outsidetrilinterp(CONST_FLOAT tcoord[MAX_ATOMS][SPACE], 
		  	CONST_FLOAT charge[MAX_ATOMS], 
		  	CONST_INT   type[MAX_ATOMS], 
		  	CONST_INT   total_atoms, 
		  	CONST_FLOAT map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS], 
		  	CONST_FLOAT inv_spacing, 
		        // float elec[MAX_ATOMS],
			      // float emap[MAX_ATOMS],
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
    double	 x, y, z;
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
     
	    // etotal += (elec[i] = map[z][y][x][ElecMap] * charge[i]) + 
		      // (emap[i] = map[z][y][x][AtomType]); 
	    etotal += (map[z][y][x][ElecMap] * charge[i]) + 
		      (map[i] = map[z][y][x][AtomType]); 

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

	    // etotal += (elec[i] = e * charge[i]) + (emap[i] = m); 
	    etotal += (e * charge[i]) + (m); 

#endif /* not MINPOINT */

	} /* inside grid */
    }/*for  0 <= i < total_atoms*/

    return( (float)etotal );
}

// added 15-jan-2001

/*----------------------------------------------------------------------------*/

float template_trilinterp( CONST_FLOAT tcoord[MAX_ATOMS][SPACE], 
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

    return((float)sqrt(etotal / (float)total_atoms ));
} // float template_trilinterp( CONST_FLOAT tcoord[MAX_ATOMS][SPACE], 

/*----------------------------------------------------------------------------*/

float outside_templ_trilinterp(CONST_FLOAT tcoord[MAX_ATOMS][SPACE], 
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

    return( (float)sqrt(etotal / (float)total_atoms));
}

/*----------------------------------------------------------------------------*/

float byatom_template_trilinterp( CONST_FLOAT tcoord[MAX_ATOMS][SPACE], 
		  	                      CONST_FLOAT charge[MAX_ATOMS], 
		  	                      CONST_INT   type[MAX_ATOMS], 
		  	                      CONST_INT   total_atoms, 
		  	                      CONST_FLOAT map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS], 
		  	                      CONST_FLOAT inv_spacing, 
		                          float elec[MAX_ATOMS], 
		                          float emap[MAX_ATOMS], 
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

    return( (float)sqrt(etotal / (float)total_atoms));
}/* EOF */
