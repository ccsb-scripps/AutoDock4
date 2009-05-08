/*

 $Id: trilinterp.cc,v 1.15 2009/05/08 23:02:18 rhuey Exp $

 AutoDock  

Copyright (C) 2009 The Scripps Research Institute. All rights reserved.

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
#include "trilinterp.h"

/* linear interpolation from l (when a=0) to h (when a=1)*/
/* (equal to (a*h)+((1-a)*l) )*/
#define LERP(a,l,h)	((l)+(((h)-(l))*(a)))

extern int ElecMap;
extern int DesolvMap;

#ifdef DEBUG
#include <stdio.h>
extern FILE *logFile;
#endif

Real trilinterp( 

 const int first_atom, // loop begins at this atom  for (i=first_atom;
 const int last_atom, // loop ends at this atom - 1       i<last_atom; i++)
 const Real tcoord[MAX_ATOMS][SPACE], // temporary coordinates
 const Real charge[MAX_ATOMS], // partial atomic charges
 const Real abs_charge[MAX_ATOMS], // absolute magnitude of partial charges
 const int   type[MAX_ATOMS], // atom type of each atom
 #include "map_declare.h"
 GridMapSetInfo *info, // info->lo[X],info->lo[Y],info->lo[Z],    minimum coordinates in x,y,z
 int some_atoms_outside_grid, // boolean
 int ignore_inter[MAX_ATOMS], // array of booleans, says to ignore computation intermolecular energies per atom
 Real elec[MAX_ATOMS], // set if not NULL - electrostatic energies, atom by atom
 Real emap[MAX_ATOMS],  // set if not NULL - intermolecular energies
 Real *p_elec_total, // set if not NULL - total electrostatic energy
 Real *p_emap_total // set if not NULL - total intermolecular energy
 )

/******************************************************************************/
/*      Name: trilinterp                                                      */
/*  Function: Trilinear interpolation of interaction energies from map[]      */
/*            using the coordinates in tcoord[].                              */
/*Copyright (C) 2009 The Scripps Research Institute. All rights reserved. */
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
    register double elec_total=0, emap_total=0;
    register int i;               /* i-th atom */

    // for (i=0; i<total_atoms; i++) {
    for (i=first_atom; i<last_atom; i++) {
        register double e, m, d; 
        register double u,   v,   w;
        register double p0u, p0v, p0w;
        register double p1u, p1v, p1w;
        register int AtomType;        /* atom type */
        register int u0,  v0,  w0;
        register int u1,  v1,  w1;

        if (ignore_inter[i]) {
            if (elec != NULL) elec[i] = 0;
            if (emap != NULL) emap[i] = 0;
            continue;
        }

        if (some_atoms_outside_grid) {
            register double x,y,z;
            x = tcoord[i][X];
            y = tcoord[i][Y];
            z = tcoord[i][Z];
            if (is_out_grid_info(x,y,z)) {
                register double epenalty;
                x -= info->center[X];
                y -= info->center[Y];
                z -= info->center[Z];
                // sqhypotenuse(x,y,z) is the square of the distance from grid's centre to atom
                epenalty = sqhypotenuse(x,y,z) * ENERGYPENALTY;
                if (elec != NULL) elec[i] = epenalty;
                if (emap != NULL) emap[i] = epenalty;
                elec_total += epenalty;
                emap_total += epenalty;
                continue;
            }
        }

        AtomType = type[i];

        u1  = (u0 = (int) (u = ((double)tcoord[i][X]-(double)info->lo[X]) * (double)info->inv_spacing)) + 1;
        p1u = 1.0L - (p0u = u - (double) u0);

        v1  = (v0 = (int) (v = ((double)tcoord[i][Y]-(double)info->lo[Y]) * (double)info->inv_spacing)) + 1;
        p1v = 1.0L - (p0v = v - (double) v0);

        w1  = (w0 = (int) (w = ((double)tcoord[i][Z]-(double)info->lo[Z]) * (double)info->inv_spacing)) + 1;
        p1w = 1.0L - (p0w = w - (double) w0);

#ifdef MINPOINT
        register int ix,iy,iz;                      /*MINPOINT*/
        ix = (p0u < p1u)? u0 : u1;				    /*MINPOINT*/
        iy = (p0v < p1v)? v0 : v1;				    /*MINPOINT*/
        iz = (p0w < p1w)? w0 : w1;				    /*MINPOINT*/

        e = map[iz][iy][ix][ElecMap];               /*MINPOINT*/
        m = map[iz][iy][ix][AtomType]; 	            /*MINPOINT*/
        d = map[iz][iy][ix][DesolvMap]; 	        /*MINPOINT*/
#else
        e = m = d = 0.0L;

        e += p1u * p1v * p1w * map[ w0 ][ v0 ][ u0 ][ElecMap];
        m += p1u * p1v * p1w * map[ w0 ][ v0 ][ u0 ][AtomType];
        d += p1u * p1v * p1w * map[ w0 ][ v0 ][ u0 ][DesolvMap];

        d += p0u * p1v * p1w * map[ w0 ][ v0 ][ u1 ][DesolvMap];
        m += p0u * p1v * p1w * map[ w0 ][ v0 ][ u1 ][AtomType];
        e += p0u * p1v * p1w * map[ w0 ][ v0 ][ u1 ][ElecMap];

        e += p1u * p0v * p1w * map[ w0 ][ v1 ][ u0 ][ElecMap];
        m += p1u * p0v * p1w * map[ w0 ][ v1 ][ u0 ][AtomType];
        d += p1u * p0v * p1w * map[ w0 ][ v1 ][ u0 ][DesolvMap];

        d += p0u * p0v * p1w * map[ w0 ][ v1 ][ u1 ][DesolvMap];
        m += p0u * p0v * p1w * map[ w0 ][ v1 ][ u1 ][AtomType];
        e += p0u * p0v * p1w * map[ w0 ][ v1 ][ u1 ][ElecMap];

        e += p1u * p1v * p0w * map[ w1 ][ v0 ][ u0 ][ElecMap];
        m += p1u * p1v * p0w * map[ w1 ][ v0 ][ u0 ][AtomType];
        d += p1u * p1v * p0w * map[ w1 ][ v0 ][ u0 ][DesolvMap];

        d += p0u * p1v * p0w * map[ w1 ][ v0 ][ u1 ][DesolvMap];
        m += p0u * p1v * p0w * map[ w1 ][ v0 ][ u1 ][AtomType];
        e += p0u * p1v * p0w * map[ w1 ][ v0 ][ u1 ][ElecMap];

        e += p1u * p0v * p0w * map[ w1 ][ v1 ][ u0 ][ElecMap];
        m += p1u * p0v * p0w * map[ w1 ][ v1 ][ u0 ][AtomType];
        d += p1u * p0v * p0w * map[ w1 ][ v1 ][ u0 ][DesolvMap];

        d += p0u * p0v * p0w * map[ w1 ][ v1 ][ u1 ][DesolvMap];
        m += p0u * p0v * p0w * map[ w1 ][ v1 ][ u1 ][AtomType];
        e += p0u * p0v * p0w * map[ w1 ][ v1 ][ u1 ][ElecMap];
#endif /* not MINPOINT */

        elec_total += e * charge[i];
        emap_total += m + d * abs_charge[i]; 

        if (elec != NULL) elec[i] = e * charge[i];
        if (emap != NULL) emap[i] = m + d * abs_charge[i];

    } // for (i=first_atom; i<last_atom; i++)

    if (p_elec_total != NULL) *p_elec_total = elec_total;
    if (p_emap_total != NULL) *p_emap_total = emap_total;

    return( (Real)elec_total+emap_total );
}

/*----------------------------------------------------------------------------*/

/* EOF */
