/*

 $Id: autocomm.h,v 1.19 2009/05/08 23:02:10 rhuey Exp $

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

/* autocomm.h */

#ifndef _AUTOCOMM
#define _AUTOCOMM

#include <sys/types.h>
#include <time.h>
/* include stdio to pick up definition of FILENAME_MAX and possibly PATH_MAX */
#include <stdio.h>

/*******************************************************************************
**      Name: autocomm.h                                                      **
**  Function: Defines Constants, common to both AUTOGRID & AUTODOCK...        **
**Copyright (C) 2009 The Scripps Research Institute. All rights reserved.
**----------------------------------------------------------------------------**
**    Author: Garrett Matthew Morris, The Scripps Research Institute          **
**      Date: 02/28/1995                                                      **
**----------------------------------------------------------------------------**
**    Inputs: none                                                            **
**   Returns: nothing                                                         **
**   Globals: all defines                                                     **
**----------------------------------------------------------------------------**
** Modification Record                                                        **
** Date     Inits   Comments                                                  **
** 02/28/95 GMM     This header was added.                                    **
*******************************************************************************/

#ifndef COMMON_STUFF
#define COMMON_STUFF

/*
** Constants,
*/

#define FALSE        0      /* Logical constant                               */
#define TRUE         1      /* Logical constant                               */

#define PI	         3.14159265358979323846   /* Mathematical constant, pi */
#define TWOPI	     6.28318530717958647692
#define HALF_PI      1.57079632679489661923

#define X            0      /* x-coordinate                                   */
#define Y            1      /* y-coordinate                                   */
#define Z            2      /* z-coordinate                                   */
#define XYZ          3      /* Dimensions of Cartesian Space                  */
#define SPACE        3      /* Dimensions of Cartesian Space                  */

#define APPROX_ZERO  1.0E-6 /* To avoid division-by-zero errors...            */
#define BIG          1.0E12 /* Very large constant                            */
#define MAX_CHARS    128    /* Number of characters in atom data & filenames  */
#define MAX_LINES    256    /* Number of lines in parameter file              */
#ifndef PATH_MAX
#define PATH_MAX     FILENAME_MAX
#endif

#ifdef USE_XCODE
#define LINE_LEN     140    /* Line length in characters                      */
#else
#define LINE_LEN     256    /* Line length in characters                      */
#endif

#if defined( USE_XCODE )
/* The stacksize limit within Xcode forces us to use smaller grids */
#define MAX_GRID_PTS 61     	/* Maximum number of grid points in 1 dimension */
#elif defined( __CYGWIN__ ) 
#define MAX_GRID_PTS 64		/* Maximum number of grid points in 1 dimension */
#else
#define MAX_GRID_PTS 128	/* Maximum number of grid points in 1 dimension */
				/* MAX_GRID_PTS 128 causes a SIGSEGV on Cygwin */
#endif

#define	EINTCLAMP    100000. /* Clamp pairwise internal energies (kcal/mol )  */

#define MAX_MAPS_PAD 0       // Use this to pad MAX_MAPS to a power of 2, for presumably-faster memory access
#define NUM_NON_VDW_MAPS 2   // Number of electrostatic and desolvation maps
#define MAX_ATOM_TYPES (16 - NUM_NON_VDW_MAPS)    /* Maximum number of atom types set to keep MAX_MAPS a power of 2 */
#define MAX_MAPS (MAX_ATOM_TYPES + NUM_NON_VDW_MAPS + MAX_MAPS_PAD) /* Maximum number of energy maps        */
                            /* 0,1,2,... are for atomic interactions          */
                            /* last two are for electrostatics and desolvation */

#define VECLENMAX    16     /* For AVS fld files...                           */

// Legacy definitions:
#define COVALENTTYPE 'Z'
#define COVALENTTYPE2 'Y'

#define CARBON		0
#define NITROGEN	1
#define OXYGEN		2
#define SULPHUR		3
#define HYDROGEN	4
#define UNKNOWN		5
#define METAL		6
#define COVALENT 7
#define COVALENT2 8
// end Legacy definitions


#define UnderLine "________________________________________________________________________________\n\n"

/*
** Common Macros...
*/

#define pr              (void) fprintf
#define pr_2x           print_2x
#define prStr           (void) sprintf
#define flushLog        (void) fflush(logFile)

#define dist(x1,y1,z1,x2,y2,z2,r) _dx=((x2)-(x1)),_dy=((y2)-(y1)),_dz=((z2)-(z1)),r=sqrt(_dx*_dx + _dy*_dy + _dz*_dz)

/*
** New types...
*/


#include "typedefs.h"


typedef char Boole;


typedef struct AtomDesc {

	Real crd[XYZ];
	Real q;
	int   type;

	} AtomDesc;


/*
** Note the following differing definitions of "times" and "time":-
**
** Arch. times()				time()
** ----- ----------------------------------	--------------------------
** Sun	 clock_t times(struct tms *buffer);	time_t time(time_t *tloc);
** SGI	 clock_t times(struct tms *buffer);	time_t time(time_t *tloc);
** HP	 clock_t times(struct tms *buffer);	time_t time(time_t *tloc);
** Alpha time_t  times(struct tms *buffer);	time_t time(time_t *tloc);
** ----- ----------------------------------	--------------------------
**	 Clock					time_t
**
** Arch	 srand48()												localtime()
** ----- -------------------------------									-----------------------------------------------
** Sun	 void srand48(long seedval);										struct tm *localtime(const time_t *clock);
** SGI	 void srand48 (long seedval);										struct tm *localtime(const time_t *clock);
** HP	 void srand48(long int seedval);									struct tm *localtime(const time_t *timer);
** Alpha void srand48 (long seed_val);										struct tm *localtime(const time_t *timer );
**
** timesys and timesyshms used to use Clock, should use time_t
**
*/

#ifdef __alpha
#define Clock time_t
#else
#define Clock clock_t
#endif /* #ifdef __alpha */


#endif

/*
 * assert that quaternions are OK
 */
#include <assert.h> // for assert in assertQuatOK
#include <math.h> // for sqrt in assertQuatOK

#define ONE_MINUS_EPSILON 0.999
#define ONE_PLUS_EPSILON 1.001

/*
 * void assertQuatOK( const Quat q )
 * {
 *     register double mag4 = hypotenuse4( q.x, q.y, q.z, q.w );
 *     assert((mag4 > ONE_MINUS_EPSILON) && (mag4 < ONE_PLUS_EPSILON));
 * }
 */
#define assertQuatOK( q ) {fflush(logFile);register double aQOK_mag4 = hypotenuse4( (q).x, (q).y, (q).z, (q).w ); assert((aQOK_mag4 > ONE_MINUS_EPSILON) && (aQOK_mag4 < ONE_PLUS_EPSILON)); }


#endif /*_AUTOCOMM*/

/*
** EOF
*/
