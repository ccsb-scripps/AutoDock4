/*

 $Id: readmap.cc,v 1.9 2007/04/27 06:01:51 garrett Exp $

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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <sys/times.h>
#include <time.h>
#include <math.h>
#include "readmap.h"

extern char dock_param_fn[];
extern char *programname;
extern int ignore_errors;
extern int ElecMap;
extern FILE *logFile;
extern int debug;

char mapf2c(Real);

Statistics readmap( char           line[LINE_LEN],
                    int            outlev,
                    Clock          jobStart,
                    struct tms     tmsJobStart,
                    Boole          B_charMap,
                    Boole          *P_B_HaveMap, 
                    int            num_maps, 
                    GridMapSetInfo *info,
                    Real           map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS],
                    // double *maps 
                    char           map_type
                  )

{
    FILE *map_file;

    char FileName[MAX_CHARS];
    char FldFileName[MAX_CHARS];
    char GpfName[MAX_CHARS];
    char ExtGpfName[MAX_CHARS];
    char message[LINE_LEN];
    char mmFileName[MAX_CHARS];
    char xyz_str[4];
    char C_mapValue;
    char map_line[LINE_LEN];
    char inputline[LINE_LEN];
    char atom_type_name[MAX_CHARS];

    Real cen[SPACE];
    Real spacing = 0.;
    double map_max;
    double map_min;
    double map_total;

    int indpf = 0;
    int nel[SPACE];
    int nv=0;
    int nvExpected = 0;

    register int xyz = 0;
    register int i = 0;
    register int j = 0;
    register int k = 0;

    struct tms tms_jobEnd;
    struct tms tms_loadEnd;
    struct tms tms_loadStart;

    Clock jobEnd;
    Clock loadEnd;
    Clock loadStart;

    Statistics map_stats;

    strcpy( xyz_str, "xyz\0" );

    //maps->atom_type = num_maps;

    /*
    \  ATOMIC AFFINITY or ELECTROSTATIC GRID MAP
     \  Read in active site grid map...
      \____________________________________________________________
     */

    (void) sscanf( line, "%*s %s", FileName );
    if ( openFile( FileName, "r", &map_file, jobStart,tmsJobStart,TRUE )) {
        *P_B_HaveMap = TRUE;
        if (debug > 0) {
            for (i=0; i < info->num_atom_types; i++) {
                (void) fprintf(logFile, "info->atom_type_name[%d] = \"%s\"\n", i, info->atom_type_name[i] );
            }
        }

        if (map_type == 'e') {
            strcpy(atom_type_name, "e\0");
        } else if (map_type == 'd') {
            strcpy(atom_type_name, "d\0");
        } else {
            strcpy(atom_type_name, info->atom_type_name[num_maps]);
        }
        pr( logFile, "Opened Grid Map %d (%s):\t\t\t\t%s\n", num_maps+1, atom_type_name, FileName );

        if (!ignore_errors) {
            pr( logFile, "Checking header information.\n" );
        }
         /*
         \ Check header lines of grid map... 
         /
         \ :Line 1  GRID_PARAMETER_FILE 
        */
        if (fgets(inputline, LINE_LEN, map_file) == NULL) {
            warn_bad_file( FileName,"Could not read GRID_PARAMETER_FILE line." );
        } else {
            (void) sscanf(inputline, "%*s %s", GpfName);
            if ( strindex( dock_param_fn, ".dpf" ) == -1) {
                pr_2x( stderr, logFile,"Can't find \".dpf\" in the dock-parameter filename.\n\n" );
                pr_2x( stderr, logFile,"AutoDock needs the extension of the grid parameter file to be \".gpf\"\nand that of the docking parameter file to be \".dpf\".\n\n" );
            } else {
                /*
                \ replace ".dpf" with ".gpf".
                */
                indpf = strindex( dock_param_fn, "dpf" );
                strcpy(ExtGpfName, dock_param_fn);
                ExtGpfName[ indpf ] = 'g';
            }
        } /* endif */
         /*
         \ :Line 2  GRID_DATA_FILE 
        */
        if (fgets(inputline, LINE_LEN, map_file) == NULL) {
            warn_bad_file( FileName,"Could not read \".fld\" GRID_DATA_FILE line." );
        } else {
            (void) sscanf(inputline, "%*s %s", FldFileName);
            if (!ignore_errors) {
                check_header_line( FldFileName, info->FN_gdfld );
            } /* endif */
        } /* endif */
         /*
         \ :Line 3  MACROMOLECULE 
        */
        if (fgets(inputline, LINE_LEN, map_file) == NULL) {
            warn_bad_file( FileName,"Could not read MACROMOLECULE line." );
        } else {
            (void) sscanf(inputline,"%*s %s", mmFileName);
            check_header_line( mmFileName, info->FN_receptor );
        } /* endif */
         /*
         \ :Line 4  SPACING 
        */
        if (fgets(inputline, LINE_LEN, map_file) == NULL) {
            warn_bad_file( FileName,"Could not read SPACING line." );
        } else {
            (void) sscanf(inputline,"%*s " FDFMT, &spacing);
            check_header_float(spacing, info->spacing, "grid point spacing", FileName );
        } /* endif */
         /*
         \ :Line 5  NELEMENTS 
        */
        if (fgets(inputline, LINE_LEN, map_file) == NULL) {
            warn_bad_file( FileName,"Could not read NELEMENTS line." );
        } else {
            (void) sscanf(inputline,"%*s %d %d %d", &nel[X], &nel[Y], &nel[Z]);
            for (xyz = 0;  xyz < SPACE;  xyz++) {
                //maps->num_points[xyz] = nel[xyz];
                //maps->num_points1[xyz] = nel[xyz] + 1;
                check_header_int( nel[xyz], info->num_points[xyz], xyz_str[xyz], FileName );
            } /* xyz */
        } /* endif */
         /* 
         \ :Line 6  CENTER
        */
        if (fgets(inputline, LINE_LEN, map_file) == NULL) {
            warn_bad_file( FileName,"Could not read CENTER line." );
        } else {
            (void) sscanf(inputline,"%*s " FDFMT3, &cen[X], &cen[Y], &cen[Z]);
            for (xyz = 0;  xyz < SPACE;  xyz++) {
                //maps->center[xyz] = cen[xyz];
                check_header_float(cen[xyz], info->center[xyz], "grid-map center", FileName );
            } /* xyz */
        } /* endif */
    } /* endif */
    flushLog;

    /*
    \   Now find the extrema of the grid-map energies,
     \  While reading in the values...
      \____________________________________________________________
     */

    map_max = -BIG;
    map_min =  BIG;
    map_total = 0.;
    nvExpected = info->num_points1[X] * info->num_points1[Y] * info->num_points1[Z];
    nv = 0;

    pr( logFile, "Number of grid points expected in  x-dimension:  %d\n", info->num_points1[X] );
    pr( logFile, "Number of grid points expected in  y-dimension:  %d\n", info->num_points1[Y] );
    pr( logFile, "Number of grid points expected in  z-dimension:  %d\n", info->num_points1[Z] );
    pr( logFile, "Looking for %d energies from Grid Map %d... \n", nvExpected, num_maps+1 );
    flushLog;
    loadStart = times( &tms_loadStart );

    for ( k = 0;  k < info->num_points1[Z];  k++) {
        for ( j = 0;  j < info->num_points1[Y];  j++) {
            for ( i = 0;  i < info->num_points1[X];  i++) {
                if (B_charMap) {
                    if (fgets(map_line, LINE_LEN, map_file) != NULL) { /*new*/
                        if (sscanf( map_line,  "%c",  &C_mapValue ) != 1) continue;
                        map[k][j][i][num_maps] = mapc2f(C_mapValue);
                        nv++;
                    }
                } else {
                    if (fgets( map_line, LINE_LEN, map_file) != NULL) { /*new*/
                        if (sscanf( map_line,  FDFMT,  &map[k][j][i][num_maps] ) != 1) continue;
                        nv++;
                    }
                }
                map_max = max( map_max, map[k][j][i][num_maps] );
                map_min = min( map_min, map[k][j][i][num_maps] );
                map_total += map[k][j][i][num_maps];
            }
        }
    }

    map_stats.number = nv;
    map_stats.minimum = map_min;
    map_stats.maximum = map_max;
    map_stats.mean = map_total / (double) nv;

    /*
    if (map_stats.number > 1) {
        double deviation = 0.;
        double sum_squares = 0.;
        for ( k = 0;  k < info->num_points1[Z];  k++) {
            for ( j = 0;  j < info->num_points1[Y];  j++) {
                for ( i = 0;  i < info->num_points1[X];  i++) {
                    deviation = map[k][j][i][num_maps] - map_stats.mean;
                    sum_squares += deviation * deviation;
                }
            }
        }
        map_stats.standard_deviation = sqrt(sum_squares / (map_stats.number - 1));
    } else {
        map_stats.standard_deviation = 0.;
    }
    */

    pr( logFile, "Closing file.\n" );
    fclose( map_file );
    pr( logFile, "%d energies found for map %d\n", nv, num_maps+1 );
    if (map_type == 'e') {
        pr( logFile, "Minimum electrostatic potential = %.2f,  maximum electrostatic potential = %.2f\n\n", map_min, map_max );
    } else {
        pr( logFile, "Minimum energy = %.2f,  maximum energy = %.2f\n\n", map_min, map_max );
    }
    pr( logFile, "Time taken (s): " );

    loadEnd = times( &tms_loadEnd );
    timesys( loadEnd - loadStart, &tms_loadStart, &tms_loadEnd );

    pr( logFile, "\n" );

    if (nv != nvExpected ) {
        prStr( message, "\n%s: wrong number of values read in. Check grid map!\n\n", programname  );
        pr_2x( stderr, logFile, message );

        jobEnd = times( &tms_jobEnd );
        timesys( jobEnd - jobStart, &tmsJobStart, &tms_jobEnd );
        pr_2x( logFile, stderr, UnderLine );

        /* END PROGRAM */
        exit(-1);
    } 

    flushLog;

    return map_stats;
}

Real mapc2f(char numin)
{
    Real numout;
    if (numin == 0) {
        numout = 0.;
    } else if (numin > 0) {
        numout = numin * 10.;
    } else {
        numout = numin /10.;
    }
    return numout;
}


/*
    char mapf2c(Real numin)
    {
        char numout;
        if (numin == 0.) {
            numout = 0;
        } else if ((-12.8 < numin) && (numin < 0.)) {
            numout = numin * 10.;
        } else if ((0. < numin) && (numin < 1280.)) {
            numout = numin / 10.;
        } else if (numin >= 1280.) {
            numout = 127;
        } else {
            numout = -128;
        }
        return numout;
    }
*/

/* EOF */
