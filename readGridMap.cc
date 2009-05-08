/*

 $Id: readGridMap.cc,v 1.8 2009/05/08 23:02:16 rhuey Exp $

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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <sys/times.h>
#include <time.h>
#include "readmap.h"

extern char dock_param_fn[];
extern char *programname;
extern int ignore_errors;
extern int ElecMap;
extern FILE *logFile;
extern int debug;

char mapf2c(Real);

void readmap( Boole *P_B_HaveMap, 
             int *P_imap, 
             int *num_atom_types, 
             Real *P_ExtSpacing, 
             char ligand_atom_types[MAX_MAPS][3],
             char *ExtFldFileName,
             int ExtGridPts1[SPACE],
             int ExtGridPts[SPACE],
             Clock jobStart,
             char line[LINE_LEN],
             char *ExtMacromolFileName,
                #include "map_declare.h"
             Real MapCenter[SPACE],
             Real MapMax[MAX_MAPS],
             Real MapMin[MAX_MAPS],
             struct tms tmsJobStart,
             Boole B_charMap,
             int outlev,
             GridMap grid_map)

{
    FILE *mapFilePtr;

    char FileName[PATH_MAX];
    char FldFileName[PATH_MAX];
    char GpfName[PATH_MAX];
    char ExtGpfName[PATH_MAX];
    char message[LINE_LEN];
    char mmFileName[PATH_MAX];
    char xyz_str[4];
    char C_mapValue;
    char mapline[LINE_LEN];
    char inputline[LINE_LEN];
    char atom_type_name[MAX_CHARS];
    char map_type = '?';

    Real cen[SPACE];
    Real spacing = 0.;

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

    strcpy( xyz_str, "xyz\0" );


    /*
    \  ATOMIC AFFINITY or ELECTROSTATIC GRID MAP
     \  Read in active site grid map...
      \____________________________________________________________
     */

    (void) sscanf( line, "%*s %s", FileName );
    if ( openFile( FileName, "r", &mapFilePtr, jobStart,tmsJobStart,TRUE )) {
        *P_B_HaveMap = TRUE;
        if (debug > 0) {
            for (i=0; i < (*num_atom_types); i++) {
                (void) fprintf(logFile, "ligand_atom_types[%d] = \"%s\"\n", i, ligand_atom_types[i] );
            }
        }
        if ((*P_imap) == (*num_atom_types)) {
            strcpy(atom_type_name, "e\0");
            map_type = 'e';
        } else if ( (*P_imap) == ((*num_atom_types)+1)) {
            strcpy(atom_type_name, "d\0");
            map_type = 'd';
        } else {
            strcpy(atom_type_name, ligand_atom_types[*P_imap]);
        }
        pr( logFile, "Opened Grid Map %d (%s):\t\t\t\t%s\n", (*P_imap)+1, atom_type_name, FileName );
        if (!ignore_errors) {
            pr( logFile, "Checking header information.\n" );
        }
         /*
         \ Check header lines of grid map... 
         /
         \ :Line 1  GRID_PARAMETER_FILE 
        */
        if (fgets(inputline, LINE_LEN, mapFilePtr) == NULL) {
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
        if (fgets(inputline, LINE_LEN, mapFilePtr) == NULL) {
            warn_bad_file( FileName,"Could not read \".fld\" GRID_DATA_FILE line." );
        } else {
            (void) sscanf(inputline, "%*s %s", FldFileName);
            if (!ignore_errors) {
                check_header_line( FldFileName, ExtFldFileName );
            } /* endif */
        } /* endif */
         /*
         \ :Line 3  MACROMOLECULE 
        */
        if (fgets(inputline, LINE_LEN, mapFilePtr) == NULL) {
            warn_bad_file( FileName,"Could not read MACROMOLECULE line." );
        } else {
            (void) sscanf(inputline,"%*s %s", mmFileName);
            check_header_line( mmFileName, ExtMacromolFileName );
        } /* endif */
         /*
         \ :Line 4  SPACING 
        */
        if (fgets(inputline, LINE_LEN, mapFilePtr) == NULL) {
            warn_bad_file( FileName,"Could not read SPACING line." );
        } else {
            (void) sscanf(inputline,"%*s " FDFMT, &spacing);
            check_header_float(spacing, *P_ExtSpacing, "grid point spacing", FileName );
        } /* endif */
         /*
         \ :Line 5  NELEMENTS 
        */
        if (fgets(inputline, LINE_LEN, mapFilePtr) == NULL) {
            warn_bad_file( FileName,"Could not read NELEMENTS line." );
        } else {
            (void) sscanf(inputline,"%*s %d %d %d", &nel[X], &nel[Y], &nel[Z]);
            for (xyz = 0;  xyz < SPACE;  xyz++) {
                check_header_int( nel[xyz], ExtGridPts[xyz], xyz_str[xyz], FileName );
            } /* xyz */
        } /* endif */
         /* 
         \ :Line 6  CENTER
        */
        if (fgets(inputline, LINE_LEN, mapFilePtr) == NULL) {
            warn_bad_file( FileName,"Could not read CENTER line." );
        } else {
            (void) sscanf(inputline,"%*s " FDFMT3, &cen[X], &cen[Y], &cen[Z]);
            for (xyz = 0;  xyz < SPACE;  xyz++) {
                check_header_float(cen[xyz], MapCenter[xyz], "grid-map center", FileName );
            } /* xyz */
        } /* endif */
    } /* endif */
    flushLog;
    /*
    \   Now find the extrema of the grid-map energies,
     \  While reading in the values...
      \____________________________________________________________
     */
    MapMax[*P_imap] = -BIG;
    MapMin[*P_imap] =  BIG;
    nvExpected = ExtGridPts1[X] * ExtGridPts1[Y] * ExtGridPts1[Z];
    nv = 0;
    pr( logFile, "Number of grid points expected in  x-dimension:  %d\n", ExtGridPts1[X] );
    pr( logFile, "Number of grid points expected in  y-dimension:  %d\n", ExtGridPts1[Y] );
    pr( logFile, "Number of grid points expected in  z-dimension:  %d\n", ExtGridPts1[Z] );
    pr( logFile, "Looking for %d energies from Grid Map %d... \n", nvExpected, (*P_imap)+1 );
    flushLog;
    loadStart = times( &tms_loadStart );
    for ( k = 0;  k < ExtGridPts1[Z];  k++) {
        for ( j = 0;  j < ExtGridPts1[Y];  j++) {
            for ( i = 0;  i < ExtGridPts1[X];  i++) {
                if (B_charMap) {
                    if (fgets(mapline, LINE_LEN, mapFilePtr) != NULL) { /*new*/
                        if (sscanf( mapline,  "%c",  &C_mapValue ) != 1) continue;
                        map[k][j][i][*P_imap] = mapc2f(C_mapValue);
                        nv++;
                    }
                } else {
                    if (fgets( mapline, LINE_LEN, mapFilePtr) != NULL) { /*new*/
                        double v;
                        if (sscanf( mapline,  "%lf",  &v ) != 1) continue;
                        map[k][j][i][*P_imap]  = v;
                        nv++;
                    }
                }
                MapMax[*P_imap] = max( MapMax[*P_imap], map[k][j][i][*P_imap] );
                MapMin[*P_imap] = min( MapMin[*P_imap], map[k][j][i][*P_imap] );
            }
        }
    }
    pr( logFile, "Closing file.\n" );
    fclose( mapFilePtr );
    pr( logFile, "%d energies found for map %d\n", nv, (*P_imap)+1 );
    if (map_type == 'e') {
        pr( logFile, "Minimum electrostatic potential = %.2f,  maximum electrostatic potential = %.2f\n\n", MapMin[*P_imap], MapMax[*P_imap] );
    } else {
        pr( logFile, "Minimum energy = %.2f,  maximum energy = %.2f\n\n", MapMin[*P_imap], MapMax[*P_imap] );
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

        exit(-1);
    } /* END PROGRAM */

    ++(*P_imap);

    flushLog;
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
