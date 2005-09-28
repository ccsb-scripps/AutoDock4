/*

 $Id: readmap.cc,v 1.4 2005/09/28 22:54:21 garrett Exp $

*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/* readmap.cc */

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

char mapf2c(FloatOrDouble);

void readmap( char           line[LINE_LEN],
              int            outlev,
 
              Clock          jobStart,
              struct tms     tmsJobStart,
        
              Boole          B_charMap,

              Boole          *P_B_HaveMap, 
              int            *P_imap, 
 
              GridMapSetInfo *info,
              FloatOrDouble map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS]
              // double *maps 
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
    char map_type = '?';

    FloatOrDouble cen[SPACE];
    FloatOrDouble spacing = 0.;
    double max[MAX_MAPS];
    double min[MAX_MAPS];

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

    //maps->atom_type = *P_imap;

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
        if ((*P_imap) == info->num_atom_types) {
            strcpy(atom_type_name, "e\0");
            map_type = 'e';
        } else if ( (*P_imap) == (info->num_atom_types + 1)) {
            strcpy(atom_type_name, "d\0");
            map_type = 'd';
        } else {
            strcpy(atom_type_name, info->atom_type_name[*P_imap]);
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
            #ifdef USE_DOUBLE
                (void) sscanf(inputline,"%*s %lf", &spacing);
            #else
                (void) sscanf(inputline,"%*s %f", &spacing);
            #endif
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
            #ifdef USE_DOUBLE
                (void) sscanf(inputline,"%*s %lf %lf %lf", &cen[X], &cen[Y], &cen[Z]);
            #else
                (void) sscanf(inputline,"%*s %f %f %f", &cen[X], &cen[Y], &cen[Z]);
            #endif
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
    max[*P_imap] = -BIG;
    min[*P_imap] =  BIG;
    nvExpected = info->num_points1[X] * info->num_points1[Y] * info->num_points1[Z];
    nv = 0;
    pr( logFile, "Number of grid points expected in  x-dimension:  %d\n", info->num_points1[X] );
    pr( logFile, "Number of grid points expected in  y-dimension:  %d\n", info->num_points1[Y] );
    pr( logFile, "Number of grid points expected in  z-dimension:  %d\n", info->num_points1[Z] );
    pr( logFile, "Looking for %d energies from Grid Map %d... \n", nvExpected, (*P_imap)+1 );
    flushLog;
    loadStart = times( &tms_loadStart );
    for ( k = 0;  k < info->num_points1[Z];  k++) {
        for ( j = 0;  j < info->num_points1[Y];  j++) {
            for ( i = 0;  i < info->num_points1[X];  i++) {
                if (B_charMap) {
                    if (fgets(map_line, LINE_LEN, map_file) != NULL) {
                        (void) sscanf( map_line,  "%c",  &C_mapValue );
                        // maps->map[k][j][i][*P_imap] = mapc2f(C_mapValue);
                        map[k][j][i][*P_imap] = mapc2f(C_mapValue);
                        nv++;
                    }
                } else {
                    if (fgets( map_line, LINE_LEN, map_file) != NULL) {
                        // (void) sscanf( map_line,  "%lf",  &maps->map[k][j][i][*P_imap] );
#ifdef USE_DOUBLE
                        (void) sscanf( map_line,  "%lf",  &map[k][j][i][*P_imap] );
#else
                        (void) sscanf( map_line,  "%f",  &map[k][j][i][*P_imap] );
#endif
                        nv++;
                    }
                }
                //max[*P_imap] = max( max[*P_imap], maps->map[k][j][i][*P_imap] );
                //min[*P_imap] = min( min[*P_imap], maps->map[k][j][i][*P_imap] );
                max[*P_imap] = max( max[*P_imap], map[k][j][i][*P_imap] );
                min[*P_imap] = min( min[*P_imap], map[k][j][i][*P_imap] );
            }
            /* nv += info->num_points1[X]; */
        }
    }
    pr( logFile, "Closing file.\n" );
    fclose( map_file );
    pr( logFile, "%d energies found for map %d\n", nv, (*P_imap)+1 );
    if (map_type == 'e') {
        pr( logFile, "Minimum electrostatic potential = %.2f,  maximum electrostatic potential = %.2f\n\n", min[*P_imap], max[*P_imap] );
    } else {
        pr( logFile, "Minimum energy = %.2f,  maximum energy = %.2f\n\n", min[*P_imap], max[*P_imap] );
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

FloatOrDouble mapc2f(char numin)
{
    FloatOrDouble numout;
    if (numin == 0) {
        numout = 0.;
    } else if (numin > 0) {
        numout = numin * 10.;
    } else {
        numout = numin /10.;
    }
    return numout;
}

void scale_map(
        double weight,
        FloatOrDouble map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS]
        )
{
    
}


/*
    char mapf2c(FloatOrDouble numin)
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
