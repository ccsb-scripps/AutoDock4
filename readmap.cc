/*

 $Id: readmap.cc,v 1.2 2003/02/26 01:35:43 garrett Exp $

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

char mapf2c(FloatOrDouble);

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
void readmap( Boole *P_B_HaveMap, 
             int *P_Imap, 
             int *P_NumAtmMaps, 
             FloatOrDouble *P_ExtSpacing, 
             char AtmTypStr[ATOM_MAPS], 
             char ExtFldFileName[MAX_CHARS],
             int ExtGridPts1[SPACE],
             int ExtGridPts[SPACE],
             Clock jobStart,
             char line[LINE_LEN],
             char ExtMacromolFileName[MAX_CHARS],
             FloatOrDouble map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS],
             FloatOrDouble MapCenter[SPACE],
             FloatOrDouble MapMax[MAX_MAPS],
             FloatOrDouble MapMin[MAX_MAPS],
             struct tms tmsJobStart,
             Boole B_charMap )

{
    FILE *mapFilePtr;

    char FileName[MAX_CHARS];
    char FldFileName[MAX_CHARS];
    char GpfName[MAX_CHARS];
    char ExtGpfName[MAX_CHARS];
    char message[LINE_LEN];
    char mmFileName[MAX_CHARS];
    char xyz_str[4];
    char C_mapValue;
    char mapline[LINE_LEN];
    char inputline[LINE_LEN];

    FloatOrDouble cen[SPACE];
    FloatOrDouble spacing = 0.;

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
        pr( logFile, "Opened Grid Map %d (%c):\t\t\t\t%s\n", 
            (*P_Imap)+1, ((*P_Imap)==(*P_NumAtmMaps))?'e':AtmTypStr[*P_Imap], FileName );
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
            if (!ignore_errors) {
                check_header_line(GpfName, ExtGpfName );
            } /* endif */
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
            #ifdef USE_DOUBLE
                (void) sscanf(inputline,"%*s %lf", &spacing);
            #else
                (void) sscanf(inputline,"%*s %f", &spacing);
            #endif
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
            #ifdef USE_DOUBLE
                (void) sscanf(inputline,"%*s %lf %lf %lf", &cen[X], &cen[Y], &cen[Z]);
            #else
                (void) sscanf(inputline,"%*s %f %f %f", &cen[X], &cen[Y], &cen[Z]);
            #endif
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
    MapMax[*P_Imap] = -BIG;
    MapMin[*P_Imap] =  BIG;
    nvExpected = ExtGridPts1[X] * ExtGridPts1[Y] * ExtGridPts1[Z];
    nv = 0;
    pr( logFile, "Number of grid points expected in  x-dimension:  %d\n", ExtGridPts1[X] );
    pr( logFile, "Number of grid points expected in  y-dimension:  %d\n", ExtGridPts1[Y] );
    pr( logFile, "Number of grid points expected in  z-dimension:  %d\n", ExtGridPts1[Z] );
    pr( logFile, "Looking for %d energies from Grid Map %d... \n", nvExpected, (*P_Imap)+1 );
    flushLog;
    loadStart = times( &tms_loadStart );
    for ( k = 0;  k < ExtGridPts1[Z];  k++) {
        for ( j = 0;  j < ExtGridPts1[Y];  j++) {
            for ( i = 0;  i < ExtGridPts1[X];  i++) {
                if (B_charMap) {
                    if (fgets(mapline, LINE_LEN, mapFilePtr) != NULL) { /*new*/
                        (void) sscanf( mapline,  "%c",  &C_mapValue );
                        map[k][j][i][*P_Imap] = mapc2f(C_mapValue);
                        nv++;
                    }
                } else {
                    if (fgets( mapline, LINE_LEN, mapFilePtr) != NULL) { /*new*/
                        #ifdef USE_DOUBLE
                            (void) sscanf( mapline,  "%lf",  &map[k][j][i][*P_Imap] );
                        #else
                            (void) sscanf( mapline,  "%f",  &map[k][j][i][*P_Imap] );
                        #endif
                        nv++;
                    }
                }
                MapMax[*P_Imap] = max( MapMax[*P_Imap], map[k][j][i][*P_Imap] );
                MapMin[*P_Imap] = min( MapMin[*P_Imap], map[k][j][i][*P_Imap] );
            }
            /* nv += ExtGridPts1[X]; */ /*new*/
        }
    }
    pr( logFile, "Closing file.\n" );
    fclose( mapFilePtr );
    pr( logFile, "%d energies found for map %d\n", nv, (*P_Imap)+1 );
    pr( logFile, "Minimum energy = %.2f,  maximum energy = %.2f\n\n", MapMin[*P_Imap], MapMax[*P_Imap] );
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

    ++(*P_Imap);

    flushLog;
}
/* EOF */
