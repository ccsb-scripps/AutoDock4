/* readfield.cc */


    #include <stdlib.h>
    #include <string.h>
    #include <stdio.h>
    #include <time.h>
    #include <sys/times.h>
    #include <sys/types.h>
    #include "readfield.h"


extern FILE *logFile;

void readfield( float *P_F_invSpacing, 
	     float *P_F_spacing,
	     char FN_gdfld[MAX_CHARS],
	     char FN_gpf[MAX_CHARS],
	     int I_gridpts1[SPACE],
	     int I_gridpts[SPACE],
	     float *P_F_xhi,
	     float *P_F_yhi,
	     float *P_F_zhi,
	     Clock jobStart,
	     char line[LINE_LEN],
	     float *P_F_xlo,
	     float *P_F_ylo,
	     float *P_F_zlo,
	     char FN_macromol[MAX_CHARS],
	     float F_mapCenter[SPACE],
	     struct tms tms_jobStart )

{
    FILE *fldFile;
    char rec9[9], inputline[LINE_LEN];
    float F_localSpacing;
    register int RI_XYZ = 0;

    /*
    ** GRID_DATA_FILE
    ** Read the (AVS-format) grid data file, .fld
    ** _____________________________________________________________
    **/
    (void) sscanf( line, "%*s %s", FN_gdfld);
     
    if ( openFile( FN_gdfld, "r", &fldFile, jobStart,tms_jobStart,TRUE )) {
	pr( logFile, "Opening Grid Map Dimensions file:\t\t%s\n\n", FN_gdfld);
    }

    /*
    ** Skip over the AVS-readable .fld  header comments, until '#SPACING'...
    */
    while( fgets(line, LINE_LEN, fldFile) != NULL ) {
	(void) sscanf(line,"%s", rec9);
	if (equal(rec9,"#SPACING", 8)) {
	    (void) sscanf(line,"%*s %f", &F_localSpacing);
	    *P_F_spacing = F_localSpacing;
	    break;
	}
    } /* endwhile */
    *P_F_invSpacing = 1. / (*P_F_spacing);
    pr( logFile, "Grid Point Spacing =\t\t\t\t%.3f Angstroms\n\n", *P_F_spacing);
    /*
    ** #NELEMENTS 
    */
    (void) fgets(inputline, LINE_LEN, fldFile);
    (void) sscanf(inputline,"%*s %d %d %d", &I_gridpts[X], &I_gridpts[Y], &I_gridpts[Z]);
    pr( logFile, "Even Number of User-specified Grid Points =\t%d x-points\n\t\t\t\t\t\t%d y-points\n\t\t\t\t\t\t%d z-points\n\n", I_gridpts[X],I_gridpts[Y],I_gridpts[Z]);
    for (RI_XYZ = 0;  RI_XYZ < SPACE;  RI_XYZ++) {
	I_gridpts1[RI_XYZ] = I_gridpts[RI_XYZ] + 1;
    } /* RI_XYZ */
    pr( logFile, "Adding the Central Grid Point makes:\t\t%d x-points\n\t\t\t\t\t\t%d y-points\n\t\t\t\t\t\t%d z-points\n\n", I_gridpts1[X], I_gridpts1[Y], I_gridpts1[Z]);
    if ( (I_gridpts[X] <= 0)||(I_gridpts[Y] <= 0)||(I_gridpts[Z] <= 0) ) {
	stop("insufficient grid points." );
	exit( -1 );
    } else if ((I_gridpts[X] > MAX_GRID_PTS)||(I_gridpts[Y] > MAX_GRID_PTS)||(I_gridpts[Z] > MAX_GRID_PTS)) {
	stop("too many grid points." );
	exit( -1 );
    }

    /*
    ** #CENTER 
    */
    (void) fgets(inputline, LINE_LEN, fldFile);
    (void) sscanf(inputline,"%*s %f %f %f", &F_mapCenter[X], &F_mapCenter[Y], &F_mapCenter[Z]);
    pr( logFile, "Coordinates of Central Grid Point of Maps =\t(%.3f, %.3f, %.3f)\n\n", F_mapCenter[X],  F_mapCenter[Y],  F_mapCenter[Z]);

    /*
    ** #MACROMOLECULE 
    */
    (void) fgets(inputline, LINE_LEN, fldFile);
    (void) sscanf(inputline,"%*s %s", FN_macromol);
    pr( logFile, "Macromolecule file used to create Grid Maps =\t%s\n\n", FN_macromol);

    /*
    ** #GRID_PARAMETER_FILE 
    */
    (void) fgets(inputline, LINE_LEN, fldFile);
    (void) sscanf(inputline,"%*s %s", FN_gpf);
    pr( logFile, "Grid Parameter file used to create Grid Maps =\t%s\n\n", FN_gpf);
    /*
    ** Close Grid-dimensions data file 
    */
    fclose(fldFile);

    /*
    ** Determine the dimensions of the grids,
    */
    *P_F_xlo = F_mapCenter[X] - (float)(I_gridpts[X]/2) * (*P_F_spacing);
    *P_F_ylo = F_mapCenter[Y] - (float)(I_gridpts[Y]/2) * (*P_F_spacing);
    *P_F_zlo = F_mapCenter[Z] - (float)(I_gridpts[Z]/2) * (*P_F_spacing);

    *P_F_xhi = F_mapCenter[X] + (float)(I_gridpts[X]/2) * (*P_F_spacing);
    *P_F_yhi = F_mapCenter[Y] + (float)(I_gridpts[Y]/2) * (*P_F_spacing);
    *P_F_zhi = F_mapCenter[Z] + (float)(I_gridpts[Z]/2) * (*P_F_spacing);

    pr( logFile, "Minimum coordinates in grid = (%.3f, %.3f, %.3f)\n",   *P_F_xlo, *P_F_ylo, *P_F_zlo);
    pr( logFile, "Maximum coordinates in grid = (%.3f, %.3f, %.3f)\n\n", *P_F_xhi, *P_F_yhi, *P_F_zhi);
    fflush( logFile );
}
/*
** EOF
*/
