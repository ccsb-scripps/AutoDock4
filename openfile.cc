/*

 $Id: openfile.cc,v 1.2.2.1 2005/03/01 00:01:29 gillet Exp $

*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/* openfile.cc */

#include <stdio.h>
#include <stdlib.h>
/* the BOINC API header file */
#ifdef BOINC
#include "diagnostics.h"
#include "boinc_api.h" 
#include "filesys.h" 		// boinc_fopen(), etc... */
#endif
#include "openfile.h"


extern char *programname;
extern FILE *logFile;
/*----------------------------------------------------------------------------*/
/* fopen rewrite to either use BOINC api or normal system call */
FILE *ad_fopen(const char *path, const char *mode)
{
  FILE *filep;
#ifdef BOINC
  int rc;
  char resolved_name[512];
  rc = boinc_resolve_filename(path, resolved_name, sizeof(resolved_name));
  if (rc){
      fprintf(stderr, "BOINC_ERROR: cannot open filename.%s\n",path);
      boinc_finish(rc);    /* back to BOINC core */
    }
    // Then open the file with boinc_fopen() not just fopen()
    filep = boinc_fopen(resolved_name, mode);
#else
    filep = fopen(path,mode);
#endif
    return filep;
}

/*----------------------------------------------------------------------------*/
int openfile( char filename[MAX_CHARS],
	      char mode[],
	      FILE **fp )

{
	if ( (*fp = ad_fopen(filename, mode)) == NULL ) {
		fprintf(stderr, "\n%s: I'm sorry; I can't find or open \"%s\"\n", programname, filename);
		fprintf(logFile,"\n%s: I'm sorry; I can't find or open \"%s\"\n", programname, filename);
		return( FALSE );
	}
	return( TRUE );
}

/*----------------------------------------------------------------------------*/
int openFile( char       filename[MAX_CHARS],
	      char       mode[],
	      FILE       **fp,
	      Clock      start,
	      struct tms tms_start,
	      Boole      mayExit)

{
    Clock  jobEnd;
    struct tms tms_jobEnd;

    if ( (*fp = ad_fopen(filename, mode)) == NULL ) {
	fprintf(stderr, "\n%s: I'm sorry; I can't find or open \"%s\"\n", programname, filename);
	fprintf(logFile,"\n%s: I'm sorry; I can't find or open \"%s\"\n", programname, filename);

	jobEnd = times( &tms_jobEnd );
	timesys( jobEnd - start, &tms_start, &tms_jobEnd );
	pr_2x( logFile, stderr, UnderLine );

	if (mayExit) {
	    exit(-1); /* END PROGRAM */
	} else {
	    return( FALSE );
	}
    }
    return( TRUE );
}
/*----------------------------------------------------------------------------*/
