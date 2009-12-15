/*

 $Id: setflags.cc,v 1.19 2009/12/15 06:15:03 mp Exp $

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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "setflags.h"
#include "openfile.h"
#include "version.h"
#include "banner.h"
#include "strindex.h"


extern FILE *parFile;
extern FILE *logFile;
extern FILE *stateFile;
extern int  write_stateFile;
extern char *programname;

extern char dock_param_fn[];
extern char AutoDockHelp[];
extern int  debug;
extern int  ignore_errors;
extern int  parse_tors_mode;
extern int  keepresnum;


int setflags( int argc, char ** argv, const char * version_num)

/*
** naming convention: 
** var   => var is an integer variable;
** var => var is a pointer to a pointer to a character variable
*/

/******************************************************************************/
/*      Name: setflags                                                        */
/*  Function: read flags from argv; return argindex of first non arg.   */
/*Copyright (C) 2009 The Scripps Research Institute. All rights reserved. */
/*----------------------------------------------------------------------------*/
/*    Author: Garrett Matthew Morris, TSRI.                                   */
/*            (Adapted from code supplied by Bruce Duncan, TSRI.)             */
/*      Date: 02/02/94                                                        */
/*----------------------------------------------------------------------------*/
/*    Inputs: argc, argv                                                 */
/*   Returns: argindex                                                      */
/*   Globals: *parFile;				                              */
/*            *logFile;					                      */
/*            *programname;						      */
/*            dock_param_fn[];						      */
/*            ignore_errors;						      */
/*            command_mode;						      */
/*----------------------------------------------------------------------------*/
/* Modification Record                                                        */
/* Date     Inits   Comments                                                  */
/* 06/11/92 GMM     Modified for Autodock flags:                              */
/*                  -c = Command mode;                                        */
/*                  -p = Parameter filename;                                  */
/*                  -l = Log filename;                                        */
/*                  -o = Use old PDBq format (q in columns 55-61)             */
/*                  -d = Debug mode;                                          */
/* 01/22/93 GMM     -i = Ignore header-checking                               */
/* 06/14/93 GMM     -t = Parse the PDBq file to check torsions, then stop.    */
/* 02/02/94 GMM     -k = Keep original residue number in output PDB clusters. */
/******************************************************************************/

{
    int argindex;
/*----------------------------------------------------------------------------*/
/* Initialize                                                                 */
/*----------------------------------------------------------------------------*/
    argindex = 1;
    programname = argv[0];
    parFile = stdin;
    logFile = stdout;
    char logFileName[PATH_MAX+2];
    static char * p_logFileName = "stdout"; // change with -l <NAME> or defaults
     // to parFile name with last 3 chars changed from "dpf" to "dlg"
    /*
     * see autoglobal.h for initialization of debug, keepresnum and logicals...
     */
    if (argc==1) { //No arguments provided
        usage(stdout, "AutoDock");
        exit(0);
    }
/*----------------------------------------------------------------------------*/
/* Loop over arguments                                                        */
/*----------------------------------------------------------------------------*/
    while((argc > 1) && (argv[1][0] == '-')){
        if (argv[1][1] == '-') argv[1]++;

        switch(argv[1][1]){
#ifdef FOO
        case 'n':
            ncount = atoi(argv[2]);
            argv++;
            argc--;
            argindex++;
            break;
#endif
        case 'd':
            debug++;
            break;
        case 'u':
        case 'h':
            usage(stdout, "AutoDock");
            exit(0);
            break;
        case 'i':
            ignore_errors = TRUE;
            break;
        case 'k':
            keepresnum--;
            break;
        case 'C':
            //show copyright
            show_copyright(stdout);
            show_warranty(stdout);
            exit(0);
            break;
        case 'c':
            //command_mode removed with 4.1 release spring 2009, mp + rh
            fprintf(stderr, "\n%s: command mode is not supported in this version of autodock\n", programname );
            break;
        case 'l':
	    p_logFileName = argv[2];
            argv++;
            argc--;
            argindex++;
            break;
        case 's':
            if ( (stateFile = ad_fopen(argv[2], "w")) == NULL ) {
#ifdef DEBUG
                fprintf(stderr, "\n State file name = %s\n", argv[2]); 
#endif /* DEBUG */
                fprintf(stderr, "\n%s: can't create state file %s\n", programname, argv[2]);
                fprintf(stderr, "\n%s: Unsuccessful Completion.\n\n", programname);
                return(-1);
            }
            else{
                fprintf(stateFile, "<?xml version=\"1.0\" ?>\n");
                fprintf(stateFile, "<autodock>\n");
                fprintf(stateFile, "\t<version>%s</version>\n", version_num);
                fprintf(stateFile, "\t<autogrid_version>%s</autogrid_version>\n", version_num);
                fprintf(stateFile, "\t<output_xml_version>%5.2f</output_xml_version>\n", OUTPUT_XML_VERSION);
                write_stateFile = TRUE;
            }
            argv++;
            argc--;
            argindex++;
            break;    
        case 'p':
            snprintf(dock_param_fn, PATH_MAX -1, "%s", argv[2] );
            if ( strindex( dock_param_fn, ".dpf") != (int) strlen(dock_param_fn) - 4){
                fprintf(stderr, "\n AutoDock needs the extension of the docking parameter file to be \".dpf\"");
                fprintf(stderr, "\n%s: Unsuccessful Completion.\n\n", programname);
                return(-1);
            }

            if ( (parFile = ad_fopen(dock_param_fn, "r")) == NULL ) {
#ifdef DEBUG
                fprintf(stderr, "\n Parameter file name = %s\n", dock_param_fn);
#endif /* DEBUG */
                fprintf(stderr, "\n%s: can't find or open parameter file %s\n", programname, dock_param_fn);
                fprintf(stderr, "\n%s: Unsuccessful Completion.\n\n", programname);
                return(-1);
            }
            argv++;
            argc--;
            argindex++;
            break;
        case 't':
            parse_tors_mode = TRUE;
            break;
        case 'v':
            fprintf(stdout, "AutoDock %-8s\n", version_num);
            fprintf(stdout, " Copyright (C) 2009 The Scripps Research Institute.\n");
            fprintf(stdout, " License GPLv2+: GNU GPL version 2 or later <http://gnu.org/licenses/gpl.html>\n");
            fprintf(stdout, " This is free software: you are free to change and redistribute it.\n");
            fprintf(stdout, " There is NO WARRANTY, to the extent permitted by law.\n");
            exit(0);
            break;
        default:
            fprintf(stderr, "%s: unknown switch \"-%c\".  \n", programname, argv[1][1]);
            usage(stderr, programname);
            return(-1);
            /* break; */
        }
        argindex++;
        argc--;
        argv++;
    }
    // set docking log file name from parameter file name if
    // a "-p <DPF>" appeared but no "-l <DLG>" appeared
    if ( parFile != stdin  && 0==strcmp(p_logFileName, "stdout")) {
            strncpy(logFileName, dock_param_fn, strlen(dock_param_fn)-4);
	    logFileName[strlen(dock_param_fn)] = '\0';
            strcat(logFileName, ".dlg");
	    }
    else snprintf(logFileName, sizeof logFileName, "%s", p_logFileName);

    if ( (logFile = ad_fopen(logFileName, "w")) == NULL ) {
#ifdef DEBUG
                fprintf(stderr, "\n Log file name = %s\n", logFileName); 
#endif /* DEBUG */
                fprintf(stderr, "\n%s: can't create log file %s\n", programname, logFileName);
                fprintf(stderr, "\n%s: Unsuccessful Completion.\n\n", programname);
                return(-1);
            }
            setlinebuf(logFile); // to ensure output even if crash
    return(argindex);
}
/* EOF */
