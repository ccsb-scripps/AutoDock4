/*

 $Id: setflags.cc,v 1.3 2004/02/12 04:32:16 garrett Exp $

*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/* setflags.cc */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "setflags.h"


extern FILE *parFile;
extern FILE *logFile;
extern char *programname;

extern char dock_param_fn[];
extern char AutoDockHelp[];
extern int  debug;
extern int  ignore_errors;
extern int  command_mode;
extern int  oldpdbq;
extern int  parse_tors_mode;
extern int  keepresnum;


int setflags( int I_argc, char * const PPC_argv[])

/*
** naming convention: 
** I_var   => var is an integer variable;
** PPC_var => var is a pointer to a pointer to a character variable
*/

/******************************************************************************/
/*      Name: setflags                                                        */
/*  Function: read flags from PPC_argv; return I_argindex of first non arg.   */
/* Copyright: (C) Garrett Matthew Morris, TSRI.                               */
/*----------------------------------------------------------------------------*/
/*    Author: Garrett Matthew Morris, TSRI.                                   */
/*            (Adapted from code supplied by Bruce Duncan, TSRI.)             */
/*      Date: 02/02/94                                                        */
/*----------------------------------------------------------------------------*/
/*    Inputs: I_argc,PPC_argv                                                 */
/*   Returns: I_argindex                                                      */
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
    int I_argindex;
/*----------------------------------------------------------------------------*/
/* Initialize                                                                 */
/*----------------------------------------------------------------------------*/
    I_argindex = 1;
    programname = PPC_argv[0];
    parFile = stdin;
    logFile = stdout;
    /*
     * see autoglobal.h for initialization of debug, keepresnum and logicals...
     */
/*----------------------------------------------------------------------------*/
/* Loop over arguments                                                        */
/*----------------------------------------------------------------------------*/
    while((I_argc > 1) && (PPC_argv[1][0] == '-')){
        switch(PPC_argv[1][1]){
#ifdef FOO
        case 'n':
            ncount = atoi(PPC_argv[2]);
            PPC_argv++;
            I_argc--;
            I_argindex++;
            break;
#endif
        case 'd':
            debug++;
            break;
        case 'u':
            usage();
	    exit(0);
            break;
        case 'o':
            oldpdbq = TRUE;
            break;
        case 'i':
            ignore_errors = TRUE;
            break;
        case 'k':
            keepresnum--;
            break;
        case 'c':
            command_mode = TRUE;
            break;
        case 'l':
            if ( (logFile = fopen(PPC_argv[2], "w")) == NULL ) {
#ifdef DEBUG
                fprintf(stderr,"\n Log file name = %s\n",PPC_argv[2]); 
#endif /* DEBUG */
                fprintf(stderr, "\n%s: can't create log file %s\n", programname, PPC_argv[2]);
                fprintf(stderr, "\n%s: Unsuccessful Completion.\n\n", programname);
                return(-1);
            }
            PPC_argv++;
            I_argc--;
            I_argindex++;
            break;
        case 'p':
            strcpy(dock_param_fn, PPC_argv[2]);
            if ( (parFile = fopen(PPC_argv[2], "r")) == NULL ) {
#ifdef DEBUG
                fprintf(stderr,"\n Parameter file name = %s\n",PPC_argv[2]);
#endif /* DEBUG */
                fprintf(stderr, "\n%s: can't find or open parameter file %s\n", programname, PPC_argv[2]);
                fprintf(stderr, "\n%s: Unsuccessful Completion.\n\n", programname);
                return(-1);
            }
            PPC_argv++;
            I_argc--;
            I_argindex++;
            break;
        case 't':
            parse_tors_mode = TRUE;
            break;
        default:
            fprintf(stderr,"%s: unknown switch \"-%c\".  Usage:\n",programname,PPC_argv[1][1]);
            fprintf(stderr,"%s %s\n",programname,AutoDockHelp);
            return(-1);
            /* break; */
        }
        I_argindex++;
        I_argc--;
        PPC_argv++;
    }
    return(I_argindex);
}
/* EOF */
