/*

 $Id: stop.cc,v 1.5 2009/06/10 00:09:09 rhuey Exp $

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
#include <string>
#include "stop.h"

using std::string;

extern char *programname;
extern FILE *logFile;


/*----------------------------------------------------------------------------*/
void stop(const char *reason)
/*----------------------------------------------------------------------------*/
{

    if (logFile == stdout) {
	fprintf( logFile, "%s: ERROR: %s\n", programname, reason);
	fprintf( logFile, "%s: Aborting...\n\n", programname);
	fprintf( logFile, "%s: Unsuccessful Completion.\n\n", programname);
	fflush(  logFile  );
    } else {
        string pn=programname;
        string r=reason;
	string message =  pn + ": ERROR: "+r+"\n";
	print_2x( logFile, stderr, message.c_str() );
	message = pn + ": Aborting...\n\n";
	print_2x( logFile, stderr, message.c_str() );
	message = pn + ": Unsuccessful Completion.\n\n";
	print_2x( logFile, stderr, message.c_str() );
	fflush(logFile);
	fflush(stderr);
    }

    exit(-1);
}
/* EOF */
