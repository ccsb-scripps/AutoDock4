/*

 $Id: banner.cc,v 1.3 2004/02/12 04:32:14 garrett Exp $

*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/* banner.cc */

    #include <stdio.h>
    #include "banner.h"

extern FILE *logFile;

void banner( FloatOrDouble version_num )

{

/*----------------------------------------------------------------------------*/
/* Output banner...                                                           */
/*----------------------------------------------------------------------------*/

    (void) fprintf(logFile, "      ________________________________________________________________\n");
    (void) fprintf(logFile, "\n");
    (void) fprintf(logFile, "__________//___________________________/////___________________/____________\n");
    (void) fprintf(logFile, "_________/__/__________________________/____/__________________/____________\n");
    (void) fprintf(logFile, "________/____/___________/_____________/_____/_________________/____________\n");
    (void) fprintf(logFile, "________/____/__/_____/_/////___/////__/_____/__/////___/////__/___/________\n");
    (void) fprintf(logFile, "_______/______/_/_____/__/_____/_____/_/_____/_/_____/_/_____/_/_//_________\n");
    (void) fprintf(logFile, "_______////////_/_____/__/_____/_____/_/_____/_/_____/_/_______//_/_________\n");
    (void) fprintf(logFile, "_______/______/_/____//__/___/_/_____/_/____/__/_____/_/_____/_/___/________\n");
    (void) fprintf(logFile, "_______/______/__////_/___///___/////__/////____/////___/////__/____/_______\n\n");
    (void) fprintf(logFile, "      ________________________________________________________________\n");
    (void) fprintf(logFile, "\n");
    (void) fprintf(logFile, "                                ______\n");
    (void) fprintf(logFile, "                               /      \\\n");
    (void) fprintf(logFile, "                              /        \\\n");
    (void) fprintf(logFile, "                             /          \\\n");
    (void) fprintf(logFile, "                             \\    /\\    /\n");
    (void) fprintf(logFile, "                              \\  /  \\  /\n");
    (void) fprintf(logFile, "                               \\/ /\\ \\/\n");
    (void) fprintf(logFile, "                                 /  \\\n");
    (void) fprintf(logFile, "                                /____\\\n");
    (void) fprintf(logFile, "\n");
    (void) fprintf(logFile, "                  ______________________________________ \n");
    (void) fprintf(logFile, "                 |                                      |\n");
    (void) fprintf(logFile, "                 |            AutoDock %3.2f             |\n", version_num );
    (void) fprintf(logFile, "                 |                                      |\n");
    (void) fprintf(logFile, "                 |            (c) 1991-2004             |\n");
    (void) fprintf(logFile, "                 |    The Scripps Research Institute    |\n");
    (void) fprintf(logFile, "                 |                                      |\n");
    (void) fprintf(logFile, "                 |       Garrett M. Morris, TSRI        |\n");
    (void) fprintf(logFile, "                 |         garrett@scripps.edu          |\n");
    (void) fprintf(logFile, "                 |                                      |\n");
    (void) fprintf(logFile, "                 |       David S. Goodsell, TSRI        |\n");
    (void) fprintf(logFile, "                 |         goodsell@scripps.edu         |\n");
    (void) fprintf(logFile, "                 |                                      |\n");
    (void) fprintf(logFile, "                 |        Arthur J. Olson, TSRI         |\n");
    (void) fprintf(logFile, "                 |           olson@scripps.edu          |\n");
    (void) fprintf(logFile, "                 |______________________________________|\n");
    (void) fprintf(logFile, "\n");
    (void) fprintf(logFile, "                  ______________________________________ \n");
    (void) fprintf(logFile, "                 |                                      |\n");
    (void) fprintf(logFile, "                 | Automated Docking of Flexible Ligand |\n");
    (void) fprintf(logFile, "                 |     to Macromolecular Receptor       |\n");
    (void) fprintf(logFile, "                 |______________________________________|\n");
    (void) fprintf(logFile, "\n\n\n");

}
/*----------------------------------------------------------------------------*/
/* EOF.                                                                       */
/*----------------------------------------------------------------------------*/
