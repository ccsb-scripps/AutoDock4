/* dpftypes.cc */

    #include <stdio.h>
    #include <stdlib.h>
    #include <string.h>
    #include "dpftypes.h"


extern int ElecMap;
extern FILE *logFile;


void dpftypes( int *P_Htype, 
	       int *P_N_all_maps, 
	       int *P_N_atm_maps, 
	       char S_atomType[ATOM_MAPS], 
	       char S_line[LINE_LEN] )

{
    register int i;
    char S_errorMessage[LINE_LEN];
    char FN_ligand[MAX_CHARS];

    /*
    \  ATOM_TYPE_NAMES
     \______________________________________________________________
    */
    sscanf( S_line, "%*s %s", S_atomType);

    ElecMap = *P_N_atm_maps = strlen( S_atomType );
    *P_N_all_maps = *P_N_atm_maps + 1;
    /*|
      |  e.g.        S_atomType   *P_N_atm_maps   *P_N_all_maps
      |                CNOH                4                5
      |       C-arrays:0123,
      |    electrostatics: 4, or *P_N_atm_maps
      |                     =strlen( S_atomType );
      |*/
    if ( *P_N_atm_maps > ATOM_MAPS ) {
      prStr( S_errorMessage, "ERROR: %d atom types declared in \"%s\", S_line \"%s\"; maximum allowed is %d.\nChange \"ATOM_MAPS\" in \"autocomm.h\"", *P_N_atm_maps, FN_ligand, S_atomType, ATOM_MAPS );
      stop( S_errorMessage );
      exit( -1 );
    }
    for ( i = 0; i < *P_N_atm_maps; i++) {
      if ( S_atomType[i] == 'H' ) {
	*P_Htype = i;
      }
    }    /* i */
    pr( logFile, "Atom descriptors for types 1-%d = ", *P_N_atm_maps );
    for (i=0; i < *P_N_atm_maps; i++) {
      pr( logFile, "%c%s", S_atomType[i], (i==(*P_N_atm_maps - 1))?" ":", " );
    }
    pr( logFile, " respectively.\n\n" );
    pr( logFile, "Electrostatic energies will be stored in Grid Map %d\n", ElecMap+1 );
}
/* EOF */
