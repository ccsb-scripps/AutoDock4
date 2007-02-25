
#ifndef QTRANSFORM
#define QTRANSFORM

#include <stdio.h>
#include "constants.h"
#include "structs.h"
#include "writePDBQT.h"
#include "qmultiply.h"
#include "torNorVec.h"

void qtransform( const Coord T,
	 	 const Quat  q,
                 Real tcoord[MAX_ATOMS][SPACE],
		 const int   natom);

void reorient( FILE *logFile, 
               const int true_ligand_atoms, 
               char atomstuff[MAX_ATOMS][MAX_CHARS],
               Real crdpdb[MAX_ATOMS][SPACE],  // original PDB coordinates from input
               Real charge[MAX_ATOMS],
               int type[MAX_ATOMS],
               ParameterEntry  parameterArray[MAX_MAPS],
               Quat q_reorient,
               Coord origin,
               const int ntor,
               int tlist[MAX_TORS][MAX_ATOMS],
               Real vt[MAX_TORS][SPACE],
               Molecule *ptr_ligand,
               const int debug );
#endif
