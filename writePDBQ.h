#ifdef WRITEPDBQSTATE
/* writePDBQ.h */
#endif /* WRITEPDBQSTATE */

#ifndef WRITEPDBQSTATE

#ifndef WRITEPDBQ
#define WRITEPDBQ
#include "constants.h"
#include "printEnergies.h"
void  writePDBQ( int   irun,
                    char  smFileName[MAX_CHARS],
		    char  dpfFN[MAX_CHARS],
                    float sml_center[SPACE],
                    State stat,
                    int   ntor,
                    float eintra,
                    float einter,
                    int   natom,
                    char  atomstuff[MAX_ATOMS][MAX_CHARS],
                    float crd[MAX_ATOMS][SPACE],
                    float emap[MAX_ATOMS],
                    float elec[MAX_ATOMS],
                    float charge[MAX_ATOMS],
                    int ligand_is_inhibitor,
		    float torsFreeEnergy);
#endif

#else /* WRITEPDBQSTATE */

/* writePDBQ.h */

#ifndef PRINTDOCKED
#define PRINTDOCKED

#include "constants.h"
#include "printEnergies.h"
#include "trilinterp.h"
#include "eintcal.h"
#include "cnv_state_to_coords.h"
#include "stateLibrary.h"

void writeStateOfPDBQ(  int   irun,
		    char  smFileName[MAX_CHARS],
		    char  dpfFN[MAX_CHARS],
		    float sml_center[SPACE],
		    State (*Ptr_state),
		    int   ntor,
		    float (*Ptr_eintra),
		    float (*Ptr_einter),
		    int   natom,
		    char  atomstuff[MAX_ATOMS][MAX_CHARS],
		    float crd[MAX_ATOMS][SPACE],
		    float emap[MAX_ATOMS],
		    float elec[MAX_ATOMS],
		    float charge[MAX_ATOMS],
        int ligand_is_inhibitor,
		    float torsFreeEnergy,
                    float vt[MAX_TORS][SPACE],
                    int   tlist[MAX_TORS][MAX_ATOMS],
                    float crdpdb[MAX_ATOMS][SPACE],
                    int   nonbondlist[MAX_NONBONDS][2],
                    float e_internal[NEINT][ATOM_MAPS][ATOM_MAPS],
                    int   type[MAX_ATOMS],
                    int   Nnb,
                    Boole B_calcIntElec,
                    float q1q2[MAX_NONBONDS],
                    float map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS],
                    float inv_spacing,
                    float xlo,
                    float ylo,
                    float zlo,
                    float xhi,
                    float yhi,
                    float zhi,
                    Boole B_template,
                    float template_energy[MAX_ATOMS],
                    float template_stddev[MAX_ATOMS]
		    );

#endif

#endif /* WRITEPDBQSTATE */

void writeMolAsPDBQ(Molecule *mol, FILE *output);

