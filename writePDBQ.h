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
                 FloatOrDouble sml_center[SPACE],
                 State stat,
                 int   ntor,
                 FloatOrDouble eintra,
                 FloatOrDouble einter,
                 int   natom,
                 char  atomstuff[MAX_ATOMS][MAX_CHARS],
                 FloatOrDouble crd[MAX_ATOMS][SPACE],
                 FloatOrDouble emap[MAX_ATOMS],
                 FloatOrDouble elec[MAX_ATOMS],
                 FloatOrDouble charge[MAX_ATOMS],
                 FloatOrDouble abs_charge[MAX_ATOMS],
                 FloatOrDouble qsp_abs_charge[MAX_ATOMS],
                 int ligand_is_inhibitor,
                 FloatOrDouble torsFreeEnergy,
                 int outlev,
                 int   ignore_inter[MAX_ATOMS],
                 const Boole         B_include_1_4_interactions,
                 const FloatOrDouble scale_1_4,
                 const FloatOrDouble sol_fn[NEINT],
                 const ParameterEntry parameterArray[MAX_MAPS]
                 );
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

void writeStateOfPDBQ(int irun,FourByteLong seed[2],
                    char  smFileName[MAX_CHARS],
                    char  dpfFN[MAX_CHARS],
                    FloatOrDouble sml_center[SPACE],
                    State (*Ptr_state),
                    int   ntor,
                    FloatOrDouble (*Ptr_eintra),
                    FloatOrDouble (*Ptr_einter),
                    int   natom,
                    char  atomstuff[MAX_ATOMS][MAX_CHARS],
                    FloatOrDouble crd[MAX_ATOMS][SPACE],
                    FloatOrDouble emap[MAX_ATOMS],
                    FloatOrDouble elec[MAX_ATOMS],
                    FloatOrDouble charge[MAX_ATOMS],
                    FloatOrDouble abs_charge[MAX_ATOMS],
                    FloatOrDouble qsp_abs_charge[MAX_ATOMS],
                    int ligand_is_inhibitor,
                    FloatOrDouble torsFreeEnergy,
                    FloatOrDouble vt[MAX_TORS][SPACE],
                    int   tlist[MAX_TORS][MAX_ATOMS],
                    FloatOrDouble crdpdb[MAX_ATOMS][SPACE],
                    int   nonbondlist[MAX_NONBONDS][MAX_NBDATA],
                    FloatOrDouble e_internal[NEINT][ATOM_MAPS][ATOM_MAPS],
                    int   type[MAX_ATOMS],
                    int   Nnb,
                    Boole B_calcIntElec,
                    FloatOrDouble q1q2[MAX_NONBONDS],
                    FloatOrDouble map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS],
                    FloatOrDouble inv_spacing,
                    FloatOrDouble xlo,
                    FloatOrDouble ylo,
                    FloatOrDouble zlo,
                    FloatOrDouble xhi,
                    FloatOrDouble yhi,
                    FloatOrDouble zhi,
                    Boole B_template,
                    FloatOrDouble template_energy[MAX_ATOMS],
                    FloatOrDouble template_stddev[MAX_ATOMS],
                    int outlev,
                    int   ignore_inter[MAX_ATOMS],
                    const Boole         B_include_1_4_interactions,
                    const FloatOrDouble scale_1_4,
                    const FloatOrDouble sol_fn[NEINT],
                    const ParameterEntry parameterArray[MAX_MAPS]
		    );

#endif

#endif /* WRITEPDBQSTATE */
extern FILE *stateFile;
extern int write_stateFile;


void writeMolAsPDBQ(Molecule *mol, FILE *output);

