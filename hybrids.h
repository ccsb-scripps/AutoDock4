/*******************************************************************************
	Due to the fact that *.c files are compiled differently than *.cc
	(because of the mangled C++ names), this is a separate header file
	that serves much the same purpose as other prototypes headers.

			2/6/96  rsh
    
        Can get rid of the individual #if-#endif brackets around the fcns.
*******************************************************************************/
#include "constants.h"
#ifndef _STRUCTS_H
#include "structs.h"
#endif
#include "gs.h"
#include "ls.h"
#include "support.h"

#ifndef CALL_GLSS
#define CALL_GLSS

State call_glss(Global_Search *global_method, Local_Search *local_method, 
		State now, unsigned int num_evals, unsigned int pop_size, 
		int outlev, unsigned int extOutputEveryNgens, Molecule *mol,
		Boole B_RandomTran0, Boole B_RandomQuat0, Boole B_RandomDihe0,
        GridMapSetInfo *info, char FN_pop_file[MAX_CHARS] );


Representation **generate_R(int num_torsions, GridMapSetInfo *info );

Representation **generate_R_quaternion(int num_torsions, GridMapSetInfo *info );

Genotype generate_Gtype(int num_torsions, GridMapSetInfo *info );

Phenotype generate_Ptype(int num_torsions, GridMapSetInfo *info );

Individual random_ind(int num_torsions, GridMapSetInfo *info );

#endif


#ifndef CALL_GLSS_TORS
#define CALL_GLSS_TORS

State call_glss_tors(Global_Search *global_method, Local_Search *local_method, 
		State now, unsigned int num_evals, unsigned int pop_size, 
		int outlev, unsigned int extOutputEveryNgens, Molecule *mol,
		Boole B_RandomDihe0,
        GridMapSetInfo *info, char FN_pop_file[MAX_CHARS] );

Representation **generate_R_tors(int num_torsions, GridMapSetInfo *info );

Genotype generate_Gtype_tors(int num_torsions, GridMapSetInfo *info );

Phenotype generate_Ptype_tors(int num_torsions, GridMapSetInfo *info );

Individual random_ind_tors(int num_torsions, GridMapSetInfo *info );

#endif

#ifndef CALL_LS
#define CALL_LS

State call_ls(Local_Search *local_method, State now, unsigned int pop_size, Molecule *mol);

#endif


#ifndef CALL_GS
#define CALL_GS

State call_gs(Global_Search *global_method, State now, unsigned int num_evals, unsigned int pop_size,
              Molecule *mol,
              int extOutputEveryNgens,
              GridMapSetInfo *info );

#endif


#ifndef MMM
#define MMM

void minmeanmax( FILE *fp, Population &pop, int num_its, GridMapSetInfo *info );

#endif
