/* AutoDock
 $Id: main.cc,v 1.114 2009/12/15 06:21:01 mp Exp $

**  Function: Performs Automated Docking of Small Molecule into Macromolecule
**Copyright (C) 2009 The Scripps Research Institute. All rights reserved.
** All Rights Reserved.
**____________________________________________________________________________
** Primary Authors: 
**            Garrett Matthew Morris, C/C++ version 
**            David Goodsell, Original FORTRAN version 1.0
**                                       e-mail: goodsell@scripps.edu
**
** Other Contributors: see file AUTHORS
**            
**            Laboratory of Arthur J. Olson
**            The Scripps Research Institute
**            Department of Molecular Biology, MB5
**            10550 North Torrey Pines Road
**            La Jolla, CA 92037.
**____________________________________________________________________________
**    Inputs: Control file, Small Molecule PDBQT file, macromolecular grid map
**            files.
**   Returns: Autodock Log File, includes docked conformation clusters (PDBQT)
**____________________________________________________________________________
** Modification Record (pre-CVS)
** Date     Inits   Comments
** 09/06/95 RSH     Added code to handle GA/LS stuff
**          DSG     Quaternion rotations
**          DSG     Generates torsions from annotated pdb list
**          DSG     Generates internal energies
**          DSG     Performs a limited Cluster analysis of conformations
** 05/07/92 GMM     C translation
** 05/14/92 GMM     Time-dependent seed in random-number generation
** 10/29/92 GMM     Application Visualization System (AVS) readable grid
**                  display file input.
**                  [AVS is a trademark of Stardent Computer Inc.]
**                  Also added the 'total_charge' check.
** 11/19/93 GMM     #ifdef NOSQRT, with non-square-rooting acceleration.
** 09/26/94 GMM     Cluster analysis now outputs RMS deviations.
** 09/28/94 GMM     Modularized code.
** 10/02/94 GMM     Distance constraints added, for Ed Moret. Accelerated.
** 09/06/95 RSH     Incorporation of GA/SW tokens
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
*******************************************************************************/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>
#include <sys/types.h> // time_t time(time_t *tloc);
#include <time.h>      // time_t time(time_t *tloc);
#include <sys/times.h>
#include <stdlib.h>
#include <string>
using std::string;

#include <sys/param.h>
#include <ctype.h> // tolower
#include <unistd.h> // sysconf

/* the BOINC API header file */
#ifdef BOINC
#include "diagnostics.h"
#include "boinc_api.h"
#include "filesys.h"                 // boinc_fopen(), etc... */
#endif

#include "coliny.h"
#include "hybrids.h"
#include "ranlib.h"
#include "gs.h"
#include "ls.h"
#include "rep.h"
#include "support.h"
#include "distdepdiel.h"
#include "calculateEnergies.h"
#include "conformation_sampler.h"
#include "main.h"
#include "alea.h"
#include "call_cpso.h"
#include "dimLibrary.h"


extern int debug;
extern int keepresnum;
extern Real idct;
extern Eval evaluate;
extern Linear_FE_Model AD4;
extern Real nb_group_energy[3]; ///< total energy of each nonbond group (intra-ligand, inter, and intra-receptor)
extern int Nnb_array[3];  ///< number of nonbonds in the ligand, intermolecular and receptor groups

static const char* const ident[] = {ident[1], "@(#)$Id: main.cc,v 1.114 2009/12/15 06:21:01 mp Exp $"};


int sel_prop_count = 0;
static Boole B_found_about_keyword = FALSE; //set false by 'move' true by 'about'
static Boole B_found_elecmap = FALSE;
static Boole B_found_desolvmap = FALSE;
static void exit_if_missing_elecmap_desolvmap_about(string  keyword); // see bottom of main.cc

//------------------------------- PSO -Work Variables declaration 

// State Structure Variable DECLARATION
State sNew[S_max];
int D; //search space dimension 
int nb_eval; // Current number of Swarm evaluations
int eval_max = 0; //Max number of Swarm iterations
int S ; // Swarm size
int S_factor=30; // Swarm size
double xmin[D_max], xmax[D_max]; // Intervals defining the search space
double Vmin[D_max], Vmax[D_max]; // Intervals defining the search space

int main (int argc, char ** argv)


{

//   MAX_MAPS
//

//   MAX_GRID_PTS & MAX_MAPS
//
static MapType map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS];
// double *map;  // Use this with malloc...

GridMapSetInfo *info;  // this information is from the AVS field file

//   MAX_ATOMS
//
char atomstuff[MAX_ATOMS][MAX_CHARS];
char pdbaname[MAX_ATOMS][5];
Real crdpdb[MAX_ATOMS][SPACE];  // original PDB coordinates from input
Real crdreo[MAX_ATOMS][SPACE];  // reoriented coordinates
Real crd[MAX_ATOMS][SPACE];     // current coordinates
Real charge[MAX_ATOMS];
Real abs_charge[MAX_ATOMS];
Real qsp_abs_charge[MAX_ATOMS];
Real elec[MAX_ATOMS];
Real emap[MAX_ATOMS];
int type[MAX_ATOMS];
int bond_index[MAX_ATOMS];
int ignore_inter[MAX_ATOMS];
Atom atoms[MAX_ATOMS];

//   MAX_TORS
//
int  tlist[MAX_TORS][MAX_ATOMS];
Real vt[MAX_TORS][SPACE];
Real F_TorConRange[MAX_TORS][MAX_TOR_CON][2];
unsigned short  US_TorE[MAX_TORS];
Boole B_isTorConstrained[MAX_TORS];
int N_con[MAX_TORS];
unsigned short US_torProfile[MAX_TORS][NTORDIVS];

//   MAX_NONBONDS
//
NonbondParam *nonbondlist = new NonbondParam[MAX_NONBONDS];

//   LINE_LEN
//
char error_message[LINE_LEN];
char message[LINE_LEN];
char line[LINE_LEN];
char torfmt[LINE_LEN];
char param[2][LINE_LEN];
char rms_atoms_cmd[LINE_LEN];
char c_mode_str[LINE_LEN];
char confsampler_type[LINE_LEN];
char autodock_parameter_version[LINE_LEN]; //eg 4.1.1

// filename max length is taken from system include file
char FN_parameter_library[PATH_MAX];
char FN_ligand[PATH_MAX];
char FN_flexres[PATH_MAX];
char FN_rms_ref_crds[PATH_MAX];
char FN_clus[PATH_MAX];
char FN_watch[PATH_MAX];
char dummy_FN_ligand[PATH_MAX];
char FN_pop_file[PATH_MAX];
char FN_trj[PATH_MAX];
//char FN_gdfld[PATH_MAX]  // now part of the GridMapSetInfo structure;
//char FN_gpf[PATH_MAX];  // now part of the GridMapSetInfo structure
//char FN_receptor[PATH_MAX];// now part of the GridMapSetInfo structure

//   MAX_CHARS
char hostnm[MAX_CHARS];

//   MAX_RECORDS
//
char PDBQT_record[MAX_RECORDS][LINE_LEN];

//   SPACE
//
Real lig_center[SPACE];
//Real map_center[SPACE];

//   MAX_RUNS
//
Real econf[MAX_RUNS];  // this is the list of energies printed in the histogram in "analysis"
State sHist[MAX_RUNS];  /* qtnHist[MAX_RUNS][QUAT],torHist[MAX_RUNS][MAX_TORS];*/
State sUnbound; // State of the unbound ligand's conformation
State sUnbound_ext; // State of the unbound ligand's conformation after extended-conformation search
//// UNCOMMENT if using Step 2 in unbound calculation ---> State sUnbound_ls; // State of the unbound ligand's conformation after a local search
State sUnbound_ad; // State of the unbound ligand's conformation after an AutoDock search

char out_acc_rej = '?';
char timeSeedIsSet[2];
char selminpar = 'm';
char S_contype[8];

//   MAX_ATOM_TYPES
//
char            *ligand_atom_type_ptrs[MAX_ATOM_TYPES]; /* array of ptrs used to parse input line of atom type names */
ParameterEntry parameterArray[MAX_ATOM_TYPES]; // input  nonbond and desolvation parameters
static ParameterEntry * foundParameter;

Real cA;
Real cB;
Real cA_unbound = 392586.8;  // repulsive
Real cB_unbound = 0.0; // attractive
Real epsij;
Real F_A;
Real F_Aova;
Real F_tor;
Real F_torPref;
Real F_torHWdth;
Real Rij;
Real sqlower;
Real squpper;
Real tmpconst;

// Distance-dependence in Desolvation Term
const double sigma = 3.6L;
const double qsolpar = 0.01097L;

// ELECSCALE converts between CGS units and SI units;
// see, e.g. p 254, "Molecular Modeling and Simulation", by Tamar Schlick, Springer.
//
// Units of ELECSCALE are (Kcal/mol ) * (Angstrom / esu^2)
// and this allows us to use distances in  Angstroms and charges in esu...
const Real ELECSCALE = 332.06363;

// const Real ELECSCALE = 83.0159075;   this ELECSCALE (corresponding to eps(r) = 1/4r) gives -7.13 kcal/mol for 1pgp Tests/test_autodock4.py

// i
double Ri, epsi, Ri_hb, epsi_hb;
hbond_type hbondi;

// j
double Rj, epsj, Rj_hb, epsj_hb;
hbond_type hbondj;

Real scale_1_4 = 0.5;
Real c=0.0;
Real clus_rms_tol = 0.0;
Real e0max = BIG;
Real eintra = 0.0;  // sum of intramolecular energy for the ligand plus that of the protein
Real einter = 0.0; // intermolecular energy between the ligand and the protein
Real etotal = 0.0;
Real AD3_FE_coeff_estat   = 1.000;
Real qtwFac = 1.0;
Real qtwStep0 = DegreesToRadians( 5.0 );
Real qtwStepFinal = DegreesToRadians( 0.5 );
Real maxrad = -1.0;
Real r2sum=0.0;
// Real RJ = 8.31441;     // in J/K/mol, Gas Constant, Atkins Phys.Chem., 2/e
// Real Rcal = 1.9871917; // in cal/K/mol, Gas Constant, RJ/4.184
// Real T0K = 273.15;        // 0 degrees Celsius, in K
Real RTreduc = 1.0;
// Real spacing = 0.0;// now part of the GridMapSetInfo structure
Real RT0 = 616.0;
Real RTFac = 0.95;
Real torsdoffac = 0.3113;
Real torsFreeEnergy = 0.0;
Real torFac = 1.0;
Real torStep0 = DegreesToRadians( 5.0 );
Real torStepFinal = DegreesToRadians( 1.0 );
Real trnFac = 1.0;
Real trnStep0 = 0.2;
Real trnStepFinal = 0.2;
Real WallEnergy = 1.0e8; /* Energy barrier beyond walls of gridmaps. */
//  The GA Stuff
Real m_rate = 0.02;
Real c_rate = 0.80;
Real alpha = 0;
Real beta = 1;
Real search_freq = 0.06;
Real unbound_internal_FE = 0.0;
Real unbound_ext_internal_FE = 0.0;
Real unbound_ad_internal_FE = 0.0;
Real emap_total = 0.;
Real elec_total = 0.;
Real charge_total = 0.;
Real etot = 0.;
//  The LS Stuff
Real rho = 1.0;
Real lb_rho = 0.01;
Real *rho_ptr = NULL;
Real *lb_rho_ptr = NULL;
Real psw_trans_scale = 1.0;//1 angstrom
Real psw_rot_scale = 0.05;  //about 3 degrees, we think
Real psw_tors_scale = 0.1; //about 6 degrees

Unbound_Model ad4_unbound_model = Unbound_Default; 
//NOT Same_As_Bound so user can specify in dpf

EnergyBreakdown eb;

initialise_energy_breakdown(&eb, 0, 0);

unsigned int outputEveryNgens = 100;

Boole B_atom_types_found = FALSE;
Boole B_isGaussTorCon = FALSE;
Boole B_constrain_dist;
Boole B_either = FALSE;
Boole B_calcIntElec = FALSE;
Boole B_calcIntElec_saved = FALSE;
Boole B_write_trj = FALSE;
Boole B_watch = FALSE;
Boole B_acconly = FALSE;
Boole B_cluster_mode = FALSE;
Boole B_havemap = FALSE;
Boole B_havenbp = FALSE;
Boole B_haveCharges;
Boole B_linear_schedule = FALSE;
Boole B_qtwReduc = FALSE;
Boole B_selectmin = FALSE;
Boole B_symmetry_flag = TRUE;
Boole B_tempChange = TRUE;
Boole B_torReduc = FALSE;
Boole B_trnReduc = FALSE;
Boole B_write_all_clusmem = FALSE;
Boole B_ShowTorE = FALSE;
Boole B_RandomTran0 = FALSE;
Boole B_RandomQuat0 = FALSE;
Boole B_RandomDihe0 = FALSE;
Boole B_CalcTrnRF = FALSE;
Boole B_CalcQtwRF = FALSE;
Boole B_CalcTorRF = FALSE;
Boole B_charMap = FALSE;
Boole B_include_1_4_interactions = FALSE;  // This was the default behaviour in previous AutoDock versions (1 to 3).
Boole B_found_move_keyword = FALSE;
Boole B_found_ligand_types = FALSE;
Boole B_found_autodock_parameter_version = FALSE;
Boole B_use_non_bond_cutoff = TRUE;
Boole B_have_flexible_residues = FALSE;  // if the receptor has flexible residues, this will be set to TRUE
Boole B_rms_atoms_ligand_only = TRUE;  // cluster on the ligand atoms only
Boole B_reorient_random = FALSE; // if true, create a new random orientation before docking
int atm1=0;
int atm2=0;
int a1=0;
int a2=0;
int atomC1;
int atomC2;
int dpf_keyword = -1;
//int gridpts1[SPACE];  // now part of the GridMapSetInfo structure
//int gridpts[SPACE];  // now part of the GridMapSetInfo structure
int n_heavy_atoms_in_ligand = 0;
int ncycles = -1;
int iCon=0;
int indcom = 0;
int ligand_is_inhibitor = 1;
int ltorfmt = 4;
int nruns = 0;
int nstepmax = -1;
int naccmax = 0;
int natom = 0;
int nconf = 0;
int ncycm1 = 1;
int ndihed = 0;
int nlig = 0;
int nres = 0;
int nmol = 0;
int Nnb = 0;
int nrejmax = 0;
int ntor = 0;
int ntor1 = 1;
int ntor_ligand = 0;
int ntorsdof = 0;
int num_maps = 0;
int num_atom_types = 0;
int nval = 0;
int outlev = -1;
int retval = 0;
int trj_end_cyc = 0;
int trj_begin_cyc = 0;
int trj_freq = 0;
int xA = 12;
int xB = 6;
int xA_unbound = 12;
int xB_unbound = 6;
int I_tor;
int I_torBarrier;
int MaxRetries = 1000; /* Default maximum number of retries for ligand init. */
int OutputEveryNTests = 1000;
int NumLocalTests = 10;
int maxTests = 10000;
int parameter_library_found = 0;
/* int beg; */
/* int end; */
/* int imol = 0; */
int outside = FALSE;
int atoms_outside = FALSE;
// unsigned int min_evals_unbound =  250000;
unsigned int max_evals_unbound = 1000000;
int saved_sInit_ntor = 0;
int confsampler_samples = 0;

unsigned short US_energy;
unsigned short US_tD;
unsigned short US_torBarrier = TORBARMAX;
unsigned short US_min = TORBARMAX;

register int i = 0;
register int j = 0;
register int k = 0;
//register int m = 0;
register int xyz = 0;
int j1 = 1;

State sInit;            /* Real qtn0[QUAT], tor0[MAX_TORS]; */

Quat q_reorient;

Molecule ligand;        /* ligand */

static Real F_A_from;
static Real F_A_to;
static Real F_lnH;
static Real F_W;
static Real F_hW;
static FourByteLong clktck = 0;

#ifndef VERSION
static string version_num = "4.2.2";
#else
static string version_num = VERSION;
#endif

struct tms tms_jobStart;
struct tms tms_gaStart;
struct tms tms_gaEnd;
struct tms tms_psoStartClock;
struct tms tms_psoEndClock;

Clock  jobStart;
Clock  gaStart;
Clock  gaEnd;

time_t time_seed;

EnergyTables *ad_energy_tables;  // Holds vdw+Hb, desolvation & dielectric lookup tables
EnergyTables *unbound_energy_tables;  // Use for computing unbound energy & conformation

Statistics map_stats;

//  The GA Stuff
FourByteLong seed[2];
unsigned int pop_size = 200;
unsigned int num_generations = 0;  //  Don't terminate on the basis of number of generations
unsigned int num_evals = 250000;
unsigned int num_evals_unbound = num_evals;
unsigned int max_its = 300;
unsigned int max_succ = 4;
unsigned int max_fail = 4;
int window_size = 10;
int low = 0;
int high = 100;
int elitism = 1;

// For Branch Crossover Mode
int end_of_branch[MAX_TORS];

Selection_Mode s_mode = Proportional;
Real linear_ranking_selection_probability_ratio = 2.0;
Xover_Mode c_mode = TwoPt;  //  can be: OnePt, TwoPt, Uniform or Arithmetic
Worst_Mode w_mode = AverageOfN;
EvalMode e_mode = Normal_Eval;
Global_Search *GlobalSearchMethod = NULL;
Local_Search *LocalSearchMethod = NULL;
//// Local_Search *UnboundLocalSearchMethod = NULL;

//---------PSO Variable Declaration 
// Local Variables 
//

// Confidence Coefficients
//double w = 0 ; // First confidence coefficient
//double w_start = 0.9 ; // First confidence coefficient-Start
//double w_end = 0.4 ; // First confidence coefficient-End
double c1 = 0; // second confidence coefficient
double c2 = 0; // second confidence coefficient
int K ; // Max number of particles informed by a given one 
int n_exec, n_exec_max; //Nbs of executions 

//Boole B_notinbox; //For checking in Grid Box
// SSM-PSO
//double mc = 0.3; // momentum constant for SSM-PSO

Clock  psoStartClock;
Clock  psoEndClock;


info = (GridMapSetInfo *) malloc( sizeof(GridMapSetInfo) );
ad_energy_tables = (EnergyTables *) malloc( sizeof(EnergyTables) );
unbound_energy_tables = (EnergyTables *) malloc( sizeof(EnergyTables) );

// Create a coordinate at the origin:
Coord origin;
origin.x = 0.;
origin.y = 0.;
origin.z = 0.;

//______________________________________________________________________________
/*
** Get the time at the start of the run...
*/

jobStart = times( &tms_jobStart );


//_____________________________________________________________________________
/*
** Boinc initialization
*/
#ifdef BOINC
    int flags = 0;
    int rc;
    flags =
      BOINC_DIAG_DUMPCALLSTACKENABLED |
      BOINC_DIAG_HEAPCHECKENABLED |
      BOINC_DIAG_REDIRECTSTDERR |
      BOINC_DIAG_REDIRECTSTDOUT ;
    boinc_init_diagnostics(flags);

#ifdef BOINCCOMPOUND
    BOINC_OPTIONS options;
    options.main_program = false;
    options.check_heartbeat = false; // monitor does check heartbeat
    options.handle_trickle_ups = false;
    options.handle_trickle_downs = false;
    options.handle_process_control = false;
    options.send_status_msgs = true;// only the worker programs (i.e. model) sends status msgs
    options.direct_process_action = true;// monitor handles suspend/quit, but app/model doesn't
    // Initialisation of Boinc
    rc =  boinc_init_options(options); //return 0 for success
    if( rc ){
      fprintf(stderr,"BOINC_ERROR: boinc_init_options() failed \n");
      exit(rc);
    }

#else
    // All BOINC applications must initialize the BOINC interface:
    rc = boinc_init();
    if (rc){
      fprintf(stderr, "BOINC_ERROR: boinc_init() failed.\n");
      exit(rc);
    }
#endif
#endif


//______________________________________________________________________________
/*
** Parse the arguments in the command line...
*/

if ( setflags(argc,argv,version_num.c_str()) == -1) {
    exit(-1);
} /* END PROGRAM */


//______________________________________________________________________________
/*
** Initialize torsion arrays and constants.
*/

(void) strcpy( torfmt, "%*s" ); /* len(torfmt) is 3 chars */

for (j = 0;  j < MAX_ATOMS;  j++ ) {
    type[j] = 0;
    ignore_inter[j] = 0;
}

for (i = 0; i  < MAX_TORS;  i++ ) {
    for (j = 0;  j < MAX_ATOMS;  j++ ) {
        tlist[i][j] = 0;
    }
}

for (i = 0; i  < MAX_TORS;  i++ ) {
    B_isTorConstrained[i] = 0;
    US_torProfile[i][0] = 0;
    N_con[i] = 0;
}

initialiseState( &sInit );
initialiseState( &(ligand.S) );

initialiseQuat( &q_reorient );

F_W      =  360.0 / NTORDIVS;
F_hW     =  F_W  / 2.0;
F_A_from = -360.0 + F_hW;
F_A_to   =  360.0 + F_hW;

for (k = 0; k < MAX_RUNS; k++) {
    for (i = 0; i  < MAX_TORS;  i++ ) {
        sHist[k].tor[i] = 0.0;
    }
}

for (i = 0; i < MAX_TORS;  i++ ) {
    if ( (ltorfmt += 4) > LINE_LEN ) {
        prStr( error_message, "%s:  ERROR: MAX_TORS = %d torsions declared in \"constants.h\";\n\t LINE_LEN = %d, Therefore you must change \"LINE_LEN\" to exceed %d...\n", programname, MAX_TORS, LINE_LEN, 4+4*MAX_TORS );
        stop( error_message );
        exit( -1 );
    } else {
        (void) strcat( torfmt, " %lf" );  /* add on 4 chars  for each new torsion... */
    }
} /* len(torfmt) is 3+4*MAX_TORS chars */

for (j = 0; j < MAX_NONBONDS; j++) {
    nonbondlist[j].a1 = nonbondlist[j].a2 = 0;
}

for (j = 0; j < MAX_RUNS; j++) {
    econf[j] = 0.0;
}

B_constrain_dist = B_haveCharges = FALSE;
ntor1 = ntor = atomC1 = atomC2 = 0;
sqlower = squpper = 0.0;

timeSeedIsSet[0] = 'F';
timeSeedIsSet[1] = 'F';

if (clktck == 0) {        /* fetch clock ticks per second first time */
    if ( (clktck = sysconf(_SC_CLK_TCK)) < (FourByteLong)0L) {
        stop("\"sysconf(_SC_CLK_TCK)\" command failed in \"main.c\"\n");
        exit( -1 );
    } else {
        idct = (Real)1.0 / (Real)clktck;
        if (debug) {
            pr(logFile, "N.B. debug is on and set to %d\n\n", debug);
            pr(logFile, "\n\nFYI:  Number of clock ticks per second = %d\n", (int)clktck);
            pr(logFile, "FYI:  Elapsed time per clock tick = %.3e milli-seconds\n\n\n\n", idct * 1000. );
        }
    }
}

(void) strcpy(FN_rms_ref_crds,"unspecified filename\0");


//______________________________________________________________________________
/*
** log(x): compute the natural (base e) logarithm of x,
*/

F_lnH = ((Real)log(0.5));

//______________________________________________________________________________
/*
** Determine output level before we ouput anything.
** We must parse the entire DPF -- silently -- for any outlev settings
** or flexible residues file specification
*/

while( fgets(line, LINE_LEN, parFile) != NULL ) { /* PARSING-DPF parFile */
    dpf_keyword = parse_dpf_line( line );

    switch( dpf_keyword ) {
    case DPF_OUTLEV:
        /*
        **  outlev
        **  Output level,
        */
        retval = sscanf( line, "%*s %d", &outlev );
        switch ( outlev ) {
        case -1:
            outputEveryNgens = (unsigned int) OUTLEV0_GENS;
            break;
        case 0:
            outputEveryNgens = (unsigned int) OUTLEV0_GENS;
            break;
        case 1:
            outputEveryNgens = (unsigned int) OUTLEV1_GENS;
            break;
        case 2:
        default:
            outputEveryNgens = (unsigned int) OUTLEV2_GENS;
            break;
        }
        break;

    case DPF_FLEXRES:
        // The DPF specifies a flexible residues file
        // -- set a flag
        // -- get the filename
        B_have_flexible_residues = TRUE;
        (void) sscanf( line, "%*s %s", FN_flexres );
        break;

    default:
        break;
    } // switch( dpf_keyword )
} // while

// Rewind DPF, so we can resume normal parsing
(void) rewind( parFile );

//______________________________________________________________________________
/*
** Output banner...
*/

banner( version_num.c_str() );

(void) fprintf(logFile, "                           $Revision: 1.114 $\n\n");
(void) fprintf(logFile, "                   Compiled on %s at %s\n\n\n", __DATE__, __TIME__);


//______________________________________________________________________________
/*
** Print the time and date when the log file was created...
*/

pr( logFile, "This file was created at:\t\t\t" );
printdate( logFile, 1 );

(void) strcpy(hostnm, "unknown host");

if (gethostname( hostnm, sizeof hostnm ) == 0) {
    pr( logFile, "                   using:\t\t\t\"%s\"\n", hostnm );
}

pr( logFile, "\nNOTE: \"rus\" stands for:\n\n      r = Real, wall-clock or elapsed time;\n      u = User or cpu-usage time;\n      s = System time\n\nAll timings are in seconds, unless otherwise stated.\n\n\n" );

//______________________________________________________________________________

(void) fprintf(logFile, "\n      ________________________________________________________________\n\n");
(void) fprintf(logFile, "                   SETTING UP DEFAULT PARAMETER LIBRARY\n");
(void) fprintf(logFile, "      ________________________________________________________________\n\n\n");

//______________________________________________________________________________
//
// Read in default parameters
//
setup_parameter_library(outlev, "default Unbound_Same_As_Bound", Unbound_Same_As_Bound);

//
// Compute the look-up table for the distance-dependent dielectric function
//
(void) fprintf(logFile, "\n\nPreparing Energy Tables for Bound Calculation:\n\n");
setup_distdepdiel(outlev, ad_energy_tables);
(void) fprintf(logFile, "Preparing Energy Tables for Unbound Calculation:\n\n");
setup_distdepdiel(outlev, unbound_energy_tables);

//______________________________________________________________________________

(void) fprintf(logFile, "\n      ___________________________________________________\n\n");
(void) fprintf(logFile,   "             PARSING INPUT DOCKING PARAMETER FILE\n");
(void) fprintf(logFile,   "      ___________________________________________________\n\n");

//______________________________________________________________________________
/*
** (Note: "dock_param_fn" set in "setflags.c"...)
*/
pr( logFile, "Docking parameter file (DPF) used for this docking:\t\t%s\n\n", dock_param_fn );

//______________________________________________________________________________
/*
** Start reading in the DPF parameter/run-control file,
*/

while( fgets(line, LINE_LEN, parFile) != NULL ) { /* PARSING-DPF parFile */
    // "line" is a string containing the current line of the input DPF.

    dpf_keyword = parse_dpf_line( line );

    switch( dpf_keyword ) {
        case -1:
            sprintf( error_message,
               "DPF> %s\n%s: ERROR: Unrecognized keyword in docking parameter file.\n",
               line, programname );
            stop( error_message );

            break;

        case DPF_NULL:
        case DPF_COMMENT:
            pr( logFile, "DPF> %s", line );
            (void) fflush(logFile);
            break;

        default:
            pr( logFile, "\n\nDPF> %s\n", line );
            indcom = strindex( line, "#" );
            if (indcom != -1) {
                line[ indcom ] = '\0'; /* Truncate "line" at the comment */
            }
            (void) fflush(logFile);
            break;
    } /* switch */

    switch( dpf_keyword ) {

//______________________________________________________________________________

    case DPF_NULL:
    case DPF_COMMENT:
        break;

//______________________________________________________________________________

    case DPF_PARAMETER_VERSION:
        /*
        ** autodock_parameter_version string
        **
        **
        **
        ** initial implementation ignores value of string
        */

        B_found_autodock_parameter_version = 1==sscanf( line, "%*s %s", autodock_parameter_version );
        pr( logFile, "\n\tAutodock parameter version %s.\n", autodock_parameter_version );
        (void) fflush(logFile);

        break;

/*____________________________________________________________________________*/

    case DPF_OUTLEV:
        /*
        **  outlev
        **  Output level,
        */
        retval = sscanf( line, "%*s %d", &outlev );
        switch ( outlev ) {
        case -1:
            pr( logFile, "Output Level = -1.  ONLY STATE VARIABLES OUTPUT, NO COORDINATES.\n" );
            outputEveryNgens = (unsigned int) OUTLEV0_GENS;
            break;
        case 0:
            pr( logFile, "Output Level = 0.\n" );
            outputEveryNgens = (unsigned int) OUTLEV0_GENS;
            break;
        case 1:
            pr( logFile, "Output Level = 1.  MINIMUM OUTPUT DURING DOCKING.\n" );
            outputEveryNgens = (unsigned int) OUTLEV1_GENS;
            break;
        case 2:
        default:
            pr( logFile, "Output Level = 2.  FULL OUTPUT DURING DOCKING.\n" );
            outputEveryNgens = (unsigned int) OUTLEV2_GENS;
            break;
        }
        if(outputEveryNgens>0) pr( logFile, "\n\tOutput population statistics every %u generations.\n", outputEveryNgens );
        else pr( logFile, "\n\tNever output generation-based population statistics.\n");
        (void) fflush(logFile);
        break;

/*____________________________________________________________________________*/

    case DPF_PARAMETER_LIBRARY:
        /*
        ** parameter_file AD4_parameters.dat
        **  or
        ** parameter_library AD4_parameters.dat
        **
        ** initial implementation based on hsearch was suggested by Mike Pique
        */

        parameter_library_found = sscanf( line, "%*s %s", FN_parameter_library );
        (void) fflush(logFile);

        read_parameter_library(FN_parameter_library, outlev);

        break;

/*____________________________________________________________________________*/

    case DPF_INCLUDE_1_4_INTERACTIONS:
        /*
         * include_1_4_interactions 0.5
         *
         * Set the Boolean variable, B_include_1_4_interactions, to TRUE.
         *
         * NOTE:  You must use this command _before_ the "move ligand.pdbqt"
         *        command, since "include_1_4_interactions" affects how the Ligand
         *        PDBQT specified by the "move" command will be interpreted.
         */
        if (B_found_move_keyword == TRUE) {
            // If we have found the move keyword already, warn the user
            // that this command ("include_1_4_interactions 0.5") should have
            // been given before this!
            pr(logFile, "\nWARNING:  This command will be ignored.\n\nYou must put this command _before_ the \"move ligand.pdbqt\" command, since this command affects how the PDBQT file will be interpreted.\n\n");
        }
        (void) sscanf( line, "%*s " FDFMT, &scale_1_4 );
        B_include_1_4_interactions = TRUE;
        print_1_4_message(logFile, B_include_1_4_interactions, scale_1_4);
        break;

//______________________________________________________________________________

    case DPF_INTELEC:
        /*
        **  intelec
        **  Calculate internal electrostatic energies...
        */
        B_calcIntElec = TRUE;
        if (outlev >= 0) {
            pr( logFile, "Electrostatic energies will be calculated for all non-bonds between moving atoms.\n\n");
        }
        retval = sscanf( line, "%*s " FDFMT, &AD3_FE_coeff_estat );
        if (retval == 1) {
            if (outlev >= 0) {
                pr(logFile, "NOTE!  Internal electrostatics will NOT be scaled by the factor specified by this command,  %.4f -- the coefficient set by this command is ignored in AutoDock 4;\n", AD3_FE_coeff_estat);
                pr(logFile, "       the coefficient that will actually be used should be set in the parameter library file.\n");
                pr(logFile, "       The coefficient for the electrostatic energy term is %.4f", AD4.coeff_estat);
                if (parameter_library_found == 1) {
                    pr( logFile, " as specified in parameter library \"%s\".\n", FN_parameter_library );
                } else {
                    pr( logFile, ", the factory default value.\n");
                }
            }
        } else {
            AD3_FE_coeff_estat = 1.0;
        }

        (void) fflush(logFile);
        break;

//______________________________________________________________________________

    case DPF_SEED:
        /*
        **  seed
        **  Set the random-number gerator's seed value,
        */
        retval = (int)sscanf( line, "%*s %s %s", param[0], param[1]);
        timeSeedIsSet[0] = 'F';
        timeSeedIsSet[1] = 'F';
        pr(logFile, "%d seed%c found.\n", retval, ((retval==1)? ' ' : 's'));
        for (j=0; j<retval; j++) {
            for (i=0; i<(int)strlen(param[j]); i++) {
                param[j][i] = (char)tolower( (int)param[j][i] );
            }
            pr(logFile, "argument \"%s\" found\n", param[j]);
        }
        if ((retval==2) || (retval==1)) {
            for (i=0; i<retval ; i++ ) {
                if (equal(param[i], "tim", 3)) {
                    timeSeedIsSet[i] = 'T';
                    seed[i] = (FourByteLong)time( &time_seed );
                    seed_random(seed[i]);
                    pr(logFile,"Random number generator was seeded with the current time, value = %ld\n",seed[i]);
                } else if (equal(param[i], "pid", 3)) {
                    timeSeedIsSet[i] = 'F';
                    seed[i] = getpid();
                    seed_random(seed[i]);
                    pr(logFile,"Random number generator was seeded with the process ID, value   = %ld\n",seed[i]);
                } else {
                    timeSeedIsSet[i] = 'F';
                    seed[i] = atol(param[i]);
                    seed_random(seed[i]);
                    pr(logFile,"Random number generator was seeded with the user-specified value  %ld\n",seed[i]);
                }
            }/*i*/
            pr(logFile, "\n");
            if (retval==2) {
                setall(seed[0], seed[1]);
                initgn(-1);  // Reinitializes the state of the current random number generator
                pr(logFile,"Portable random number generator was seeded with the user-specified values  %ld, %ld\n", seed[0], seed[1]);
            }
        } else {
            pr(logFile, "Error encountered reading seeds!\n");
        }

        (void) fflush(logFile);
        break;

/*____________________________________________________________________________*/

    case DPF_LIGAND_TYPES:
        /*
         *  Read in the ligand atom type names, e.g.
         *
         *  ligand_types C HD OA P               # ligand atom type names
         *
         *  The order of the arguments is the index that will
         *  be used for look up in the grid maps, "map_index".
         */

        //  Use the function "parsetypes" to read in the atom types;
        //
        //  The array "ligand_atom_type_ptrs" is returned, having been filled with pointers
        //  to the beginning of each "atom type word" (not atom type characters);
        //  In AutoDock 4, an atom type can be either 1 or 2 characters long.
        //
        num_atom_types = parsetypes(line, ligand_atom_type_ptrs, MAX_ATOM_TYPES);
        if (num_atom_types<0){
            prStr( error_message, "%s:  ERROR! Too many atom types have been found: maximum is %d; we cannot continue !\n\n", programname, MAX_ATOM_TYPES );
            pr_2x( logFile, stderr, error_message );
            exit(-1);
        }



        B_found_ligand_types = TRUE;

        // This is not necessary if we increment num_maps one-at-a-time as we read each atom map in
        // num_maps += num_atom_types;
        info->num_atom_types = num_atom_types;

        for (i=0; i<num_atom_types; i++) {
            strcpy(info->atom_type_name[i], ligand_atom_type_ptrs[i]);
#ifdef DEBUG
            (void) fprintf(logFile, "%d %s ->%s\n",i, ligand_atom_type_ptrs[i], info->atom_type_name[i]);
#endif
        }

        if (num_atom_types > 0) {
            B_atom_types_found = TRUE;
        } else {
            prStr( error_message, "%s:  ERROR! No atom types have been found; we cannot continue without this information!\n\n", programname );
            pr_2x( logFile, stderr, error_message );
            prStr( error_message, "%s:  ERROR! Are you trying to use an AutoDock 3 DPF with AutoDock 4?\n\n", programname );
            pr_2x( logFile, stderr, error_message );
            exit(-1);
        }

        if (debug > 0) {
            for (i=0; i<num_atom_types; i++) {
                (void) fprintf(logFile, "info->atom_type_name[%d] = %s\n", i, info->atom_type_name[i] );
            }
        }

        // For all ligand atom types... set up the map_index
        // "ligand_types"
        for (i=0; i<num_atom_types; i++) {
            foundParameter = apm_find(info->atom_type_name[i]);
            if (foundParameter != NULL ) {
                // Not NULL means we have found this atom type's parameters.
                // Set the ParameterEntry's "map_index" member to the
                // 0-based index it had in the list of ligand types supplied in the DPF "types" line:
                foundParameter->map_index = i;
                parameterArray[i] = *(foundParameter);
                if (outlev > 0) {
                    (void) fprintf( logFile, "Parameters found for ligand type \"%s\" (grid map index = %d, weighted well depth, epsilon = %6.4f)", foundParameter->autogrid_type, foundParameter->map_index, foundParameter->epsij );
                    if (parameter_library_found == 1) {
                        pr( logFile, " in parameter library \"%s\".\n", FN_parameter_library );
                    } else {
                        pr( logFile, "\n");
                    }
                }
            } else {
                // We could not find this parameter -- return error here
                prStr( error_message,"%s: ERROR:  Unknown ligand atom type \"%s\"; add parameters for it to the parameter library first!\n", programname, info->atom_type_name[i]);
                pr_2x( logFile, stderr, error_message );
                if (parameter_library_found == 1) {
                    prStr( error_message,"%s:         Edit the parameter library file \"%s\" and try again.\n", programname, FN_parameter_library );
                    pr_2x( logFile, stderr, error_message );
                }
                exit(-1);
            } // if / else apm_find
        } // for i
        pr( logFile, "\n\n");

        (void) fflush( logFile);

        // Calculate the internal energy table

        // loop over atom types, i, from 1 to number of atom types
        for (i=0; i<num_atom_types; i++) {

            //  Find internal energy parameters, i.e.  epsilon and r-equilibrium values...
            //  Lennard-Jones and Hydrogen Bond Potentials

            Ri = parameterArray[i].Rij;
            epsi = parameterArray[i].epsij;
            Ri_hb = parameterArray[i].Rij_hb;
            epsi_hb = parameterArray[i].epsij_hb;
            hbondi = parameterArray[i].hbond;

            // loop over atom types, j, from i to number of atom types
            for (j=i; j<num_atom_types; j++) {

                //  Find internal energy parameters, i.e.  epsilon and r-equilibrium values...
                //  Lennard-Jones and Hydrogen Bond Potentials

                Rj = parameterArray[j].Rij;
                epsj = parameterArray[j].epsij;
                Rj_hb = parameterArray[j].Rij_hb;
                epsj_hb = parameterArray[j].epsij_hb;
                hbondj = parameterArray[j].hbond;

                // we need to determine the correct xA and xB exponents
                xA = 12; // for both LJ, 12-6 and HB, 12-10, xA is 12
                xB =  6; // assume we have LJ, 12-6

                if ( ((hbondi == DS) || (hbondi == D1)) && ((hbondj == AS) || (hbondj == A1) || (hbondj == A2)) ) {
                    // i is a donor and j is an acceptor.
                    // i is a hydrogen, j is a heteroatom
                    // we need to calculate the arithmetic mean of Ri_hb and Rj_hb  // not in this Universe...  :-(
                    //Rij = arithmetic_mean(Ri_hb, Rj_hb);// not in this Universe...  :-(
                    Rij = Rj_hb;
                    // we need to calculate the geometric mean of epsi_hb and epsj_hb  // not in this Universe...  :-(
                    //epsij = geometric_mean(epsi_hb, epsj_hb);// not in this Universe...  :-(
                    epsij = epsj_hb;
                    xB = 10;
                } else if ( ((hbondi == AS) || (hbondi == A1) || (hbondi == A2)) && ((hbondj == DS) || (hbondj == D1))) {
                    // i is an acceptor and j is a donor.
                    // i is a heteroatom, j is a hydrogen
                    // we need to calculate the arithmetic mean of Ri_hb and Rj_hb// not in this Universe...  :-(
                    //Rij = arithmetic_mean(Ri_hb, Rj_hb);// not in this Universe...  :-(
                    Rij = Ri_hb;
                    // we need to calculate the geometric mean of epsi_hb and epsj_hb// not in this Universe...  :-(
                    //epsij = geometric_mean(epsi_hb, epsj_hb);// not in this Universe...  :-(
                    epsij = epsi_hb;
                    xB = 10;
                } else {
                    // we need to calculate the arithmetic mean of Ri and Rj
                    Rij = arithmetic_mean(Ri, Rj);
                    // we need to calculate the geometric mean of epsi and epsj
                    epsij = geometric_mean(epsi, epsj);
                }

                /* Check that the Rij is reasonable */
                if ((Rij < RIJ_MIN) || (Rij > RIJ_MAX)) {
                    (void) fprintf( logFile,
                    "WARNING: pairwise distance, Rij, %.2f, is not a very reasonable value for the equilibrium separation of two atoms! (%.2f Angstroms <= Rij <= %.2f Angstroms)\n\n", Rij, RIJ_MIN, RIJ_MAX);
                    (void) fprintf( logFile, "Perhaps you meant to use \"intnbp_coeffs\" instead of \"intnbp_r_eps\"?\n\n");
                    /* gmm commented out for dave goodsell, mutable atoms
                     * exit(-1); */
                }
                /* Check that the epsij is reasonable */
                if ((epsij < EPSIJ_MIN) || (epsij > EPSIJ_MAX)) {
                    (void) fprintf( logFile,
                    "WARNING: well-depth, epsilon_ij, %.2f, is not a very reasonable value for the equilibrium potential energy of two atoms! (%.2f kcal/mol <= epsilon_ij <= %.2f kcal/mol)\n\n", epsij, EPSIJ_MIN, EPSIJ_MAX);
                    (void) fprintf( logFile, "Perhaps you meant to use \"intnbp_coeffs\" instead of \"intnbp_r_eps\"?\n\n");
                    /* gmm commented out for dave goodsell, mutable atoms
                     * exit(-1); */
                }
                /* Defend against division by zero... */
                if (xA != xB) {
                    cA = (tmpconst = epsij / (Real)(xA - xB)) * pow( (double)Rij, (double)xA ) * (Real)xB;
                    cB = tmpconst * pow( (double)Rij, (double)xB ) * (Real)xA;
                    pr(logFile, "\nCalculating internal non-bonded interaction energies for docking calculation;\n");
                    intnbtable( &B_havenbp, a1, a2, info, cA, cB, xA, xB, AD4.coeff_desolv, sigma, ad_energy_tables, BOUND_CALCULATION );
                    pr(logFile, "\nCalculating internal non-bonded interaction energies for unbound conformation calculation;\n");
                    intnbtable( &B_havenbp, a1, a2, info, cA_unbound, cB_unbound, xA_unbound, xB_unbound, AD4.coeff_desolv, sigma, unbound_energy_tables, UNBOUND_CALCULATION );
                    // Increment the atom type numbers, a1 and a2, for the internal non-bond table
                    a2++;
                    if (a2 >= info->num_atom_types) {
                        a1++;
                        a2 = a1;
                    }

                } else {
                    pr(logFile,"WARNING: Exponents must be different, to avoid division by zero!\n\tAborting...\n");
                    exit(-1);
                }
                (void) fflush(logFile);

            } // for j
        } // for i
        break;

//______________________________________________________________________________

    case DPF_FLD:
        /*
        ** fld
        ** GRID_DATA_FILE
        ** Read the (AVS-format) grid data file, .fld
        */
        // TO DO: add outlev
        readfield( info, line, jobStart, tms_jobStart );
        num_maps = 0;

        /*
        // Dynamically allocate memory for the maps
        map = NewGridMapSet(info);

        if (map == NULL) {
            prStr(error_message, "%s:  Sorry, there is not enough memory to store the grid maps.  Please use smaller maps and/or fewer atom types.\n", programname);
            stop(error_message);
            exit(1);
        }
        // Initialise the maps
        for (i=0; i<num_map_values; i++) {
            map[i] = 0.0L;
        }
        */
        (void) fflush(logFile);
        break;

//______________________________________________________________________________

    case DPF_MAP:
        /*
        ** map
        ** ATOMIC AFFINITY, ELECTROSTATIC POTENTIAL OR DESOLVATION ENERGY GRID MAP
        ** Read in active site grid map...
        */
        B_charMap = FALSE;
        if (B_atom_types_found == TRUE) {
            // Read in the AutoGrid atomic affinity map
            // map_index could be incremented here if we had the atom_type stored in each map...
            map_stats = readmap( line, outlev, jobStart, tms_jobStart, B_charMap, &B_havemap, num_maps, info, map, 'a' );
            pr(logFile, "Min= %.3lf Mean= %.3lf Max= %.3lf\n\n",
                    map_stats.minimum, map_stats.mean, map_stats.maximum);
            num_maps++;

        } else {
            prStr( error_message, "%s:  ERROR! No atom types have been found; we cannot continue without this information!\n\n", programname );
            pr_2x( logFile, stderr, error_message );
            prStr( error_message, "%s:  ERROR! Are you trying to use an AutoDock 3 DPF with AutoDock 4?\n\n", programname );
            pr_2x( logFile, stderr, error_message );
            exit(-1);
        }
        (void) fflush(logFile);
        break;

//______________________________________________________________________________

    case DPF_ELECMAP:
        /*
         *  elecmap file.e.map
         */
        map_stats = readmap( line, outlev, jobStart, tms_jobStart, B_charMap, &B_havemap, num_maps, info, map, 'e' );
        pr(logFile, "Min= %.3lf Mean= %.3lf Max= %.3lf\n\n",
                map_stats.minimum, map_stats.mean, map_stats.maximum);
        ElecMap = num_maps;
        B_found_elecmap = TRUE;
        num_maps++;
        (void) fflush(logFile);
        break;

//______________________________________________________________________________

    case DPF_DESOLVMAP:
        /*
         *  desolvmap file.d.map
         */
        map_stats = readmap( line, outlev, jobStart, tms_jobStart, B_charMap, &B_havemap, num_maps, info, map, 'd' );
        pr(logFile, "Min= %.3lf Mean= %.3lf Max= %.3lf\n\n",
                map_stats.minimum, map_stats.mean, map_stats.maximum);
        DesolvMap = num_maps;
        B_found_desolvmap = TRUE;
        num_maps++;
        (void) fflush(logFile);
        break;

//______________________________________________________________________________

    case DPF_CHARMAP:
        /*
        ** charmap
        ** ATOMIC AFFINITY, ELECTROSTATIC POTENTIAL OR DESOLVATION ENERGY GRID MAP
        ** Read in active site grid map...
        */
        B_charMap = TRUE;
        if (B_atom_types_found == TRUE) {
            // map_index could be incremented here if we had the atom_type stored in each map...
            map_stats = readmap( line, outlev, jobStart, tms_jobStart, B_charMap, &B_havemap, num_maps, info, map, 'c' );
            pr(logFile, "Min= %.3lf Mean= %.3lf Max= %.3lf\n\n",
                    map_stats.minimum, map_stats.mean, map_stats.maximum);
            num_maps++;
        } else {
            prStr( error_message, "%s:  ERROR! No atom types have been found; we cannot continue without this information!\n\n", programname );
            pr_2x( logFile, stderr, error_message );
            prStr( error_message, "%s:  ERROR! Are you trying to use an AutoDock 3 DPF with AutoDock 4?\n\n", programname );
            pr_2x( logFile, stderr, error_message );
            exit(-1);
        }
        (void) fflush(logFile);
        break;

//______________________________________________________________________________

    case DPF_MOVE:
        /*
        ** move ligand_file.pdbqt
        ** Specify the movable ligand,
        */
        //
        // Initialisations that must be done before reading in a new ligand...
        //
        if (num_maps != num_atom_types + NUM_NON_VDW_MAPS) { //  dsolv map and elec map
            prStr(error_message, "\n\nMISSING MAP ERROR:\nnumber of maps %d does not match number expected for %d ligand types. \nUnable to continue.\n", num_maps, num_atom_types);
            stop(error_message);
            exit(-1);
        }

        nconf = 0;
        for (k = 0; k < MAX_RUNS; k++) {
            for (i = 0; i  < MAX_TORS;  i++ ) {
                sHist[k].tor[i] = 0.0;
            }
            econf[k] = 0.0;
        }
        for (j = 0;  j < MAX_ATOMS;  j++ ) {
            type[j] = 0;
            ignore_inter[j] = 0;
        }
        for (i = 0; i  < MAX_TORS;  i++ ) {
            for (j = 0;  j < MAX_ATOMS;  j++ ) {
                tlist[i][j] = 0;
            }
            B_isTorConstrained[i] = 0;
            US_torProfile[i][0] = 0;
            N_con[i] = 0;
        }
        for (j = 0; j < MAX_NONBONDS; j++) {
            nonbondlist[j].a1 = nonbondlist[j].a2 = 0;
        }
        for (j=0; j<3; j++) {
            Nnb_array[j] = 0;
            nb_group_energy[j] = 0.0;
        }
        initialiseState( &sInit );
        initialiseState( &(ligand.S) );
        initialiseQuat( &q_reorient );
        B_constrain_dist = B_haveCharges = FALSE;
        ntor1 = ntor = atomC1 = atomC2 = 0;
        ntor_ligand = 0;
        ntorsdof = 0;
        sqlower = squpper = 0.0;
        strcpy( FN_pop_file, "");  // means don't print pop_file
        Nnb = 0;
        ligand_is_inhibitor = 1;
        initialise_energy_breakdown(&eb, 0, 0);
        //
        // end of initialization
        //

        // this is the DPF_MOVE section...
        B_found_move_keyword = TRUE;
        B_found_about_keyword = FALSE; //set false by 'move' true by 'about'

        print_1_4_message(logFile, B_include_1_4_interactions, scale_1_4);

        natom=0;
        ligand = readPDBQT( line,
                            num_atom_types,
                            &natom,
                            crdpdb, crdreo, charge, &B_haveCharges,
                            type, bond_index,
                            pdbaname, FN_ligand, FN_flexres, B_have_flexible_residues, atomstuff, &n_heavy_atoms_in_ligand,
                            &B_constrain_dist, &atomC1, &atomC2,
                            &sqlower, &squpper,
                            &ntor1, &ntor, &ntor_ligand,
                            tlist, vt,
                            &Nnb, nonbondlist,
                            jobStart, tms_jobStart, hostnm, &ntorsdof, outlev,
                            ignore_inter,
                            B_include_1_4_interactions,
                            atoms, PDBQT_record, end_of_branch );

        // pre-calculate some values we will need later in computing the desolvation energy
        //
        for (i=0; i<natom; i++) {
            abs_charge[i] = fabs(charge[i]);
            qsp_abs_charge[i] = qsolpar * abs_charge[i];
        }
        pr(logFile, "Number of atoms in ligand:  %d\n\n", true_ligand_atoms);
        pr(logFile, "Number of non-hydrogen atoms in ligand:  %d\n\n", n_heavy_atoms_in_ligand);

        pr(logFile, "Number of vibrational degrees of freedom of ligand:  %d\n\n\n", (3 * true_ligand_atoms) - 6 );
        pr( logFile, "Number of torsional degrees of freedom = %d\n", ntorsdof);

        torsFreeEnergy = (Real)ntorsdof * AD4.coeff_tors;

        pr( logFile, "Estimated loss of torsional free energy upon binding = %+.4f kcal/mol\n\n", torsFreeEnergy);

        for (i=0;i<natom;i++) {
            if (ignore_inter[i] == 1) {
                pr(logFile, "Special Boundary Conditions:\n");
                pr(logFile, "____________________________\n\n");
                pr(logFile, "AutoDock will ignore the following atoms in the input PDBQT file \nin intermolecular energy calculations:\n");
                pr(logFile, "\n(This is because these residue atoms are at the boundary between \nflexible and rigid, and since they cannot move \nthey will not affect the total energy.)\n\n");
                break;
            }
        }
        for (i=0;i<natom;i++) {
            if (ignore_inter[i] == 1) {
                pr(logFile, "Atom number %d:  %s\n", i+1, atomstuff[i] );
            }
        }
        pr(logFile, "\n");

        if (!B_haveCharges) {
            pr( logFile, "%s: WARNING! No partial atomic charges have been supplied yet.\n\n",programname);
        } else {
            if (Nnb > 0) {
            pr(logFile,"Calculating the product of the partial atomic charges, q1*q2, for all %d non-bonded pairs...\n\n", Nnb);
            pr(logFile," -- Scaled q1*q2 means multiplied by both  %.1lf (for conversion later on to kcal/mol)\n", (double)ELECSCALE);
            pr(logFile,"    and by the AD4 FF electrostatics coefficient, %.4lf\n\n", (double)AD4.coeff_estat);
            if (outlev >= 0) {
                pr(logFile,"Non-bonded                           Scaled\n");
                pr(logFile,"   Pair     Atom1-Atom2    q1*q2      q1*q2\n");
                pr(logFile,"__________  ___________  _________  _________\n");
                for (i = 0;  i < Nnb;  i++) {
                    atm1 = nonbondlist[i].a1;
                    atm2 = nonbondlist[i].a2;
                    int t1 = nonbondlist[i].t1;
                    int t2 = nonbondlist[i].t2;
                    nonbondlist[i].desolv =
                           ( parameterArray[t2].vol * (parameterArray[t1].solpar + qsp_abs_charge[atm1])
                           + parameterArray[t1].vol * (parameterArray[t2].solpar + qsp_abs_charge[atm2]) );
                    nonbondlist[i].q1q2 = charge[atm1] * charge[atm2];
                    pr(logFile,"   %4d     %5d-%-5d    %5.2f",i+1,atm1+1,atm2+1,nonbondlist[i].q1q2);
                    nonbondlist[i].q1q2 *= ELECSCALE * AD4.coeff_estat;
                    pr(logFile,"     %6.2f\n",nonbondlist[i].q1q2);
                }
                pr(logFile,"\n");
                } // if outlev > 0
            } // if NNb > 0
        } // else

        sInit.ntor = ligand.S.ntor;
        ++nmol;
        ++nlig;

        (void) fflush(logFile);
        break;

/*____________________________________________________________________________*/

    case DPF_FLEXRES:
        /*
         * flexible_residues file.pdbqt
         */
        (void) fflush(logFile);
        pr(logFile, "\nThe flexible residues will be read in from \"%s\".\n", FN_flexres);
        break;


#ifdef USING_COLINY
/*____________________________________________________________________________*/

    case DPF_COLINY:
    {
        //ostdiostream fstr(logFile);
        //ostdiostream fstr(logFile->_file);
        //CommonIO::set_streams(&fstr,&fstr,&cin);

        struct tms tms_colinyStart;
        struct tms tms_colinyEnd;

        Clock  colinyStart;
        Clock  colinyEnd;

        int coliny_seed;
        char algname[LINE_LEN];
        char nruns_str[LINE_LEN];
        (void) sscanf(line, "%*s %s %d", algname, &nruns);
        (void) sscanf(line, "%*s %s %s", algname, nruns_str);

        if (strcmp(algname,"help")==0) {
            std::vector<double> initvec;
            coliny_init(algname, "", 0);
            prStr(error_message, "%s:  ERROR:  no optimizer type specified.", programname);
            stop(error_message);
            exit(-1);
        }
        else if (strcmp(nruns_str,"help")==0) {
            std::vector<double> initvec;
            coliny_init(algname, nruns_str, 0);
            prStr(error_message, "%s:  ERROR:  no optimizer type specified.", programname);
            stop(error_message);
            exit(-1);
        }


            if (nruns>MAX_RUNS) {
                prStr(error_message, "%s:  ERROR:  %d runs requested, but only dimensioned for %d.\nChange \"MAX_RUNS\" in \"constants.h\".", programname, nruns, MAX_RUNS);
                stop(error_message);
                exit(-1);
            }
            exit_if_missing_elecmap_desolvmap_about("coliny");

            evaluate.setup( crd, charge, abs_charge, qsp_abs_charge, type, natom,
                            map, elec, emap, nonbondlist, ad_energy_tables,
                            Nnb, B_calcIntElec, B_isGaussTorCon, B_isTorConstrained, B_ShowTorE,
                            US_TorE, US_torProfile,
                            vt, tlist,
                            crdpdb, crdreo, sInit, ligand, ignore_inter, B_include_1_4_interactions, scale_1_4,
                            unbound_internal_FE, info,
                            B_use_non_bond_cutoff, B_have_flexible_residues);

            evaluate.compute_intermol_energy(TRUE);

            char domain[1024];
            // NOTE: Coliny enforces the bound constraints, but since the
            // torsion angles are periodic, we simply prevent the optimizer
            // from going too far.
            if (sInit.ntor > 0) {
                sprintf(domain,"[%f,%f] [%f,%f] [%f,%f] [-1000.0,1000.0]^3 [-3.1416,3.1416] [-3.1416,3.1416]^%d",(double)info->lo[X], (double)info->hi[X], (double)info->lo[Y], (double)info->hi[Y], (double)info->lo[Z], (double)info->hi[Z], sInit.ntor);
                //sprintf(domain,"[%f,%f] [%f,%f] [%f,%f] [-1.0,1.1]^3 [-6.2832,12.5664] [-6.2832,12.5664]^%d",(double)info->lo[X], (double)info->hi[X], (double)info->lo[Y], (double)info->hi[Y], (double)info->lo[Z], (double)info->hi[Z], sInit.ntor);
            } else {
                sprintf(domain,"[%f,%f] [%f,%f] [%f,%f] [-1000.0,1000.0]^3 [-3.1416,3.1416]",(double)info->lo[X], (double)info->hi[X], (double)info->lo[Y], (double)info->hi[Y], (double)info->lo[Z], (double)info->hi[Z]);
            }
            pr(logFile, "Number of Coliny %s dockings = %d run%c\n", algname, nruns, (nruns>1)?'s':' ');
            pr(logFile, "Search Domain: %s\n", domain);

            //
            // COLINY-SPECIFIC LOGIC - BEGIN
            //

            try {

                std::vector<double> initvec, finalpt;
                // set up initial point
                initvec.resize(7+sInit.ntor);
                initvec[0] = sInit.T.x;
                initvec[1] = sInit.T.y;
                initvec[2] = sInit.T.z;
                /*
                 * axis-angle (nx,ny,nz,ang) suffers from bias
                initvec[3] = sInit.Q.nx;
                initvec[4] = sInit.Q.ny;
                initvec[5] = sInit.Q.nz;
                initvec[6] = DegreesToRadians( sInit.Q.ang );
                */
                sInit.Q = convertRotToQuat( sInit.Q );
                initvec[3] = sInit.Q.x;
                initvec[4] = sInit.Q.y;
                initvec[5] = sInit.Q.z;
                initvec[6] = sInit.Q.w;
                for (j=0; j < sInit.ntor ; j++) {
                  initvec[j+7] = DegreesToRadians(sInit.tor[j]);
                }
                coliny_init(algname, domain, sInit.ntor+7);

                for (j=0; j<nruns; j++) {
                  fprintf( logFile, "\n\n\tBEGINNING Coliny %s DOCKING\n",algname);
                  pr(logFile, "\nDoing %s run:  %d/%d.\n", algname, j+1, nruns);

                  //coliny uses a single seed
                  coliny_seed = seed[0]+seed[1]+j;
                  pr(logFile, "Seed: %d [%ld+%ld+%d]\n", coliny_seed, seed[0], seed[1], j);
                  //pr(logFile, "Seeds:  %ld %ld\n", seed[0], seed[1]);
                  (void) fflush(logFile);

                  colinyStart = times(&tms_colinyStart);

                  finalpt.resize( initvec.size() );
                  int neval, niters;
                  coliny_minimize( coliny_seed, initvec, finalpt, neval, niters );
                  //fstr.flush();

                  make_state_from_rep( (double *)&(finalpt[0]), int(finalpt.size()), &sHist[nconf]);

                  pr(logFile, "\nTotal Num Evals: %d\n", neval);
                  printState(logFile, sHist[nconf], 2);

                  colinyEnd = times(&tms_colinyEnd);
                  pr(logFile, "Time taken for this %s run:\n", algname);
                  timesyshms(colinyEnd-colinyStart, &tms_colinyStart, &tms_colinyEnd);
                  pr(logFile, "\n");

                  pr(logFile, "Total number of Energy Evaluations: %d\n", (int)evaluate.evals() );
                  //pr(logFile, "Total number of Iterations:        %d\n", (int)niters);

                  pr(logFile, "\nFinal docked state:\n");
                  pr( logFile, UnderLine );
                  pr( logFile, "\n\n\tFINAL Coliny %s DOCKED STATE\n",algname );
                  pr( logFile,     "\t____________________________________\n\n\n" );
                  (void) fflush(logFile);

                  writePDBQT( j, seed, FN_ligand, dock_param_fn, lig_center,
                              sHist[nconf], ntor, &eintra, &einter, natom, atomstuff,
                              crd, emap, elec,
                              charge, abs_charge, qsp_abs_charge,
                              ligand_is_inhibitor,
                              torsFreeEnergy,
                              vt, tlist, crdpdb, nonbondlist,
                              ad_energy_tables,
                              type, Nnb, B_calcIntElec,
                              map,
                              outlev,
                              ignore_inter,
                              B_include_1_4_interactions, scale_1_4, parameterArray, unbound_internal_FE,
                              info, DOCKED, PDBQT_record, B_use_non_bond_cutoff, B_have_flexible_residues, ad4_unbound_model);

                  // See also "calculateEnergies.cc", switch(ad4_unbound_model)
                  if (ad4_unbound_model == Unbound_Same_As_Bound) {
                      // Update the unbound internal energy, setting it to the current internal energy
                      unbound_internal_FE = eintra;
                  }
                  econf[nconf] = eintra + einter + torsFreeEnergy - unbound_internal_FE;
                  evaluate.reset();

                  ++nconf;

                } // Next run
                if(write_stateFile){
                  fprintf(stateFile,"\t</runs>\n");
                  (void) fflush(stateFile);
                }
                (void) fflush(logFile);
            }
            catch (std::exception& err) {
              (void)fprintf(logFile, "Caught Exception: %s\n", err.what());
              exit(1);
            }

    }
    break;
#endif


//______________________________________________________________________________

    case DPF_ABOUT:
        /*
        **  about
        **  Rotation center for current ligand,
        */
        (void) sscanf( line, "%*s " FDFMT3, &lig_center[X], &lig_center[Y], &lig_center[Z]);
        pr( logFile, "Small molecule center of rotation =\t" );
        pr( logFile, "(%+.3f, %+.3f, %+.3f)\n\n", lig_center[X], lig_center[Y], lig_center[Z]);
        B_found_about_keyword = TRUE; //set false by 'move' true by 'about'
        /*
        **  Center the ligand,
        */
        if ( nmol == 0 ) {
            pr( logFile, "Must specify a ligand PDBQT file, using the \"move\" command.\n");
        } else {
            if (outlev >= 0) {
                pr( logFile, "Translating small molecule by:\t" );
                pr( logFile, "(%+.3f, %+.3f, %+.3f)\n\n", -lig_center[X], -lig_center[Y], -lig_center[Z]);
            }
            /*
            **  Zero-out on central point...
            */
            maxrad = -1.0;
            for ( i=0; i<true_ligand_atoms; i++ ) { /*new, gmm, 6-23-1998*/
                r2sum=0.0;
                for (xyz = 0;  xyz < SPACE;  xyz++) {
                    c = crd[i][xyz] = (crdpdb[i][xyz] -= lig_center[xyz]);
                    r2sum += c*c;
                } /* xyz */
                maxrad = max(maxrad,sqrt(r2sum));
            } /* i */
            if (outlev >= 0) {
                pr( logFile, "Furthest ligand atom from \"about\" center is %.3f Angstroms (maxrad).\n\n",maxrad);
            }
        }
        (void) fflush(logFile);
        break;

/*____________________________________________________________________________*/

    case DPF_REORIENT:
        /*
         *  reorient random
         *      # applies a random rotation to the input ligand
         * -OR-
         *  reorient standard
         *      # moves the ligand such that
         *      # the first three atoms lie parallel to the xy-plane, and
         *      # the first two atoms lie parallel to the x-axis
         * -OR-
         *  reorient <axis-x> <axis-y> <axis-z> <angle>
         *      # applies the specified rotation to the input ligand
         */
        (void) sscanf( line, "%*s %s", param[0] );
        { // Parse the reorient command
            for (i=0; i<6; i++) {
                param[0][i] = (char)tolower( (int)param[0][i] );
            }
            if (equal(param[0],"random",6)) {
                // reorient random
                B_reorient_random = TRUE; // create a new random orientation before docking

                create_random_orientation( &q_reorient );

            } else if (equal(param[0],"standard",8)) {
                { // reorient standard
                B_reorient_random = FALSE; // do not create a new random orientation before docking

                if (true_ligand_atoms >= 3 ) {
                    // Move the ligand such that
                    // the first three atoms lie parallel to the xy-plane, and
                    // the first two atoms lie parallel to the x-axis
                    Vector vec_01,     // vector between ligand atoms 0 and 1
                           vec_12,     // vector between ligand atoms 1 and 2
                           vec_normal, // vector perpendicular to plane of vec_01 and vec_12
                           vec_x_axis, // vector along the X-axis
                           vec_z_axis, // vector along the Z-axis
                           vec_reorient_axis; // vector describing the axis about which to reorient
                    // Set the X and Z axes:
                    vec_x_axis[X] = 1.;
                    vec_x_axis[Y] = 0.;
                    vec_x_axis[Z] = 0.;
                    vec_z_axis[X] = 0.;
                    vec_z_axis[Y] = 0.;
                    vec_z_axis[Z] = 1.;
                    for (xyz = 0;  xyz < SPACE;  xyz++) {
                        vec_01[xyz] = (double)( crdpdb[1][xyz] - crdpdb[0][xyz] );
                        vec_12[xyz] = (double)( crdpdb[2][xyz] - crdpdb[1][xyz] );
                    }
                    // Compute the normal to vec_01 and vec_12
                    Cross_product( vec_normal, vec_01, vec_12 );
                    Print_vector( logFile, "vec_01", vec_01 );
                    Print_vector( logFile, "vec_12", vec_12 );
                    Print_vector( logFile, "vec_normal", vec_normal );
                    Print_vector( logFile, "vec_z_axis", vec_z_axis );
                    // Compute the angle between vec_01 and vec_12
                    double angle_012 = 0.;
                    angle_012 = Angle_between( vec_01, vec_12 );
                    pr( logFile, "Angle between vectors 01 and 12 = %.2f degrees\n", RadiansToDegrees( angle_012 ) );
                    if ( ( fabs(angle_012) < APPROX_ZERO ) || ( ( fabs(angle_012) > (PI - APPROX_ZERO) ) && ( fabs(angle_012) < (PI + APPROX_ZERO) ) ) ) {
                        // angle is too small or "too linear" to align the molecule into the xy-plane
                        pr( logFile, "%s:  WARNING!  The angle between the first three atoms is not suitable (%6.3f degrees) to align them with the xy-plane.\n", programname, RadiansToDegrees( angle_012 ) );
                    } else {
                        // Calculate angle between vec_normal and the z-axis
                        double angle_n1z = 0.;  // Angle between vec_normal and the z-axis
                        angle_n1z = Angle_between( vec_normal, vec_z_axis );
                        pr( logFile, "Angle between vec_normal and vec_z_axis = %.2f degrees\n", RadiansToDegrees( angle_n1z ) );
                        //
                        // We need to rotate the molecule about the normal to vec_normal and vec_z_axis
                        Cross_product( vec_reorient_axis, vec_normal, vec_z_axis );
                        //
                        // Set the rotation axis for reorientation
                        q_reorient.nx = vec_reorient_axis[X];
                        q_reorient.ny = vec_reorient_axis[Y];
                        q_reorient.nz = vec_reorient_axis[Z];
                        //
                        // Normalise the vector defining the axis of rotation:
                        q_reorient = normRot( q_reorient );
                        //
                        // Set the angle for reorientation of the first 3 atoms
                        // into the xy-plane
                        q_reorient.ang = -angle_n1z;
                        //
                        // Convert the rotation-about-axis components (nx,ny,nz,ang)
                        // to a rotation-quaternion (x,y,z,w):
                        q_reorient = convertRotToQuat( q_reorient );

                        // Rotate ligand into the xy-plane...
                        // qtransform( origin, q_reorient, crdreo, true_ligand_atoms );
                        qtransform( origin, q_reorient, crdpdb, true_ligand_atoms );

                        // Compute the updated vec_01, the vector between atom 0 and atom 1,
                        // since the preceding "qtransform" changed the coordinates.
                        for (xyz = 0;  xyz < SPACE;  xyz++) {
                            // vec_01[xyz] = (double)( crdreo[1][xyz] - crdreo[0][xyz] );
                            vec_01[xyz] = (double)( crdpdb[1][xyz] - crdpdb[0][xyz] );
                        }
                        //
                        // Compute the angle between vec_01 and the x-axis:
                        double angle_01x = 0.;
                        angle_01x = Angle_between( vec_01, vec_x_axis );
                        //
                        pr( logFile, "Angle between vector 01 and the x-axis = %.2f degrees\n", RadiansToDegrees( angle_01x ) );
                        //
                        // The rotation axis to rotate the first two atoms, 0 and 1,
                        // to be parallel to the x-axis, will be
                        // perpendicular to the xy-plane, i.e. the z-axis,
                        // since the molecule's first 3 atoms are now in the xy-plane.
                        q_reorient.nx = vec_z_axis[X];
                        q_reorient.ny = vec_z_axis[Y];
                        q_reorient.nz = vec_z_axis[Z];
                        //
                        // Set the rotation angle:
                        q_reorient.ang = angle_01x;
                        //
                        // Build the quaternion from the axis-angle rotation values:
                        q_reorient = convertRotToQuat( q_reorient );
                    } // angle_012 is appropriate to align into xy-plane

                } else {
                    prStr( error_message, "%s: ERROR! Insufficient atoms in the ligand.  There must be at least three atoms in the ligand to use this command.\n", programname );
                    stop( error_message );
                    exit( -1 );
                }
                } // reorient standard
            } else {
                { // reorient <nx> <ny> <nz> <angle>
                    B_reorient_random = FALSE; // do not create a new random orientation before docking

                    // Read the specified initial orientation for the ligand
                    retval = (int)sscanf( line,"%*s %lf %lf %lf %lf", &(q_reorient.nx), &(q_reorient.ny), &(q_reorient.nz), &(q_reorient.ang) );
                    if ( retval == 4 ) {
                        // Normalise the vector defining the axis of rotation:
                        q_reorient = normRot( q_reorient );
                        // Make sure angle is in radians, and ranges from -PI to PI
                        q_reorient.ang = DegreesToRadians( q_reorient.ang ); // convert from degrees to radians
                        q_reorient.ang = ModRad( q_reorient.ang ); // wrap to range (0, 2*PI) using modulo 2*PI
                        q_reorient.ang = WrpRad( q_reorient.ang ); // wrap to range (-PI, PI)
                        pr( logFile, "After normalising the vector, and converting the angle to radians, the axis-angle rotation becomes ((%.3f, %.3f, %.3f), %.2f radians)\n",
                                q_reorient.nx, q_reorient.ny, q_reorient.ny, q_reorient.ang);
                        // Convert the rotation-about-axis components (nx,ny,nz,ang)
                        // to a rotation-quaternion (x,y,z,w):
                        q_reorient = convertRotToQuat( q_reorient );
                    } else {
                        prStr( error_message, "%s: ERROR! Please specify the vector and rotation angle using four real numbers.\n", programname );
                        stop( error_message );
                        exit( -1 );
                    }
                } // reorient <nx> <ny> <nz> <angle>
            } // endif
        } // end parsing reorient command line

        // reorient( logFile, true_ligand_atoms, atomstuff, crdreo, charge, type,
        reorient( logFile, true_ligand_atoms, atomstuff, crdpdb, charge, type,
                  parameterArray, q_reorient, origin, ntor, tlist, vt, &ligand, debug );

        (void) fflush(logFile);
        break;


//______________________________________________________________________________

    case DPF_TRAN0:
        /*
        **  tran0
        **  Initial_translation,
        */
        (void) sscanf( line, "%*s %s", param[0]);
        for (i=0; i<6; i++) {
            param[0][i] = (char)tolower( (int)param[0][i] );
        }
        if (equal(param[0],"random",6)) {
            B_RandomTran0 = TRUE;
            ligand.S.T.x = sInit.T.x = random_range( info->lo[X], info->hi[X] );
            ligand.S.T.y = sInit.T.y = random_range( info->lo[Y], info->hi[Y] );
            ligand.S.T.z = sInit.T.z = random_range( info->lo[Z], info->hi[Z] );
        } else {
            B_RandomTran0 = FALSE;
            (void) sscanf( line,"%*s %lf %lf %lf", &(sInit.T.x), &(sInit.T.y), &(sInit.T.z));
            ligand.S.T.x = sInit.T.x;
            ligand.S.T.y = sInit.T.y;
            ligand.S.T.z = sInit.T.z;
        }
        if (outlev >= 0) {
            pr( logFile, "Initial translation =\t\t\t(%.3f, %.3f, %.3f) Angstroms\n", sInit.T.x, sInit.T.y, sInit.T.z );
        }
        (void) fflush(logFile);
        break;

//______________________________________________________________________________

    case DPF_QUAT0:
    case DPF_AXISANGLE0:
    case DPF_QUATERNION0:
        /*
         * Handles both axisangle0 and quaternion0
         *
         *  axisangle0 1. 0. 0. 0.
         *  axisangle0 random
         *  ( quat0 <--- deprecated )
         *  Initial_quaternion, specified as an axis and angle
         *
         *  quaternion0 0. 0. 0. 1.
         *  quaternion0 random
         *  Initial_quaternion, specified as the four components (qx, qy, qz, qw)
         */
        {
        // Local Block...
        double a, b, c, d;
        (void) sscanf( line, "%*s %s", param[0]);
        for (i=0; i<6; i++) {
            param[0][i] = (char)tolower( (int)param[0][i] );
        }
        if (equal(param[0],"random",6)) {
            // Make a random initial quaternion,
            // and set the boolean B_RandomQuat0 to true,
            // so we can generate random quaternions in population-based methods.
            B_RandomQuat0 = TRUE;
            create_random_orientation( &(sInit.Q) );
            if (outlev >= 0) {
                pr( logFile, "Each run will begin with a new, random initial orientation.\n");
            }
        } else {
            // Read in the user-defined axis-angle values for the initial quaternion
            // and set the boolean B_RandomQuat0 to false,
            B_RandomQuat0 = FALSE;
            (void) sscanf( line, "%*s %lf %lf %lf %lf", &a, &b, &c, &d);
            sInit.Q = (dpf_keyword == DPF_QUATERNION0) ?
                      quatComponentsToQuat(a,b,c,d) :
                      axisDegreeToQuat(a,b,c,d);
        }
        ligand.S.Q = sInit.Q;
        if (outlev >= 0) {
            if (dpf_keyword == DPF_QUATERNION0) {
                pr( logFile, "Initial quaternion,  (x,y,z,w) =\t( %.3f, %.3f, %.3f, %.3f ),\n", sInit.Q.x, sInit.Q.y, sInit.Q.z, sInit.Q.w);
            } else {
                //if (dpf_keyword == DPF_QUAT0 && B_RandomQuat0)  {
                    //pr( logFile, "WARNING quat0 command is obsolete. Now use quaternion0 instead\n");
                //}
                if (dpf_keyword == DPF_QUAT0 && !B_RandomQuat0)  {
                    pr( logFile, "WARNING quat0 command is obsolete. Now use axisangle0 instead\n");
                }
                if (!B_RandomQuat0) {
                    pr( logFile, "Initial axis-angle,  (nx,ny,nz,ang) =\t( %.3f, %.3f, %.3f, %.1f deg ),\n", a, b, c, d );
                }
                pr( logFile, "Initial quaternion,  (x,y,z,w) =\t( %.3f, %.3f, %.3f, %.3f ),\n", sInit.Q.x, sInit.Q.y, sInit.Q.z, sInit.Q.w);
            }
#ifdef DEBUG
            pr( logFile, "Initial Quaternion sInit.Q:\n\n");
            printQuat( logFile, sInit.Q );
            pr( logFile, "Initial Quaternion ligand.S.Q:\n\n");
            printQuat( logFile, ligand.S.Q );
#endif
        }
        (void) fflush(logFile);
        }
        break;

//______________________________________________________________________________

    case DPF_NDIHE:
        /*
        **  ndihe
        **  Number of dihedral angles to be specified by "dihe0"
        */
        (void) sscanf( line, "%*s %d", &ndihed );
        if ( nmol == 0 ) {
            if (outlev >= 0) {
                pr( logFile, "Must specify a ligand PDBQT file, using the \"move\" command.\n");
            }
        } else {
            if (outlev >= 0) {
                pr( logFile, "%s: WARNING!  The \"ndihe\" command is no longer supported.  The number of torsions in the PDBQT file(s) is the number that will be used (i.e. %d)\n", programname, ntor);
            }
            if ( ndihed != ntor ) {
                pr( logFile, "%s: WARNING!  You requested %d torsions, but I found %d in PDBQT-file specifications.\n", programname, ndihed, ntor );
            } /* if */
        }
        (void) fflush(logFile);
        break;

//______________________________________________________________________________

    case DPF_DIHE0:
        /*
        **  dihe0
        **  Initial dihedral angles, input in degrees,
        */
        (void) sscanf( line, "%*s %s", param[0]);
        for (i=0; i<6; i++) {
            param[0][i] = (char)tolower( (int)param[0][i] );
        }
        if (equal(param[0],"random",6)) {
            B_RandomDihe0 = TRUE;
            sInit.ntor = nval = ntor;
            for ( i=0; i<nval; i++ ) {
                sInit.tor[i] = random_range( -180.0, 180.0 );
            }
        } else {
            B_RandomDihe0 = FALSE;
            retval = (int)sscanf( line, torfmt, TOR_ARG_LIST );
            if (retval == 0) {
                pr( logFile, "WARNING!  AutoDock could not read any torsions!\n" );
            } else if (retval == EOF) {
                pr( logFile, "WARNING!  End of file encountered while reading dihe0 line\n");
            } else if (retval < ntor) {
                pr( logFile, "WARNING!  Only %d initial torsion angles were detected on input line.\n",retval);
                pr( logFile, "WARNING!  I am sorry, the number of torsions detected in the PDBQT files was %d torsions.\n", ntor);
            } else {
                if (outlev >= 0) {
                    pr( logFile, "%d initial torsion angles were detected on input line.\n", retval );
                }
            }
            nval = retval;
        }
        if (nval != ntor) {
            pr( logFile, "%s: WARNING!  The number of torsions specified (%d) does not match the number found in the PDBQT file (i.e. %d)\n", programname, nval, ntor);
        }
        for ( i=0; i<nval; i++ ) {
            if (outlev >= 0) {
                pr( logFile, "\tInitial torsion %2d = %7.2f deg\n", (i+1), sInit.tor[i] ); /* sInit.tor is in degrees */
                /* Convert sInit.tor[i] into radians */
            }
            ligand.S.tor[i] = sInit.tor[i] = DegreesToRadians( sInit.tor[i] ); /* sInit.tor is now in radians  Added:05-01-95 */
        }
        (void) fflush(logFile);
        break;

//______________________________________________________________________________

    case DPF_TSTEP:
        /*
        **  tstep
        **  Translation_step,
        */
        retval = (int)sscanf( line, "%*s " FDFMT2, &trnStep0, &trnStepFinal );
        if (retval == 0) {
            pr( logFile, "WARNING!  AutoDock could not read any arguments!\n" );
        } else if (retval == EOF) {
            pr( logFile, "WARNING!  End of file encountered!\n");
        } else if (retval > 0) {
            pr( logFile, "Initial cycle, maximum translation step = +/- %-.1f Angstroms\n", trnStep0);
        }
        if (retval == 2) {
            B_CalcTrnRF = TRUE;
            if (outlev >= 0) {
                pr( logFile, "Final cycle,   maximum translation step = +/- %-.1f Angstroms\n", trnStepFinal);
                pr( logFile, "Reduction factor will be calculated when number of cycles has been read in.\n");
            }
        }
        (void) fflush(logFile);
        break;

//______________________________________________________________________________

    case DPF_QSTEP:
        /*
        **  qstep
        **  Quaternion_step,
        */
        retval = (int)sscanf( line, "%*s " FDFMT2, &qtwStep0, &qtwStepFinal );
        if (retval == 0) {
            pr( logFile, "WARNING!  AutoDock could not read any arguments!\n" );
        } else if (retval == EOF) {
            pr( logFile, "WARNING!  End of file encountered!\n");
        } else if (retval > 0) {
            if (outlev >= 0) {
                pr( logFile, "Initial cycle, maximum quaternion angle step = +/- %-.1f deg\n", qtwStep0);
            }
            /* convert to radians */
            qtwStep0 = DegreesToRadians( qtwStep0 );
        }
        if (retval == 2) {
            B_CalcQtwRF = TRUE;
            if (outlev >= 0) {
                pr( logFile, "Final cycle,   maximum quaternion angle step = +/- %-.1f deg\n", qtwStepFinal);
                pr( logFile, "Reduction factor will be calculated when number of cycles has been read in.\n");
            }
            /* convert to radians */
            qtwStepFinal = DegreesToRadians( qtwStepFinal );
        }
        (void) fflush(logFile);
        break;

//______________________________________________________________________________

    case DPF_DSTEP:
        /*
        **  dstep
        **  Torsion_step,
        */
        retval = (int)sscanf( line, "%*s " FDFMT2, &torStep0, &torStepFinal );
        if (retval == 0) {
            pr( logFile, "WARNING!  AutoDock could not read any arguments!\n" );
        } else if (retval == EOF) {
            pr( logFile, "WARNING!  End of file encountered!\n");
        } else if (retval > 0) {
            if (outlev >= 0) {
                pr( logFile, "Initial cycle, maximum torsion angle step = +/- %-.1f deg\n", torStep0);
            }
            /* convert to radians */
            torStep0 = DegreesToRadians( torStep0 );
        }
        if (retval == 2) {
            B_CalcTorRF = TRUE;
            if (outlev >= 0) {
                pr( logFile, "Final cycle,   maximum torsion angle step = +/- %-.1f deg\n", torStepFinal);
                pr( logFile, "Reduction factor will be calculated when number of cycles has been read in.\n");
            }
            /* convert to radians */
            torStepFinal = DegreesToRadians( torStepFinal );
        }
        (void) fflush(logFile);
        break;

//______________________________________________________________________________

    case DPF_TRNRF:
        /*
        **  trnrf
        **  Translation reduction factor,
        */
        (void) sscanf( line, "%*s " FDFMT, &trnFac );
        if (outlev >= 0) {
            pr( logFile, "Reduction factor for translations =\t%-.3f /cycle\n", trnFac );
        }
        B_trnReduc = (trnFac != 1.);
        (void) fflush(logFile);
        break;

//______________________________________________________________________________

    case DPF_QUARF:
        /*
        **  quarf
        **  Quaternion reduction factor,
        */
        (void) sscanf( line, "%*s " FDFMT, &qtwFac );
        if (outlev >= 0) {
            pr( logFile, "Reduction factor for quaternion angle =\t%-.3f /cycle\n", qtwFac );
        }
        B_qtwReduc = (qtwFac != 1.);
        (void) fflush(logFile);
        break;

//______________________________________________________________________________

    case DPF_DIHRF:
        /*
        **  dihrf
        **  Torsion reduction factor,
        */
        (void) sscanf( line, "%*s " FDFMT, &torFac );
        if (outlev >= 0) {
            pr( logFile, "Reduction factor for torsion angles =\t%-.3f /cycle\n", torFac );
        }
        B_torReduc = (torFac != 1.);
        (void) fflush(logFile);
        break;

//______________________________________________________________________________

    case DPF_FLEX:
        /*
        **  flex
        **  Flexible side-chains, cannot translate:
        */
        nmol++;
        nres++;
        (void) fflush(logFile);
        break;

//______________________________________________________________________________

    case DPF_INTNBP_REQM_EPS:
    case DPF_INTNBP_COEFFS:
        /*
        **  intnbp_r_eps
        **  Read internal energy parameters:
        **  Lennard-Jones and Hydrogen Bond Potentials,
        **  DPF_INTNBP_REQM_EPS: Using epsilon and r-equilibrium values...
        **  DPF_INTNBP_COEFFS: Using coefficients...
        */
        (void) sscanf( line, "%*s " FDFMT2 " %d %d %s %s", &Rij, &epsij, &xA, &xB, param[0], param[1] );
        if ( dpf_keyword == DPF_INTNBP_REQM_EPS ) {
        /* check that the Rij is reasonable */
	/* SF ...but only if there are no G-atoms. */        /* SF RING CLOSURE */

	if ((Rij <= 2.0 ) && (epsij >= EPSIJ_MAX )) {    /* RING CLOSURE */
 	     (void) fprintf( logFile, "Ring closure distance potential found for atom type %s :\n    Equilibrium distance   = %.2f Angstroms \n    Equilibrium potential  = %.6f Kcal/mol\n    Pseudo-LJ coefficients = %d-%d \n\n", param[1] , Rij, epsij, xA, xB); /* SF RING CLOSURE */
			}   /* SF RING CLOSURE */
	else { /* SF RING CLOSURE */

	        if ((Rij < RIJ_MIN) || (Rij > RIJ_MAX)) {
        	    (void) fprintf( logFile,
	            "WARNING: pairwise distance, Rij, %.2f, is not a very reasonable value for the equilibrium separation of two atoms! (%.2f Angstroms <= Rij <= %.2f Angstroms)\n\n", Rij, RIJ_MIN, RIJ_MAX);
	            (void) fprintf( logFile, "Perhaps you meant to use \"intnbp_coeffs\" instead of \"intnbp_r_eps\"?\n\n");
	            /* GMM COMMENTED OUT FOR DAVE GOODSELL, MUTABLE ATOMS
	             * exit(-1); */
	     	     }
	        /* check that the epsij is reasonable */
	        if ((epsij < EPSIJ_MIN) || (epsij > EPSIJ_MAX)) {
	            (void) fprintf( logFile,
	            "WARNING: well-depth, epsilon_ij, %.2f, is not a very reasonable value for the equilibrium potential energy of two atoms! (%.2f kcal/mol <= epsilon_ij <= %.2f kcal/mol)\n\n", epsij, EPSIJ_MIN, EPSIJ_MAX);
	            (void) fprintf( logFile, "Perhaps you meant to use \"intnbp_coeffs\" instead of \"intnbp_r_eps\"? \n\n");
	            /* GMM COMMENTED OUT FOR DAVE GOODSELL, MUTABLE ATOMS
	             * exit(-1); */
	        }
	        } /* RING CLOSURE */

	     }

        /* Defend against division by zero... */
        if (xA != xB) {
            if ( dpf_keyword == DPF_INTNBP_REQM_EPS ) {
            // Calculate the coefficients from Rij and epsij
            cA = (tmpconst = epsij / (Real)(xA - xB)) * pow( (double)Rij, (double)xA ) * (Real)xB;
            cB = tmpconst * pow( (double)Rij, (double)xB ) * (Real)xA;
            } else {
            cA = Rij;
            cB = epsij;
            }

            int a[2]; /* atom types of this interaction pair */
            for (int i=0;i<2;i++) {
                foundParameter = apm_find(param[i]);
                if ( NULL == foundParameter ) {
                    prStr( error_message,"%s: ERROR:  Unknown ligand atom type \"%s\"; add parameters for it to the parameter library first!\n", programname, param[i]);
                    exit(-1);
                }
                else a[i] = foundParameter->map_index;
            }
            pr(logFile, "\nCalculating internal non-bonded interaction energies for docking calculation;\n");
            intnbtable( &B_havenbp, a[0], a[1], info, cA, cB, xA, xB, AD4.coeff_desolv, sigma, ad_energy_tables, BOUND_CALCULATION );
           // pr(logFile, "\nCalculating internal non-bonded interaction energies for unbound conformation calculation;\n");
            //intnbtable( &B_havenbp, a[0], a[1], info, cA_unbound, cB_unbound, xA_unbound, xB_unbound, AD4.coeff_desolv, sigma, unbound_energy_tables, UNBOUND_CALCULATION );
        } else {
            pr(logFile,"WARNING: Exponents must be different, to avoid division by zero!\n\tAborting...\n");
            exit(-1);
        }
        (void) fflush(logFile);
        break;

//______________________________________________________________________________


    case DPF_UNBOUND_INTNBP_COEFFS:
        /*
        **  unbound_intnbp_coeffs
        **  Read internal energy parameters for unbound extended state calculation:
        */
        (void) sscanf( line, "%*s " FDFMT2 " %d %d", &cA_unbound, &cB_unbound, &xA_unbound, &xB_unbound );

        pr(logFile, "\nSetting the internal non-bonded interaction energy parameters for the\nunbound docking calculation, E = %.1f / r^%d - %.1f / r^%d\n\n", cA_unbound, xA_unbound, cB_unbound, xB_unbound);
        (void) fflush(logFile);
        break;

//______________________________________________________________________________

    case DPF_RT0:
        /*
        **  rt0
        **  Initial Temperature,
        */
        (void) sscanf( line, "%*s " FDFMT, &RT0 );
        if (RT0 <= 0.) {
            pr( logFile, "\nWARNING!  Negative temperatures not allowed! Will default to RT = 616 cal mol.\n" );
            RT0 = 616.0;
        }
        if (outlev >= 0) {
            pr( logFile, "\n\t\tTEMPERATURE SCHEDULE INFORMATION\n" );
            pr( logFile, "\t\t________________________________\n\n" );
            pr( logFile, "               -1 -1                 -1 -1\n" );
            pr( logFile, "R = %5.3f J mol  K    = %5.3f cal mol  K  \n\n", RJ, Rcal );
            pr( logFile, "                                        -1\n" );
            pr( logFile, "Initial R*Temperature = %8.2f cal mol\n", RT0 );
            pr( logFile, "      (=> Temperature = %8.2f K\n", RT0/Rcal );
            pr( logFile, "                   or = %8.2f C)\n\n", RT0/Rcal - T0K );
        }
        break;

//______________________________________________________________________________

    case DPF_RTRF:
        /*
        **  rtrf
        **  Temperature reduction factor,
        */
        (void) sscanf( line, "%*s " FDFMT, &RTFac);
        if (outlev >= 0) {
            pr( logFile, "R*Temperature reduction factor = %8.2f\t/cycle\n", RTFac );
        }
        if (RTFac >= 1.) {
            stop("Cooling is impossible with a reduction\n\tfactor greater than or equal to 1.0!" );
            exit( -1 );
        } else if (RTFac == 0.0 ) {
            stop("Cooling is impossible with a ZERO reduction factor!" );
            exit( -1 );
        } else if (RTFac < 0.0 ) {
            stop("Cooling is impossible with a NEGATIVE reduction factor!" );
            exit( -1 );
        }
        (void) fflush(logFile);
        B_tempChange = ( RTFac != 1.0 );
        break;

//______________________________________________________________________________

    case DPF_RUNS:
        /*
        **  runs
        **  Number of docking runs,
        */
        (void) sscanf( line, "%*s %d", &nruns );
        if ( nruns > MAX_RUNS ) {
            prStr( error_message, "%s:  ERROR: %d runs were requested, but AutoDock is only dimensioned for %d.\nChange \"MAX_RUNS\" in \"constants.h\".", programname, nruns, MAX_RUNS);
            stop( error_message );
            exit( -1 );
        }
        pr( logFile, "Number of runs =\t\t\t\t%8d run%c\n", nruns, (nruns > 1)?'s':' ');
        (void) fflush(logFile);
        break;

//______________________________________________________________________________

    case DPF_CYCLES:
        /*
        **  cycles
        **  Number of constant temperature SA cycles,
        */
        (void) sscanf( line, "%*s %d", &ncycles );
        if (ncycles < 0) {
            pr( logFile, "WARNING!  Negative number of cycles found!  Using default value.\n");
            ncycles = 50;
        }
        pr( logFile, "Maximum number of cycles =\t\t\t%8d cycles\n\n", ncycles);
        if (B_linear_schedule) {
            RTreduc = RT0 / ncycles;
            if (outlev >= 0) {
                pr( logFile, "\nA linear temperature reduction schedule was requested...\n" );
                pr( logFile, "Annealing temperature will be reduced by %.3f cal mol per cycle.\n\n", RTreduc );
            }
        }
        (void) fflush(logFile);
        break;

//______________________________________________________________________________

    case DPF_ACCS:
        /*
        **  accs
        **  Maximum number of steps accepted,
        */
        (void) sscanf( line, "%*s %d", &naccmax );
        if (naccmax < 0) {
            naccmax = 100;
            pr( logFile, "WARNING!  Negative number of accepted moves found!  Using default value.\n");
        }
        if (outlev >= 0) {
            pr( logFile, "Maximum number accepted per cycle =\t\t%8d steps\n", naccmax);
        }
        if (nrejmax != 0) {
            nstepmax = naccmax + nrejmax;
            if (outlev >= 0) {
                pr( logFile, "                                           \t_________\n" );
                pr( logFile, "Maximum possible number of steps per cycle =\t%8d\tsteps\n\n", nstepmax);
            }
        }
        (void) fflush(logFile);
        break;

//______________________________________________________________________________

    case DPF_REJS:
        /*
        **  rejs
        **  Maximum number of steps rejected,
        */
        (void) sscanf( line, "%*s %d", &nrejmax );
        if (nrejmax < 0) {
            nrejmax = 100;
            pr( logFile, "WARNING!  Negative number of rejected moves found!  Using default value.\n");
        }
        if (outlev >= 0) {
            pr( logFile, "Maximum number rejected per cycle =\t\t%8d steps\n", nrejmax);
        }
        if (naccmax != 0) {
            nstepmax = naccmax + nrejmax;
            if (outlev >= 0) {
                pr( logFile, "                                           \t_________\n" );
                pr( logFile, "Maximum possible number of steps per cycle =\t%8d steps\n\n", nstepmax);
            }
        }
        (void) fflush(logFile);
        break;

//______________________________________________________________________________

    case DPF_SELECT:
        /*
        **  select
        **  Select either minimum or last state from previous cycle,
        */
        (void) sscanf( line, "%*s %c", &selminpar );
        B_selectmin = (selminpar == 'm');
        if ( B_selectmin ) {
            if (outlev >= 0) {
                pr( logFile, "%s will begin each new cycle\nwith the state of minimum energy from the previous annealing cycle.\n", programname);
            }
        } else {
            if (outlev >= 0) {
                pr( logFile, "%s will begin each new cycle\nwith the last state from the previous annealing cycle.\n", programname);
            }
        }
        (void) fflush(logFile);
        break;

//______________________________________________________________________________

    case DPF_RMSTOL:
        /*
        **  rmstol
        **  Cluster tolerance,
        */
        (void) sscanf( line, "%*s " FDFMT, &clus_rms_tol);
        if (outlev >= 0) {
            pr( logFile, "Maximum RMS tolerance for conformational cluster analysis = %.2f Angstroms\n", clus_rms_tol);
        }
        (void) fflush(logFile);
        break;

//______________________________________________________________________________

    case DPF_RMSREF:
        /*
        **  rmsref
        **  RMS Reference Coordinates:
        */
        (void) sscanf( line, "%*s %s", FN_rms_ref_crds);
        if (outlev >= 0) {
            pr( logFile, "RMS reference coordinates will taken from \"%s\"\n", FN_rms_ref_crds );
        }
        (void) fflush(logFile);
        break;

//______________________________________________________________________________

    case DPF_RMSATOMS:
        /*
        **  rmsatoms ligand_only
        **  rmsatoms all
        **
        **  Set the atoms to compute the RMSD values for cluster analysis
        **  either "ligand_only" (the default) or "all" moving atoms (ligand + receptor)
        */
        retval = sscanf( line, "%*s %s", rms_atoms_cmd);
        if (retval != 1) {
            pr( logFile, "%s:  ERROR: please specify an argument (either \"ligand_only\" or \"all\").  By default, only the ligand atoms will be used for the cluster analysis.\n", programname );
            B_rms_atoms_ligand_only = TRUE;  // cluster on the ligand atoms only
        } else {
            if ( strncmp( rms_atoms_cmd, "ligand_only", 11 ) == 0 ) {
                if (outlev >= 0) {
                    pr( logFile, "RMS clustering will be performed on the ligand atoms only.\n" );
                }
                B_rms_atoms_ligand_only = TRUE;  // cluster on the ligand atoms only
            } else if ( strncmp( rms_atoms_cmd, "all", 3 ) == 0 ) {
                if (outlev >= 0) {
                    pr( logFile, "RMS clustering will be performed on the moving atoms of the receptor plus all the ligand atoms.\n" );
                }
                B_rms_atoms_ligand_only = FALSE;  // cluster on the ligand atoms plus moving receptor atoms
            } else {
                if (outlev >= 0) {
                    pr( logFile, "RMS clustering will be performed on the ligand atoms only.\n" );
                }
                B_rms_atoms_ligand_only = TRUE;  // cluster on the ligand atoms only
            }
        }
        (void) fflush(logFile);
        break;

//______________________________________________________________________________

    case DPF_TRJFRQ:
        /*
        **  trjfrq
        **  Trajectory frequency,
        */
        (void) sscanf( line, "%*s %d", &trj_freq);
        B_write_trj = (trj_freq > 0);
        if (outlev >= 0) {
            pr( logFile, UnderLine );
            pr( logFile, "\t\tTRAJECTORY INFORMATION\n" );
            pr( logFile, "\t\t______________________\n\n\n" );
        }
        if (B_write_trj) {
            if (outlev >= 0) {
                pr( logFile, "Output frequency for trajectory frames =\tevery %d step%s\n", trj_freq, (trj_freq > 1)?"s.":"." );
            }
        } else {
            if (outlev >= 0) {
                pr( logFile, "No trajectory of states will be written.\n\n" );
                pr( logFile, "Subsequent \"trjbeg\", \"trjend\", \"trjout\" and \"trjsel\" parameters will be ignored.\n\n" );
            }
        }
        (void) fflush(logFile);
        break;

//______________________________________________________________________________

    case DPF_TRJBEG:
        /*
        **  trjbeg
        **  Trajectory begin cycle,
        */
        (void) sscanf( line, "%*s %d", &trj_begin_cyc );
        if (outlev >= 0) {
            pr( logFile, "Begin outputting trajectory of states at cycle:\t%d\n", trj_begin_cyc );
        }
        if (trj_begin_cyc < 0) {
            trj_begin_cyc = 0;
        } else if (trj_begin_cyc > ncycles) {
            trj_begin_cyc = trj_end_cyc = ncycles;
        }
        --trj_begin_cyc;
        (void) fflush(logFile);
        break;

//______________________________________________________________________________

    case DPF_TRJEND:
        /*
        **  trjend
        **  Trajectory end cycle,
        */
        (void) sscanf( line, "%*s %d", &trj_end_cyc );
        if (outlev >= 0) {
            pr( logFile, "Cease outputting trajectory of states at cycle:\t%d\n", trj_end_cyc );
        }
        if (trj_end_cyc > ncycles) {
            trj_end_cyc = ncycles;
        } else if (trj_end_cyc < 0) {
            trj_end_cyc = 1;
        }
        --trj_end_cyc;
        (void) fflush(logFile);
        break;

//______________________________________________________________________________

    case DPF_TRJOUT:
        /*
        **  trjout
        **  Trajectory file,
        */
        (void) sscanf( line, "%*s %s", FN_trj );
        if (outlev >= 0) {
            pr( logFile, "\nWrite trajectory of state variables to file: \"%s\"\n", FN_trj);
        }
        (void) fflush(logFile);
        break;

//______________________________________________________________________________

    case DPF_TRJSEL:
        /*
        **  trjsel
        **  Trajectory select,
        */
        (void) sscanf( line, "%*s %c", &out_acc_rej );
        B_acconly = (out_acc_rej == 'A');
        B_either  = (out_acc_rej == 'E');
        if (B_acconly) {
            if (outlev >= 0) {
                pr( logFile, "Output *accepted* states only.\n" );
            }
        } else if (B_either) {
            if (outlev >= 0) {
                pr( logFile, "Output *either* accepted or rejected states.\n" );
            }
        } else {
            pr( logFile, "WARNING: Missing or unknown accepted/rejected output flag.\n" );
        }
        (void) fflush(logFile);
        break;

//______________________________________________________________________________

    case DPF_EXTNRG:
        /*
        **  extnrg
        **  Wall Energy,
        */
        (void) sscanf( line, "%*s " FDFMT, &WallEnergy );
        if (outlev >= 0) {
            pr( logFile, "External grid energy (beyond grid map walls) = %.2f\n\n", WallEnergy );
        }
        (void) fflush(logFile);
        break;

//______________________________________________________________________________

    case DPF_CLUSTER:
        /*
        **  cluster
        **  Cluster mode,
        */
        (void) sscanf( line, "%*s %s", FN_clus );
        B_cluster_mode = TRUE;
        if (!B_found_about_keyword){
                prStr(error_message, "%s:  ERROR:  no \"about\" command has been specified!\n", programname);
                stop(error_message);
                exit(-1);
        }
        if (outlev >= 0) {
            pr( logFile, "Cluster mode is now set.\n\n" );
        }
        clmode( num_atom_types, clus_rms_tol,
                hostnm, jobStart, tms_jobStart,
                B_write_all_clusmem, FN_clus, crdpdb, lig_center,
                B_symmetry_flag, FN_rms_ref_crds );
        (void) fflush(logFile);
        break;

//______________________________________________________________________________

    case DPF_CLUSALL:
        /*
        ** write_all_clusmem
        ** Write all cluster members...
        */
        B_write_all_clusmem = TRUE;
        if (outlev >= 0) {
            pr( logFile, "All members of each cluster will be written out after the clustering histogram.\n(This is instead of outputting just the lowest energy member in each.)\n\n" );
        }
        break;

//______________________________________________________________________________

    case DPF_RMSNOSYM:
        /*
        **  rmsnosym
        **  Calculate RMS values in the normal way,
        **  ignoring any atom-type equivalences...
        */
        B_symmetry_flag = FALSE;
        if (outlev >= 0) {
            pr( logFile, "Symmetry will be ignored in RMS calculations.\n\n" );
        }
        (void) fflush(logFile);
        break;

//______________________________________________________________________________

    case DPF_SCHEDLIN:
        /*
        **  linear_schedule
        **  Use a linear (arithmetic) temperature
        **  reduction schedule.  This is necessary for
        **  more accurate entropy estimations...
        */
        B_linear_schedule = TRUE;
        if (outlev >= 0) {
            pr( logFile, "A linear temperature reduction schedule will be used...\n\n" );
        }
        if (ncycles == -1) {
            pr( logFile, "\nWARNING!  Please specify the number of cycles first!\n\n" );
        } else {
            RTreduc = RT0 / ncycles;
            if (outlev >= 0) {
                pr( logFile, "Annealing temperature will be reduced by %.3f cal mol per cycle.\n\n", RTreduc );
            }
        }
        (void) fflush(logFile);
        break;

//______________________________________________________________________________

    case DPF_WATCH:
        /*
        **  watch
        **  for watching a job's progress PDBQT file in AVS,
        */
        (void) sscanf( line, "%*s %s", FN_watch);
        if (B_write_trj) {
            pr(logFile,"\nAutoDock will create the watch-file \"%s\", for real-time monitoring of runs.\n\n", FN_watch);
            pr(logFile,"\nThe watch-file will be updated every %d moves, in accordance with the trajectory parameters..\n\n", trj_freq);
            B_watch = TRUE;
        } else {
            pr(logFile,"\nYou must set \"trjfrq\" to be greater than zero. No watch-file will be created.\n\n");
            B_watch = FALSE;
        }
        (void) fflush(logFile);
        break;

//______________________________________________________________________________

    case DPF_GAUSSTORCON:
    case DPF_HARDTORCON:
        /*
        ** "gausstorcon" Add Gaussian torsion contraints,
        ** "hardtorcon"  Add Hard torsion contraints,
        */
        (void) sscanf( line, "%*s %d " FDFMT2, &I_tor, &F_torPref, &F_torHWdth);
        if (I_tor <= 0) {
            pr( logFile, "\nTorsion IDs less than 1 (%d) are not allowed!\n\n", I_tor);
        } else if (I_tor > ntor) {
            pr( logFile, "\nRequested torsion ID (%d) is larger than the number of torsions found (%d)!\n\n", I_tor, ntor);
        } else { /* torsion-ID accepted */
            --I_tor;    /* Because humans start at 1, and C at 0... */
            if ( B_isTorConstrained[I_tor] == 0 ) {

                if (dpf_keyword ==  DPF_GAUSSTORCON) {
                    B_isGaussTorCon = TRUE;
                    B_isTorConstrained[I_tor] = 1;
                    /*
                    ** Initialize... Torsion Energy Profile...
                    ** Set energies at every torsion division
                    ** to the user-defined (maximum) barrier energy,
                    */
                    for (US_tD = 0;  US_tD < NTORDIVS;  US_tD++) {
                        US_torProfile[I_tor][US_tD] = US_torBarrier;
                    }
                } else {
                    /*
                    ** DPF_HARDTORCON
                    */
                    B_isTorConstrained[I_tor] = 2;
                }
            }
            if (dpf_keyword ==  DPF_GAUSSTORCON) {
                (void) strcpy( S_contype, " half-" );
            } else {
                (void) strcpy( S_contype, " " );
            }
                /*
            ** F_torPref ranges from -180.0 to +180.0 degrees...
            */
            F_torPref = WrpDeg(ModDeg(F_torPref));
            if (F_torHWdth < 0.) {
                pr(logFile,"\nI'm sorry, negative%swidths (%.1f) are not allowed. I will use the default (%.1f) instead.\n\n", S_contype, F_torHWdth, DEFHWDTH);
                F_torHWdth = DEFHWDTH;
            } else if (F_torHWdth > 90.) {
                pr(logFile,"\nI'm sorry, your requested%swidth (%.1f) is too large. I will use the default (%.1f) instead.\n\n", S_contype, F_torHWdth, DEFHWDTH);
                F_torHWdth = DEFHWDTH;
            }
            pr( logFile, "For torsion %d, Adding a constrained-torsion zone centered on %.1f degrees;\n%swidth = %.1f degrees.\n\n", 1+I_tor, F_torPref, S_contype, F_torHWdth);

            if (dpf_keyword == DPF_GAUSSTORCON) {
                /*
                ** Calculate the torsion energy profile;
                ** combine this with previous profile without
                ** losing any information.
                */
                for (F_A = F_A_from;  F_A <= F_A_to;  F_A += F_W) {
                    F_Aova = (F_A - F_torPref) / F_torHWdth;
                    US_energy = (unsigned short) (((Real)US_torBarrier) * (1.0 - exp(F_lnH * F_Aova*F_Aova)));
                    /*
                    ** if F_A(<-180.or>180), wrap to -180to180,
                    */
                    F_tor = WrpDeg(ModDeg(F_A));
                    /*
                    ** Convert from F_tor to US_tD
                    */
                    US_tD = (unsigned short) ((F_tor - F_hW + 180.)/F_W);
                    US_torProfile[I_tor][US_tD] = min(US_energy,US_torProfile[I_tor][US_tD]);
                }/* for F_A */
                /*
                ** Ensure that the lowest point(s) in the profile are
                ** zero...
                */
                US_min = TORBARMAX;
                for (US_tD = 0;  US_tD < NTORDIVS;  US_tD++) {
                    US_min = min(US_min,US_torProfile[I_tor][US_tD]);
                }
                for (US_tD = 0;  US_tD < NTORDIVS;  US_tD++) {
                    US_torProfile[I_tor][US_tD] -= US_min;
                }
            } else { /*DPF_HARDTORCON*/

                iCon = N_con[I_tor] + 1;
                if (iCon < MAX_TOR_CON) {
                    F_TorConRange[I_tor][N_con[I_tor]][LOWER] = F_torPref - 0.5* F_torHWdth;
                    F_TorConRange[I_tor][N_con[I_tor]][UPPER] = F_torPref + 0.5* F_torHWdth;
                    N_con[I_tor] = iCon;
                } else {
                    pr(logFile,"\n\n I'm sorry, you can only have %d (=MAX_TOR_CON) torsion constraints.\nIf you need more, change the \"#define MAX_TOR_CON\" line in \"constants.h\".\n\n",MAX_TOR_CON);
                }/* Still room to add another constraint. */
            } /*DPF_HARDTORCON*/
        }/* torsion-ID accepted */
        (void) fflush(logFile);
        break;

//______________________________________________________________________________

    case DPF_BARRIER:
        /*
        **  barrier
        **  Define torsion-barrier energy...
        */
        (void) sscanf( line, "%*s %d", &I_torBarrier);
        US_torBarrier = (unsigned short)I_torBarrier;
        US_torBarrier = min(US_torBarrier, TORBARMAX);
        pr(logFile,"\nTorsion barrier energy is set to %uhd\n\n", US_torBarrier);
        (void) fflush(logFile);
        break;

//______________________________________________________________________________

    case DPF_SHOWTORPEN:
        /*
        **  showtorpen
        **  Show torsion's penalty energy.
        */
        B_ShowTorE = TRUE;
        pr(logFile,"\nConstrained torsion penalty energies will be stored during docking, and output after each run\n\n");
        (void) fflush(logFile);
        break;

//______________________________________________________________________________

    case DPF_E0MAX:
        /*
        **  e0max
        **  Set maximum initial energy,
        */
        retval = sscanf( line, "%*s " FDFMT " %d", &e0max, &MaxRetries );
        if (retval == 0) {
            pr( logFile, "Could not read any arguments!\n" );
        } else if (retval == EOF) {
            pr( logFile, "End of file encountered!\n");
        } else if (retval == 1) {
            pr(logFile,"Using the default maximum number of retries for initialization, %d retries.\n\n", MaxRetries);
        } else if (retval == 2) {
            pr(logFile,"Using user-specified maximum number of retries for initialization, %d retries.\n\n", MaxRetries);
        }
        if (e0max < 0.) {
            e0max = 1000.0;
        }
        pr(logFile,"\nIf the initial energy is greater than e0max, %.3f,\nthen a new, random initial state will be created.\n\n",e0max);
        (void) fflush(logFile);
        break;

//______________________________________________________________________________

    case DPF_SIMANNEAL:
        /*
        ** simanneal
        */
        /*
        ** Calculate reduction factor based on initial and final step values,
        ** and number of cycles...
        */
            ncycm1 = ncycles - 1;
            if (ncycm1 < 0) {
                ncycm1 = 1;
            }
            if (B_CalcTrnRF) {
                trnFac = RedFac(trnStep0, trnStepFinal, ncycm1);
                pr( logFile, "Calculated reduction factor for translations     = %-.3f /cycle\n", trnFac);
                B_trnReduc = (trnFac != 1.);
            }
            if (B_CalcQtwRF) {
                qtwFac = RedFac(qtwStep0, qtwStepFinal, ncycm1);
                pr( logFile, "Calculated reduction factor for quaternion angle = %-.3f /cycle\n", qtwFac );
                B_qtwReduc = (qtwFac != 1.);
            }
            if (B_CalcTorRF) {
                torFac    = RedFac(torStep0, torStepFinal, ncycm1);
                pr( logFile, "Calculated reduction factor for torsion angles   = %-.3f /cycle\n", torFac );
                B_torReduc = (torFac != 1.);
            }
            pr(logFile, "\n");
            /*
            ** Number of ligands read in...
            */
            if (nlig > 0) {
                pr( logFile, "\nTotal number of ligands read in by the DPF \"move\" command = %d\n\n", nlig );
            }
            if (nres > 0) {
                pr( logFile, "\nTotal number of residues read in by the DPF \"flex\" command = %d\n\n", nres );
            }
            pr(logFile, "\n");

            if (!B_found_about_keyword){
                    prStr(error_message, "%s:  ERROR:  no \"about\" command has been specified!\n", programname);
                    stop(error_message);
                    exit(-1);
            }

            if (B_havenbp) {
                nbe( info, ad_energy_tables, num_atom_types );
            }
            if (B_cluster_mode) {
                clmode( num_atom_types, clus_rms_tol,
                        hostnm, jobStart, tms_jobStart,
                        B_write_all_clusmem, FN_clus, crdpdb, lig_center,
                        B_symmetry_flag, FN_rms_ref_crds );
            }
            for (j = 0; j < MAX_RUNS; j++) {
                econf[j] = torsFreeEnergy;
            }
            if (ad4_unbound_model==Unbound_Default) ad4_unbound_model = Unbound_Same_As_Bound;
            pr(logFile, "Unbound model to be used is %s.\n", report_parameter_library());
            /* ___________________________________________________________________
            **
            ** Begin the automated docking simulation,
            ** ___________________________________________________________________
            */
            simanneal( &nconf, Nnb, WallEnergy, atomstuff, charge, abs_charge, qsp_abs_charge, B_calcIntElec,
                        crd, crdpdb, dock_param_fn,
                        ad_energy_tables,
                        econf, B_either,
                        elec, emap,
                        ncycles, nruns, jobStart,
                        map,
                        naccmax, natom, nonbondlist, nrejmax, ntor1, ntor, outlev,
                        sInit, sHist,   qtwFac, B_qtwReduc, qtwStep0,
                        B_selectmin, FN_ligand, lig_center, RT0, B_tempChange, RTFac,
                        tms_jobStart, tlist, torFac, B_torReduc, torStep0,
                        FN_trj, trj_end_cyc, trj_begin_cyc, trj_freq, trnFac,
                        B_trnReduc, trnStep0, type, vt, B_write_trj,
                        B_constrain_dist, atomC1, atomC2, sqlower, squpper,
                        B_linear_schedule, RTreduc,
                        /*maxrad,*/
                        B_watch, FN_watch,
                        B_isGaussTorCon, US_torProfile, B_isTorConstrained,
                        B_ShowTorE, US_TorE, F_TorConRange, N_con,
                        B_RandomTran0, B_RandomQuat0, B_RandomDihe0,
                        e0max, torsFreeEnergy, MaxRetries, ligand_is_inhibitor,
                        ignore_inter,
                        B_include_1_4_interactions, scale_1_4,
                        parameterArray, unbound_internal_FE,
                        info, B_use_non_bond_cutoff,
                        B_have_flexible_residues,
                        PDBQT_record,
                        ad4_unbound_model
                        );

            (void) fflush(logFile);
        break;

//______________________________________________________________________________

    case DPF_SET_GA:

      if (GlobalSearchMethod != NULL) {
          pr(logFile, "Deleting the previous settings for the Genetic Algorithm.\n");
          (void) fflush(logFile);
          delete GlobalSearchMethod;
          GlobalSearchMethod = NULL;
      }

      if (debug > 0) {
        if(outputEveryNgens>0) pr( logFile, "\n\tOutput population statistics every %u generations.\n", outputEveryNgens );
        else pr( logFile, "\n\tNever output generation-based population statistics.\n");
      }
      GlobalSearchMethod = new Genetic_Algorithm(e_mode, s_mode, c_mode, w_mode, elitism, c_rate, m_rate,
                                                 window_size, num_generations, outputEveryNgens );
      ((Genetic_Algorithm *)GlobalSearchMethod)->mutation_values( low, high, alpha, beta, trnStep0, qtwStep0, torStep0  );
      ((Genetic_Algorithm *)GlobalSearchMethod)->initialize(pop_size, 7+sInit.ntor);

      if (s_mode==LinearRanking){
        (void)((Genetic_Algorithm *)GlobalSearchMethod)->set_linear_ranking_selection_probability_ratio(linear_ranking_selection_probability_ratio);
          pr( logFile, "\n\tSet linear_ranking_selection_probability_ratio to %f.\n", linear_ranking_selection_probability_ratio );
      }
      

      (void) fflush(logFile);
      break;
//______________________________________________________________________________

    case DPF_SET_SW1:

      if (LocalSearchMethod != NULL) {
          pr(logFile, "Deleting the previous settings for the local search Solis-Wets algorithm (SW1 object).\n");
          delete LocalSearchMethod;
          LocalSearchMethod = NULL;
      }

      pr(logFile, "Creating a new Local Search object using the Solis-Wets algorithm (SW1) with the current settings.\n\n");
      LocalSearchMethod = new Solis_Wets1(7+sInit.ntor, max_its, max_succ, max_fail, rho, lb_rho, 2.0, 0.5, search_freq);

      (void) fflush(logFile);
      break;
//______________________________________________________________________________

    case DPF_SET_PSW1:

      if (LocalSearchMethod != NULL) {
          pr(logFile, "Deleting the previous settings for the local search pseudo-Solis-Wets algorithm (pSW1 object).\n");
          delete LocalSearchMethod;
          LocalSearchMethod = NULL;
      }

      pr(logFile, "Creating a new Local Search object using the pseudo-Solis-Wets algorithm (pSW1) with the current settings.\n\n");

      //  Allocate space for the variable rho's
      rho_ptr = new Real[7+sInit.ntor];
      lb_rho_ptr = new Real[7+sInit.ntor];

      //  Initialize the rho's corresponding to the translation
      //  0,1,2   x,y,z
      //  3,4,5,6 qx,qy,qz,qw
      //  7,...   tor1
//These scale values can be changed in the dpf, officially unsupported 4/2009
//Real psw_trans_scale = 1.0;//1 angstrom
//Real psw_rot_scale = 0.05;  //about 3 degrees, we think
//Real psw_tors_scale = 0.1; //about 6 degrees

      for (j=0; j<3; j++) {
         // j=0,1,2
         rho_ptr[j] = rho * psw_trans_scale;// formerly trnStep0;
         lb_rho_ptr[j] = lb_rho * psw_trans_scale; //once trnStepFinal;
      }

      //  Initialize the rho's corresponding to the quaterion
      for (; j<7; j++) {
         // j=3,4,5,6
         rho_ptr[j] = rho * psw_rot_scale;// formerly qtwStep0;
         lb_rho_ptr[j] = lb_rho * psw_rot_scale; //once qtwStepFinal;
      }

      //  Initialize the rho's corresponding to the torsions
      for (; j<7+sInit.ntor; j++) {
         // j=7,...
         rho_ptr[j] = rho * psw_tors_scale;// formerly torStep0;
         lb_rho_ptr[j] = lb_rho * psw_tors_scale;//formerly torStepFinal;
      }

      LocalSearchMethod = new Pseudo_Solis_Wets1(7+sInit.ntor, max_its, max_succ, max_fail, 2.0, 0.5, search_freq, rho_ptr, lb_rho_ptr);

      (void) fflush(logFile);
      break;

//______________________________________________________________________________

    case DPF_SET_PATTERN:

      if (LocalSearchMethod != NULL) {
          pr(logFile, "Deleting the previous settings for the local search Pattern Search algorithm (PS object).\n");
          delete LocalSearchMethod;
          LocalSearchMethod = NULL;
      }


      pr(logFile, "Creating a new Local Search object using the Pattern Search algorithm (PS) with the current settings.\n\n");
      LocalSearchMethod = new Pattern_Search(7+sInit.ntor, max_succ, rho, lb_rho, 2.0, 0.5, search_freq);

      (void) fflush(logFile);
      break;

//______________________________________________________________________________

    case DPF_GALS:
        (void) fflush( logFile );
        /*
        ** Genetic Algorithm-Local search,  a.k.a. Lamarckian Genetic Algorithm
        */
            (void) sscanf( line, "%*s %d",&nruns );
            if ( nruns > MAX_RUNS ) {
                prStr( error_message, "%s:  ERROR: %d runs requested, but only dimensioned for %d.\nChange \"MAX_RUNS\" in \"constants.h\".", programname, nruns, MAX_RUNS);
                stop( error_message );
                exit( -1 );
            } else if ((GlobalSearchMethod==NULL)||(LocalSearchMethod==NULL)) {
                prStr(error_message, "%s:  ERROR:  You must use \"set_ga\" to allocate both Global Optimization object AND Local Optimization object.\n", programname);
                stop(error_message);
                exit(-1);
            }
            exit_if_missing_elecmap_desolvmap_about("gals");
            pr( logFile, "Number of requested LGA dockings = %d run%c\n", nruns, (nruns > 1)?'s':' ');
            if (ad4_unbound_model==Unbound_Default) ad4_unbound_model = Unbound_Same_As_Bound;
            pr(logFile, "Unbound model to be used is %s.\n", report_parameter_library());

#ifdef DEBUG
            pr( logFile, "\nAbout to call evaluate.setup(), sInit:\n\n");
            printState( logFile, sInit, 2 );
#endif

            evaluate.setup( crd, charge, abs_charge, qsp_abs_charge, type, natom, map,
                            elec, emap, nonbondlist, ad_energy_tables, Nnb,
                            B_calcIntElec, B_isGaussTorCon, B_isTorConstrained,
                            B_ShowTorE, US_TorE, US_torProfile, vt, tlist, crdpdb, crdreo, sInit, ligand,
                            ignore_inter,
                            B_include_1_4_interactions, scale_1_4,
                            unbound_internal_FE, info, B_use_non_bond_cutoff, B_have_flexible_residues);
                            //parameterArray, unbound_internal_FE, info, B_use_non_bond_cutoff, B_have_flexible_residues);

            evaluate.compute_intermol_energy(TRUE);

            if(write_stateFile){
              fprintf(stateFile,"\t<run_requested>%d</run_requested>\n",nruns);
              fprintf(stateFile,"\t<runs>\n");
            }
            for (j=0; j<nruns; j++) {
                j1=j+1;

                (void) fprintf( logFile, "\n\n\tBEGINNING LAMARCKIAN GENETIC ALGORITHM DOCKING\n");
                (void) fflush( logFile );
                pr( logFile, "\nRun:\t%d / %d\n", j1, nruns );

                pr(logFile, "Date:\t");
                printdate( logFile, 2 );
                (void) fflush( logFile );

                gaStart = times( &tms_gaStart );

                if (B_reorient_random == TRUE) {
                    // create a new random orientation before docking
                    create_random_orientation( &q_reorient );
                    // reorient the ligand
                    reorient( logFile, true_ligand_atoms, atomstuff, crdpdb, charge, type,
                              parameterArray, q_reorient, origin, ntor, tlist, vt, &ligand, debug );
                    // update the evaluate object
                    evaluate.update_crds( crdpdb, vt );
                }

                //  Can get rid of the following line
                ((Genetic_Algorithm *)GlobalSearchMethod)->initialize(pop_size, 7+sInit.ntor);

                // Reiterate output level...
                pr(logFile, "Output level is set to %d.\n\n", outlev);

                // Start Lamarckian GA run -- Bound simulation
                sHist[nconf] = call_glss( GlobalSearchMethod, LocalSearchMethod,
                                          sInit,
                                          num_evals, pop_size,
                                          outlev,
                                          outputEveryNgens, &ligand,
                                          B_RandomTran0, B_RandomQuat0, B_RandomDihe0,
                                          info, FN_pop_file, end_of_branch );
                // State of best individual at end of GA-LS run is returned.
                // Finished Lamarckian GA run

                gaEnd = times( &tms_gaEnd );
                pr( logFile, "\nRun completed;  time taken for this run:\n");
                timesyshms( gaEnd - gaStart, &tms_gaStart, &tms_gaEnd );
                pr( logFile, "\n");
                printdate( logFile, 1 );
                (void) fflush( logFile );

                pr(logFile, "Total number of Energy Evaluations: %lu\n", evaluate.evals() );
                pr(logFile, "Total number of Generations:        %u\n", ((Genetic_Algorithm *)GlobalSearchMethod)->num_generations());

                pr( logFile, "\n\n\tFINAL LAMARCKIAN GENETIC ALGORITHM DOCKED STATE\n" );
                pr( logFile,     "\t_______________________________________________\n\n\n" );

                writePDBQT( j, seed,  FN_ligand, dock_param_fn, lig_center,
                            sHist[nconf], ntor, &eintra, &einter, natom, atomstuff,
                            crd, emap, elec,
                            charge, abs_charge, qsp_abs_charge,
                            ligand_is_inhibitor,
                            torsFreeEnergy,
                            vt, tlist, crdpdb, nonbondlist,
                            ad_energy_tables,
                            type, Nnb, B_calcIntElec,
                            map,
                            outlev, ignore_inter,
                            B_include_1_4_interactions, scale_1_4, parameterArray, unbound_internal_FE,
                            info, DOCKED, PDBQT_record, B_use_non_bond_cutoff, B_have_flexible_residues, ad4_unbound_model);

                // See also "calculateEnergies.cc", switch(ad4_unbound_model)
                if (ad4_unbound_model == Unbound_Same_As_Bound) {
                    // Update the unbound internal energy, setting it to the current internal energy
                    unbound_internal_FE = eintra;
                }
                econf[nconf] = eintra + einter + torsFreeEnergy - unbound_internal_FE;

                ++nconf;

                pr( logFile, UnderLine );
            } // Next LGA run
            if(write_stateFile){
               fprintf(stateFile,"\t</runs>\n");
               (void) fflush(stateFile);
            }
            (void) fflush(logFile);
        break;

//______________________________________________________________________________

    case DPF_LS:
       // do_local_only
       (void) sscanf(line, "%*s %d", &nruns);

            if ( nruns > MAX_RUNS ) {

               prStr( error_message, "%s:  ERROR: %d runs requested, but only dimensioned for %d.\nChange \"MAX_RUNS\" in \"constants.h\".", programname, nruns, MAX_RUNS);
               stop( error_message );
               exit( -1 );

           } else if (LocalSearchMethod==NULL) {

               prStr(error_message, "%s:  ERROR:  You must use \"set_sw1\", \"set_psw1\" or \"set_pattern\" to create a Local Optimization object.\n", programname);
               stop(error_message);
               exit(-1);
            }
           exit_if_missing_elecmap_desolvmap_about("ls");
           pr( logFile, "Number of Local Search (LS) only dockings = %d run%c\n", nruns, (nruns > 1)?'s':' ');
           if (ad4_unbound_model==Unbound_Default) ad4_unbound_model = Unbound_Same_As_Bound;
           pr(logFile, "Unbound model to be used is %s.\n", report_parameter_library());
           evaluate.setup( crd, charge, abs_charge, qsp_abs_charge, type, natom, map,
                           elec, emap,
                           nonbondlist,
                           ad_energy_tables,
                           Nnb, B_calcIntElec, B_isGaussTorCon,B_isTorConstrained,
                           B_ShowTorE, US_TorE, US_torProfile, vt, tlist, crdpdb, crdreo, sInit, ligand,
                           ignore_inter,
                           B_include_1_4_interactions, scale_1_4,
                           unbound_internal_FE, info, B_use_non_bond_cutoff, B_have_flexible_residues);
                           //parameterArray, unbound_internal_FE, info, B_use_non_bond_cutoff, B_have_flexible_residues);

            evaluate.compute_intermol_energy(TRUE);

           if(write_stateFile){
             fprintf(stateFile,"\t<run_requested>%d</run_requested>\n",nruns);
             fprintf(stateFile,"\t<runs>\n");
           }

           for (j=0; j<nruns; j++) {

               (void) fprintf( logFile, "\n\n\tBEGINNING SOLIS & WETS LOCAL SEARCH DOCKING\n");
               pr( logFile, "\nRun:\t%d / %d\n", j+1, nruns );
               (void) fflush( logFile );

               pr(logFile, "Date:\t");
               printdate( logFile, 2 );
               (void) fflush( logFile );

               gaStart = times( &tms_gaStart );

               sHist[nconf] = call_ls(LocalSearchMethod, sInit, pop_size, &ligand);

               pr(logFile, "There were %lu Energy Evaluations.\n", evaluate.evals());

               gaEnd = times( &tms_gaEnd );
               pr( logFile, "Time taken for this Local Search (LS) run:\n");
               timesyshms( gaEnd - gaStart, &tms_gaStart, &tms_gaEnd );
               pr( logFile, "\n");
               (void) fflush( logFile );

               pr( logFile, "\n\n\tFINAL LOCAL SEARCH DOCKED STATE\n" );
               pr( logFile,     "\t_______________________________\n\n\n" );

               writePDBQT( j, seed, FN_ligand, dock_param_fn, lig_center,
                           sHist[nconf], ntor, &eintra, &einter, natom, atomstuff,
                           crd, emap, elec,
                           charge, abs_charge, qsp_abs_charge,
                           ligand_is_inhibitor,
                           torsFreeEnergy,
                           vt, tlist, crdpdb, nonbondlist,
                           ad_energy_tables,
                           type, Nnb, B_calcIntElec,
                           map,
                           outlev,
                           ignore_inter,
                           B_include_1_4_interactions, scale_1_4, parameterArray, unbound_internal_FE,
                           info, DOCKED, PDBQT_record, B_use_non_bond_cutoff, B_have_flexible_residues, ad4_unbound_model);

               // See also "calculateEnergies.cc", switch(ad4_unbound_model)
               if (ad4_unbound_model == Unbound_Same_As_Bound) {
                   // Update the unbound internal energy, setting it to the current internal energy
                   unbound_internal_FE = eintra;
               }
               econf[nconf] = eintra + einter + torsFreeEnergy - unbound_internal_FE;

               ++nconf;

               pr( logFile, UnderLine );

           } // Next run
           if(write_stateFile){
             fprintf(stateFile,"\t</runs>\n");
             (void) fflush(stateFile);
           }
           (void) fflush(logFile);
       break;

//______________________________________________________________________________

    case DPF_GS:
      (void) sscanf(line, "%*s %d", &nruns);

          if (nruns>MAX_RUNS) {

              prStr(error_message, "%s:  ERROR:  %d runs requested, but only dimensioned for %d.\nChange \"MAX_RUNS\" in \"constants.h\".", programname, nruns, MAX_RUNS);
              stop(error_message);
              exit(-1);

          } else if (GlobalSearchMethod==NULL) {
              prStr(error_message, "%s:  ERROR:  You must use \"set_ga\" to allocate a Global Optimization object.\n", programname);
              stop(error_message);
              exit(-1);
           }

          exit_if_missing_elecmap_desolvmap_about("gs");

          pr(logFile, "Number of Genetic Algorithm (GA) only dockings = %d run%c\n", nruns, (nruns>1)?'s':' ');
          if (ad4_unbound_model==Unbound_Default) ad4_unbound_model = Unbound_Same_As_Bound;
          pr(logFile, "Unbound model to be used is %s.\n", report_parameter_library());


          evaluate.setup( crd, charge, abs_charge, qsp_abs_charge, type, natom, map,
                          elec, emap,
                          nonbondlist,
                          ad_energy_tables,
                          Nnb, B_calcIntElec, B_isGaussTorCon,B_isTorConstrained,
                          B_ShowTorE, US_TorE, US_torProfile, vt, tlist, crdpdb, crdreo, sInit, ligand,
                          ignore_inter,
                          B_include_1_4_interactions, scale_1_4,
                          unbound_internal_FE, info, B_use_non_bond_cutoff, B_have_flexible_residues);
                          //parameterArray, unbound_internal_FE, info, B_use_non_bond_cutoff, B_have_flexible_residues);

            evaluate.compute_intermol_energy(TRUE);

          if(write_stateFile){
            fprintf(stateFile,"\t<run_requested>%d</run_requested>\n",nruns);
            fprintf(stateFile,"\t<runs>\n");
          }
          for (j=0; j<nruns; j++) {

              fprintf( logFile, "\n\n\tBEGINNING GENETIC ALGORITHM DOCKING\n");
              pr(logFile, "\nDoing Genetic Algorithm run:  %d/%d.\n", j+1, nruns);
              (void) fflush( logFile );

              pr(logFile, "Date:\t");
              printdate( logFile, 2 );
              (void) fflush( logFile );

              gaStart = times(&tms_gaStart);

              sHist[nconf] = call_gs( GlobalSearchMethod, sInit, num_evals, pop_size,
                                      &ligand, outputEveryNgens, info, end_of_branch);

              pr(logFile, "\nFinal docked state:\n");
              printState(logFile, sHist[nconf], 2);

              gaEnd = times(&tms_gaEnd);
              pr(logFile, "Time taken for this Genetic Algorithm (GA) run:\n");
              timesyshms(gaEnd-gaStart, &tms_gaStart, &tms_gaEnd);
              pr(logFile, "\n");
              (void) fflush(logFile);

              pr(logFile, "Total number of Energy Evaluations: %lu\n", evaluate.evals() );
              pr(logFile, "Total number of Generations:        %u\n", ((Genetic_Algorithm *)GlobalSearchMethod)->num_generations());


              pr( logFile, "\n\n\tFINAL GENETIC ALGORITHM DOCKED STATE\n" );
              pr( logFile,     "\t____________________________________\n\n\n" );

              writePDBQT( j, seed, FN_ligand, dock_param_fn, lig_center,
                          sHist[nconf], ntor, &eintra, &einter, natom, atomstuff,
                          crd, emap, elec,
                          charge, abs_charge, qsp_abs_charge,
                          ligand_is_inhibitor,
                          torsFreeEnergy,
                          vt, tlist, crdpdb, nonbondlist,
                          ad_energy_tables,
                          type, Nnb, B_calcIntElec,
                          map,
                          outlev,
                          ignore_inter,
                          B_include_1_4_interactions, scale_1_4, parameterArray, unbound_internal_FE,
                          info, DOCKED, PDBQT_record, B_use_non_bond_cutoff, B_have_flexible_residues, ad4_unbound_model);

              // See also "calculateEnergies.cc", switch(ad4_unbound_model)
              if (ad4_unbound_model == Unbound_Same_As_Bound) {
                  // Update the unbound internal energy, setting it to the current internal energy
                  unbound_internal_FE = eintra;
              }
              econf[nconf] = eintra + einter + torsFreeEnergy - unbound_internal_FE;

              ++nconf;

              pr( logFile, UnderLine );

          } // Next run
          if(write_stateFile){
            fprintf(stateFile,"\t</runs>\n");
            (void) fflush(stateFile);
          }
          (void) fflush(logFile);
      break;

//______________________________________________________________________________

    case GA_pop_size:
       (void) sscanf(line, "%*s %u", &pop_size);
       pr(logFile, "A population of %u individuals will be used\n", pop_size);
        (void) fflush(logFile);
       break;

//______________________________________________________________________________

    case GA_num_generations:
       (void) sscanf(line, "%*s %u", &num_generations);
       pr(logFile, "The GA will run for at most %u generations.\n", num_generations);
        (void) fflush(logFile);
       break;

//______________________________________________________________________________

    case GA_num_evals:
       (void) sscanf(line, "%*s %u", &num_evals);
       pr(logFile, "There will be at most %u function evaluations used.\n", num_evals);
        (void) fflush(logFile);
       break;

//______________________________________________________________________________

    case GA_window_size:
       (void) sscanf(line, "%*s %d", &window_size);
       pr(logFile, "The GA's selection window is %d generations.\n", window_size);
        (void) fflush(logFile);
       break;

//______________________________________________________________________________

    case GA_low:
       (void) sscanf(line, "%*s %d", &low);
       pr(logFile, "Setting low to %d.\n", low);
        (void) fflush(logFile);
       break;

//______________________________________________________________________________

    case GA_high:
       (void) sscanf(line, "%*s %d", &high);
       pr(logFile, "Setting high to %d.\n", high);
        (void) fflush(logFile);
       break;

//______________________________________________________________________________

    case GA_elitism:
       (void) sscanf(line, "%*s %d", &elitism);
       pr(logFile, "The %d best will be preserved each generation.\n", elitism);
        (void) fflush(logFile);
       break;

//______________________________________________________________________________

    case GA_mutation_rate:
       (void) sscanf(line, "%*s " FDFMT, &m_rate);
       pr(logFile, "The mutation rate is %f.\n", m_rate);
        (void) fflush(logFile);
      // if m_rate is out of range, make_table will fail 
      if (m_rate < 0 || m_rate > 1) {
	prStr(error_message, "mutation rate must be within range 0 to 1 (inclusive).\n");
	stop(error_message);
	}
       break;

//______________________________________________________________________________

    case GA_crossover_rate:
       (void) sscanf(line, "%*s " FDFMT, &c_rate);
       pr(logFile, "The crossover rate is %f.\n", c_rate);
        (void) fflush(logFile);
       break;

//______________________________________________________________________________

    case GA_Cauchy_alpha:
       (void) sscanf(line, "%*s " FDFMT, &alpha);
       pr(logFile, "The alpha parameter (for the Cauchy distribution) is being set to %f.\n",
          alpha);
        (void) fflush(logFile);
       break;

//______________________________________________________________________________

    case GA_Cauchy_beta:
       (void) sscanf(line, "%*s " FDFMT, &beta);
       pr(logFile, "The beta parameter (for the Cauchy distribution) is being set to %f.\n",
          beta);
        (void) fflush(logFile);
       break;

//______________________________________________________________________________

    case SW_max_its:
       (void) sscanf(line, "%*s %u", &max_its);
       pr(logFile, "Solis & Wets algorithms will perform at most %u iterations.\n", max_its);
        (void) fflush(logFile);
       break;

//______________________________________________________________________________

    case SW_max_succ:
       (void) sscanf(line, "%*s %u", &max_succ);
       pr(logFile, "Solis & Wets algorithms expand rho every %u in a row successes.\n", max_succ);
        (void) fflush(logFile);
      break;

//______________________________________________________________________________

    case SW_max_fail:
       (void) sscanf(line, "%*s %u", &max_fail);
       pr(logFile, "Solis & Wets algorithms contract rho every %u in a row failures.\n", max_fail);
        (void) fflush(logFile);
      break;

//______________________________________________________________________________

    case SW_rho:
       (void) sscanf(line, "%*s " FDFMT, &rho);
       pr(logFile, "rho is set to %f.\n", rho);
        (void) fflush(logFile);
      break;


//______________________________________________________________________________

    case SW_lb_rho:
        (void) sscanf(line, "%*s " FDFMT, &lb_rho);
        pr(logFile, "rho will never get smaller than %f.\n", lb_rho);
        (void) fflush(logFile);
        break;

//______________________________________________________________________________

    case PSW_TRANS_SCALE:
        (void) sscanf(line, "%*s " FDFMT, &psw_trans_scale);
        pr(logFile, "psw_trans_scale is set to %f.\n", psw_trans_scale);
        (void) fflush(logFile);
        break;

//______________________________________________________________________________

    case PSW_ROT_SCALE:
        (void) sscanf(line, "%*s " FDFMT, &psw_rot_scale);
        pr(logFile, "psw_rot_scale is set to %f.\n", psw_rot_scale);
        (void) fflush(logFile);
        break;

//______________________________________________________________________________

    case PSW_TORS_SCALE:
        (void) sscanf(line, "%*s " FDFMT, &psw_tors_scale);
        pr(logFile, "psw_tors_scale is set to %f.\n", psw_tors_scale);
        (void) fflush(logFile);
        break;

//______________________________________________________________________________


    case LS_search_freq:
        (void) sscanf(line, "%*s " FDFMT, &search_freq);
        pr(logFile, "Local search will be performed with frequency %f.\n", search_freq);
        (void) fflush(logFile);
        break;

//______________________________________________________________________________
   
   case PSO_c1:
        (void) sscanf(line, "%*s %lf", &c1);
        pr(logFile, "PSO will be performed with the First Confidence Coefficient  (C1) %lf.\n", c1);
        (void) fflush(logFile);
        break;

//______________________________________________________________________________

   case PSO_c2:
        (void) sscanf(line, "%*s %lf", &c2);
        pr(logFile, "PSO will be performed with the Second Confidence Coefficient (C2) %lf.\n", c2);
        (void) fflush(logFile);
        break;

//______________________________________________________________________________

   case PSO_k:
        (void) sscanf(line, "%*s %d", &K);
        pr(logFile, "Max number of particles informed by a given one = %d.\n", K);
        (void) fflush(logFile);
        break;

//______________________________________________________________________________

   case PSO_swarm_moves:
        (void) sscanf(line, "%*s %d", &eval_max);
        pr(logFile, "There will be %d swarm Moves.\n", eval_max);
        (void) fflush(logFile);
        break;

//______________________________________________________________________________

   case PSO_swarm_size_factor:
        (void) sscanf(line, "%*s %d", &S_factor);
        pr(logFile, "There will be %d Swarm Size Factor.\n", S_factor);
        (void) fflush(logFile);
        break;

//______________________________________________________________________________

   case PSO_n_exec:
        (void) sscanf(line, "%*s %d", &n_exec_max);
        pr(logFile, "Number of requested PSO runs = %d.\n", n_exec_max);
        (void) fflush(logFile);
        break;

//______________________________________________________________________________

   case DPF_PSO_CONSTRICTION:

	/* 
	 * PSO Work 
	 */
        pr( logFile, "\nTotal number of torsions in the ligand = %d \n", ntor);
        D = 7 + ntor; //Dimension D 
        pr( logFile, "\nTotal number of dimensions is equal to the number of Degrees of Freedom = %d \n", D);
        pr( logFile, "Number of requested PSO dockings = %d run%c\n", n_exec_max, (n_exec_max > 1)?'s':' ');
        //No. of Particles in the Swarm 
        S = S_factor + (int) (2 * sqrt(D)); //S - Swarm Size (~somewhat analogous to pop_size)
        pr( logFile, "\nTotal number in swarm size = %d \n", S);
        
        //Checking the No. of Runs and No. of Particles 
        if (n_exec_max > R_max ) //set to 1024 in constants.h
        {
		    prStr( error_message, "ERROR: %d runs requested, but only dimensioned for %d.\nChange \"R_max\" in \"constants.h\".", nruns, R_max);
            stop( error_message );
            exit( -1 );
        } else if ( S  > S_max ) 
        {
		    prStr( error_message, "ERROR: %d Particles initialised in Swarm, but only dimensioned for %d.\nChange \"S_max\" in \"constants.h\".", S, S_max);
            stop( error_message );
            exit( -1 );
        }

        if (outlev > 3) {
            pr( logFile, "\ncalling initialiseDimension: info->lo=%f %f %f , info->hi=%f %f %f,  \n", info->lo[0], info->lo[1], info->lo[2], info->hi[0], info->hi[1], info->hi[1]);
	    };
        initialiseDimension(info, xmin, xmax, D);
        if (outlev > 3) {
            pr( logFile, "\nAFTER initDim trans: xmin=%f %f %f , xmax=%f %f %f,  \n", xmin[0], xmin[1], xmin[2], xmax[0], xmax[1], xmax[1]);
            pr( logFile, "\nAFTER initDim quat: xmin=%f %f %f %f, xmax=%f %f %f %f,  \n", xmin[3], xmin[4], xmin[5], xmin[6], xmax[3], xmax[4], xmax[5],xmax[6]);
        };

        evaluate.setup( crd, charge, abs_charge, qsp_abs_charge, type, natom, map,
                        elec, emap, nonbondlist, ad_energy_tables, Nnb,
                        B_calcIntElec, B_isGaussTorCon, B_isTorConstrained,
                        B_ShowTorE, US_TorE, US_torProfile, vt, tlist, crdpdb, crdreo, sInit, ligand,
                        ignore_inter, B_include_1_4_interactions, scale_1_4,
                        unbound_internal_FE, info, B_use_non_bond_cutoff, B_have_flexible_residues);
                        //parameterArray, unbound_internal_FE, info, B_use_non_bond_cutoff, B_have_flexible_residues);

        evaluate.compute_intermol_energy(TRUE);

//	//Initialise the Dimension of xmax, xmin using the function
//initialiseDimension(xlo, xhi, ylo, yhi, zlo, zhi, xmin, xmax, D);


//    evaluate.setup(crd, charge, type, natom, map, inv_spacing, 
//              elec, emap,
//              xlo, xhi, ylo, yhi, zlo, zhi, nonbondlist, e_internal, Nnb, 
//              B_calcIntElec, q1q2, B_isGaussTorCon, B_isTorConstrained,
//              B_ShowTorE, US_TorE, US_torProfile, vt, tlist, crdpdb, sInit, mol, 
//              B_template, template_energy, template_stddev);

        for (n_exec = 0; n_exec < n_exec_max; n_exec++)
        {
            j1 = n_exec + 1;
            (void) fprintf( logFile, "\n\n\tBEGINNING PARTICLE SWARM OPTIMIZATION\n");
            (void) fflush( logFile );
            pr(logFile, "\nRun:\t%d / %d\n", j1, n_exec_max );
            pr(logFile, "\nWorking with Constriction PSO\n");
            pr(logFile, "Date:\t");
            printdate( logFile, 2 );
		    (void) fflush( logFile );
            psoStartClock = times(&tms_psoStartClock);
// Reiterate output level...
//pr(logFile, "Output level is set to %d.... about to invoke call_cpso:n_exec=%d, S=%d, D=%d, 
//   *xmin= %f, *xmax=%f, eval_max=%d, K=%f\n\n", outlev, n_exec, S, D, *xmin, *xmax, eval_max, K);
//Start Particle Swarm Optimization Run	               
//sHist[n_exec] = call_cpso(n_exec, sInit, S, D, xmin, xmax, eval_max, K, c1, c2, outlev,
//					7+sInit.ntor, max_its, max_succ, max_fail, 2.0, 0.5, search_freq, rho_ptr, lb_rho_ptr);	
// swarmsize factor S is analogous to pop_size
            sHist[n_exec] = call_cpso(LocalSearchMethod, sInit, n_exec, S, D, xmin, xmax, eval_max, K, c1, c2, outlev);
            //Finished Particle Swarm Optimization Run
            psoEndClock = times(&tms_psoEndClock);
            pr(logFile, "\nRun completed. Time taken for this PSO run:\n");
            timesyshms(psoEndClock-psoStartClock, &tms_psoStartClock, &tms_psoEndClock);
            pr(logFile, "\n");
            printdate( logFile, 1 );
            (void) fflush(logFile);
            pr(logFile, "\nTotal number of Energy Evaluations: %lu\n", evaluate.evals() );
            pr(logFile, "Total number of Swarm moves 	  %d\n", eval_max);
            pr( logFile, "\n\n\tFINAL PARTICLE SWARM OPTIMIZATION (varCPSO-ls) DOCKED STATE\n" );
            pr( logFile,     "\t_________________________________________________________________\n\n\n" );
//Printing the Run in PDBQ Format 
//		    writeStateOfPDBQ(n_exec , FN_ligand, dock_param_fn, lig_center, 
//        	            &(sHist[n_exec]), ntor, &eintra, &einter, natom, atomstuff, 
//                	    crd, emap, elec, charge, 
//	                    ligand_is_inhibitor,
//        	            torsFreeEnergy,
//                	    vt, tlist, crdpdb, nonbondlist, e_internal,
//	                    type, Nnb, B_calcIntElec, q1q2,
//        	            map, inv_spacing, xlo, ylo, zlo, xhi, yhi, zhi,
//                	    B_template, template_energy, template_stddev,
//	                    outlev);

            writePDBQT( n_exec, seed,  FN_ligand, dock_param_fn, lig_center,
                       sHist[nconf], ntor, &eintra, &einter, natom, atomstuff,
                       crd, emap, elec, charge, 
                       abs_charge, qsp_abs_charge,
                       ligand_is_inhibitor,
                       torsFreeEnergy,
                       vt, tlist, crdpdb, nonbondlist,
                       ad_energy_tables,
                       type, Nnb, B_calcIntElec, map,
                       outlev, ignore_inter, B_include_1_4_interactions, 
                       scale_1_4, parameterArray, unbound_internal_FE,
                       info, DOCKED, PDBQT_record, B_use_non_bond_cutoff, //@@info
                       B_have_flexible_residues, ad4_unbound_model);
	
            econf[nconf] = eintra + einter; // new2
            ++nconf;
            pr( logFile, UnderLine );
        }
        /* PSO Work*/
        break;

//______________________________________________________________________________

    case DPF_ANALYSIS:
        /*
        ** analysis
        */
        /* _____________________________________________________________________
        **
        ** Perform Cluster analysis on results of docking,
        ** _____________________________________________________________________
        */
            analysis( Nnb, atomstuff, charge, abs_charge, qsp_abs_charge, B_calcIntElec, clus_rms_tol,
                      crdpdb, ad_energy_tables, map, econf, nruns,
                      natom, nonbondlist, nconf, ntor, sHist, FN_ligand,
                      lig_center, B_symmetry_flag, tlist, type, vt, FN_rms_ref_crds,
                      torsFreeEnergy, B_write_all_clusmem, ligand_is_inhibitor,
                      outlev,
                      ignore_inter, B_include_1_4_interactions, scale_1_4,
                      unbound_internal_FE,
                      info, B_use_non_bond_cutoff, B_have_flexible_residues,
                      B_rms_atoms_ligand_only, ad4_unbound_model);

            (void) fflush(logFile);
        break;

//______________________________________________________________________________

    case DPF_TORSDOF:
        /*
        ** torsdof %d %f
        */
        retval = sscanf( line, "%*s %d " FDFMT, &ntorsdof, &torsdoffac );
        if (retval == 2) {
            pr( logFile, "WARNING:  The torsional DOF coefficient is now read in from the parameter file; the value specified here (%.4lf) will be ignored.\n\n", (double)torsdoffac);
        }
        pr( logFile, "Number of torsional degrees of freedom = %d\n", ntorsdof);
        pr( logFile, "Free energy coefficient for torsional degrees of freedom = %.4f", AD4.coeff_tors);
        if (parameter_library_found == 1) {
            pr( logFile, " as specified in parameter library \"%s\".\n\n", FN_parameter_library );
        } else {
            pr( logFile, ", the factory default value.\n\n");
        }

        torsFreeEnergy = (Real)ntorsdof * AD4.coeff_tors;

        pr( logFile, "Estimated loss of torsional free energy upon binding = %+.4f kcal/mol\n\n", torsFreeEnergy);
        (void) fflush(logFile);
        break;

//______________________________________________________________________________

    case DPF_INVESTIGATE:
        /*
        ** Bin energies by RMSD from reference structure
        **
        ** investigate 100000 1000000 100
        */
        (void) sscanf( line, "%*s %d %d %d", &OutputEveryNTests, &maxTests, &NumLocalTests );
        (void) fprintf( logFile, "OutputEveryNTests= %d\n", OutputEveryNTests);
        (void) fprintf( logFile, "maxTests= %d\n", maxTests );
        (void) fprintf( logFile, "NumLocalTests= %d\n\n", NumLocalTests );
        (void) investigate( Nnb, charge, abs_charge, qsp_abs_charge, B_calcIntElec,
                            crd, crdpdb, ad_energy_tables,
                            maxTests,
                            map, natom, nonbondlist, ntor,
                            outlev, tlist, type, vt, B_isGaussTorCon, US_torProfile,
                            B_isTorConstrained, B_ShowTorE, US_TorE,
                            F_TorConRange, N_con, B_symmetry_flag, FN_rms_ref_crds,
                            OutputEveryNTests, NumLocalTests, trnStep0, torStep0,
                            ignore_inter,
                            B_include_1_4_interactions, scale_1_4,
                            unbound_internal_FE,
                            info, B_use_non_bond_cutoff, B_have_flexible_residues );

        (void) fflush(logFile);
        break;

//______________________________________________________________________________

    case DPF_LIG_NOT_INHIB:
        /*
        ** ligand_is_not_inhibitor
        */
        ligand_is_inhibitor = 0;
        pr( logFile, "\nThis ligand is not an inhibitor, so dissociation constants (Kd) will be calculated, not inhibition constants (Ki).\n\n" );
        (void) fflush(logFile);
        break;

/*____________________________________________________________________________*/

    case DPF_UNBOUND_MODEL:
        /*
        **  unbound_model { extended [energy <FLOAT>]| compact | bound }
        **    extended is alias for "compute_unbound_extended" token
        */
        char unbound_model_type[LINE_LEN];
        (void) sscanf( line, "%*s %s", unbound_model_type );

        if (equal( unbound_model_type, "bound", 5 )
        || equal( unbound_model_type, "same_as_bound", 13 )
        || equal( unbound_model_type, "unbound_same_as_bound", 21 )) {
            if (ad4_unbound_model != Unbound_Same_As_Bound)  // default for Autodock 4.1
                setup_parameter_library(outlev, "Unbound_Same_As_Bound", Unbound_Same_As_Bound);
            ad4_unbound_model = Unbound_Same_As_Bound;
        } else if (equal( unbound_model_type, "extended", 8 )) {
            if (ad4_unbound_model != Unbound_Default) { //illegal to set extended after other
                pr( logFile, "%s:  ERROR:  Setting unbound model type twice: \"%s\" .\n",
                    programname, unbound_model_type );
                stop("");
            }
            if ( (1== sscanf( line, "%*s extended energy " FDFMT, &unbound_internal_FE ))){
                ad4_unbound_model = Extended;
                setup_parameter_library(outlev, "unbound_extended", ad4_unbound_model);
            }
            else goto process_DPF_COMPUTE_UNBOUND_EXTENDED; // case DPF_COMPUTE_UNBOUND_EXTENDED below
        } else if (equal( unbound_model_type, "compact", 6 )) {
            ad4_unbound_model = Compact;
        } else {
            // note that "User" is not acceptable in dpf file
            pr( logFile, "%s:  ERROR:  Unrecognized unbound model type \"%s\" .\n",
                    programname, unbound_model_type );
            stop("");
        }
        (void) fflush(logFile);
        break;

/*____________________________________________________________________________*/

    case DPF_UNBOUND:
        /*
         * unbound FLOAT
         * unbound energy FLOAT
         */
        if (ad4_unbound_model != Unbound_Default && ad4_unbound_model!= User) { //illegal to set user after other
            pr( logFile, "%s:  ERROR:  Setting unbound model type twice!\n",
                programname );
            stop("");
        }
        if ((1!= sscanf( line, "%*s " FDFMT, &unbound_internal_FE ))
        && (1!= sscanf( line, "%*s energy" FDFMT, &unbound_internal_FE ))){
            pr( logFile, "%s:  ERROR:  Non-numeric unbound model energy \"%s\" .\n",
                    programname, line);
            stop("");
        }
        pr(logFile, "The internal energy of the unbound state was set to %+.3lf kcal/mol\n", unbound_internal_FE);
        ad4_unbound_model = User;
        pr(logFile, "The unbound ligand energy model was set to User\n\n");
        (void) fflush(logFile);
        break;

/*____________________________________________________________________________*/
    case DPF_COMPUTE_UNBOUND_EXTENDED:
        /*
         *  compute_unbound_extended
         */
    process_DPF_COMPUTE_UNBOUND_EXTENDED:
        if (ntor > 0) {
           (void) sprintf( message, "%s: WARNING: Using autodock4.0 unbound extended model in autodock4.2!\n", programname );
            print_2x( logFile, stderr, message );
            if (ad4_unbound_model != Unbound_Default) { //illegal to set extended after other
                pr( logFile, "%s:  ERROR:  Setting unbound model type twice!\n",
                    programname );
                stop("");
            }
            ad4_unbound_model = Extended;
            setup_parameter_library(outlev, "unbound_extended", ad4_unbound_model);

            pr(logFile, "Computing the energy of the unbound state of the ligand,\ngiven the torsion tree defined in the ligand file.\n\n");
            (void) fflush( logFile );

            // The initial goal is to obtain an extended conformation of the ligand.

            // Step 0 // {
            //
            // Set termination criteria for unbound calculations
            // -------------------------------------------------
            //
            // Set the maximum number of energy evaluations for finding the unbound conformation
            // if num_evals is less than this, then a shorter unbound docking will be performed
            max_evals_unbound = 1000000; // 1 million
            num_evals_unbound = num_evals > max_evals_unbound ?  max_evals_unbound :  num_evals;
            // end of Step 0 // }

            // Step 1 // {
            //
            // Run a hybrid global-local search using the unbound energy tables (set to be repulsive-only)
            // -------------------------------------------------------------------------------------------
            //
            //  *  Turn off the use of the non-bond cutoff
            //  *  Turn off internal electrostatics
            //  *  Turn off intermolecular energy calculations
            // TODO Do not translate or rotate the ligand in unbound searches
            //
            /*
             *  Genetic Algorithm-Local search,  a.k.a.
             *  Lamarckian Genetic Algorithm
             */
            if ((GlobalSearchMethod==NULL)||(LocalSearchMethod==NULL)) {
                prStr(error_message, "%s:  ERROR:  You must use \"set_ga\" to allocate both Global Optimization object AND Local Optimization object.\n", programname);
                stop(error_message);
                exit(-1);
            }
            exit_if_missing_elecmap_desolvmap_about("compute_unbound_extended");

            //
            // Do not use a non-bond cutoff, this helps to produce the "most" extended conformation
            // especially with long inhibitors
            B_use_non_bond_cutoff = FALSE;
            //
            // Save the current value of B_calcIntElec, so we can restore it later.
            B_calcIntElec_saved = B_calcIntElec;
            //
            // Set the calculation of internal electrostatics to FALSE:
            // B_calcIntElec = FALSE;
            //
            // Assume the unbound state of the receptor is the same as the input coordinates from the
            // flexible residues file.  This means we must not change the rotatable bonds in the
            // flexible residues of the receptor during the unbound extended search.
            // We can turn off rotation of the flexres by setting ntor to ntor_ligand.
            // Save the current value of "ntor" in the "sInit" state variable, set it to number of torsions
            // in the ligand for the unbound extended search, then restore it afterwards.
            saved_sInit_ntor = sInit.ntor;
            sInit.ntor = ntor_ligand;
            //
            // Use the repulsive unbound energy tables, "unbound_energy_tables",
            // to drive the molecule into an extended conformation
            evaluate.setup( crd, charge, abs_charge, qsp_abs_charge, type, natom, map,
                            elec, emap, nonbondlist, unbound_energy_tables, Nnb,
                            B_calcIntElec, B_isGaussTorCon, B_isTorConstrained,
                            B_ShowTorE, US_TorE, US_torProfile, vt, tlist, crdpdb, crdreo, sInit, ligand,
                            ignore_inter,
                            B_include_1_4_interactions, scale_1_4,
                            unbound_internal_FE, info, B_use_non_bond_cutoff, B_have_flexible_residues);
                            //parameterArray, unbound_internal_FE, info, B_use_non_bond_cutoff, B_have_flexible_residues);
            //
            // Turn off computing the intermolecular energy, we will only consider the intramolecular energy
            // to determine the unbound state of the flexible molecule:
            evaluate.compute_intermol_energy(FALSE);
            //
            (void) fprintf( logFile, "\n\n\tBEGINNING COMPUTATION OF UNBOUND EXTENDED STATE USING LGA\n");
            (void) fprintf( logFile,     "\t_________________________________________________________\n\n\n");
            (void) fflush( logFile );
            //
            pr(logFile, "Date:\t");
            printdate( logFile, 2 );
            (void) fflush( logFile );
            gaStart = times( &tms_gaStart );
            //  Can get rid of the following line
            ((Genetic_Algorithm *)GlobalSearchMethod)->initialize(pop_size, 7+sInit.ntor);
            //
            // Start Lamarckian GA run searching only torsions -- Unbound simulation
            // sUnbound_ext = call_glss_tors( GlobalSearchMethod, LocalSearchMethod,
            sUnbound_ext = call_glss( GlobalSearchMethod, LocalSearchMethod,
                                      sInit,
                                      num_evals_unbound, pop_size,
                                      outlev,
                                      outputEveryNgens, &ligand,
                                      // B_RandomDihe0, // use this line with call_glss_tors()
                                      B_RandomTran0, B_RandomQuat0, B_RandomDihe0,
                                      info, FN_pop_file, end_of_branch );
            // State of best individual at end of GA-LS run, sUnbound_ext, is returned.
            // Finished Lamarckian GA run
            gaEnd = times( &tms_gaEnd );
            pr( logFile, "\nFinished Lamarckian Genetic Algorithm (LGA), time taken:\n");
            timesyshms( gaEnd - gaStart, &tms_gaStart, &tms_gaEnd );
            pr( logFile, "\n");
            printdate( logFile, 1 );
            (void) fflush( logFile );
            pr(logFile, "\nTotal number of Energy Evaluations: %lu\n", evaluate.evals() );
            pr(logFile, "Total number of Generations:        %u\n", ((Genetic_Algorithm *)GlobalSearchMethod)->num_generations());
            // end of Step 1 // }

            // Step 2 // {
            //
            // Do a short local search using the standard internal energy tables
            // -----------------------------------------------------------------
            //
            // turn on internal electrostatics
            // but keep intermolecular energy calculations off
            //
            // Turn on calculation of internal electrostatics:
            //// B_calcIntElec = TRUE;
            //
            // Use the standard AutoDock energy tables to compute the internal energy
            // Use this value to set unbound_internal_FE
            //// evaluate.setup( crd, charge, abs_charge, qsp_abs_charge, type, natom, map,
                            //// elec, emap, nonbondlist, ad_energy_tables, Nnb,
                            //// B_calcIntElec, B_isGaussTorCon, B_isTorConstrained,
                            //// B_ShowTorE, US_TorE, US_torProfile, vt, tlist, crdpdb, crdreo, sInit, ligand,
                            //// ignore_inter,
                            //// B_include_1_4_interactions, scale_1_4,
                            //// unbound_internal_FE, info, B_use_non_bond_cutoff, B_have_flexible_residues);
            //
            // --- Start Local Search ---
            //// pr( logFile, "\nPerforming local search using standard AutoDock scoring function\n" );
            //// pr( logFile, "\nUsing UnboundLocalSearchMethod = new Solis_Wets1(7+sInit.ntor, 300, 4, 4, 1., 0.01, 2., 0.5, 1.);\n\n" );
            // Create a local search object
            // * Use an initial rho value of 0.1 (default is set in DPF by "sw_rho 1.0")
            //   to ensure smaller, 'more local' steps.
            // * Use a search frequency of 1.0 (default is set in DPF by "ls_search_freq 0.06")
            //// unsigned int ls_pop_size = 150;
            // max_its = 300
            // max_succ = 4
            // max_fail = 4
            // rho = 1.
            // lb_rho = 0.01
            // expansion = 2.
            // contraction = 0.5
            // search_freq = 1.
            //// UnboundLocalSearchMethod = new Solis_Wets1(7+sInit.ntor, 300, 4, 4, 1., 0.01, 2., 0.5, 1.);
            // Perform a local search, using the standard AutoDock 4 scoring function
            //// sUnbound_ls = call_ls( UnboundLocalSearchMethod, sUnbound_ext, ls_pop_size, &ligand );
            //// // sUnbound_ext = sUnbound_ls; // if you want to update sUnbound_ext to be sUnbound_ls...
            // --- Finished Local Search ---
            // end of Step 2 // }

            // Step 3 // {
            //
            // Restore the AutoDock 4 force field for docking
            // ----------------------------------------------
            //
            // Remember to turn on the use of the non-bond cutoff
            B_use_non_bond_cutoff = TRUE;
            //
            // Restore the setting for calculation of internal electrostatics to the saved value:
            B_calcIntElec = B_calcIntElec_saved;
            //
            // Restore the number of torsions in the state variables "sInit" and "sUnbound_ext"
            sInit.ntor = saved_sInit_ntor;
            sUnbound_ext.ntor = saved_sInit_ntor;
            //
            // Use the standard AutoDock energy tables to compute the internal energy
            // Use this value to set unbound_internal_FE
            evaluate.setup( crd, charge, abs_charge, qsp_abs_charge, type, natom, map,
                            elec, emap, nonbondlist, ad_energy_tables, Nnb,
                            B_calcIntElec, B_isGaussTorCon, B_isTorConstrained,
                            B_ShowTorE, US_TorE, US_torProfile, vt, tlist, crdpdb, crdreo, sInit, ligand,
                            ignore_inter,
                            B_include_1_4_interactions, scale_1_4,
                            unbound_internal_FE, info, B_use_non_bond_cutoff, B_have_flexible_residues);
                            //parameterArray, unbound_internal_FE, info, B_use_non_bond_cutoff, B_have_flexible_residues);
            // end of Step 3 // }

            // Step 4 // {
            //
            // Compute the energy of the unbound extended state
            // ------------------------------------------------
            //
            // Convert from unbound state to unbound coordinates
            cnv_state_to_coords( sUnbound_ext, vt, tlist, sUnbound_ext.ntor, crdpdb, crd, natom );
            //
            // Calculate the unbound internal energy using the standard AutoDock energy function
            (void) eintcalPrint(nonbondlist, ad_energy_tables, crd, Nnb, B_calcIntElec, B_include_1_4_interactions, scale_1_4, qsp_abs_charge, B_use_non_bond_cutoff, B_have_flexible_residues);
            //
            // eintcal() and eintcalPrint() set the values of the hideously-global nb_group_energy[]
            unbound_ext_internal_FE = nb_group_energy[INTRA_LIGAND] + nb_group_energy[INTRA_RECEPTOR];
            //
            pr(logFile, "\n\nThe internal energy of the unbound extended state was computed to be %+.3lf kcal/mol\n\n", unbound_ext_internal_FE);
            // end of Step 4 // }

            // Step 5 // {
            //
            // Decide whether to use extended or AutoDock state for unbound state
            // ------------------------------------------------------------------
            //
            if (unbound_ext_internal_FE > 0.0) {
                // Unbound extended state has an internal energy that is positive

                // Step 5.1 // {
                //
                // Repeat Step 1 with the standard AutoDock internal energy potentials
                //
                // Run a hybrid global-local search using the autodock energy tables
                // -----------------------------------------------------------------
                //
                //  *  Turn off the use of the non-bond cutoff
                //  *  Turn off internal electrostatics
                //  *  Turn off intermolecular energy calculations
                // TODO Do not translate or rotate the ligand in unbound searches
                //
                /*
                 *  Genetic Algorithm-Local search,  a.k.a.
                 *  Lamarckian Genetic Algorithm
                 */
                (void) fprintf( logFile, "\n\n\tBEGINNING COMPUTATION OF UNBOUND AUTODOCK STATE USING LGA\n");
                (void) fprintf( logFile,     "\t_________________________________________________________\n\n\n");
                (void) fflush( logFile );
                //
                pr(logFile, "Date:\t");
                printdate( logFile, 2 );
                (void) fflush( logFile );
                gaStart = times( &tms_gaStart );
                //  Can get rid of the following line
                ((Genetic_Algorithm *)GlobalSearchMethod)->initialize(pop_size, 7+sInit.ntor);
                //
                // Start Lamarckian GA run searching only torsions -- Unbound simulation
                // sUnbound_ad = call_glss_tors( GlobalSearchMethod, LocalSearchMethod,
                sUnbound_ad = call_glss( GlobalSearchMethod, LocalSearchMethod,
                                         sInit,
                                         num_evals_unbound, pop_size,
                                         outlev,
                                         outputEveryNgens, &ligand,
                                         B_RandomTran0, B_RandomQuat0, B_RandomDihe0,
                                         info, FN_pop_file, end_of_branch );
                // State of best individual at end of GA-LS run, sUnbound_ad, is returned.
                // Finished Lamarckian GA run
                gaEnd = times( &tms_gaEnd );
                pr( logFile, "\nFinished Lamarckian Genetic Algorithm (LGA), time taken:\n");
                timesyshms( gaEnd - gaStart, &tms_gaStart, &tms_gaEnd );
                pr( logFile, "\n");
                printdate( logFile, 1 );
                (void) fflush( logFile );
                pr(logFile, "\nTotal number of Energy Evaluations: %lu\n", evaluate.evals() );
                pr(logFile, "Total number of Generations:        %u\n", ((Genetic_Algorithm *)GlobalSearchMethod)->num_generations());
                // Restore the number of torsions in the state variable "sUnbound_ad"
                sUnbound_ad.ntor = saved_sInit_ntor;
                // end of Step 5.1 // }

                // Step 5.2 // {
                //
                // Compute the energy of the unbound AutoDock state
                // ------------------------------------------------
                //
                // Convert from unbound state to unbound coordinates
                cnv_state_to_coords( sUnbound_ad, vt, tlist, sUnbound_ad.ntor, crdpdb, crd, natom );
                //
                // Calculate the unbound internal energy using the standard AutoDock energy function
                (void) eintcalPrint(nonbondlist, ad_energy_tables, crd, Nnb, B_calcIntElec, B_include_1_4_interactions, scale_1_4, qsp_abs_charge, B_use_non_bond_cutoff, B_have_flexible_residues);
                //
                // eintcal() and eintcalPrint() set the values of the hideously-global nb_group_energy[]
                unbound_ad_internal_FE = nb_group_energy[INTRA_LIGAND] + nb_group_energy[INTRA_RECEPTOR];
                //
                pr(logFile, "\n\nThe internal energy of the unbound AutoDock state was computed to be %+.3lf kcal/mol\n\n", unbound_ad_internal_FE);
                // end of Step 5.2 // }

                if (unbound_ad_internal_FE < unbound_ext_internal_FE) {
                    pr(logFile, "NOTE:   The AutoDock internal energy of the \"extended\" state was higher\nNOTE:   than that of the state obtained by searching using the AutoDock internal\nNOTE:   energy potentials.\nNOTE:   The unbound state was set to the AutoDock optimum state, not the \"extended\" state.\n\n");
                    unbound_internal_FE = unbound_ad_internal_FE;
                    sUnbound = sUnbound_ad;
                } else {
                    pr(logFile, "NOTE:   Although the AutoDock internal energy of the \"extended\" state was positive, it was lower\nNOTE:   than that of the state obtained by searching using the AutoDock internal\nNOTE:   energy potentials.\nNOTE:   The unbound state was set to the \"extended\" state.\n\n");
                    unbound_internal_FE = unbound_ext_internal_FE;
                    sUnbound = sUnbound_ext;
                }
            } else {
                // Unbound extended state has an internal energy that is negative
                unbound_internal_FE = unbound_ext_internal_FE;
                sUnbound = sUnbound_ext;
                pr(logFile, "NOTE:   The AutoDock internal energy of the \"extended\" state was negative.\n\nNOTE:   The unbound state was set to the \"extended\" state.\n\n");
            }
            //
            pr(logFile, "\n\nThe internal energy of the unbound state was set to %+.3lf kcal/mol\n\n", unbound_internal_FE);
            // end of Step 5 // }

            // Step 6 // {
            //
            // Convert from unbound state to unbound coordinates
            cnv_state_to_coords( sUnbound, vt, tlist, sUnbound.ntor, crdpdb, crd, natom );
            // end of Step 6 // }

            // Step 7 // {
            //
            // Output the coordinates of the unbound state
            pr( logFile, "\n\n\tFINAL UNBOUND STATE\n" );
            pr( logFile,     "\t___________________\n\n\n" );
            //
            writePDBQT( -1, seed,  FN_ligand, dock_param_fn, lig_center,
                        sUnbound, ntor, &eintra, &einter, natom, atomstuff,
                        crd, emap, elec,
                        charge, abs_charge, qsp_abs_charge,
                        ligand_is_inhibitor,
                        torsFreeEnergy,
                        vt, tlist, crdpdb, nonbondlist,
                        ad_energy_tables,
                        type, Nnb, B_calcIntElec,
                        map,
                        outlev,
                        ignore_inter,
                        B_include_1_4_interactions, scale_1_4, parameterArray, unbound_internal_FE,
                        info, UNBOUND, PDBQT_record, B_use_non_bond_cutoff, B_have_flexible_residues, ad4_unbound_model);
            // end of Step 7 // }

            // Step 8 // {
            //
            // Remember to reset the energy evaluator back to computing the intermolecular energy between
            // the flexible and the rigid molecules.
            evaluate.compute_intermol_energy(TRUE);
            // end of Step 8 // }

        } else {
            pr(logFile, "NOTE:  AutoDock cannot compute the energy of the unbound state, since the ligand is rigid.\n\n");
            pr(logFile, "NOTE:  Use the \"unbound energy\" command to set the energy of the unbound state, if known from a previous calculation where the ligand was treated as flexible.\n\n");
            unbound_internal_FE = 0.0L;
            ad4_unbound_model = User;
            pr(logFile, "\n\nThe internal energy of the unbound state was set to %+.3lf kcal/mol\n\n", unbound_internal_FE);
        }

        pr( logFile, UnderLine );
        (void) fflush(logFile);
        break;

/*____________________________________________________________________________*/

    case DPF_EPDB:
        /*
         *  epdb
         *
         *  Computes the energy of the ligand specified by the "move lig.pdbqt" command.
         *  Return the energy of the Small Molecule.
         *  FN_ligand must be in   PDBQT-format;
         *  flag can be:-
         *  0 = NEW, or   PDBQT-71, and
         *  1 = OLD, or   PDBQT-55 (old PDBq format).
         */
        outside = FALSE;
        atoms_outside = FALSE;
        eintra = 0.0L;
        einter = 0.0L;
        etotal = 0.0L;

        pr(logFile, "WARNING This command, \"epdb\", currently computes the energy of the ligand specified by the \"move lig.pdbqt\" command.\n");
        retval = sscanf(line, "%*s %s", dummy_FN_ligand);
        if (retval >= 1) {
            pr(logFile, "WARNING  -- it will not read in the PDBQT file specified on the \"epdb\" command line.\n");
        }

        /*
        (void) sscanf(line, "%*s %s", FN_ligand);
        pr(logFile, "epdb %s\n\n", FN_ligand);
        natom = 0;
        print_1_4_message(logFile, B_include_1_4_interactions, scale_1_4);
        //
        ligand = readPDBQT( line,
                            num_atom_types,
                            &natom,
                            crdpdb, crdreo, charge, &B_haveCharges,
                            type, bond_index,
                            pdbaname, FN_ligand, FN_flexres, B_have_flexible_residues, atomstuff, &n_heavy_atoms_in_ligand,
                            &B_constrain_dist, &atomC1, &atomC2,
                            &sqlower, &squpper,
                            &ntor1, &ntor, &ntor_ligand,
                            tlist, vt,
                            &Nnb, nonbondlist,
                            jobStart, tms_jobStart, hostnm, &ntorsdof, outlev,
                            ignore_inter,
                            B_include_1_4_interactions,
                            atoms, PDBQT_record, end_of_branch );

        //
        // pre-calculate some values we will need later in computing the desolvation energy
        //
        for (i=0; i<natom; i++) {
            abs_charge[i] = fabs(charge[i]);
            qsp_abs_charge[i] = qsolpar * abs_charge[i];
        }
        */
        exit_if_missing_elecmap_desolvmap_about("epdb");

        // to restore the original coordinates, we must temporarily undo the centering transformation
        for ( i=0; i<true_ligand_atoms; i++ ) {
            for (xyz = 0;  xyz < SPACE;  xyz++) {
                crdpdb[i][xyz] += lig_center[xyz];
            }
        }
        // determine if any atoms are outside the grid box
        atoms_outside = FALSE;
        for (i=0; i<natom; i++) {
            outside = is_out_grid_info(crdpdb[i][X], crdpdb[i][Y], crdpdb[i][Z]);
            if (outside) {
                atoms_outside = TRUE;
                (void) sprintf( message, "%s: WARNING: Atom %d (%.3f, %.3f, %.3f) is outside the grid!\n", programname, i+1, crdpdb[i][X], crdpdb[i][Y], crdpdb[i][Z] );
                print_2x( logFile, stderr, message );
                outside = FALSE; /* Reset outside */
            }
        }
        pr(logFile, "Number of \"true\" ligand atoms:  %d\n", true_ligand_atoms);
        //
        for (i=0;i<natom;i++) {
            if (ignore_inter[i] == 1) {
                pr(logFile, "Special Boundary Conditions:\n");
                pr(logFile, "____________________________\n\n");
                pr(logFile, "AutoDock will ignore the following atoms in the input PDBQT file \nin intermolecular energy calculations:\n");
                pr(logFile, "\n(This is because these residue atoms are at the boundary between \nflexible and rigid, and since they cannot move \nthey will not affect the total energy.)\n\n");
                break;
            }
        }
        for (i=0;i<natom;i++) {
            if (ignore_inter[i] == 1) {
                pr(logFile, "Atom number %d:  %s\n", i+1, atomstuff[i] );
            }
        }
        pr(logFile, "\n");

        sInit.ntor = ligand.S.ntor;

        // save any currently-computed unbound internal FE
        //unbound_internal_FE_saved = unbound_internal_FE;
        //ad4_unbound_model_saved = ad4_unbound_model;
        //ad4_unbound_model = User;

        // Initialise to zero, since we may call "epdb" more than once in a single DPF
        // Set the unbound free energy -- assume it is zero, since this is a "single-point energy calculation"
        unbound_internal_FE = 0.0;

        // Calculate the internal energy
        if (ntor > 0) {
            (void) eintcalPrint(nonbondlist, ad_energy_tables, crdpdb, Nnb, B_calcIntElec, B_include_1_4_interactions, scale_1_4, qsp_abs_charge, B_use_non_bond_cutoff, B_have_flexible_residues);
        }

        pr(logFile, "Unbound model to be used is %s.\n", report_parameter_library());
        // calculate the energy breakdown for the input coordinates, "crdpdb"
        eb = calculateBindingEnergies( natom, ntor, unbound_internal_FE, torsFreeEnergy, B_have_flexible_residues,
                                crdpdb, charge, abs_charge, type, map, info, outside,
                                ignore_inter, elec, emap, &elec_total, &emap_total,
                                nonbondlist, ad_energy_tables, Nnb, B_calcIntElec,
                                B_include_1_4_interactions, scale_1_4, qsp_abs_charge, 
                                B_use_non_bond_cutoff, ad4_unbound_model);

        pr(logFile, "\n\n\t\tIntermolecular Energy Analysis\n");
        pr(logFile,     "\t\t==============================\n\n");
        pr(logFile, "Atom  vdW+Hb+Elec  vdW+Hbond  Electrosta  Partial          Coordinates         \n");
        pr(logFile, "Type    Energy      Energy    tic Energy  Charge      x         y         z    \n");
        pr(logFile, "____  __________  __________  __________  _______  ________  ________  ________\n");
        //          "1234  0123456789  0123456789  0123456789  1234567  12345678  12345678  12345678"

        charge_total = 0.;
        etot = 0.;
        for (i = 0;  i < natom;  i++) {
            etot = emap[i] + elec[i];
            pr(logFile, "%4d  %10.2f  %10.2f  %10.2f  %7.3f  %8.4f  %8.4f  %8.4f\n", (type[i]+1), etot, emap[i], elec[i], charge[i], crdpdb[i][X], crdpdb[i][Y], crdpdb[i][Z]);
            charge_total += charge[i];
        } /*i*/
        pr(logFile, "      __________  __________  __________  _______\n");
        pr(logFile, "Total %10.2lf  %10.2lf  %10.2lf  %7.3lf\n",        (double)(emap_total + elec_total), (double)emap_total, (double)elec_total, (double)charge_total);
        pr(logFile, "      __________  __________  __________  _______\n");
        pr(logFile, "      vdW+Hb+Elec  vdW+Hbond  Electrosta  Partial\n");
        pr(logFile, "        Energy      Energy    tic Energy  Charge\n\n");

        pr(logFile, "Total Intermolecular Interaction Energy   = %+.3lf kcal/mol\n", (double)(eb.e_inter + eb.e_intra) );
        pr(logFile, "Total Intermolecular vdW + Hbond Energy   = %+.3lf kcal/mol\n", (double)emap_total);
        pr(logFile, "Total Intermolecular Electrostatic Energy = %+.3lf kcal/mol\n\n\n", (double)elec_total);

        printEnergies( &eb, "epdb: USER    ", ligand_is_inhibitor, emap_total, elec_total, B_have_flexible_residues, ad4_unbound_model );
        pr(logFile, "\n");

        // remember to re-center the ligand
        for ( i=0; i<true_ligand_atoms; i++ ) {
            for (xyz = 0;  xyz < SPACE;  xyz++) {
                crdpdb[i][xyz] -= lig_center[xyz];
            }
        }

        // restore the saved unbound internal FE
        //unbound_internal_FE = unbound_internal_FE_saved;
        //ad4_unbound_model = ad4_unbound_model_saved;

        (void) fflush(logFile);
        break;

/*____________________________________________________________________________*/

    case DPF_TERMINATION:
        /*
         *  ga_termination energy 0.1
         *  ga_termination evals 25000  // the best energy did not change in this time
         *  ga_termination time 120 s
         */
        /*
        (void) sscanf( line, "%*s %d", &i );
        */
        (void) fflush(logFile);
        break;

/*____________________________________________________________________________*/

    case GA_CROSSOVER_MODE:
        /*
         * ga_crossover_mode OnePt
         * ga_crossover_mode TwoPt
         * ga_crossover_mode Uniform
         * ga_crossover_mode Arithmetic
         *
         * Xover_Mode c_mode = OnePt;  //  can be: OnePt, TwoPt, Uniform or Arithmetic
         */
        (void) sscanf( line, "%*s %s", c_mode_str );
        if (strcmp(c_mode_str, "onept") == 0) {
            c_mode = OnePt;
            pr(logFile, "One-point crossover will be used in GA and LGA searches.\n");
        } else if (strcmp(c_mode_str, "twopt") == 0) {
            c_mode = TwoPt;
            pr(logFile, "Two-point crossover will be used in GA and LGA searches.\n");
        } else if (strcmp(c_mode_str, "uniform") == 0) {
            c_mode = Uniform;
            pr(logFile, "Uniform crossover will be used in GA and LGA searches.\n");
        } else if (strcmp(c_mode_str, "arithmetic") == 0) {
            c_mode = Arithmetic;
            pr(logFile, "Arithmetic crossover will be used in GA and LGA searches.\n");
        } else if (strcmp(c_mode_str, "branch") == 0) {
            c_mode = Branch;
            pr(logFile, "Branch crossover will be used in GA and LGA searches.\n");
        } else {
            c_mode = TwoPt; // default
        }
        (void) fflush(logFile);
        break;

/*____________________________________________________________________________*/

    case GA_TOURNAMENT_SELECTION:
        /*
         * ga_tournament_selection
         */
        
	// M Pique - does not appear to work so disabling for now (Oct 2009)
        //s_mode = Tournament;
        //pr(logFile, "Tournament selection will be used in GA and LGA searches.\n");
        prStr( error_message, "%s:  ERROR! Tournament selection is not yet implemented! Instead proportional selection will be used in GA and LGA searches\n\n", programname);
        pr_2x( logFile, stderr, error_message );
        s_mode = Proportional;
        (void) fflush(logFile);
        break;

/*____________________________________________________________________________*/

    case GA_LINEAR_RANKING_SELECTION:
        /*
         * ga_linear_ranking_selection [probability_ratio Real]
         */
        
        (void) sscanf( line, "%*s %s", c_mode_str );
        if (strcmp(c_mode_str, "probability_ratio") == 0){ 
            (void) sscanf( line, "%*s %*s %f", &linear_ranking_selection_probability_ratio);
        }

        s_mode = LinearRanking;
        pr(logFile, "Linear ranking selection will be used in GA and LGA searches.\n");
        (void) fflush(logFile);
        break;

/*____________________________________________________________________________*/

    case GA_PROPORTIONAL_SELECTION:
        /*
         * ga_proportional_selection 
         */
        
        s_mode = Proportional;
        pr(logFile, "Proportional selection will be used in GA and LGA searches.\n");
        (void) fflush(logFile);
        break;

/*____________________________________________________________________________*/

    case GA_BOLTZMAN_SELECTION:
        /*
         * ga_boltzman_selection 
         */
        
        prStr( error_message, "%s:  ERROR! Boltzman selection is not yet implemented! Instead proportional selection will be used in GA and LGA searches\n\n", programname);
        pr_2x( logFile, stderr, error_message );
        s_mode = Proportional;
        (void) fflush(logFile);
        break;

/*____________________________________________________________________________*/


    case DPF_POPFILE:
        /*
         *  output_pop_file
         *
         *  Used to write out the population to a file at the end of
         *  every GA.
         */
        (void) sscanf( line, "%*s %s", FN_pop_file);
        (void) fflush(logFile);
        pr( logFile, "The population will be written to the file \"%s\" at the end of every generation.\n", FN_pop_file);
        break;

 /*____________________________________________________________________________*/

     case DPF_CONFSAMPLER:
        /*
         * confsampler
         *
         * Scan a region around conformations saved in sHist array.
         *
         */

        (void) sscanf( line, "%*s %s %d", confsampler_type, &confsampler_samples);
        pr( logFile, "Scanning local regions around each docked conformation.\n");

        exit_if_missing_elecmap_desolvmap_about("confsampler");

        if (strcmp(confsampler_type, "systematic") == 0) {
            systematic_conformation_sampler(sHist, nconf, vt, crdpdb, tlist, lig_center, natom, type, info);
        }

        else {
            random_conformation_sampler(sHist, nconf, confsampler_samples, vt, crdpdb, tlist, lig_center, natom, type, info);
        }
        (void) fflush(logFile);

        break;

/*____________________________________________________________________________*/

    case DPF_COPYRIGHT:
        /*
         * 'copyright' to show the Gnu GPL copyright
         */
        show_copyright(logFile);
        (void) fflush(logFile);
        break;

/*____________________________________________________________________________*/

    case DPF_WARRANTY:
        /*
         * 'warranty' to show the Gnu GPL warranty
         */
        show_warranty(logFile);
        (void) fflush(logFile);
        break;

/*_12yy_______________________________________________________________________*/

    case DPF_:
        /*
         *
         */
        /*
        (void) sscanf( line, "%*s %d", &i );
        (void) fflush(logFile);
        break;
        */

//______________________________________________________________________________

    default:
        /*
        **  Do nothing...
        */
        break;

//______________________________________________________________________________

    } /* switch( dpf_keyword ) */

} /* while PARSING-DPF parFile */

/* __________________________________________________________________________
**
** Close the docking parameter file...
** __________________________________________________________________________
*/
pr( logFile, ">>> Closing the docking parameter file (DPF)...\n\n" );
//pr( logFile, UnderLine );
(void) fclose( parFile );


//______________________________________________________________________________
/*
** Print the time and date when the docking has finished...
*/

pr( logFile, "This docking finished at:\t\t\t" );
printdate( logFile, 1 );
pr( logFile, "\n\n\n" );

success( hostnm, jobStart, tms_jobStart );

 if(write_stateFile){
   fprintf(stateFile,"</autodock>\n");
   (void) fclose( stateFile );
 }
(void) fclose( logFile );

// delete arrays
delete []nonbondlist;


//________________________________________________________________________________
/*
** End of Boinc
*/
#ifdef BOINCCOMPOUND
 boinc_fraction_done(1.);
#endif
#ifdef BOINC
    boinc_finish(0);       /* should not return */
#endif

return 0;

} /* END OF PROGRAM */

/* AutoDock main private utility functions
*/
static void exit_if_missing_elecmap_desolvmap_about(string keyword)
{

    char error_message[LINE_LEN];

    if (!B_found_elecmap) {
         prStr(error_message, "%s:  %s command: no \"elecmap\" command has been specified!\n", programname, keyword.c_str());
         stop(error_message);
         exit(-1);
    } else if (!B_found_desolvmap) {
         prStr(error_message, "%s:  %s command: no \"desolvmap\" command has been specified!\n", programname, keyword.c_str());
         stop(error_message);
         exit(-1);
    } else if (!B_found_about_keyword){
         prStr(error_message, "%s:  %s command: no \"about\" command has been specified!\n", programname, keyword.c_str());
         stop(error_message);
         exit(-1);
    }
}


#ifdef BOINC
/*  Dummy graphics API entry points.
 *  This app does not do graphics, but it still must provide these callbacks.
 */

void app_graphics_render(int xs, int ys, double time_of_day) {}
void app_graphics_reread_prefs(){}
void boinc_app_mouse_move(int x, int y, bool left, bool middle, bool right ){}
void boinc_app_mouse_button(int x, int y, int which, bool is_down){}
void boinc_app_key_press(int wParam, int lParam){}
void boinc_app_key_release(int wParam, int lParam){}
#endif

/* EOF */
