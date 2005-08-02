/*

 $Id: main.cc,v 1.18 2005/08/02 23:02:54 garrett Exp $

*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/* main.cc */

// possibly unnecessary // #include <iostream.h>
#include <math.h>

/* the BOINC API header file */
#ifdef BOINC
#include "diagnostics.h"
#include "boinc_api.h" 
#include "filesys.h" 		// boinc_fopen(), etc... */
#endif

#include "coliny.h"
#include "hybrids.h"
#include "ranlib.h"
#include "gs.h"
#include "ls.h"
#include "rep.h"
#include "support.h"

#include <sys/types.h> // time_t time(time_t *tloc);
#include <time.h>      // time_t time(time_t *tloc);
#include <sys/times.h>
#include <stdlib.h>
#include <string.h>
#include <sys/param.h>
#include <ctype.h> // tolower
#include <unistd.h> // sysconf

#include "main.h"

#ifdef sun
    extern "C"
    {
        /* Needed on Sun */
        int gethostname(char *name,int namelen);
    }
#endif

extern int debug;
extern int keepresnum;
extern FloatOrDouble idct;
extern Eval evaluate;
extern Linear_FE_Model AD4;
extern int oldpdbq;


int sel_prop_count = 0;


int main (int argc, char * const argv[], char * const envp[])

/*******************************************************************************
**      Name: main  (AutoDock)                                                **
**  Function: Performs Automated Docking of Small Molecule into Macromolecule **
** Copyright: (C) 1994-2005 TSRI, Arthur J. Olson's Labortatory.              **
**____________________________________________________________________________**
**   Authors: Garrett Matthew Morris, Current C/C++ version 4.0               **
**                                       e-mail: garrett@scripps.edu          **
**                                                                            **
**            David Goodsell, Orignal FORTRAN version 1.0                     **
**                                       e-mail: goodsell@scripps.edu         **
**                                                                            **
**            The Scripps Research Institute                                  **
**            Department of Molecular Biology, MB5                            **
**            10550 North Torrey Pines Road                                   **
**            La Jolla, CA 92037.                                             **
**                                                                            **
**      Date: 02/10/04                                                        **
**____________________________________________________________________________**
**    Inputs: Control file, Small Molecule PDBQT file, macromolecular grid map**
**            files.                                                          **
**   Returns: Autodock Log File, includes docked conformation clusters (PDBQT)**
**   Globals: X, Y, Z, SPACE, MAX_GRID_PTS, NEINT, MAX_ATOMS,                 **
**            MAX_MAPS, ATOM_MAPS, A_DIVISOR, LINE_LEN,TRUE,FALSE             **
**            programname, command_mode,                                      **
**            command_in_fp, command_out_fp, parFile, logFile.                **
**____________________________________________________________________________**
** Modification Record                                                        **
** Date     Inits   Comments                                                  **
** 09/06/95 RSH     Added code to handle GA/LS stuff                          **
**          DSG     Quaternion rotations                                      **
**          DSG     Generates torsions from annotated pdb list                **
**          DSG     Generates internal energies                               **
**          DSG     Performs a limited Cluster analysis of conformations      **
** 05/07/92 GMM     C translation                                             **
** 05/14/92 GMM     Time-dependent seed in random-number generation           **
** 10/29/92 GMM     Application Visualization System (AVS) readable grid      **
**                  display file input.                                       **
**                  [AVS is a trademark of Stardent Computer Inc.]            **
** 04/17/93 GMM     '-o' oldpdbq flag for q in columns 55-61;                 **
**                  the default PDBQ format is now columns 71-76;             **
**                  Also added the 'total_charge' check.                      **
** 11/19/93 GMM     #ifdef NOSQRT, with non-square-rooting acceleration.      **
** 09/26/94 GMM     Cluster analysis now outputs RMS deviations.              **
** 09/28/94 GMM     Modularized code.                                         **
** 10/02/94 GMM     Distance constraints added, for Ed Moret. Accelerated.    **
** 09/06/95 RSH     Incorporation of GA/SW tokens                             **
*******************************************************************************/

{
char atm_typ_str[ATOM_MAPS]; //  "atm_typ_str" used to serve a similar role to "ligand_atom_types" in AD3.
char *ligand_atom_type_ptrs[MAX_MAPS]; /* array of ptrs used to parse input line of atom type names */
char atomstuff[MAX_ATOMS][MAX_CHARS];
char error_message[LINE_LEN];
char message[LINE_LEN];
char FN_clus[MAX_CHARS];
char FN_watch[MAX_CHARS];
char FN_gdfld[MAX_CHARS];
char FN_gpf[MAX_CHARS];
char FN_ligand[MAX_CHARS];
char dummy_FN_ligand[MAX_CHARS];
char FN_template_energy_file[MAX_CHARS];
char FN_trj[MAX_CHARS];
char FN_receptor[MAX_CHARS];
char FN_rms_ref_crds[MAX_CHARS];
char FN_parameter_library[MAX_CHARS];
char hostnm[MAX_CHARS];
char line[LINE_LEN];
char line_template[LINE_LEN];
char out_acc_rej = '?';
char param[2][MAX_CHARS];
char timeSeedIsSet[2];
char pdbaname[MAX_ATOMS][5];
char selminpar = 'm';
char S_contype[8];
char torfmt[LINE_LEN];
char ligand_atom_types[MAX_MAPS][3];

static ParameterEntry * foundParameter;
ParameterEntry parameterArray[MAX_MAPS];

FILE *template_energy_file;

FloatOrDouble crdpdb[MAX_ATOMS][SPACE];
FloatOrDouble crd[MAX_ATOMS][SPACE];
FloatOrDouble lig_center[SPACE];
FloatOrDouble map_center[SPACE];
FloatOrDouble vt[MAX_TORS][SPACE];
FloatOrDouble cA;
FloatOrDouble cB;
FloatOrDouble charge[MAX_ATOMS];
FloatOrDouble abs_charge[MAX_ATOMS];
FloatOrDouble qsp_abs_charge[MAX_ATOMS];
FloatOrDouble econf[MAX_RUNS];
FloatOrDouble elec[MAX_ATOMS];
FloatOrDouble emap[MAX_ATOMS];
FloatOrDouble epsij;
FloatOrDouble F_A;
FloatOrDouble F_Aova;
FloatOrDouble F_tor;
FloatOrDouble F_TorConRange[MAX_TORS][MAX_TOR_CON][2];
FloatOrDouble F_torPref;
FloatOrDouble F_torHWdth;
FloatOrDouble mapmax[MAX_MAPS];
FloatOrDouble mapmin[MAX_MAPS];
FloatOrDouble q1q2[MAX_NONBONDS];
FloatOrDouble Rij;
FloatOrDouble sqlower;
FloatOrDouble squpper;
FloatOrDouble template_energy[MAX_ATOMS]; // template energy value for each atom
FloatOrDouble template_stddev[MAX_ATOMS]; // and standard deviation of this energy
FloatOrDouble tmpconst;
FloatOrDouble xhi;
FloatOrDouble xlo;
FloatOrDouble yhi;
FloatOrDouble ylo;
FloatOrDouble zhi;
FloatOrDouble zlo;

// Distance-dependence in Desolvation Term
const double sigma = 3.6L;
const double qsolpar = 0.01097L;

// i
double Ri, epsi, Ri_hb, epsi_hb;
hbond_type hbondi;
// j
double Rj, epsj, Rj_hb, epsj_hb;
hbond_type hbondj;

#ifdef USE_DOUBLE
FloatOrDouble scale_1_4 = 0.5L;
FloatOrDouble c=0.0L;
FloatOrDouble clus_rms_tol = 0.0L;
FloatOrDouble e0max = BIG;
FloatOrDouble eintra = 0.0L;
FloatOrDouble einter = 0.0L;
FloatOrDouble etotal = 0.0L;
FloatOrDouble AD3_FE_coeff_estat   = 1.000L;
FloatOrDouble inv_spacing = 0.0L;
FloatOrDouble qtwFac = 1.0L;
FloatOrDouble qtwStep0 = 5.0L;
FloatOrDouble qtwStepFinal = 5.0L;
FloatOrDouble maxrad = -1.0L;
FloatOrDouble r2sum=0.0L;
FloatOrDouble RJ = 8.31441L;     // in J/K/mol, Gas Constant, Atkins Phys.Chem., 2/e
FloatOrDouble Rcal = 1.9871917L; // in cal/K/mol, Gas Constant, RJ/4.184
FloatOrDouble T0K = 273.15L;     // 0 degrees Celsius, in K
FloatOrDouble RTreduc = 1.0L;
FloatOrDouble spacing = 0.0L;
FloatOrDouble RT0 = 616.0L;
FloatOrDouble RTFac = 0.95L;
FloatOrDouble torsdoffac = 0.3113L;
FloatOrDouble torsFreeEnergy = 0.0L;
FloatOrDouble torFac = 1.0L;
FloatOrDouble torStep0 = 5.0L;
FloatOrDouble torStepFinal = 5.0L;
FloatOrDouble trnFac = 1.0L;
FloatOrDouble trnStep0 = 0.2L;
FloatOrDouble trnStepFinal = 0.2L;
FloatOrDouble WallEnergy = 1.0e8L; /* Energy barrier beyond walls of gridmaps. */
//  The GA Stuff
FloatOrDouble m_rate = 0.02L;
FloatOrDouble c_rate = 0.80L;
FloatOrDouble alpha = 0.0L;
FloatOrDouble beta = 1.0L;
FloatOrDouble search_freq = 0.06L;
FloatOrDouble rho = 1.0L;
FloatOrDouble lb_rho = 0.01L;
FloatOrDouble *rho_ptr = NULL;
FloatOrDouble *lb_rho_ptr = NULL;
FloatOrDouble unbound_internal_FE = 0.0L;
FloatOrDouble emap_total = 0.L;
FloatOrDouble elec_total = 0.L;
FloatOrDouble charge_total = 0.L;
FloatOrDouble etot = 0.L;

#else

FloatOrDouble scale_1_4 = 0.5;
FloatOrDouble c=0.0;
FloatOrDouble clus_rms_tol = 0.0;
FloatOrDouble e0max = BIG;
FloatOrDouble eintra = 0.0;
FloatOrDouble einter = 0.0;
FloatOrDouble etotal = 0.0;
FloatOrDouble AD3_FE_coeff_estat   = 1.000;
FloatOrDouble inv_spacing = 0.0;
FloatOrDouble qtwFac = 1.0;
FloatOrDouble qtwStep0 = 5.0;
FloatOrDouble qtwStepFinal = 5.0;
FloatOrDouble maxrad = -1.0;
FloatOrDouble r2sum=0.0;
FloatOrDouble RJ = 8.31441;     // in J/K/mol, Gas Constant, Atkins Phys.Chem., 2/e
FloatOrDouble Rcal = 1.9871917; // in cal/K/mol, Gas Constant, RJ/4.184
FloatOrDouble T0K = 273.15;        // 0 degrees Celsius, in K
FloatOrDouble RTreduc = 1.0;
FloatOrDouble spacing = 0.0;
FloatOrDouble RT0 = 616.0;
FloatOrDouble RTFac = 0.95;
FloatOrDouble torsdoffac = 0.3113;
FloatOrDouble torsFreeEnergy = 0.0;
FloatOrDouble torFac = 1.0;
FloatOrDouble torStep0 = 5.0;
FloatOrDouble torStepFinal = 5.0;
FloatOrDouble trnFac = 1.0;
FloatOrDouble trnStep0 = 0.2;
FloatOrDouble trnStepFinal = 0.2;
FloatOrDouble WallEnergy = 1.0e8; /* Energy barrier beyond walls of gridmaps. */
//  The GA Stuff
FloatOrDouble m_rate = 0.02;
FloatOrDouble c_rate = 0.80;
FloatOrDouble alpha = 0;
FloatOrDouble beta = 1;
FloatOrDouble search_freq = 0.06;
FloatOrDouble rho = 1.0;
FloatOrDouble lb_rho = 0.01;
FloatOrDouble *rho_ptr = NULL;
FloatOrDouble *lb_rho_ptr = NULL;
FloatOrDouble unbound_internal_FE = 0.0;
FloatOrDouble emap_total = 0.;
FloatOrDouble elec_total = 0.;
FloatOrDouble charge_total = 0.;
FloatOrDouble etot = 0.;

#endif

unsigned short US_TorE[MAX_TORS];
unsigned int outputEveryNgens = 100;

Boole B_atom_types_found = FALSE;
Boole B_isGaussTorCon = FALSE;
Boole B_constrain_dist;
Boole B_either = FALSE;
Boole B_calcIntElec = FALSE;
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
Boole B_template = FALSE;
Boole B_torReduc = FALSE;
Boole B_trnReduc = FALSE;
Boole B_write_all_clusmem = FALSE;
Boole B_isTorConstrained[MAX_TORS];
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

int atm1=0;
int atm2=0;
int a1=0;
int a2=0;
int atomC1;
int atomC2;
int curatm=0;
int dpf_keyword = -1;
int dummy_oldpdbq;
int gridpts1[SPACE];
int gridpts[SPACE];
int Htype = 0;
int ncycles = -1;
int iCon=0;
int imap=0;
int indcom = 0;
int ligand_is_inhibitor = 1;
int ltorfmt = 4;
int nruns = 0;
int nstepmax = -1;
int naccmax = 0;
int natom = 0;
int nonbondlist[MAX_NONBONDS][MAX_NBDATA];
int nconf = 0;
int ncycm1 = 1;
int ndihed = 0;
int nlig = 0;
int nres = 0;
int nmol = 0;
int Nnb = 0;
int nrejmax = 0;
int ntor;
int ntor1;
int ntorsdof = 0;
int num_all_maps = 0;
int num_atom_types = 0;
int nval = 0;
int outlev = -1;
int retval = 0;
int status = 0;
int tlist[MAX_TORS][MAX_ATOMS];
int trj_end_cyc = 0;
int trj_begin_cyc = 0;
int trj_freq = 0;
int type[MAX_ATOMS];
int bond_index[MAX_ATOMS];
int xA;
int xB;
int I_tor;
int I_torBarrier;
int N_con[MAX_TORS];
int MaxRetries = 1000; /* Default maximum number of retries for ligand init. */
int OutputEveryNTests = 1000;
int NumLocalTests = 10;
int maxTests = 10000;
int ignore_inter[MAX_ATOMS];
int par_lib_FN_found = 0;
/* int beg; */
/* int end; */
/* int imol = 0; */
int outside = FALSE;
int atoms_outside = FALSE;

unsigned short US_energy;
unsigned short US_tD;
unsigned short US_torBarrier = TORBARMAX;
unsigned short US_torProfile[MAX_TORS][NTORDIVS];
unsigned short US_min = TORBARMAX;

register int i = 0;
register int j = 0;
int j1 = 1;
register int k = 0;
register int xyz = 0;


State sInit;                             /* FloatOrDouble qtn0[QUAT], tor0[MAX_TORS]; */
State sHist[MAX_RUNS];  /*qtnHist[MAX_RUNS][QUAT],torHist[MAX_RUNS][MAX_TORS];*/

Molecule mol;        /* ligand */

static FloatOrDouble e_internal[NEINT][ATOM_MAPS][ATOM_MAPS];
static FloatOrDouble sol_fn[NEINT];  // internal desolvation energy function
static FloatOrDouble F_A_from;
static FloatOrDouble F_A_to;
static FloatOrDouble F_lnH;
static FloatOrDouble F_W;
static FloatOrDouble F_hW;
static FloatOrDouble map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS];
static FourByteLong clktck = 0;

static FloatOrDouble ELECSCALE = 83.015909;      /* relative dielectric epsilon = 4*r */
/* static FloatOrDouble ELECSCALE = 332.06363;   // relative dielectric epsilon = r   */
static FloatOrDouble version = 4.00;

struct tms tms_jobStart;
struct tms tms_gaStart;
struct tms tms_gaEnd;

Clock  jobStart;
Clock  gaStart;
Clock  gaEnd;

time_t time_seed;

//  The GA Stuff
FourByteLong seed[2];
unsigned int pop_size = 50;
unsigned int num_generations = 0;  //  Don't terminate on the basis of number of generations
unsigned int num_evals = 150000;
unsigned int max_its = 30;
unsigned int max_succ = 4;
unsigned int max_fail = 4;
int window_size = 10;
int low = 0;
int high = 100;
int elitism = 1;


Selection_Mode s_mode = Proportional;
Xover_Mode c_mode = TwoPt;
Worst_Mode w_mode = AverageOfN;
EvalMode e_mode = Normal_Eval;
Global_Search *GlobalSearchMethod = NULL;
Local_Search *LocalSearchMethod = NULL;

// GridMap grid_map;

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

if ( setflags(argc,argv) == -1) {
    exit(-1);
} /* END PROGRAM */

//______________________________________________________________________________
/*
** Initialize torsion arrays and constants.
*/

(void) strcpy( torfmt, "%*s\0" ); /* len(torfmt) is 4 chars */

for (j = 0;  j < MAX_ATOMS;  j++ ) {
    type[j] = 0;
    template_energy[j] = 0.0;
    template_stddev[j] = 1.0;
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
initialiseState( &(mol.S) );

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
        prStr( error_message, "ERROR: MAX_TORS = %d torsions declared in \"constants.h\";\n\t LINE_LEN = %d, Therefore you must change \"LINE_LEN\" to exceed %d...\n", MAX_TORS, LINE_LEN, 4+4*MAX_TORS );
        stop( error_message );
        exit( -1 );
    } else {
        (void) strcat( torfmt, " %lf\0" );  /* add on 4 chars  for each new torsion... */
    }
} /* len(torfmt) is 4+4*MAX_TORS chars */

for (j = 0; j < MAX_NONBONDS; j++) {
    nonbondlist[j][ATM1] = nonbondlist[j][ATM2] = 0;
}

for (j = 0; j < MAX_RUNS; j++) {
    // isort[j] = j;
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
        idct = (FloatOrDouble)1.0 / (FloatOrDouble)clktck;
        if (debug) {
            pr(logFile, "N.B. debug is on and set to %d\n\n", debug);
            pr(logFile, "\n\nFYI:  Number of clock ticks per second = %d\n", (int)clktck);
            pr(logFile, "FYI:  Elapsed time per clock tick = %.3e milli-seconds\n\n\n\n", idct * 1000. );
        }
    }
}

(void) strcpy(FN_rms_ref_crds,"unspecified filename\0");

//______________________________________________________________________________
//
// Read in default parameters

// read_parameter_library(FN_parameter_library, outlev);

//______________________________________________________________________________
/*
** log(x): compute the natural (base e) logarithm of x,
*/

F_lnH = ((FloatOrDouble)log(0.5));

//______________________________________________________________________________
/*
** Output banner...
*/

banner( version );

//______________________________________________________________________________
/*
** Print the time and date when the file was created...
*/

pr( logFile, "This file was created at:\t\t\t" );
printdate( logFile, 1 );

(void) strcpy(hostnm, "unknown host\0");

if (gethostname( hostnm, MAX_CHARS ) == 0) {
    pr( logFile, "                   using:\t\t\t\"%s\"\n", hostnm );
}

pr( logFile, "\nNOTE: \"rus\" stands for:\n\n      r = Real, wall-clock or elapsed time;\n      u = User or cpu-usage time;\n      s = System time\n\nAll timings are in seconds, unless otherwise stated.\n\n\n" );

//______________________________________________________________________________
/*
** (Note: "dock_param_fn" set in "setflags.c"...)
*/
pr( logFile, "Parameter file used for this docking:\t\t%s\n\n", dock_param_fn );


//______________________________________________________________________________
/*
** Start reading in the DPF parameter/run-control file,
*/

while( fgets(line, LINE_LEN, parFile) != NULL ) { /* PARSING-DPF parFile */
    // "line" is a string containing the current line of the input DPF.

    dpf_keyword = parse_dpf_line( line );

    switch( dpf_keyword ) {
        case -1:
            pr( logFile, "DPF> %s", line );
            prStr( error_message,"%s: WARNING: Unrecognized keyword in docking parameter file.\n", programname );
            pr_2x( logFile, stderr, error_message );
            continue;
            /* break; */

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
            pr( logFile, "Output Level = 0.  NO OUTPUT DURING DOCKING.\n" );
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
        pr( logFile, "\n\tOutput every %u generations.\n", outputEveryNgens );
        (void) fflush(logFile);
        break;

//______________________________________________________________________________

    case DPF_TYPES:
        /*
        ** types
        **       - for backwards compatibility with AutoDock 3 and earlier.
        ** ATOM_TYPE_NAMES
        */
        for (i=0; i<ATOM_MAPS-1; i++) {
            atm_typ_str[i] = '?'; 
        }
        atm_typ_str[ATOM_MAPS-1] = '\0';
        // This is deprecated 
        prStr( error_message, "%s: WARNING:  In AutoDock 4, the new command, \"ligand_types\" is preferable to the \"types\" command.\n\n", programname);
        pr_2x( logFile, stderr, error_message );
        dpftypes( &Htype, &num_all_maps, &num_atom_types, atm_typ_str, line );
        if (num_atom_types > 0) {
            B_atom_types_found = TRUE;
        } else {
            prStr( error_message, "%s:  ERROR!  No atom types have been found; we cannot continue without this information!\n\n", programname );
            pr_2x( logFile, stderr, error_message );
            exit(-1);
        }
        (void) fflush(logFile);
        break;

//______________________________________________________________________________

    case DPF_FLD:
        /*
        ** fld
        ** GRID_DATA_FILE
        ** Read the (AVS-format) grid data file, .fld
        */
        // TO DO: add outlev
        readfield( &inv_spacing, &spacing, FN_gdfld,
            FN_gpf, gridpts1, gridpts, &xhi,&yhi,&zhi,
            jobStart, line, &xlo,&ylo,&zlo, FN_receptor,
            map_center, tms_jobStart );
        (void) fflush(logFile);
        break;

//______________________________________________________________________________

    case DPF_MAP:
        /*
        ** map
        ** ATOMIC AFFINITY or ELECTROSTATIC GRID MAP
        ** Read in active site grid map...
        */
        B_charMap = FALSE;
        // "imap" is incremented each time we call "readmap".
        // map_index could be incremented here if we had the atom_type stored in each map...
        if (B_atom_types_found == TRUE) {
            readmap( &B_havemap, &imap, &num_atom_types, &spacing,
                     ligand_atom_types, FN_gdfld, gridpts1, gridpts,
                     jobStart, line, FN_receptor, map, map_center,
                     mapmax, mapmin, tms_jobStart, B_charMap, outlev);
                     //grid_map
                   //);
        } else {
            prStr( error_message, "%s:  ERROR! No atom types have been found; we cannot continue without this information!\n\n", programname );
            pr_2x( logFile, stderr, error_message );
            exit(-1);
        }
        (void) fflush(logFile);
        break;

//______________________________________________________________________________

    case DPF_CHARMAP:
        /*
        ** charmap
        ** ATOMIC AFFINITY or ELECTROSTATIC GRID MAP
        ** Read in active site grid map...
        */
        B_charMap = TRUE;
        if (B_atom_types_found == TRUE) {
            readmap( &B_havemap, &imap, &num_atom_types, &spacing,
                     ligand_atom_types, FN_gdfld, gridpts1, gridpts,
                     jobStart, line, FN_receptor, map, map_center,
                     mapmax, mapmin, tms_jobStart, B_charMap, outlev);
                     // grid_map
                   // );
        } else {
            prStr( error_message, "%s:  ERROR! No atom types have been found; we cannot continue without this information!\n\n", programname );
            pr_2x( logFile, stderr, error_message );
            exit(-1);
        }
        (void) fflush(logFile);
        break;

//______________________________________________________________________________

    case DPF_MOVE:
        /*
        ** move
        ** Movable ligand,
        */
        B_found_move_keyword = TRUE;

        print_1_4_message(logFile, B_include_1_4_interactions, scale_1_4);
 
        natom=0;
        mol = readPDBQT(line,
                        num_atom_types,
                        &natom,
                        crdpdb, charge, &B_haveCharges,
                        type, bond_index,
                        pdbaname, FN_ligand, atomstuff, Htype,
                        &B_constrain_dist, &atomC1, &atomC2,
                        &sqlower, &squpper,
                        &ntor1, &ntor, tlist, vt,
                        &Nnb, nonbondlist,
                        jobStart, tms_jobStart, hostnm, &ntorsdof, outlev,
                        ignore_inter,
                        B_include_1_4_interactions);
 
        // pre-calculate some values we will need later in computing the desolvation energy
        //
        for (i=0; i<natom; i++) {
            abs_charge[i] = fabs(charge[i]);
            qsp_abs_charge[i] = qsolpar * abs_charge[i];
        }
        pr(logFile, "Number of \"true\" ligand atoms:  %d\n", true_ligand_atoms);
 
        for (i=0;i<natom;i++) {
            if (ignore_inter[i] == 1) {
                pr(logFile, "Special Boundary Conditions:\n");
                pr(logFile, "----------------------------\n\n");
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
            pr(logFile,"Calculating the product of the partial atomic charges q1*q2 for all %d non-bonded pairs...\n\n\n",Nnb);
            if (outlev >= 0) {
                pr(logFile,"Non-bonded                           Scaled\n");
                pr(logFile,"   Pair     Atom1-Atom2    q1*q2      q1*q2\n");
                pr(logFile,"__________  ___________  _________  _________\n");
                for (i = 0;  i < Nnb;  i++) {
                    atm1 = nonbondlist[i][ATM1];
                    atm2 = nonbondlist[i][ATM2];
                    q1q2[i] = charge[atm1] * charge[atm2];
                    pr(logFile,"   %4d     %5d-%-5d    %5.2f",i+1,atm1+1,atm2+1,q1q2[i]);
                    q1q2[i] *= ELECSCALE * AD4.coeff_estat;
                    pr(logFile,"    %5.2f\n",q1q2[i]);
                }
                pr(logFile,"\n");
            }
        }

        sInit.ntor = mol.S.ntor;
        ++nmol;
        ++nlig;

        (void) fflush(logFile);
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

        char algname[64];
        char nruns_str[64];
        (void) sscanf(line, "%*s %s %d", algname, &nruns);
        (void) sscanf(line, "%*s %s %s", algname, nruns_str);

        if (strcmp(algname,"help")==0) {
            utilib::BasicArray<double> initvec;
            coliny_init(algname, "", initvec);
            prStr(error_message, "ERROR:  no optimizer type specified.");
            stop(error_message);
            exit(-1);
        }
        else if (strcmp(nruns_str,"help")==0) {
            utilib::BasicArray<double> initvec;
            coliny_init(algname, nruns_str, initvec);
            prStr(error_message, "ERROR:  no optimizer type specified.");
            stop(error_message);
            exit(-1);
        }

        if (!command_mode) {
            if (nruns>MAX_RUNS) {
                prStr(error_message, "ERROR:  %d runs requested, but only dimensioned for %d.\nChange \"MAX_RUNS\" in \"constants.h\".", nruns, MAX_RUNS);
                stop(error_message);
                exit(-1);
            }

            evaluate.setup(crd, charge, abs_charge, qsp_abs_charge, type, natom, map, inv_spacing,
                           elec, emap,
                           xlo, xhi, ylo, yhi, zlo, zhi, nonbondlist,
                           e_internal, Nnb, B_calcIntElec, q1q2, B_isGaussTorCon,
                           B_isTorConstrained, B_ShowTorE, US_TorE, US_torProfile,
                           vt, tlist, crdpdb, sInit, mol,
                           B_template, template_energy, template_stddev,
                           ignore_inter,
                           B_include_1_4_interactions, scale_1_4, sol_fn, parameterArray, unbound_internal_FE);

            char domain[1024];
            // NOTE: Coliny enforces the bound constraints, but since the
            // torsion angles are periodic, we simply prevent the optimizer
            // from going too far.
            if (sInit.ntor > 0) {
                sprintf(domain,"[%f,%f] [%f,%f] [%f,%f] [-1.0,1.1]^3 [-6.2832,12.5664] [-6.2832,12.5664]^%d",(double)xlo, (double)xhi, (double)ylo, (double)yhi, (double)zlo, (double)zhi, sInit.ntor);
            } else {
                sprintf(domain,"[%f,%f] [%f,%f] [%f,%f] [-1.0,1.1]^3 [-3.1416,3.1416]",(double)xlo, (double)xhi, (double)ylo, (double)yhi, (double)zlo, (double)zhi);
            }
            ucout << domain << endl;
            ucout << Flush;

            pr(logFile, "Number of Coliny %s dockings = %d run%c\n", algname, nruns, (nruns>1)?'s':' ');

            //
            // COLINY-SPECIFIC LOGIC - BEGIN
            //

            try {

                utilib::BasicArray<double> initvec, finalpt;
                coliny_init(algname, domain, initvec);

                for (j=0; j<nruns; j++) {

                  fprintf( logFile, "\n\n\tBEGINNING Coliny %s DOCKING\n",algname);
                  pr(logFile, "\nDoing %s run:  %d/%d.\n", algname, j+1, nruns);

                  if (timeSeedIsSet[0] == 'T') {
                      seed[0] = (FourByteLong)time( &time_seed );
                  }
                  if (timeSeedIsSet[1] == 'T') {
                      seed[1] = (FourByteLong)time( &time_seed );
                  }
                  setall(seed[0], seed[1]);
                  initgn(-1); // Reinitializes the state of the current generator

                  pr(logFile, "Seeds:  %ld %ld\n", seed[0], seed[1]);
                  (void) fflush(logFile);

                  colinyStart = times(&tms_colinyStart);

                  finalpt.resize( initvec.size() );
                  int neval, niters;
                  coliny_minimize( seed[0]+seed[1]*j+j, initvec, finalpt, neval, niters );
                  //fstr.flush();

                  make_state_from_rep( (double *)&finalpt, int(finalpt.size()), &sHist[nconf]);

                  pr(logFile, "\nFinal docked state:\n");
                  pr(logFile, "\nTotal Num Evals: %d\n", neval);
                  printState(logFile, sHist[nconf], 2);

                  colinyEnd = times(&tms_colinyEnd);
                  pr(logFile, "Time taken for this %s run:\n", algname);
                  timesyshms(colinyEnd-colinyStart, &tms_colinyStart, &tms_colinyEnd);
                  pr(logFile, "\n");
                  (void) fflush(logFile);

                  pr(logFile, "Total number of Energy Evaluations: %d\n", (int)evaluate.evals() );
                  pr(logFile, "Total number of Iterations:        %d\n", (int)niters);

                  pr( logFile, UnderLine );
                  pr( logFile, "\n\n\tFINAL Coliny %s DOCKED STATE\n",algname );
                  pr( logFile,     "\t____________________________________\n\n\n" );

                  writeStateOfPDBQ( j, seed, FN_ligand, dock_param_fn, lig_center,
                        &(sHist[nconf]), ntor, &eintra, &einter, natom, atomstuff,
                        crd, emap, elec, 
                        charge, abs_charge, qsp_abs_charge,
                        ligand_is_inhibitor,
                        torsFreeEnergy,
                        vt, tlist, crdpdb, nonbondlist, e_internal,
                        type, Nnb, B_calcIntElec, q1q2,
                        map, inv_spacing, xlo, ylo, zlo, xhi, yhi, zhi,
                        B_template, template_energy, template_stddev,
                        outlev,
                        ignore_inter,
                        B_include_1_4_interactions, scale_1_4, sol_fn, parameterArray, unbound_internal_FE);

		  econf[nconf] = eintra + einter; // new2
		  
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

        } else {
            (void)fprintf(logFile, "NOTE: Command mode has been set, so optimization cannot be performed.\n\n");
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
        #ifdef USE_DOUBLE
            (void) sscanf( line, "%*s %lf %lf %lf", &lig_center[X], &lig_center[Y], &lig_center[Z]);
        #else
            (void) sscanf( line, "%*s %f %f %f", &lig_center[X], &lig_center[Y], &lig_center[Z]);
        #endif
        pr( logFile, "Small molecule center of rotation =\t" );
        pr( logFile, "(%+.3f, %+.3f, %+.3f)\n\n", lig_center[X], lig_center[Y], lig_center[Z]);
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
            mol.S.T.x = sInit.T.x = random_range( xlo, xhi );
            mol.S.T.y = sInit.T.y = random_range( ylo, yhi );
            mol.S.T.z = sInit.T.z = random_range( zlo, zhi );
        } else {
            B_RandomTran0 = FALSE;
            (void) sscanf( line,"%*s %lf %lf %lf", &(sInit.T.x), &(sInit.T.y), &(sInit.T.z));
            mol.S.T.x = sInit.T.x;
            mol.S.T.y = sInit.T.y;
            mol.S.T.z = sInit.T.z;
            if (outlev >= 0) {
                pr( logFile, "Initial translation =\t\t\t(%.3f, %.3f, %.3f) Angstroms\n", sInit.T.x, sInit.T.y, sInit.T.z );
            }
            (void) fflush(logFile);
        }
        break;

//______________________________________________________________________________

    case DPF_QUAT0:
        /*
        **  quat0
        **  Initial_quaternion,
        */
        (void) sscanf( line, "%*s %s", param[0]);
        for (i=0; i<6; i++) {
            param[0][i] = (char)tolower( (int)param[0][i] );
        }
        if (equal(param[0],"random",6)) {
            B_RandomQuat0 = TRUE;
            sInit.Q.nx  = random_range( -1.0,   1.0 );
            sInit.Q.ny  = random_range( -1.0,   1.0 );
            sInit.Q.nz  = random_range( -1.0,   1.0 );
            sInit.Q.ang = random_range(  0.0, 360.0 );

            if (outlev >= 0) {
                pr( logFile, "Each run will begin with a new, random initial quaternion.\n");
                }
        } else {
            B_RandomQuat0 = FALSE;

        }
        (void) sscanf( line, "%*s %lf %lf %lf %lf", &sInit.Q.nx, &sInit.Q.ny, &sInit.Q.nz, &sInit.Q.ang);
        if (outlev >= 0) {
            pr( logFile, "Initial quaternion,  q = [(x,y,z),w] =\t[ (%.3f, %.3f, %.3f), %.1f deg ],\n", sInit.Q.nx, sInit.Q.ny, sInit.Q.nz, sInit.Q.ang);
        }
        mol.S.Q.nx  = sInit.Q.nx;
        mol.S.Q.ny  = sInit.Q.ny;
        mol.S.Q.nz  = sInit.Q.nz;
        mol.S.Q.ang = sInit.Q.ang;

        mol.S.Q.ang = sInit.Q.ang = Rad( sInit.Q.ang ); /*convert to radians*/
        mkUnitQuat( &sInit.Q );
        mkUnitQuat( &(mol.S.Q) );

        if (outlev >= 0) {
            pr( logFile, "Quaternion Vector Normalized to:  %.3f %.3f %.3f\n", sInit.Q.nx, sInit.Q.ny, sInit.Q.nz );
        }
        (void) fflush(logFile);
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
                pr( logFile, "Number of torsions = %d\n", ndihed);
            }
            if ( ndihed != ntor ) {
                pr( logFile, "WARNING! You requested %d torsions, but I found %d in PDBQT-file specifications.\n", ndihed, ntor );
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
            sInit.ntor = nval = ndihed;
            for ( i=0; i<nval; i++ ) {
                sInit.tor[i] = random_range( -180.0, 180.0 );
            }
        } else {
            B_RandomDihe0 = FALSE;
            retval = (int)sscanf( line, torfmt, TOR_ARG_LIST );
            if (retval == 0) {
                pr( logFile, "WARNING!  Could not read any torsions!\n" );
            } else if (retval == EOF) {
                pr( logFile, "WARNING!  End of file encountered while reading dihe0 line\n");
            } else if (retval < ndihed) {
                pr( logFile, "WARNING!  Only %d initial torsion angles were detected on input line.\n",retval);
                pr( logFile, "WARNING!  I am sorry, you set 'ndihe', the number of dihedrals, to %d torsions.\n", ndihed);
            } else {
                if (outlev >= 0) {
                    pr( logFile, "%d initial torsion angles were detected on input line.\n", retval );
                }
            }
            nval = retval;
        }
        for ( i=0; i<nval; i++ ) {
            if (outlev >= 0) {
                pr( logFile, "\tInitial torsion %2d = %7.2f deg\n", (i+1), sInit.tor[i] ); /* sInit.tor is in degrees */
                /* Convert sInit.tor[i] into radians */
            }
            mol.S.tor[i] = sInit.tor[i] = Rad( sInit.tor[i] ); /* sInit.tor is now in radians  Added:05-01-95 */
        }
        (void) fflush(logFile);
        break;

//______________________________________________________________________________

    case DPF_TSTEP:
        /*
        **  tstep
        **  Translation_step,
        */
        #ifdef USE_DOUBLE
            retval = (int)sscanf( line, "%*s %lf %lf", &trnStep0, &trnStepFinal );
        #else
            retval = (int)sscanf( line, "%*s %f %f", &trnStep0, &trnStepFinal );
        #endif
        if (retval == 0) {
            pr( logFile, "WARNING!  Could not read any arguments!\n" );
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
        #ifdef USE_DOUBLE
            retval = (int)sscanf( line, "%*s %lf %lf", &qtwStep0, &qtwStepFinal );
        #else
            retval = (int)sscanf( line, "%*s %f %f", &qtwStep0, &qtwStepFinal );
        #endif
        if (retval == 0) {
            pr( logFile, "WARNING!  Could not read any arguments!\n" );
        } else if (retval == EOF) {
            pr( logFile, "WARNING!  End of file encountered!\n");
        } else if (retval > 0) {
            if (outlev >= 0) {
                pr( logFile, "Initial cycle, maximum quaternion angle step = +/- %-.1f deg\n", qtwStep0);
            }
            /* convert to radians */
            qtwStep0 = Rad( qtwStep0 );
        }
        if (retval == 2) {
            B_CalcQtwRF = TRUE;
            if (outlev >= 0) {
                pr( logFile, "Final cycle,   maximum quaternion angle step = +/- %-.1f deg\n", qtwStepFinal);
                pr( logFile, "Reduction factor will be calculated when number of cycles has been read in.\n");
            }
            /* convert to radians */
            qtwStepFinal = Rad( qtwStepFinal );
        }
        (void) fflush(logFile);
        break;

//______________________________________________________________________________

    case DPF_DSTEP:
        /*
        **  dstep
        **  Torsion_step,
        */
        #ifdef USE_DOUBLE
            retval = (int)sscanf( line, "%*s %lf %lf", &torStep0, &torStepFinal );
        #else
            retval = (int)sscanf( line, "%*s %f %f", &torStep0, &torStepFinal );
        #endif
        if (retval == 0) {
            pr( logFile, "WARNING!  Could not read any arguments!\n" );
        } else if (retval == EOF) {
            pr( logFile, "WARNING!  End of file encountered!\n");
        } else if (retval > 0) {
            if (outlev >= 0) {
                pr( logFile, "Initial cycle, maximum torsion angle step = +/- %-.1f deg\n", torStep0);
            }
            /* convert to radians */
            torStep0 = Rad( torStep0 );
        }
        if (retval == 2) {
            B_CalcTorRF = TRUE;
            if (outlev >= 0) {
                pr( logFile, "Final cycle,   maximum torsion angle step = +/- %-.1f deg\n", torStepFinal);
                pr( logFile, "Reduction factor will be calculated when number of cycles has been read in.\n");
            }
            /* convert to radians */
            torStepFinal = Rad( torStepFinal );
        }
        (void) fflush(logFile);
        break;

//______________________________________________________________________________

    case DPF_TRNRF:
        /*
        **  trnrf
        **  Translation reduction factor,
        */
        #ifdef USE_DOUBLE
            (void) sscanf( line, "%*s %lf", &trnFac );
        #else
            (void) sscanf( line, "%*s %f", &trnFac );
        #endif
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
        #ifdef USE_DOUBLE
            (void) sscanf( line, "%*s %lf", &qtwFac );
        #else
            (void) sscanf( line, "%*s %f", &qtwFac );
        #endif
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
        #ifdef USE_DOUBLE
            (void) sscanf( line, "%*s %lf", &torFac );
        #else
            (void) sscanf( line, "%*s %f", &torFac );
        #endif
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
        /*
        **  intnbp_r_eps
        **  Read internal energy parameters:
        **  Lennard-Jones and Hydrogen Bond Potentials,
        **  Using epsilon and r-equilibrium values...
        */
        #ifdef USE_DOUBLE
            (void) sscanf( line, "%*s %lf %lf %d %d", &Rij, &epsij, &xA, &xB );
        #else
            (void) sscanf( line, "%*s %f %f %d %d", &Rij, &epsij, &xA, &xB );
        #endif
        /* check that the Rij is reasonable */
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
            (void) fprintf( logFile, "Perhaps you meant to use \"intnbp_coeffs\" instead of \"intnbp_r_eps\"?\n\n");
            /* GMM COMMENTED OUT FOR DAVE GOODSELL, MUTABLE ATOMS
             * exit(-1); */
        }
        /* Defend against division by zero... */
        if (xA != xB) {
            cA = (tmpconst = epsij / (FloatOrDouble)(xA - xB)) * pow( (double)Rij, (double)xA ) * (FloatOrDouble)xB;
            cB = tmpconst * pow( (double)Rij, (double)xB ) * (FloatOrDouble)xA;
            intnbtable( &B_havenbp, &a1, &a2, num_atom_types, ligand_atom_types, cA, cB, xA, xB, e_internal, sol_fn, AD4.coeff_desolv, sigma );
        } else {
            pr(logFile,"WARNING: Exponents must be different, to avoid division by zero!\n\tAborting...\n");
            exit(-1);
        }
        (void) fflush(logFile);
        break;

//______________________________________________________________________________

    case DPF_INTNBP_COEFFS:
        /*
        **  intnbp_coeffs
        **  Read internal energy parameters:
        **  Lennard-Jones and Hydrogen Bond Potentials,
        **  Using coefficients...
        */
        #ifdef USE_DOUBLE
            (void) sscanf( line, "%*s %lf %lf %d %d", &cA, &cB, &xA, &xB );
        #else
            (void) sscanf( line, "%*s %f %f %d %d", &cA, &cB, &xA, &xB );
        #endif

        /* Defend against division by zero... */
        if (xA != xB) {
            intnbtable( &B_havenbp, &a1, &a2, num_atom_types, ligand_atom_types, cA, cB, xA, xB, e_internal, sol_fn, AD4.coeff_desolv, sigma );
        } else {
            pr(logFile,"WARNING: Exponents must be different. Aborting...\n");
            exit(-1);
        }

        (void) fflush(logFile);
        break;

//______________________________________________________________________________

    case DPF_RT0:
        /*
        **  rt0
        **  Initial Temperature,
        */
        #ifdef USE_DOUBLE
            (void) sscanf( line, "%*s %lf", &RT0 );
        #else
            (void) sscanf( line, "%*s %f", &RT0 );
        #endif
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
        #ifdef USE_DOUBLE
            (void) sscanf( line, "%*s %lf", &RTFac);
        #else
            (void) sscanf( line, "%*s %f", &RTFac);
        #endif
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
            prStr( error_message, "ERROR: %d runs were requested, but AutoDock is only dimensioned for %d.\nChange \"MAX_RUNS\" in \"constants.h\".", nruns, MAX_RUNS);
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
        (void) sscanf( line, "%*s %1s", &selminpar );
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
        #ifdef USE_DOUBLE
            (void) sscanf( line, "%*s %lf", &clus_rms_tol);
        #else
            (void) sscanf( line, "%*s %f", &clus_rms_tol);
        #endif
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
        (void) sscanf( line, "%*s %1s", &out_acc_rej );
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
        #ifdef USE_DOUBLE
            (void) sscanf( line, "%*s %lf", &WallEnergy );
        #else
            (void) sscanf( line, "%*s %f", &WallEnergy );
        #endif
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

    case DPF_INTELEC:
        /*
        **  intelec
        **  Calculate internal electrostatic energies...
        */
        B_calcIntElec = TRUE;
        if (outlev >= 0) {
            pr( logFile, "Internal electrostatic energies will be calculated.\n\n");
        }
        #ifdef USE_DOUBLE
            retval = sscanf( line, "%*s %lf", &AD3_FE_coeff_estat );
        #else
            retval = sscanf( line, "%*s %f", &AD3_FE_coeff_estat );
        #endif
        if (retval == 1) {
            if (outlev >= 0) {
                pr(logFile, "WARNING!  Internal electrostatics will NOT be scaled by a factor of %.4f -- this coefficient is ignored in AutoDock 4;\n", AD3_FE_coeff_estat);
                pr(logFile, "          the coefficient that will actually be used should be set in the parameter library file, \"%s\".\n", FN_parameter_library);
                pr(logFile, "          In this case, the coefficient for the electrostatics term is %.4f\n", AD4.coeff_estat);
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
        #ifdef USE_DOUBLE
            (void) sscanf( line, "%*s %d %lf %lf", &I_tor, &F_torPref, &F_torHWdth);
        #else
            (void) sscanf( line, "%*s %d %f %f", &I_tor, &F_torPref, &F_torHWdth);
        #endif
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
            F_torPref = Wrp(ModDeg(F_torPref));
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
                    US_energy = (unsigned short) (((FloatOrDouble)US_torBarrier) * (1.0 - exp(F_lnH * F_Aova*F_Aova)));
                    /*
                    ** if F_A(<-180.or>180), wrap to -180to180,
                    */
                    F_tor = Wrp(ModDeg(F_A));
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
        #ifdef USE_DOUBLE
            retval = sscanf( line, "%*s %lf %d", &e0max, &MaxRetries );
        #else
            retval = sscanf( line, "%*s %f %d", &e0max, &MaxRetries );
        #endif
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
        if (!command_mode) {
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
            if (B_havenbp) {
                nbe( atm_typ_str, e_internal, num_atom_types );
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
            /* ___________________________________________________________________
            **
            ** Begin the automated docking simulation,
            ** ___________________________________________________________________
            */
            simanneal(&nconf, Nnb, WallEnergy, atomstuff, charge, abs_charge, qsp_abs_charge, B_calcIntElec,
                    q1q2,crd,crdpdb,dock_param_fn,e_internal,econf,B_either,
                    elec,emap,xhi,yhi,zhi,ncycles,inv_spacing,nruns,jobStart,
                    xlo,ylo,zlo,map,
                    naccmax, natom, nonbondlist, nrejmax, ntor1, ntor, outlev,
                    sInit, sHist,   qtwFac, B_qtwReduc, qtwStep0,
                    B_selectmin,FN_ligand,lig_center,RT0,B_tempChange,RTFac,
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
                    sol_fn, parameterArray, unbound_internal_FE);
            (void) fflush(logFile);
        } else {
            (void)fprintf(logFile, "NOTE: Command mode has been set, so simulated annealing cannot be performed.\n\n");
        }
        break;

//______________________________________________________________________________

    case DPF_SET_GA:

      if (GlobalSearchMethod != NULL) {
          pr(logFile, "Deleting the previous settings for the Genetic Algorithm.\n");
          delete GlobalSearchMethod;
      }

      pr(logFile, "Creating a new Genetic Algorithm object with the current parameters.\n\n");
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
      if (debug > 0) {
          pr( logFile, "\n\tOutput every %u generations.\n", outputEveryNgens );
      }
      GlobalSearchMethod = new Genetic_Algorithm(e_mode, s_mode, c_mode, w_mode, elitism, c_rate, m_rate, window_size, num_generations, outputEveryNgens);
      ((Genetic_Algorithm *)GlobalSearchMethod)->mutation_values(low, high, alpha, beta);
      ((Genetic_Algorithm *)GlobalSearchMethod)->initialize(pop_size, 7+sInit.ntor);

      (void) fflush(logFile);
      break;
//______________________________________________________________________________

    case DPF_SET_SW1:

      if (LocalSearchMethod != NULL) {
          pr(logFile, "Deleting the previous settings for the local search Solis-Wets algorithm (SW1 object).\n");
          delete LocalSearchMethod;
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
      }

      pr(logFile, "Creating a new Local Search object using the pseudo-Solis-Wets algorithm (pSW1) with the current settings.\n\n");

      //  Allocate space for the variable rho's
      rho_ptr = new FloatOrDouble[7+sInit.ntor];
      lb_rho_ptr = new FloatOrDouble[7+sInit.ntor];

      //  Initialize the rho's corresponding to the translation
      for (j=0; j<3; j++) {
         rho_ptr[j] = trnStep0;
         lb_rho_ptr[j] = trnStepFinal;
      }

      //  Initialize the rho's corresponding to the quaterion
      for (; j<7; j++) {
         rho_ptr[j] = qtwStep0;
         lb_rho_ptr[j] = qtwStepFinal;
      }

      //  Initialize the rho's corresponding to the torsions
      for (; j<7+sInit.ntor; j++) {
         rho_ptr[j] = torStep0;
         lb_rho_ptr[j] = torStepFinal;
      }

      LocalSearchMethod = new Pseudo_Solis_Wets1(7+sInit.ntor, max_its, max_succ, max_fail, 2.0, 0.5, search_freq, rho_ptr, lb_rho_ptr);

      (void) fflush(logFile);
      break;

//______________________________________________________________________________

    case DPF_GALS:
        (void) fflush( logFile );
        /*
        ** Genetic Algorithm-Local search
        */
        if (!command_mode) {
            (void) sscanf( line, "%*s %d", &nruns );

            if ( nruns > MAX_RUNS ) {

                prStr( error_message, "ERROR: %d runs requested, but only dimensioned for %d.\nChange \"MAX_RUNS\" in \"constants.h\".", nruns, MAX_RUNS);
                stop( error_message );
                exit( -1 );

            } else if ((GlobalSearchMethod==NULL)||(LocalSearchMethod==NULL)) {

                prStr(error_message, "ERROR:  Need to use \"set_ga\" to allocate both Global Optimization object AND Local Optimization object.\n");
                stop(error_message);
                exit(-1);

            }

            pr( logFile, "Number of requested LGA dockings = %d run%c\n", nruns, (nruns > 1)?'s':' ');

            evaluate.setup(crd, charge, abs_charge, qsp_abs_charge, type, natom, map, inv_spacing,
              elec, emap,
              xlo, xhi, ylo, yhi, zlo, zhi, nonbondlist, e_internal, Nnb,
              B_calcIntElec, q1q2, B_isGaussTorCon, B_isTorConstrained,
              B_ShowTorE, US_TorE, US_torProfile, vt, tlist, crdpdb, sInit, mol,
              B_template, template_energy, template_stddev,
              ignore_inter,
              B_include_1_4_interactions, scale_1_4, sol_fn, parameterArray, unbound_internal_FE);

	    if(write_stateFile){
	      fprintf(stateFile,"\t<run_requested>%d</run_requested>\n",nruns);
	      fprintf(stateFile,"\t<runs>\n");
	    }
            for (j=0; j<nruns; j++) {
                j1 = j + 1;

                (void) fprintf( logFile, "\n\n\tBEGINNING LAMARCKIAN GENETIC ALGORITHM DOCKING\n");
                (void) fflush( logFile );
                pr( logFile, "\nRun:\t%d / %d\n", j1, nruns );

                // Update time-based RNG seeds...
                if (timeSeedIsSet[0] == 'T') {
                    pr(logFile, "Updating First Time-dependent Seed Now.\n");
                    seed[0] = (FourByteLong)time( &time_seed );
                }
                if (timeSeedIsSet[1] == 'T') {
                    pr(logFile, "Updating Second Time-dependent Seed Now.\n");
                    seed[1] = (FourByteLong)time( &time_seed );
                }
                setall(seed[0], seed[1]);
                initgn(-1); // Reinitializes the state of the current random number generator
                pr(logFile, "Seeds:\t%ld %ld\nDate:\t", seed[0], seed[1]);
                printdate( logFile, 2 );
                (void) fflush( logFile );

                gaStart = times( &tms_gaStart );

                //  Can get rid of the following line
                ((Genetic_Algorithm *)GlobalSearchMethod)->initialize(pop_size, 7+sInit.ntor);

                // Reiterate output level...
                pr(logFile, "Output level is set to %d.\n\n", outlev);

                // Start Lamarckian GA run
                sHist[nconf] = call_glss( GlobalSearchMethod, LocalSearchMethod, sInit,
                                          num_evals, pop_size, xlo, xhi,
                                          ylo, yhi, zlo, zhi, outlev,
                                          outputEveryNgens, &mol,
                                          B_template, B_RandomTran0,
                                                           B_RandomQuat0, B_RandomDihe0);
                                          // State of best individual at end
                                          // of GA-LS run is returned.
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


                writeStateOfPDBQ( j,seed,  FN_ligand, dock_param_fn, lig_center,
                    &(sHist[nconf]), ntor, &eintra, &einter, natom, atomstuff,
                    crd, emap, elec, 
                    charge, abs_charge, qsp_abs_charge, 
                    ligand_is_inhibitor,
                    torsFreeEnergy,
                    vt, tlist, crdpdb, nonbondlist, e_internal,
                    type, Nnb, B_calcIntElec, q1q2,
                    map, inv_spacing, xlo, ylo, zlo, xhi, yhi, zhi,
                    B_template, template_energy, template_stddev,
                    outlev,
                    ignore_inter,
                    B_include_1_4_interactions, scale_1_4, sol_fn, parameterArray, unbound_internal_FE);

                econf[nconf] = eintra + einter; // new2

                ++nconf;

                pr( logFile, UnderLine );
            } // Next LGA run
	    if(write_stateFile){
	       fprintf(stateFile,"\t</runs>\n");
	       (void) fflush(stateFile);
	    }
            (void) fflush(logFile);
        } else {
            (void)fprintf(logFile, "NOTE: Command mode has been set.  Sorry, genetic algorithm-local search cannot be performed.\n\n");
        }
        break;

//______________________________________________________________________________

    case DPF_LS:
       (void) sscanf(line, "%*s %d", &nruns);

        if (!command_mode) {
            if ( nruns > MAX_RUNS ) {

               prStr( error_message, "ERROR: %d runs requested, but only dimensioned for %d.\nChange \"MAX_RUNS\" in \"constants.h\".", nruns, MAX_RUNS);
               stop( error_message );
               exit( -1 );

           } else if (LocalSearchMethod==NULL) {

               prStr(error_message, "ERROR:  Need to use \"set_sw1\" or \"set_psw1\" to allocate a Local Optimization object.\n");
               stop(error_message);
               exit(-1);
           }

           pr( logFile, "Number of Local Search (LS) only dockings = %d run%c\n", nruns, (nruns > 1)?'s':' ');


           evaluate.setup(crd, charge, abs_charge, qsp_abs_charge, type, natom, map, inv_spacing,
              elec, emap,
              xlo, xhi, ylo, yhi, zlo, zhi, nonbondlist,
              e_internal, Nnb, B_calcIntElec, q1q2,B_isGaussTorCon,B_isTorConstrained,
              B_ShowTorE, US_TorE, US_torProfile, vt, tlist, crdpdb, sInit, mol,
              B_template, template_energy, template_stddev,
              ignore_inter,
              B_include_1_4_interactions, scale_1_4, sol_fn, parameterArray, 
              unbound_internal_FE);

	   if(write_stateFile){
	     fprintf(stateFile,"\t<run_requested>%d</run_requested>\n",nruns);
	     fprintf(stateFile,"\t<runs>\n");
	   }

           for (j=0; j<nruns; j++) {

               pr( logFile, "\nDoing Local Search run: %d / %d.\n", j+1, nruns );
               if (timeSeedIsSet[0] == 'T') {
                   seed[0] = (FourByteLong)time( &time_seed );
               }
               if (timeSeedIsSet[1] == 'T') {
                   seed[1] = (FourByteLong)time( &time_seed );
               }
               setall(seed[0], seed[1]);
               initgn(-1); // Reinitializes the state of the current generator

               pr(logFile, "Seeds: %ld %ld\n", seed[0], seed[1]);
               (void) fflush( logFile );
               gaStart = times( &tms_gaStart );

               sHist[nconf] = call_ls(LocalSearchMethod, sInit, pop_size, &mol);

               pr(logFile, "There were %lu Energy Evaluations.\n", evaluate.evals());

               gaEnd = times( &tms_gaEnd );
               pr( logFile, "Time taken for this Local Search (LS) run:\n");
               timesyshms( gaEnd - gaStart, &tms_gaStart, &tms_gaEnd );
               pr( logFile, "\n");
               (void) fflush( logFile );

               pr( logFile, "\n\n\tFINAL LOCAL SEARCH DOCKED STATE\n" );
               pr( logFile,     "\t_______________________________\n\n\n" );

               writeStateOfPDBQ( j,seed, FN_ligand, dock_param_fn, lig_center,
                    &(sHist[nconf]), ntor, &eintra, &einter, natom, atomstuff,
                    crd, emap, elec, 
                    charge, abs_charge, qsp_abs_charge, 
                    ligand_is_inhibitor,
                    torsFreeEnergy,
                    vt, tlist, crdpdb, nonbondlist, e_internal,
                    type, Nnb, B_calcIntElec, q1q2,
                    map, inv_spacing, xlo, ylo, zlo, xhi, yhi, zhi,
                    B_template, template_energy, template_stddev,
                    outlev,
                    ignore_inter,
                    B_include_1_4_interactions, scale_1_4, sol_fn, parameterArray, unbound_internal_FE);

               econf[nconf] = eintra + einter; // new2

               ++nconf;

               pr( logFile, UnderLine );

           } // Next run
	   if(write_stateFile){
	     fprintf(stateFile,"\t</runs>\n");
	     (void) fflush(stateFile);
	   }
           (void) fflush(logFile);
       } else {
            (void)fprintf(logFile, "NOTE: Command mode has been set, so local search cannot be performed.\n\n");
       }
       break;

//______________________________________________________________________________

    case DPF_GS:
      (void) sscanf(line, "%*s %d", &nruns);

      if (!command_mode) {
          if (nruns>MAX_RUNS) {

              prStr(error_message, "ERROR:  %d runs requested, but only dimensioned for %d.\nChange \"MAX_RUNS\" in \"constants.h\".", nruns, MAX_RUNS);
              stop(error_message);
              exit(-1);

          } else if (GlobalSearchMethod==NULL) {

              prStr(error_message, "ERROR:  Need to use \"set_ga\" to allocate a Global Optimization object.\n");
              stop(error_message);
              exit(-1);
          }

          pr(logFile, "Number of Genetic Algorithm (GA) only dockings = %d run%c\n", nruns, (nruns>1)?'s':' ');


          evaluate.setup(crd, charge, abs_charge, qsp_abs_charge, type, natom, map, inv_spacing,
             elec, emap,
             xlo, xhi, ylo, yhi, zlo, zhi, nonbondlist,
             e_internal, Nnb, B_calcIntElec, q1q2, B_isGaussTorCon,B_isTorConstrained,
             B_ShowTorE, US_TorE, US_torProfile, vt, tlist, crdpdb, sInit, mol,
             B_template, template_energy, template_stddev,
             ignore_inter,
             B_include_1_4_interactions, scale_1_4, sol_fn, parameterArray, unbound_internal_FE);

	  if(write_stateFile){
	    fprintf(stateFile,"\t<run_requested>%d</run_requested>\n",nruns);
	    fprintf(stateFile,"\t<runs>\n");
	  }
          for (j=0; j<nruns; j++) {

              fprintf( logFile, "\n\n\tBEGINNING GENETIC ALGORITHM DOCKING\n");
              pr(logFile, "\nDoing Genetic Algorithm run:  %d/%d.\n", j+1, nruns);

              if (timeSeedIsSet[0] == 'T') {
                  seed[0] = (FourByteLong)time( &time_seed );
              }
              if (timeSeedIsSet[1] == 'T') {
                  seed[1] = (FourByteLong)time( &time_seed );
              }
              setall(seed[0], seed[1]);
              initgn(-1); // Reinitializes the state of the current generator

              pr(logFile, "Seeds:  %ld %ld\n", seed[0], seed[1]);
              (void) fflush(logFile);

              gaStart = times(&tms_gaStart);

              sHist[nconf] = call_gs(GlobalSearchMethod, sInit, num_evals, 
                      pop_size, xlo, xhi, ylo, yhi, zlo, zhi, &mol, outputEveryNgens);

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

              writeStateOfPDBQ( j, seed, FN_ligand, dock_param_fn, lig_center,
                    &(sHist[nconf]), ntor, &eintra, &einter, natom, atomstuff,
                    crd, emap, elec, 
                    charge, abs_charge, qsp_abs_charge, 
                    ligand_is_inhibitor,
                    torsFreeEnergy,
                    vt, tlist, crdpdb, nonbondlist, e_internal,
                    type, Nnb, B_calcIntElec, q1q2,
                    map, inv_spacing, xlo, ylo, zlo, xhi, yhi, zhi,
                    B_template, template_energy, template_stddev,
                    outlev,
                    ignore_inter,
                    B_include_1_4_interactions, scale_1_4, sol_fn, parameterArray, unbound_internal_FE);

              econf[nconf] = eintra + einter; // new2

              ++nconf;

              pr( logFile, UnderLine );

          } // Next run
	  if(write_stateFile){
	    fprintf(stateFile,"\t</runs>\n");
	    (void) fflush(stateFile);
	  }
          (void) fflush(logFile);
      } else {
            (void)fprintf(logFile, "NOTE: Command mode has been set, so global search cannot be performed.\n\n");
      }
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
       #ifdef USE_DOUBLE
           (void) sscanf(line, "%*s %lf", &m_rate);
       #else
           (void) sscanf(line, "%*s %f", &m_rate);
       #endif
       pr(logFile, "The mutation rate is %f.\n", m_rate);
        (void) fflush(logFile);
       break;

//______________________________________________________________________________

    case GA_crossover_rate:
       #ifdef USE_DOUBLE
           (void) sscanf(line, "%*s %lf", &c_rate);
       #else
           (void) sscanf(line, "%*s %f", &c_rate);
       #endif
       pr(logFile, "The crossover rate is %f.\n", c_rate);
        (void) fflush(logFile);
       break;

//______________________________________________________________________________

    case GA_Cauchy_alpha:
       #ifdef USE_DOUBLE
           (void) sscanf(line, "%*s %lf", &alpha);
       #else
           (void) sscanf(line, "%*s %f", &alpha);
       #endif
       pr(logFile, "The alpha parameter (for the Cauchy distribution) is being set to %f.\n",
          alpha);
        (void) fflush(logFile);
       break;

//______________________________________________________________________________

    case GA_Cauchy_beta:
       #ifdef USE_DOUBLE
           (void) sscanf(line, "%*s %lf", &beta);
       #else
           (void) sscanf(line, "%*s %f", &beta);
       #endif
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
       #ifdef USE_DOUBLE
           (void) sscanf(line, "%*s %lf", &rho);
       #else
           (void) sscanf(line, "%*s %f", &rho);
       #endif
       pr(logFile, "rho is set to %f.\n", rho);
        (void) fflush(logFile);
      break;

//______________________________________________________________________________

    case SW_lb_rho:
        #ifdef USE_DOUBLE
            (void) sscanf(line, "%*s %lf", &lb_rho);
        #else
            (void) sscanf(line, "%*s %f", &lb_rho);
        #endif
        pr(logFile, "rho will never get smaller than %f.\n", lb_rho);
        (void) fflush(logFile);
        break;

//______________________________________________________________________________

    case LS_search_freq:
        #ifdef USE_DOUBLE
            (void) sscanf(line, "%*s %lf", &search_freq);
        #else
            (void) sscanf(line, "%*s %f", &search_freq);
        #endif
        pr(logFile, "Local search will be performed with frequency %f.\n", search_freq);
        (void) fflush(logFile);
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
        if (!command_mode) {
            analysis( Nnb, atomstuff, charge, abs_charge, qsp_abs_charge, B_calcIntElec, q1q2, clus_rms_tol,
                  crdpdb, e_internal, inv_spacing, map, econf, nruns,
                  xlo,ylo,zlo,
                  natom, nonbondlist, nconf, ntor, sHist, FN_ligand,
                  lig_center, B_symmetry_flag, tlist, type, vt, FN_rms_ref_crds,
                  torsFreeEnergy, B_write_all_clusmem, ligand_is_inhibitor,
                  B_template, template_energy, template_stddev, outlev,
                  ignore_inter, B_include_1_4_interactions, scale_1_4, 
                  sol_fn, parameterArray, unbound_internal_FE);
            (void) fflush(logFile);
        } else {
            (void)fprintf(logFile, "NOTE: Command mode has been set, so cluster analysis cannot be performed.\n\n");
        }
        break;

//______________________________________________________________________________

    case DPF_TORSDOF:
        /*
        ** torsdof %d %f
        */
        #ifdef USE_DOUBLE
            retval = sscanf( line, "%*s %d %lf", &ntorsdof, &torsdoffac );
            if (retval == 2) {
                pr( logFile, "WARNING:  The torsional DOF coefficient is now read in from the parameter file; the value specified here (%.4lf) will be ignored.\n\n", torsdoffac);
            }
        #else
            retval = sscanf( line, "%*s %d %f", &ntorsdof, &torsdoffac );
            if (retval == 2) {
                pr( logFile, "WARNING:  The torsional DOF coefficient is now read in from the parameter file; the value specified here (%.4f) will be ignored.\n\n", torsdoffac);
            }
        #endif
        pr( logFile, "Number of torsional degrees of freedom = %d\n", ntorsdof);
        pr( logFile, "Note: this must exclude any torsions involving -OH and -NH2 groups.\n\n");
        pr( logFile, "Free energy coefficient for torsional degrees of freedom = %.4f  (from \"%s\")\n\n", AD4.coeff_tors, FN_parameter_library);

        torsFreeEnergy = (FloatOrDouble)ntorsdof * AD4.coeff_tors;

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
        (void) investigate( Nnb, charge, abs_charge, qsp_abs_charge, B_calcIntElec, q1q2,
                crd, crdpdb, e_internal, xhi, yhi, zhi,
                inv_spacing, maxTests, xlo, ylo, zlo,
                map, natom, nonbondlist, ntor,
                outlev, tlist, type, vt, B_isGaussTorCon, US_torProfile,
                B_isTorConstrained, B_ShowTorE, US_TorE,
                F_TorConRange, N_con, B_symmetry_flag, FN_rms_ref_crds,
                OutputEveryNTests, NumLocalTests, trnStep0, torStep0,
                ignore_inter,
                B_include_1_4_interactions, scale_1_4, 
                sol_fn, parameterArray, unbound_internal_FE);

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

    case DPF_TEMPL_ENERGY:
        /*
        ** Read in the template energy file
        ** Perform a template docking...
        */
        B_template = TRUE;
        pr( logFile, "\nAutoDock will use a template scoring function to perform docking.\n\n");
        (void) sscanf( line, "%*s %s", FN_template_energy_file);
        pr( logFile, "Template energy file \"%s\" will be read in.  Expecting %d pairs of values\n\n", FN_template_energy_file, natom);
        if ((template_energy_file = ad_fopen(FN_template_energy_file, "r")) == NULL) {
            pr(logFile, "\n%s: ERROR:  I'm sorry, I cannot find or open \"%s\".\n\n", programname, FN_template_energy_file);
        }
        curatm = 0;
        while (fgets(line_template, LINE_LEN, template_energy_file) != NULL) {
            #ifdef USE_DOUBLE
                retval = (int)sscanf(line_template, "%lf %lf", &template_energy[curatm], &template_stddev[curatm]);
            #else
                retval = (int)sscanf(line_template, "%f %f", &template_energy[curatm], &template_stddev[curatm]);
            #endif
            if (retval != 2) {
                pr(logFile, "\nWARNING: AutoDock expects the template energy file to have two values on each line: the energy and the standard deviation.  %d values were found on line %d.\n\n", retval, curatm+1);
            }
            if (template_stddev[curatm] == 0.0) {
                pr(logFile, "\nWARNING: the template energy standard deviation for atom %d is zero! This is not allowed, and will be set to 1.0.\n\n", curatm+1);
                template_stddev[curatm] = 1.0;
            }
            curatm++;
        }
        if (curatm != natom) {
            pr(logFile, "\nWARNING: The number of values in the template energy file (%d) does not match the number of atoms in the input PDBQT file (%d).\n", curatm, natom);
        }
        (void) fclose(template_energy_file);
        (void) fflush(logFile);
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
        if (B_found_move_keyword == TRUE) {  // If we have found the move keyword already, warn the user that his command should be given first!
            pr(logFile, "\nWARNING:  This command will be ignored.\n\nYou must put this command _before_ the \"move ligand.pdbqt\" command, since this command affects how the PDBQT file will be interpreted.\n\n");
        }
        #ifdef USE_DOUBLE
            (void) sscanf( line, "%*s %lf", &scale_1_4 );
        #else
            (void) sscanf( line, "%*s %f", &scale_1_4 );
        #endif
        B_include_1_4_interactions = TRUE;
        print_1_4_message(logFile, B_include_1_4_interactions, scale_1_4);
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

        par_lib_FN_found = sscanf( line, "%*s %s", FN_parameter_library );
        (void) fflush(logFile);

        read_parameter_library(FN_parameter_library, outlev);

        break;

/*____________________________________________________________________________*/

    case DPF_LIGAND_TYPES:
        /*
         *  Read in the ligand atom type names
         *
         *  ligand_types C HD OA P               # ligand atom type names
         *
         *  The order of the arguments is the index that will 
         *  be used for look up in the grid maps, "map_index".
         */
        
        //  Use "parsetypes" to read in the atom types;
        //  
        //  The array "ligand_atom_type_ptrs" is returned, having been filled with pointers
        //  to the beginning of each "atom type word" (not atom type characters);
        //  an atom type can be either 1 or 2 characters long.
        //
        //  "atm_typ_str" used to serve a similar role to "ligand_atom_types" in AD3.
        num_atom_types = parsetypes(line, ligand_atom_type_ptrs, MAX_ATOM_TYPES);

        for (i=0; i<num_atom_types; i++) {
            strcpy(ligand_atom_types[i], ligand_atom_type_ptrs[i]);
#ifdef DEBUG
            (void) fprintf(logFile, "%d %s ->%s\n",i, ligand_atom_type_ptrs[i], ligand_atom_types[i]);
#endif
        }

        if (num_atom_types > 0) {
            B_atom_types_found = TRUE;
        } else {
            prStr( error_message, "%s:  ERROR! No atom types have been found; we cannot continue without this information!\n\n", programname );
            pr_2x( logFile, stderr, error_message );
            exit(-1);
        }

        if (debug > 0) {
            for (i=0; i<num_atom_types; i++) {
                (void) fprintf(logFile, "ligand_atom_types[%d] = %s\n", i, ligand_atom_types[i] );
            }
        }

        // For all ligand atom types... set up the map_index
        // "ligand_types"
        for (i=0; i<num_atom_types; i++) {
            foundParameter = apm_find(ligand_atom_types[i]);
            if (foundParameter != NULL ) {
                // We have found this atom type's parameters
                // Now we can set the ParameterEntry's map_index to the
                // 0-based index it had in the list of ligand types supplied
                // in this DPF line:
                foundParameter->map_index = i;
                parameterArray[i] = *(foundParameter);
                if (outlev > 0) {
                    (void) fprintf( logFile, "Parameters found for ligand type \"%s\" (grid map index = %d, weighted epsilon = %6.4f)", foundParameter->autogrid_type, foundParameter->map_index, foundParameter->epsij );
                    if (par_lib_FN_found == 1) {
                        pr( logFile, " in parameter library \"%s\".\n", FN_parameter_library );
                    } else {
                        pr( logFile, "\n");
                    }
                }
            } else {
                // We could not find this parameter -- return error here
                prStr( error_message,"%s: ERROR:  Unknown ligand atom type \"%s\"; add parameters for it to the parameter library first!\n", programname, ligand_atom_types[i]);
                pr_2x( logFile, stderr, error_message );
                if (par_lib_FN_found == 1) {
                    prStr( error_message,"%s:         Edit the parameter library file \"%s\" and try again.\n", programname, FN_parameter_library );
                    pr_2x( logFile, stderr, error_message );
                }
                exit(-1);
            } // if / else apm_find
        } // for i
        pr( logFile, "\n\n");

        ElecMap = num_atom_types;
        DesolvMap = num_atom_types + 1;

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
                    cA = (tmpconst = epsij / (FloatOrDouble)(xA - xB)) * pow( (double)Rij, (double)xA ) * (FloatOrDouble)xB;
                    cB = tmpconst * pow( (double)Rij, (double)xB ) * (FloatOrDouble)xA;
                    intnbtable( &B_havenbp, &a1, &a2, num_atom_types, ligand_atom_types, cA, cB, xA, xB, e_internal, sol_fn, AD4.coeff_desolv, sigma );
                } else {
                    pr(logFile,"WARNING: Exponents must be different, to avoid division by zero!\n\tAborting...\n");
                    exit(-1);
                }
                (void) fflush(logFile);

            } // for j
        } // for i
        break;

/*____________________________________________________________________________*/

    case DPF_UNBOUND:
        /*
         * unbound 0.0
         */
#ifdef USE_DOUBLE
        (void) sscanf( line, "%*s %lf", &unbound_internal_FE );
        pr(logFile, "\nThe internal free energy of the unbound state is %+.3lf kcal/mol\n\n", unbound_internal_FE);
#else
        (void) sscanf( line, "%*s %f", &unbound_internal_FE );
        pr(logFile, "\nThe internal free energy of the unbound state is %+.3f kcal/mol\n\n", unbound_internal_FE);
#endif
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
        retval = sscanf(line, "%*s %s %d", dummy_FN_ligand, &dummy_oldpdbq);
        if (retval >= 1) {
            pr(logFile, "WARNING  -- it will not read in the PDBQT file specified on the \"epdb\" command line.\n");
        }

        /*
        (void) sscanf(line, "%*s %s %d", FN_ligand, &oldpdbq);
        pr(logFile, "epdb %s %d\n\n", FN_ligand, oldpdbq);

        natom = 0;

        print_1_4_message(logFile, B_include_1_4_interactions, scale_1_4);
        //
        mol = readPDBQT(line,
                        num_atom_types,
                        &natom,
                        crdpdb, charge, &B_haveCharges,
                        type, bond_index,
                        pdbaname, FN_ligand, atomstuff, Htype,
                        &B_constrain_dist, &atomC1, &atomC2,
                        &sqlower, &squpper,
                        &ntor1, &ntor, tlist, vt,
                        &Nnb, nonbondlist,
                        jobStart, tms_jobStart, hostnm, &ntorsdof, outlev,
                        ignore_inter,
                        B_include_1_4_interactions);
        //
        // pre-calculate some values we will need later in computing the desolvation energy
        //
        for (i=0; i<natom; i++) {
            abs_charge[i] = fabs(charge[i]);
            qsp_abs_charge[i] = qsolpar * abs_charge[i];
        }
        */

        for ( i=0; i<true_ligand_atoms; i++ ) {
            for (xyz = 0;  xyz < SPACE;  xyz++) {
                crdpdb[i][xyz] += lig_center[xyz];
            }
        }
        // determine if any atoms are outside the grid box
        atoms_outside = FALSE;
        for (i=0; i<natom; i++) {
            outside = is_out_grid(crdpdb[i][X], crdpdb[i][Y], crdpdb[i][Z]);
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
                pr(logFile, "----------------------------\n\n");
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
        //
        sInit.ntor = mol.S.ntor;
        //
        if (ntor > 0) {
            eintra = eintcalPrint(nonbondlist, e_internal, crdpdb, Nnb, B_calcIntElec, q1q2, B_include_1_4_interactions, scale_1_4, qsp_abs_charge, sol_fn, parameterArray, unbound_internal_FE);
        } else {
            eintra = 0.0L;
        }
        if (outside) {
            etotal = (einter = outsidetrilinterp4byatom(crdpdb, charge, abs_charge, type, natom, map, inv_spacing, elec,emap, xlo,ylo,zlo, xhi,yhi,zhi, map_center[0], map_center[1], map_center[2], ignore_inter)) + eintra; // gmm 2001.11.07
        } else {
            etotal = (einter = trilinterp4(crdpdb, charge, abs_charge, type, natom, map, inv_spacing, elec, emap, xlo, ylo, zlo, ignore_inter)) + eintra;
        }
        pr(logFile, "\n\n\t\tIntermolecular Energy Analysis\n");
        pr(logFile,     "\t\t==============================\n\n\n");
        pr(logFile, "Atom  NB.+ Elec.  Non-bonded  Electrosta  Partial          Coordinates         \n");
        pr(logFile, "Type    Energy      Energy    tic Energy  Charge      x         y         z    \n");
        pr(logFile, "----  ----------  ----------  ----------  -------  --------  --------  --------\n");
        /*          "1234  0123456789  0123456789  0123456789  1234567  12345678  12345678  12345678"*/
        /*          "----  ----------  ----------  ----------  -------  --------  --------  --------"*/
        emap_total = 0.;
        elec_total = 0.;
        charge_total = 0.;
        etot = 0.;

        for (i = 0;  i < natom;  i++) {
            etot = emap[i] + elec[i];
            pr(logFile, "%4d  %10.2f  %10.2f  %10.2f  %7.3f  %8.4f  %8.4f  %8.4f\n", (type[i]+1), etot, emap[i], elec[i], charge[i], crdpdb[i][X], crdpdb[i][Y], crdpdb[i][Z]);
            emap_total += emap[i];
            elec_total += elec[i];
            charge_total += charge[i];
        } /*i*/
        pr(logFile, "      ----------  ----------  ----------  -------\n");
#ifdef USE_DOUBLE
        pr(logFile, "Total %10.2lf  %10.2lf  %10.2lf  %7.3lf\n\n", (emap_total + elec_total), emap_total, elec_total, charge_total);
        pr(logFile, "    E_total                          = %.2lf kcal/mol\n\n", etotal);
        pr(logFile, "    E_intermolecular_atomic-affinity = %.2lf kcal/mol\n",   emap_total);
        pr(logFile, "    E_intermolecular_electrostatic   = %.2lf kcal/mol\n\n", elec_total);
#else
        pr(logFile, "Total %10.2f  %10.2f  %10.2f  %7.3f\n\n", (emap_total + elec_total), emap_total, elec_total, charge_total);
        pr(logFile, "    E_total                          = %.2f kcal/mol\n\n", etotal);
        pr(logFile, "    E_intermolecular_atomic-affinity = %.2f kcal/mol\n",   emap_total);
        pr(logFile, "    E_intermolecular_electrostatic   = %.2f kcal/mol\n\n", elec_total);
#endif
        printEnergies(einter, eintra, torsFreeEnergy, "epdb: USER    ", ligand_is_inhibitor, elec_total, emap_total);
        pr(logFile, "\n");
        for ( i=0; i<true_ligand_atoms; i++ ) {
            for (xyz = 0;  xyz < SPACE;  xyz++) {
                crdpdb[i][xyz] -= lig_center[xyz];
            }
        }
        unbound_internal_FE = eintra;
        (void) fflush(logFile);
        break;

/*____________________________________________________________________________*/

    case DPF_TERMINATION:
        /*
         *  termination_criterion energy 0.1
         *  termination_criterion evals 25000
         *  termination_criterion time 120 s
         */
        /*
        (void) sscanf( line, "%*s %d", &i );
        */
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


/* _________________________________________________________________________
**
** If in command-mode, set the command file-pointers to standard i/o,
** _________________________________________________________________________
*/
if (command_mode) {
    status = cmdmode( natom,jobStart,tms_jobStart,
                      xlo,ylo,zlo, xhi,yhi,zhi, inv_spacing,
                      map, e_internal, WallEnergy, vt, tlist, ntor,
                      Nnb, nonbondlist, atomstuff, crdpdb,
                      hostnm, type, charge, abs_charge, qsp_abs_charge, B_calcIntElec, q1q2,
                      atm_typ_str, torsFreeEnergy,
                      ligand_is_inhibitor, map_center,
                      ignore_inter,
                      B_include_1_4_interactions, scale_1_4, 
                      sol_fn, parameterArray, unbound_internal_FE);
    exit( status );  /* "command_mode" exits here... */
}

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

/* EOF */

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
