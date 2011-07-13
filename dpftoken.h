/*

 $Id: dpftoken.h,v 1.33 2011/07/13 05:08:26 mp Exp $

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


/******************************************************************************
 *      Name: dpftoken.h                                                      *
 *  Function: Define tokens for parsing DPFs (docking parameter files)        *
 *Copyright (C) 2009 The Scripps Research Institute. All rights reserved.
 *----------------------------------------------------------------------------*
 *    Author: Garrett Matthew Morris, The Scripps Research Institute          *
 *      Date: 02/28/1995                                                      *
 *----------------------------------------------------------------------------*
 *    Inputs: none                                                            *
 *   Returns: nothing                                                         *
 *   Globals: none                                                            *
 *----------------------------------------------------------------------------*
 * Modification Record                                                        *
 * Date     Inits   Comments                                                  *
 * 09/06/95 RSH     GA/SW tokens added                                        *
 * 02/28/95 GMM     This header added                                         *
 ******************************************************************************/

#ifndef DPF_TOKENS
#define DPF_TOKENS

enum DpfTokens {
  DPF_UNKNOWN                  = -1,
  DPF_NULL                     =  0,
  DPF_COMMENT                  =  1,
  DPF_BLANK_LINE               =  2,
// 2 // (DPF_TYPES was removed, since the "types" command is no longer supported in AD4
  DPF_FLD                      =  3,
  DPF_MAP                      =  4,
  DPF_MOVE                     =  5,
  DPF_ABOUT                    =  6,
  DPF_TRAN0                    =  7,
  DPF_AXISANGLE0               =  8,
  DPF_NDIHE                    =  9,
  DPF_DIHE0                    = 10,
  DPF_TSTEP                    = 11,
  DPF_QSTEP                    = 12,
  DPF_DSTEP                    = 13,
  DPF_TRNRF                    = 14,
  DPF_QUARF                    = 15,
  DPF_DIHRF                    = 16,
  DPF_FLEX                     = 17,
  DPF_INTNBP_COEFFS            = 18,
  DPF_RT0	               = 19,
  DPF_RTRF                     = 20,
  DPF_RUNS                     = 21,
  DPF_CYCLES                   = 22,
  DPF_ACCS                     = 23,
  DPF_REJS                     = 24,
  DPF_SELECT                   = 25,
  DPF_OUTLEV                   = 26,
  DPF_RMSTOL                   = 27,
  DPF_TRJFRQ                   = 28,
  DPF_TRJBEG                   = 29,
  DPF_TRJEND                   = 30,
  DPF_TRJOUT                   = 31,
  DPF_TRJSEL                   = 32,
  DPF_EXTNRG                   = 33,
  DPF_NEWCRD                   = 34,
  DPF_CLUSTER                  = 35,
  DPF_CLUSALL                  = 36,
  DPF_RMSNOSYM                 = 37,
  DPF_SCHEDLIN                 = 38,
  DPF_RMSREF                   = 39,
  DPF_INTELEC                  = 40,
  DPF_SEED                     = 41,
  DPF_INTNBP_REQM_EPS          = 42,
  DPF_WATCH                    = 43,
  DPF_GAUSSTORCON              = 44,
  DPF_BARRIER                  = 45,
  DPF_SHOWTORPEN               = 46,
  DPF_HARDTORCON               = 47,
  DPF_E0MAX                    = 48,
  DPF_CHARMAP                  = 49,
  DPF_RAMP_VDW_REPULSION       = 50,
  DPF_SIMANNEAL	               = 51,
  DPF_GALS                     = 52,
  DPF_SET_GA                   = 53,
  DPF_SET_SW1                  = 54,
  DPF_SET_PSW1                 = 55,
  GA_pop_size                  = 56,
  GA_num_generations           = 57,
  GA_num_evals                 = 58,
  GA_window_size               = 59,
  GA_low                       = 60,
  GA_high	               = 61,
  GA_elitism                   = 62,
  GA_mutation_rate             = 63,
  GA_crossover_rate            = 64,
  GA_Cauchy_alpha              = 65,
  GA_Cauchy_beta               = 66,
  SW_max_its                   = 67,
  SW_max_succ                  = 68,
  SW_max_fail                  = 69,
  SW_rho                       = 70,
  SW_lb_rho                    = 71,
  LS_search_freq               = 72,
  DPF_LS                       = 73,
  DPF_GS                       = 74,
  DPF_ANALYSIS	               = 75,
  DPF_TORSDOF	               = 76,
  DPF_INVESTIGATE	       = 77,
  DPF_LIG_NOT_INHIB            = 78,
  DPF_TEMPL_ENERGY             = 79,
  DPF_HOLLOW_OUT               = 80,
  DPF_COLINY	               = 81,
  DPF_INCLUDE_1_4_INTERACTIONS = 82,
  DPF_PARAMETER_LIBRARY	       = 83,
  DPF_RECEPTOR_TYPES	       = 84,
  DPF_LIGAND_TYPES	       = 85,
  DPF_UNBOUND	               = 86,
  DPF_EPDB	               = 87,
  DPF_TERMINATION	       = 88,
  GA_CROSSOVER_MODE	       = 89,
  DPF_POPFILE                  = 90,
  DPF_SET_PATTERN              = 91,
  DPF_COMPUTE_UNBOUND_EXTENDED = 92,
  DPF_FLEXRES                  = 93,
  DPF_ELECMAP                  = 94,
  DPF_DESOLVMAP                = 95,
  DPF_UNBOUND_INTNBP_COEFFS    = 96,
  DPF_RMSATOMS                 = 97,
  DPF_CONFSAMPLER              = 98,
  DPF_REORIENT                 = 99,
  DPF_QUATERNION0              =100,
  DPF_COPYRIGHT                =101,
  DPF_WARRANTY                 =102,
  DPF_QUAT0	                   =103,
  DPF_PARAMETER_VERSION        =104,
  DPF_UNBOUND_MODEL            =105,
  PSW_TRANS_SCALE              =106,
  PSW_ROT_SCALE                =107,
  PSW_TORS_SCALE               =108,
  GA_PROPORTIONAL_SELECTION    =109,
  GA_TOURNAMENT_SELECTION      =110,
  GA_BOLTZMAN_SELECTION        =111,
  PSO_W_START	                   , // the initializers are unnecessary so...
  PSO_W_END	                   ,
  PSO_TVMAX                    ,
  PSO_QVMAX                    ,
  PSO_RVMAX                    ,
  PSO_C1 		       ,
  PSO_C2 		       ,
  PSO_NEIGHBORS_DYNAMIC	       ,
  PSO_RANDOM_BY_DIMENSION      ,
  PSO_ADAPTIVE_VELOCITY       ,
  PSO_K                       ,
  DPF_PSO_CONSTRICTION	      , //this is the no longer supported PSO implemented in autodock so far (5/2011)
  PSO_STAGE2CONSTRICTION	,
  PSO_INTERPOLATE_AS_SCALARS	,
  DPF_PSO_SSM		        ,
  DPF_PARSWARMOPT              ,  //this is the current PSO implemented in autodock (7/2011)
  GA_LINEAR_RANKING_SELECTION  ,
  DPF_RMSMODE                  ,
  DPF_SCALE_EINTERMOL          ,
  DPF_OUTPUT_POP_STATS         ,
  DPF_OUTPUT_RESNUM_AS         ,
  DPF_SMOOTH                   ,
};
#endif
/* EOF */
