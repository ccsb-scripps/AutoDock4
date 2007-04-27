/*

 $Id: dpftoken.h,v 1.16 2007/04/27 06:01:48 garrett Exp $

 AutoDock 

 Copyright (C) 1989-2007,  Garrett M. Morris, David S. Goodsell, Ruth Huey, Arthur J. Olson, 
 All Rights Reserved.

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
 * Copyright: (C) Garrett Matthew Morris, TSRI                                *
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

#define	DPF_		-1
#define	DPF_NULL	 0
#define	DPF_COMMENT	 1
// 2 // (DPF_TYPES was removed, since the "types" command is no longer supported in AD4
#define	DPF_FLD		 3
#define	DPF_MAP		 4
#define	DPF_MOVE	 5
#define	DPF_ABOUT	 6
#define	DPF_TRAN0	 7
#define	DPF_AXISANGLE0	 8
#define	DPF_NDIHE	 9
#define	DPF_DIHE0	10
#define	DPF_TSTEP	11
#define	DPF_QSTEP	12
#define	DPF_DSTEP	13
#define	DPF_TRNRF	14
#define	DPF_QUARF	15
#define	DPF_DIHRF	16
#define	DPF_FLEX	17
#define	DPF_INTNBP_COEFFS	18
#define	DPF_RT0		19
#define	DPF_RTRF	20
#define	DPF_RUNS	21
#define	DPF_CYCLES	22
#define	DPF_ACCS	23
#define	DPF_REJS	24
#define	DPF_SELECT	25
#define	DPF_OUTLEV	26
#define	DPF_RMSTOL	27
#define	DPF_TRJFRQ	28
#define	DPF_TRJBEG	29
#define	DPF_TRJEND	30
#define	DPF_TRJOUT	31
#define	DPF_TRJSEL	32
#define	DPF_EXTNRG	33
#define	DPF_NEWCRD	34
#define	DPF_CLUSTER	35
#define	DPF_CLUSALL	36
#define	DPF_RMSNOSYM	37
#define	DPF_SCHEDLIN	38
#define	DPF_RMSREF	39
#define	DPF_INTELEC	40
#define	DPF_SEED	41
#define	DPF_INTNBP_REQM_EPS	42
#define	DPF_WATCH	43
#define	DPF_GAUSSTORCON	44
#define	DPF_BARRIER	45
#define	DPF_SHOWTORPEN	46
#define	DPF_HARDTORCON	47
#define	DPF_E0MAX	48
#define	DPF_CHARMAP	49
#define	DPF_RAMP_VDW_REPULSION 50
#define	DPF_SIMANNEAL	51
#define	DPF_GALS	52
#define DPF_SET_GA	53
#define DPF_SET_SW1	54
#define DPF_SET_PSW1	55
#define GA_pop_size	56
#define GA_num_generations	57
#define GA_num_evals	58
#define GA_window_size	59
#define GA_low		60
#define GA_high		61
#define GA_elitism	62
#define GA_mutation_rate	63
#define GA_crossover_rate	64
#define GA_Cauchy_alpha	65
#define GA_Cauchy_beta	66
#define SW_max_its	67
#define SW_max_succ	68
#define SW_max_fail	69
#define SW_rho		70
#define SW_lb_rho	71
#define LS_search_freq	72
#define DPF_LS		73
#define DPF_GS		74
#define	DPF_ANALYSIS	75
#define	DPF_TORSDOF	76
#define	DPF_INVESTIGATE	77
#define DPF_LIG_NOT_INHIB 78
#define DPF_TEMPL_ENERGY 79
#define DPF_HOLLOW_OUT 80
#define DPF_COLINY	81
#define DPF_INCLUDE_1_4_INTERACTIONS	82
#define DPF_PARAMETER_LIBRARY	83
#define DPF_RECEPTOR_TYPES	    84
#define DPF_LIGAND_TYPES	    85
#define DPF_UNBOUND	    86
#define DPF_EPDB	    87
#define DPF_TERMINATION	    88
#define GA_CROSSOVER_MODE	    89
#define DPF_POPFILE 90
#define DPF_SET_PATTERN 91
#define DPF_COMPUTE_UNBOUND_EXTENDED	    92
#define DPF_FLEXRES 93
#define DPF_ELECMAP 94
#define DPF_DESOLVMAP 95
#define	DPF_UNBOUND_INTNBP_COEFFS	96
#define	DPF_RMSATOMS	97
#define DPF_CONFSAMPLER	98
#define DPF_REORIENT	99
#define	DPF_QUATERNION0	 100
#define DPF_COPYRIGHT 101
#define DPF_WARRANTY 102

#endif
/* EOF */
