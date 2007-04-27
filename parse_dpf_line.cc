/*

 $Id: parse_dpf_line.cc,v 1.18 2007/04/27 06:01:50 garrett Exp $

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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include "parse_dpf_line.h"


#define NUM_LEXEMES_AUTODOCK 113 // this is the length of the tokentable of AutoDock-related lexemes 
#define NUM_LEXEMES_COLINY     1 // this is the length of the tokentable of Coliny-related lexemes 


int parse_dpf_line( char line[LINE_LEN] )

/******************************************************************************/
/*      Name: parse_dpf_line                                                  */
/*  Function: Parse the docking parameter file line                           */
/* Copyright: (C) 1994, TSRI                                                  */
/*----------------------------------------------------------------------------*/
/*    Author: Garrett Morris, The Scripps Research Institute                  */
/*      Date: 19/05/94                                                        */
/*----------------------------------------------------------------------------*/
/*    Inputs: line                                                            */
/*   Returns: integer token describing the keyword found.                     */
/*   Globals: none.                                                           */
/*----------------------------------------------------------------------------*/
/* Modification Record                                                        */
/* Date     Inits   Comments                                                  */
/* 06/09/95 RSH     Changed to an array implementation                        */
/* 19/05/94 GMM     Entered code.                                             */
/******************************************************************************/

{
    int j, i, token = DPF_;               /* return -1 if nothing is recognized. */
    char c[LINE_LEN];

    // tokentablesize should be set to the length of the tokentable,
    //
#if defined(USING_COLINY)
    const int tokentablesize = NUM_LEXEMES_AUTODOCK + NUM_LEXEMES_COLINY;
#else
    const int tokentablesize = NUM_LEXEMES_AUTODOCK;
#endif

    const struct {
       char *lexeme;
       int tokenvalue;
    } tokentable[] = {{"ligand", DPF_MOVE},  // 1
                      {"fld", DPF_FLD}, // 2
                      {"map", DPF_MAP}, // 3
                      {"move", DPF_MOVE}, // 4
                      {"about", DPF_ABOUT}, // 5
                      {"tran0", DPF_TRAN0}, // 6
                      {"quat0", DPF_AXISANGLE0}, // 7
                      {"ndihe", DPF_NDIHE}, // 8
                      {"dihe0", DPF_DIHE0}, // 9
                      {"torsdof", DPF_TORSDOF}, // 10
                      {"tstep", DPF_TSTEP}, // 11
                      {"qstep", DPF_QSTEP}, // 12
                      {"dstep", DPF_DSTEP}, // 13
                      {"trnrf", DPF_TRNRF}, // 14
                      {"quarf", DPF_QUARF}, // 15
                      {"dihrf", DPF_DIHRF}, // 16
                      {"flex", DPF_FLEX}, // 17
                      {"intnbp_coeffs", DPF_INTNBP_COEFFS}, // 18
                      {"rt0", DPF_RT0}, // 19
                      {"rtrf", DPF_RTRF}, // 20
                      {"runs", DPF_RUNS}, // 21
                      {"cycles", DPF_CYCLES}, // 22
                      {"accs", DPF_ACCS}, // 23
                      {"rejs", DPF_REJS}, // 24
                      {"select", DPF_SELECT}, // 25
                      {"outlev", DPF_OUTLEV}, // 26
                      {"rmstol", DPF_RMSTOL}, // 27
                      {"trjfrq", DPF_TRJFRQ}, // 28
                      {"trjbeg", DPF_TRJBEG}, // 29
                      {"trjend", DPF_TRJEND}, // 30
                      {"trjout", DPF_TRJOUT}, // 31
                      {"trjsel", DPF_TRJSEL}, // 32
                      {"extnrg", DPF_EXTNRG}, // 33
                      {"newcrd", DPF_NEWCRD}, // 34
                      {"cluster", DPF_CLUSTER}, // 35
                      {"write_all", DPF_CLUSALL}, // 36
                      {"write_all_cluster_members", DPF_CLUSALL}, // 37
                      {"charmap", DPF_CHARMAP}, // 38
                      {"rmsnosym", DPF_RMSNOSYM}, // 39
                      {"rmsref", DPF_RMSREF}, // 40
                      {"watch", DPF_WATCH}, // 41
                      {"linear_schedule", DPF_SCHEDLIN}, // 42
                      {"schedule_linear", DPF_SCHEDLIN}, // 43
                      {"linsched", DPF_SCHEDLIN}, // 44
                      {"schedlin", DPF_SCHEDLIN}, // 45
                      {"intelec", DPF_INTELEC}, // 46
                      {"seed", DPF_SEED}, // 47
                      {"e0max", DPF_E0MAX}, // 48
                      {"simanneal", DPF_SIMANNEAL}, // 49
                      {"hardtorcon", DPF_HARDTORCON}, // 50
                      {"intnbp_r_eps", DPF_INTNBP_REQM_EPS}, // 51
                      {"gausstorcon", DPF_GAUSSTORCON}, // 52
                      {"barrier", DPF_BARRIER}, // 53
                      {"showtorpen", DPF_SHOWTORPEN}, // 54
                      {"ga_run", DPF_GALS}, // 55
                      {"gals_run", DPF_GALS}, // 56
                      {"do_gals", DPF_GALS}, // 57
                      {"set_ga", DPF_SET_GA}, // 58
                      {"set_sw1", DPF_SET_SW1}, // 59
                      {"set_psw1", DPF_SET_PSW1}, // 60
                      {"analysis", DPF_ANALYSIS}, // 61
                      {"ga_pop_size", GA_pop_size}, // 62
                      {"ga_num_generations", GA_num_generations}, // 63
                      {"ga_num_evals", GA_num_evals}, // 64
                      {"ga_window_size", GA_window_size}, // 65
                      {"ga_low", GA_low}, // 66
                      {"ga_high", GA_high}, // 67
                      {"ga_elitism", GA_elitism}, // 68
                      {"ga_mutation_rate", GA_mutation_rate}, // 69
                      {"ga_crossover_rate", GA_crossover_rate}, // 70
                      {"ga_cauchy_alpha", GA_Cauchy_alpha}, // 71
                      {"ga_cauchy_beta", GA_Cauchy_beta}, // 72
                      {"sw_max_its", SW_max_its}, // 73
                      {"sw_max_succ", SW_max_succ}, // 74
                      {"sw_max_fail", SW_max_fail}, // 75
                      {"sw_rho", SW_rho}, // 76
                      {"sw_lb_rho", SW_lb_rho}, // 77
                      {"do_local_only", DPF_LS}, // 78
                      {"ls_run", DPF_LS}, // 79
                      {"do_global_only", DPF_GS}, // 80
                      {"ga_only_run", DPF_GS}, // 81
                      {"ls_search_freq", LS_search_freq}, // 82
                      {"bin_energies_by_rmsd", DPF_INVESTIGATE}, // 83
                      {"investigate", DPF_INVESTIGATE}, // 84
              {"ligand_is_not_inhibitor", DPF_LIG_NOT_INHIB}, // 85
              {"template", DPF_TEMPL_ENERGY}, // 86
              {"template_energy_file", DPF_TEMPL_ENERGY}, // 87
              {"include_1_4_interactions", DPF_INCLUDE_1_4_INTERACTIONS}, // 88
              {"parameter_library", DPF_PARAMETER_LIBRARY}, // 89
              {"parameter_file", DPF_PARAMETER_LIBRARY} // 90
              , {"receptor_types", DPF_RECEPTOR_TYPES}  // 91
              , {"ligand_types", DPF_LIGAND_TYPES}      // 92
              , {"unbound", DPF_UNBOUND}      // 93
              , {"epdb", DPF_EPDB}      // 94
              , {"ga_termination_criterion", DPF_TERMINATION}      // 95
              , {"ga_termination", DPF_TERMINATION}      // 96
              , {"ga_crossover_mode", GA_CROSSOVER_MODE}      // 97
              , {"output_pop_file", DPF_POPFILE}      // 98
              , {"set_pattern", DPF_SET_PATTERN}      // 99
              , {"compute_unbound_extended", DPF_COMPUTE_UNBOUND_EXTENDED} // 100
              , {"set_unbound_energy", DPF_UNBOUND}      // 101
              , {"flexible_residues", DPF_FLEXRES} // 102
              , {"flexres", DPF_FLEXRES} // 103
              , {"elecmap", DPF_ELECMAP} // 104
              , {"desolvmap", DPF_DESOLVMAP} // 105
              , {"unbound_intnbp_coeffs", DPF_UNBOUND_INTNBP_COEFFS} // 106
              , {"rmsatoms", DPF_RMSATOMS} // 107
              , {"confsampler", DPF_CONFSAMPLER} // 108
              , {"reorient", DPF_REORIENT} // 109
              , {"axisangle0", DPF_AXISANGLE0} // 110
              , {"quaternion0", DPF_QUATERNION0} // 111
              , {"copyright", DPF_COPYRIGHT} // 112
              , {"warranty", DPF_WARRANTY} // 113
			   // Remember to define NUM_LEXEMES_AUTODOCK earlier

#if defined(USING_COLINY)
              , {"coliny", DPF_COLINY}  // 1 
               // Remember to define NUM_LEXEMES_COLINY earlier
#endif
              };

    c[0] = '\0';
    for (j=0; ((line[j]!='\0')&&(line[j]!=' ')&&(line[j]!='\t')&&(line[j]!='\n')); j++) {
        /*  Ignore case */
        c[j] = (char)tolower((int)line[j]);
        /*(void)fprintf(stderr,"%c",c[j]);*/
    }
    /*(void)fprintf(stderr,"/n,j = %d\n",j);*/

    /*  Recognize one character tokens  */

    if ((c[0]=='\n') || (c[0]=='\0')) {
        token = DPF_NULL;
    } else if (c[0]=='#') {
        token = DPF_COMMENT;
    }

    /*  Recognize token strings  */

    for (i=0;  (i < tokentablesize) && (token == DPF_);  i++) {
        /*(void)fprintf(stderr,"i = %d, tokentable[i].lexeme = %s, tokentable[i].value = %d, c = %s\n",i,tokentable[i].lexeme,tokentable[i].tokenvalue,c);*/
        if (strncasecmp(tokentable[i].lexeme, c, j) == 0) {
            token = tokentable[i].tokenvalue;
        }
    }
    return(token);
}
/* EOF */
