/*

 $Id: autoglobal.h,v 1.12 2007/04/27 06:01:47 garrett Exp $

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

#ifndef _AUTOGLOBAL
#define _AUTOGLOBAL

#include <sys/types.h>
#include <string.h>
#include <stdio.h>

#include "structs.h"

/******************************************************************************/
/*      Name: autoglobal.h                                                    */
/*  Function: Global variables for Autodock modules.                          */
/* Copyright: (C) TSRI                                                        */
/*----------------------------------------------------------------------------*/
/*    Author: Garrett M. Morris                                               */
/*                                                                            */
/*            e-mail: garrett@scripps.edu				      */
/*                                                                            */
/*            The Scripps Research Institute                                  */
/*            Department of Molecular Biology, MB5                            */
/*            10666 North Torrey Pines Road                                   */
/*            La Jolla, CA 92037.                                             */
/*                                                                            */
/*      Date: 03/18/93                                                        */
/*----------------------------------------------------------------------------*/
/*    Inputs: None.                                                           */
/*   Returns: Nothing.                                                        */
/*   Globals: programname, AutoGridHelp, AutoDockHelp, command_mode,          */
/*            command_in_fp, command_out_fp, GPF, logFile.                    */
/*----------------------------------------------------------------------------*/
/* Modification Record                                                        */
/* Date     Inits   Comments                                                  */
/* 03/18/93 GMM     Created.                                                  */
/******************************************************************************/


/*----------------------------------------------------------------------------*/
/* Globals,                                                                   */
/*----------------------------------------------------------------------------*/

char    *programname;
char    AutoDockHelp[] = "           -p parameter_filename\n           -l log_filename\n           -o (Use old PDBQ format, charge q in columns 55-61)\n           -k (Keep original residue numbers)\n           -i (Ignore header-checking)\n           -t (Parse the PDBQ file to check torsions, then stop.)\n           -c < command_file (Command mode, by file)\n           -c | control_program (Command mode, by control_program)\n\n";

char    AutoGridHelp[] = "-p parameter_filename\n-l log_filename\n-o (old PDBQ format)\n-d (increment debug level)\n-u (display this message)\n";

char    dock_param_fn[MAX_CHARS];
char    grid_param_fn[MAX_CHARS];

int     command_mode = FALSE;
int     debug = 0;
int	    ElecMap = 0;
int	    DesolvMap = 0;
int     ignore_errors = FALSE;
int     keepresnum = 1;
int     parse_tors_mode = FALSE;
int	    true_ligand_atoms = 0;
int     write_stateFile = FALSE;
// For energy breakdown of non-bonded interactions
int     Nnb_array[3] = {0};    // number of nonbonds in the ligand, intermolecular and receptor groups

Real	idct = 1.0;
// For energy breakdown of non-bonded interactions
Real    nb_group_energy[3] = {0.0};  // total energy of each nonbond group (intra-ligand, inter, and intra-receptor)

FILE    *command_in_fp;
FILE    *command_out_fp;
FILE    *parFile;
FILE    *GPF;
FILE    *logFile;
FILE    *stateFile;

Linear_FE_Model AD3;
Linear_FE_Model AD4_wrt_3;
Linear_FE_Model AD4;

#endif /*_AUTOGLOBAL*/
/* EOF */
