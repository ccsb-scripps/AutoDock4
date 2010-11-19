/*

 $Id: PDBQT_tokens.h,v 1.2.2.1 2010/11/19 20:09:29 rhuey Exp $

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

/* pdbqtokens.h */

/******************************************************************************
 *      Name: pdbqtokens.h                                                    *
 *  Function: Defines the tokens for PDBQ parsing.                            *
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
 * 02/28/95 GMM     This header added                                         *
 * 09/27/07 Huameng  add token for multiple ligands							  *
 ******************************************************************************/


#define PDBQ_UNRECOGNIZED    (-1)      /* PDBQ-file token for unrecognized symbols. */
#define PDBQ_NULL    0      /* PDBQ-file token for '\n' or '\0'. */
#define PDBQ_ROOT    1      /* PDBQ-file token for "ROOT" */
#define PDBQ_ENDROOT 2      /* PDBQ-file token for "ENDROOT". */
#define PDBQ_ATOM    3      /* PDBQ-file token for "ATOM". */
#define PDBQ_HETATM  4      /* PDBQ-file token for "HETATM". */
#define PDBQ_TORS    5      /* PDBQ-file token for "TORS". */
#define PDBQ_BRANCH  6      /* PDBQ-file token for "BRANCH". */
#define PDBQ_ENDBRANCH    7 /* PDBQ-file token for "ENDBRANCH". */
#define PDBQ_ENDTORS 8      /* PDBQ-file token for "ENDTORS". */
#define PDBQ_REMARK  9      /* PDBQ-file token for "REMARK". */
#define PDBQ_CONSTRAINT 10  /* PDBQ-file token for "CONSTRAINT". */
#define PDBQ_BEGIN_RES 11   /* PDBQ-file token for "BEGIN_RES". */
#define PDBQ_END_RES 12     /* PDBQ-file token for "END_RES". */
#define PDBQ_TORSDOF 14     /* PDBQ-file token for "TORSDOF" - torsional degrees of freedom. */
#define PDBQ_CONECT  15     /* PDBQ-file token for "CONECT" - PDB connectivity record. */

#define PDBQ_BEGIN_LIG 16   /* PDBQ-file token for "BEGIN_LIG (ligand)". */
#define PDBQ_END_LIG   17   /* PDBQ-file token for "END_LIG (ligand)". */
