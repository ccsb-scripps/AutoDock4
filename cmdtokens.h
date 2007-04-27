/*

 $Id: cmdtokens.h,v 1.2 2007/04/27 06:01:48 garrett Exp $

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
 *      Name: cmdtokens.h                                                     *
 *  Function: Define the tokens for parsing commands in AUTODOCK command mode.*
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
 * 09/06/95 RSH     GA-LS constants added                                     *
 * 02/28/95 GMM     This header added                                         *
 ******************************************************************************/


#define COM_NULL     0      /* Command mode token for '\n' or '\0'. */
#define COM_STOP     1      /* Command mode token for "STOP". */
#define COM_EVAL     2      /* Command mode token for "EVALUATE". */
#define COM_OUTC     3      /* Command mode token for "OUTCOORDS". */
#define COM_OUTE     4      /* Command mode token for "OUTENERGY_ATOMIC". */
#define COM_TRJ      5      /* Command mode token for "TRAJECTORY". */
#define COM_EPDB     6      /* Command mode token for "ENERGY_PDB". */
