/* cmdtokens.h */

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
