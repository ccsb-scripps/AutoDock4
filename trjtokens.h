/* trjtokens.h */

/******************************************************************************
 *      Name: trjtokens.h                                                     *
 *  Function: Defines tokens for parsing AUTODOCK trajectory files.           *
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
 ******************************************************************************/


#define TRJ_NULL     0      /* Trajectory-file token for '\n' or '\0'. */
#define TRJ_NTOR     1      /* Trajectory-file token for "ntorsions" */
#define TRJ_RUN      2      /* Trajectory-file token for "run". */
#define TRJ_CYCLE    3      /* Trajectory-file token for "cycle". */
#define TRJ_STATE    4      /* Trajectory-file token for "state". */
#define TRJ_TEMP     5      /* Trajectory-file token for "temp". */
