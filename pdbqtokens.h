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
#define PDBQ_TORSDOF 14      /* PDBQ-file token for "TORSDOF". - torsional degrees of freedom */
