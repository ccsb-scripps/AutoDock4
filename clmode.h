/*

 $Id: clmode.h,v 1.7 2009/05/08 23:02:11 rhuey Exp $

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


#ifndef CLMODE
#define CLMODE
#include "constants.h"
#include "strindex.h"
#include "readPDBQT.h"
#include "get_atom_type.h"
#include "getpdbcrds.h"
#include "sort_enrg.h"
#include "cluster_analysis.h"
#include "prClusterHist.h"
#include "bestpdb.h"
#include "success.h"
#include "qmultiply.h"
#include "openfile.h"
void  clmode( int   num_atm_maps, 
              Real clus_rms_tol, 
              char  hostnm[MAX_CHARS], 
              Clock jobStart,
              struct tms tms_jobStart, 
              Boole B_write_all_clusmem, 
              char  *clusFN, 
              Real crdpdb[MAX_ATOMS][SPACE], 
              Real sml_center[SPACE], 
              Boole B_symmetry_flag,
              char  *rms_ref_crds );
#endif
