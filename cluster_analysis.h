/*

 $Id: cluster_analysis.h,v 1.4 2007/04/27 06:01:48 garrett Exp $

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

#ifndef CLUSTER_ANALYSIS
#define CLUSTER_ANALYSIS
#include "constants.h"
#include "getrms.h"


int  cluster_analysis( Real clus_rms_tol, 
                       int   cluster[MAX_RUNS][MAX_RUNS], 
                       int   num_in_clus[MAX_RUNS], 
                       int   isort[MAX_RUNS], 
                       int   nconf, 
                       int   natom, 
                       int   type[MAX_ATOMS], 
                       Real crd[MAX_RUNS][MAX_ATOMS][SPACE], 
                       Real crdpdb[MAX_ATOMS][SPACE], 
                       Real sml_center[SPACE], 
                       Real clu_rms[MAX_RUNS][MAX_RUNS], 
                       Boole B_symmetry_flag,
                       Real ref_crds[MAX_ATOMS][SPACE],
                       int   ref_natoms,
                       Real ref_rms[MAX_RUNS]);
#endif
