/*

 $Id: prClusterHist.h,v 1.5 2009/05/08 23:02:15 rhuey Exp $

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

#ifndef PRCLUSTERHIST
#define PRCLUSTERHIST

#include "constants.h"

void  prClusterHist(int   ncluster,
                    int   irunmax,
                    Real clus_rms_tol,
                    int   num_in_clu[MAX_RUNS],
                    int   cluster[MAX_RUNS][MAX_RUNS],
                    Real econf[MAX_RUNS],
                    Real clu_rms[MAX_RUNS][MAX_RUNS],
		    Real ref_rms[MAX_RUNS]);
#endif
