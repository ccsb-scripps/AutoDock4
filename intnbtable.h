/*

 $Id: intnbtable.h,v 1.8 2007/04/27 06:01:49 garrett Exp $

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

#ifndef INTNBTABLE
#define INTNBTABLE
#include "constants.h"
#include "timesys.h"
#include "structs.h"

void intnbtable(Boole *P_B_havenbp,
                int   a1,
                int   a2,
                GridMapSetInfo *info,
                Real cA,
                Real cB,
                int   xA,
                int   xB,
                double coeff_desolv,
                double sigma,
                EnergyTables *ad_tables,
                Boole B_is_unbound_calculation );
#endif
