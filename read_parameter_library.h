/*

 $Id: read_parameter_library.h,v 1.9 2009/06/10 00:09:09 rhuey Exp $

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

#ifndef _READ_PARAMETER_LIBRARY
#define _READ_PARAMETER_LIBRARY

#include "autocomm.h"

void read_parameter_library(
        char *FN_parameter_library,
        int outlev
        );

void setup_parameter_library(
        int outlev,
        const char * model_text,
        Unbound_Model unbound_model
        );

char * report_parameter_library();

void setup_distdepdiel( int outlev, 
                        EnergyTables *ptr_ad_energy_tables  // Holds vdw+Hb, desolvation & dielectric lookup tables
                      );


#endif
