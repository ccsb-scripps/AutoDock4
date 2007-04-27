/*

 $Id: coliny.h,v 1.4 2007/04/27 06:01:48 garrett Exp $

 AutoDock 

 Copyright (C) 1989-2007,  William Hart, Garrett M. Morris, David S. Goodsell, Ruth Huey, Arthur J. Olson, 
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

//
// coliny.h
//
// Interface to Coliny optimizers
//

#ifndef __coliny_h
#define __coliny_h

#if defined(USING_COLINY)

#include <vector>

//
// Initialize the 'algname' coliny optimizer over 'domain'
//
void coliny_init(char* algname, char* domain);

//
// Perform minimization with a given seed and initial point. Return
// summary statistics
//
void coliny_minimize(int seed, std::vector<double>& initpt,
				std::vector<double>& finalpt,
				int& neval, int& niters);

#endif

#endif
