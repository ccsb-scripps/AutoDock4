/*

 $Id: mapping.cc,v 1.3 2007/04/27 06:01:49 garrett Exp $

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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/********************************************************************
     These are the user defined functions that perform the
     mapping between Genotype and Phenotype and its inverse

				rsh 9/95
********************************************************************/
#include "support.h"

extern FILE *logFile;

//  This should be made more efficient.  As it is now, we (de facto) AlwaysEval!!!!
Phenotype Individual::mapping(void)
{

#ifdef DEBUG
   (void)fprintf(logFile, "mapping.cc/Phenotype Individual::mapping(void)\n");
#endif /* DEBUG */

   phenotyp.write(*genotyp.vread(0), 0);
   phenotyp.write(*genotyp.vread(1), 1);
   phenotyp.write(*genotyp.vread(2), 2);
   phenotyp.write(*genotyp.vread(3), 3);
   phenotyp.write(*genotyp.vread(4), 4);
   value(Normal_Eval);

   return(phenotyp);
}

Genotype Individual::inverse_mapping(void)
{

#ifdef DEBUG
   (void)fprintf(logFile, "mapping.cc/Genotype Individual::inverse_mapping(void)\n");
#endif /* DEBUG */

   genotyp.write(*phenotyp.vread(0), 0);
   genotyp.write(*phenotyp.vread(1), 1);
   genotyp.write(*phenotyp.vread(2), 2);
   genotyp.write(*phenotyp.vread(3), 3);
   genotyp.write(*phenotyp.vread(4), 4);

   return(genotyp);
}
