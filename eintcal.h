/*
 $Id: eintcal.h,v 1.6.6.1 2005/10/10 16:35:52 alther Exp $
*/

#ifndef EINTCAL_H
#define EINTCAL_H

#include "constants.h"
#include "structs.h"
#include "typedefs.h"

FloatOrDouble  eintcal(const int            nonbondlist[MAX_NONBONDS][MAX_NBDATA],
                       const EnergyTables   *ad_energy_tables,
                       const FloatOrDouble  tcoord[MAX_ATOMS][SPACE],
                       const int            Nnb,
                       const Boole          B_calcIntElec,
                       const FloatOrDouble  q1q2[MAX_NONBONDS],
                       const Boole          B_include_1_4_interactions,
                       const FloatOrDouble  scale_1_4,
                       const FloatOrDouble  qsp_abs_charge[MAX_ATOMS],
                       const ParameterEntry parameterArray[MAX_MAPS],
                       const FloatOrDouble  unbound_internal_FE
                      );

#endif        // EINTCAL_H
