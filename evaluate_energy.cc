/*

 $Id: evaluate_energy.cc,v 1.6 2005/09/28 22:54:20 garrett Exp $

*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/* evaluate_energy.cc */

    #include "structs.h"
    #include "evaluate_energy.h"


extern FILE *logFile;

FloatOrDouble evaluate_energy( 

    FloatOrDouble crd[MAX_ATOMS][SPACE],
    FloatOrDouble charge[MAX_ATOMS],
    FloatOrDouble abs_charge[MAX_ATOMS],
    FloatOrDouble qsp_abs_charge[MAX_ATOMS],
    int   type[MAX_ATOMS],
    int   natom,
    FloatOrDouble map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS],
    int   nonbondlist[MAX_NONBONDS][MAX_NBDATA],

    EnergyTables *ptr_ad_energy_tables,
    
    int   Nnb,
    Boole B_calcIntElec,
    FloatOrDouble q1q2[MAX_NONBONDS],
    Boole B_isGaussTorCon,
    Boole B_isTorConstrained[MAX_TORS],
    State now,
    Boole B_ShowTorE,
    unsigned short US_TorE[MAX_TORS],
    unsigned short US_torProfile[MAX_TORS][NTORDIVS],
    const Boole         B_include_1_4_interactions,
    const FloatOrDouble scale_1_4,
    const ParameterEntry parameterArray[MAX_MAPS],
    const FloatOrDouble unbound_internal_FE,
    GridMapSetInfo *info

   )

{
    FloatOrDouble e = 0.;
    int   I_tor = 0;
    int   indx = 0;

    e = quicktrilinterp( crd, charge, abs_charge, type, natom, map, info);

    /* pr(logFile,"e(tril)=%10.2f,  ",e); / *###*/

    e += eintcal( nonbondlist, ptr_ad_energy_tables, crd, Nnb, B_calcIntElec, q1q2, B_include_1_4_interactions, scale_1_4, qsp_abs_charge, parameterArray, unbound_internal_FE);

    /* pr(logFile,"e(eintcal)=%10.2f\n",e); / *###*/

    /*
    ** FASTER, LESS ACCURATE  METHOD:
    ** (alternative to above line).
    **
    ** Only bother to calculate internal energy if
    ** trilinterp energy is low enough.
    **
    ** if ( (e = quicktrilinterp( crd, charge, abs_charge, type, natom,
    ** map, inv_spacing, xlo,ylo,zlo)) <
    ** ENERGY_CUTOFF) {
    **   e += (*Addr_eintra = eintcal( nonbondlist,ptr_ad_energy_tables,crd,Nnb,B_calcIntElec,q1q2, B_include_1_4_interactions, scale_1_4, abs_charge, parameterArray, unbound_internal_FE));
    ** }
    */

    if (B_isGaussTorCon) {
        for (I_tor = 0; I_tor <= now.ntor; I_tor++) {
            if (B_isTorConstrained[I_tor] == 1) {
                indx = Rad2Div( now.tor[I_tor] );
                if (B_ShowTorE) {
                    e += (FloatOrDouble)(US_TorE[I_tor] = US_torProfile[I_tor][indx]);
                } else {
                    e += (FloatOrDouble)US_torProfile[I_tor][indx];
                }
            }/*if*/
        }/*I_tor*/
    }/*if*/

    return e;
}
/* EOF */