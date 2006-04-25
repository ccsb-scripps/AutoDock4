/*

 $Id: evaluate_energy.cc,v 1.10 2006/04/25 22:32:09 garrett Exp $

*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/* evaluate_energy.cc */

    #include "structs.h"
    #include "evaluate_energy.h"


extern FILE *logFile;

Real evaluate_energy( 

    Real crd[MAX_ATOMS][SPACE],
    Real charge[MAX_ATOMS],
    Real abs_charge[MAX_ATOMS],
    Real qsp_abs_charge[MAX_ATOMS],
    int   type[MAX_ATOMS],
    int   natom,
    Real map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS],
    int   **nonbondlist,

    EnergyTables *ptr_ad_energy_tables,
    
    int   Nnb,
    Boole B_calcIntElec,
    Real q1q2[MAX_NONBONDS],
    Boole B_isGaussTorCon,
    Boole B_isTorConstrained[MAX_TORS],
    State now,
    Boole B_ShowTorE,
    unsigned short US_TorE[MAX_TORS],
    unsigned short US_torProfile[MAX_TORS][NTORDIVS],
    const Boole         B_include_1_4_interactions,
    const Real scale_1_4,
    const ParameterEntry parameterArray[MAX_MAPS],
    const Real unbound_internal_FE,
    GridMapSetInfo *info

   )

{
    Real e = 0.;
    int   I_tor = 0;
    int   indx = 0;

    e = trilinterp( crd, charge, abs_charge, type, natom, map, 
            info, ALL_ATOMS_INSIDE_GRID, NULL_IGNORE_INTERMOL, 
            NULL_ELEC, NULL_EVDW,
            NULL_ELEC_TOTAL, NULL_EVDW_TOTAL);

    e += eintcal( nonbondlist, ptr_ad_energy_tables, crd, Nnb, B_calcIntElec, q1q2, B_include_1_4_interactions, scale_1_4, qsp_abs_charge, parameterArray) - unbound_internal_FE;

    if (B_isGaussTorCon) {
        for (I_tor = 0; I_tor <= now.ntor; I_tor++) {
            if (B_isTorConstrained[I_tor] == 1) {
                indx = Rad2Div( now.tor[I_tor] );
                if (B_ShowTorE) {
                    e += (Real)(US_TorE[I_tor] = US_torProfile[I_tor][indx]);
                } else {
                    e += (Real)US_torProfile[I_tor][indx];
                }
            }/*if*/
        }/*I_tor*/
    }/*if*/

    return e;
}
/* EOF */
