#ifndef EVALUATE_ENERGY
#define EVALUATE_ENERGY

#include <stdio.h>
#include "constants.h"
#include "trilinterp.h"
#include "eintcal.h"

float evaluate_energy( 
    float crd[MAX_ATOMS][SPACE],
    float charge[MAX_ATOMS],
    int   type[MAX_ATOMS],
    int   natom,
    float map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS],
    float inv_spacing,
    float xlo,
    float ylo,
    float zlo,
    int   nonbondlist[MAX_NONBONDS][2],
    float e_internal[NEINT][ATOM_MAPS][ATOM_MAPS],
    int   Nnb,
    Boole B_calcIntElec,
    float q1q2[MAX_NONBONDS],
    Boole B_isGaussTorCon,
    Boole B_isTorConstrained[MAX_TORS],
    State now,
    Boole B_ShowTorE,
    unsigned short US_TorE[MAX_TORS],
    unsigned short US_torProfile[MAX_TORS][NTORDIVS]);

/*
float evaluate_energy( 
    float crd[MAX_ATOMS][SPACE],
    float charge[MAX_ATOMS],
    int   type[MAX_ATOMS],
    int   natom,
    float map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS],
    float inv_spacing,
    float xlo,
    float ylo,
    float zlo,
/ *    float *Addr_eintra, * /
    int   nonbondlist[MAX_NONBONDS][2],
    float e_internal[NEINT][ATOM_MAPS][ATOM_MAPS],
    int   Nnb,
    Boole B_calcIntElec,
    float q1q2[MAX_NONBONDS],
    Boole B_isGaussTorCon,
    / * int   ntor, * /
    Boole B_isTorConstrained[MAX_TORS],
    State now, / * WAS: float torNow[MAX_TORS], * /
    Boole B_ShowTorE,
    unsigned short US_TorE[MAX_TORS],
    unsigned short US_torProfile[MAX_TORS][NTORDIVS]);
    */

#endif
