/********************************************************************
    The header file for the eval class

                                rsh 09/06/95
********************************************************************/
#ifndef _EVAL_H
#define _EVAL_H

#include "autocomm.h"
#include "constants.h"
#include "grid.h"
#include "rep.h"
#include "structs.h"
#include "typedefs.h"


#if defined(USING_COLINY)
void make_state_from_rep(double *x, int n, State *now);
#endif

void make_state_from_rep(Representation **rep, State *stateNow);

class Eval
{
   private:
      UnsignedFourByteLong num_evals;
      int natom, Nnb;
      GridMapSetInfo *info;
      FloatOrDouble eval_elec[MAX_ATOMS]; // gmm added 21-Jan-1998, for writePDBQState
      FloatOrDouble eval_emap[MAX_ATOMS]; // gmm added 21-Jan-1998, for writePDBQState
      Boole B_calcIntElec, B_isGaussTorCon, B_ShowTorE;
      State stateNow;
      unsigned short *US_TorE, (*US_torProfile)[NTORDIVS];
      int *type, (*nonbondlist)[MAX_NBDATA], (*tlist)[MAX_ATOMS];
      FloatOrDouble *q1q2, *charge, *abs_charge, *qsp_abs_charge;
      FloatOrDouble (*crd)[SPACE], (*vt)[SPACE], (*crdpdb)[SPACE];
      EnergyTables *ptr_ad_energy_tables;
      FloatOrDouble (*map)[MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS];
      Boole *B_isTorConstrained;
      Molecule mol;
      Boole B_template; // Use the template-docking scoring function if true 15-jan-2001
      // FloatOrDouble template_energy[MAX_ATOMS]; // atomic template values
      // FloatOrDouble template_stddev[MAX_ATOMS]; // atomic template values
      FloatOrDouble *template_energy; // atomic template values
      FloatOrDouble *template_stddev; // atomic template values
      int ignore_inter[MAX_ATOMS]; // gmm 2002-05-21, for CA, CB in flexible sidechains
      Boole         B_include_1_4_interactions; // gmm 2005-01-8, for scaling 1-4 nonbonds
      FloatOrDouble scale_1_4;                  // gmm 2005-01-8, for scaling 1-4 nonbonds
      ParameterEntry *parameterArray;
      FloatOrDouble  unbound_internal_FE;

   public:
      Eval(void);
      void setup(FloatOrDouble init_crd[MAX_ATOMS][SPACE],
          FloatOrDouble  init_charge[MAX_ATOMS],
          FloatOrDouble  init_abs_charge[MAX_ATOMS],
          FloatOrDouble  init_qsp_abs_charge[MAX_ATOMS],
          int            init_type[MAX_ATOMS], int init_natom,
          FloatOrDouble  init_map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS],

          FloatOrDouble  init_elec[MAX_ATOMS], // gmm added 21-Jan-1998, for writePDBQState
          FloatOrDouble  init_emap[MAX_ATOMS], // gmm added 21-Jan-1998, for writePDBQState

          int            init_nonbondlist[MAX_NONBONDS][MAX_NBDATA],
          EnergyTables   *init_ptr_ad_energy_tables,
          int init_Nnb,
          Boole          init_B_calcIntElec, FloatOrDouble init_q1q2[MAX_NONBONDS],
          Boole          init_B_isGaussTorCon, Boole init_B_isTorConstrained[MAX_TORS],
          Boole          init_B_ShowTorE, unsigned short init_US_TorE[MAX_TORS],
          unsigned short init_US_torProfile[MAX_TORS][NTORDIVS],
          FloatOrDouble  init_vt[MAX_TORS][SPACE], int init_tlist[MAX_TORS][MAX_ATOMS],
          FloatOrDouble  init_crdpdb[MAX_ATOMS][SPACE], State stateInit, Molecule molInit,
          Boole          init_B_template,
          FloatOrDouble  init_template_energy[MAX_ATOMS],
          FloatOrDouble  init_template_stddev[MAX_ATOMS],
          int            init_ignore_inter[MAX_ATOMS],
          Boole          init_B_include_1_4_interactions, // gmm 2005-01-8, for scaling 1-4 nonbonds
          FloatOrDouble  init_scale_1_4,                   // gmm 2005-01-8, for scaling 1-4 nonbonds
          ParameterEntry init_parameterArray[MAX_MAPS],
          FloatOrDouble  init_unbound_internal_FE,
          GridMapSetInfo *init_info
          );

      double operator()(Representation **);
#if defined(USING_COLINY)
      double operator()(double*, int);
#endif
      double eval();			// WEH - a basic change that facilitates the use of Coliny
      UnsignedFourByteLong evals(void);
      void reset(void);
      int write(FILE *out_file, Representation **rep);
};

inline Eval::Eval(void)
: num_evals(0)
{
}

inline void Eval::setup(FloatOrDouble init_crd[MAX_ATOMS][SPACE],
                        FloatOrDouble init_charge[MAX_ATOMS],
                        FloatOrDouble init_abs_charge[MAX_ATOMS],
                        FloatOrDouble init_qsp_abs_charge[MAX_ATOMS],
                        int init_type[MAX_ATOMS],
                        int init_natom,
                        FloatOrDouble init_map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS],

                        FloatOrDouble init_elec[MAX_ATOMS], // gmm added 21-Jan-1998, for writePDBQState
                        FloatOrDouble init_emap[MAX_ATOMS], // gmm added 21-Jan-1998, for writePDBQState
                        int init_nonbondlist[MAX_NONBONDS][MAX_NBDATA],
                        EnergyTables   *init_ptr_ad_energy_tables,
                        int init_Nnb,
                        Boole init_B_calcIntElec, FloatOrDouble init_q1q2[MAX_NONBONDS],
                        Boole init_B_isGaussTorCon,
                        Boole init_B_isTorConstrained[MAX_TORS],
                        Boole init_B_ShowTorE,
                        unsigned short init_US_TorE[MAX_TORS],
                        unsigned short init_US_torProfile[MAX_TORS][NTORDIVS],
                        FloatOrDouble init_vt[MAX_TORS][SPACE],
                        int init_tlist[MAX_TORS][MAX_ATOMS],
                        FloatOrDouble init_crdpdb[MAX_ATOMS][SPACE],
                        State stateInit,
                        Molecule molInit,

                        Boole init_B_template,
                        FloatOrDouble init_template_energy[MAX_ATOMS],
                        FloatOrDouble init_template_stddev[MAX_ATOMS],

                        int   init_ignore_inter[MAX_ATOMS],

                        Boole         init_B_include_1_4_interactions,
                        FloatOrDouble init_scale_1_4,

                        ParameterEntry init_parameterArray[MAX_MAPS],

                        FloatOrDouble init_unbound_internal_FE,
                        GridMapSetInfo *init_info
                       )

{
    register int i;

    crd = init_crd;
    charge = init_charge;
    abs_charge = init_abs_charge;
    qsp_abs_charge = init_qsp_abs_charge;
    type = init_type;
    natom = init_natom;
    map = init_map;

    nonbondlist = init_nonbondlist;
    ptr_ad_energy_tables = init_ptr_ad_energy_tables;
    Nnb = init_Nnb;
    B_calcIntElec = init_B_calcIntElec;
    q1q2 = init_q1q2;
    B_isGaussTorCon = init_B_isGaussTorCon;
    B_isTorConstrained = init_B_isTorConstrained;
    B_ShowTorE = init_B_ShowTorE;
    US_TorE = init_US_TorE;
    US_torProfile = init_US_torProfile;
    vt = init_vt;
    tlist = init_tlist;
    crdpdb = init_crdpdb;
    stateNow = stateInit;
    num_evals = 0;
    for (i=0; i<MAX_ATOMS; i++) {
       init_elec[i] = init_emap[i] = 0.0;
       ignore_inter[i] = init_ignore_inter[i];
    }
    mol = molInit;
    B_template = init_B_template;
    template_energy = init_template_energy;
    template_stddev = init_template_stddev;

    B_include_1_4_interactions = init_B_include_1_4_interactions;
    scale_1_4 = init_scale_1_4;

    parameterArray = init_parameterArray;

    unbound_internal_FE = init_unbound_internal_FE;

    info = init_info;
}

inline UnsignedFourByteLong Eval::evals(void)
{
   return(num_evals);
}

inline void Eval::reset(void)
{
   num_evals = 0;
}

#endif
