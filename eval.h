/*

 $Id: eval.h,v 1.22 2009/10/02 22:12:38 rhuey Exp $

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

/********************************************************************
    The header file for the eval class

                                rsh 09/06/95
********************************************************************/
#ifndef _EVAL_H
#define _EVAL_H

#include <stdio.h>
#include "structs.h"
#include "rep.h"
#include "qmultiply.h"
#include "cnv_state_to_coords.h"
#include "trilinterp.h"
#include "eintcal.h"
#include "energy.h"

#ifdef DEBUG
extern FILE *logFile;
#endif

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
      Real eval_elec[MAX_ATOMS]; // gmm added 21-Jan-1998, for writePDBQState
      Real eval_emap[MAX_ATOMS]; // gmm added 21-Jan-1998, for writePDBQState
      Boole B_calcIntElec, B_isGaussTorCon, B_ShowTorE;
      State stateNow;
      unsigned short *US_TorE, (*US_torProfile)[NTORDIVS];
      int *type, (*tlist)[MAX_ATOMS];
      NonbondParam *nonbondlist;
      Real *charge, *abs_charge, *qsp_abs_charge;
      Real (*crd)[SPACE], (*vt)[SPACE], (*crdpdb)[SPACE], (*crdreo)[SPACE];
      EnergyTables *ptr_ad_energy_tables;
      Real (*map)[MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS];
      Boole *B_isTorConstrained;
      Molecule mol;
      int ignore_inter[MAX_ATOMS]; // gmm 2002-05-21, for CA, CB in flexible sidechains
      Boole         B_include_1_4_interactions; // gmm 2005-01-8, for scaling 1-4 nonbonds
      Real scale_1_4;                  // gmm 2005-01-8, for scaling 1-4 nonbonds
      //ParameterEntry *parameterArray;
      Real  unbound_internal_FE;
      Boole B_compute_intermol_energy; // use for computing unbound state
      Boole B_use_non_bond_cutoff;  // set this to FALSE if we are computing unbound extended conformations
      Boole B_have_flexible_residues;

   public:
      Eval(void);
      void setup( Real init_crd[MAX_ATOMS][SPACE],
                  Real  init_charge[MAX_ATOMS],
                  Real  init_abs_charge[MAX_ATOMS],
                  Real  init_qsp_abs_charge[MAX_ATOMS],
                  int            init_type[MAX_ATOMS], int init_natom,
                  Real  init_map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS],

                  Real  init_elec[MAX_ATOMS], // gmm added 21-Jan-1998, for writePDBQState
                  Real  init_emap[MAX_ATOMS], // gmm added 21-Jan-1998, for writePDBQState

                  NonbondParam *init_nonbondlist,
                  EnergyTables   *init_ptr_ad_energy_tables,
                  int init_Nnb,
                  Boole          init_B_calcIntElec,
                  Boole          init_B_isGaussTorCon, Boole init_B_isTorConstrained[MAX_TORS],
                  Boole          init_B_ShowTorE, unsigned short init_US_TorE[MAX_TORS],
                  unsigned short init_US_torProfile[MAX_TORS][NTORDIVS],
                  Real  init_vt[MAX_TORS][SPACE], int init_tlist[MAX_TORS][MAX_ATOMS],
                  Real  init_crdpdb[MAX_ATOMS][SPACE], 
                  Real  init_crdreo[MAX_ATOMS][SPACE], 
                  State stateInit, Molecule molInit,
                  int            init_ignore_inter[MAX_ATOMS],
                  Boole          init_B_include_1_4_interactions, // gmm 2005-01-8, for scaling 1-4 nonbonds
                  Real  init_scale_1_4,                   // gmm 2005-01-8, for scaling 1-4 nonbonds
                  //ParameterEntry init_parameterArray[MAX_ATOM_TYPES], // input  nonbond and desolvation parameters
                  Real  init_unbound_internal_FE,
                  GridMapSetInfo *init_info,
                  Boole  init_B_use_non_bond_cutoff,  // set this to FALSE if we are computing unbound extended conformations
                  Boole  init_B_have_flexible_residues
                  );
      void update_crds( Real init_crdreo[MAX_ATOMS][SPACE], 
                        Real init_vt[MAX_TORS][SPACE] );

      double operator()(Representation **);
      double operator()(Representation **, int); // GMM - allows calculation of a particular term of the total energy
#if defined(USING_COLINY)
      double operator()(double*, int);
#endif
      double eval();    // WEH - a basic change that facilitates the use of Coliny
      double eval(int); // GMM - allows calculation of a particular term of the total energy
      UnsignedFourByteLong evals(void);
      void reset(void);
      int write(FILE *out_file, Representation **rep);
      double evalpso(State *state);
      void compute_intermol_energy(Boole init_B_compute_intermol_energy); // for computing unbound state
};

inline Eval::Eval(void)
: num_evals(0)
{
}

inline void Eval::setup(Real init_crd[MAX_ATOMS][SPACE],
                        Real init_charge[MAX_ATOMS],
                        Real init_abs_charge[MAX_ATOMS],
                        Real init_qsp_abs_charge[MAX_ATOMS],
                        int init_type[MAX_ATOMS],
                        int init_natom,
                        Real init_map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS],

                        Real init_elec[MAX_ATOMS], // gmm added 21-Jan-1998, for writePDBQState
                        Real init_emap[MAX_ATOMS], // gmm added 21-Jan-1998, for writePDBQState
                        NonbondParam *init_nonbondlist,
                        EnergyTables   *init_ptr_ad_energy_tables,
                        int init_Nnb,
                        Boole init_B_calcIntElec, 
                        Boole init_B_isGaussTorCon,
                        Boole init_B_isTorConstrained[MAX_TORS],
                        Boole init_B_ShowTorE,
                        unsigned short init_US_TorE[MAX_TORS],
                        unsigned short init_US_torProfile[MAX_TORS][NTORDIVS],
                        Real init_vt[MAX_TORS][SPACE],
                        int init_tlist[MAX_TORS][MAX_ATOMS],
                        Real init_crdpdb[MAX_ATOMS][SPACE],
                        Real init_crdreo[MAX_ATOMS][SPACE],
                        State stateInit,
                        Molecule molInit,

                        int init_ignore_inter[MAX_ATOMS],

                        Boole init_B_include_1_4_interactions,
                        Real init_scale_1_4,

                        //ParameterEntry init_parameterArray[MAX_ATOM_TYPES], // input  nonbond and desolvation parameters

                        Real init_unbound_internal_FE,
                        GridMapSetInfo *init_info,
                        Boole init_B_use_non_bond_cutoff,  // set this to FALSE if we are computing unbound extended conformations
                        Boole init_B_have_flexible_residues
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
    B_isGaussTorCon = init_B_isGaussTorCon;
    B_isTorConstrained = init_B_isTorConstrained;
    B_ShowTorE = init_B_ShowTorE;
    US_TorE = init_US_TorE;
    US_torProfile = init_US_torProfile;
    vt = init_vt;
    tlist = init_tlist;
    crdpdb = init_crdpdb;
    crdreo = init_crdreo;
    // set all of the components of the State, one at a time...
    copyState( &stateNow, stateInit );
#ifdef DEBUG
    pr(logFile, "\n\nstateNow:\n");
    printState( logFile, stateNow, 2 );
#endif
    num_evals = 0;
    for (i=0; i<MAX_ATOMS; i++) {
       init_elec[i] = init_emap[i] = 0.0;
       ignore_inter[i] = init_ignore_inter[i];
    }
    mol = molInit;

    B_include_1_4_interactions = init_B_include_1_4_interactions;
    scale_1_4 = init_scale_1_4;

    //parameterArray = init_parameterArray;

    unbound_internal_FE = init_unbound_internal_FE;

    info = init_info;
    B_compute_intermol_energy = TRUE; // default is "Yes, calculate the intermolecular energy".

    B_use_non_bond_cutoff = init_B_use_non_bond_cutoff;  // set this to FALSE if we are computing unbound extended conformations
    B_have_flexible_residues = init_B_have_flexible_residues;
}

inline void Eval::update_crds( Real init_crdreo[MAX_ATOMS][SPACE], 
                               Real init_vt[MAX_TORS][SPACE] )
{
    crdreo = init_crdreo;
    vt = init_vt;
}

inline void Eval::compute_intermol_energy(Boole init_B_compute_intermol_energy)
    // For computing the conformation and the internal energy of the unbound state.
{
    B_compute_intermol_energy = init_B_compute_intermol_energy;
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
