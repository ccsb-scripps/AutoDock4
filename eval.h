/********************************************************************
    The header file for the eval class

                                rsh 09/06/95
********************************************************************/
#ifndef _EVAL_H
#define _EVAL_H

#include "structs.h"
#include "rep.h"

#include <stdio.h>
#include "qmultiply.h"
#include "cnv_state_to_coords.h"
#include "trilinterp.h"
#include "eintcal.h"
#include "energy.h"

void make_state_from_rep(Representation **rep, State *stateNow);

class Eval
{
   private:
      UnsignedFourByteLong num_evals;
      int natom, Nnb;
      float inv_spacing, xlo, xhi, ylo, yhi, zlo, zhi;  // gmm added new private members xhi,yhi,zhi
      float xcen, ycen, zcen; // gmm added 14-Jan-1998, center of grid
      float eval_elec[MAX_ATOMS]; // gmm added 21-Jan-1998, for writePDBQState
      float eval_emap[MAX_ATOMS]; // gmm added 21-Jan-1998, for writePDBQState
      Boole B_calcIntElec, B_isGaussTorCon, B_ShowTorE;
      State stateNow;
      unsigned short *US_TorE, (*US_torProfile)[NTORDIVS]; 
      int *type, (*nonbondlist)[2], (*tlist)[MAX_ATOMS];
//      float (*q1q2), *Addr_eintra, (*charge);
      float *q1q2, *charge;
      float (*crd)[SPACE], (*vt)[SPACE], (*crdpdb)[SPACE]; 
      float (*e_internal)[ATOM_MAPS][ATOM_MAPS]; 
      float (*map)[MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS];
      Boole *B_isTorConstrained; 
      Molecule mol;
      Boole B_template; // Use the template-docking scoring function if true 15-jan-2001
      // float template_energy[MAX_ATOMS]; // atomic template values
      // float template_stddev[MAX_ATOMS]; // atomic template values
      float *template_energy; // atomic template values
      float *template_stddev; // atomic template values
   
   public:
      Eval(void);
      void setup(float init_crd[MAX_ATOMS][SPACE], float init_charge[MAX_ATOMS], 
            int   init_type[MAX_ATOMS], int init_natom, 
            float init_map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS], 
            float init_inv_spacing, 
            float init_elec[MAX_ATOMS], // gmm added 21-Jan-1998, for writePDBQState
            float init_emap[MAX_ATOMS], // gmm added 21-Jan-1998, for writePDBQState
            float init_xlo, float init_xhi, // added by GMM
            float init_ylo, float init_yhi, // xhi,yhi,zhi
            float init_zlo, float init_zhi, // ...
//            float init_Addr_eintra, int init_nonbondlist[MAX_NONBONDS][2], 
            int init_nonbondlist[MAX_NONBONDS][2], 
            float init_e_internal[NEINT][ATOM_MAPS][ATOM_MAPS], int init_Nnb, 
            Boole init_B_calcIntElec, float init_q1q2[MAX_NONBONDS], 
            Boole init_B_isGaussTorCon, Boole init_B_isTorConstrained[MAX_TORS], 
            Boole init_B_ShowTorE, unsigned short init_US_TorE[MAX_TORS], 
            unsigned short init_US_torProfile[MAX_TORS][NTORDIVS], 
            float init_vt[MAX_TORS][SPACE], int init_tlist[MAX_TORS][MAX_ATOMS], 
            float init_crdpdb[MAX_ATOMS][SPACE], State stateInit, Molecule molInit,
            Boole init_B_template, 
            float init_template_energy[MAX_ATOMS], 
            float init_template_stddev[MAX_ATOMS]);
      double operator()(Representation **);
      UnsignedFourByteLong evals(void);
      void reset(void);
      int write(FILE *out_file, Representation **rep);
};

inline Eval::Eval(void)
: num_evals(0)
{
}

inline void Eval::setup(float init_crd[MAX_ATOMS][SPACE], float init_charge[MAX_ATOMS], 
            int init_type[MAX_ATOMS], int init_natom, 
            float init_map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS], 
            float init_inv_spacing, 
            float init_elec[MAX_ATOMS], // gmm added 21-Jan-1998, for writePDBQState
            float init_emap[MAX_ATOMS], // gmm added 21-Jan-1998, for writePDBQState
            float init_xlo, float init_xhi,   // gmm
            float init_ylo, float init_yhi,   // gmm
            float init_zlo, float init_zhi,   // gmm
//            float init_Addr_eintra, int init_nonbondlist[MAX_NONBONDS][2], 
            int init_nonbondlist[MAX_NONBONDS][2], 
            float init_e_internal[NEINT][ATOM_MAPS][ATOM_MAPS], int init_Nnb, 
            Boole init_B_calcIntElec, float init_q1q2[MAX_NONBONDS], 
            Boole init_B_isGaussTorCon, Boole init_B_isTorConstrained[MAX_TORS], 
            Boole init_B_ShowTorE, unsigned short init_US_TorE[MAX_TORS], 
            unsigned short init_US_torProfile[MAX_TORS][NTORDIVS], 
            float init_vt[MAX_TORS][SPACE], int init_tlist[MAX_TORS][MAX_ATOMS], 
            float init_crdpdb[MAX_ATOMS][SPACE], State stateInit,
            Molecule molInit, 
            Boole init_B_template, 
            float init_template_energy[MAX_ATOMS], 
            float init_template_stddev[MAX_ATOMS])
{
   register int i;

   crd = init_crd; 
   charge = init_charge; 
   type = init_type; 
   natom = init_natom; 
   map = init_map; 
   inv_spacing = init_inv_spacing; 
   xlo = init_xlo; 
   xhi = init_xhi; // gmm
   xcen = 0.5*(xhi+xlo); // gmm 14-jan-98
   ylo = init_ylo; 
   yhi = init_yhi; // gmm
   ycen = 0.5*(yhi+ylo); // gmm 14-jan-98
   zlo = init_zlo; 
   zhi = init_zhi; // gmm
   zcen = 0.5*(zhi+zlo); // gmm 14-jan-98
//   Addr_eintra = init_Addr_eintra; 
   nonbondlist = init_nonbondlist; 
   e_internal = init_e_internal; 
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
   }
   mol = molInit;
   B_template = init_B_template;
   template_energy = init_template_energy;
   template_stddev = init_template_stddev;
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
