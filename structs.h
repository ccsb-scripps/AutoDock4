/* structs.h */
#ifndef _STRUCTS_H
#define _STRUCTS_H

#include "constants.h"
#include "typedefs.h"

/* *****************************************************************************
 *      Name: structs.h                                                       *
 *  Function: Defines structures used in Molecular Applications.              *
 * Copyright: (C) Garrett Matthew Morris, TSRI                                *
 *----------------------------------------------------------------------------*
 *    Author: Garrett Matthew Morris, The Scripps Research Institute          *
 *      Date: SEP/07/1995                                                     *
 *----------------------------------------------------------------------------*
 *    Inputs: none                                                            *
 *   Returns: nothing                                                         *
 *   Globals: none                                                            *
 *----------------------------------------------------------------------------*
 * Modification Record                                                        *
 * Date     Inits   Comments                                                  *
 * 02/28/95 GMM     This header added                                         *
 ***************************************************************************** */

/* ____________________________________________________________________________ */

#ifdef USE_INT_AS_LONG
typedef int FourByteLong;
typedef unsigned int UnsignedFourByteLong;
#else
typedef long FourByteLong;
typedef unsigned long UnsignedFourByteLong;
#endif

/* ____________________________________________________________________________ */

typedef struct coord
{
  double x;			/* Cartesian x-coordinate */
  double y;			/* Cartesian y-coordinate */
  double z;			/* Cartesian z-coordinate */
} Coord;

/* ____________________________________________________________________________ */

typedef struct quat
{
  double nx;			/* unit vector's x-component */
  double ny;			/* unit vector's y-component */
  double nz;			/* unit vector's z-component */
  double ang;			/* angle of rotation about unit-vector */
  double x;			/* quaternion's x-component */
  double y;			/* quaternion's y-component */
  double z;			/* quaternion's z-component */
  double w;			/* quaternion's w-component */
  double qmag;			/* quaternion's 4-D magnitude */
} Quat;

/* ____________________________________________________________________________ */

typedef struct energy
{
  double total;			/* total energy */
  double intra;			/* intramolecular energy, a.k.a. "internal" energy */
  double inter;			/* intermolecular energy */
  double FE;			/* estimated Free Energy of binding */
} Energy;

/* ____________________________________________________________________________ */

typedef struct state
{
  Coord T;			/* coordinates of center of molecule */
  Quat Q;			/* rigid-body orientation */
  double tor[MAX_TORS];		/* torsion angles in radians */
  int ntor;			/* number of torsions in molecule */
  int hasEnergy;		/* if 0, this state has an undefined energy */
  Energy e;			/* energy structure */
} State;

/* ____________________________________________________________________________ */

typedef struct molecule
{
  FloatOrDouble crdpdb[MAX_ATOMS][SPACE];	/* original coordinates of atoms */
  FloatOrDouble crd[MAX_ATOMS][SPACE];	/* current coordinates of atoms */
  char atomstr[MAX_ATOMS][MAX_CHARS];	/* strings describing atoms, from PDB file, cols,1-30. */
  int natom;			/* number of atoms in molecule */
  FloatOrDouble vt[MAX_TORS][SPACE];	/* vectors  of torsions */
  int tlist[MAX_TORS][MAX_ATOMS];	/* torsion list of movable atoms */
  State S;			/* state of molecule */
} Molecule;

/* ____________________________________________________________________________ */

typedef struct rotamer
{
  double tor[MAX_TORS_IN_ROTAMER];	/* torsion angles in radians */
  int ntor;			/* number of torsions */
} Rotamer;

/* ____________________________________________________________________________ */

typedef struct atom
{
  Coord pt;			/* transformed point */
  Coord pt0;			/* untransformed point, original PDB coords */
  double q;			/* partial atomic charge */
  Boole Bq;			/* TRUE if the atom has a charge */
  int type;			/* atom type */
  Boole isH;			/* TRUE if atom is a hydrogen */
  int id;			/* serial ID */
  char name[5];			/* PDB atom name; formerly "pdbaname" */
  char str[MAX_CHARS];		/* PDB atom string; formerly "atomstuff" */
  int nnb;			/* number of non-bonds for this atom */
} Atom;
/* ____________________________________________________________________________ */

typedef struct bond
{
  Atom *atom1;
  Atom *atom2;
  double bondLength;
  Coord bondVector;
} Bond;
/* ____________________________________________________________________________ */

typedef struct pair_id
{
  Atom *atom1;			/* pointer to one atom in pair */
  Atom *atom2;			/* pointer to other atom */
} PairID;
/* ____________________________________________________________________________ */

typedef struct dist_constraint
{
  PairID bond;			/* two atoms defining distance constraint */
  double length;		/* current bond length */
  double lower;			/* lower bound on distance */
  double upper;			/* upper bound on distance */
} DisCon;

/* ____________________________________________________________________________ */

typedef struct torsion
{
  PairID rotbnd;		/* atom serial-IDs of rotatable bond */
  int nmoved;			/* number of atoms moved by this */
  int IDmove[MAX_ATOMS];	/* atom serial-IDs of atoms moved by this */
  Coord vt;			/* bond-vector of rotatable bond */
} Torsion;

/* ______________________________________________________________________________
** Molecular fragments, side-chains, even entire ligands */

typedef struct group
{
  int natom;			/* Number of atoms in fragment */
  Atom atm[MAX_ATOMS];		/* Atom data */
  int ntor;			/* Number of torsions in fragment */
  Torsion tors[MAX_TORS];	/* Torsion data */
  int nnb;			/* Number of non-bonds in fragment */
  PairID nbs[MAX_NONBONDS];	/* Non-bond data */
  Boole B_constrain;		/* TRUE if any distance constraints */
  DisCon distcon;		/* Distance constraint data */
  char pdbqfilnam[MAX_CHARS];	/* PDBQ filename holding these data */
} Group;

/* ______________________________________________________________________________
** Grid Map */

typedef struct grid_map
{
    double map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS];
    Coord center;
    Coord max;
    Coord min;
    double spacing;
    int atom_type;
    int num_points[3];
    int num_points_1[3];
} GridMap;

/* ______________________________________________________________________________
** Parameter Dictionary */

#include <search.h>

#define MAX_NUM_AUTOGRID_TYPES 100
#define MAX_LEN_AUTOGRID_TYPE 7

enum hbond_type
{ NON, DS, D1, AS, A1, A2 };	/* hbonding character: */

typedef struct parameter_entry
{				// was "parm_info" in earlier AutoGrid 4 code
  char autogrid_type[MAX_LEN_AUTOGRID_TYPE + 1];	/* autogrid_type is based on babel_types assigned by PyBabel */
  double Rij;			/* Lennard-Jones equilibrium separation */
  double epsij;			/* Lennard-Jones energy well-depth */
  double vol;			/* solvation volume */
  double solpar;		/* solvation parameter */
  hbond_type hbond;		/* hbonding character: 
				   NON: none, 
				   DS: spherical donor 
				   D1: directional donor
				   AS: spherical acceptor
				   A1: acceptor of 1 directional hbond
				   A2: acceptor of 2 directional hbonds */
  double Rij_hb;		/* 12-10 Lennard-Jones equilibrium separation */
  double epsij_hb;		/* 12-10 Lennard-Jones energy well-depth */
  int rec_index;		/* used to set up receptor atom_types */
  int map_index;		/* used to set up map atom_types */
  int bond_index;		/* used to set up bonds; corresponds to the enum in mdist.h */
} ParameterEntry;

typedef struct linear_FE_model
{
    double coeff_vdW;                 // Free energy coefficient for van der Waals term
    double coeff_hbond;               // Free energy coefficient for H-bonding term
    double coeff_estat;               // Free energy coefficient for electrostatics term
    double coeff_desolv;              // Free energy coefficient for desolvation term
    double coeff_tors;                // Free energy coefficient for torsional term

    double stderr_vdW;                // Free energy standard error for van der Waals term
    double stderr_hbond;              // Free energy standard error for H-bonding term
    double stderr_estat;              // Free energy standard error for electrostatics term
    double stderr_desolv;             // Free energy standard error for desolvation term
    double stderr_tors;               // Free energy standard error for torsional term
} Linear_FE_Model;

#endif
/* EOF */
