/*

 $Id: conformation_sampler.h,v 1.7 2007/04/27 06:01:48 garrett Exp $

 AutoDock 

 Copyright (C) 1989-2007,  Max Chang, Rik Belew, Garrett M. Morris, David S. Goodsell, Ruth Huey, Arthur J. Olson, 
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

#include "constants.h"
#include "structs.h"
#include "support.h"
#include "getrms.h"
#include "autocomm.h"
#include "cnv_state_to_coords.h"
#include "getrms.h"

#ifndef _CONFORMATION_SAMPLER_H
#define _CONFORMATION_SAMPLER_H

#define BASE_DIMENSIONS 7
#define NUM_BINS 20
#define BIN_SIZE 0.1

class ConformationSampler {
	public:
		State base_state, probe_state;
		Individual base_ind, probe_ind;
		Phenotype base_point, probe_point;
		Real base_axis_angle[4];
        Quat base_q;
		Real base_crd[MAX_ATOMS][SPACE]; //probe_crd?;
		Real base_energy, total_energy, total_favorable_energy;
		Real min_energy, min_energy_rmsd;
		Real Boltzmann_sum, Boltzmann_diff_sum;
		int dimensionality, evals, favorable_evals;
		Real temp_rotation_angle;

		Real min_values[BASE_DIMENSIONS-1];
		Real max_values[BASE_DIMENSIONS-1];
		
		int bin_count[NUM_BINS];
		int bin_count_favorable[NUM_BINS];
		Real bin_total_energy[NUM_BINS];
		Real bin_total_favorable_energy[NUM_BINS];
		Real bin_min_energy[NUM_BINS];
		Real bin_max_energy[NUM_BINS];
		Real bin_Boltzmann_sum[NUM_BINS];

		ConformationSampler(State);
		~ConformationSampler(void);

		void random_sample(void);
		void random_sample(int);
		void systematic_search(int index);
		Real current_energy(void);
		Real current_rmsd(void);
		Real reference_rmsd(void);
		Real fraction_favorable(void);
		Real average_favorable_energy(void);
		Real energy_volume(void);
		Real RK_entropy(void);
		void output_statistics(void);
		Real partition_function(void);
		Real partition_function(int bin);
		Real entropy_estimate(void);
	private:
		Real normalized_volume(void);
		Real normalized_Boltzmann(void);
		Real configurational_integral(void);
		void update_bounds(void);
};

void systematic_conformation_sampler(State hist[MAX_RUNS], int nconf, Real init_vt[MAX_TORS][SPACE], Real init_crdpdb[MAX_ATOMS][SPACE], int init_tlist[MAX_TORS][MAX_ATOMS], Real init_lig_center[SPACE], int init_natom, int init_type[MAX_ATOMS], GridMapSetInfo *init_info);
void random_conformation_sampler(State hist[MAX_RUNS], int nconf, int num_samples, Real init_vt[MAX_TORS][SPACE], Real init_crdpdb[MAX_ATOMS][SPACE], int init_tlist[MAX_TORS][MAX_ATOMS], Real init_lig_center[SPACE], int init_natom, int init_type[MAX_ATOMS], GridMapSetInfo *init_info);
Individual set_ind(GridMapSetInfo *info, State state);
void raaEuler(Real raa[4], Real euler[3]);
void testMatrix(void);
void raaMatrix(Real raa[4], Real matrix[3][3]);
void matrixraa(Real matrix[3][3], Real raa[4]);
void multiplyraa(Real raa1[4], Real raa2[4], Real raa_result[4]);
void matrixMultiply(Real m1[3][3], Real m2[3][3], Real result[3][3]);
void rand_axis(Real axis[4], Real angle);
void setup_reference_coordinates(void);

#endif
