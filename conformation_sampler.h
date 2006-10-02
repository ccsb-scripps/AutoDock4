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
#define NUM_BINS 10
#define BIN_SIZE 0.1

class ConformationSampler {
	public:
		State base_state, probe_state;
		Individual base_ind, probe_ind;
		Phenotype base_point, probe_point;
		Real base_axis_angle[4];
		Real base_crd[MAX_ATOMS][SPACE]; //probe_crd?;
		Real base_energy, total_energy, total_favorable_energy;
		Real min_energy, min_energy_rmsd;
		Real Boltzmann_sum;
		int dimensionality, evals, favorable_evals;
		Real temp_rotation_angle;
		Real rotation_angles[3];

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
void setup_reference_coordinates(void);

#endif
