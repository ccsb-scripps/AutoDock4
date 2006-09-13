#include "conformation_sampler.h"
#include "hybrids.h"
#include "ranlib.h"
#include <math.h>

#define AVOGADRO 6.022e23f
#define RK_CONSTANT 0.0019872065 // constant in entropy calculation
#define TEMP 298 // temperature in entropy calculation

#define RMSD_SYMMETRY TRUE
#define TRAN_STEP 0.03 // size of translation steps (x,y,z)
#define ROT_ANG_STEP 0.025 // size of step for rotation angle
#define TOR_ANG_STEP 0.03 // size of step for torsion angles
#define RHO 0.2 // parameter for random sampling

#define DEFAULT_RANDOM_SAMPLES 10000
#define DEFAULT_INCREMENTAL_STEPS 5 // note that this is steps up and down, i.e. 
									// X +/- 3
#define ROTATION_ANGLE_INDEX 6
#define Z_TRANSLATON_INDEX 2
#define X_ROTATION_INDEX 3
#define Y_ROTATION_INDEX 4
#define Z_ROTATION_INDEX 5

#define ICO_X 0.525731112119133606
#define ICO_Y 0.850650808352039932

extern class Eval evaluate;
extern FILE *logFile;

Real (*vt)[SPACE], (*crdpdb)[SPACE];
int (*tlist)[MAX_ATOMS];
Real *lig_center;
int natom;
int *type;
GridMapSetInfo *info;

Real crd[MAX_ATOMS][SPACE];
Real ref_crd[MAX_ATOMS][SPACE];

const Real vertices[12][3] = {{-ICO_X, 0., ICO_Y}, {ICO_X, 0., ICO_Y}, {-ICO_X, 0., -ICO_Y}, {ICO_X, 0., -ICO_Y},
                              {0., ICO_Y, ICO_X}, {0., ICO_Y, -ICO_X}, {0., -ICO_Y, ICO_X}, {0., -ICO_Y, -ICO_X},
                              {ICO_Y, ICO_X, 0.}, {-ICO_Y, ICO_X, 0.}, {ICO_Y, -ICO_X, 0.}, {-ICO_Y, -ICO_X, 0.}};

ConformationSampler::ConformationSampler(State init_state) {
	base_state = init_state;
	base_ind = set_ind(info, init_state);
	base_point = base_ind.phenotyp;	
	base_energy = base_point.evaluate(Normal_Eval);
	cnv_state_to_coords(init_state, vt, tlist, init_state.ntor, crdpdb, base_crd, natom);
	
	dimensionality = BASE_DIMENSIONS + init_state.ntor;
	evals = 0;
	favorable_evals = 0;
	total_energy = 0.0;
	total_favorable_energy = 0.0;
	min_energy = base_energy;
	min_energy_rmsd = 0.0;
	
	// set up the temp variables
	probe_state = base_state;
	probe_ind = base_ind;
	probe_point = base_point;

	// store axis/angle representation in an array
	for (unsigned int i=0; i < 4; i++) {
		base_axis_angle[i] = probe_point.gread(i+X_ROTATION_INDEX).real;
	}
	
	// initialize bounds
	for (int i=0; i < 3; i++) {
		max_values[i] = min_values[i] = base_point.gread(i).real;
		max_values[i+3] = min_values[i+3] = 0.0; // Eulerian angles set to 0
		
		rotation_angles[i] = 0.0;
	}
	
	// reset bins
	for (int i=0; i < NUM_BINS; i++) {
		bin_total_energy[i] = 0.0;
		bin_total_favorable_energy[i] = 0.0;
		bin_count[i] = 0;
		bin_count_favorable[i] = 0;
		bin_min_energy[i] = 0.0;
		bin_max_energy[i] = base_energy;
	}
}

void ConformationSampler::random_sample(void) {
	random_sample(1);
}

void ConformationSampler::random_sample(int num_samples) {
	Real multiplier;
	
	for (int sample=0; sample < num_samples; sample++) {
		probe_point = base_point;
		multiplier = ranf();
		//multiplier = genunf(0.0, 1.25);

		// perturb translation and torsion angles randomly
		for (unsigned int i=0; i < (unsigned int) dimensionality; i++) {
			if (i >= X_ROTATION_INDEX && i <= Z_ROTATION_INDEX) continue;
			probe_point.write(probe_point.gread(i).real + gennor(0.0, multiplier*RHO) , i);
			//probe_point.write(probe_point.gread(i).real + genunf(-1.0 * RHO, RHO) , i);
		}
		
		// adjust current rotation randomly using multiplication
		Real random_axis_angle[4];
		for (int i=0; i < 4; i++) {
			random_axis_angle[i] = gennor(0.0, multiplier*RHO);
			//random_axis_angle[i] = genunf(-1.0 * RHO, RHO);
		}
		
		Real new_axis_angle[4];
		multiplyraa(base_axis_angle, random_axis_angle, new_axis_angle);
		for (unsigned int i=0; i < 4; i++) {
			probe_point.write(new_axis_angle[i], i+X_ROTATION_INDEX);
		}

		current_energy();
	}
}

ConformationSampler::~ConformationSampler(void) {
}

// NOTE: currently, the torsional free energy penalty is not included.
// Since this is an entropic term, I believe we can ignore it in this analysis.
Real ConformationSampler::current_energy(void) {
	evals++;
	Real energy = probe_point.evaluate(Normal_Eval);
	Real rmsd = current_rmsd();
	
	total_energy += energy;
	
	// store information on minimum energy conformation
	if (energy < min_energy) {
		min_energy = energy;
		min_energy_rmsd = rmsd;
	}
	
	int bin = (int)(rmsd/BIN_SIZE);
	if (bin < NUM_BINS) {
		bin_count[bin]++;
		bin_total_energy[bin] += energy;
		if (energy < bin_min_energy[bin]) bin_min_energy[bin] = energy;
		if (energy > bin_max_energy[bin]) bin_max_energy[bin] = energy;
	}
	
	if (energy < 0.0) {
		favorable_evals++;
		total_favorable_energy += energy;
		update_bounds();
		
		if (bin < NUM_BINS) {
			bin_count_favorable[bin]++;
			bin_total_favorable_energy[bin] += energy;
		}
	}
	return energy;
}

Real ConformationSampler::current_rmsd(void) {
	probe_ind.phenotyp = probe_point;
	probe_ind.inverse_mapping();
	probe_state = probe_ind.state(base_state.ntor);
	cnv_state_to_coords(probe_state, vt, tlist, probe_state.ntor, crdpdb, crd, natom);
	return getrms(crd, base_crd, RMSD_SYMMETRY, natom, type);
}

Real ConformationSampler::reference_rmsd(void) {
	return getrms(base_crd, ref_crd, RMSD_SYMMETRY, natom, type);
}

void ConformationSampler::update_bounds(void) {
	Real euler[3];
	Real raa[4];
	Real current_value;
	
	// set up axis-angle array
	for (int i=0; i < 4; i++) {
		raa[i] = probe_point.gread(i+3).real;
	}
	raaEuler(raa, euler);
	
	// check existing translation bounds
	for (int i=0; i < 3; i++) {
		current_value = probe_point.gread(i).real;
		if (current_value < min_values[i]) min_values[i] = current_value;
		else if (current_value > max_values[i]) max_values[i] = current_value; 
	}
	
	// check rotation bounds
	for (int i=0; i < 3; i++) {
		if (euler[i] < min_values[i+3]) min_values[i+3] = euler[i];
		else if (euler[i] > max_values[i+3]) max_values[i+3] = euler[i];
	}
}

void ConformationSampler::systematic_search(int index) {
	
	// for rotation axes, rotate using the pre-defined vertices 
	if (index <= Z_ROTATION_INDEX && index >= X_ROTATION_INDEX) {
		
		Real temp_axis_angle[4];
		Real new_axis_angle[4];
		
		for (int i=0; i < 12; i++) {
			
			// set up rotation
			temp_axis_angle[0] = vertices[i][0];
			temp_axis_angle[1] = vertices[i][1];
			temp_axis_angle[2] = vertices[i][2];
			temp_axis_angle[3] = temp_rotation_angle;
			multiplyraa(base_axis_angle, temp_axis_angle, new_axis_angle);
			
			probe_point.write(new_axis_angle[0], X_ROTATION_INDEX);
			probe_point.write(new_axis_angle[1], Y_ROTATION_INDEX);
			probe_point.write(new_axis_angle[2], Z_ROTATION_INDEX);
			probe_point.write(new_axis_angle[3], ROTATION_ANGLE_INDEX);
			
			systematic_search(Z_TRANSLATON_INDEX); // go to translation
		}
	}
	
	// translation, rotation angle, and torsion angles
	// step through 
	else {
		int num_steps = DEFAULT_INCREMENTAL_STEPS; // steps up or down
		Real start, step_size;
			
		// set step sizes for different dimensions
		if (index <= Z_TRANSLATON_INDEX) step_size = TRAN_STEP;
		else if (index == ROTATION_ANGLE_INDEX) step_size = ROT_ANG_STEP;
		else step_size = TOR_ANG_STEP;
		
		// for the rotation angle, use different bounds in order to avoid
		// symmetry problems
		if (index == ROTATION_ANGLE_INDEX) start = -2*num_steps*step_size;
		else start = base_point.gread(index).real - num_steps * step_size;
		
		for (int current = 0; current <= 2 * num_steps; current++) {

			if (index == ROTATION_ANGLE_INDEX) {
				temp_rotation_angle = start + current*step_size;
			}
			else {
				probe_point.write(start + current * step_size, index);
			}		
			
			if (index == 0) {
				(void)current_energy();
			}
			else {				
				// check if the rotation angle is 0, to avoid shifting axis unnecessarily
				if (index == ROTATION_ANGLE_INDEX && current == 2 * num_steps) {
					probe_point.write(base_point.gread(X_ROTATION_INDEX), X_ROTATION_INDEX);
					probe_point.write(base_point.gread(Y_ROTATION_INDEX), Y_ROTATION_INDEX);
					probe_point.write(base_point.gread(Z_ROTATION_INDEX), Z_ROTATION_INDEX);
					probe_point.write(base_point.gread(ROTATION_ANGLE_INDEX), ROTATION_ANGLE_INDEX);
					systematic_search(Z_TRANSLATON_INDEX); // skip to translation
					//current_energy();// DEBUGGING
				}
				else {
					systematic_search(index-1);
				}
			}
		} 
	}
}

Real ConformationSampler::fraction_favorable(void) {
	return favorable_evals/evals;
}

Real ConformationSampler::average_favorable_energy(void) {
	if (favorable_evals == 0) return 0;
	else return total_favorable_energy/favorable_evals;
}

Real ConformationSampler::energy_volume(void) {
	return total_favorable_energy/evals;
}

Real ConformationSampler::configurational_integral(void) {
	Real Vb = 1.0;
	for (int i=0; i < 6; i++) {
		Vb *= (max_values[i]-min_values[i]);
	}
	return Vb;
}

/*
 * estimate entropy, as described by Ruvinsky and Kozintsev
 * return (0.0019872065)*(298)*math.log(Vb * 6.02 * 10**23 / (8*math.pi*math.pi))
 */

Real ConformationSampler::RK_entropy(void) {
	return RK_CONSTANT * TEMP * log(configurational_integral() * AVOGADRO/ (8 * PI * PI));
}

void ConformationSampler::output_statistics(void) {
	fprintf(logFile, "Conformation starting energy: %.3f\n", base_energy);
	fprintf(logFile, "RMSD from reference state: %.3f\n", reference_rmsd());
	fprintf(logFile, "Fraction of favorable evaluations: %.3f\n", (Real)favorable_evals/evals);
	fprintf(logFile, "Average favorable energy: %.3f\n", total_favorable_energy/favorable_evals);
	fprintf(logFile, "Estimated energy volume: %.3f\n", total_favorable_energy/evals);
	fprintf(logFile, "Minimum energy found: %.3f (%.3f A from starting point)\n", min_energy, min_energy_rmsd);
	//fprintf(logFile, "Bins in local region.\n");
	fprintf(logFile, "\nRMSD       #     fraction  Volume    Avg. (-)   Min E     Max E\n");
	for (int i=0; i < NUM_BINS; i++) {
		fprintf(logFile, "%.1f    %7d    %2.3f    %2.3f    %2.3f    %2.3f    %2.3f\n", (i+1)*BIN_SIZE, bin_count[i], (Real)bin_count_favorable[i]/bin_count[i], bin_total_favorable_energy[i]/bin_count[i], bin_total_favorable_energy[i]/bin_count_favorable[i], bin_min_energy[i], bin_max_energy[i]);
	}
	fprintf(logFile, "%d evaluations.\n\n", evals);
}

void systematic_conformation_sampler(State hist[MAX_RUNS], int nconf, Real init_vt[MAX_TORS][SPACE], Real init_crdpdb[MAX_ATOMS][SPACE], int init_tlist[MAX_TORS][MAX_ATOMS], Real init_lig_center[SPACE], int init_natom, int init_type[MAX_ATOMS], GridMapSetInfo *init_info) {
	vt = init_vt;
	crdpdb = init_crdpdb;
	tlist = init_tlist;
	lig_center = init_lig_center;
	natom = init_natom;
	type = init_type;
	info = init_info;
	
	setup_reference_coordinates();
	
	fprintf(logFile, "Initiating a systematic search.\n");
	for (int i=0; i < nconf; i++) {
		fprintf(logFile, "\nConformation %d:\n", i+1);
		State base_state = hist[i];
		ConformationSampler CS(base_state);
		//CS.systematic_search(CS.dimensionality-1);
		CS.systematic_search(BASE_DIMENSIONS-1);
		CS.output_statistics();
	}
	fprintf(logFile,"\n\n");
}

void random_conformation_sampler(State hist[MAX_RUNS], int nconf, int num_samples, Real init_vt[MAX_TORS][SPACE], Real init_crdpdb[MAX_ATOMS][SPACE], int init_tlist[MAX_TORS][MAX_ATOMS], Real init_lig_center[SPACE], int init_natom, int init_type[MAX_ATOMS], GridMapSetInfo *init_info) {
	vt = init_vt;
	crdpdb = init_crdpdb;
	tlist = init_tlist;
	lig_center = init_lig_center;
	natom = init_natom;
	type = init_type;
	info = init_info;
	
	setup_reference_coordinates();
	
	if (num_samples == 0) num_samples = DEFAULT_RANDOM_SAMPLES;
	
	fprintf(logFile, "Initiating a random search using %d samples near each conformation.\n", num_samples);
	for (int i=0; i < nconf; i++) {
		fprintf(logFile, "\nConformation %d:\n", i+1);
		State base_state = hist[i];
		ConformationSampler CS(base_state);
		CS.random_sample(num_samples);
		CS.output_statistics();
	}
	
	fprintf(logFile,"\n\n");
}


/* copied (and slightly modified) from non-included code in call_glss.cc */
Individual set_ind(GridMapSetInfo *info, State state)
{
   Genotype temp_Gtype;
   Phenotype temp_Ptype;
   int i;

   temp_Gtype = generate_Gtype(state.ntor, info);
   temp_Ptype = generate_Ptype(state.ntor, info);

   // use the state to generate a Genotype
   temp_Gtype.write(state.T.x, 0);
   temp_Gtype.write(state.T.y, 1);
   temp_Gtype.write(state.T.z, 2);
   temp_Gtype.write(state.Q.nx, 3);
   temp_Gtype.write(state.Q.ny, 4);
   temp_Gtype.write(state.Q.nz, 5);
   temp_Gtype.write(state.Q.ang, 6);
   for (i=0;i<state.ntor; i++) {
       temp_Gtype.write(state.tor[i], 7+i);
   };

   Individual temp(temp_Gtype, temp_Ptype);   

   // use mapping to generate a Phenotype
   temp.phenotyp =  temp.mapping();

   return(temp);
}

void raaEuler(Real raa[4], Real euler[3]) {
	Real s = sin(raa[3]);
	Real c = cos(raa[3]);
	Real t = 1.0 - c;
	
	// check for singularities
	if (raa[0]*raa[1]*t + raa[2]*s > 0.998) {
		euler[0] = 0.0;
		euler[1] = atan2(raa[0]*sin(raa[3]/2), cos(raa[3]/2));
		euler[2] = PI/2;
	}
	
	else if (raa[0]*raa[1]*t + raa[2]*s < -0.998) {
		euler[0] = 0.0;
		euler[1] = -atan2(raa[0]*sin(raa[3]/2), cos(raa[3]/2));
		euler[2] = -PI/2;
	}
	
	euler[0] = atan2(raa[0]*s - raa[1]*raa[2]*t , 1 - (raa[0]*raa[0] + raa[2]*raa[2])*t);
	euler[1] = atan2(raa[1]*s - raa[0]*raa[2]*t , 1 - (raa[1]*raa[1] + raa[2]*raa[2])*t);
	euler[2] = asin(raa[0]*raa[1]*t + raa[2]*s);
}

void raaMatrix(Real raa[4], Real matrix[3][3]) {
	Real angle_cos = cos(raa[3]);
	Real angle_sin = sin(raa[3]);
	Real t = 1.0 - angle_cos;
	
	// make sure that input vecotr is a unit vector
	Real length = hypotenuse(raa[0], raa[1], raa[2]);
	raa[0] /= length;
	raa[1] /= length;
	raa[2] /= length;
	
	matrix[0][0] = angle_cos + raa[0]*raa[0]*t;
	matrix[1][1] = angle_cos + raa[1]*raa[1]*t;
	matrix[2][2] = angle_cos + raa[2]*raa[2]*t;
	
	Real tmp1 = raa[0]*raa[1]*t;
	Real tmp2 = raa[2]*angle_sin;
	matrix[1][0] = tmp1 + tmp2;
	matrix[0][1] = tmp1 - tmp2;
	
	tmp1 = raa[0]*raa[2]*t;
	tmp2 = raa[1]*angle_sin;
	matrix[2][0] = tmp1 - tmp2;
	matrix[0][2] = tmp1 + tmp2;
	
	tmp1 = raa[1]*raa[2]*t;
	tmp2 = raa[0]*angle_sin;
	matrix[2][1] = tmp1 + tmp2;
	matrix[1][2] = tmp1 - tmp2;
}

void matrixraa(Real matrix[3][3], Real raa[4]) {
	Real length = hypotenuse(matrix[2][1] - matrix[1][2], matrix[2][0] - matrix[0][2], matrix[1][0] - matrix[0][1]);
	
	// need to check acos() parameter to avoid values out of range
	Real cosine = (matrix[0][0] + matrix[1][1] + matrix[2][2] - 1)/2;
	if (cosine > 1.0) raa[3] = 0;
	else if (cosine < -1.0) raa[3] = PI;
	else raa[3] = acos(cosine);
	
	raa[0] = (matrix[2][1] - matrix[1][2])/length;
	raa[1] = (matrix[0][2] - matrix[2][0])/length;
	raa[2] = (matrix[1][0] - matrix[0][1])/length;
}

void multiplyraa(Real raa1[4], Real raa2[4], Real raa_result[4]) {
	Real matrix1[3][3];
	Real matrix2[3][3];
	Real result_matrix[3][3];
	
	raaMatrix(raa1, matrix1);
	raaMatrix(raa2, matrix2);
	matrixMultiply(matrix1, matrix2, result_matrix);
	matrixraa(result_matrix, raa_result);
}

void matrixMultiply(Real m1[3][3], Real m2[3][3], Real result[3][3]) {
	result[0][0] = m1[0][0]*m2[0][0] + m1[0][1]*m2[1][0] + m1[0][2]*m2[2][0];
	result[0][1] = m1[0][0]*m2[0][1] + m1[0][1]*m2[1][1] + m1[0][2]*m2[2][1];
	result[0][2] = m1[0][0]*m2[0][2] + m1[0][1]*m2[1][2] + m1[0][2]*m2[2][2];
	result[1][0] = m1[1][0]*m2[0][0] + m1[1][1]*m2[1][0] + m1[1][2]*m2[2][0];
	result[1][1] = m1[1][0]*m2[0][1] + m1[1][1]*m2[1][1] + m1[1][2]*m2[2][1];
	result[1][2] = m1[1][0]*m2[0][2] + m1[1][1]*m2[1][2] + m1[1][2]*m2[2][2];
	result[2][0] = m1[2][0]*m2[0][0] + m1[2][1]*m2[1][0] + m1[2][2]*m2[2][0];
	result[2][1] = m1[2][0]*m2[0][1] + m1[2][1]*m2[1][1] + m1[2][2]*m2[2][1];
	result[2][2] = m1[2][0]*m2[0][2] + m1[2][1]*m2[1][2] + m1[2][2]*m2[2][2];
}

void setup_reference_coordinates(void) {
	for (int i = 0;  i < natom;  i++) {
		ref_crd[i][0] = lig_center[0] + crdpdb[i][0];
		ref_crd[i][1] = lig_center[1] + crdpdb[i][1];
		ref_crd[i][2] = lig_center[2] + crdpdb[i][2];
	}
}