#ifndef COPYDIMENSION
#define COPYDIMENSION

#include "structs.h"
#include "constants.h"
#include "alea.h"


//void copyDimension( State *S, Position R);
//void copyState2Dimension( Position *R, State S);
void make_state_from_position(State *S, Position P);
void copy_state_to_position(Position *R , State S);
void initialiseDimension(GridMapSetInfo *info, double *xmin, double *xmax, int D);
void initializeStateDimensions(GridMapSetInfo *info, double *xmin, double *xmax, 
				double pso_qstep, double pso_dstep, int d_freedom);
void initialiseParticle(int s, int D, Position *Xi, Velocity *Vi, double *xmin, double *xmax, double *Vmin, double *Vmax);
void swarmActivity(int S, int D, Position *Xi, int nb_eval, int outlev);

#endif
