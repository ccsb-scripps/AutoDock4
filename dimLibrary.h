#ifndef COPYDIMENSION
#define COPYDIMENSION

//#include "main.h"
#include "structs.h"
#include "constants.h"
#include "alea.h"

//void initialiseState( State *S );
//void copyState( State *destination, State  source);

void copyDimension( State *S, Position R);
void copyState2Dimension( Position *R, State S);
//void copyState2Dimension( Position R, State *S);
void initialiseDimension(GridMapSetInfo * info, double *xmin, double *xmax, int D);
//void initialiseDimension(float xlo, float xhi, float ylo, float yhi, float zlo, float zhi, double *xmin, double *xmax, int D);
void initialiseParticle(int s, int D, Position *Xi, Velocity *Vi, double *xmin, double *xmax, double *Vmin, double *Vmax);
void swarmActivity(int S, int D, Position *Xi, int nb_eval, int outlev);

/*void printState( FILE *fp,
		 State state, 
		 int detail );

void writeState( FILE *fp, 
		 State state );

int checkState(State *D);

Molecule copyStateToMolecule(State *source, Molecule *mol);*/

#endif
