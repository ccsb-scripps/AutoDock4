#ifndef COPYSTATE
#define COPYSTATE

#include "structs.h"
#include "constants.h"

void initialiseState( State *S );

void initialiseQuat( Quat *Q );

void copyState( State *destination,
		State  source);

void printState( FILE *fp,
		 State state, 
		 int detail );

void writeState( FILE *fp, 
		 State state );

int checkState(State *D);

Molecule copyStateToMolecule(State *source, Molecule *mol);
#endif
