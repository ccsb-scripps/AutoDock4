/*
 * mdist.h
 *
 * Peter Reilly et al.
 *
 */

#include "autocomm.h"

void mdist();

enum {C,N,O,H,XX,P,S};
	
double mindist[ATOM_MAPS][ATOM_MAPS];
double maxdist[ATOM_MAPS][ATOM_MAPS];

void mdist() {

	register int i,j;

    // Zero all the mindist and maxdist elements.
	for (i=0; i<ATOM_MAPS; i++) {
		for (j=0; j<ATOM_MAPS; j++) {
			mindist[i][j] = 0.0L;
			maxdist[i][j] = 0.0L;
		}
	}

    // Set all the mindist and maxdist elements to the defaults for AutoDock versions 1 - 3...
	for (i=0; i<ATOM_MAPS; i++) {
		for (j=0; j<ATOM_MAPS; j++) {
			mindist[i][j] = 0.9;
			maxdist[i][j] = 2.1;
		}
	}

    /*
     * These values were contributed by Peter Reilly et al.:
     *
    mindist[C][C] = 1.32;
    maxdist[C][C] = 1.545;
    mindist[C][N] = 1.32;
    maxdist[C][N] = 1.39;
    mindist[C][O] = 1.20;
    maxdist[C][O] = 1.43;
    mindist[C][XX] = 1.07;
    maxdist[C][XX] = 1.15;
    mindist[C][S] = 1.80;
    maxdist[C][S] = 1.84;
    mindist[N][H] = 0.99;
    maxdist[N][H] = 1.03;
    mindist[N][XX] = 0.99;
    maxdist[N][XX] = 1.03;
    mindist[O][H] = 0.94;
    maxdist[O][H] = 0.98;
    mindist[O][P] = 1.47;
    maxdist[O][P] = 1.63;
    mindist[H][S] = 1.316;
    maxdist[H][S] = 1.356;
    mindist[XX][S] = 1.316;
    maxdist[XX][S] = 1.356;
    mindist[S][S] = 2.018;
    maxdist[S][S] = 2.058;
     */
    mindist[C][H] = 1.07;
    maxdist[C][H] = 1.15;
    mindist[H][C] = 1.07;
    maxdist[H][C] = 1.15;

    mindist[N][H] = 0.99;
    maxdist[N][H] = 1.10; // maxdist[N][H] = 1.03;
    mindist[H][N] = 0.99;
    maxdist[H][N] = 1.10; // maxdist[H][N] = 1.03;

    mindist[O][H] = 0.94;
    maxdist[O][H] = 1.10; // maxdist[O][H] = 0.98;
    mindist[H][O] = 0.94;
    maxdist[H][O] = 1.10; // maxdist[H][O] = 0.98;

    mindist[S][H] = 1.316;
    maxdist[S][H] = 1.356;
    mindist[H][S] = 1.316;
    maxdist[H][S] = 1.356;

}; 
