#ifndef TORSION
#define TORSION

void  torsion( State now,
               Real crd[MAX_ATOMS][SPACE], 
	       Real vt[MAX_TORS][SPACE], 
               int   tlist[MAX_TORS][MAX_ATOMS], 
               int   ntor );
#endif
