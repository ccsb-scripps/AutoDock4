#ifndef TORSION
#define TORSION

void  torsion( State now,
               float crd[MAX_ATOMS][SPACE], 
	       float vt[MAX_TORS][SPACE], 
               int   tlist[MAX_TORS][MAX_ATOMS], 
               int   ntor );
#endif
