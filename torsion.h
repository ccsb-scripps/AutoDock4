#ifndef TORSION
#define TORSION

void  torsion( State now,
               FloatOrDouble crd[MAX_ATOMS][SPACE], 
	       FloatOrDouble vt[MAX_TORS][SPACE], 
               int   tlist[MAX_TORS][MAX_ATOMS], 
               int   ntor );
#endif
