#ifndef TORSION
#define TORSION

void  torsion( const State now,
               Real crd[MAX_ATOMS][SPACE], 
               const Real vt[MAX_TORS][SPACE], 
               const int tlist[MAX_TORS][MAX_ATOMS], 
               const int ntor );
#endif
