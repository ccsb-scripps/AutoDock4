#ifndef CHANGESTATE
#define CHANGESTATE

State  changeState( State last,      /* ...must be a normalized quaternion! */
		    FloatOrDouble trnStep,
		    FloatOrDouble torStep,
		    int   ntor,
		    FloatOrDouble F_TorConRange[MAX_TORS][MAX_TOR_CON][2],
		    int   N_con[MAX_TORS]);
#endif
