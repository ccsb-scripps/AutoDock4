#ifndef CHANGESTATE
#define CHANGESTATE

State  changeState( State last,      /* ...must be a normalized quaternion! */
		    float trnStep,
		    float torStep,
		    int   ntor,
		    float F_TorConRange[MAX_TORS][MAX_TOR_CON][2],
		    int   N_con[MAX_TORS]);
#endif
