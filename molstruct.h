/* structure for molecules; side-chains, etc. */

typedef struct molecule {

	int   natom;
	double untrnfm_crdpdb[ MAX_ATOMS ][ NTRN ];
	double trnsfmd_crdpdb[ MAX_ATOMS ][ NTRN ];
	double charge[ MAX_ATOMS ];
	int   type[ MAX_ATOMS ];
	char  pdbaname[ MAX_ATOMS ][ 5 ];
	char  atomstuff[ MAX_ATOMS ][ MAX_CHARS ];

	Boole B_haveCharges;
	char  pdbqFileName[ MAX_CHARS ];
	int   Htype;
	Boole B_constrain;
	int   atomC1;
	int   atomC2;
	double sqlower;
	double squpper;

	int   ntor1;
	int   ntor;
	int   tlist[ MAX_TORS ][ MAX_ATOMS ];
	double vt[ MAX_TORS ][ NTRN ];

	int   Nnb;
	int   Nnbonds[ MAX_ATOMS ];
	int   nonbondlist[MAX_NONBONDS][MAX_NBDATA];
} Molecule;
