/* evaluate_energy.cc */

    #include "structs.h"
    #include "evaluate_energy.h"


extern FILE *logFile;

float evaluate_energy( 

    float crd[MAX_ATOMS][SPACE],
    float charge[MAX_ATOMS],
    int   type[MAX_ATOMS],
    int   natom,
    float map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS],
    float inv_spacing,
    float xlo,
    float ylo,
    float zlo,
    int   nonbondlist[MAX_NONBONDS][2],
    float e_internal[NEINT][ATOM_MAPS][ATOM_MAPS],
    int   Nnb,
    Boole B_calcIntElec,
    float q1q2[MAX_NONBONDS],
    Boole B_isGaussTorCon,
    Boole B_isTorConstrained[MAX_TORS],
    State now,
    Boole B_ShowTorE,
    unsigned short US_TorE[MAX_TORS],
    unsigned short US_torProfile[MAX_TORS][NTORDIVS])

{
    float e = 0.;
    int   I_tor = 0;
    int   indx = 0;

    e = quicktrilinterp( crd, charge, type, natom, map,
			 	   inv_spacing, xlo, ylo, zlo);

    /* pr(logFile,"e(tril)=%10.2f,  ",e); / *###*/

/*    e += (*Addr_eintra = eintcal( nonbondlist, e_internal, crd, type, Nnb, B_calcIntElec, q1q2) ); */

    e += eintcal( nonbondlist, e_internal, crd, type, Nnb, B_calcIntElec, q1q2);

    /* pr(logFile,"e(eintcal)=%10.2f\n",e); / *###*/

    /*
    ** FASTER, LESS ACCURATE  METHOD:
    ** (alternative to above line).
    **
    ** Only bother to calculate internal energy if
    ** trilinterp energy is low enough.
    **
    ** if ( (e = quicktrilinterp( crd, charge, type, natom,
    ** map, inv_spacing, xlo,ylo,zlo)) <
    ** ENERGY_CUTOFF) {
    **   e += (*Addr_eintra = eintcal( nonbondlist,e_internal,
    **   crd,type,Nnb,B_calcIntElec,q1q2 ));
    ** }
    */

    if (B_isGaussTorCon) {
	for (I_tor = 0; I_tor <= now.ntor; I_tor++) {
	    if (B_isTorConstrained[I_tor] == 1) {
		indx = Rad2Div( now.tor[I_tor] );
		if (B_ShowTorE) {
		    e += (float)(US_TorE[I_tor] = US_torProfile[I_tor][indx]);
		} else {
		    e += (float)US_torProfile[I_tor][indx];
		}
	    }/*if*/
	}/*I_tor*/
    }/*if*/

    return e;
}
/* EOF */
