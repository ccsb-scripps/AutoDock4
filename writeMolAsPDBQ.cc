/*

 $Id: writeMolAsPDBQ.cc,v 1.1.6.1 2005/10/10 16:51:41 alther Exp $

*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/* writeMolAsPDBQ.cc */

#include <stdio.h>
#include <string.h>
#include "structs.h"
#include "stateLibrary.h"
#include "writeMolAsPDBQ.h"

extern int      keepresnum;
extern FILE    *logFile;

void
writeMolAsPDBQ(Molecule * mol, FILE * output)
/*
 * write out the Molecule "mol" as a PDBQ formatted file to the file pointed
 * to by "output".
 */
{
	int             i;
	char            rec14[14], rec9[9];

	(void) fprintf(output, "\n");
	(void) fprintf(output, "DIAGNOSTIC: MODEL     %4d\n", 1);

	(void) fprintf(output, "DIAGNOSTIC: USER    NEWDPF tran0 %f %f %f\n", mol->S.T.x, mol->S.T.y, mol->S.T.z);
	(void) fprintf(output, "DIAGNOSTIC: USER    NEWDPF quat0 %f %f %f %f\n", mol->S.Q.nx, mol->S.Q.ny, mol->S.Q.nz, ((mol->S.Q.ang) * 57.29577951));
	if (mol->S.ntor > 0) {
		(void) fprintf(output, "DIAGNOSTIC: USER    NEWDPF ndihe %d\n", mol->S.ntor);
		(void) fprintf(output, "DIAGNOSTIC: USER    NEWDPF dihe0 ");
		for (i = 0; i < mol->S.ntor; i++) {
			(void) fprintf(output, "%.2f ", ((mol->S.tor[i]) * 57.29577951));
		}
		(void) fprintf(output, "\n");
	}
	(void) fprintf(output, "DIAGNOSTIC: USER                              x       y       z     vdW  Elec       q\n");
	(void) fprintf(output, "DIAGNOSTIC: USER                           _______ _______ _______ _____ _____    ______\n");
	if (keepresnum > 0) {
		for (i = 0; i < mol->natom; i++) {
			strncpy(rec14, &(mol->atomstr[i][13]), (size_t) 13);
			(void) fprintf(output, FORMAT_PDBQ_ATOM_RESSTR, "DIAGNOSTIC: ", i + 1, rec14, mol->crd[i][0], mol->crd[i][1], mol->crd[i][2], 0.0, 0.0, 0.0);
			(void) fprintf(output, "\n");
		}
	} else {
		for (i = 0; i < mol->natom; i++) {
			strncpy(rec9, &(mol->atomstr[i][13]), (size_t) 8);
			(void) fprintf(output, FORMAT_PDBQ_ATOM_RESNUM, "DIAGNOSTIC: ", i + 1, rec9, 1, mol->crd[i][0], mol->crd[i][1], mol->crd[i][2], 0.0, 0.0, 0.0);
			(void) fprintf(output, "\n");
		}
	}
	(void) fprintf(output, "DIAGNOSTIC: TER\n");
	(void) fprintf(output, "DIAGNOSTIC: ENDMDL\n");
	(void) fprintf(output, "_______________________________________________________________________________\n\n");
	(void) fflush(output);
}

/* EOF */
