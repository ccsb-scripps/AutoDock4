/*

 $Id: writePDBQ.cc,v 1.12 2006/04/25 22:33:33 garrett Exp $

*/

/* writePDBQ.cc */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>
#include <stdio.h>
#include <string.h>
#include "structs.h"
#include "stateLibrary.h"
#include "writePDBQ.h"

extern int      keepresnum;
extern FILE    *logFile;

extern int      write_stateFile;
extern FILE    *stateFile;

void
writePDBQ(int irun,
	  char smFileName[MAX_CHARS],
	  char dpfFN[MAX_CHARS],
	  Real sml_center[SPACE],
	  State state,
	  int ntor,
	  Real eintra,
	  Real einter,
	  int natom,
	  char atomstuff[MAX_ATOMS][MAX_CHARS],
	  Real crd[MAX_ATOMS][SPACE],
	  Real emap[MAX_ATOMS],
	  Real elec[MAX_ATOMS],
	  Real charge[MAX_ATOMS],
	  Real abs_charge[MAX_ATOMS],
	  Real qsp_abs_charge[MAX_ATOMS],
	  int ligand_is_inhibitor,
	  Real torsFreeEnergy,
	  int outlev,
	  int ignore_inter[MAX_ATOMS],
	  const Boole B_include_1_4_interactions,
	  const Real scale_1_4,

	  const ParameterEntry parameterArray[MAX_MAPS],
	  const Real unbound_internal_FE)
{
	int             i = 0;

	Real   emap_total = 0.0L;
	Real   elec_total = 0.0L;
	Real   MaxValue = 99.99L;

	char            rec14[14], rec9[9];

	if (outlev > -1) {
		//output of coordinates // gmm 2001 - 11 - 01

        pr(logFile, "\n");
		pr(logFile, "DOCKED: MODEL     %4d\n", irun + 1);
		pr(logFile, "DOCKED: USER    Run = %d\n", irun + 1);
		pr(logFile, "DOCKED: USER    DPF = %s\n", dpfFN);

		printEnergies(einter, eintra, torsFreeEnergy, "DOCKED: USER    ", ligand_is_inhibitor, emap_total, elec_total, unbound_internal_FE);

		(void) fprintf(logFile, "DOCKED: USER    NEWDPF move %s\n", smFileName);
		(void) fprintf(logFile, "DOCKED: USER    NEWDPF about %f %f %f\n", sml_center[X], sml_center[Y], sml_center[Z]);
		(void) fprintf(logFile, "DOCKED: USER    NEWDPF tran0 %f %f %f\n", state.T.x, state.T.y, state.T.z);
		(void) fprintf(logFile, "DOCKED: USER    NEWDPF quat0 %f %f %f %f\n", state.Q.nx, state.Q.ny, state.Q.nz, Deg(WrpRad(ModRad(state.Q.ang))));
		if (ntor > 0) {
			(void) fprintf(logFile, "DOCKED: USER    NEWDPF ndihe %d\n", ntor);
			(void) fprintf(logFile, "DOCKED: USER    NEWDPF dihe0 ");
			for (i = 0; i < ntor; i++) {
				(void) fprintf(logFile, "%.2f ", Deg(state.tor[i]));
			}
			(void) fprintf(logFile, "\n");

		}
		/* write state file */
		if (write_stateFile) {
			pr(stateFile, "\t\t<move>%s</move>\n", smFileName);
			pr(stateFile, "\t\t<about>%f %f %f</about>\n", sml_center[X], sml_center[Y], sml_center[Z]);

			pr(stateFile, "\t\t<tran0>%f %f %f</tran0>\n", state.T.x, state.T.y, state.T.z);
			pr(stateFile, "\t\t<quat0>%f %f %f %f</quat0>\n", state.Q.nx, state.Q.ny, state.Q.nz, Deg(WrpRad(ModRad(state.Q.ang))));

			if (ntor > 0) {
				pr(stateFile, "\t\t<ndihe>%d</ndihe>\n", ntor);
				pr(stateFile, "\t\t<dihe0>");
				for (i = 0; i < ntor; i++) {
					(void) fprintf(stateFile, "%.2f ", Deg(WrpRad(ModRad(state.tor[i]))));
				}
				(void) fprintf(stateFile, "\n");
				pr(stateFile, "</dihe0>\n");
			}
			pr(stateFile, "\t</run>\n");
		}		/* end write state file */
		(void) fprintf(logFile, "DOCKED: USER                              x       y       z     vdW  Elec       q\n");
		(void) fprintf(logFile, "DOCKED: USER                           _______ _______ _______ _____ _____    ______\n");
		if (keepresnum > 0) {
			for (i = 0; i < natom; i++) {
				strncpy(rec14, &atomstuff[i][13], (size_t) 13);
				(void) fprintf(logFile, FORMAT_PDBQ_ATOM_RESSTR, "DOCKED: ", i + 1, rec14, crd[i][X], crd[i][Y], crd[i][Z], min(emap[i], MaxValue), min(elec[i], MaxValue), charge[i]);
				(void) fprintf(logFile, "\n");
			}
		} else {
			for (i = 0; i < natom; i++) {
				strncpy(rec9, &atomstuff[i][13], (size_t) 8);
				(void) fprintf(logFile, FORMAT_PDBQ_ATOM_RESNUM, "DOCKED: ", i + 1, rec9, irun + 1, crd[i][X], crd[i][Y], crd[i][Z], min(emap[i], MaxValue), min(elec[i], MaxValue), charge[i]);
				(void) fprintf(logFile, "\n");
			}
		}
		(void) fprintf(logFile, "DOCKED: TER\n");
		(void) fprintf(logFile, "DOCKED: ENDMDL\n");
		//(void) fprintf(logFile, UnderLine);
		(void) fflush(logFile);
	}
} //end of function, void
/* EOF */
