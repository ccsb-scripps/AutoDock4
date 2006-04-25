/*

 $Id: cluster_analysis.cc,v 1.3 2006/04/25 22:31:55 garrett Exp $

*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/* cluster_analysis.cc */


#include <math.h>

    #include "cluster_analysis.h"


int cluster_analysis( Real clus_rms_tol, 
		      int cluster[MAX_RUNS][MAX_RUNS], 
		      int num_in_clus[MAX_RUNS], 
		      int isort[MAX_RUNS], 
		      int nconf, 
		      int natom, 
		      int type[MAX_ATOMS],
		      Real crd[MAX_RUNS][MAX_ATOMS][SPACE], 
		      Real crdpdb[MAX_ATOMS][SPACE], 
		      Real sml_center[SPACE], 
		      Real clu_rms[MAX_RUNS][MAX_RUNS], 
		      Boole B_symmetry_flag,
		      Real ref_crds[MAX_ATOMS][SPACE],
		      int ref_natoms,
		      Real ref_rms[MAX_RUNS])
{
/* __________________________________________________________________________
  | Cluster Analysis                                                         |
  |__________________________________________________________________________|
  |  If conformations are within clus_rms_tol, scored as the same.           |
  |  Compares atoms with same type, not atom name, to find propellers.       |
  |__________________________________________________________________________|*/

    register int compare = 0,
                 i = 0;

    int	  nClusters = 1,
	  thisconf = 0,
	  new_conf = FALSE;

    Real rms = 0.;

    if (ref_natoms == -1) {

/* No reference coordinates were defined, so we must
   un-center the original PDBQ coordinates, */

	for (i = 0;  i < natom;  i++) {
		ref_crds[i][0] = sml_center[0] + crdpdb[i][0];
		ref_crds[i][1] = sml_center[1] + crdpdb[i][1];
		ref_crds[i][2] = sml_center[2] + crdpdb[i][2];
	}/*i*/
    }

/* Assign the index of the lowest energy to 0,0 in "cluster" */

    thisconf = cluster[0][0] = isort[0]; 

/* Set number in 0-th cluster to 1 */
/* Also initialize total number of clusters to 1 */

    num_in_clus[0] = nClusters = 1;	 
    clu_rms[0][0]  = getrms(crd[thisconf], ref_crds,
			    B_symmetry_flag, natom, type);

/* Go through *all* conformations... */
    for ( i=0; i<nconf; i++) {

/* Calculate the RMSD to the reference structure: */
	ref_rms[i] = getrms(crd[i], ref_crds, B_symmetry_flag, natom, type);
    }

/* Go through all conformations *except* 0-th... */
    for ( i=1; i<nconf; i++) {	

/* "thisconf" is the index of the energy-sorted i-th conformation: */
        thisconf = isort[i];

/* Assume this is a new conformation until proven otherwise... */
        new_conf = TRUE;

	for ( compare=0; compare<nClusters; compare++ ) {

	    rms = getrms(crd[thisconf], crd[cluster[compare][0]],
			 B_symmetry_flag,natom,type);

/* Check rms; if greater than tolerance, */
            if ( rms > clus_rms_tol ) {

                continue; /* to compare next conformation, */

            } else {

/* Otherwise, add this conformation to the current cluster, */
		cluster[compare][ num_in_clus[compare] ] = thisconf;
		clu_rms[compare][ num_in_clus[compare] ] = rms;
		num_in_clus[compare]++;
		new_conf = FALSE;
		break; /* Go on to next i iteration... */
	    }
        } /* next comparison... */

        if (new_conf) {

/* Start a new cluster...  */
            cluster[ nClusters ][0] = thisconf;
	    clu_rms[ nClusters ][0] = getrms(crd[thisconf], ref_crds,
					     B_symmetry_flag,natom,type);
            num_in_clus[ nClusters ] = 1;

/* Increment the number of clusters... */
            nClusters++;

        } /* endif */

    } /* next i... */

    return nClusters;
}
/* EOF */
