/*

 $Id: prClusterHist.cc,v 1.4 2006/04/25 22:32:44 garrett Exp $

*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/* prClusterHist.cc */

    #include <stdio.h>
    #include "prClusterHist.h"


extern FILE *logFile;

void prClusterHist( int ncluster,
		    int irunmax,
		    Real clus_rms_tol,
		    int num_in_clu[MAX_RUNS],
		    int cluster[MAX_RUNS][MAX_RUNS],
		    Real econf[MAX_RUNS],
		    Real clu_rms[MAX_RUNS][MAX_RUNS],
		    Real ref_rms[MAX_RUNS])
{
    Real         etot = 0.,
		  eavg = 0.;

    register int  Rank=0,
	          j=0;

    int           num_multi_mem_clu=0,
		  ConfNum=0,
		  Rank1=1;

    (void)fprintf( logFile, "\n\n" );
    (void)fprintf( logFile, "________________________________________________________________________________\n\n" );

    (void)fprintf( logFile, "Number of distinct conformational clusters found = %d,  out of %d runs,\nUsing an rmsd-tolerance of %.1f A\n\n\n", ncluster, irunmax, (double)clus_rms_tol );
 
    (void)fprintf( logFile, "\tCLUSTERING HISTOGRAM\n" );
    (void)fprintf( logFile, "\t____________________\n\n");

    (void)fprintf( logFile, "________________________________________________________________________________\n");
    (void)fprintf( logFile, "     |           |     |           |     |                                    \n" );
    (void)fprintf( logFile, "Clus | Lowest    | Run | Mean      | Num | Histogram                          \n" );
    (void)fprintf( logFile, "-ter | Binding   |     | Binding   | in  |                                    \n" );
    (void)fprintf( logFile, "Rank | Energy    |     | Energy    | Clus|    5    10   15   20   25   30   35\n" );
    (void)fprintf( logFile, "_____|___________|_____|___________|_____|____:____|____:____|____:____|____:___");

    for (Rank = 0;  Rank < ncluster;  Rank++) {

	Rank1 = Rank + 1;

	(void)fprintf( logFile, "\n%4d |", Rank1);
	if (econf[cluster[Rank][0]] > 999999.99) {
	    (void)fprintf( logFile, "%+10.2e", econf[cluster[Rank][0]]);
	} else {
	    (void)fprintf( logFile, "%+10.2f", econf[cluster[Rank][0]]);
	}
	(void)fprintf( logFile, " |%4d |", cluster[Rank][0]+1);

        if (num_in_clu[Rank] > 1) {
	    /* Calculate average energy in cluster */
	    etot = 0.;
	    for (j = 0;  j < num_in_clu[Rank]; j++ ) {
		etot += econf[ cluster[Rank][j] ];
	    }
	    eavg = etot / (Real)num_in_clu[Rank];
            num_multi_mem_clu++;
	    if (eavg > 999999.99) {
		(void)fprintf( logFile, "%+10.2e |", eavg );
	    } else {
		(void)fprintf( logFile, "%+10.2f |", eavg );
	    }
	} else {
	    if (econf[cluster[Rank][0]] > 999999.99) {
		(void)fprintf( logFile, "%+10.2e |", econf[cluster[Rank][0]]);
	    } else {
		(void)fprintf( logFile, "%+10.2f |", econf[cluster[Rank][0]]);
	    }
	}
	(void)fprintf( logFile, "%4d |", num_in_clu[Rank] );

        for (j=0;  j<num_in_clu[Rank]; j++) {
	    (void)fprintf( logFile, "#" );
	}/*j*/

    }/*Rank*/

    (void)fprintf( logFile, "\n" );
    (void)fprintf( logFile, "_____|___________|_____|___________|_____|______________________________________\n");
    (void)fprintf( logFile, "\n" );

    if (num_multi_mem_clu > 0) {
        (void)fprintf( logFile, "\nNumber of multi-member conformational clusters found = %d, out of %d runs.\n\n", num_multi_mem_clu, irunmax );
    }

    (void)fprintf( logFile, "\n\n" );
    (void)fprintf( logFile, "\tRMSD TABLE\n" );
    (void)fprintf( logFile, "\t__________\n\n");

    (void)fprintf( logFile, "_____________________________________________________________________\n");
    (void)fprintf( logFile, "     |      |      |           |         |                 |\n");
    (void)fprintf( logFile, "Rank | Sub- | Run  | Binding   | Cluster | Reference       | Grep\n");
    (void)fprintf( logFile, "     | Rank |      | Energy    | RMSD    | RMSD            | Pattern\n");
    (void)fprintf( logFile, "_____|______|______|___________|_________|_________________|___________\n" );

    for (Rank = 0;  Rank < ncluster;  Rank++) {
	Rank1 = Rank + 1;
	ConfNum = 0;
        for (j=0;  j<num_in_clu[Rank]; j++) {
	    ++ConfNum;
	    if (econf[cluster[Rank][j]] > 999999.99) {
		(void)fprintf( logFile, "%4d   %4d   %4d  %+10.2e  %8.2f  %8.2f           RANKING\n", Rank1, ConfNum,  cluster[Rank][j]+1, econf[cluster[Rank][j]], ((j==0)?(0.):(clu_rms[Rank][j])), ref_rms[cluster[Rank][j]] );
	    } else {
		(void)fprintf( logFile, "%4d   %4d   %4d  %+10.2f  %8.2f  %8.2f           RANKING\n", Rank1, ConfNum,  cluster[Rank][j]+1, econf[cluster[Rank][j]], ((j==0)?(0.):(clu_rms[Rank][j])), ref_rms[cluster[Rank][j]] );
	    }
        }   /*j*/
        //(void)fprintf( logFile, ".......................................................................\n" );
    }/*Rank*/
    (void)fprintf( logFile, "_______________________________________________________________________\n\n");

    //(void)fprintf( logFile, "________________________________________________________________________________\n\n" );

    fflush( logFile );
}
/* EOF */
	/*
	** kend = num_in_clu[Rank] - 1;
        ** (void)fprintf( logFile, "________________________________________________________________________________\n" );
        ** (void)fprintf( logFile, "ClusterRank=%d\t", Rank1 );
        ** (void)fprintf( logFile, "NumberOfMembers= %d:\n", num_in_clu[Rank] );
	*/
	/*
        ** for (j=0;  j<num_in_clu[Rank]; j+=OUTNUMCLUST) {
            ** kmax = min(num_in_clu[Rank],(j + OUTNUMCLUST));
            ** (void)fprintf( logFile, "\nClusterRank=%d\tNum=         \t", Rank1 );
            ** for (k = j; k < kmax;  k++ ) {
                ** (void)fprintf( logFile, "%8d%s", ConfNum, (k == kend)?";":"," );
		** ++ConfNum;
            ** }
            ** (void)fprintf( logFile, "\nClusterRank=%d\tRun=         \t", Rank1 );
            ** for (k = j; k < kmax;  k++ ) {
                ** (void)fprintf( logFile, "%8d%s", 1+cluster[Rank][k], (k == kend)?";":"," );
            ** }
            ** (void)fprintf( logFile, "\nClusterRank=%d\tEnergy=      \t", Rank1 );
            ** for (k = j; k < kmax; k++ ) {
                ** (void)fprintf( logFile, "%+8.2f%s", (double)econf[cluster[Rank][k]], (k == kend)?";":","  );
            ** }
            ** (void)fprintf( logFile, "\nClusterRank=%d\tClusterRMS=  \t", Rank1 );
            ** for (k = j; k < kmax; k++ ) {
                ** (void)fprintf( logFile, "%8.2f%s", (double)((k==0)?(0.):(clu_rms[Rank][k])), (k == kend)?";":","  );
            ** }
            ** (void)fprintf( logFile, "\nClusterRank=%d\tReferenceRMS=\t", Rank1 );
            ** for (k = j; k < kmax; k++ ) {
                ** (void)fprintf( logFile, "%8.2f%s", (double)ref_rms[cluster[Rank][k]], (k == kend)?";":","  );
            ** }
            ** (void)fprintf( logFile, "\n" );
	    */
