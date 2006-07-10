/*

 $Id: prClusterHist.cc,v 1.5 2006/07/10 22:51:46 garrett Exp $

*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/* prClusterHist.cc */

#include <math.h>
#include <stdio.h>
#include "prClusterHist.h"
#include "constants.h"


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
    Real          etot = 0.,
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

    double U_internal_energy = 0.0;
    double U_internal_energy_summation = 0.0;
    double Q_partition_function = 0.0;
    double A_free_energy = 0.0;
    double S_entropy = 0.0;
    double this_energy = 0.0;
    double RT = Rcal * TK;

    // Loop over all the clusters
    for (Rank = 0;  Rank < ncluster;  Rank++) {
        Rank1 = Rank + 1;

        // Print the Cluster Rank
        (void)fprintf( logFile, "\n%4d |", Rank1);

        // Print the Lowest Binding Energy
        if (econf[cluster[Rank][0]] > 999999.99) {
            (void)fprintf( logFile, "%+10.2e", econf[cluster[Rank][0]]);
        } else {
            (void)fprintf( logFile, "%+10.2f", econf[cluster[Rank][0]]);
        }

        // Print the Run
        (void)fprintf( logFile, " |%4d |", cluster[Rank][0]+1);

        // Print the Mean Binding Energy
        if (num_in_clu[Rank] > 1) {
            /* Calculate average energy in cluster */
            etot = 0.;
            for (j = 0;  j < num_in_clu[Rank]; j++ ) {
                this_energy = econf[ cluster[Rank][j] ];

                U_internal_energy_summation += this_energy * exp( -this_energy / RT );
                Q_partition_function += exp( -this_energy / RT );

                etot += this_energy;
            }
            eavg = etot / (Real)num_in_clu[Rank];
            num_multi_mem_clu++;
            if (eavg > 999999.99) {
                (void)fprintf( logFile, "%+10.2e |", eavg );
            } else {
                (void)fprintf( logFile, "%+10.2f |", eavg );
            }
        } else {
            this_energy = econf[ cluster[Rank][0] ];

            U_internal_energy_summation += this_energy * exp( -this_energy / RT );
            Q_partition_function += exp( -this_energy / RT );

            if (this_energy > 999999.99) {
                (void)fprintf( logFile, "%+10.2e |", this_energy);
            } else {
                (void)fprintf( logFile, "%+10.2f |", this_energy);
            }
        }

        // Print the Number in Cluster
        (void)fprintf( logFile, "%4d |", num_in_clu[Rank] );

        // Print the Histogram
        // Print a '#' symbol for each member in this cluster
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

    // Finish the calculation of the internal energy, U, which depends on the partition function, Q:
    U_internal_energy = (1. / Q_partition_function) * U_internal_energy_summation;
    // Calculate the free energy, A
    A_free_energy = -RT * log( Q_partition_function );
    // Calculate the entropy, S
    S_entropy = -( A_free_energy - U_internal_energy ) / TK;

    // Print the Statistical Thermodyamics values 
    //
    (void)fprintf( logFile, "\n\n\tSTATISTICAL MECHANICAL ANALYSIS\n" );
    (void)fprintf( logFile, "\t_______________________________\n" );
    (void)fprintf( logFile, "\n\n" );
    (void)fprintf( logFile, "Internal (Boltzmann-weighted) energy, U = %.2f kcal/mol at Temperature, T = %.2f K\n", U_internal_energy, TK );
    (void)fprintf( logFile, "Partition function, Q = %.2f kcal/mol at Temperature, T = %.2f K\n", Q_partition_function, TK );
    (void)fprintf( logFile, "Free energy, A = %.2f kcal/mol at Temperature, T = %.2f K\n", A_free_energy, TK );
    (void)fprintf( logFile, "Entropy, S = %.2f kcal/mol/K at Temperature, T = %.2f K\n", S_entropy, TK );
    (void)fprintf( logFile, "\n" );
    (void)fprintf( logFile, "_______________________________________________________________________\n\n");

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
