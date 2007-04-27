/*

 $Id: cmdmode.cc,v 1.20 2007/04/27 06:01:48 garrett Exp $

 AutoDock 

 Copyright (C) 1989-2007,  Garrett M. Morris, David S. Goodsell, Ruth Huey, Arthur J. Olson, 
 All Rights Reserved.

 AutoDock is a Trade Mark of The Scripps Research Institute.

 This program is free software; you can redistribute it and/or
 modify it under the terms of the GNU General Public License
 as published by the Free Software Foundation; either version 2
 of the License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

 */


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/times.h>
#include <time.h>
#include <string.h>
#include <ctype.h>
#include "cmdmode.h"
#include "cmdtokens.h"
#include "trjtokens.h"
#include "eintcal.h"


extern FILE *logFile;
extern FILE *command_in_fp;
extern FILE *command_out_fp;

extern char *programname;
extern char dock_param_fn[];
extern char AutoDockHelp[];

extern int ignore_errors;
extern int keepresnum;
extern int debug;
extern int parse_tors_mode;

int cmdmode(int   natom,
             Clock jobStart,
             struct tms tms_jobStart,
             Real map[MAX_GRID_PTS][MAX_GRID_PTS][MAX_GRID_PTS][MAX_MAPS],

                    EnergyTables *ptr_ad_energy_tables,

             Real WallEnergy,
             Real vt[MAX_TORS][SPACE],
             int   tlist[MAX_TORS][MAX_ATOMS],
             int   ntor,
             int   Nnb,
             NonbondParam *nonbondlist,
             char  atomstuff[MAX_ATOMS][MAX_CHARS],
             Real crdpdb[MAX_ATOMS][SPACE],
             char  hostnm[MAX_CHARS],
             int   type[MAX_ATOMS],
             Real charge[MAX_ATOMS],
             Real abs_charge[MAX_ATOMS],
             Real qsp_abs_charge[MAX_ATOMS],
             Boole B_calcIntElec,
             char  atm_typ_str[ATOM_MAPS],
             Real torsFreeEnergy,
             int ligand_is_inhibitor,
             int ignore_inter[MAX_ATOMS],
             const Boole         B_include_1_4_interactions,
             const Real scale_1_4,
             const ParameterEntry parameterArray[MAX_MAPS],
             const Real unbound_internal_FE,

             GridMapSetInfo *info,
             Boole B_have_flexible_residues,
             Boole B_use_non_bond_cutoff
            )

{
    char message[LINE_LEN],
         command[LINE_LEN],
         trjline[LINE_LEN],
         filename[MAX_CHARS],
         line[LINE_LEN],
         rec5[5],
         pdbaname[MAX_ATOMS][5],
         rec8[MAX_ATOMS][9], /* rec[X] gives (X+1) elements...   */
         rec14[MAX_ATOMS][15],
         trjFileName[MAX_CHARS],
         lastmove = '?';

    int  com_id = 0,
         nframes = 0,
         indpf = 0,
         offset[VECLENMAX],
         nat = 0,
         outside = FALSE,
         ntor_old = 0,
         keyword_id = -1,
         irun = -1,
         icycle = -1,
         nstep = 0,
         veclen = 0,
         movecode = -1;

    register int   i = 0,
         ii = 0;

    Clock  jobEnd;

    struct tms tms_jobEnd;

    Real eintra = 0.,
          einter = 0.,
          etotal = 0.,
          etot   = 0.,
          crd[MAX_ATOMS][SPACE],
          elec[MAX_ATOMS],
          emap[MAX_ATOMS],
          T = 0.,
          charge_total = 0.,
          emap_total = 0.,
          elec_total = 0.,
          E = 0.,
          Eint   = 0.;

    FILE *pdbFile, *trjFile;

    State S;

    ParameterEntry thisparm;

    EnergyBreakdown eb;

    for (i = 0;  i < natom;  i++) {
        strncpy(rec8[i], &atomstuff[i][13], (size_t)8);
        strncpy(rec14[i], &atomstuff[i][6], (size_t)14);
    }

/*
    Set up the command file-pointers to standard i/o,
*/
    set_cmd_io_std();

    pr_2x(logFile, stderr, UnderLine);
    prStr(message, "%s: NOW IN COMMAND MODE.\n", programname);
    pr_2x(logFile, stderr, message);
    fflush(logFile);
    fflush(stderr);

/*
    Now read in the Commands...
*/
    while ((fgets(command, LINE_LEN, command_in_fp)) != NULL) {
/* 
        Parse the command line...
*/
        com_id = parse_com_line(command);

        if (com_id == -1) {
            pr(stderr, "%s: ERROR: Unrecognized command: \"%c%c%c%c\"\n",
               programname, command[0], command[1], command[2], command[3]);
            continue;
        } /* endif */

/* 
        Act on that command...
*/
        switch(com_id) {
/*
            ____________________________________________________________________
*/
            case COM_NULL:
                pr(logFile, "%s\n", command);
                break;
/*
            ____________________________________________________________________
*/
            case COM_STOP:
/*
                Write out AVS field-file suitable for "animating" a
                trajectory file.
*/
                pr(logFile, "\n\nTotal number of atoms = %d\n", natom);
                pr(logFile, "Total number of trajectory frames = %d\n\n", nframes);

                indpf = strindex(dock_param_fn, ".dpf");
                strncpy(filename, dock_param_fn, (size_t)indpf);
                filename[ indpf ] = '\0';
                strcat(filename, ".dlg.pdb\0");

                offset[0]=5;
                offset[1]=6;
                offset[2]=7;
                offset[3]=4;
                offset[4]=1;

                print_avsfld(logFile, 5, natom, nframes, offset, 8,
                                  "x y z ResidueNum AtomNum", filename);

                pr(logFile, "COMMAND: stop\n\n");

                success(hostnm, jobStart, tms_jobStart);

                return 0;
/*
            ____________________________________________________________________
*/

            case COM_EPDB:
/*
                EPDB filename flag
                Return the energy of the Small Molecule.
                filename must be in PDBQ-format;
                flag can be:-
*/
                sscanf(command, "%*s %s", filename);
                pr(logFile, "COMMAND: epdb %s\n\n", filename);
 
                nat = 0;
                eintra = einter = etotal = 0.;
                outside = FALSE;
 
                if (openFile(filename, "r", &pdbFile, jobStart, tms_jobStart, FALSE)) {
                    while ((fgets(line, LINE_LEN, pdbFile)) != NULL) {
                        for (ii = 0; ii < 4; ii++) {
                            rec5[ii] = (char) tolower((int)line[ii]);
                        }
                        if (equal(rec5, "atom", 4) || equal(rec5, "heta", 4)) {
                            int serial;
                            readPDBQTLine(line, &serial, crd[nat], &charge[nat], &thisparm);
                            strncpy(pdbaname[natom], &line[12], (size_t)4);
                            type[nat]=get_atom_type(pdbaname[natom]);
                            if (type[nat] == -1) {
                                jobEnd = times(&tms_jobEnd);
                                timesys(jobEnd - jobStart, &tms_jobStart, &tms_jobEnd);
                                pr_2x(logFile, stderr, UnderLine);
                                return -1;
                            } /* endif */
                            outside = is_out_grid_info(crd[nat][X], crd[nat][Y], crd[nat][Z]);
                            nat++;
                            if (outside) {
                                (void) sprintf( message, "%s: WARNING: Atom %d is outside the grid!\n(%s)\n", programname, nat, line);
                                print_2x( logFile, stderr, message );
                                /* Reset outside */
                                outside = FALSE;
                            }
                        } /* endif atom or hetatm */
                    } /* endwhile */
                    fclose(pdbFile);
                    natom = nat;
                    if (ntor > 0) {
                        eintra = eintcalPrint(nonbondlist, ptr_ad_energy_tables, crd, Nnb, B_calcIntElec, B_include_1_4_interactions, scale_1_4, abs_charge, parameterArray, B_use_non_bond_cutoff, B_have_flexible_residues) - unbound_internal_FE;
                    } else {
                        eintra = 0.0 - unbound_internal_FE;
                    }
                    pr(logFile, "\n\n\t\tIntermolecular Energy Analysis\n");
                    pr(logFile,     "\t\t==============================\n\n\n");
                    pr(logFile, "Atom  NB.+ Elec.  Non-bonded  Electrosta  Partial          Coordinates         \n");
                    pr(logFile, "Type    Energy      Energy    tic Energy  Charge      x         y         z    \n");
                    pr(logFile, "____  __________  __________  __________  _______  ________  ________  ________\n");
                    /*            1234  0123456789  0123456789  0123456789  1234567  12345678  12345678  12345678"*/
                    /*           " ---  ----------  ----------  ----------  -------  --------  --------  --------"*/
                    emap_total = 0.;
                    elec_total = 0.;
                    charge_total = 0.;
                    for (i = 0;  i < natom;  i++) {
                        etot = emap[i] + elec[i];
                        pr(logFile, "%4d  %10.2f  %10.2f  %10.2f  %7.3f  %8.4f  %8.4f  %8.4f\n",
                                (type[i]+1), etot, emap[i], elec[i], charge[i], crd[i][X], crd[i][Y], crd[i][Z]);
                        emap_total += emap[i];
                        elec_total += elec[i];
                        charge_total += charge[i];
                    } /*i*/
                    pr(logFile, "      __________  __________  __________  _______\n\n");
                    pr(logFile, "Total %10.2f  %10.2f  %10.2f  %7.3f\n\n",
                        (emap_total + elec_total), emap_total, elec_total, charge_total);
                
                    pr(command_out_fp, "%.2f\n", etotal);
                    pr(logFile, "    E_intermolecular_atomic-affinity = %.2f kcal/mol\n", emap_total);
                    pr(logFile, "    E_intermolecular_electrostatic   = %.2f kcal/mol\n", elec_total);
                    printEnergies( &eb, "epdb: USER    ", ligand_is_inhibitor, emap_total, elec_total, B_have_flexible_residues);
                    pr(logFile, "\n");
                    fflush(logFile);
                    fflush(stderr);
                } /* endif */
                break;

/*
            ____________________________________________________________________
*/
            case COM_EVAL:
/*
                EVALUATE x y z  nx ny nz angle [in deg]
                tor1
                tor2
                ...
*/
                sscanf(command, "%*s %lf %lf %lf %lf %lf %lf %lf", 
                    &(S.T.x), &(S.T.y), &(S.T.z),  
                    &(S.Q.nx), &(S.Q.ny), &(S.Q.nz),  &(S.Q.ang));
                pr(logFile, "COMMAND: eval %lf %lf %lf %lf %lf %lf %lf\n         ", 
                    S.T.x, S.T.y, S.T.z, S.Q.nx, S.Q.ny, S.Q.nz, S.Q.ang);
                /**/
                S.Q.ang = DegreesToRadians(S.Q.ang);
                mkUnitQuat(&(S.Q));
                /**/
                for (i=0; i<ntor; i++) {
                    fgets(command, LINE_LEN, command_in_fp);
                    sscanf(command, "%lf", &(S.tor[i])); /* S.tor in degrees */
                    pr(logFile, "%lf ", S.tor[i]); /* S.tor in degrees */
                    S.tor[i] = DegreesToRadians(S.tor[i]); /* S.tor converted to radians */
                }
                cnv_state_to_coords(S,  vt, tlist, ntor,  crdpdb, crd, natom);
                if (ntor > 0) {
                    eintra = eintcalPrint(nonbondlist, ptr_ad_energy_tables, crd, Nnb, B_calcIntElec, B_include_1_4_interactions, scale_1_4, abs_charge, parameterArray, B_use_non_bond_cutoff, B_have_flexible_residues) - unbound_internal_FE;
                } else {
                    eintra = 0.0 - unbound_internal_FE;
                }
                outside = FALSE;
                for (i = 0;  i < natom;  i++) {
                    outside = is_out_grid_info(crd[i][X], crd[i][Y], crd[i][Z]);
                    if (outside) {
                        // break; // gmm 2001.11.07
                        (void) sprintf( message, "%s: WARNING: Atom %d is outside the grid!\n", programname, i+1);
                        print_2x( logFile, stderr, message );
                        /* Reset outside */
                        outside = FALSE;
                    }
                } /*i*/
                prStr(message, "%f\n", etotal);
                pr_2x(command_out_fp, logFile, message);
                fflush(logFile);
                fflush(stderr);
                break;
/*
            ____________________________________________________________________
*/
            case COM_OUTC:
                pr(logFile, "COMMAND: outc\n\n");
                pr(logFile, "USER    Total Non-bonded  Energy of Complex = %+8.2f\n", etotal);
                for (i = 0;  i < natom;  i++) {
                    pr(command_out_fp, "%f %f %f\n",  crd[i][X], crd[i][Y], crd[i][Z]);
                    pr(logFile, FORMAT_PDBQ_ATOM_RESNUM, "", i+1, rec8[i], nframes, crd[i][X], crd[i][Y], crd[i][Z], 0.0, 0.0, 0.0);
                    pr(logFile,"\n");
                } /*i*/
                strcpy(message, "TER\n\0");
                pr_2x(command_out_fp, logFile, message);
                fflush(logFile);
                fflush(stderr);
                nframes++;
                break;
/*
            ____________________________________________________________________
*/
            case COM_OUTE:
                pr(logFile, "COMMAND: oute\n\n");
                printEnergies( &eb, "oute: USER    ", ligand_is_inhibitor, emap_total, elec_total, B_have_flexible_residues);
/*                 prStr(message, "USER    Total Internal Energy of Small Molecule = %.2f\n", eintra); */
/*                 pr_2x(command_out_fp, logFile, message); */
/*                 prStr(message, "USER    Total Docked Energy of Complex = %.2f\n", etotal); */
/*                 pr_2x(command_out_fp, logFile, message); */
/*                 prStr(message, "USER    Predicted Free Energy of Binding = %.2f\n", einter + torsFreeEnergy); */
/*                 pr_2x(command_out_fp, logFile, message); */
                for (i = 0;  i < natom;  i++) {
                    pr(command_out_fp, "%.14s  %10.2f%10.2f\n", rec14[i], emap[i], elec[i]);
                    pr(logFile,        "%.14s  %10.2f%10.2f\n", rec14[i], emap[i], elec[i]);
                } /*i*/
                fflush(logFile);
                fflush(stderr);
                break;
/*
            ____________________________________________________________________
*/
            case COM_TRJ:
/*
                Convert a trajectory, from .trj file into PDB...
*/
                sscanf(command, "%*s %s", trjFileName);
                pr(logFile, "COMMAND: traj %s\n\n", trjFileName);
                if (!openFile(trjFileName, "r", &trjFile,jobStart,tms_jobStart,FALSE)) {
                    return -1;
                } else {
                    ntor_old = ntor;
                    T = 0.;

                    while ((fgets(trjline, LINE_LEN, trjFile)) != NULL) {
                        keyword_id = parse_trj_line(trjline);

                        if (keyword_id == -1) {
                            pr(stderr, "%s: ERROR: Unrecognized keyword in Trajectory file: \"%s\"\n",
                                programname, trjline);
                            continue;
                        }

                        switch(keyword_id) {
/*
                        ________________________________________________________
*/
                            case TRJ_NULL:
                                break;
/*
                        ________________________________________________________
*/
                            case TRJ_NTOR:
                                sscanf(trjline, "%*s %d", &ntor);
                                if (ntor != ntor_old) {
                                    prStr(message, "\n%s: mis-matched number of torsions in trajectory file %s - find correct docking-parameter file and substrate-PDBQ file.\n", programname, trjFileName);
                                    pr_2x(stderr, logFile, message);
                                }
                                break;
/*
                        ________________________________________________________
*/
                            case TRJ_RUN:
                                sscanf(trjline, "%*s %d", &irun);
                                break;
/*
                        ________________________________________________________
*/
                            case TRJ_CYCLE:
                                sscanf(trjline, "%*s %d", &icycle);
                                break;
/*
                        ________________________________________________________
*/
                            case TRJ_TEMP:
                                sscanf(trjline, "%*s " FDFMT, &T);
                                break;
/*
                        ________________________________________________________
*/
                            case TRJ_STATE:
                                if (input_state(&S, trjFile, trjline, ntor, 
                                    &nstep, &E, &Eint, &lastmove) != (int)0) {
                                    /*...input_state ensures tor is in radians*/
                                    cnv_state_to_coords(S, vt, tlist, ntor, crdpdb, crd, natom);
                                    outside = FALSE;
                                    for (i = 0;  i < natom;  i++) {
                                        outside = is_out_grid_info(crd[i][X], crd[i][Y], crd[i][Z]);
                                        if (outside) {
                                            // break; // gmm 2001.11.07
                                            (void) sprintf( message, "%s: WARNING: Atom %d is outside the grid!\n", programname, i+1);
                                            print_2x( logFile, stderr, message );
                                            /* Reset outside */
                                            outside = FALSE;
                                        }
                                    } /*i*/
                                    pr(logFile, "USER   Run %d  Cycle %d  Step %d  Temp %.2f K  %c  Etot %.2f  Eint   %.2f\n", irun, icycle, nstep, T, lastmove, E, Eint);
                                    switch (lastmove) {
                                        case 'A':
                                            movecode = 1;
                                            break;
                                        case 'a':
                                            movecode = 2;
                                            break;
                                        case 'e':
                                            movecode = 3;
                                            break;
                                        case 'R':
                                            movecode = 4;
                                            break;
                                        default:
                                            movecode = -1;
                                            break;
                                    } /* endswitch */
                                    for (i = 0;  i < natom;  i++) {
                                        /*
                                        Write coordinates and energy to command output...

                                        lastmove = AaeR,
                                        movecode = 1234
                                        */
                                        pr(command_out_fp, "%.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.2f %d %d %d %d\n", crd[i][X], crd[i][Y], crd[i][Z],  emap[i], elec[i], emap[i]+elec[i],  Eint, etotal, T, movecode, nstep, icycle, irun);
                                        pr(logFile, FORMAT_PDBQ_ATOM_RESNUM, "", i+1, rec8[i], nstep, crd[i][X], crd[i][Y], crd[i][Z], 0.0, 0.0, 0.0);
                                        pr(logFile,"\n");
                                    } /*i*/
                                    pr(logFile, "TER\n");
                                    fflush(logFile);
                                    fflush(stderr);
                                    nframes++;
                                } /* endif inputstate */
                                break;
/*
**                        ________________________________________________________
*/
                           default:
                                pr(stderr, "%s: ERROR: Unrecognized keyword: '%s'\n", programname, trjline);
                                break;
/*
**                        ________________________________________________________
*/
                        }

                    } /*  end of while -- EOF of trjFile  */
                    fclose(trjFile);
                }
                break;
/*
**            ____________________________________________________________________
*/
            default:
                pr(stderr, "%s: ERROR: Unrecognized command: '%s'\n", programname, line);
                break;
/*
**            ____________________________________________________________________
*/
        } /* End of switch. */
    } /* End of while. */
/*

Write out an AVS-readable field file, for input of trajectory file.

*/
    indpf = strindex(dock_param_fn, ".dpf");
    strncpy(filename, dock_param_fn, (size_t)indpf);
    filename[ indpf ] = '\0';
    strcat(filename, ".dlg.pdb\0");

    for (i=0, veclen = 13; i<veclen; i++) {
        offset[i] = i;
    }

    print_avsfld(logFile, veclen, natom, nframes, offset, 13, 
        "x y z E_atom E_elec E_atom_elec E_internal E_total Temp MoveCode Step Cycle Run", filename);

    jobEnd = times(&tms_jobEnd);
    timesys(jobEnd - jobStart, &tms_jobStart, &tms_jobEnd);
    pr_2x(logFile, stderr, UnderLine);

    return 0;
}
/* EOF */
