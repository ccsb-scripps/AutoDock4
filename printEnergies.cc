/*

 $Id: printEnergies.cc,v 1.10 2006/04/17 06:04:08 garrett Exp $

*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/* printEnergies.cc */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "printEnergies.h"

extern FILE *logFile;
extern FILE *stateFile;
extern int write_stateFile;

#define print1000(file, x) pr(file,  ((fabs((x)) >= 0.0) && ((fabs(x)) <= 1000.)) ? "%+7.2f" : "%+11.2e" , (x));

void print1000_no_sign(FILE* file, double x) {
    pr(file,  ((fabs((x)) >= 0.01) && ((fabs(x)) <= 1000.)) ? "%7.2f" : "%11.2e" , (x));
}

void print_molar(FILE* file, double x) {
    // 1e-3  <= x < 1     mM millimolar
    // 1e-6  <= x < 1e-3  uM micromolar
    // 1e-9  <= x < 1e-6  nM nanomolar
    // 1e-12 <= x < 1e-9  pM picomolar
    // 1e-15 <= x < 1e-12 fM femtomolar
    // 1e-18 <= x < 1e-15 aM attomolar
    // 1e-21 <= x < 1e-18 zM zeptomolar
    // 1e-24 <= x < 1e-21 yM yottomolar
    //          x < 1e-24    sub-yottomolar
    if ((fabs((x)) > 1e-3) && ((fabs(x)) <= 1.)) {
        pr(file, "%7.2f mM (millimolar)", x*1e3);
    } else if ((fabs((x)) > 1e-6) && ((fabs(x)) <= 1e-3)) {
        pr(file, "%7.2f uM (micromolar)", x*1e6);
    } else if ((fabs((x)) > 1e-9) && ((fabs(x)) <= 1e-6)) {
        pr(file, "%7.2f nM (nanomolar)", x*1e9);
    } else if ((fabs((x)) > 1e-12) && ((fabs(x)) <= 1e-9)) {
        pr(file, "%7.2f pM (picomolar)", x*1e12);
    } else if ((fabs((x)) > 1e-15) && ((fabs(x)) <= 1e-12)) {
        pr(file, "%7.2f fM (femtomolar)", x*1e15);
    } else if ((fabs((x)) > 1e-18) && ((fabs(x)) <= 1e-15)) {
        pr(file, "%7.2f aM (attomolar)", x*1e18);
    } else if ((fabs((x)) > 1e-21) && ((fabs(x)) <= 1e-18)) {
        pr(file, "%7.2f zM (zeptomolar)", x*1e21);
    } else if ((fabs((x)) > 1e-24) && ((fabs(x)) <= 1e-21)) {
        pr(file, "%7.2f yM (yottomolar)", x*1e24);
    } else {
        pr(file, "%11.2e M (molar)", x);
    }
}

void printEnergies( 
        FloatOrDouble einter,
        FloatOrDouble eintra,
        FloatOrDouble torsFreeEnergy,
        char  *prefixString,
        int ligand_is_inhibitor,
        FloatOrDouble emap_total,
        FloatOrDouble elec_total,
        FloatOrDouble unbound_internal_FE
        )
{
    FloatOrDouble deltaG = 0.0;
    FloatOrDouble Ki = 1.0;
    // FloatOrDouble RJ = 8.31441;  // in J/K/mol, Gas Constant, Atkins Phys.Chem., 2/e
    FloatOrDouble Rcal = 1.9871917; // in cal/K/mol, Gas Constant, RJ/4.184
    FloatOrDouble TK = 298.15;      // Room temperature, in K

    // equilibrium:   E  +  I  <=>    EI
    // binding:       E  +  I   ->    EI         K(binding),      Kb
    // dissociation:     EI     ->  E  +  I      K(dissociation), Kd
    //
    //                            1
    //         K(binding) = ---------------
    //                      K(dissociation)
    // so:
    //      ln K(binding) = -ln K(dissociation)
    //              ln Kb = -ln Kd
    // Ki = dissociation constant of the enzyme-inhibitor complex = Kd
    //      [E][I]
    // Ki = ------
    //       [EI]
    // so:
    //              ln Kb = -ln Ki
    // deltaG(binding)    = -R*T*ln Kb
    // deltaG(inhibition) =  R*T*ln Ki
    //
    // Binding and Inhibition occur in opposite directions, so we 
    // lose the minus-sign:  deltaG = R*T*lnKi,  _not_ -R*T*lnKi
    // => deltaG/(R*T) = lnKi
    // => Ki = exp(deltaG/(R*T))
    deltaG = einter + eintra + torsFreeEnergy - unbound_internal_FE;
    if (deltaG < 0.0) {
        Ki = exp((deltaG*1000.)/(Rcal*TK));
    }

    if (strncmp(prefixString, "UNBOUND", 7) != 0 ) {
        pr( logFile, "%sEstimated Free Energy of Binding    = ", prefixString);
        print1000(logFile, deltaG);
        pr( logFile, " kcal/mol  [=(1)+(2)+(3)-(4)]\n");

        if (deltaG < 0.0) {
            if (ligand_is_inhibitor == 1) {
                pr( logFile, "%sEstimated Inhibition Constant, Ki   = ", prefixString);
            } else {
                pr( logFile, "%sEstimated Dissociation Constant, Kd = ", prefixString);
            }
            // print1000_no_sign(logFile, Ki);
            print_molar(logFile, Ki);
            pr( logFile, "  [Temperature = %.2f K]\n", TK);
        }

        pr( logFile, "%s\n", prefixString);
    }

    pr( logFile, "%s(1) Final Intermolecular Energy     = ", prefixString);
    print1000(logFile, einter);
    pr( logFile, " kcal/mol\n");
    pr( logFile, "%s    vdW + Hbond + desolv Energy     = ", prefixString);
    print1000(logFile, emap_total);
    pr( logFile, " kcal/mol\n");
    pr( logFile, "%s    Electrostatic Energy            = ", prefixString);
    print1000(logFile, elec_total);
    pr( logFile, " kcal/mol\n");

    pr( logFile, "%s(2) Final Internal Energy           = ", prefixString);
    print1000(logFile, eintra);
    pr( logFile, " kcal/mol\n");

    pr( logFile, "%s(3) Torsional Free Energy           = ", prefixString);
    print1000(logFile, torsFreeEnergy);
    pr( logFile, " kcal/mol\n");

    pr( logFile, "%s(4) Unbound System's Energy         = ", prefixString);
    print1000(logFile, unbound_internal_FE);
    pr( logFile, " kcal/mol\n");

    pr( logFile, "%s\n", prefixString);
    pr( logFile, "%s\n", prefixString);
}

void printStateEnergies( FloatOrDouble einter, FloatOrDouble eintra, FloatOrDouble torsFreeEnergy, char  *prefixString, int ligand_is_inhibitor, FloatOrDouble unbound_internal_FE )
{
    FloatOrDouble deltaG = 0.0;
    FloatOrDouble Ki = 1.0;
    // FloatOrDouble RJ = 8.31441;  // in J/K/mol, Gas Constant, Atkins Phys.Chem., 2/e
    FloatOrDouble Rcal = 1.9871917; // in cal/K/mol, Gas Constant, RJ/4.184
    FloatOrDouble TK = 298.15;      // Room temperature, in K

    // equilibrium:   E  +  I  <=>    EI
    // binding:       E  +  I   ->    EI         K(binding),      Kb
    // dissociation:     EI     ->  E  +  I      K(dissociation), Kd
    //
    //                            1
    //         K(binding) = ---------------
    //                      K(dissociation)
    // so:
    //      ln K(binding) = -ln K(dissociation)
    //              ln Kb = -ln Kd
    // Ki = dissociation constant of the enzyme-inhibitor complex = Kd
    //      [E][I]
    // Ki = ------
    //       [EI]
    // so:
    //              ln Kb = -ln Ki
    // deltaG(binding)    = -R*T*ln Kb
    // deltaG(inhibition) =  R*T*ln Ki
    //
    // Binding and Inhibition occur in opposite directions, so we 
    // lose the minus-sign:  deltaG = R*T*lnKi,  _not_ -R*T*lnKi
    // => deltaG/(R*T) = lnKi
    // => Ki = exp(deltaG/(R*T))
    deltaG = einter + eintra + torsFreeEnergy - unbound_internal_FE;
    if (deltaG < 0.0) {
        Ki = exp((deltaG*1000.)/(Rcal*TK));
    }

    pr(stateFile, "\t\t<free_NRG_binding>");
    print1000(stateFile, deltaG);
    pr(stateFile, "</free_NRG_binding>\n");

    if (deltaG < 0.0) {
        if (ligand_is_inhibitor == 1) {
            pr(stateFile, "\t\t<Ki>");
            print1000_no_sign(stateFile, Ki);
            pr(stateFile, "</Ki>\n");
        } else {
            pr(stateFile, "\t\t<Kd>");
            print1000_no_sign(stateFile, Ki);
            pr(stateFile, "</Kd>\n");
        }
        pr(stateFile, "\t\t<Temp>%.2f</Temp>\n", TK); //temperature in K
    } 

    pr(stateFile, "\t\t<final_intermol_NRG>");
    print1000(stateFile, einter);
    pr(stateFile, "</final_intermol_NRG>\n");

    pr(stateFile, "\t\t<internal_ligand_NRG>");
    print1000(stateFile, eintra);
    pr(stateFile, "</internal_ligand_NRG>\n");

    pr(stateFile, "\t\t<torsonial_free_NRG>");
    print1000(stateFile, torsFreeEnergy);
    pr(stateFile, "</torsonial_free_NRG>\n"); 
}
