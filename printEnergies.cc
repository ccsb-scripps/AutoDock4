/*

 $Id: printEnergies.cc,v 1.3 2005/03/11 02:11:30 garrett Exp $

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

#define print1000(x) pr(logFile,  ((fabs((x)) >= 0.001) && ((fabs(x)) <= 1000.)) ? "%+7.2f" : "%+11.2e" , (x));

void print1000_no_sign(double x) {
    pr(logFile,  ((fabs((x)) >= 0.01) && ((fabs(x)) <= 1000.)) ? "%7.2f" : "%11.2e" , (x));
}

void printEnergies( FloatOrDouble einter, FloatOrDouble eintra, FloatOrDouble torsFreeEnergy, char  *prefixString, int ligand_is_inhibitor )
{
    FloatOrDouble deltaG = 0.0;
    FloatOrDouble Ki = 1.0;
    FloatOrDouble edocked=0.0;
    // FloatOrDouble RJ = 8.31441;  // in J/K/mol, Gas Constant, Atkins Phys.Chem., 2/e
    FloatOrDouble Rcal = 1.9871917; // in cal/K/mol, Gas Constant, RJ/4.184
    FloatOrDouble TK = 298.15;      // Room temperature, in K

    edocked = einter + eintra;
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
    deltaG = einter + torsFreeEnergy;
    if (deltaG < 0.0) {
        Ki = exp((deltaG*1000.)/(Rcal*TK));
    }

    pr( logFile, "%sEstimated Free Energy of Binding    = ", prefixString);
    print1000(deltaG);
    pr( logFile, " kcal/mol  [=(1)+(3)]\n");

    if (deltaG < 0.0) {
        if (ligand_is_inhibitor == 1) {
            pr( logFile, "%sEstimated Inhibition Constant, Ki   = ", prefixString);
        } else {
            pr( logFile, "%sEstimated Dissociation Constant, Kd = ", prefixString);
        }
        print1000_no_sign(Ki);
        pr( logFile, "       [Temperature = %.2f K]\n", TK);
    }

    pr( logFile, "%s\n", prefixString);

    pr( logFile, "%sFinal Docked Energy                 = ", prefixString);
    print1000(edocked);
    pr( logFile, " kcal/mol  [=(1)+(2)]\n");

    pr( logFile, "%s\n", prefixString);

    pr( logFile, "%s(1) Final Intermolecular Energy     = ", prefixString);
    print1000(einter);
    pr( logFile, " kcal/mol\n");

    pr( logFile, "%s(2) Final Internal Energy of Ligand = ", prefixString);
    print1000(eintra);
    pr( logFile, " kcal/mol\n");

    pr( logFile, "%s(3) Torsional Free Energy           = ", prefixString);
    print1000(torsFreeEnergy);
    pr( logFile, " kcal/mol\n");

    pr( logFile, "%s\n", prefixString);

}
